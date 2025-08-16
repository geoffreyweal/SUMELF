"""
supporting_methods.py, Geoffrey Weal, 5/7/2022

This script contains methods used by the python scripts in this folder.
"""
from math import pi
import numpy as np

from ase import Atom, Atoms
from ase.constraints import FixAtoms, FixBondLengths
from ase.optimize import FIRE

from SUMELF.SUMELF.general_methods.geometry_methods                                                                                                                                    import get_unit_vector
from SUMELF.SUMELF.general_methods.distance_methods                                                                                                                                    import get_distance
from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods                                                                                                               import get_bond_lengths
from SUMELF.SUMELF.general_methods.geometry_methods                                                                                                                                    import get_rotation_matrix_around_arbitrary_axis
from SUMELF.SUMELF.general_methods.unit_cell_methods                                                                                                                                   import get_cell_corner_points
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.get_new_bonding_vector_methods.get_best_bonding_vector_with_no_neighbours import get_best_bonding_vector_with_no_neighbours
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.get_new_bonding_vector_methods.get_best_bonding_vector_with_one_neighbour import get_best_bonding_vector_with_one_neighbour
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.interactions_methods.obtain_H_bonding_vectors                             import is_H_bond_acceptor
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.get_new_bonding_vector_methods.supporting_methods                         import project_u_onto_v, record_hydrogen_positions

from SUMELF.SUMELF.calculators.Coulomb import Coulomb

same_atom_distance_tolerance = 0.05
angle_tolerance = 15.0 * (pi/180.0)
def get_two_best_bonding_vector_with_no_neighbour(molecule, donor_index, total_no_of_neighbours, no_of_pairs_of_lone_electrons, wrapped_crystal, crystal_graph, acceptor_info, element_to_attach=None):
	"""
	This method will provide a optimal bonding unit vector for hydrogen bonding system that requires us to obtain two missing hydrogens, where the H-bond donor is bound to no atoms.

	Parameters
	----------
	molecule : ase.Atoms
		This is the molecule of interest in the crystal.
	donor_index : int
		This is the index of the H-bond donor
	total_no_of_neighbours : int
		This is the total number of neighbours that we want to surround the central atom
	no_of_pairs_of_lone_electrons : int
		This is the number of lone pairs of electrons surrounding the central atom.
	wrapped_crystal : ase.Atoms
		This is the crystal object.
	crystal_graph : networkx.Graph
		This is the networkx graph for the full crystal.
	acceptor_info : list of ints
		These are the indices of all the H-bond acceptors in the crystal, and their positions.
	element_to_attach : str.
		This is the element you are using this method to attach a new atom to your molecule. If None given, their will be an issue given.

	Returns
	-------
	new_bonding_unit_vector : numpy.array
		This is the newly created bonding unit vector
	"""

	# First, get the initial potential for the first unit vector
	new_initial_bonding_unit_vector_1 = get_best_bonding_vector_with_no_neighbours(molecule, donor_index, wrapped_crystal, acceptor_info)
	
	# Second, get the initial potential for the second unit vector
	new_initial_bonding_unit_vector_2 = get_best_bonding_vector_with_one_neighbour(molecule, donor_index, new_initial_bonding_unit_vector_1, total_no_of_neighbours, no_of_pairs_of_lone_electrons, wrapped_crystal, crystal_graph, acceptor_info, element_to_attach=element_to_attach)

	# Third, create the hydrogen bonding system
	H_bonding_system, bond_length = get_hydrogen_bonding_system(molecule, donor_index, total_no_of_neighbours, wrapped_crystal, crystal_graph, new_initial_bonding_unit_vector_1, new_initial_bonding_unit_vector_2, element_to_attach=element_to_attach)

	# Fourth, obtain the optimised bonding unit vectors.
	new_initial_bonding_unit_vector_1, new_initial_bonding_unit_vector_2 = optimise_hydrogen_bonding_system(H_bonding_system, bond_length)

	# Fifth, return the bonding unit vectors.
	return new_initial_bonding_unit_vector_1, new_initial_bonding_unit_vector_2

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_hydrogen_bonding_system(molecule, donor_index, total_no_of_neighbours, wrapped_crystal, crystal_graph, new_initial_bonding_unit_vector_1, new_initial_bonding_unit_vector_2, element_to_attach=None):
	"""
	This method will create the hydrogen bonding system for performing a coulombic local optimisation upon the hydrogens in the system.

	Parameters
	----------
	molecule : ase.Atoms
		This is the molecule of interest in the crystal.
	donor_index : int
		This is the index of the H-bond donor
	total_no_of_neighbours : int.
		This is the number of atoms surrounding the H-bond donor atom.
	wrapped_crystal : ase.Atoms
		This is the crystal object.
	crystal_graph : networkx.Graph
		This is the networkx graph for the full crystal.
	new_initial_bonding_unit_vector_1 : numpy.array
		This is the first bonding unit vector obtained using the get_best_bonding_vector_with_no_neighbours method. This is an initial guess position for this bonding vector.
	new_initial_bonding_unit_vector_2 : numpy.array
		This is the second bonding unit vector obtained using the get_best_bonding_vector_with_one_neighbour method. This is an initial guess position for this bonding vector.
	element_to_attach : str.
		This is the element you are using this method to attach a new atom to your molecule. If None given, their will be an issue given.

	Returns
	-------
	H_bonding_system : ase.Atoms
		This is the H-bonding system.
	bond_length : float
		This is the bond length between the H-bond donor and the hydrogen.

	"""

	# Before beginning, get the bond length betwen atoms and element_to_attach
	bond_lengths = get_bond_lengths(element_to_attach)

	# First, obtain the position, elements, and charges for the wrapped crystal.
	crystal_wrapped_elements  = wrapped_crystal.get_chemical_symbols()
	crystal_wrapped_positions = wrapped_crystal.get_positions()
	crystal_wrapped_charges = wrapped_crystal.get_initial_charges()

	# Second, obtain the element for the donor atom.
	donor_atom_element = molecule[donor_index].symbol

	# Third, get the position for the donor atom.
	donor_atom_position = molecule[donor_index].position
	donor_atom_charge = molecule[donor_index].charge

	# Fourth, determine the bond length between the newly created H and the donor atom
	bond_length_between_H_and_donor_atom = bond_lengths[donor_atom_element][total_no_of_neighbours]

	# Fifth, determine the positions of the the cell corners for the crystal
	cell_corner_points = get_cell_corner_points(crystal_cell_lattice=wrapped_crystal.get_cell(),super_cell_reach=2)

	# Sixth, create the H-bonding system.
	hydrogen_bonding_system = Atoms()
	for acceptor_index in range(len(wrapped_crystal)):

		# 6.1: If the third atom is a hydrogen, move on to the nex
		if not is_H_bond_acceptor(acceptor_index, crystal_wrapped_elements, crystal_wrapped_charges):
			continue

		# 6.2: obtain the position of this third atom. 
		original_acceptor_atom_position = crystal_wrapped_positions[acceptor_index]

		# 6.3: Try out all the positions from trnslations about the original unit cell by one unit cell length in each i,j,k direction. 
		for corner_position in cell_corner_points:

			# 6.3.1: Get the acceptor_atom_position after translation.
			acceptor_atom_position = original_acceptor_atom_position + corner_position

			# 6.3.2: Get the distance from the newly plced hydrogen to the acceptor atom (assuming linear direction).
			H_to_H_acceptor_distance = get_distance(acceptor_atom_position, donor_atom_position) - bond_length_between_H_and_donor_atom

			# 6.3.3: If the following is true, you are probably comparing the current atom to itself. Continue on.
			if H_to_H_acceptor_distance <= same_atom_distance_tolerance:
				continue

			# 6.3.4: If the hydrogen and the H acceptor atom are within the max h_bnding distance, report this info (assume linear direction of bond from H-donor to H to H-acceptor).
			if H_to_H_acceptor_distance <= 3.5:

				# 6.3.4.1: Get the indices of the hydrogens bound to the H-bonding acceptor.
				neighbours = crystal_graph[acceptor_index]
				neighbours_that_are_hydrogens = [neighbour_index for neighbour_index in neighbours if (wrapped_crystal[neighbour_index].symbol == 'H')]
				neighbours_that_are_hydrogens = list(set(neighbours_that_are_hydrogens))

				# 6.3.4.2: Get the positions of the hydrogen bound to this H-bond acceptor.
				neighbours_that_are_hydrogens_positions = [(wrapped_crystal[hydrogen_index].position + corner_position, wrapped_crystal[hydrogen_index].charge) for hydrogen_index in neighbours_that_are_hydrogens]

				# 6.3.4.2: Add the atom to the hydrogen bonding system.
				atom = Atom(symbol=wrapped_crystal[acceptor_index].symbol,position=acceptor_atom_position,charge=wrapped_crystal[acceptor_index].charge-0.5)
				hydrogen_bonding_system.append(atom)
				for neighbour_that_is_a_hydrogen_position, charge in neighbours_that_are_hydrogens_positions:
					atom = Atom(symbol='H',position=neighbour_that_is_a_hydrogen_position,charge=charge+1.0)
					hydrogen_bonding_system.append(atom)

	# ------------------------------------------------------------------------------------------------------------------------------
	# Seventh, create the hydrogen donor system

	# 7.1: Get the bond length for the H-donor --- H bond.
	bond_length = bond_lengths[donor_atom_element][total_no_of_neighbours]

	# 7.2, create the H-binding donor system
	hydrogen_donor_system = Atoms()

	# 7.3: Add the H-donor to it
	hydrogen_donor_system.append(Atom(symbol=donor_atom_element,position=donor_atom_position,charge=donor_atom_charge-0.5))

	# 7.4: Obtain the hydrogen postions and add them to the hydrogen_donor_system
	for new_initial_bonding_unit_vector in [new_initial_bonding_unit_vector_1, new_initial_bonding_unit_vector_2]:
		hydrogen_position = donor_atom_position + (bond_length*new_initial_bonding_unit_vector)
		hydrogen_donor_system.append(Atom(symbol='H',position=hydrogen_position,charge=+1.0))

	# ------------------------------------------------------------------------------------------------------------------------------
	# Eighth, create the full hydrogen bonding system that ontains both the hydrogen donor system and the hydrogen acceptors
	H_bonding_system = hydrogen_bonding_system.copy() + hydrogen_donor_system.copy()

	# Ninth, return the system
	return H_bonding_system, bond_length

def optimise_hydrogen_bonding_system(H_bonding_system, bond_length): 
	"""
	This method is designed to optimise the hydrogens in the hydrogen bonding system.

	Parameters
	----------
	H_bonding_system : ase.Atoms
		This is the H-bonding system.
	bond_length : float
		This is the bond length between the H-bond donor and the hydrogen.

	Returns
	-------
	system : ase.Atoms
		This is the newly created bonding unit vector.
	"""

	# First, determine the number of fixed atoms in the system
	no_of_stationary_atoms = len(H_bonding_system)-3

	# Second, fix all the atoms that are not the hydrogens we want to optimise.
	fixed_atoms = FixAtoms(indices=list(range(no_of_stationary_atoms+1)))

	# Third, provide the bonds that need to be fixed to preserve H-donor -- H bond lengths and H -- H-donor -- H angle.
	pairs = [(no_of_stationary_atoms,no_of_stationary_atoms+1), (no_of_stationary_atoms,no_of_stationary_atoms+2), (no_of_stationary_atoms+1,no_of_stationary_atoms+2)]
	fixed_bonds = FixBondLengths(pairs)

	# Fourth, contrain the required atoms, bond lengths, and bond angles. 
	H_bonding_system.set_constraint([fixed_bonds, fixed_atoms])

	# Fifth, set the coulombic calculator to the H-bonding system
	H_bonding_system.calc = Coulomb()

	# Sixth, perform an optimisation upon the 
	optimisation = FIRE(H_bonding_system, logfile=None) #trajectory='H2O.traj', logfile=None)
	optimisation.run(fmax=0.001, steps=100)

	# Seventh, obtain the positions of the H-bond donor and the hydrogens bound to the H-bond donor after optimisation.
	H_bond_donor = H_bonding_system[-3].position
	H1_donor     = H_bonding_system[-2].position
	H2_donor     = H_bonding_system[-1].position

	# Eighth, obtain the new bonding unit vector after optimisation.
	new_initial_bonding_unit_vector_1 = get_unit_vector(H1_donor - H_bond_donor)
	new_initial_bonding_unit_vector_2 = get_unit_vector(H2_donor - H_bond_donor)

	# Ninth, return new_initial_bonding_unit_vector_1 and new_initial_bonding_unit_vector_2
	return new_initial_bonding_unit_vector_1, new_initial_bonding_unit_vector_2

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------







