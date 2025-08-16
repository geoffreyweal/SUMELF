"""
obtain_H_bonding_vectors.py, Geoffrey Weal, 6/7/2022

This script is designed to identify new H-bonding interactions.
"""
from math import pi
import numpy as np
import logging

from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods import get_bond_lengths
from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods import get_bond_angle
from SUMELF.SUMELF.general_methods.unit_cell_methods                     import get_cell_corner_points
from SUMELF.SUMELF.general_methods.angle_methods                         import get_angle
from SUMELF.SUMELF.general_methods.geometry_methods                      import get_unit_vector
from SUMELF.SUMELF.general_methods.distance_methods                      import get_distance
from SUMELF.SUMELF.general_methods.geometry_methods                      import get_rotation_matrix_around_arbitrary_axis
from SUMELF.SUMELF.general_methods.unit_cell_methods                     import centre_molecule_in_cell

hydrogen_bond_donor    = ['N', 'O', 'F'] # 'C',
hydrogen_bond_acceptor = ['N', 'O', 'F']
bond_lengths = get_bond_lengths('H')

def obtain_H_bonding_vectors(original_molecule, molecule_graph, index_to_attach_Hs_to, total_no_of_neighbours, no_of_pairs_of_lone_electrons, crystal, crystal_graph):
	"""
	This method is designed to identify any hydrogens that if created could form hydrogen bonding, and provide the bonding vector for forming hydrogen bonding.

	Parameters
	----------
	molecule : ase.Atoms
		This is the molecule you want to check if its hydrogens are involved in hydrogen-bonding.
	molecule_graph : networkx.Graph
		This is the graph of the molecule.
	index_to_attach_Hs_to : int 
		This is the index to check its neighbouring hydrogens for hydrogen-bonding
	total_no_of_neighbours : int
		This is the total number of atoms surrounding the donor atom.
	no_of_pairs_of_lone_electrons : int
		This is the number of pairs of lone electrons about the donor atom.
	crystal : ase.Atoms
		This is the full crystal object.
	crystal_graph : networkx.Graph
		This is the networkx graph for the full crystal.

	Returns
	-------
	hydrogen_bonding_atoms : dict
		This gives the information about the possible hydrogen bonding acceptors and their positions
	"""

	# First, obtain the element for the central atom.
	central_atom_element = original_molecule[index_to_attach_Hs_to].symbol

	# Second, only continue with this method if index_to_attach_Hs_to is a hydrogen_bond_donor atom.
	if central_atom_element not in hydrogen_bond_donor:
		return []

	# Third, move the molecule as close into the origin unit cell as possible.
	molecule = original_molecule.copy()
	molecule.set_cell(crystal.get_cell())
	molecule.set_pbc(True)
	centre_molecule_in_cell(molecule,crystal_cell_lattice=crystal.get_cell(),move_molecule=True)

	# Fourth, determine the positions of the the cell corners for the crystal
	cell_corner_points = get_cell_corner_points(crystal_cell_lattice=crystal.get_cell(),super_cell_reach=2)

	# Fifth, make a copy of the crystal and wrap it.Record the position of the wrapped crystal
	wrapped_crystal = crystal.copy()
	wrapped_crystal.wrap()
	crystal_wrapped_elements  = wrapped_crystal.get_chemical_symbols()
	crystal_wrapped_positions = wrapped_crystal.get_positions()
	crystal_wrapped_charges = wrapped_crystal.get_initial_charges()

	# Sixth, get the position for the central atom.
	central_atom_position = molecule[index_to_attach_Hs_to].position

	# Seventh, determine the bond length between the newly created H and the central atom
	bond_length_between_H_and_central_atom = bond_lengths[central_atom_element][total_no_of_neighbours]

	# Eighth, look through all the atoms in the crystal to see if the introduction of a hydrogen would be involved in hydrogen bonding. 
	hydrogen_bonding_details = []
	for acceptor_index in range(len(wrapped_crystal)):

		# 8.2.1: If the third atom is a hydrogen, move on to the next atom in wrapped_crystal
		if not is_H_bond_acceptor(acceptor_index, crystal_wrapped_elements, crystal_wrapped_charges):
			continue

		# 8.2.2: Obtain the position of this third atom. 
		original_acceptor_atom_position = crystal_wrapped_positions[acceptor_index]

		# 8.2.3: Try out all the positions from translations about the original unit cell by one unit cell length in each i,j,k direction. 
		for corner_position in cell_corner_points:

			# 8.2.3.1a: Get the acceptor_atom_position after translation.
			acceptor_atom_position = original_acceptor_atom_position + corner_position

			# 8.2.3.1b: Check that this acceptor_atom_position is not the same or about as central_atom_position. This likely indicates the two atoms at these positions are the same atom.
			distance_BETWEEN_acceptor_atom_position_AND_central_atom_position = get_distance(acceptor_atom_position, central_atom_position)
			if distance_BETWEEN_acceptor_atom_position_AND_central_atom_position < 0.1:
				continue

			# 8.2.3.2: Get the distance from the newly plced hydrogen to the acceptor atom (assuming linear direction).
			H_to_H_acceptor_distance = distance_BETWEEN_acceptor_atom_position_AND_central_atom_position - bond_length_between_H_and_central_atom
			# Make sure that the distance is not less than 0. Good sign that something has gone wrong.
			if H_to_H_acceptor_distance <= 0.0:
				import pdb; pdb.set_trace()
				raise Exception('Huh?')

			# 8.2.3.3: If the hydrogen and the H acceptor atom are within the max h_bnding distance, report this info (assume linear direction of bond from H-donor to H to H-acceptor).
			if H_to_H_acceptor_distance <= 3.5:

				# 8.2.3.3.1: Get the indices of the hydrogens bound to the H-bonding acceptor.
				neighbours = crystal_graph[acceptor_index]
				neighbours_that_are_hydrogens = [neighbour_index for neighbour_index in neighbours if (wrapped_crystal[neighbour_index].symbol == 'H')]
				neighbours_that_are_hydrogens = list(set(neighbours_that_are_hydrogens))

				# 8.2.3.3.2: Get the positions of the hydrogen bound to this H-bond acceptor.
				neighbours_that_are_hydrogens_positions = [(wrapped_crystal[hydrogen_index].position + corner_position) for hydrogen_index in neighbours_that_are_hydrogens]

				# 8.2.3.3.3: Add details to hydrogen_bonding_details
				hydrogen_bonding_details.append((acceptor_index, acceptor_atom_position, neighbours_that_are_hydrogens_positions))

	# Ninth, return the hydrogen bonding information
	return hydrogen_bonding_details

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

def is_H_bond_acceptor(acceptor_index, crystal_wrapped_elements, crystal_wrapped_charges):
	"""
	This method will determine if the atom of interest is a hydrogen bond acceptor

	Parameters
	----------
	acceptor_index : int
		This is the atom index for the atom that we are determining if it is a H-bond acceptor.
	crystal_wrapped_elements : list of str.
		These are the elements in the crystal.
	crystal_wrapped_charges : list of float
		These are the charges in the crystal.

	Returns
	-------
	True if the atom is a H-bond acceptor, False if not. 
	"""

	# First, determine if the element for this atom is a possible H-bonding element
	if crystal_wrapped_elements[acceptor_index] not in hydrogen_bond_acceptor:
		return False 

	# Second, determine if the atom has a positive, negative, or neutral charge
	if crystal_wrapped_charges[acceptor_index] < 0:
		return False

	# If you have got to this point, the atom is a H-bond acceptor.
	return True

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

