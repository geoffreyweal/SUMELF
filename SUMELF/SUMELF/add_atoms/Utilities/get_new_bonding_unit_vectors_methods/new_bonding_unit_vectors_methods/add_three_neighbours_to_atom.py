"""
add_three_neighbours_to_atom.py, Geoffrey Weal, 5/7/2022

This script contains a method for determining the bonding vectors for adding three new hydrogen atoms to the central atom.
"""
import numpy as np
from math import pi
from random import uniform

from SUMELF.SUMELF.general_methods.geometry_methods                      import get_unit_vector
from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods import get_bond_angle
from SUMELF.SUMELF.general_methods.geometry_methods                      import get_rotation_matrix_around_arbitrary_axis

from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.add_two_neighbours_to_atom_with_one_current_neighbour  import add_two_neighbours_to_atom_with_one_current_neighbour
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.add_two_neighbours_to_atom_with_two_current_neighbours import add_two_neighbours_to_atom_with_two_current_neighbours
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.get_new_bonding_vector                                 import get_new_bonding_vector

minimal_angle = 30.0 * (pi/180.0)
def add_three_neighbours_to_atom(bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, no_of_Hs_to_attach, number_of_lone_pairs_of_electrons, crystal, crystal_graph, new_bonding_unit_vectors_recently_been_created, element_to_attach=None, logger=None):
	"""
	This method will return the bonding vectors for adding three new hydrogens to the atom.

	Note that:

				neighbours_to_add = no_of_Hs_to_attach + number_of_lone_pairs_of_electrons

				total_no_of_neighbours = no_of_Hs_to_attach + no_of_neighbouring_atoms_before_adding_hydrogens + number_of_lone_pairs_of_electrons # (this may exclude the number_of_lone_pairs_of_electrons value)
				

	The types of system that might use this method are:

	no_of_Hs_to_attach | no_of_neighbouring_atoms_before_adding_hydrogens | number_of_lone_pairs_of_electrons
	-------------------|--------------------------------------------------|----------------------------------
			1 		   |                         1                        |                 2                
			2 		   |                         1                        |                 1                
			2 		   |                         0                        |                 1                
			3 		   |                         1                        |                 0                
			3 		   |                         0                        |                 0                

	Parameters
	----------
	bonding_unit_vectors : list of numpy.array 
		These are the bonding unit vectors for the current neighbours that are already bonded to the atom of focus.
	molecule : ase.Atoms
		This is the molecule you want to remove the aliphatic carbons to.
	molecule_graph : networkx.Graph
		This is the graph of this molecule.
	index_to_attach_Hs_to : int.
		This is the index that you want to add hydrogens to.
	atom_element : str. 
		This is the element for the central atom.
	no_of_Hs_to_attach : int.
		This is the number of hydrogens that will be attached to the central atom (number of new bonding unit vectors to obtain).
	number_of_lone_pairs_of_electrons : int.
		These are the number of pairs of lone electrons about the central atom.
	crystal : ase.Atoms
		This is the ase.Atoms object of the crystal.
	crystal_graph : networkx.graph
		This is the graph of the associated crystal
	element_to_attach : str.
		This is the element you are using this method to attach a new atom to your molecule. If None given, their will be an issue given.

	Returns
	-------
	new_bonding_unit_vectors : list of numpy.array 
		These are the bonding unit vectors to add new hydrogens to in the current atom. 
	"""

	# Before beginning, obtain the number of atoms already surrounding the central atom.
	no_of_neighbouring_atoms_before_adding_hydrogens = len(bonding_unit_vectors)
	assert no_of_neighbouring_atoms_before_adding_hydrogens in [0, 1]

	# After this, obtain the total number of neighbours, as well as the data required from the data_for_getting_new_bonding_vector list.
	total_no_of_neighbours = no_of_neighbouring_atoms_before_adding_hydrogens + no_of_Hs_to_attach

	if no_of_neighbouring_atoms_before_adding_hydrogens == 0:
		# First, if there are no bonding_unit_vectors for this atom, perform this task:

		raise Exception('checked add_three_neighbours_to_atom. need to test practically.')

		# 1.1: Place the first hydrogen in a random place.
		new_hydrogen_1_bond_vector = get_unit_vector(next(get_new_bonding_vector(molecule, molecule_graph, index_to_attach_Hs_to, total_no_of_neighbours, number_of_lone_pairs_of_electrons, crystal, crystal_graph, default=(0,0,1), new_bonding_unit_vectors_recently_been_created=new_bonding_unit_vectors_recently_been_created, element_to_attach=element_to_attach, logger=logger)))
		
		# 1.2: We have obtained one new hydrogen bonding vector, reduce the no_of_Hs_to_attach by one
		new_no_of_Hs_to_attach = no_of_Hs_to_attach - 1

		# 1.3: Only perform the next step if you have more hydrogens to attach to molecule (new_no_of_Hs_to_attach >= 1)
		if new_no_of_Hs_to_attach >= 1:

			# 1.3.1: Obtain positions for the other two neighbours (be them atoms or a lone pair of electrons)
			new_hydrogen_bond_vectors = add_two_neighbours_to_atom_with_one_current_neighbour([new_hydrogen_1_bond_vector], molecule, molecule_graph, index_to_attach_Hs_to, atom_element, new_no_of_Hs_to_attach, number_of_lone_pairs_of_electrons, new_bonding_unit_vectors_recently_been_created + [new_hydrogen_1_bond_vector], element_to_attach=element_to_attach, logger=logger)

	elif no_of_neighbouring_atoms_before_adding_hydrogens == 1:
		# Second, if there is one bonding vector in the bonding_unit_vectors list for this atom, perform this task:

		# 2.1: Get the single bonding_unit_vector from bonding_unit_vectors
		bonding_unit_vector = bonding_unit_vectors[0]

		# 2.2: Get a randomly generated perpendicular vector to the line from the origin to bonding_unit_vector
		for a_vector in get_new_bonding_vector(molecule, molecule_graph, index_to_attach_Hs_to, total_no_of_neighbours, number_of_lone_pairs_of_electrons, crystal, crystal_graph, default=[(1,0,0),(0,1,0)], new_bonding_unit_vectors_recently_been_created=new_bonding_unit_vectors_recently_been_created, element_to_attach=element_to_attach, logger=logger):
			a_unit_vector = get_unit_vector(a_vector)
			if np.arccos(np.dot(bonding_unit_vector,a_unit_vector)) > minimal_angle:
				perpendicular_vector = get_unit_vector(np.cross(bonding_unit_vector,a_unit_vector))
				break

		# 2.3: Obtain the matrix for rotation 120 degree around the perpendicular_vector axis
		new_bonding_angle = get_bond_angle(atom_element, no_of_Hs_to_attach+no_of_neighbouring_atoms_before_adding_hydrogens, number_of_lone_pairs_of_electrons, element_to_attach=element_to_attach)
		rotation_matrix_around_perpendicular_vector = get_rotation_matrix_around_arbitrary_axis(perpendicular_vector, new_bonding_angle)

		# 2.4: Rotate bonding_unit_vector around perpendicular_vector by 120 degrees
		new_hydrogen_1_bond_vector = rotation_matrix_around_perpendicular_vector @ bonding_unit_vector

		# 2.5: We have obtained one new hydrogen bonding vector, reduce the no_of_Hs_to_attach by one
		new_no_of_Hs_to_attach = no_of_Hs_to_attach - 1

		# Only perform the next step if you have more hydrogens to attach to molecule (new_no_of_Hs_to_attach >= 1).
		# * If you are at this point, you have two atoms bonded to the central atom, and want to add two new atoms, 
		#   or a new atom and a lone pair of electrons. 
		if new_no_of_Hs_to_attach >= 1:

			if not new_no_of_Hs_to_attach == 2:
				raise Exception('Error: new_no_of_Hs_to_attach shoudl probably be equal to 2. new_no_of_Hs_to_attach: '+str(new_no_of_Hs_to_attach))

			# 2.6: Get the two points that are equally spaced between bonding_unit_vector and new_hydrogen_1_bond_vector
			current_set_of_bonding_vectors = [bonding_unit_vector, new_hydrogen_1_bond_vector]
			new_hydrogen_bond_vectors = add_two_neighbours_to_atom_with_two_current_neighbours(current_set_of_bonding_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, new_no_of_Hs_to_attach, number_of_lone_pairs_of_electrons, new_bonding_unit_vectors_recently_been_created + [new_hydrogen_1_bond_vector], element_to_attach=element_to_attach, logger=logger)

	else:
		# Third, these method are only designed to have 4 atoms around an atom at most. 
		raise Exception('Huh')

	# Fourth, extract the bonding vectors from the results of add_two_neighbours_to_atom_with_one_current_neighbour/add_two_neighbours_to_atom_with_one_current_neighbour
	if   new_no_of_Hs_to_attach == 0:	
		pass
	elif new_no_of_Hs_to_attach == 1:	
		new_hydrogen_2_bond_vector, = new_hydrogen_bond_vectors
	elif new_no_of_Hs_to_attach == 2:	
		new_hydrogen_2_bond_vector, new_hydrogen_3_bond_vector = new_hydrogen_bond_vectors
	else:
		raise Exception('Huh')

	# Fifth, return the bonding vectors for the hydrogens, dont worry about returning the bonding vectors for the lone pairs of electrons.
	if   no_of_Hs_to_attach == 1: 
		return [get_unit_vector(new_hydrogen_1_bond_vector)]
	elif no_of_Hs_to_attach == 2: 
		return [get_unit_vector(new_hydrogen_1_bond_vector), get_unit_vector(new_hydrogen_2_bond_vector)]
	elif no_of_Hs_to_attach == 3: 
		return [get_unit_vector(new_hydrogen_1_bond_vector), get_unit_vector(new_hydrogen_2_bond_vector), get_unit_vector(new_hydrogen_3_bond_vector)]
	else:
		raise Exception('Huh?')






