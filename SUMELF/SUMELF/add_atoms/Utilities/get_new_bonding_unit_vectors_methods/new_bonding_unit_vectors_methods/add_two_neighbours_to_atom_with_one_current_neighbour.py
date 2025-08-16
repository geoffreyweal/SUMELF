"""
add_two_neighbours_to_atom_with_one_current_neighbour.py, Geoffrey Weal, 5/7/2022

This script contains a method for determining the bonding vectors for adding two new hydrogen atoms to the central atom that already contains one neighbour.
"""
from SUMELF.SUMELF.general_methods.geometry_methods                      import get_unit_vector
from SUMELF.SUMELF.general_methods.geometry_methods                      import planeFit
from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods import get_bond_angle
from SUMELF.SUMELF.general_methods.geometry_methods                      import get_rotation_matrix_around_arbitrary_axis

def add_two_neighbours_to_atom_with_one_current_neighbour(bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, no_of_Hs_to_attach, number_of_lone_pairs_of_electrons, new_bonding_unit_vectors_recently_been_created, element_to_attach=None, logger=None):
	"""
	This method will return the bonding vectors for adding two new hydrogens to the atom.

	Note that:

				neighbours_to_add = no_of_Hs_to_attach + number_of_lone_pairs_of_electrons

				total_no_of_neighbours = no_of_Hs_to_attach + no_of_neighbouring_atoms_before_adding_hydrogens + number_of_lone_pairs_of_electrons # (this may exclude the number_of_lone_pairs_of_electrons value)
				

	The types of system that might use this method are:

	no_of_Hs_to_attach | no_of_neighbouring_atoms_before_adding_hydrogens | number_of_lone_pairs_of_electrons
	-------------------|--------------------------------------------------|----------------------------------
			1 		   |                         1                        |                 1                
			2 		   |                         1                        |                 0                

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
	element_to_attach : str.
		This is the element you are using this method to attach a new atom to your molecule. If None given, their will be an issue given.

	Returns
	-------
	new_bonding_unit_vectors : list of numpy.array 
		These are the bonding unit vectors to add new hydrogens to in the current atom. 
	"""

	# Before beginning, obtain the number of atoms already surrounding the central atom.
	no_of_neighbouring_atoms_before_adding_hydrogens = len(bonding_unit_vectors)
	assert no_of_neighbouring_atoms_before_adding_hydrogens == 1

	# First, determine the index of the neighbouring atom attached to the central atom
	neighbouring_atom_index = list(molecule_graph[index_to_attach_Hs_to])[0]

	# Second, obtain the indices of atoms bound to the neighbouring atom
	neighbours_of_neighbouring_atom_indices = list(molecule_graph[neighbouring_atom_index])

	# Third, get all the positions of all the atoms in the plane of the likely double bond. 
	plane_positions = [molecule[neighbouring_atom_index].position] + [molecule[index].position for index in neighbours_of_neighbouring_atom_indices]

	# Fourth, obtain the plane for the atoms in the double bond and neighbours
	_, normal = planeFit(plane_positions)

	# Fifth, obtain the angle to rotate the hydrogen to add around the normal axis by.
	new_bonding_angle = get_bond_angle(atom_element, no_of_Hs_to_attach+no_of_neighbouring_atoms_before_adding_hydrogens, number_of_lone_pairs_of_electrons, element_to_attach=element_to_attach)
	if new_bonding_angle == 'flat':
		from ase.visualize import view
		view(molecule)
		print('index_to_attach_Hs_to = '+str(index_to_attach_Hs_to))
		import pdb; pdb.set_trace()
		raise Exception('Flat has come up. Check to make sure the angle should be right')

	# Sixth, obtain the matrix to rotate the hydrogen to add around the normal axis.
	rotation_matrix_around_neighbouring_bonding_unit_vector = get_rotation_matrix_around_arbitrary_axis(normal, new_bonding_angle)
	if no_of_Hs_to_attach == 2:
		negative_rotation_matrix_around_neighbouring_bonding_unit_vector = get_rotation_matrix_around_arbitrary_axis(normal, -new_bonding_angle)

	# Seventh, obtain the unit vector for the neighbouring atom
	neighbouring_bonding_unit_vector = get_unit_vector(bonding_unit_vectors[0])

	# Eighth, obtain the positions for bonding the two new hydrogen atoms too, which are obtained by rotating the neighbouring_bonding_unit_vector by 120 and 240 degree around the normal vector. 
	new_hydrogen_1_bond_vector = rotation_matrix_around_neighbouring_bonding_unit_vector @ neighbouring_bonding_unit_vector
	if no_of_Hs_to_attach == 2:
		new_hydrogen_2_bond_vector = negative_rotation_matrix_around_neighbouring_bonding_unit_vector @ neighbouring_bonding_unit_vector

	# Ninth, return the bonding vectors for the hydrogens, dont worry about returning the bonding vectors for the lone pairs of electrons.
	if   no_of_Hs_to_attach == 1: 
		return [get_unit_vector(new_hydrogen_1_bond_vector)]
	elif no_of_Hs_to_attach == 2: 
		return [get_unit_vector(new_hydrogen_1_bond_vector), get_unit_vector(new_hydrogen_2_bond_vector)]
	else:
		raise Exception('Huh?')






