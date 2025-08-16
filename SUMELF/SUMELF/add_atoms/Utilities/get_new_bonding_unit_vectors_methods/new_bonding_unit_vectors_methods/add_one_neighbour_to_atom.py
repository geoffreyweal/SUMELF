"""
add_one_neighbour_to_atom.py, Geoffrey Weal, 5/7/2022

This script contains a method for determining the bonding vectors for adding one new hydrogen atoms to the central atom.
"""
from SUMELF.SUMELF.general_methods.geometry_methods import get_unit_vector

def add_one_neighbour_to_atom(bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, no_of_Hs_to_attach, number_of_lone_pairs_of_electrons, new_bonding_unit_vectors_recently_been_created, element_to_attach=None, logger=None):
	"""
	This method will return the bonding vectors for adding one new hydrogens to the atom.

	Note that:

				neighbours_to_add = no_of_Hs_to_attach + number_of_lone_pairs_of_electrons

				total_no_of_neighbours = no_of_Hs_to_attach + no_of_neighbouring_atoms_before_adding_hydrogens + number_of_lone_pairs_of_electrons # (this may exclude the number_of_lone_pairs_of_electrons value)
				

	The types of system that might use this method are:

	no_of_Hs_to_attach | no_of_neighbouring_atoms_before_adding_hydrogens | number_of_lone_pairs_of_electrons
	-------------------|--------------------------------------------------|----------------------------------
			1 		   |                         3                        |                 0                
			1 		   |                         2                        |                 0                
			1 		   |                         1                        |                 0                

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
	assert no_of_neighbouring_atoms_before_adding_hydrogens in [1, 2, 3]

	# First, get the central point across all bonding vectors
	centre_point = sum(bonding_unit_vectors)/len(bonding_unit_vectors)

	# Second, get the new_bonding_unit_vector, which is the unit vector of middle_point pointing in the opposite direction. 
	new_bonding_unit_vector = get_unit_vector(-centre_point)

	# Third, return new_bonding_unit_vector as a list. 
	return [new_bonding_unit_vector]




