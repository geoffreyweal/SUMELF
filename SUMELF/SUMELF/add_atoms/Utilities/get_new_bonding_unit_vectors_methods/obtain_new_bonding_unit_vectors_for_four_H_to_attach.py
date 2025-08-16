"""
obtain_new_bonding_unit_vectors_for_four_H_to_attach.py, Geoffrey Weal, 4/7/2022

This script is designed to obtain the new bonding_unit_vectors for adding four hydrogen atoms about an atom. 
"""
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.add_four_neighbours_to_atom import add_four_neighbours_to_atom

def obtain_new_bonding_unit_vectors_for_four_H_to_attach(bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, number_of_lone_pairs_of_electrons, crystal, crystal_graph, element_to_attach=None, logger=None):
	"""
	This method is designed to obtain the new bonding_unit_vectors for adding four hydrogen atoms about an atom. 

	Note that the methods used are based on the equation:

				neighbours = no_of_neighbouring_atoms_before_adding_hydrogens + no_of_Hs_to_attach + number_of_lone_pairs_of_electrons

	Parameters
	----------
	bonding_unit_vectors : list of numpy.array 
		These are the bonding unit vectors for the current neighbours that are already bonded to the atom of focus.
	molecule : ase.Atoms
		This is the molecule you want to remove the aliphatic carbons to.
	molecule_graph : networkx.Graph
		This is the graph of this molecule.
	index_to_attach_Hs_to :list 
		This is the atom to attach new hydrogens to.
	atom_element : str.
		This is the element of the atom we will be adding hydrogens to.
	number_of_lone_pairs_of_electrons : int
		This is the number of pairs of lone electrons in the atom we want to attach hydrogen about. 
	crystal : ase.Atoms
		This is the ase.Atoms object of the crystal.
	crystal_graph : networkx.graph
		This is the graph of the associated crystal
	element_to_attach : str.
		This is the element you are using this method to attach a new atom to your molecule. If None given, their will be an issue given.
	logger : logging
		This is the logger to write information and warnings about the running of the program to.

	Returns
	-------
	new_bonding_unit_vectors : list of numpy.array 
		These are the bonding unit vectors to add four new hydrogens to in the current atom. 
	"""

	# Before beginning, note the number of new hydrogens to add to the system
	no_of_Hs_to_attach = 4

	# First, determine the number of neighbouring atoms about the central atom before adding hydrogens.
	no_of_neighbouring_atoms_before_adding_hydrogens = len(bonding_unit_vectors)

	# Second, make sure that no_of_neighbouring_atoms_before_adding_hydrogens == 0 and number_of_lone_pairs_of_electrons == 0, 
	# as we are only dealing with atoms that can have up to four atoms about them.
	if not (no_of_neighbouring_atoms_before_adding_hydrogens == 0 and number_of_lone_pairs_of_electrons == 0):
		raise Exception('Warning: Atom is a bit weird. no_of_neighbouring_atoms_before_adding_hydrogens = '+str(no_of_neighbouring_atoms_before_adding_hydrogens)+', number_of_lone_pairs_of_electrons = '+str(number_of_lone_pairs_of_electrons))

	# Third, there are four hydrogens being added to an atom with two other atoms around it. We need the atom to have tetrahedral geometry. 
	new_bonding_unit_vectors = add_four_neighbours_to_atom(bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, no_of_Hs_to_attach, number_of_lone_pairs_of_electrons, crystal, crystal_graph, [], element_to_attach=element_to_attach)
	
	# Fourth, check that four new_bonding_unit_vector is given in the new_bonding_unit_vectors list
	if not (isinstance(new_bonding_unit_vectors,list) and (len(new_bonding_unit_vectors) == 4)):
		raise Exception('Error: Only one new_bonding_unit_vector should be given in the new_bonding_unit_vectors list. len(new_bonding_unit_vectors) = '+str(len(new_bonding_unit_vectors)))

	# Fifth, return the new bonding unit vectors for placing a hydrogen around the central atom
	return new_bonding_unit_vectors



