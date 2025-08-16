"""
get_new_bonding_unit_vectors.py, Geoffrey Weal, 4/7/2022

This script is designed to determine the new bonding unit vectors for adding hydrogens about atoms. 
"""
from SUMELF.SUMELF.general_methods.general_molecules_methods                                                                      import get_number_of_lone_pairs_of_electron_pairs
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.obtain_new_bonding_unit_vectors_for_one_H_to_attach   import obtain_new_bonding_unit_vectors_for_one_H_to_attach
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.obtain_new_bonding_unit_vectors_for_two_H_to_attach   import obtain_new_bonding_unit_vectors_for_two_H_to_attach
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.obtain_new_bonding_unit_vectors_for_three_H_to_attach import obtain_new_bonding_unit_vectors_for_three_H_to_attach
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.obtain_new_bonding_unit_vectors_for_four_H_to_attach  import obtain_new_bonding_unit_vectors_for_four_H_to_attach

Atoms_set_up_for_adding_hydrogens_to = ['B', 'C', 'N', 'O', 'S']
def get_new_bonding_unit_vectors(no_of_Hs_to_attach, bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, number_of_lone_pairs_of_electrons, crystal, crystal_graph, element_to_attach=None, logger=None):
	"""
	This method will return the bonding vectors for adding new hydrogens to the atom.

	Parameters
	----------
	no_of_Hs_to_attach : int
		This is the number of hydrogens we want to obtain bonding unit vectors for. These could be new hydrogens being added, or vectors to update the positions of existing hydrogen atoms. 
	bonding_unit_vectors : list of numpy.array 
		These are the bonding unit vectors for the current neighbours that are already bonded to the atom of focus.
	molecule : ase.Atoms
		This is the molecule you want to remove the aliphatic carbons to.
	molecule_graph : networkx.Graph
		This is the graph of this molecule.
	index_to_attach_Hs_to :list 
		This is the atom to attach new hydrogens to.
	number_of_lone_pairs_of_electrons : int.
		This is the number of lone pairs of electrons about the molecule.
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
		These are the bonding unit vectors to add new hydrogens to in the current atom. 
	"""

	# First, obtain the element for this atom.
	atom_element = molecule[index_to_attach_Hs_to].symbol

	# Second, this method has only been designed for for certain elements:
	if atom_element not in Atoms_set_up_for_adding_hydrogens_to:
		to_string  = 'Issue with adding hydrogens to molecule. This program has only been designed to add hydrogens to: '+str(Atoms_set_up_for_adding_hydrogens_to)+'\n'
		to_string += 'The atom you are trying to add hydrogens to is '+str(atom_element)+'\n'
		to_string += 'Check this out.'
		import pdb; pdb.set_trace()
		raise Exception(to_string)

	# Third, obtain the new bonding unit vectors for attaching new hydrogens to the atom of interest.
	if   no_of_Hs_to_attach == 0:
		import pdb; pdb.set_trace()
		raise Exception('Warning: You are not adding any hydrogens to an atom in this crystal.')
	elif no_of_Hs_to_attach == 1:
		new_bonding_unit_vectors = obtain_new_bonding_unit_vectors_for_one_H_to_attach  (bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, number_of_lone_pairs_of_electrons, crystal, crystal_graph, element_to_attach=element_to_attach, logger=logger)
	elif no_of_Hs_to_attach == 2:
		new_bonding_unit_vectors = obtain_new_bonding_unit_vectors_for_two_H_to_attach  (bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, number_of_lone_pairs_of_electrons, crystal, crystal_graph, element_to_attach=element_to_attach, logger=logger)
	elif no_of_Hs_to_attach == 3:
		new_bonding_unit_vectors = obtain_new_bonding_unit_vectors_for_three_H_to_attach(bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, number_of_lone_pairs_of_electrons, crystal, crystal_graph, element_to_attach=element_to_attach, logger=logger)
	elif no_of_Hs_to_attach == 4:
		new_bonding_unit_vectors = obtain_new_bonding_unit_vectors_for_four_H_to_attach (bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, number_of_lone_pairs_of_electrons, crystal, crystal_graph, element_to_attach=element_to_attach, logger=logger)
	else:
		raise Exception('Warning: This method has only been designed to attach up to 4 hydrogens to atoms to date.')

	# Fourth, check that we have obtain the expected number of new bonding unit vectors.
	#         * This should be the same as no_of_Hs_to_attach
	if not (len(new_bonding_unit_vectors) == no_of_Hs_to_attach):
		toString  = 'Error: The number of new bonding unit vectors is not the same as no_of_Hs_to_attach.\n'
		toString += 'We expect these two values to be the same.\n'
		toString += 'no_of_Hs_to_attach = '+str(no_of_Hs_to_attach)+'\n'
		toString += 'len(new_bonding_unit_vectors) = '+str(len(new_bonding_unit_vectors))+'\n'
		toString += 'new_bonding_unit_vectors = '+str(new_bonding_unit_vectors)+'\n'
		toString += 'Check this, as this might indicate a programming error.'
		raise Exception(toString)

	# Fifth, return new_bonding_unit_vectors
	return new_bonding_unit_vectors

















'''








	indices_of_neighbouring_hydrogens_in_original_molecule = tuple(index for index in molecule_graph[index_to_attach_Hs_to] if (molecule[index].symbol == 'H'))

	indices_of_hydrogens_that_will_be_modified = [index for index, position in bonding_unit_vectors if index in indices_of_neighbouring_hydrogens_in_original_molecule]




	import pdb; pdb.set_trace()


	molecule_graph[index_to_attach_Hs_to]



	# First, determine the number of atoms neighbouring atom index_to_attach_Hs_to after hydrogen have been attached
	total_no_of_neighbours_after_adding_H_atoms = len(bonding_unit_vectors) + total_no_of_Hs_to_attach




	# Second, determine the number of hydrogens currently attached to this atom currently.
	no_of_neighbouring_H_atoms_before_adding_hydrogens = len([index for index in molecule_graph[index_to_attach_Hs_to] if (molecule[index].symbol == 'H')])

	# Third, make sure that total_no_of_Hs_to_attach > no_of_neighbouring_H_atoms_before_adding_hydrogens.
	if not (total_no_of_Hs_to_attach >= no_of_neighbouring_H_atoms_before_adding_hydrogens):
		raise Exception('Issue: It is required that total_no_of_Hs_to_attach >= no_of_neighbouring_H_atoms_before_adding_hydrogens')

	# Fourth, determine the number of hydrogens to attach to this atom.
	no_of_Hs_to_attach = total_no_of_Hs_to_attach - no_of_neighbouring_H_atoms_before_adding_hydrogens

	# Eighth, make sure that hydrogens are being added to the system, otherwise why are we using this method.
	if no_of_Hs_to_attach == 0:
		import pdb; pdb.set_trace()
		raise Exception('Warning: You are not adding any hydrogens to an atom in this crystal.')

	# Sixth, determine how many neighbours your atom will have after you attach hydrogens to it.
	total_no_of_neighbours_after_adding_H_atoms = no_of_neighbouring_atoms_before_adding_hydrogens + no_of_Hs_to_attach

	# Seventh, this method will only add hydrogens to atoms so long as that atom will have less than or equal to 4 neighbouring atoms in total.
	if total_no_of_neighbours_after_adding_H_atoms > 4:
		from ase.visualize import view
		view(molecule)
		print('Can not add hydrogens to an atoms that will have more than 4 atoms around it in total. index: '+str(index_to_attach_Hs_to))
		print(index_to_attach_Hs_to, total_no_of_Hs_to_attach, hybridisation, formal_charge)
		import pdb; pdb.set_trace()
		exit()



	import pdb; pdb.set_trace()

	# Eighth, determine the number of valence electrons attached to the index_to_attach_Hs_to atom.
	number_of_lone_pairs_of_electrons = get_number_of_lone_pairs_of_electron_pairs(hybridisation, total_no_of_neighbours_after_adding_H_atoms, formal_charge, index_to_attach_Hs_to, molecule, molecule_graph, logger=logger)

	# Ninth, obtain the new_bonding_unit_vectors
	data_for_getting_new_bonding_vector = [crystal, crystal_graph]
	new_bonding_unit_vectors = obtain_new_bonding_unit_vectors(bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, no_of_Hs_to_attach, number_of_lone_pairs_of_electrons, data_for_getting_new_bonding_vector, element_to_attach=element_to_attach, logger=logger)

	# Tenth, return the new_bonding_unit_vectors
	return new_bonding_unit_vectors

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------

def obtain_new_bonding_unit_vectors(bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, no_of_Hs_to_attach, number_of_lone_pairs_of_electrons, data_for_getting_new_bonding_vector, element_to_attach=None, logger=None):
	"""
	This method is designed to obtain the new bonding_unit_vectors for adding hydrogen atoms about an atom. 

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
	no_of_Hs_to_attach : int
		This is the number of hydrogens that we want to attach to this atom
	number_of_lone_pairs_of_electrons : int
		This is the number of pairs of lone electrons in the atom we want to attach hydrogen about. 
	data_for_getting_new_bonding_vector : list
		This is a collection of information needed for obtaining the new bonding vector: [crystal, crystal_graph, hybridisation_of_atoms_in_molecule] 
	element_to_attach : str.
		This is the element you are using this method to attach a new atom to your molecule. If None given, their will be an issue given.
	logger : logging
		This is the logger to write information and warnings about the running of the program to.
	
	Returns
	-------
	new_bonding_unit_vectors : list of numpy.array 
		These are the bonding unit vectors to add new hydrogens to in the current atom. 
	"""

	# First, obtain the element for the atom of interest
	atom_element = molecule[index_to_attach_Hs_to].symbol

	data_for_getting_new_bonding_vector = [crystal, crystal_graph]

	# Second, obtain the new bonding unit vectors for attaching new hydrogens to the atom of interest.
	if no_of_Hs_to_attach == 0:
		import pdb; pdb.set_trace()
		raise Exception('Warning: You are not adding any hydrogens to an atom in this crystal.')
	if no_of_Hs_to_attach == 1:
		new_bonding_unit_vectors = obtain_new_bonding_unit_vectors_for_one_H_to_attach  (bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, number_of_lone_pairs_of_electrons, data_for_getting_new_bonding_vector, element_to_attach=element_to_attach, logger=logger)
	elif no_of_Hs_to_attach == 2:
		new_bonding_unit_vectors = obtain_new_bonding_unit_vectors_for_two_H_to_attach  (bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, number_of_lone_pairs_of_electrons, data_for_getting_new_bonding_vector, element_to_attach=element_to_attach, logger=logger)
	elif no_of_Hs_to_attach == 3:
		new_bonding_unit_vectors = obtain_new_bonding_unit_vectors_for_three_H_to_attach(bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, number_of_lone_pairs_of_electrons, data_for_getting_new_bonding_vector, element_to_attach=element_to_attach, logger=logger)
	elif no_of_Hs_to_attach == 4:
		new_bonding_unit_vectors = obtain_new_bonding_unit_vectors_for_four_H_to_attach (bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, number_of_lone_pairs_of_electrons, data_for_getting_new_bonding_vector, element_to_attach=element_to_attach, logger=logger)
	else:
		raise Exception('Warning: This method has only been designed to attach up to 4 hydrogens to atoms to date.')

	# Third, return new_bonding_unit_vectors
	return new_bonding_unit_vectors

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------


'''

