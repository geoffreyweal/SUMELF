"""
supporting_methods.py, Geoffrey Weal, 5/7/2022

This script contains methods used by the python scripts in this folder.
"""
import numpy as np
import logging

from ase.data import atomic_numbers

from SUMELF.SUMELF.general_methods.geometry_methods                                                                                                                                       import get_unit_vector
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.interactions_methods.obtain_H_bonding_vectors                                import obtain_H_bonding_vectors
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.interactions_methods.obtain_H_pi_interaction_vectors                         import obtain_H_pi_interaction_vectors
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.interactions_methods.obtain_ideal_methyl_like_hydrogen_positions             import obtain_ideal_methyl_like_hydrogen_positions

from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.get_new_bonding_vector_methods.get_best_bonding_vector_with_no_neighbours    import get_best_bonding_vector_with_no_neighbours
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.get_new_bonding_vector_methods.get_best_bonding_vector_with_one_neighbour    import get_best_bonding_vector_with_one_neighbour
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.get_new_bonding_vector_methods.get_two_best_bonding_vector_with_no_neighbour import get_two_best_bonding_vector_with_no_neighbour

def get_new_bonding_vector(molecule, molecule_graph, index_to_attach_Hs_to, total_no_of_neighbours, no_of_pairs_of_lone_electrons, crystal, crystal_graph, default=None, new_bonding_unit_vectors_recently_been_created=[], element_to_attach=None, logger=None):
	"""
	This method is designed to determine a new bonding unit vector for adding the next hydrogen to the molecule. 

	Parameters
	----------
	molecule : ase.Atoms
		This is the molecule you want to check if its hydrogens are involved in hydrogen-bonding.
	molecule_graph : networkx.Graph
		This is the graph of the molecule.
	index_to_attach_Hs_to : int.
		This is the index to check its neighbouring hydrogens for hydrogen-bonding
	total_no_of_neighbours : list of int.
		This is the total number of atoms about the index_to_attach_Hs_to atom.
	no_of_pairs_of_lone_electrons : list of int.
		This is the total number of pairs of lone electrons about the index_to_attach_Hs_to atom.
	new_bonding_unit_vectors_recently_been_created : list of numpy.arrays
		These are the new_bonding_unit_vectors that have been created during the current process of obtaining the new bonding vectors.
	crystal ; ase.Atoms
		This is the full crystal object.
	crystal_graph : networkx.Graph
		This is the networkx graph for the full crystal.
	default : list, numpy.array, None
		These are all the vectors to try first if no hydrogen bonding or H-pi aromatic ring interactions are found. 
	element_to_attach : str.
		This is the element you are using this method to attach a new atom to your molecule. If None given, their will be an issue given.

	Yields
	------
	A unit vector : numpy.array
	"""

	# First, obtain all the hydrogen bonding information if index_to_attach_Hs_to is set as a H-bonding donor.
	if element_to_attach == 'H':
		hydrogen_bonding_details = obtain_H_bonding_vectors(molecule, molecule_graph, index_to_attach_Hs_to, total_no_of_neighbours, no_of_pairs_of_lone_electrons, crystal, crystal_graph)
	else:
		hydrogen_bonding_details = []

	# Second, wrap the crystal in order to make sure you can determine the position of atoms in the origin unit cell.
	wrapped_crystal = crystal.copy()
	wrapped_crystal.wrap()

	# ------------------------------------------------------------------------------------------------------------------------------

	# Third, if their are hydrogen_bonding_details or H_pi_interaction_details found, use these to obtain the new_bonding_vector.
	if not len(hydrogen_bonding_details) == 0:

		# 3.1: Write data to the logger if a logger is given.
		if logger is not None:
			logger.info('Hydrogen: '+str(hydrogen_bonding_details))

		# 3.2: Obtain the indices of the atoms that neighbour index_to_attach_Hs_to
		neighbour_indices = list(molecule_graph[index_to_attach_Hs_to])

		# 3.3: Obtain the total number of neighbours, including new bonding unit vectors
		recently_made_new_bonding_unit_vectors = list(new_bonding_unit_vectors_recently_been_created)
		total_number_of_neighbour_indices_and_new_bonding_unit_vectors = len(neighbour_indices+recently_made_new_bonding_unit_vectors)

		# 3.4: Obtain the new bonding unit vectors.
		if   total_number_of_neighbour_indices_and_new_bonding_unit_vectors == 0:
			if logger is not None:
				logger.info('The get_best_bonding_vector_with_no_neighbours method was used to obtain H-bonding vectors')
			yield get_best_bonding_vector_with_no_neighbours(molecule, index_to_attach_Hs_to, wrapped_crystal, hydrogen_bonding_details, element_to_attach=element_to_attach)
		elif total_number_of_neighbour_indices_and_new_bonding_unit_vectors == 1:
			if   len(recently_made_new_bonding_unit_vectors) == 0:
				if logger is not None:
					logger.info('The get_best_bonding_vector_with_one_neighbour method was used to obtain H-bonding vectors')
				yield get_best_bonding_vector_with_one_neighbour(molecule, index_to_attach_Hs_to, neighbour_indices[0], total_no_of_neighbours, no_of_pairs_of_lone_electrons, wrapped_crystal, crystal_graph, hydrogen_bonding_details, element_to_attach=element_to_attach)
			elif len(recently_made_new_bonding_unit_vectors) == 1:
				if logger is not None:
					logger.info('The get_two_best_bonding_vector_with_no_neighbour method was used to obtain H-bonding vectors')
				new_bonding_unit_vectors = get_two_best_bonding_vector_with_no_neighbour(molecule, index_to_attach_Hs_to, total_no_of_neighbours, no_of_pairs_of_lone_electrons, wrapped_crystal, crystal_graph, hydrogen_bonding_details, element_to_attach=element_to_attach)
				for new_bonding_unit_vector in new_bonding_unit_vectors:
					yield new_bonding_unit_vector
			else:
				raise Exception('No method has been set up for this point.')
		else:
			raise Exception('No method has been set up for this point.')

	# ------------------------------------------------------------------------------------------------------------------------------
	# Fourth, obtain all the H-pi aromatic ring interaction details that would involve a hydrogen bound to the index_to_attach_Hs_to atom.
	#         * If their is competition between whether a hydrogen bond or has a H-pi aromatic ring, we assume the  hydrogen bond is most dominant.

	# 4.1: Obtain the vectors for a hydrogen being near a pi bond
	new_bonding_unit_vectors_from_H_pi_interactions = obtain_H_pi_interaction_vectors(molecule, molecule_graph, index_to_attach_Hs_to, total_no_of_neighbours, no_of_pairs_of_lone_electrons, crystal, crystal_graph, element_to_attach=element_to_attach)
	
	# 4.2: Print information to the logger.
	if (len(new_bonding_unit_vectors_from_H_pi_interactions) > 0) and (logger is not None):
		logger.info('Pi: '+str(new_bonding_unit_vectors_from_H_pi_interactions))

	# 4.3: Yield all the unit vectors for adding a hydrogen to a molecule near a pi interaction.
	for new_bonding_unit_vector in new_bonding_unit_vectors_from_H_pi_interactions:
		yield new_bonding_unit_vector

	# ------------------------------------------------------------------------------------------------------------------------------
	# Fifth, obtain the positions of hydrogens such that they are onset of atoms from the neighbours.

	# 5.1: Obtain the bond vectors for atoms that are attached in a methyl-like fashion. 
	#      * This will give the vector required for a staggered conformation. 
	new_bonding_unit_vectors_from_between_neighbour_atoms_neighbours = obtain_ideal_methyl_like_hydrogen_positions(molecule, molecule_graph, index_to_attach_Hs_to, total_no_of_neighbours, no_of_pairs_of_lone_electrons, crystal, crystal_graph)

	# 5.2: Print information to the logger.
	if (len(new_bonding_unit_vectors_from_between_neighbour_atoms_neighbours) > 0) and (logger is not None):
		logger.info('Methyl-like bond unit vectors: '+str(new_bonding_unit_vectors_from_between_neighbour_atoms_neighbours))	

	# 5.3: Yield all the unit vectors for adding a hydrogen to an atom in a molecule in a methyl-like array.
	for new_bonding_unit_vector in new_bonding_unit_vectors_from_between_neighbour_atoms_neighbours:
		yield new_bonding_unit_vector

	# ------------------------------------------------------------------------------------------------------------------------------
	# Sixth, if default vectors are given, return them

	# 6.1: Obtain the default vectors
	if default is None:	
		defaults = []
	elif isinstance(default[0],(int,float)):
		if logger is not None:
			logger.info('Default vectors being used')
		defaults = [default]
	else:
		if logger is not None:
			logger.info('Default vectors being used')
		defaults = default

	# 6.2: Yield the default vectors.
	for vector in defaults:
		yield get_unit_vector(np.array(vector))

	# ------------------------------------------------------------------------------------------------------------------------------
	# Seventh, yield randomly generated unit vectors.

	# 7.2: Print information to the logger.
	if logger is not None:
		logger.info('Yielding randomly generated unit vectors')

	# 7.2: If at this point, return a randomly generated unit vector.
	for _ in range(100):
		yield get_unit_vector(np.array((uniform(-1.0,1.0), uniform(-1.0,1.0), uniform(-1.0,1.0))))
	else:
		raise Exception('Huh? There has been repetition in a neverending for loop.')

# ----------------------------------------------------------------------------------------------------------------------------------




