"""
get_bonding_unit_vectors.py, Geoffrey Weal, 6/7/2022

This script is designed to determine the current bonding unit vectors about the atom of interest, which is index_to_attach_Hs_to
"""
from SUMELF.SUMELF.general_methods.geometry_methods import get_unit_vector

def get_bonding_unit_vectors(molecule, molecule_graph, index_to_attach_Hs_to):
	"""
	This method will obtain the unit vectors for all neighbours attached to index_to_attach_Hs_to.

	Parameters
	----------
	molecule : ase.Atoms
		This is the molecule you want to remove the aliphatic carbons to.
	molecule_graph : networkx.Graph
		This is the graph of this molecule.
	index_to_attach_Hs_to :list 
		This is the atom to attach new hydrogens to.

	Returns
	-------
	bonding_unit_vectors : list of numpy.arrays
		These are the unit vectors that point in the directions of existing bonds to neighbours around the central atom.
	"""

	# First, determine the position of the atom to bind new hydrogens to.
	atom_position = molecule[index_to_attach_Hs_to].position

	# Second, initialise the list to add bonding unit vectors to. These vector are from atom index_to_attach_Hs_to to each of its neighbours. 
	bonding_unit_vectors = []

	# Third, obtain the binding unit vectors for neighbouring atoms attached to this atoms.
	for neighbour_index in molecule_graph.adj[index_to_attach_Hs_to].keys():

		# 3.1: Obtain the position of the neighbour atom
		neighbour_position = molecule[neighbour_index].position

		# 3.2: Obtain the bonding vector from the current atom and its neighbour
		bonding_vector = neighbour_position - atom_position

		# 3.3: Convert bonding_vector into a unit vector
		bonding_unit_vector = get_unit_vector(bonding_vector)

		# 3.4: Append bonding_unit_vector to bonding_unit_vectors.
		bonding_unit_vectors.append((neighbour_index, bonding_unit_vector))

	# Fourth, return the bonding_unit_vectors
	return bonding_unit_vectors