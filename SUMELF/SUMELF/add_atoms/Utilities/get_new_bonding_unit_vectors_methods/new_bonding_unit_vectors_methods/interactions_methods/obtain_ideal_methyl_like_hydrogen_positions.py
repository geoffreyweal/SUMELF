"""
obtain_ideal_methyl_like_hydrogen_positions.py, Geoffrey Weal, 27/4/24

This method is designed to provide the vectors for creating a methyl group with H atoms in a staggered conformation.
"""
import numpy as np
from ase.data import atomic_numbers

from SUMELF.SUMELF.general_methods.geometry_methods                      import get_unit_vector
from SUMELF.SUMELF.general_methods.geometry_methods                      import rotate_vector_around_axis

def obtain_ideal_methyl_like_hydrogen_positions(molecule, molecule_graph, index_to_attach_Hs_to, total_no_of_neighbours, no_of_pairs_of_lone_electrons, crystal, crystal_graph):
	"""
	This method is designed to provide the vectors for creating a methyl group with H atoms in a staggered conformation.

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

	# First, obtain the number of neighbours about the atom of interest.
	neighbouring_indices = tuple(molecule_graph[index_to_attach_Hs_to].keys())

	# Second, determine if this atom contains 4 neighbours.
	if total_no_of_neighbours + no_of_pairs_of_lone_electrons < 4:
		return []

	# Third, determine if the atom is currently only bonded to one atom. 
	if len(neighbouring_indices) != 1:
		return []

	# Fourth, obtain the index of the only neighbour bond to index_to_attach_Hs_to.
	neighbour_index = neighbouring_indices[0]

	# Fifth, obtain the neighbours attached to neighbour_index.
	neighbours_of_neighbour_indices = tuple([atom_index for atom_index in molecule_graph[neighbouring_indices[0]].keys() if (atom_index != index_to_attach_Hs_to)])

	# Sixth, determine the elements of neighbours_of_neighbour_indices.
	elements_of_nns = [molecule[atom_index].symbol for atom_index in neighbours_of_neighbour_indices]

	# Seventh, obtain the atomic number of neighbours_of_neighbour_indices. 
	atomic_numbers_of_nns = [atomic_numbers[symbol] for symbol in elements_of_nns]

	# Eighth, sort the neighbours by their atomic number
	sorted_neighbours_of_neighbour_indices = [atom_index for atom_index, atomic_number in sorted(zip(neighbours_of_neighbour_indices, atomic_numbers_of_nns), key=lambda x:(-x[1], x[0]))]

	# Ninth, initialise the list to hold the new bonding vectors.
	new_bond_vectors = []

	# Tenth, for each neighbour in sorted_neighbours_of_neighbour_indices. 
	for nn_index in sorted_neighbours_of_neighbour_indices:

		# 10.1: Obtain the bonding vector from nn_index to neighbour_index
		nn_bond_unit_vector = get_unit_vector(molecule[neighbour_index].position - molecule[nn_index].position)

		# 10.2: Obtain the bonding vector from nn_index to index_to_attach_Hs_to
		from_neighbour_to_focus_atom_bond_unit_vector = get_unit_vector(molecule[index_to_attach_Hs_to].position - molecule[neighbour_index].position)

		# 10.3: Obtain the rotation vector for rotating vectors from from_neighbour_to_focus_atom_bond_unit_vector upwards
		rotation_vector = -get_unit_vector(np.cross(nn_bond_unit_vector, from_neighbour_to_focus_atom_bond_unit_vector))

		# 10.4: Obtain the bonding vector for adding a new hydrogen. 
		new_bond_vector = get_unit_vector(rotate_vector_around_axis(from_neighbour_to_focus_atom_bond_unit_vector, np.radians(180-109.5), rotation_vector))

		# 10.5: Add new_bond_vector to new_bond_vectors
		new_bond_vectors.append(new_bond_vector)





	'''
	from ase import Atom

	debugging_molecule = molecule.copy()

	central_atom_position = molecule[index_to_attach_Hs_to].position

	for new_bonding_unit_vector in new_bond_vectors:

		added_position = central_atom_position + new_bonding_unit_vector

		debugging_molecule.append(Atom('X', added_position))

	from ase.visualize import view

	view([molecule, debugging_molecule])

	import pdb; pdb.set_trace()

	'''







	# Eleventh, return the bonding vectors
	return new_bond_vectors

	