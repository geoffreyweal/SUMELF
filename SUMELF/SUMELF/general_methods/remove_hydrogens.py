"""
remove_hydrogens.py, Geoffrey Weal, 19/7/22

This script is designed to hold methods useful for removing hydrogens from molecules.
"""
import numpy as np
import networkx as nx
from copy import deepcopy

def remove_hydrogens(molecule, graph=None, remove_deuterium=True, give_no_H_atoms_on_each_atom=False):
	"""
	This method will remove all hydrogens from your molecule.

	Parameters
	----------
	molecule : ase.Atoms
		This is the molecule you to remove all hydrogens from.
	graph : nx.Graph
		This is the graph that is associated with the molecule. Give as None if you do not want to give a graph for this molecule. 
	remove_deuterium : bool.
		If true, deuterium atoms will also be removed. 
	give_no_H_atoms_on_each_atom : bool.
		This boolean will indicate if the user wants to return a list that presents all the hydrogens attached to each atom in the molecule.

	Returns
	-------
	copied_molecule : ase.Atoms
		This is the molecule with all hydrogens removed.
	copied_graph : nx.Graph
		This is the graph of the molecule with all hydrogens removed.
	"""

	# First, take note of indices that have been removed, and make a copy of the indices and molecule object.
	remove_indices = []
	original_indices = list(range(len(molecule)))
	copied_molecule = molecule.copy()
	if give_no_H_atoms_on_each_atom:
		no_of_H_on_atoms_in_molecule = []

	# Second, take the AddedOrModifiedAtomList list from the copied_molecule object.
	if 'AddedOrModifiedAtomList' in copied_molecule.arrays:
		raise Exception('Check, as AddedOrModifiedAtomList may have been removed from protocols.')

	# Third, remove any atoms that are hydrogens.
	for atom_index in range(len(molecule)-1,-1,-1):
		symbol =  molecule[atom_index].symbol
		if symbol == 'H' or (remove_deuterium and symbol == 'D'):
			del copied_molecule[atom_index]
			del original_indices[atom_index]
			remove_indices.append(atom_index)

	# Fourth, determine what should be returned.
	return_graph = not (graph is None)

	# Fifth, copy and update the graph with no hydrogens if desired:
	if return_graph:
		copied_graph = deepcopy(graph)
		add_hydrogen_no_to_graph(copied_graph, remove_indices)
		copied_graph.remove_nodes_from(remove_indices)
		original_to_new_indices = {original_indices[atom_index]: atom_index for atom_index in range(len(original_indices))}
		copied_graph = nx.relabel_nodes(copied_graph, original_to_new_indices)

	# Sixth, obtain the number of hydrogens attached to each atom in the molecule.
	if give_no_H_atoms_on_each_atom:
		no_of_H_on_atoms_in_molecule = [copied_graph.nodes[no]['H'] for no in copied_graph.nodes]

	# Seventh, return molecules and graphs of molecules without hydrogens attached to them, and no_of_H variables based on inputs.
	if return_graph and give_no_H_atoms_on_each_atom:
		return copied_molecule, copied_graph, no_of_H_on_atoms_in_molecule
	elif return_graph and (not give_no_H_atoms_on_each_atom):
		return copied_molecule, copied_graph
	elif (not return_graph) and give_no_H_atoms_on_each_atom:
		return copied_molecule, no_of_H_on_atoms_in_molecule
	else:
		return copied_molecule

# ------------------------------------------------------------------------------------------------------------------------

def add_hydrogen_no_to_graph(copied_graph, remove_indices):
	"""
	This method is designed to label atoms in the graph with the number of hydrogen atoms they are attached to.
	"""

	# First, set up the dict that records which atoms are bound to hydrogens in the molecule.
	no_of_hydrogens_attached_to_atoms_in_molecule = {}

	# Second, set a boolean to record if no_of_neighbouring_non_cord_H is given in the molecule graphs. 
	all_does_contain_no_of_neighbouring_non_cord_H_details     = []
	all_does_not_contain_no_of_neighbouring_non_cord_H_details = []

	# Third, record hydrogens attached to each atom in the molecule. 
	for node_index, node_features in copied_graph.nodes.items():

		# 3.1: If the node is a hydrogen atom, ignore it.
		if node_index in remove_indices:
			continue

		# 3.2: Get the indices of neighbouring atoms bound to the node.
		neighbour_indices = list(copied_graph[node_index].keys())

		# 3.3: Determine if any of the neighbouring atoms are hydrogens.
		hydrogens_attached_to_atom = tuple(set(neighbour_indices) & set(remove_indices))

		# 3.4: Determine the number of hydrogen atoms bound to the node.
		no_of_hydrogens_attached_to_atom = len(hydrogens_attached_to_atom)

		# 3.5: Record if 'no_of_neighbouring_non_cord_H' was in the node_features for this node in copied_graph
		if 'no_of_neighbouring_non_cord_H' in node_features:
			all_does_contain_no_of_neighbouring_non_cord_H_details.append((copied_graph.name, node_index))
		else:
			all_does_not_contain_no_of_neighbouring_non_cord_H_details.append((copied_graph.name, node_index))

		# 3.6: If 'no_of_neighbouring_non_cord_H' is given, add it to no_of_hydrogens_attached_to_atom
		if 'no_of_neighbouring_non_cord_H' in node_features:
			no_of_hydrogens_attached_to_atom += node_features['no_of_neighbouring_non_cord_H']
			del node_features['no_of_neighbouring_non_cord_H']

		# 3.7: Record the number of hydrogen atoms bound to this node.
		no_of_hydrogens_attached_to_atoms_in_molecule[node_index] = no_of_hydrogens_attached_to_atom

	# Fourth, make sure that there is consistancy across all atoms for all molecules for the 'no_of_neighbouring_non_cord_H' feature.
	#         * We want either all atom to have the 'no_of_neighbouring_non_cord_H' feature, or none of them have this feature. 
	if (len(all_does_contain_no_of_neighbouring_non_cord_H_details) != 0) and (len(all_does_not_contain_no_of_neighbouring_non_cord_H_details) != 0):
		toString  = 'Error: There is inconsistancy for the "no_of_neighbouring_non_cord_H" feature. Some atoms have them, some do not.'+'\n'
		tostring += 'Nodes containing "no_of_neighbouring_non_cord_H":     '+str(all_does_contain_no_of_neighbouring_non_cord_H_details)+'\n'
		tostring += 'Nodes not containing "no_of_neighbouring_non_cord_H": '+str(all_does_not_contain_no_of_neighbouring_non_cord_H_details)+'\n'
		tostring += 'Check this, as this indicates a programming error.'
		raise Exception(tostring)

	# Fifth, record the number of atoms bound to nodes in this molecule.
	nx.set_node_attributes(copied_graph, no_of_hydrogens_attached_to_atoms_in_molecule, 'H')
	
# ------------------------------------------------------------------------------------------------------------------------


