"""
is_solvent.py, Geoffrey Weal, 14/2/24

This method is designed to determine if a molecule is a solvent or not.
"""
import networkx as nx
from copy import deepcopy
import networkx.algorithms.isomorphism as iso
from SUMELF.SUMELF.general_methods.is_solvent_methods.load_graphs_of_known_solvents import load_graphs_of_known_solvents

nm = iso.categorical_node_match(attr=['E', 'H'], default=[None, 0])
def is_solvent(non_hydrogen_graph):
	"""
	This method is designed to determine if a molecule is a solvent or not.

	Parameters
	----------
	non_hydrogen_graph : networkx.Graph
		This is the graph of the molecule, with hydrogens attached to nodes of the "heavy" atoms in the graph. If they are not, the method will first remove hydrogen nodes and attach them to the "heavy" atom as node variables. 

	Returns
	-------
	True if the molecule is a solvent. False if not. 
	"""

	# First, if the solvent_graphs global variable does not currently exist, load it into memory.
	if "solvent_graphs" not in globals():
		global solvent_graphs
		solvent_graphs = load_graphs_of_known_solvents()

	# Second, make a copy of the graphs
	non_hydrogen_graph = deepcopy(non_hydrogen_graph)

	# Third, convert the graph into one that does not contain hydrogen nodes, but hydrogens are attached to the "heavy" atoms nodes of graphs as node variables. 
	#        * This method will also add hydrogens from "no_of_neighbouring_non_cord_H" to the hydrogens list. 
	non_hydrogen_graph = remove_hydrogen_nodes(non_hydrogen_graph)

	# Fourth, for each solvent_graph in the solvent_graphs list.
	for solvent_graph in solvent_graphs:

		# 4.1: If non_hydrogen_graph is the same as the solvent_graph, non_hydrogen_graph is a solvent. 
		if nx.is_isomorphic(non_hydrogen_graph, solvent_graph, node_match=nm):

			# 4.3: non_hydrogen_graph is the same as solvent_graph, return True
 			return True

	# Fifth, non_hydrogen_graph does not match any solvent in the solvent_graphs list, return False.
	return False

def remove_hydrogen_nodes(original_graph):
	"""
	This method is designed to remove the hydrogens as nodes from the original graph and adding them to the "heavy" atom nodes in the graph as node variables.

	This method will also add atoms from 'no_of_neighbouring_non_cord_H' to the "heavy" atoms number of hydrogen variable. 

	Parameters
	----------
	original_graph : networkx.Graph
		This is the original graph for the molecule. 

	Returns
	-------
	non_hydrogen_graph : networkx.Graph
		This is the graph to remove hydrogens as nodes from, and add them and hydrogens from 'no_of_neighbouring_non_cord_H' to the "heavy" atom nodes in graph as node variables.  
	"""

	# First, initalise the atoms and hydrogen_indices lists
	atoms = []
	hydrogen_indices = []

	# Second, for each atom (node) in original_graph
	for atom_index, atom_details in original_graph.nodes.items():

		# 2.1: Determine what the key is in atom_details that indicates what the element is. 
		if   'element' in atom_details.keys():
			raise Exception('Error: element is old, should be E')
		elif 'E'       in atom_details.keys():
			element_key = 'E'
		else:
			raise Exception('Error: no element key')

		# 2.2: Get the element of the atom
		element = atom_details[element_key]

		# 2.3: If the atom is a hydrogen, record the indices and move on to the next atom in original_graph
		if element in ['H', 'D', 'T']:
			hydrogen_indices.append(atom_index)
			continue

		# 2.4: Initialise the int to record the number of hydrogens bound to the atom.
		number_of_bonded_hydrogens = 0

		# 2.5: For each atom bound to atom_index.
		for atom_bonded_to in original_graph[atom_index]:

			# 2.5.1: Obtain the details of the bonded atom
			neighbouring_atom_details = original_graph.nodes[atom_bonded_to]

			# 2.5.2: Get the element of the bonded atom
			neighbouring_atom_element = neighbouring_atom_details[element_key]

			# 2.5.3: If the bonded atom is a hydrogen, increment number_of_bonded_hydrogens
			if neighbouring_atom_element in ['H', 'D', 'T']:
				number_of_bonded_hydrogens += 1

		# 2.6: If 'no_of_neighbouring_non_cord_H' is in atom_details, add it to number_of_bonded_hydrogens
		if 'no_of_neighbouring_non_cord_H' in atom_details:
			number_of_bonded_hydrogens += atom_details['no_of_neighbouring_non_cord_H']

		# 2.7: Add the element to the atoms list.
		atoms.append((atom_index, {'E': element, 'H': number_of_bonded_hydrogens}))

	# Third, initialise the bonds list.
	bonds = []

	# Fourth, for each bond (edge) in original_graph.
	for atom1_index, atom2_index in original_graph.edges.keys():

		# 4.1: If either atom1_index or atom2_index is in hydrogen_indices, it is a hydrogen, so dont include it.
		if (atom1_index in hydrogen_indices) or (atom2_index in hydrogen_indices):
			continue

		# 4.2: Create the bond tuple
		bond = tuple(sorted([atom1_index, atom2_index]))

		# 4.3: Add bond to the bonds list. 
		bonds.append(bond)

	# Fifth, initialise the networkx.Graph object
	non_hydrogen_graph = nx.Graph()

	# Sixth, add the atoms to the graph as nodes
	non_hydrogen_graph.add_nodes_from(atoms)

	# Seventh, add the bonds to the graph as edges
	non_hydrogen_graph.add_edges_from(bonds)

	# Eighth, remap the solvent graph so that the indices are consecutive. 
	remap_atoms = {old_atom_index: new_atom_index for (new_atom_index, old_atom_index) in enumerate(sorted(non_hydrogen_graph.nodes.keys()))}

	# Ninth, remap the nodes in solvent_graph with consecutive indices. 
	non_hydrogen_graph = nx.relabel_nodes(non_hydrogen_graph, remap_atoms)

	# Tenth, return non_hydrogen_graph
	return non_hydrogen_graph

def contains_hydrogen(non_hydrogen_graph):
	"""
	This method is designed to determine if the graph has hydrogens as nodes.

	Parameters
	----------
	non_hydrogen_graph : networkx.Graph
		This is the graph to check if it has hydrogens as nodes.

	Returns
	-------
	True if the graph contains hydrogens as nodes. False if not. 
	"""

	# First, for each node in the non_hydrogen_graph.nodes nodes list. 
	for atom_index, atom_details in non_hydrogen_graph.nodes.items():

		# Second, determine what the key is in atom_details that indicates what the element is. 
		if   'element' in atom_details.keys():
			raise Exception('Error: element is old, should be E')
		elif 'E'       in atom_details.keys():
			element_key = 'E'
		else:
			raise Exception('Error: no element key')

		# Third, if the atom is a hydrogen, return True
		if atom_details[element_key] in ['H', 'D', 'T']:
			return True

	# Fourth, did not find a hydrogen node, so return false. 
	return False

# ========================================================================================================



