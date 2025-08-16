"""
Geoffrey Weal, add_graph_to_ASE_Atoms_object.py, 25/1/2024

This script is designed to add information about the atoms and bonds for a molecule from a networkx object into the ASE object of the molecule itself. 
"""
import numpy as np
from copy import deepcopy
from itertools import permutations, product

def add_graph_to_ASE_Atoms_object(molecule, molecule_graph):
	"""
	This script is designed to add information about the atoms and bonds for a molecule from a networkx object into the ASE object of the molecule itself. 

	Parameters
	----------
	molecule : ase.Atoms
		This is the molecule of interest to add graph details to
	molecule_graph : networkx.graph
		This is the graph that resembles the molecule given by the ASE object. 
	"""

	# ============================================================================================================
	# First, copy graph node details to the atoms in the ASE Atoms object for the molecule

	# 1.1: Get all the atom property names for nodes in molecule_graph
	all_atom_property_names = set()
	for node_properties in molecule_graph.nodes.values():
		all_atom_property_names.update(node_properties.keys())

	# 1.2: Check that all atom properties are found in all atoms (nodes) in molecule_graph
	nodes_that_have_all_properties = []
	node_property_name_errors = {}
	for node_index, node_properties in molecule_graph.nodes.items(): 
		node_property_names = set(node_properties.keys())
		difference = all_atom_property_names - node_property_names
		if len(difference) > 0:
			node_property_name_errors[node_index] = list(difference)
		else:
			nodes_that_have_all_properties.append(node_index)
	if len(node_property_name_errors) > 0:
		to_string  = 'Error: Some atoms (nodes) in the graph do not contain expected property names.\n'
		to_string += f'Nodes with missing properties: {sorted(node_property_name_errors.items())}\n'
		to_string += f'Nodes with no missing properties: {sorted(nodes_that_have_all_properties)}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 1.3: Remove 'E' from node properties (if it is given). This is the element of the atom, which is already given in the ase.Atoms object.
	for prefix, suffix in product(['E', 'element', 'number', 'position', 'charge', 'initial_charge'], ['', 's']):
		atom_property = prefix+suffix
		if atom_property in all_atom_property_names:
			all_atom_property_names.remove(atom_property)

	# 1.4: Initialse the properties to record
	atom_properties = {property_name: [] for property_name in all_atom_property_names}

	# 1.5: For each node in the graph.
	for index, (node_index, node_properties) in enumerate(molecule_graph.nodes.items()):

		# 1.5.1: Check that the index ordering is correct for this list.
		#        * The indices in molecule_graph should be consecutive from 0 to len(molecule_graph)
		if not index == node_index:
			raise Exception('Huh?')

		# 1.5.2: Copy the information for the node from the graph into a list in atom_properties
		for atom_property in atom_properties:
			atom_properties[atom_property].append(node_properties[atom_property])

	# 1.6: Add each property in atom_properties into the ASE Atoms object for molecule.
	for atom_property in atom_properties:
		molecule.set_array(atom_property, np.array(atom_properties[atom_property]))

	# ============================================================================================================
	# Second, move NeighboursList to the end of the array list in the ASE atoms object.
	#    * This is show NeighboursList is given at the end of each atom in the saved xyz file

	# 2.1: Initialise the NeighboursList
	NeighboursList = []

	# 2.2: Collect all the neighbours for all atoms in the molecule
	for atom_index in molecule_graph.nodes.keys():

		# 2.2.1: Initialse the bonded_indices that will store the indices of all the atoms bonded to atom_index
		bonded_indices = []

		# 2.2.2: For each bond involving atom_index
		for (bond_atom1_index, bond_atom2_index) in molecule_graph.edges(atom_index):

			# 2.2.3: Check that the indices in the bond are not the same (I.e. the atom is not bonded to itself).
			if bond_atom1_index == bond_atom2_index:
				raise Exception('Error: Atom bonded to itself: '+str((bond_atom1_index, bond_atom2_index)))

			# 2.2.4: Obtain the other atom in the bond
			if   atom_index == bond_atom1_index:
				other_atom_index = bond_atom2_index
			elif atom_index == bond_atom2_index:
				other_atom_index = bond_atom1_index
			else:
				raise Exception('Error: atom_index not in bond. atom_index: '+str(atom_index)+'; bond: '+str((bond_atom1_index, bond_atom2_index)))

			# 2.2.5: Add the bonded atom as a neighbour of atom_index to bonded_indices
			bonded_indices.append(other_atom_index)
				
		# 2.2.6: Sort bonded_indices
		bonded_indices.sort()

		# 2.2.7: Convert the list into a string
		if len(bonded_indices) == 0:
			bonded_indices = '-'
		else:
			bonded_indices = ','.join([str(index) for index in bonded_indices])

		# 2.2.8: Add bonded_indices to NeighboursList
		NeighboursList.append(bonded_indices)

	# 2.3: Convert list to numpy array
	NeighboursList = np.array(NeighboursList)

	if 'hybridisation' in molecule.arrays.keys():
		hybridisation = molecule.get_array('hybridisation')
		del molecule.arrays['hybridisation']
		molecule.set_array('hybridisation', hybridisation)

	# 2.4: Add NeighboursList array to molecule
	molecule.set_array('NeighboursList',NeighboursList)

	# ============================================================================================================
	# Third, copy graph edge details to the ASE Atoms object as bonds for the molecule

	# 3.1: Get all the bond property names for edges in molecule_graph
	all_bond_property_names = set()
	for edge_properties in molecule_graph.edges.values():
		all_bond_property_names.update(edge_properties.keys())

	# 3.2: Check that all bond properties are found in all bonds (edges) in molecule_graph
	edge_property_name_errors = {}
	for edge_index, edge_properties in molecule_graph.edges.items(): 
		edge_property_names = set(edge_properties.keys())
		difference = all_bond_property_names - edge_property_names
		if len(difference) > 0:
			edge_property_name_errors[node_index] = list(difference)
	if len(edge_property_name_errors) > 0:
		raise Exception('Error: Some atoms (nodes) in the graph do not contain expected property names: '+str(edge_property_name_errors))

	# 3.1: Initialse the properties to record
	bond_properties = {property_name: [] for property_name in all_bond_property_names}

	# 3.2: For each edge in the graph.
	for index, (edge_indices, edge_properties) in enumerate(sorted(molecule_graph.edges.items())):

		# 3.2.1: Check that this edge exists as a bond in the molecule object by inspecting the NeighboursList
		for index1, index2 in permutations(edge_indices):
			neighbours = eval('['+NeighboursList[index1]+']')
			if not index2 in neighbours:
				raise Exception('Huh?')

		# 3.2.2: Copy the information for the edge from the graph into a list in bond_properties
		for bond_property in bond_properties:
			bond_properties[bond_property].append((edge_indices, edge_properties[bond_property]))

	# 3.3: Add each property in bond_properties into the ASE Atoms object for molecule (if there are properties in the bond_properties dictionary) 
	if len(bond_properties) > 0:

		# 3.3.1: Initialise a list that contains the bond properties that have issues. 
		bond_properties_with_issues = []

		# 3.3.2: Create a sorted list of all the bonds in the molecule from the molecule_graph
		molecule_edges = sorted(molecule_graph.edges())

		# 3.3.3: For each bond property in bond_properties
		for bond_property_key, bond_property_values in bond_properties.items():

			# 3.3.3.1: Obtain the bonds (edges) for all the bonds in bond_properties[bond_property_key]
			bond_property_edges = sorted([bond_atom_indices for bond_atom_indices, value in bond_property_values])

			# 3.3.3.2: If bond_property_edges does not match molecule_edges, there is a problem. 
			if not bond_property_edges == molecule_edges:

				# 3.3.3.3: Check the bonds that differ between bond_property_edges and molecule_edges
				bonds_not_found_in_molecule_edges = bond_property_edges ^ molecule_edges

				# 3.3.3.4: Add problems to bond_properties_with_issues
				bond_properties_with_issues.append((bond_property_key, bonds_not_found_in_molecule_edges))

		# 3.3.4: Report any inconsistencies with bond_properties
		if len(bond_properties_with_issues) > 0:
			to_string  = 'Error: Some bonds given in BondProperties do not match the bonds in the NeighboursList.\n'
			to_string += 'Bond properties in BondProperties that are inconsistent with NeighboursList (protperty name: inconsistent bonds):\n'
			for bond_property_key, bonds_not_found_in_molecule_edges in bond_properties_with_issues:
				to_string += str(bond_property_key)+': '+str(list(bonds_not_found_in_molecule_edges))+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

		# 3.3.5: Add bond_properties to molecule.info
		molecule.info['BondProperties'] = str(bond_properties)

	# ============================================================================================================

