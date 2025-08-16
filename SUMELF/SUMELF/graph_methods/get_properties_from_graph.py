"""
get_property_names_from_graph.py, Geoffrey Weal, 26/1/24

This method is designed to obtain all the names of node/edge properties from a networkx Graph object.
"""
from copy import deepcopy

def get_node_property_names_from_graph(graph):
	"""
	This method is designed to obtain all the names of node properties from a networkx Graph object.

	Parameters
	----------
	graph : networkx.graph or None
		This is the graph to get node property names from.

	Returns
	-------
	node_property_names : list
		These are the name of all the properties in the graph object (except for E (element))
	"""

	# First, check that all the node in the graph are consistent (contain all the same properties)
	#        This method also returns the name of the properties in the node.
	node_property_names = check_all_nodes_have_same_property_names(graph)

	# Second, return node_property_names
	return node_property_names

# ================================================================================================================

def get_node_properties_from_graph(graph, property_name, remove_node_property_from_graph=False):
	"""
	This method will return a list of all the values of a property for all nodes in the graph

	Parameters
	----------
	graph : networkx.graph or None
		This is the graph to get node property names from.
	property_name : str.
		This is the name of the properties you would like to obtain.

	Returns
	-------
	node_property : dict
		These are the name of all the properties in the graph object (except for E (element))
	"""

	# First, check that each node in graph contains property_name 
	check_all_nodes_have_property(graph, property_name)

	# Second, initialise dictionary for recording property_name value for each node in the graph
	property_values = {}

	# Third, for each node in the graph
	for node_index, node_properties in graph.nodes.items():

		# 3.1: property_name should be in node_properties. 
		#      If not, their is a programming issue as this should 
		#      have been spotted in the check_all_nodes_have_property method.
		if not (property_name in node_properties.keys()):
			raise Exception('Programming Error: Should not get here, as this error should have been spotted in check_all_nodes_have_property. Could not find "'+str(property_name)+'" in node_properties.')

		# 3.2: Get the value of the property from node_properties.
		property_values[node_index] = deepcopy(node_properties[property_name])

	# Fourth, remove property_name from the nodes in the graph if desired. 
	if remove_node_property_from_graph:
		remove_node_properties_from_graph(graph, property_name)

	# Fifth, return property_values
	return property_values

# ================================================================================================================

def add_node_properties_to_graph(graph, property_name, property_value_for_nodes):
	"""
	This method will remove a property from the nodes in the graph.

	Parameters
	----------
	graph : networkx.graph or None
		This is the graph to get node property names from.
	property_name : str.
		This is the name of the properties you would like to obtain.
	property_value_for_nodes : list
		This is the list of properties for property_name to add to each node in graph
	"""

	# First, check that each node in graph does not already contains property_name 
	#        We do not want this.
	check_all_nodes_do_not_have_property(graph, property_name)

	# Second, remove property_name from each node in the graph
	for node_index, node_properties in graph.nodes.items():

		# 2.1: property_name should not be in node_properties. 
		#      If so, their is a programming issue as this should 
		#      have been spotted in the check_all_nodes_do_not_have_property method.
		if (property_name in node_properties.keys()):
			raise Exception('Programming Error: Should not get here, as this error should have been spotted in check_all_nodes_do_not_have_property. Found "'+str(property_name)+'" in node_properties when it should not be in node_properties.')

		# 2.2: Remove property_name from each node in the graph.
		node_properties[property_name] = deepcopy(property_value_for_nodes[node_index])

# ================================================================================================================

def remove_node_properties_from_graph(graph, property_name, force_remove=False):
	"""
	This method will remove a property from the nodes in the graph.

	It is initally required that all atoms in the graph contain the property, as this should be the case before removing it.

	If you dont want this requirement, set force_remove to True

	Parameters
	----------
	graph : networkx.graph or None
		This is the graph to get node property names from.
	property_name : str.
		This is the name of the properties you would like to obtain.
	force_remove : bool.
		This indicates if you want to remove the property from the node, even if the node does not contain the the property.

	Returns
	-------
	node_property : dict
		These are the name of all the properties in the graph object (except for E (element))
	"""

	# First, check that each node in graph contains property_name 
	if not force_remove:
		check_all_nodes_have_property(graph, property_name)

	# Second, remove property_name from each node in the graph
	for node_index, node_properties in graph.nodes.items():

		# 3.1: If the node does not contain the property, either:
		#      * If force_remove is True  -> continue 
		#      * If force_remove is False -> Report a programming error. 
		if not (property_name in node_properties.keys()):
			if force_remove:
				continue
			else:
				to_string  = 'Programming Error: Should not get here, as this error should have been spotted in check_all_nodes_have_property.\n'
				to_string += f'Could not find "{property_name}" in node_properties.\n'
				to_string += 'Check this.'
				raise Exception(to_string)

		# 3.2: Remove property_name from each node in the graph.
		del node_properties[property_name]

# ================================================================================================================
# ================================================================================================================
# ================================================================================================================

def check_all_nodes_have_same_property_names(graph):
	"""
	This method is designed to check that every node in the graph contains the same properties

	Parameters
	----------
	graph : networkx.graph or None
		This is the graph to get node property names from.

	Returns
	-------
	node_property_names : list
		These are the name of all the properties in the graph object (except for E (element))
	"""

	# First, initialse the node_property_names as a set
	node_property_names = set()

	# Second, obtain all the names of properties from all nodes in the graph
	for node_properties in graph.nodes.values():
		node_property_names.update(node_properties.keys())

	# Third, remove 'E', as this is not an important property to return (it is the element of the atom, that is kept for tracking reasons. The ASE Atoms object records all of this). 
	if 'E' in node_property_names:
		node_property_names.remove('E')

	# Fourth, check that each property in node_property_names is in each node in graph
	errors = {}
	for node_index, node_properties in graph.nodes.items(): 
		node_properties = set(node_properties.keys())
		differences = node_property_names - node_properties
		if len(differences) > 0:
			errors[node_index] = sorted(differences)

	# Fifth, remove an error if some of the atoms (nodes) in the graph do not contain a property.
	if len(errors) > 0:
		raise Exception('Error. Not all atoms have the same properties. The following atoms are missing property information: '+str(errors))

	# Sixth, return node_property_names as a sorted list:
	return sorted(node_property_names)

def check_all_nodes_have_property(graph, property_name):
	"""
	This method will check that all nodes in the graph contain the property_name.

	Parameters
	----------
	graph : networkx.graph or None
		This is the graph to get node property names from.
	property_name : str.
		This is the name of the properties you would like to obtain.
	"""

	# First, initialise the list for nodes that do not have the property you want to obtain.
	could_not_find_property_for_nodes = []

	# Second, for each node in the graph
	for node_index, node_properties in graph.nodes.items():

		# Check if the node has this property name. If not, record this as an error.
		if not (property_name in node_properties.keys()):
			could_not_find_property_for_nodes.append(node_index)

	# Third, raise an issue if one/some atom(s) do not have this property
	if len(could_not_find_property_for_nodes):
		raise Exception('Error: The following nodes do not contain "'+str(property_name)+'" as a property in this graph: '+str(sorted(could_not_find_property_for_nodes)))

def check_all_nodes_do_not_have_property(graph, property_name):
	"""
	This method will check that all nodes in the graph do not contain the property_name.

	Parameters
	----------
	graph : networkx.graph or None
		This is the graph to get node property names from.
	property_name : str.
		This is the name of the properties you would like to obtain.
	"""

	# First, initialise the list for nodes that contain the property you want to obtain.
	found_property_for_nodes = []

	# Second, for each node in the graph
	for node_index, node_properties in graph.nodes.items():

		# 2.1: Check if the node has this property name. If it does, record this as an error.
		if (property_name in node_properties.keys()):
			found_property_for_nodes.append(node_index)

	# Third, raise an issue if one/some atom(s) have this property
	if len(found_property_for_nodes):
		raise Exception('Error: One or more nodes contain "'+str(property_name)+'" as a property in this graph.\nNone of the nodes should have this property.\nNodes that contain "'+str(property_name)+'" are: '+str(sorted(found_property_for_nodes)))


# ================================================================================================================
# ================================================================================================================
# ================================================================================================================


