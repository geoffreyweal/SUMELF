"""
general_methods.py, Geoffrey Weal, 17/2/22

This script includes methods for reconstructing a molecule from the crystal structure using a component assembly approach.
"""
from networkx import Graph, connected_components

from SUMELF.SUMELF.general_methods.distance_methods  import get_distance, less_than_or_equal_to_max_bondlength, get_covalent_bond_distance
from SUMELF.SUMELF.general_methods.unit_cell_methods import get_cell_corner_points

def create_molecule_using_component_assembly_approach(molecule_in_crystal, crystal_graph, crystal, mapping_atom_in_crystal_TO_atom_in_molecule_indices, take_shortest_distance=False, to_print=True, logger=None):
	"""
	This method will reconstruct a molecule from the crystal structure by creating adding 26 identical cells around the origin unit cell, and keeping only the component that reconstructs the molecule.

	Parameters
	----------
	molecule_in_crystal : ase.Atoms
		This is the molecule in the crystal that you want to reconstruct.
	crystal_graph : networkx.Graph
		This is the graph of the full crystal.
	crystal : ase.Atoms
		This is the original crystal.
	mapping_atom_in_crystal_TO_atom_in_molecule_indices : list
		This is the conversion of the indices of the atoms in the crystal --> to --> atoms in the molecule.
	take_shortest_distance : bool
		If true, take the shortest distance that was found if a distance shorter than the maximum distance between elements if found. If false, give an error if the distance between two atoms was greater than the max distance expected for those two elements to be bonded to each other. Default: False
	to_print : bool
		This tag will determine if you want to show a progress bar of atoms being added to the super_lattice, as this can take a bit fo time. Not used in this method. Default: True.
	logger : logging/None
		This contains the object for reporting information and warning messages to.

	Returns
	-------
	complete_molecule : ase.Atoms
		This is the molecule that has been attached back together again.
	molecule_graph : networkx.graph
		This is the graph associated to the completed molecule.
	"""

	# First, get information about the molecule and crystal	
	crystal_elements = crystal.get_chemical_symbols()
	crystal_cell_lattice = crystal.get_cell()

	# Second, get the points of displacement to move to lattices currounting the (0,0,0) cell
	cell_points = get_cell_corner_points(crystal_cell_lattice,super_cell_reach=1)

	# Third, determine which atoms are apart of which components of the molecule in the crystal, as well as the bonds that are disconnected in the crystal.
	# We will use the crystal_graph to make two graphs of the molecule in complete form and the molecule in component form. 
	disconnected_atom_pairs = []
	bonding_vectors = []
	molecule_component_graph_edges = []
	molecule_graph_edges = []
	for m_in_c_index1, molecule_index1 in mapping_atom_in_crystal_TO_atom_in_molecule_indices.items():
		element1 = crystal_elements[m_in_c_index1]
		for m_in_c_index2 in crystal_graph[m_in_c_index1]:
			molecule_index2 = mapping_atom_in_crystal_TO_atom_in_molecule_indices[m_in_c_index2]
			element2 = crystal_elements[m_in_c_index2]
			# Determine if a pair of bonded atoms fall off the unit cell from each other. 
			bondlength = crystal.get_distance(m_in_c_index1,m_in_c_index2,mic=False)
			if less_than_or_equal_to_max_bondlength(bondlength, element1, element2):
				molecule_component_graph_edges.append((molecule_index1, molecule_index2))
			else:
				disconnected_atom_pair = (molecule_index1, molecule_index2) if (molecule_index1 < molecule_index2) else (molecule_index2, molecule_index1)
				if not disconnected_atom_pair in disconnected_atom_pairs:
					disconnected_atom_pairs.append(disconnected_atom_pair)
					position1 = molecule_in_crystal[molecule_index1].position
					position2 = molecule_in_crystal[molecule_index2].position
					# Determine the bonding vector from atom 1 to atom 2
					displacement = get_molecule_cell_displacement(element1, element2, position1, position2, cell_points, m_in_c_index1, m_in_c_index2, molecule_in_crystal, take_shortest_distance=take_shortest_distance, logger=logger, crystal=crystal)
					bonding_vector = (position2 + displacement) - position1
					bonding_vectors.append(bonding_vector)
			# We know that these two atoms are connected, so add this edge to the molecule_graph in molecule indices
			molecule_graph_edges.append((molecule_index1, molecule_index2))

	# Fourth, create the graphs for obtaining molecules from the crystal.
	molecule_component_graph = Graph()
	molecule_component_graph.add_nodes_from(list((index, {'E': molecule_in_crystal[index].symbol}) for index in range(len(molecule_in_crystal))))
	molecule_component_graph.add_edges_from(molecule_component_graph_edges)
	del molecule_component_graph_edges
	molecule_graph = Graph() # make the molecule graph while we are doing things.
	molecule_graph.add_nodes_from(list((index, {'E': molecule_in_crystal[index].symbol}) for index in range(len(molecule_in_crystal))))
	molecule_graph.add_edges_from(molecule_graph_edges)
	del molecule_graph_edges

	# Fifth, add the bond vector from atom1 to atom2 in the molecule
	disconnected_atom_pairs = list(zip(disconnected_atom_pairs, bonding_vectors))

	# Sixth, determine the indices of all the components of the molecule in the crystal.
	molecule_components_graphs = list(connected_components(molecule_component_graph))

	# Seventh, reconnect the molecule by reconnecting the disconnected bonds.
	connected_molecule = connect_bonds_together(molecule_in_crystal, disconnected_atom_pairs, crystal_cell_lattice, molecule_components_graphs)

	# Eighth, return connected_molecule and molecule_graph
	return connected_molecule, molecule_graph

# ===========================================================================================================================================================================================

def get_molecule_cell_displacement(element1, element2, position1, position2, cell_points, m_in_c_index1, m_in_c_index2, molecule_in_crystal, take_shortest_distance=False, logger=None, crystal=None):
	"""
	This method will determine what displacement is needed to move atom 2 so that we can reconnect the bond from atom 1 to atom 2 in your crystal structure.

	Parameters
	----------
	element1: str.
		The element of the first atom.
	element2: str.
		The element of the second atom.
	position1 : numpy.array
		A 1x3 row matrix containing the x, y, and z coordinates of the atom1. Given in Å. 
	position2 : numpy.array
		A 1x3 row matrix containing the x, y, and z coordinates of the atom2. Given in Å. 
	cell_points : list
		This is a list of all the distances between identical points in each ijk section of the supercell from the origin unit cell. Given in Å. 
	molecule_in_crystal : ase.Atoms
		This is the molecule in the crystal. Used only for debugging.
	m_in_c_index1 : int
		This is one of the atoms in the crystal involved in the issue. Used only for debugging.
	m_in_c_index2 : int
		This is one of the atoms in the crystal involved in the issue. Used only for debugging.
	take_shortest_distance : bool
		If true, take the shortest distance that was found if a distance shorter than the maximum distance between elements if found. If false, give an error if the distance between two atoms was greater than the max distance expected for those two elements to be bonded to each other. Default: False
	crystal : ase.Atoms
		This is the original crystal.

	Returns
	-------
	displacement : numpy.array
		A 1x3 row matrix containing the displacement required to move atom 2 so that it is bonding distance to atom 1. Given in Å. 
	"""

	# First, initalise the shortest_displacement and shortest_bondlength variables.
	shortest_displacement = None
	shortest_bondlength = float('inf')

	# Second, for each displacement in cell_points
	for displacement in cell_points:

		# 2.1: Get the bond distance between atom1 and atom2
		bondlength = get_distance(position1, position2 + displacement)

		# 2.2: If the bondlength is less than the maximum distance between element1 and element2, 
		# then we have found a position where there is a bond between atom1 and atom2.
		# Therefore, return that displacement.
		if less_than_or_equal_to_max_bondlength(bondlength, element1, element2):
			return displacement

		# 2.3: If this distance is less than shortest_bondlength, record it in case we can not find a bond distance that means this requirement.
		if bondlength < shortest_bondlength:
			shortest_displacement = displacement
			shortest_bondlength = bondlength

	# Third, get the covalent bond distance between element1 and element2.
	covalent_bond_distance = get_covalent_bond_distance(element1, element2)

	# Fourth, if you are happy to take the shortest distance, do this.
	if take_shortest_distance:

		# 4.1: If the shortest distance between the two components is twice that of covalent_bond_distance, this is a problem and needs to be reported as a Exception.
		if shortest_bondlength > (2.0 * covalent_bond_distance):
			to_string  = 'The shortest bond distance is greater than two times the maximum covalent bond distance.\n'
			to_string += 'This potentially means something is up with your crystal structure.\n'
			to_string += 'Issue between atoms with indices: '+str(m_in_c_index1)+' and '+str(m_in_c_index2)+'.\n'
			crystal_debug = crystal.copy()
			crystal_debug[m_in_c_index1].symbol = 'Ar'
			crystal_debug[m_in_c_index2].symbol = 'Ar'
			from ase.visualize import view
			view([crystal,crystal_debug])
			to_string += 'Maximum bonds length between '+str(element1)+' and '+str(element2)+': '+str(covalent_bond_distance)+'\n'
			to_string += 'Shortest cell diplacement vector(s): '+str(shortest_displacement)+'\n'
			to_string += 'Shortest bond length '+str(shortest_bondlength)+' A.\n'
			to_string += 'Cell diplacement vector: Bond length (A)\n'
			for displacement in cell_points:
				bondlength = get_distance(position1, position2 + displacement)
				to_string += str(displacement)+': '+str(bondlength)+'\n'
			raise Exception(to_string)

		# 4.2: Report any warning issues to the logger if present.
		if logger is not None:
			to_string  = 'Had to take the shortest distance during making molecule method.\n'
			to_string += 'Maximum bonds length between '+str(element1)+' and '+str(element2)+': '+str(covalent_bond_distance)+'\n'
			to_string += 'Shortest cell diplacement vector(s): '+str(shortest_displacement)+'\n'
			to_string += 'Cell diplacement vector: Bond length (A)\n'
			for displacement in cell_points:
				bondlength = get_distance(position1, position2 + displacement)
				to_string += str(displacement)+': '+str(bondlength)+'\n'
			logger.warning(to_string)

		# 4.3: Return shortest_displacement
		return shortest_displacement

	# Fourth, if you do not want to take the shortest distance, return an error statement.
	print('Error in def connect_bonds_together, in component_assembly_approach.py')
	print('Could not determine where to place component 2 so that it connects with component 1')
	print('Check this problem out')
	print('Atom indices: '+str(m_in_c_index1)+' and '+str(m_in_c_index2))
	print('Maximum bonds length between '+str(element1)+' and '+str(element2)+': '+str(covalent_bond_distance))
	print('Cell diplacement vector: Bond length (A)')
	for displacement in cell_points:
		bondlength = get_distance(position1, position2 + displacement)
		print(str(displacement)+': '+str(bondlength))
	from ase.visualize import view
	view(molecule_in_crystal)
	import pdb; pdb.set_trace()
	print('This program will not finish without completing')
	exit()

def connect_bonds_together(molecule_in_crystal, disconnected_atom_pairs, crystal_cell_lattice, molecule_components_graphs):
	"""
	This method will reconstruct a molecule from the crystal structure by creating adding 26 identical cells around the origin unit cell, and keeping only the component that reconstructs the molecule.

	Parameters
	----------
	molecule_in_crystal : ase.Atoms
		This is the molecule in the crystal that you want to reconstruct.
	disconnected_atom_pairs : list
		This is a list of all the atoms that are disconnected in the crystal which need to be reconnected.
	crystal_cell_lattice : numpy.array
		A 3x3 row matrix giving the cell vectors in the i, j, and k basis set directions (in cartesian xyz coordinates). 
	molecule_components_graphs: list
		These are a list of indices of all the individual components of the disconnected molecule in the crystal.

	Returns
	-------
	complete_molecule : ase.Atoms
		This is the molecule that has been attached back together again.
	"""

	# First, get the points of displacement to move to lattices currounting the (0,0,0) cell
	cell_points = get_cell_corner_points(crystal_cell_lattice,super_cell_reach=1)

	# Second, reconnect each pair of atoms that are current disconnected, as given by the disconnected_atom_pairs list
	connected_molecule = molecule_in_crystal.copy()
	connected_molecule_graphs = [sorted(graph) for graph in sorted(molecule_components_graphs,key=lambda x: (len(x), -min(x)), reverse=True)]

	for (index1, index2), atom_1_to_atom_2_vector in disconnected_atom_pairs: 

		# 2.1: Determine what subgraphs (component graph) each atom in the disconnected bond comes from
		sg_index1 = None; sg_index2 = None; 
		subgraph1 = None; subgraph2 = None; 
		for mcg_index in range(len(connected_molecule_graphs)):
			subgraph = connected_molecule_graphs[mcg_index]
			if (subgraph1 is None) and (index1 in subgraph):
				subgraph1 = subgraph
				sg_index1 = mcg_index
			if (subgraph2 is None) and (index2 in subgraph):
				subgraph2 = subgraph
				sg_index2 = mcg_index
			if (subgraph1 is not None) and (subgraph2 is not None):
				break
		if (subgraph1 is None) or (subgraph2 is None):
			print('Error in def connect_bonds_together, in component_assembly_approach.py')
			print('A subgraph could not be identified, which should not happen')
			print('subgraph1: '+str(subgraph1))
			print('subgraph2: '+str(subgraph2))
			print('Check this out')
			import pdb; pdb.set_trace()
			print('This program will not finish without completing')
			exit()

		# 2.2: If sg_index1 == sg_index2, then this disconnected bond has been connected already in the
		#      process of connecting other parts of the molecule together
		if (sg_index1 == sg_index2):
			continue

		# 2.3: Using the atom_1_to_atom_2_vector vector, determine the displacement vector that is required 
		#      to move component 2 from where it is to where it needs to go to connect to component 1
		position1 = connected_molecule[index1].position
		position2_current = connected_molecule[index2].position
		position2_to_place = position1 + atom_1_to_atom_2_vector
		displacement_required_to_join_component_2_to_component_1 = position2_to_place - position2_current

		# 2.4: Move component 2 into place so that it connects with component 1.
		new_positions = connected_molecule.get_positions()
		for subgraph2_index in subgraph2:
			new_positions[subgraph2_index] += displacement_required_to_join_component_2_to_component_1
		connected_molecule.set_positions(new_positions)

		# 2.5: Connect the graphs of component 1 and 2 together, as we have now connected these two components together
		for mcg_index in sorted([sg_index1,sg_index2],reverse=True):
			del connected_molecule_graphs[mcg_index]
		connected_molecule_graphs.append(set(list(subgraph1)+list(subgraph2)))

	# Third, return the molecule that has now been reconnected. 
	return connected_molecule


