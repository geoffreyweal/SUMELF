"""
general_methods.py, Geoffrey Weal, 17/2/22

This script includes methods for reconstructing a molecule from the crystal structure using a super lattice approach.
"""
from ase import Atoms
from networkx import Graph, connected_components

from SUMELF.SUMELF.general_methods.distance_methods import less_than_or_equal_to_max_bondlength

def create_molecule_using_super_lattice_approach(molecule_in_crystal, crystal_graph, crystal, mapping_atom_in_crystal_TO_atom_in_molecule_indices, to_print=True):
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
	to_print : bool
		This tag will determine if you want to show a progress bar of atoms being added to the super_lattice, as this can take a bit fo time. Default: True.

	Returns
	-------
	complete_molecule : ase.Atoms
		This is the molecule that has been attached back together again.
	molecule_graph : networkx.graph
		This is the graph associated to the completed molecule.
	"""

	# First, get information about molecule and crystal.
	molecule_in_crystal_elements = molecule_in_crystal.get_chemical_symbols()
	molecule_in_crystal_no_of_atoms = len(molecule_in_crystal)
	crystal_cell_lattice = crystal.get_cell()
	crystal_elements = crystal.get_chemical_symbols()

	# Second, make the super-lattice containing only the molecule. This super-lattice consists of 27 lattices in total.
	cell_points = get_cell_corner_points(crystal_cell_lattice, super_cell_reach=1)
	no_of_cells = len(cell_points)
	super_molecule = Atoms()
	for displacement in cell_points:
		moved_molecule = molecule_in_crystal.copy()
		moved_molecule.set_positions(moved_molecule.get_positions() + displacement)
		super_molecule += moved_molecule

	# Third, make a graph of the super-lattice.
	super_molecule_graph = nx.Graph()
	super_molecule_graph.add_nodes_from(list((index, {'E': super_molecule[index].symbol}) for index in range(len(super_molecule))))
	molecule_graph = Graph() # make the molecule graph while we are doing things.
	molecule_graph.add_nodes_from(list((index, {'E': molecule_in_crystal[index].symbol}) for index in range(len(molecule_in_crystal))))
	indice_to_add = [ii*len(molecule_in_crystal) for ii in range(no_of_cells)]
	the_list = list(mapping_atom_in_crystal_TO_atom_in_molecule_indices.keys()) # mapping_atom_in_crystal_TO_atom_in_molecule_indices.keys() is equivalent to molecules_indices
	if to_print:
		the_list = tqdm(the_list)

	# Fourth, add edges to super_molecule_graph
	for index1 in the_list:
		element1 = crystal_elements[index1]
		m_in_c_index1 = mapping_atom_in_crystal_TO_atom_in_molecule_indices[index1]
		for index2 in crystal_graph[index1]:
			element2 = crystal_elements[index2]
			m_in_c_index2 = mapping_atom_in_crystal_TO_atom_in_molecule_indices[index2] # issue
			bondlength = super_molecule.get_distance(m_in_c_index1,m_in_c_index2,mic=True)
			if less_than_or_equal_to_max_bondlength(bondlength, element1, element2):
				molecule_graph.add_edge(molecule_index1, molecule_index2)
			for ii in range(len(indice_to_add)):
				new_m_in_c_index1 = m_in_c_index1+indice_to_add[ii]
				for jj in range(len(indice_to_add)):
					new_m_in_c_index2 = m_in_c_index2+indice_to_add[jj]
					bondlength = super_molecule.get_distance(new_m_in_c_index1,new_m_in_c_index2,mic=False)
					if less_than_or_equal_to_max_bondlength(bondlength, element1, element2):
						super_molecule_graph.add_edge(new_m_in_c_index1, new_m_in_c_index2)

	# Fifth, determine the components of the super-lattice, and record the component with the most number of atoms, as this should be a complete molecule
	sep_super_molecule_graph = list(nx.connected_components(super_molecule_graph))
	complete_molecule_graph = max(sep_super_molecule_graph,key=lambda x:len(x))
	if len(complete_molecule_graph) == len(molecule_in_crystal):
		print('Error in def create_molecule_using_super_lattice_approach, in super_lattice_approach.py')
		print('The completed component obtained from the super-lattice does not have the same number of atoms as molecule')
		print('This means there is something wrong with the program, as the completed component is not complete')
		print('Showing this "completed" component, and the molecule in the crystal.')
		from ase.visualise import view
		view([complete_molecule_graph,molecule_in_crystal])
		print('ckeck this out')
		import pdb; pdb.set_trace()
		print('This program will now stop without finishing')
		exit()

	# Sixth, obtain the indices of the completed molecule, and order then in the same way as they are ordered in the crystal.
	complete_molecule_indices = list(complete_molecule_graph)
	complete_molecule_indices.sort(key=lambda x: x%molecule_in_crystal_no_of_atoms)

	# Seventh, make the completed molecule.
	complete_molecule = Atoms()
	complete_molecule.set_cell(crystal_cell_lattice)
	for index in complete_molecule_indices:
		complete_molecule.append(super_molecule[index])

	# Seventh, make the molecule_graph from super_molecule_graph.
	# Need to work on this still
	import pdb; pdb.set_trace()
	raise Exception('Need to make a method for obtaining the molecule_graph from crystal_graph')

	# Eighth, return complete_molecule and molecule_graph
	return complete_molecule, molecule_graph


