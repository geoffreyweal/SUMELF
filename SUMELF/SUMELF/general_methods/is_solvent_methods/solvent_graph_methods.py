"""
solvent_graph_methods.py, Geoffrey Weal, 15/2/24

This script contains various methods for reading and writing graph methods for determining solvents. 
"""
import os
import networkx as nx
from importlib.util import find_spec

# These variables indicate the folders where solvent files are placed.
general_path_to_solvent_mol_files   = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'solvent_files', 'solvent_mol_files'))
general_path_to_solvent_graph_files = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'solvent_files', 'solvent_graph_files'))

def create_graph_from_mol_files():
	"""
	This method is designed to convert mol files to graphs.
	"""

	# First, determine if you have the CCDC Python API.
	ccdc_spec = find_spec("ccdc")
	if ccdc_spec is None:
		toString = ''
		toString += '\n'
		toString += '================================================'+'\n'
		toString += 'This is the SUMELF (Supporting Methods for Electronic Functions) Program'+'\n'
		toString += 'Version: '+str(__version__)+'\n'
		toString += '\n'
		toString += 'The SUMELF program requires the CCDC Python API.'+'\n'
		toString += '\n'
		toString += 'Install the CCDC Python API by following the instructions at https://downloads.ccdc.cam.ac.uk/documentation/API/installation_notes.html#installation.'+'\n'
		toString += '\n'
		toString += 'This program will exit before beginning'+'\n'
		toString += '================================================'+'\n'
		raise ImportError(toString)

	# Second, load the MoleculeReader from the ccdc program
	from ccdc.io import MoleculeReader

	# Third, check if the folder containing the mol files exists.
	if not os.path.exists(general_path_to_solvent_mol_files):
		raise Exception('Error: Could not find the folder containing solvent mol files.\nThe path to the folder that doesnt exist is: '+str(general_path_to_solvent_mol_files)+'\nCheck this.')

	# Fourth, create the folder to save graph data into if it doesnt exist.
	if not os.path.exists(general_path_to_solvent_graph_files):
		os.makedirs(general_path_to_solvent_graph_files)

	# Fifth, for each mol file in general_path_to_solvent_mol_files
	for file in os.listdir(general_path_to_solvent_mol_files):

		# Sixth, determine if the file is a mol file
		if not file.endswith('.mol'):
			continue

		# Seventh, read the mol file.
		molecule_data = MoleculeReader(general_path_to_solvent_mol_files+'/'+file)[0]

		# Eighth, create the networkx.Graph solvent object
		solvent = make_solvent_graph(molecule_data)

		# Ninth, save the graph file for this mol file into general_path_to_solvent_graph_files
		write_graph(general_path_to_solvent_graph_files+'/'+file.replace('.mol','.txt'), solvent)

def make_solvent_graph(molecule_data):
	"""
	This method is designed to create the graph file for the solvent ccdc.molecule.Molecule object. 

	Parameters
	----------
	molecule_data : ccdc.molecule.Molecule
		This is the ccdc molecules object of the solvent from the mol file

	Returns
	-------
	solvent_graph : networkx.Graph
		This is the networkx graph object of the solvent. This graph object doe not have hydrogens as nodes, but instead hydrogens are bound to the "heavy" atoms as node variables. 
	"""

	# First, obtain the information about the atoms from the molecule_data object.
	#        * Also, record the indices of hydrogens.
	atoms = {}; hydrogen_indices = []
	for atom in molecule_data.atoms:

		# 1.1: Obtain the index and element of the current atom
		index   = atom.index
		element = atom.atomic_symbol

		# 1.2: Make sure that index does not already exist in the atoms.keys() list.
		if index in atoms.keys():
			raise Exception('Error: index '+str(index)+' already in atoms dictionary. This may mean there is a duplicate index in molecule_data. Check this.')

		# 1.3: Determine if the atom is a hydrogen
		if element in ['H', 'D', 'T']:
			hydrogen_indices.append(index)
			continue

		# 1.4: Make sure that index in not in the atoms dictionary already. 
		if index in atoms.keys():
			raise Exception('Huh?')

		# 1.5: Append atom information to the atoms dictionary. 
		atoms[index] = {'E': element, 'H': 0}

	# Second, obtain the information about the bonds from the molecule_data object.
	bonds = []
	for bond in molecule_data.bonds:

		# 2.1: Get the indices of the atoms in the bond
		atom_indices_in_bond = bond.atoms
		atom_1_index = atom_indices_in_bond[0].index
		atom_2_index = atom_indices_in_bond[1].index

		# 2.2: Make sure that the atom is not bonded to itself (which would not make sense).
		if atom_1_index == atom_2_index:
			raise Exception('Huh?')

		# 2.3: Add hydrogens to "heavy" atom node variables.
		if atom_1_index in hydrogen_indices:
			atoms[atom_2_index]['H'] += 1
			continue
		if atom_2_index in hydrogen_indices:
			atoms[atom_1_index]['H'] += 1
			continue

		# 2.4: Make the bond as a tuple as (low value index, high value index) between two "heavy" atoms.
		bond = tuple(sorted([atom_1_index, atom_2_index]))

		# 2.5: Append the bond to the bonds list. 
		bonds.append(bond)

	# Third, make the solvent graph object
	solvent_graph = nx.Graph()

	# Fourth, add the atoms information to solvent_graph
	solvent_graph.add_nodes_from(atoms.items())

	# Fifth, add the bonds information to solvent_graph
	solvent_graph.add_edges_from(bonds)

	# Sixth, remap the solvent graph so that the indices are consecutive. 
	print('Check the below is working as expected.')
	import pdb; pdb.set_trace()
	remap_atoms = {old_atom_index: new_atom_index for (new_atom_index, old_atom_index) in enumerate(sorted(solvent_graph.nodes.keys()))}

	# Seventh, remap the nodes in solvent_graph with consecutive indices. 
	solvent_graph = nx.relabel_nodes(solvent_graph, remap_atoms)

	# Eighth, return solvent_graph.
	return solvent_graph

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def read_graph(path_to_graph_file):
	"""
	Read the graph file to memory.

	Parameters
	----------
	path_to_graph_file : str.
		This is the path to the graph file. 

	Returns
	-------
	system_graph : networkx.Graph
		This is the netowrkx.Graph object for this graph file. 
	"""

	# First, initialise the atoms dictionary and bonds list.
	atoms = {}
	bonds = []

	# Second, begin reading in graph file. 
	with open(path_to_graph_file, 'r') as GRAPH_FILE:

		# Third, initialise the getting_atoms boolean and atom_index. 
		getting_atoms = True
		atom_index = 0

		# Fourth, for each line in GRAPH_FILE
		for line in GRAPH_FILE:

			# Fifth, remove the end of line spaces.
			line = line.rstrip()

			# Sixth, if there is a line space in the file, this indcates we have finished reading the atoms information and now reading the bonds information. 
			if line == '':
				getting_atoms = False
				continue

			# Seventh, record the atoms or bonds information. 
			if getting_atoms:

				# ATOMS INFORMATION

				# 7.1: Read the atom information from the line, and record variables as their most appropriate type. 
				element, number_of_bound_hydrogen = line.split(':')
				number_of_bound_hydrogen = int(number_of_bound_hydrogen)

				# 7.2: Record the information about the atoms into a dictionary
				atom_info = {'E': element, 'H': number_of_bound_hydrogen}

				# 7.3: Add the atoms information to the atoms dictionary.
				atoms[atom_index] = atom_info

				# 7.4: Increment the atom index
				atom_index += 1

			else:

				# ATOMS INFORMATION

				# 7.5: Get the information about the bond, and record variables as their most appropriate type.  
				atom1_index, atom2_index = line.split()
				atom1_index = int(atom1_index)
				atom2_index = int(atom2_index)

				# 7.6: Add the bond information to the bonds list.
				bonds.append((atom1_index, atom2_index))

	# Third, initialise the graphs object
	system_graph = nx.Graph()

	# Fourth, add the atoms information to system_graph as nodes.
	system_graph.add_nodes_from(atoms.items())

	# Fifth, add the bonds information to system_graph as edges.
	system_graph.add_edges_from(bonds)

	# Sixth, return solvent_graph.
	return system_graph

def write_graph(path_to_graph_file, system_graph):
	"""
	Write the graph file to disk.

	Parameters
	----------
	path_to_graph_file : str.
		This is the path to the graph file. 
	system_graph : networkx.Graph
		This is the netowrkx.Graph object for this graph file. 
	"""

	# First, make sure that the indices of system_graph object is consecutive from 0 to len(system_graph)-1.
	if not sorted(system_graph.nodes.keys()) == list(range(len(system_graph))):
		raise Exception('Error: Graph object node indices is not consecutive from 0 to len(system_graph)-1.\nCheck this\nNode indices = '+str(sorted(system_graph.nodes.keys())))

	# Second, open the graph file.
	with open(path_to_graph_file, 'w') as GRAPH_FILE:

		# Third, record each atom in the GRAPH_FILE file. 
		for atom_index, atom_details in system_graph.nodes.items():
			GRAPH_FILE.write(str(atom_details['E'])+': '+str(atom_details['H'])+'\n')

		# Fourth, add a empty line to indicate we are now recording bond information. 
		GRAPH_FILE.write('\n')

		# Fifth, record each bond in the GRAPH_FILE file. 
		for atom1_index, atom2_index in system_graph.edges.keys():
			atom1_index, atom2_index = sorted((atom1_index, atom2_index))
			GRAPH_FILE.write(str(atom1_index)+' '+str(atom2_index)+'\n')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 



