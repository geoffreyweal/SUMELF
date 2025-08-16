"""
obtain_graph.py, Geoffrey Weal, 17/2/22

This script contains methods for obtaining a graph for a chemical system, be in a crystal structure or molecule.
"""
import sys
from networkx import Graph
import multiprocessing as mp

from tqdm import tqdm

from SUMELF.SUMELF.general_methods.distance_methods     import less_than_or_equal_to_max_bondlength
from SUMELF.SUMELF.graph_methods.NeighboursList_methods import get_neighbours_from_NeighboursList, convert_NeighboursList_to_bonds
from SUMELF.SUMELF.graph_methods.BondProperties_methods import get_neighbours_from_BondProperties

def obtain_graph(original_system,mic=True,to_print=False,name='',no_of_cpus=1,bonds_to_ignore=None,logger=None):
	"""
	This method will give a graph for a chemical system, be in a crystal structure or molecule.

	Parameters
	----------
	original_system : ase.Atoms
		The chemical system in Atomic Simulation Environment format.
	mic : bool
		Use the Minimum Image Convention to give the bond distance, even if the molecule has moved off the edge of the crystal on to the other side of the cystal due to the periodic boundary conditions.
	to_print : bool
		Print a progress bar for atoms being added to the graph.
	name : str.
		This is the name of the molecule.
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
	bonds_to_ignore : list
		This is a list of the atom pairs (given as indices) to ignore bonding between. This is given as a list of index pairs. Default: None
	logger : logging/None
		This contains the object for reporting information and warning messages to.

	Returns
	-------
	system : ase.atoms
		This is the original ASE Atoms object, but with NeighboursList and BondProperties removed to prevent confusion when changes are made with bonding indices
	system_graph : networkx.Graph
		This is the graph for the chemical system.
	"""

	# First, set up the system, graph and the bonds_to_ignore list
	system = original_system.copy()
	system_graph = Graph(name=name)
	if bonds_to_ignore is None:
		bonds_to_ignore = []

	# Second, get atom properties that are important for the graph. 
	all_atom_properties = {}
	for atom_property in system.arrays.keys():
		if atom_property in ['numbers', 'positions', 'initial_charges', 'NeighboursList']:
			continue
		all_atom_properties[atom_property] = system.arrays[atom_property]

	# Third, add atoms as node on the Graph.
	all_nodes = []
	for atom_index, atom in enumerate(system):
		atom_properties = {'E': system[atom_index].symbol}
		for atom_property_name, atom_property in all_atom_properties.items():
			atom_properties[atom_property_name] = atom_property[atom_index]
		all_nodes.append((atom_index, atom_properties))

	# Fourth, add all_nodes as nodes to system_graph
	system_graph.add_nodes_from(all_nodes)

	# Fifth, remove non-standard arrays from ase Atoms object. 
	#        This is to:
	#        * Prevent confusion between double information in the ase.Atoms object and the networkx.graph object, and
	#        * Because these pieces of information are not given to the ase.Atom object when extracting atom information from an ase.Atoms object.
	for atom_property_name in all_atom_properties.keys():
		del system.arrays[atom_property_name]

	# --------------------------------------------------------------------------------------------------

	# Sixth, determine if we have NeighboursList and BondProperties
	have_NeighboursList = 'NeighboursList' in system.arrays.keys()
	have_BondProperties = 'BondProperties' in system.info.keys()

	# Seventh, add bonds as edges to the graph.
	# 		 * Add bonds between atoms as edges in graph if the distance between atoms are less 
	# 		   than some maximum bond distance given by the less_than_or_equal_to_max_bondlength method.
	if have_NeighboursList or have_BondProperties:
		# 7.1: Have information about the bonds in your molecule. Use these to create you graphs. 

		# 7.1.1: Obtain the bonds from the NeighboursList
		if have_NeighboursList:
			NeighboursList            = get_neighbours_from_NeighboursList(system.arrays['NeighboursList'])
			bonds_from_NeighboursList = convert_NeighboursList_to_bonds(NeighboursList)
		
		# 7.1.2: Obtain the bonds from the BondProperties
		if have_BondProperties:
			BondProperties = eval(system.info['BondProperties'])
			bonds_from_BondProperties, BondProperties_bondasKey = get_neighbours_from_BondProperties(BondProperties)

		# ----------------------------------------------------------------------------------------------

		# 7.1.3: Convert NeighboursList/BondProperties into all_edges
		if   have_NeighboursList and have_BondProperties:

			# 7.1.3.1: Check that bonds_from_NeighboursList and bonds_from_BondProperties are the same (they should be).
			if not (bonds_from_NeighboursList == bonds_from_BondProperties):
				raise Exception('Error: NeighboursList and BondProperties are different: '+str(set(bonds_from_NeighboursList) ^ set(bonds_from_BondProperties)))

			# 7.1.3.2: Set all_edges = bonds_from_NeighboursList. Note: bonds_from_NeighboursList and bonds_from_BondProperties are the sane, so can choose either to be all_edges.
			all_edges = bonds_from_NeighboursList

		elif have_NeighboursList:

			# 7.1.3.3: Set all_edges = bonds_from_NeighboursList.
			all_edges = bonds_from_NeighboursList

		elif have_BondProperties:

			# 7.1.3.4: Set all_edges = bonds_from_BondProperties.
			all_edges = bonds_from_BondProperties

		else:

			# 7.1.3.5: Should not get here in the code, as "have_NeighboursList or 
			#          "have_BondProperties" should prevent the program from getting here
			to_string  = 'Error: You do not seem to have "NeighboursList" or "BondProperties" in your system.\n'
			to_string += 'Check your crystal files, as these may be missing.\n'
			to_string += 'If not, this may be a programming error\n'
			raise Exception(to_string)

		# ----------------------------------------------------------------------------------------------

		# 7.1.4: Include only edges that dont exist in the bonds_to_ignore list
		all_edges_except_bonds_to_ignore = []
		for atom1_index, atom2_index in all_edges:

			# 7.1.4.1: Check that atom1_index == atom2_index, as this would indicate that the atom is bonded to itself.
			if atom1_index == atom2_index:
				raise Exception(f'Error: atom {atom1_index} seems to be bonded to itself? This should be the case.')

			# 7.1.4.2: Check that this bond hasn't been included in bonds_to_ignore
			edge         = (atom1_index, atom2_index)
			reverse_edge = (atom2_index, atom1_index)
			if (edge in bonds_to_ignore) or (reverse_edge in bonds_to_ignore):
				continue

			# 7.1.4.3: Obtain the bond properties for this bond (if it exists in BondProperties_bondasKey)
			if have_BondProperties:
				bond_properties = BondProperties_bondasKey.get(tuple(sorted(edge)), {})
			else:
				bond_properties = {}

			# 7.1.4.4: Add bond as edge in all_edges_except_bonds_to_ignore
			all_edges_except_bonds_to_ignore.append([atom1_index, atom2_index, bond_properties])

		# ----------------------------------------------------------------------------------------------

		# 7.1.5: Add all_edges list to system_graph networkx graph.
		system_graph.add_edges_from(all_edges_except_bonds_to_ignore)

		# 7.1.6: Remove the 'NeighboursList' from system.arrays, to avoid confusion with the system_graph graph
		if have_NeighboursList:
			del system.arrays['NeighboursList']

		# 7.1.7: Remove the 'BondProperties' from system.arrays, to avoid confusion with the system_graph graph
		if have_BondProperties:
			del system.info['BondProperties']

	else:
		# 7.2: COULD NOT FIND have_NeighboursList or have_BondProperties. So will try to manually obtain the bonded between atoms based on distance between atoms in your system.

		# 7.2.1: If neighbour list is not given, get the networkx graph using this method.
		elements = system.get_chemical_symbols()
		no_of_calculations = int((len(system) * (len(system) - 1)) / 2)
		if to_print:
			print('Creating '+str(name)+' graph: ')
			print('# of atoms in '+str(name)+': '+str(len(system)))
			print('Will perform 1+2+...+'+str(len(system))+' = '+str(no_of_calculations)+' dist. measurements between atoms in '+str(name)+'.')

		# 7.2.2: Collect the edges and bond them to the system_graph in a series or parallel fashion
		if no_of_cpus > 1:
			if logger is not None:
				logger.info('Obtining the graph based on atom-atom distances using map_async in python multiprocessing')
			with mp.Manager() as manager: 
				edges = manager.list()
				pool = mp.Pool(no_of_cpus)
				if to_print:
					print('Creating graphs (Please wait until after 100%, as the process will still be running.)', file=sys.stderr)
					pool.map_async(get_connected_atoms_in_system, tqdm(get_inputs(system, elements, mic, bonds_to_ignore, edges), total=no_of_calculations, desc='Creating '+str(name)+' graph', unit='calc'))
				else:
					pool.map_async(get_connected_atoms_in_system, get_inputs(system, elements, mic, bonds_to_ignore, edges))
				pool.close()
				if logger is not None:
					logger.info('joining the pool')
				pool.join()
				if logger is not None:
					logger.info('Number of Edges to add to networkx graph: '+str(len(edges)))
				system_graph.add_edges_from(list(edges))
			if logger is not None:
				logger.info('Edges have been added to networkx graph')
		else:
			edges = []
			if to_print:
				for input_variables in tqdm(get_inputs(system, elements, mic, bonds_to_ignore, edges), total=no_of_calculations, desc='Creating '+str(name)+' graph', unit='calc'):
					get_connected_atoms_in_system(input_variables)
			else:
				for input_variables in get_inputs(system, elements, mic, bonds_to_ignore, edges):
					get_connected_atoms_in_system(input_variables)	
			system_graph.add_edges_from(edges)

	# --------------------------------------------------------------------------------------------------

	# Eight, return the graph for the system (crystal or molecule).
	#return system, system_graph
	return system_graph

# ===============================================================================================================================================

def get_inputs(system, elements, mic, bonds_to_ignore, edges):
	"""
	This method is designed to provide the inputs required for obtaining the bonding length.

	Parameters
	----------
	system : ase.Atoms
		This is the system (crystal) to obtain information about.
	elements : list of str.
		This is a list of the elements in the system (crystal)>
	mic : bool
		Use the Minimum Image Convention to give the bond distance, even if the molecule has moved off the edge of the crystal on to the other side of the cystal due to the periodic boundary conditions.
	bonds_to_ignore : list
		This is a list of the atom pairs (given as indices) to ignore bonding between. This is given as a list of index pairs. Default: None
	edges : list or manager.list
		This is the list of edges to record data to

	Returns
	-------
	atom1_index : int
		This is the index of the first atom to examine.
	atom2_index : int
		This is the index of the second atom to examine.
	element1 : str.
		This is the element of the first atom to examine.
	element2 : str.
		This is the element of the second atom to examine.
	system : ase.Atoms
		This is the system to obtain information about.
	mic : bool
		Use the Minimum Image Convention to give the bond distance, even if the molecule has moved off the edge of the crystal on to the other side of the cystal due to the periodic boundary conditions.
	"""

	# First, for the first atom.
	for atom1_index in range(len(system)):

		# Second, get the element of the first atom.
		element1 = elements[atom1_index]

		# Third, for the second atom.
		for atom2_index in range(atom1_index+1, len(system)):

			if ((atom1_index, atom2_index) in bonds_to_ignore) or ((atom2_index, atom1_index) in bonds_to_ignore):
				continue

			# Fourth, get the element of the second atom.
			element2 = elements[atom2_index]

			# Fifth, return variables for calculations.
			yield (atom1_index, atom2_index, element1, element2, system, mic, edges)

def get_connected_atoms_in_system(input_variables):
	"""
	This method is designed to determine which atoms are connected to each other in the crystal

	Parameters
	----------
	atom1_index : int
		This is the index of the first atom to examine.
	atom2_index : int
		This is the index of the second atom to examine.
	element1 : str.
		This is the element of the first atom to examine.
	element2 : str.
		This is the element of the second atom to examine.
	system : ase.Atoms
		This is the system to obtain information about.
	mic : bool
		Use the Minimum Image Convention to give the bond distance, even if the molecule has moved off the edge of the crystal on to the other side of the cystal due to the periodic boundary conditions.

	Returns
	-------
	If True, the atoms are connect. If False, atoms are not connected
	"""

	# First, obtain the input variables.
	atom1_index, atom2_index, element1, element2, system, mic, edges = input_variables

	# Second, obtain the bond length between atom 1 and atom 2
	bondlength = system.get_distance(atom1_index,atom2_index,mic=mic)

	# Third, determine if the atom is within bonding distance for the atoms elements.
	is_connected = less_than_or_equal_to_max_bondlength(bondlength, element1, element2)

	# Fourth, return the result of is_connected
	if is_connected:
		edges.append((atom1_index, atom2_index))

# ===============================================================================================================================================



