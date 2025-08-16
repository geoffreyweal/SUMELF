"""
process_crystal_with_spacegroup_kinds.py, Geoffrey Weal, 2/6/2022

This script is designed to remove aliphatic sidechains from your molecules in the crystal file if you know the space group for the crystal object.
"""

from copy import deepcopy

from networkx import connected_components

from ase import Atoms
from ase.spacegroup import Spacegroup

from SUMELF.SUMELF.process_crystal_methods.process_crystal_with_default import process_crystal
#from SUMELF.SUMELF.process_crystal_methods.reconnect_broken_molecules import reconnect_broken_molecules # Not sure if needed anymore.

def process_crystal_with_spacegroup_kinds(crystal, crystal_graph=None, bonds_to_ignore=None, take_shortest_distance=False, make_molecule_method='component_assembly_approach', no_of_cpus=1, logger=None):
	"""
	This method is designed to remove aliphatic sidechains from your molecules in the crystal file if you know the space group for the crystal object.

	Parameters
	----------
	crystal : ase.Atoms
		This is ase object of the crystal. 
	crystal_graph : networkx.graph or None
		This is the crystal_graph that has already been processed. If you want this method to process the graph for this crystal, set crystal_graph to None. Default: None
	bonds_to_ignore : list
		This is a list of the atom pairs (given as indices) to ignore bonding between. This is given as a list of index pairs. Default: None
	take_shortest_distance : bool
		If true, take the shortest distance that was found if a distance shorter than the maximum distance between elements if found. If false, give an error if the distance between two atoms was greater than the max distance expected for those two elements to be bonded to each other. Default: False
	make_molecule_method : str.
		This is the name of the method you want to use to create the molecule. See https://github.com/geoffreyweal/ECCP for more information. Default: 'component_assembly_approach'. 
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
	logger : logging/None
		This contains the object for reporting information and warning messages to.

	Returns
	-------
	molecules : list of ase.Atoms
		These are the molecules that make up the crystal (before any spacegroup symmetry is considered)
	molecule_graphs : list of networkx.Graph
		These are the graph that are associated with the molecules in the molecules list. 
	wrapped_molecules : list of ase.Atoms
		These are the molecules in the molecules list, but stilled wrapped.
	wrapped_molecule_graphs : list of networkx.Graph
		These are the networkx.Graphs of the associated wrapped molecule graphs
	SolventsList : list of ints
		This list indicates which molecules in the molecules list are solvents
	symmetry_operations : list of (numpy.array(3*3), numpy.array(3*1))
		These are the rotation/reflection and translation matrices for each symmetry operation in the espacegroup.
	"""

	raise Exception('Error: Check if working correctly')

	raise Exception('Make sure that molecules are numbered based in order of the lowest index atom of the molecule in the crystal. Also make sure atoms in each molecule are ordered like that in the crystal.')
	import pdb; pdb.set_trace()

	# First, get all the spacegroup kinds in the crystal file.
	spacegroup_kinds = crystal.arrays['spacegroup_kinds']

	# Second, separate all the atoms that can be mapped onto each other by the spacegroup and place them into each of their components. 
	components = []
	crystal_counter = {}
	for index in range(len(crystal)):

		# 2.1: Get the spacegroup kind of the current atom
		spacegroup_kind = spacegroup_kinds[index]

		# 2.2: Add more Atoms objects to components if new components need to be created.
		crystal_counter[spacegroup_kind] = crystal_counter.get(spacegroup_kind,0)+1
		while len(components) < crystal_counter[spacegroup_kind]:
			component = Atoms(cell=crystal.get_cell(), pbc=True)
			components.append(component)

		# 2.3: Add the current atom to the component is should go with.
		component_index = crystal_counter[spacegroup_kind]-1
		atom = deepcopy(crystal[index])
		atom.tag = spacegroup_kind
		components[component_index].append(atom)

	# ==================================================================================================================================

	#import pdb; pdb.set_trace()
	#raise Exception('check here') # translate_molecules_for_symmetry_operations

	# Fourth, turn all the lists in the wrapped_molecules_indices_dict dictionary into tuples
	'''
	for key in wrapped_molecules_indices_dict.keys():
		wrapped_molecules_indices_dict[key] = tuple(wrapped_molecules_indices_dict[key])
	'''

	# Fifth, append the neighbourlist to the molecules in the wrapped_molecules_dict.
	'''
	has_NeighboursList = 'NeighboursList' in crystal.arrays.keys()
	if has_NeighboursList:
		for mol_index in wrapped_molecules_dict.keys():
			wrapped_molecules_dict[mol_index].arrays['NeighboursList'] = np.array([crystal.arrays['NeighboursList'][atom_index] for atom_index in wrapped_molecules_indices_dict[mol_index]])
	'''

	translate_molecules_for_symmetry_operations = {}

	# ==================================================================================================================================

	# Third, take only the first component. 
	component = components[0]

	# Fourth, process the first component to obtain the molecules it contains.
	wrapped_molecules, wrapped_molecule_graphs = process_crystal(component, translate_molecules_for_symmetry_operations, name='component 1', make_molecule_method=make_molecule_method, no_of_cpus=no_of_cpus, bonds_to_ignore=bonds_to_ignore, logger=logger)

	# Fifth, get the symmetry operations of the spacegroup for the crystal. 
	# If the spacegroup is not given, it is assumed to be 'P 1' (i.e. all atoms in the crystal has been given in the unit cell.)
	has_SymOpts    = 'SymOpts'    in crystal.info.keys()
	has_spacegroup = 'spacegroup' in crystal.info.keys()
	if   has_SymOpts:
		symmetry_operations = get_symmetry_operations(crystal.info['SymOpts'].split())
	elif has_spacegroup:
		spacegroup = Spacegroup(spacegroup=crystal.info['spacegroup'])
		symmetry_operations = spacegroup.get_symop()
	else:
		# Assume P 1 symmetry 
		spacegroup = Spacegroup(spacegroup='P 1')
		symmetry_operations = spacegroup.get_symop()

	# Sixth, some molecules may be broken due to symmetry within the unit cell for this space group.
	# This method will reconstruct any broken molecules.
	molecules = {}
	molecule_graphs = {}
	for index in range(len(wrapped_molecules)):
		molecule, molecule_graph = reconnect_broken_molecules(wrapped_molecules[index], wrapped_molecule_graphs[index], symmetry_operations=symmetry_operations)
		molecules[index] = molecule
		molecule_graphs[index] = molecule_graph

	# Seventh, get the SolventsList from crystal file if it has it.
	SolventsList = crystal.info['SolventsList'] if ('SolventsList' in crystal.info.keys()) else []
	if not ((type(SolventsList) == type(np.array([]))) or isinstance(SolventsList,list)):
		SolventsList = np.array([int(SolventsList)])

	raise Exception('Check')

	# Eighth, return the molecules and their molecule_graphs
	return molecules, molecule_graphs, SolventsList, symmetry_operations

# ---------------------------------------------------------------------------------------------------------------------------





