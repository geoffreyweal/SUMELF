"""
process_crystal.py, Geoffrey Weal, 19/7/22

This script is designed to process the crystal from the crystal file
"""
import numpy as np
from copy import deepcopy

from ase import Atoms
from ase.io import read

from SUMELF.SUMELF.general_methods.general_molecules_methods                     import read_crystal
from SUMELF.SUMELF.graph_methods.obtain_graph                                    import obtain_graph
from SUMELF.SUMELF.graph_methods.get_properties_from_graph                       import get_node_property_names_from_graph

from SUMELF.SUMELF.process_crystal_methods.process_crystal_with_moleculeList     import process_crystal_with_moleculeList
from SUMELF.SUMELF.process_crystal_methods.process_crystal_with_spacegroup_kinds import process_crystal_with_spacegroup_kinds
from SUMELF.SUMELF.process_crystal_methods.process_crystal_with_default          import process_crystal_with_default

atoms_object = Atoms()

def process_crystal(filepath, crystal_graph=None, bonds_to_ignore=None, take_shortest_distance=False, return_list=False, no_of_cpus=1, logger=None, print_progress=True):
	"""
	This method is designed to obtain the individual molecules, as well as other quantities, from the crystal file.

	Parameters
	----------
	filepath : str. or ase.Atoms
		This is the path to the crystal file you want to process, or the ase.Atyoms object of the crystal.
	crystal_graph : networkx.graph or None
		This is the crystal_graph that has already been processed. If you want this method to process the graph for this crystal, set crystal_graph to None. Default: None
	bonds_to_ignore : list
		This is a list of the atom pairs (given as indices) to ignore bonding between. This is given as a list of index pairs. Default: None
	take_shortest_distance : bool
		If true, take the shortest distance that was found if a distance shorter than the maximum distance between elements if found. If false, give an error if the distance between two atoms was greater than the max distance expected for those two elements to be bonded to each other. Default: False
	return_list : bool.
		If true, return molecules, molecule_graphs as lists. If False, return molecules, molecule_graphs as dictionaries. Default: True
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
	logger : logging/None
		This contains the object for reporting information and warning messages to.
	print_progress : bool.
		Indicates if you want to write what the method is currently doing to the terminal. Default: True

	Returns
	-------
	molecules : list of ase.Atoms
		These are the unique molecules obtained from the crystal file (accounting for symmetry in the crystal's spacegroup).
	molecule_graphs : list of networkx.Graphs
		These are a list of all the graph s associated with each molecule in the molecules list
	solvent_molecule_names : list of ints
		This list indicates which of the molecules are likely solvents.
	symmetry_operations : list of (numpy.array(3*3), numpy.array(3*1))
		These are the rotation/reflection and translation matrices for each symmetry operation in the spacegroup/crystal.
	unit_cell_lattice : numpy.array
		These are the three vectors that describe the cell.
	"""

	# First, read in the crystal file in the ASE.
	if isinstance(filepath, str):
		crystal = read_crystal(filepath)
	elif type(filepath) == type(atoms_object):
		crystal = deepcopy(filepath)
	else:
		raise Exception('Error: The input crystal should be either the path to the crystal, or an ase.Atoms object of the crystal. crystal: '+str(crystal))

	# Second, get the graph of the crystal.
	if crystal_graph is None:
		crystal, crystal_graph = obtain_graph(crystal,name='crystal')

	# Third, determine what features are constained in the file.
	has_MoleculeList       = ('MoleculeList'     in get_node_property_names_from_graph(crystal_graph)) # if (crystal_graph is not None) else False
	has_spacegroup_kinds   = ('spacegroup_kinds' in get_node_property_names_from_graph(crystal_graph)) # if (crystal_graph is not None) else False

	# Fourth, obtain all the unique molecules from the crystal (that are not repeated due to symmetry operations of the spacegroup).
	if   has_MoleculeList:

		# 4.1: Take advantage of the given spacegroup and MoleculeList for this crystal file.
		if print_progress:
			print('Processing Crystal using MoleculeList')
		molecules, molecule_graphs, solvent_molecule_names, symmetry_operations = process_crystal_with_moleculeList(crystal, crystal_graph, bonds_to_ignore=bonds_to_ignore, take_shortest_distance=take_shortest_distance, make_molecule_method='component_assembly_approach', no_of_cpus=no_of_cpus, logger=logger, print_progress=print_progress)

	elif has_spacegroup_kinds:

		# 4.2: Take advantage of the given spacegroup and spacegroup_kinds for this crystal file.
		if print_progress:
			print('Processing Crystal using Spacegroup Kinds')
		'''
		raise Exception('This method needs checking. Update for BondProperties')
		print('Will run but need to check this')
		raise Exception('update solvent_molecule_names in process_crystal_with_default')
		raise Exception('Give molecules and graphs as a dictionary by their names, not their indices')
		'''
		molecules, molecule_graphs, solvent_molecule_names, symmetry_operations = process_crystal_with_default(crystal, crystal_graph=crystal_graph, bonds_to_ignore=bonds_to_ignore, take_shortest_distance=take_shortest_distance, make_molecule_method='component_assembly_approach', no_of_cpus=no_of_cpus, logger=logger, print_progress=print_progress)
		#molecules, molecule_graphs, solvent_molecule_names, symmetry_operations = process_crystal_with_spacegroup_kinds(crystal, crystal_graph=crystal_graph, take_shortest_distance=take_shortest_distance, 

	else:

		# 4.3: Assume that every atom in the crystal is contained in the unit cell of the input file.
		#      If a MoleculeList is given, can take advantage of this. 
		#      This method can be used whether the crystal is provided with a spacegroup or not.
		if print_progress:
			print('Processing Crystal using Default Method')
		#raise Exception('Make sure molecules are created in the order expected from the lowest atom indicxe in each molecule.')
		#raise Exception('Give molecules and graphs as a dictionary by their names, not their indices')
		molecules, molecule_graphs, solvent_molecule_names, symmetry_operations = process_crystal_with_default(crystal, crystal_graph=crystal_graph, bonds_to_ignore=bonds_to_ignore, take_shortest_distance=take_shortest_distance, make_molecule_method='component_assembly_approach', no_of_cpus=no_of_cpus, logger=logger, print_progress=print_progress)

	# Fifth, remove any solvents that are not in molecules (due to being being structurally the same due to the symmetry of the crystal).
	for index, solvent_mol_name in enumerate(solvent_molecule_names):
		if not solvent_mol_name in molecules.keys():
			del solvent_molecule_names[index]

	# Sixth, if you would like to return the molecules and molecule_graphs as lists rather than dictionaries, do this here
	if return_list:
		molecules       = [molecule       for ii, molecule       in sorted(molecules.items(),      key=lambda x:x[0])]
		molecule_graphs = [molecule_graph for ii, molecule_graph in sorted(molecule_graphs.items(),key=lambda x:x[0])]

	# Seventh, return the molecules, molecule_graphs, solvent_molecule_names, symmetry_operations, and the unit_cell_lattice
	return molecules, molecule_graphs, solvent_molecule_names, symmetry_operations, crystal.get_cell()

# ---------------------------------------------------------------------------------------------------------------------------
	

