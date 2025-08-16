"""
process_crystal_with_default.py, Geoffrey Weal, 5/6/2022

This script is designed to process the crystal given into the molecules it is made up. 
"""
import numpy as np
from copy import deepcopy

from networkx import connected_components

from ase import Atoms
from ase.spacegroup import Spacegroup

from SUMELF.SUMELF.general_methods.get_symmetry_operations import get_symmetry_operations
from SUMELF.SUMELF.graph_methods.obtain_graph              import obtain_graph
from SUMELF.SUMELF.process_crystal_methods.make_molecule   import make_molecule
from SUMELF.SUMELF.process_crystal_methods.utilities       import get_solvent_molecule_names_from_SolventsList

def process_crystal_with_default(crystal, crystal_graph=None, bonds_to_ignore=None, produce_all_molecules=False, take_shortest_distance=False, make_molecule_method='component_assembly_approach', no_of_cpus=1, logger=None, print_progress=False):
	"""
	This method is designed to process the crystal given into the molecules it is made up. 

	Note if you are using this method, we assume that the user has included every atom that the crystal contains in the unit cell, as the spacegroup is not given. 

	Parameters
	----------
	crystal : ase.Atoms
		This is ase object of the crystal. 
	crystal_graph : networkx.graph or None
		This is the crystal_graph that has already been processed. If you want this method to process the graph for this crystal, set crystal_graph to None. Default: None
	bonds_to_ignore : list
		This is a list of the atom pairs (given as indices) to ignore bonding between. This is given as a list of index pairs. Default: None
	produce_all_molecules : bool.
		This boolean indicates if you want to return all the molecules in the crystal, or just the ones that can be replicated using the symmetry operations. Default: False
	take_shortest_distance : bool
		If true, take the shortest distance that was found if a distance shorter than the maximum distance between elements if found. If false, give an error if the distance between two atoms was greater than the max distance expected for those two elements to be bonded to each other. Default: False
	make_molecule_method : str.
		This is the name of the method you want to use to create the molecule. See https://github.com/geoffreyweal/ECCP for more information. Default: 'component_assembly_approach'. 
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
	logger : logging/None
		This contains the object for reporting information and warning messages to.
	print_progress : bool.
		Indicates if you want to write what the method is currently doing to the terminal. Default: False

	Returns
	-------
	molecules : list of ase.Atoms
		These are the molecules that make up the crystal (before any spacegroup symmetry is considered).
	molecule_graphs : list of networkx.Graph
		These are the graph that are associated with the molecules in the molecules list. 
	SolventsList : list of ints
		This list indicates which molecules in the molecules list are solvents.
	symmetry_operations : list of (numpy.array(3*3), numpy.array(3*1))
		These are the rotation/reflection and translation matrices for each symmetry operation the crystal contans. 
	"""

	# First, get the symmetry operations of the spacegroup for the crystal. 
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

	# Second, get the molecules and the graphs associated with each molecule in the crystal.
	molecules, molecule_graphs = process_crystal(crystal, name='crystal', take_shortest_distance=take_shortest_distance, make_molecule_method=make_molecule_method, no_of_cpus=no_of_cpus, bonds_to_ignore=bonds_to_ignore, logger=logger)

	# Third, obtain the list of molecules that are solvent in molecules
	solvent_molecule_indices = get_solvent_molecule_names_from_SolventsList(crystal, molecules)

	# Fourth, return the new crystal structure that has had it's aliphatic sidechains removed from non-solvent molecules.
	return molecules, molecule_graphs, solvent_molecule_indices, symmetry_operations

# ---------------------------------------------------------------------------------------------------------------------------

def process_crystal(original_crystal, name='crystal', take_shortest_distance=False, make_molecule_method='component_assembly_approach', no_of_cpus=1, bonds_to_ignore=None, logger=None):
	"""
	This method is designed to to remove aliphatic sidechains from your molecules in the crystal file.

	Parameters
	----------
	original_crystal : ase.Atoms
		This is ase object of the original crystal object.
	name : str.
		This is the name of the crystal.
	take_shortest_distance : bool
		If true, take the shortest distance that was found if a distance shorter than the maximum distance between elements if found. If false, give an error if the distance between two atoms was greater than the max distance expected for those two elements to be bonded to each other. Default: False
	make_molecule_method : str.
		This is the name of the method you want to use to create the molecule. See https://github.com/geoffreyweal/ECCP for more information. Default: 'component_assembly_approach'. 
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
	bonds_to_ignore : list
		This is a list of the atom pairs (given as indices) to ignore bonding between. This is given as a list of index pairs. Default: None
	logger : logging/None
		This contains the object for reporting information and warning messages to.
	
	Returns
	-------
	molecules : list of ase.Atoms objects
		These are all the molecules obtainde from the crystal.
	molecule_graphs : list of networkx.Graph objects
		These are the molecule graphs for the associated molecules in the crystal.
	"""

	# First, openning remarks
	no_of_char_in_divides = 57
	divide_string = '.'+'-'*no_of_char_in_divides+'.'
	if logger is not None:
		logger.info('Number of atoms in '+str(name)+': '+str(len(original_crystal)))

	# Second, obtain the crystal graph
	#         * Also, update the crystal object so it doesnt contain bonding/neighbours information
	crystal, crystal_graph = obtain_graph(original_crystal,mic=True,to_print=False,name=name,no_of_cpus=no_of_cpus,bonds_to_ignore=bonds_to_ignore)

	# Third, determine each of the indicidual molecules in the crystal.
	if logger is not None:
		logger.info('Determining separate molecules')
	molecules_indices = list(connected_components(crystal_graph))
	if logger is not None:
		logger.info(divide_string)
		logger.info('No of all individual molecules: '+str(len(molecules_indices)))

	# Fourth, order the names of each molecule by the lowest index of each set in molecules_indices.
	# This helps with retaining molecule name in dynamic systems. 
	molecules_indices.sort(key=lambda x: min(x))

	#raise Exception('Make sure that molecules are numbered based in order of the lowest index atom of the molecule in the crystal. Also make sure atoms in each molecule are ordered like that in the crystal.')
	#import pdb; pdb.set_trace()

	# Fifth, create each molecule in the crystal in human-friendly form.
	molecules = {}; molecule_graphs = {}
	for index in range(len(molecules_indices)):

		# 5.1: Write to the logger if logger is present.
		if logger is not None:
			logger.info('Processing molecule '+str(index+1))

		# 5.2: Get the indices of atoms in the crystal associated with the index molecule.
		molecule_indices = list(molecules_indices[index])

		# 5.3: Make the molecule of interest, and its graph
		molecule, molecule_graph = make_molecule(molecule_indices, crystal_graph, crystal, take_shortest_distance=take_shortest_distance, make_molecule_method=make_molecule_method, logger=logger)

		# 5.4: Append molecule and its graphs to lists.
		molecules[index+1] = molecule
		molecule_graphs[index+1] = molecule_graph

	# Sixth, return data.
	return molecules, molecule_graphs

# -----------------------------------------------------------------------------------------------------------------------------
