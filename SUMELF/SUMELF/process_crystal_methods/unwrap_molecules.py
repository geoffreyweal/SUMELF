"""
unwrap_molecules.py, Geoffrey Weal, 2/6/24

This method is designed to unwrap the molecule, allowing the molecule to form as one connected unit.
"""
import sys
from tqdm import tqdm
import multiprocessing as mp

from SUMELF.SUMELF.process_crystal_methods.make_molecule import make_molecule

def unwrap_molecules(molecules, molecule_graphs, take_shortest_distance, make_molecule_method, no_of_cpus, logger, print_progress=True):
	"""
	This method is designed to unwrap the molecule, allowing the molecule to form as one connected unit.

	Parameters
	----------
	molecules : list of ase.Atoms
		These are the molecules that make up the crystal (before any spacegroup symmetry is considered).
	molecule_graphs : list of networkx.Graph
		These are the graph that are associated with the molecules in the molecules list. 
	take_shortest_distance : bool
		If true, take the shortest distance that was found if a distance shorter than the maximum distance between elements if found. If false, give an error if the distance between two atoms was greater than the max distance expected for those two elements to be bonded to each other. Default: False
	make_molecule_method : str.
		This is the name of the method you want to use to create the molecule. See https://github.com/geoffreyweal/ECCP for more information. Default: 'component_assembly_approach'. 
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
	logger : logging/None
		This contains the object for reporting information and warning messages to.
	print_progress : bool.
		Indicates if you want to write what the method is currently doing to the terminal. Default: True
	"""

	if no_of_cpus == 1: # If you want to perform this task with one cpu, perform task without multiprocessing.

		# First, initialise the dictionary to store unwrapped molecules in.
		unwrapped_molecules = {}

		# Second, initialise the progress bar.
		pbar = get_inputs(molecules, molecule_graphs, take_shortest_distance, make_molecule_method, logger, unwrapped_molecules)

		# Third, indicate if you want to print the progress to the terminal.
		if print_progress:
			pbar = tqdm(pbar, total=len(molecules), desc='Unwrapping molecules', unit='calc')

		# Fourth, unwrap each molecule in molecules.
		for input_data in pbar:
			unwrap_molecule_single_process(input_data)

		# Fifth, close pbar
		pbar.close()

	else: # If you want to run this method on multiple cpus. 

		# Sixth, create the manager to save the dictionary of unwrapped molecules to.
		with mp.Manager() as manager:

			# Seventh, create the dictionary to store unwrapped molecules to. 
			unwrapped_molecules = manager.dict()

			# Eighth, write warning message to the user.
			if print_progress:
				print('Unwrapping molecules (Please wait until after 100%, as the process will still be running.)', file=sys.stderr)

			# Ninth, initialise the progress bar.
			pbar = get_inputs(molecules, molecule_graphs, take_shortest_distance, make_molecule_method, logger, unwrapped_molecules)

			# Tenth, indicate if you want to print the progress to the terminal.
			if print_progress:
				pbar = tqdm(pbar, total=len(molecules), desc='Unwrapping molecules', unit='calc')

			# Eleventh, run the multiprocessing jobs.
			#process_map(unwrap_molecule_single_process,        get_inputs(molecules, molecule_graphs, take_shortest_distance, make_molecule_method, logger, unwrapped_molecules), total=len(molecules), desc='Unwrapping molecules', unit='calc', max_workers=no_of_cpus)
			pool = mp.Pool(no_of_cpus)
			#pool.map_async(unwrap_molecule_single_process, tqdm(get_inputs(molecules, molecule_graphs, take_shortest_distance, make_molecule_method, logger, unwrapped_molecules), total=len(molecules), desc='Unwrapping molecules', unit='calc'))
			pool.map_async(unwrap_molecule_single_process, pbar)
			pool.close()
			pool.join()

			# Eighth, save the dictionary from a multiprocessing Manager dictionary as a regular dictionary.
			unwrapped_molecules = dict(unwrapped_molecules)

	# Twelfth, return unwrapped_molecules
	return unwrapped_molecules

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def get_inputs(molecules, molecule_graphs, take_shortest_distance, make_molecule_method, logger, unwrapped_molecules):
	"""
	This generator will provide the input data for the make_molecule_single_process method

	Parameters
	----------
	molecules : list of ase.Atoms
		These are the molecules that make up the crystal (before any spacegroup symmetry is considered).
	molecule_graphs : list of networkx.Graph
		These are the graph that are associated with the molecules in the molecules list. 
	take_shortest_distance : bool
		If true, take the shortest distance that was found if a distance shorter than the maximum distance between elements if found. If false, give an error if the distance between two atoms was greater than the max distance expected for those two elements to be bonded to each other. Default: False
	make_molecule_method : str.
		This is the name of the method you want to use to create the molecule. See https://github.com/geoffreyweal/ECCP for more information. Default: 'component_assembly_approach'. 
	logger : logging/None
		This contains the object for reporting information and warning messages to.
	unwrapped_molecules : dict/multiprocessing.manager.dict
		This is the dictionary to store unwrapped molecules in.

	Returns
	-------
	mol_name :int
		This is the name of the molecule we want to unwrap.
	wrapped_molecule : ase.Atoms
		This is the ase.Atoms object for the molecule. 
	molecule_graph : networkx.Graph
		this isthe graph that are associated molecule given by wrapped_molecule
	take_shortest_distance : bool
		If true, take the shortest distance that was found if a distance shorter than the maximum distance between elements if found. If false, give an error if the distance between two atoms was greater than the max distance expected for those two elements to be bonded to each other. Default: False
	make_molecule_method : str.
		This is the name of the method you want to use to create the molecule. See https://github.com/geoffreyweal/ECCP for more information. Default: 'component_assembly_approach'. 
	logger : logging/None
		This contains the object for reporting information and warning messages to.
	unwrapped_molecules : dict/multiprocessing.manager.dict
		This is the dictionary to store unwrapped molecules in.
	"""

	# First, for each wrapped molecule in the molecules dictionary
	for mol_name, wrapped_molecule in sorted(molecules.items(), key=lambda x:x[0]):

		# Second, yield input data
		yield (mol_name, wrapped_molecule, molecule_graphs[mol_name], take_shortest_distance, make_molecule_method, logger, unwrapped_molecules)

def unwrap_molecule_single_process(input_data):
	"""
	This method is designed to provide the single cpu process for unwrapping a molecule.

	Parameters
	----------
	mol_name :int
		This is the name of the molecule we want to unwrap.
	wrapped_molecule : ase.Atoms
		This is the ase.Atoms object for the molecule. 
	molecule_graph : networkx.Graph
		this isthe graph that are associated molecule given by wrapped_molecule
	take_shortest_distance : bool
		If true, take the shortest distance that was found if a distance shorter than the maximum distance between elements if found. If false, give an error if the distance between two atoms was greater than the max distance expected for those two elements to be bonded to each other. Default: False
	make_molecule_method : str.
		This is the name of the method you want to use to create the molecule. See https://github.com/geoffreyweal/ECCP for more information. Default: 'component_assembly_approach'. 
	logger : logging/None
		This contains the object for reporting information and warning messages to.
	unwrapped_molecules : dict/multiprocessing.manager.dict
		This is the dictionary to store unwrapped molecules in.
	"""

	# First, obtain all the input variables from input_data.
	mol_name, wrapped_molecule, molecule_graph, take_shortest_distance, make_molecule_method, logger, unwrapped_molecules = input_data

	# Second, this list allow all atoms in the molecule to be processed in the unwrapped molecule.
	molecule_indices = list(range(len(wrapped_molecule)))

	# Third, this method will unwrap the molecule
	molecule, _ = make_molecule(molecule_indices, molecule_graph, wrapped_molecule, take_shortest_distance=take_shortest_distance, make_molecule_method=make_molecule_method, logger=logger)

	# Fourth, save the unwrapped molecule to molecules
	unwrapped_molecules[mol_name] = molecule

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

