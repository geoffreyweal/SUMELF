"""
EKMC_Only_Setup.py, Geoffrey Weal, 12/8/22

This script is designed to set up the files needed for runningt the excitonic kinetic Monte Carlo algorithm. 

"""
import os
from copy import deepcopy
from numpy import array_equal

from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.initial_setup_data_methods                     import save_initial_setup_data, read_initial_setup_data
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_molecules                                  import get_molecules
from SUMELF                                                                         import obtain_graph

from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.coupling_methods.get_ATC_coupling_data         import get_ATC_coupling_data
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.coupling_methods.get_EET_coupling_data         import get_EET_coupling_data
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.coupling_methods.add_zero_charges_to_molecules import add_zero_charges_to_molecules
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.coupling_methods.add_charges_to_molecules      import add_charges_to_molecules
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_RE_and_bandgap_data                        import get_RE_and_bandgap_data

from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_dimer_coupling_values                      import get_dimer_coupling_values
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_constant_rate_law_data                     import get_constant_rate_law_data

from SUMELF                                                                         import remove_folder, make_folder

def EKMC_Only_Setup(molecules_path, functional_and_basis_set, kinetic_model, short_range_couplings, long_range_couplings, short_range_rCut, long_range_rCut, rCut_mol_dist_description, reorganisation_and_bandgap_energy_details, kinetics_details, include_solvents=True, path_to_EKMC_simulations='', path_to_initial_EKMC_setup_files=None, no_of_cpus_for_setup=1):
	"""
	This program is designed to simulate the movement of an exciton through a OPV crystal system.

	Parameters
	----------
	molecules_path : str.
		The path to the molecules that make up the crystal, as obtained from the ECCP program. 
	functional_and_basis_set : str.
		This is the folder name for the functional and basis set used in calculations. This name is given from the ECCP program. 
	kinetic_model : str.
		This is the type of kinetic model you would like to use.
	short_range_couplings : dict. 
		This describes the coupling between molecules in the crystal at distances less than and including short_range_rCut. This dictionary includes the model ('model': None if you dont want to include short range coupling) and the files which will be used to obtain the coupling energies.
	long_range_couplings : dict.
		This describes the coupling between molecules in the crystal at distances between short_range_rCut and long_range_rCut. This dictionary includes the model ('model': None if you dont want to include long range coupling) and the files which will be used to obtain the coupling energies.
	short_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the short-range coupling model.
	long_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the large-range coupling model.
	rCut_mol_dist_description : str.
		This is the description given to how the distance between two molecules is obtained. 
	reorganisation_and_bandgap_energy_details : dict
		This dictionary contains all the details about where to get the data to calculate the reorganisation energy for an exciton moving from one molecule to another, and bandgap energy for a molecule. 
	kinetics_details : dict.
		This contains all the information required to obtain the rate constants (except for electronic coupling and ATC charge data).
	path_to_EKMC_simulations : str.
		This is the path to the ekmc file that contains all the information about the ekmc simulation for this crystal.
	path_to_initial_EKMC_setup_files : str.
		This is the path to the EKMC setup files that contain spatial coupling data. This info can be useful to have if you are performing EKMC simulations with repeated sims with different energetic disorders and reorganisation energies for example. 
	no_of_cpus_for_setup : int.
		This is the number of cpus used to setup the EKMC simulations.
	"""

	# First, perform the initial component of the setup.
	if (path_to_initial_EKMC_setup_files is None) or (path_to_initial_EKMC_setup_files == ''):
		molecule_names, molecules, crystal_cell_lattice, all_coupling_values = initial_setup(molecules_path, functional_and_basis_set, kinetic_model, short_range_couplings, long_range_couplings, short_range_rCut, long_range_rCut, rCut_mol_dist_description, kinetics_details, include_solvents=include_solvents, no_of_cpus_for_setup=no_of_cpus_for_setup)
	elif not os.path.exists(path_to_initial_EKMC_setup_files):
		molecule_names, molecules, crystal_cell_lattice, all_coupling_values = initial_setup(molecules_path, functional_and_basis_set, kinetic_model, short_range_couplings, long_range_couplings, short_range_rCut, long_range_rCut, rCut_mol_dist_description, kinetics_details, include_solvents=include_solvents, no_of_cpus_for_setup=no_of_cpus_for_setup)
		save_initial_setup_data(path_to_initial_EKMC_setup_files, molecule_names, molecules, crystal_cell_lattice, all_coupling_values)
	else:
		print('FOUND INITIAL SETUP DATA IN '+str(path_to_initial_EKMC_setup_files)+'. WILL USE THIS FOR GETTING MOLECULES, CRYSTAL LATTICE DATA, KINETIC DATA, AND COUPLING VALUES.')
		molecule_names, molecules, crystal_cell_lattice, all_coupling_values = read_initial_setup_data(path_to_initial_EKMC_setup_files)

	# Second, extract the non0change componenets of kinetics_details
	non_changing_lattice_kinetics_details = get_non_changing_lattice_kinetics_details(kinetic_model, kinetics_details)

	# Third, obtain the constant rate law data for the rate law you want to use. 
	constant_rate_data = get_constant_rate_law_data(kinetic_model, non_changing_lattice_kinetics_details)

	# Fourth, get the ground and excited structure energy gaps for molecules in the crystal.
	crystal_name = os.path.basename(molecules_path)
	molecule_bandgap_energy_data, dimer_reorganisation_energy_data, conformationally_equivalent_molecules = get_RE_and_bandgap_data(reorganisation_and_bandgap_energy_details, crystal_name, functional_and_basis_set, molecules_path, molecule_names)

	# Fifth, save this data to disk.
	save_KMC_data_to_disk(path_to_EKMC_simulations, molecule_names, molecules, crystal_cell_lattice, kinetic_model, non_changing_lattice_kinetics_details, molecule_bandgap_energy_data, dimer_reorganisation_energy_data, conformationally_equivalent_molecules, constant_rate_data, all_coupling_values)

# ==================================================================================================================================================================================================

def initial_setup(molecules_path, functional_and_basis_set, kinetic_model, short_range_couplings, long_range_couplings, short_range_rCut, long_range_rCut, rCut_mol_dist_description, kinetics_details, include_solvents=True, no_of_cpus_for_setup=1):
	"""
	This program is designed to simulate the movement of an exciton through a OPV crystal system. This is the data that is saved to file. 

	Parameters
	----------
	molecules_path : str.
		The path to the molecules that make up the crystal, as obtained from the ECCP program. 
	functional_and_basis_set : str.
		This is the folder name for the functional and basis set used in calculations. This name is given from the ECCP program. 
	kinetic_model : str.
		This is the type of kinetic model you would like to use.
	short_range_couplings : dict. 
		This describes the coupling between molecules in the crystal at distances less than and including short_range_rCut. This dictionary includes the model ('model': None if you dont want to include short range coupling) and the files which will be used to obtain the coupling energies.
	long_range_couplings : dict.
		This describes the coupling between molecules in the crystal at distances between short_range_rCut and long_range_rCut. This dictionary includes the model ('model': None if you dont want to include long range coupling) and the files which will be used to obtain the coupling energies.
	short_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the short-range coupling model.
	long_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the large-range coupling model.
	rCut_mol_dist_description : str.
		This is the description given to how the distance between two molecules is obtained. 
	kinetics_details : dict.
		This contains all the information required to obtain the rate constants (except for electronic coupling and ATC charge data)
	include_solvents : bool.
		This boolean indicates if you want to include solvents in your KMC model. Default: True.
	no_of_cpus_for_setup : int.
		This is the number of cpus used to setup the EKMC simulations.
	"""

	# First, get the absolute path to the folders that we will be looking at.
	crystal_name    = os.path.basename(molecules_path)
	molecules_path  = os.path.abspath(molecules_path)

	# Second, get the molecules that make up the crystal, and the crystal lattice.
	molecule_names, molecules = get_molecules(molecules_path, include_solvents=include_solvents)
	lowest_number_molecule_name = min(molecules.keys())
	if not all([array_equal(molecules[molecule_name].get_cell()[:], molecules[lowest_number_molecule_name].get_cell()[:]) for molecule_name in molecules.keys()]):
		raise Exception('Error. Not all the unit cells of the molecules in your crystal are the same. Check the unit cells of your molecules out.')
	crystal_cell_lattice = molecules[lowest_number_molecule_name].get_cell()

	# Third, get the graphs of molecules
	molecule_graphs = {}
	for molecule_name in molecules.keys():
		#import pdb; pdb.set_trace()
		#updated_molecule, molecule_graph = obtain_graph(molecules[molecule_name],mic=False,name='molecule_'+str(molecule_name))
		molecule_graph                   = obtain_graph(molecules[molecule_name],mic=False,name='molecule_'+str(molecule_name))
		molecule_graphs[molecule_name]   = deepcopy(molecule_graph)
		#molecules[molecule_name]         = updated_molecule.copy()

	# Fourth, check that rCut_mol_dist_description is one of the following, and get the distance between the molecules in the dimer.
	if rCut_mol_dist_description not in ['centre_of_mass', 'centre_of_molecule', 'nearest_atoms_method']:
		toString  = 'Error: Your input for "rCut_mol_dist_description" must be either:\n'
		toString += '\t* "centre_of_mass":       Take the dimer distance from the centre of masses of each molecule.\n'
		toString += '\t* "centre_of_molecule":   Take the dimer distance from the centre of molecules of each molecule (excluding hydrogen atoms).\n'
		toString += '\t* "nearest_atoms_method": Take the dimer distance as the closest distance between atoms in each molecule (excluding hydrogen atoms).\n'
		toString += 'Check this in your input python script and rerun.'
		raise Exception(toString)

	# Fifth, obtain the coupling data for short range coupling. 
	# The model that you can use must be either ATC or EET. 
	if not isinstance(short_range_rCut,float):
		raise Exception('Error: you need to set short_range_rCut to some float value. short_range_rCut = '+str(short_range_rCut))
	if short_range_couplings['model'] == 'ATC':
		short_range_coupling_data = get_ATC_coupling_data(short_range_couplings)
	elif short_range_couplings['model'] == 'EET':
		short_range_coupling_data = get_EET_coupling_data(short_range_couplings, crystal_name, functional_and_basis_set, molecules_path)
	else:
		toString  = 'Error: The short-range coupling model method must be either EET or ATC. '+'\n'
		toString += "short_range_couplings['model'] = "+str(short_range_couplings['model'])+'\n'
		toString += 'Check this.'+'\n'
		raise Exception(toString)

	# Sixth, obtain the coupling data for long range coupling.
	# The model that you can use must be either ATC or None. 
	if (long_range_couplings['model'] is None) or (long_range_couplings['model'].lower() == 'none'):
		# What to do if you are not using any long range coupling description
		long_range_couplings['model'] = None
		long_range_coupling_data = {}
		if long_range_rCut >= short_range_rCut:
			short_range_rCut = long_range_rCut
		long_range_rCut = short_range_rCut
	elif long_range_couplings['model'] == 'ATC':
		# What to do if you are using the ATC long range coupling description
		if not isinstance(long_range_rCut,float):
			raise Exception('Error: you need to set long_range_rCut to some float value. long_range_rCut = '+str(long_range_rCut))
		elif long_range_rCut < short_range_rCut:
			raise Exception('Error: long_range_rCut is smaller than short_range_rCut. short_range_rCut = '+str(short_range_rCut)+', long_range_rCut = '+str(long_range_rCut))
		long_range_coupling_data = get_ATC_coupling_data(long_range_couplings)
	else:
		# You need to do one of the two options above, raise exception if you have not. 
		toString  = 'Error: The long-range coupling model method must be either ATC or None. '+'\n'
		toString += "long_range_couplings['model'] = "+str(long_range_couplings['model'])+'\n'
		toString += 'Check this.'+'\n'
		raise Exception(toString)

	# Seventh, if short_range_couplings['model'] == long_range_couplings['model'], set short_range_rCut to long_range_rCut so that their are no potential future issues with setting this value
	if short_range_couplings['model'] == long_range_couplings['model']:
		if short_range_couplings['model'] == 'ATC':
			if not short_range_couplings['path_to_ATC_folder'] == long_range_couplings['path_to_ATC_folder']:
				toString  = 'Error: Your ATC inputs are not the same'+'\n'
				toString += "short_range_couplings['path_to_ATC_folder'] = "+str(short_range_couplings['path_to_ATC_folder'])+'\n'
				toString += "long_range_couplings['path_to_ATC_folder']  = "+str(long_range_couplings['path_to_ATC_folder'])+'\n'
				toString += 'Check this.'
				raise Exception(toString)

	# Eighth, collect all the short and long range coupling values together for use later in this setup program.
	coupling_data = (short_range_couplings['model'], short_range_coupling_data, long_range_couplings['model'], long_range_coupling_data)

	# Ninth, if the ATC model is used for either or both the short or long range models, add atomic transition charges to the molecules in the crystal.
	# Here, charges will be added to molecules from chg files in ATC folder. 
	# The charges added to the molecules are the atomic transition charges. 
	if (short_range_couplings['model'] == 'ATC') or (long_range_couplings['model'] == 'ATC'):
		if short_range_couplings['model'] == 'ATC':
			ATC_folder_path              = short_range_couplings['path_to_ATC_folder']
		elif long_range_couplings['model'] == 'ATC':
			ATC_folder_path              = long_range_couplings['path_to_ATC_folder']
		add_charges_to_molecules(molecules, molecule_graphs, ATC_folder_path, molecules_path, crystal_name, functional_and_basis_set, include_solvents=include_solvents)
	else:
		add_zero_charges_to_molecules(molecules, molecule_graphs)

	# Tenth, extract the non-change components of kinetics_details
	non_changing_lattice_kinetics_details = get_non_changing_lattice_kinetics_details(kinetic_model, kinetics_details)

	# Eleventh, obtain the coupling values for each dimer in the crystal, expanding to short_range_rCut and long_range_rCut
	all_coupling_values = get_dimer_coupling_values(kinetic_model, molecules, crystal_cell_lattice, coupling_data, non_changing_lattice_kinetics_details, short_range_rCut=short_range_rCut, long_range_rCut=long_range_rCut, rCut_mol_dist_description=rCut_mol_dist_description, no_of_cpus_for_setup=no_of_cpus_for_setup)

	# Twelfth, convert molecules dictionary to list corresponding to molecule_names
	molecules_list = [molecules[int(molecule_name.replace('S',''))] for molecule_name in molecule_names]

	# Thirteenth, return data.
	return molecule_names, molecules_list, crystal_cell_lattice, all_coupling_values

# ==================================================================================================================================================================================================

def get_non_changing_lattice_kinetics_details(kinetic_model, kinetics_details):
	"""
	This method is designed to store the kinetic details that do not change during setup.

	Parameters
	----------
	kinetic_model : str.
		This is the type of kinetic model you would like to use.
	kinetics_details : dict.
		This contains all the information required to obtain the rate constants (except for electronic coupling and ATC charge data).
	
	Returns
	-------
	return a dictionary containing all the kinetic details that do not change during the simulation. 
	"""

	# First, get the names of the neighbourhood dictionary keys to obtain.
	names = []
	if   kinetic_model.lower() == 'marcus':
		names += ['temperature']
	elif kinetic_model.lower() == 'mlj':
		names += ['huang_rhys_factor','uu_max','vv_max','WW','temperature']
	else:
		raise Exception('Error: kinetic_model must be either "Marcus" or "MLJ". kinetic_model = '+str(kinetic_model))
	names += ['coupling_disorder', 'energetic_disorder']

	# Second, record the kinetic detail to a new dictionary.
	non_changing_lattice_kinetics_details = {}
	for name in names:
		non_changing_lattice_kinetics_details[name] = deepcopy(kinetics_details[name])

	# Third, return non_changing_lattice_kinetics_details
	return non_changing_lattice_kinetics_details

def save_KMC_data_to_disk(path_to_KMC_setup_data, molecule_names, molecules, crystal_cell_lattice, kinetic_model, kinetics_details, molecule_bandgap_energy_data, dimer_reorganisation_energy_data, conformationally_equivalent_molecules, constant_rate_data, all_coupling_values):
	"""
	This method is designed to save the data required to simulate the exciton kmc simulation, including the electronic details of the local behaviour of each molecule in the crystal.

	Parameters
	----------
	path_to_KMC_setup_data : str.
		This is the path to the ekmc file that contains all the information about the ekmc simulation for this crystal.
	molecule_names : list of str or ints
		These are the names of the molecules in the crystal, including if it is a solvent or not.
	molecules : list of ase.Atoms
		This is the list of molecules in the crystal
	crystal_cell_lattice : ase.Cell
		This contain information about the unit cell.
	kinetic_model : str.
		This is the type of kinetic model you would like to use.
	kinetics_details : dict.
		This contains all the information required to obtain the rate constants (except for electronic coupling and ATC charge data).
	molecule_bandgap_energy_data : dict.
		These are the bandgap energies of all unique molecules in the crystal
	dimer_reorganisation_energy_data : dict.
		These are reorganisation energies for all dimers in the crystal
	conformationally_equivalent_molecules : dict. 
		This dictionary contains informatino about which conformationally equivalent molecules are assigned to which conformatinoally unique molecules (and theirfore have the same reorganisation energy data assigned to them).
	constant_rate_data : tuple
		These are the constants in the rate equation that do not change between neighbours.
	all_coupling_values : dict.
		This dictionary contains all the information about the coupling between all the dimers involving the exciton donor in the (0,0,0) unit cell in your crystal.
	"""

	# First, make the folder to place the KMC_setup_data.ekmc file in 
	#remove_folder(path_to_KMC_setup_data)
	make_folder(path_to_KMC_setup_data)

	# Second, create the placement for the molecule names list that is sorted from lowest to highest numbered name.
	molecule_list = [(molecule_name, tuple(molecule.get_center_of_mass())) for molecule_name, molecule in sorted(zip(molecule_names, molecules), key=lambda x: int(str(x[0]).replace('S','')))]
	molecule_list = [str(molecule_name)+': '+str(molecule_com) for molecule_name, molecule_com in molecule_list]

	# Third, create the KMC_setup_data.ekmc file.
	path_to_file = path_to_KMC_setup_data+'/'+f'KMC_setup_data.ekmc'
	with open(path_to_file,'w') as KMC_setup_data:
		KMC_setup_data.write('{'+', '.join(molecule_list)+'}\n')
		KMC_setup_data.write(str(tuple([tuple(xx) for xx in crystal_cell_lattice[:]]))+'\n')
		KMC_setup_data.write(str(kinetic_model)+'\n')
		KMC_setup_data.write(str(kinetics_details)+'\n')
		KMC_setup_data.write(str(molecule_bandgap_energy_data)+'\n')
		KMC_setup_data.write(str(dimer_reorganisation_energy_data)+'\n')
		KMC_setup_data.write(str(conformationally_equivalent_molecules)+'\n')
		KMC_setup_data.write(str(constant_rate_data)+'\n')
		KMC_setup_data.write(save_all_coupling_values_data(all_coupling_values)+'\n')

def save_all_coupling_values_data(all_coupling_values):
	"""
	This method is designed to order the printing of all_coupling_values.

	Parameters
	----------
	all_coupling_values : dict.
		This dictionary contains all the information about the neighbourhoods that surrounded each molecule in your crystal.

	Returns
	-------
	to_string : str.
		This is the written string that contains information about all the coupling energies and etc between molecules in all the local neighbourhoods.
	"""

	# First, create a list to place data into
	to_string = []
	for mol1 in sorted(all_coupling_values.keys()):
		to_string_mol1_data = []
		for mol2 in sorted(all_coupling_values[mol1].keys()):
			to_string_mol1_mol2_data = []
			for cell_point, rate_constant_data in sorted(all_coupling_values[mol1][mol2].items(), key=lambda cpd: tuple([order_number(cp) for cp in cpd[0]])):
				to_string_mol1_mol2_data.append(str(cell_point)+': '+str(rate_constant_data))
			to_string_mol1_mol2_data = '{' + ', '.join(to_string_mol1_mol2_data) + '}'
			to_string_mol1_data.append(str(mol2)+': '+str(to_string_mol1_mol2_data))
		to_string_mol1_data = '{' + ', '.join(to_string_mol1_data) + '}'
		to_string.append(str(mol1)+': '+str(to_string_mol1_data))
	to_string = '{' + ', '.join(to_string) + '}'

	return to_string

def order_number(value):
	if value == 0:
		return value
	elif value < 0:
		return abs(value)*2
	elif value > 0:
		return abs(value)*2 - 1

# ==================================================================================================================================================================================================



