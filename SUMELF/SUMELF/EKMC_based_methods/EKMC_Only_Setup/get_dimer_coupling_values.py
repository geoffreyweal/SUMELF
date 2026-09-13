"""
get_rate_law_data.py, Geoffrey Weal, 20/3/22

get_rate_law_data is a method that will give all the neighbours 
"""
import numpy as np
from copy import deepcopy

from tqdm import tqdm
from multiprocessing import Pool

from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_dimer_coupling_values_methods.get_dimers                     import get_dimers
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_dimer_coupling_values_methods.assign_coupling_value_to_dimer import assign_coupling_value_to_dimer

def get_dimer_coupling_values(kinetic_model, molecules, crystal_cell_lattice, coupling_data, kinetics_details, short_range_rCut, long_range_rCut, rCut_mol_dist_description, no_of_cpus_for_setup=1):
	"""
	This method will determine how molecules in the crystral neighbour each other over an extended lattice.

	This method will also gather other information, such as distances, orbital overlaps, coulomb energies, etc.

	Parameters
	----------
	kinetic_model : str.
		This is the kinetic model the user would like to use in their kinetic Monte Carlo algorithm upon the crystal.
	molecules : list of ase.Atoms
		This is a list of all the molecules in your crystal
	crystal_cell_lattice : ase.Cell
		This is the unit cell matrix
	coupling_data : tuple
		This tuple contains information about the couplings between the short and long range dimers in the crystal. 
		The information contained in this tuple are:
		short_range_couplings_model_name : str.
			This is the name of the model used for obtaining coupling between short-range dimers (distance less than or equal to rCut_neighourhood).
		short_range_coupling_data : dict.
			This is all the coupling information involving the short-range dimers (distance less than or equal to rCut_neighourhood).
		long_range_couplings_model_name : str.
			This is the name of the model used for obtaining coupling between long-range dimers (distance greater than than rCut_neighourhood).
		long_range_coupling_data : dict.
			This is all the coupling information involving the long-range dimers (distance greater than than rCut_neighourhood).
	kinetics_details : dict.
		These are the details required for performing the kinetic Monte Carlo algorithm.
	short_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the short-range coupling model.
	long_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the large-range coupling model.
	rCut_mol_dist_description : str.
		This is the description given to how the distance between two molecules is obtained. 
	no_of_cpus_for_setup : int.
		This is the number of cpus used to setup the EKMC simulations.
		
	Returns
	------- 
	all_local_neighbourhoods : dict.
		This dictionary contains all the information about the neighbourhoods that surrounded each molecule in your crystal.
	"""

	# First, print that the rates will be obtained from the energy data.
	print('-------------')
	print(('Getting rate constant values between molecules in surrounding unit cells.').upper())

	# Fifth, determine all the molecules that are within neighbourhood rCut (long_range_rCut) of each molecule in the central unit cell.
	# Also, get the rate constant values for each molecule that do not change with disorder. 
	if no_of_cpus_for_setup == 1:
		print('Obtaining rate constant information between molecules via a single CPU')
		all_dimer_coupling_values = get_dimers_coupling_values_using_single_cpu(molecules, crystal_cell_lattice, kinetic_model, kinetics_details, coupling_data, short_range_rCut, long_range_rCut, rCut_mol_dist_description)
	else:
		print('Obtaining rate constant information between molecules via multiple CPUs (cpu='+str(no_of_cpus_for_setup)+')')
		all_dimer_coupling_values = get_dimers_coupling_values_using_multi_cpu (molecules, crystal_cell_lattice, kinetic_model, kinetics_details, coupling_data, short_range_rCut, long_range_rCut, rCut_mol_dist_description, no_of_cpus_for_setup)
	
	# Seventh, return all_local_neighbourhoods
	print('-------------')
	return all_dimer_coupling_values
	
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_dimers_coupling_values_using_single_cpu(molecules, crystal_cell_lattice, kinetic_model, kinetics_details, coupling_data, short_range_rCut, long_range_rCut, rCut_mol_dist_description):
	"""
	This method will cycle through the molecules in the crystal and obtain the kinetic information needed to obtain rate constants for a EKMC simulation.

	Parameters
	----------
	molecules : list of ase.Atoms
		This is the list of molecules in the crystal to cycle through.
	crystal_cell_lattice : dict of int: numpy.arrays
		This is the 3x3 crystal lattice matrix.
	kinetic_model : str.
		This is the kinetic model the user would like to use in their kinetic Monte Carlo algorithm upon the crystal.
	kinetics_details : dict.
		These are the details required for performing the kinetic Monte Carlo algorithm.
	coupling_data : tuple
		This tuple contains information about the couplings between the short and long range dimers in the crystal. 
		The information contained in this tuple are:
		short_range_couplings_model_name : str.
			This is the name of the model used for obtaining coupling between short-range dimers (distance less than or equal to rCut_neighourhood).
		short_range_coupling_data : dict.
			This is all the coupling information involving the short-range dimers (distance less than or equal to rCut_neighourhood).
		long_range_couplings_model_name : str.
			This is the name of the model used for obtaining coupling between long-range dimers (distance greater than than rCut_neighourhood).
		long_range_coupling_data : dict.
			This is all the coupling information involving the long-range dimers (distance greater than than rCut_neighourhood).
	short_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the short-range coupling model.
	long_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the large-range coupling model.
	rCut_mol_dist_description : str.
		This is the description given to how the distance between two molecules is obtained. 
	"""

	# First, set up the progress bar
	pbar = tqdm(unit=' dimers assessed')

	# Second, set up the all_local_neighbourhoods dictionary.
	all_coupling_values = {}

	# Third, obtain all the kinetic information needed to obtain rate constants between molecule in a crystal for a EKMC simulation.
	for dimers in get_dimers(molecules, crystal_cell_lattice, long_range_rCut): 

		# 3.1: Give a description of what dimers are being looked at.
		toString = ['['+str(mol_1_name)+', '+str(mol_2_name)+', '+str(cell_point)+']' for mol_1_name, mol_2_name, _, _, cell_point in dimers]
		pbar.set_description('Getting coupling for (Mol1, Mol2, Rel Cell): '+', '.join(toString))

		# 3.2: Obtain the rate constant information data between molecule 1 and molecule 2 across various unit cells
		mol_1_name, mol_2_name, molecule1, molecule2, cell_point = dimers[0]
		coupling_energy = assign_coupling_value_to_dimer((mol_1_name, mol_2_name, molecule1, molecule2, cell_point, crystal_cell_lattice, short_range_rCut, long_range_rCut, rCut_mol_dist_description, coupling_data))

		# 3.3: Record coupling_energy in all_coupling_values for all dimers in dimers. These dimers will all have the same coupling value. 
		if not (coupling_energy == 0.0):
			for mol_1_name, mol_2_name, _, _, cell_point in dimers:
				all_coupling_values.setdefault(mol_1_name,{}).setdefault(mol_2_name,{})[cell_point] = coupling_energy

		# 3.4: Update the progress bar.
		pbar.update(len(dimers))

	# Fourth, close the progress bar. 
	pbar.close()

	# Fifth, return the all_local_neighbourhoods dictionary
	return all_coupling_values

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_dimers_coupling_values_using_multi_cpu(molecules, crystal_cell_lattice, kinetic_model, kinetics_details, coupling_data, short_range_rCut, long_range_rCut, rCut_mol_dist_description, no_of_cpus_for_setup):
	"""
	This method will cycle through the molecules in the crystal and obtain the kinetic information needed to obtain rate constants for a EKMC simulation.

	Parameters
	----------
	molecules : list of ase.Atoms
		This is the list of molecules in the crystal to cycle through.
	crystal_cell_lattice : dict of int: numpy.arrays
		This is the 3x3 crystal lattice matrix.
	kinetic_model : str.
		This is the kinetic model the user would like to use in their kinetic Monte Carlo algorithm upon the crystal.
	kinetics_details : dict.
		These are the details required for performing the kinetic Monte Carlo algorithm.
	coupling_data : tuple
		This tuple contains information about the couplings between the short and long range dimers in the crystal. 
		The information contained in this tuple are:
		short_range_couplings_model_name : str.
			This is the name of the model used for obtaining coupling between short-range dimers (distance less than or equal to rCut_neighourhood).
		short_range_coupling_data : dict.
			This is all the coupling information involving the short-range dimers (distance less than or equal to rCut_neighourhood).
		long_range_couplings_model_name : str.
			This is the name of the model used for obtaining coupling between long-range dimers (distance greater than than rCut_neighourhood).
		long_range_coupling_data : dict.
			This is all the coupling information involving the long-range dimers (distance greater than than rCut_neighourhood).
	short_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the short-range coupling model.
	long_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the large-range coupling model.
	rCut_mol_dist_description : str.
		This is the description given to how the distance between two molecules is obtained. 
	no_of_cpus_for_setup : int.
		This is the number of cpus used to setup the EKMC simulations.
	"""

	# First, set up the progress bar
	pbar = tqdm(unit=' dimers assessed')

	# Second, set up the all_coupling_values dictionary.
	all_coupling_values = {}

	# Third, obtain the generator for obtaining input values 
	input_values_generator = input_values(molecules, crystal_cell_lattice, short_range_rCut, long_range_rCut, rCut_mol_dist_description, coupling_data)

	# Fourth, obtain all the kinetic information needed to obtain rate constants between molecule in a crystal for a EKMC simulation.
	with Pool(processes=no_of_cpus_for_setup) as pool:
		for dimers, coupling_energy in pool.imap_unordered(get_dimer_coupling_values_multi, input_values_generator):

			# 4.1: Give a description of what dimers are being looked at.
			toString = ['['+str(mol_1_name)+', '+str(mol_2_name)+', '+str(cell_point)+']' for mol_1_name, mol_2_name, _, _, cell_point in dimers]
			pbar.set_description('Obtained coupling for (Mol1, Mol2, Rel Cell): '+', '.join(toString))

			# 4.2: Record coupling_energy in all_coupling_values for all dimers in dimers. These dimers will all have the same coupling value. 
			if not (coupling_energy == 0.0):
				for mol_1_name, mol_2_name, _, _, cell_point in dimers:
					all_coupling_values.setdefault(mol_1_name,{}).setdefault(mol_2_name,{})[cell_point] = coupling_energy

			# 4.3: update the progress bar.
			pbar.update(len(dimers))

	# Fifth, return the all_coupling_values dictionary
	return all_coupling_values

def input_values(molecules, crystal_cell_lattice, short_range_rCut, long_range_rCut, rCut_mol_dist_description, coupling_data):
	for dimers in get_dimers(molecules, crystal_cell_lattice, long_range_rCut):
		yield (dimers, crystal_cell_lattice, short_range_rCut, long_range_rCut, rCut_mol_dist_description, coupling_data)

def get_dimer_coupling_values_multi(input_information):
	"""
	Thid method is designed to run assign_coupling_value_to_dimer in the multiprocessing component of the get_dimers_coupling_values_using_multi_cpu method. 
	"""

	# First, extract information from input_information
	dimers, crystal_cell_lattice, short_range_rCut, long_range_rCut, rCut_mol_dist_description, coupling_data = input_information

	# Second, obtain the rate constant information data between molecule 1 and molecule 2 across various unit cells
	mol_1_name, mol_2_name, molecule1, molecule2, cell_point = dimers[0]
	coupling_energy = assign_coupling_value_to_dimer((mol_1_name, mol_2_name, molecule1, molecule2, cell_point, crystal_cell_lattice, short_range_rCut, long_range_rCut, rCut_mol_dist_description, coupling_data))

	# Third, return dimers and coupling_energy
	return (dimers, coupling_energy)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


