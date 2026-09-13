"""
assign_coupling_value_to_dimer.py, Geoffrey Weal, 25/12/22

This script is designed to return the coupling value that is appropriate to the dimer, whether the dimer is short or long range (compared to rCut_neighbour).
"""

import numpy as np
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.calculators.get_coulomb_energy import get_coulomb_energy
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_dimer_coupling_values_methods.are_molecules_within_rCut import are_molecules_within_rCut

def assign_coupling_value_to_dimer(input_information):
	"""
	This method is designed to return the coupling value that is appropriate to the dimer, whether the dimer is short or long range (compared to rCut_neighbour).

	Parameters
	----------
	molname1 : int
		This is the name of the first molecule.
	molname2 : int
		This is the name of the second molecule.
	cell_point : int
		This is the relative difference in unit cells of the second molecule in the dimer of interest relative to the first molecule (in the dimer of interest).
	crystal_cell_lattice : 3x3 np.array
		This is the matrix that descripts the unit cell for the crystal. 
	are_molecules_within_rCut : str.
		This string indicates if this dimer is within short_rCut or long_rCut
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

	Returns
	-------
	coupling_energy : float
		This is the coupling energy value between the molecules in the dimer of interest. 
	"""

	# First, separate input_information into its components, includeing the coupling_data into short-range and long-range data sets. 
	molname1, molname2, molecule1, molecule2, cell_point, crystal_cell_lattice, short_range_rCut, long_range_rCut, rCut_mol_dist_description, (short_range_couplings_model_name, short_range_coupling_data, long_range_couplings_model_name, long_range_coupling_data) = input_information

	# Second, determine if the two molecules in the dimer are within short range rCut or long range rCut of each other. 
	is_molecules_within_rCut = are_molecules_within_rCut(molecule1, molecule2, cell_point, crystal_cell_lattice, short_range_rCut, long_range_rCut, rCut_mol_dist_description)

	# Third, obtain the coupling for this dimer set based on dimer_distance vs short_range_rCut
	if is_molecules_within_rCut == 'Short Range': 

		# 4.1: Use short-range coupling method
		if short_range_couplings_model_name == 'EET':
			coupling_energy = short_range_coupling_data.get((molname1, molname2, *cell_point), None)
			if coupling_energy is None:
				raise Exception('Error. There is not EET coupling energy assigned to this dimer. Molecule1: '+str(molname1)+', Molecule2: '+str(molname2)+', cell_point:'+str(cell_point))
		elif short_range_couplings_model_name == 'ATC':
			displacement_vector = np.matmul(np.array(cell_point), crystal_cell_lattice)
			coupling_energy = get_coulomb_energy(molecule1, molecule2, displacement_vector, short_range_coupling_data['relative_permittivity'])
		else:
			toString  = 'Error: The short-range coupling model method must be either EET or ATC. '+'\n'
			toString += 'short_range_couplings_model_name = '+str(short_range_couplings_model_name)+'\n'
			toString += 'Check this.'+'\n'
			raise Exception(toString)

	elif is_molecules_within_rCut == 'Long Range': 

		# 4.2: Use long-range coupling method
		if long_range_couplings_model_name == 'ATC':
			displacement_vector = np.matmul(np.array(cell_point), crystal_cell_lattice)
			coupling_energy = get_coulomb_energy(molecule1, molecule2, displacement_vector, long_range_coupling_data['relative_permittivity'])
		elif long_range_couplings_model_name is None:
			coupling_energy = 0.0
		else:
			toString  = 'Error: The long-range coupling model method must be either ATC or None. '+'\n'
			toString += 'long_range_couplings_model_name = '+str(long_range_couplings_model_name)+'\n'
			toString += 'Check this.'+'\n'
			raise Exception(toString)
	else:

		# 4.3: The two dimers are outside of rCut, so set coupling energy to 0.0 eV
		coupling_energy = 0.0

	# Fourth, return the coupling energy between the molecules in the dimer.
	return coupling_energy
	
# ============================================================================================================================================
