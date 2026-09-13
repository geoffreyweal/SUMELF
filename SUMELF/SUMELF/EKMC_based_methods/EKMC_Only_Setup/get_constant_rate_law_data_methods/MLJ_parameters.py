"""
get_rate_constant_data.py, Geoffrey Weal, 23/3/22

This script is designed to obtain all the values that are used in obtaining rate constants that do not change between moving from unit cell to unit cell.
"""

from math import exp, pi

from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_rate_law_data_methods.get_probability_vibrostate_on_electstate_occupied import get_probability_vibrostate_on_electstate_occupied
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_rate_law_data_methods.get_franck_condon_overlap import get_franck_condon_overlap

# Constants
kB = 8.617333262145 * (10.0 ** -5.0) # eV K-1
h_bar = 6.582119569 * (10.0 ** -16.0) # eV s

def get_constant_MLJ_parameters(molecule1, molecule2, displacement_vector, coupling_energy, huang_rhys_factor, uu_max, vv_max, WW, classical_reorganisation_energy, temperature):

	raise Exception('Check')

	# Second, get the M constant.
	M_constant = get_M_constant(coupling_energy,classical_reorganisation_energy,temperature)

	# Third, get the X constant.
	X_constant = get_X_constant(classical_reorganisation_energy,temperature)

	#Fourth, get the probabilities for the likelihood of being in the ith vibrational state in a molecule.
	vib_state_occupation_probs = get_probability_vibrostate_on_electstate_occupied(uu_max,WW,temperature)

	# Fifth, get all the constant values in the double sum
	non_changing_across_lattice_data_for_uv_inputs = {}
	for uu in range(0,uu_max+1,1):
		vib_state_uu_occupation_prob = vib_state_occupation_probs[uu]
		for vv in range(0,vv_max+1,1):
			N_constant = get_N_constant(uu,vv,vib_state_uu_occupation_prob,huang_rhys_factor)
			Y_constant = get_Y_constant(uu,vv,WW,classical_reorganisation_energy,temperature)
			Z_constant = get_Z_constant(uu,vv,WW,classical_reorganisation_energy,temperature)
			non_changing_across_lattice_data_for_uv_inputs[(uu,vv)] = (N_constant,Y_constant,Z_constant)

	return (M_constant, X_constant, non_changing_across_lattice_data_for_uv_inputs)

def get_M_constant(coupling_energy,classical_reorganisation_energy,temperature):
	"""
	This method is used to obtain the M constant

	Parameters
	----------
	coupling_energy : float
		This is the total coupling energy of the dimer, given as V12. This value is given in eV.
	classical_reorganisation_energy : float
		This is the classical component of the amount of energy that is required for the molecule in the excited state to change geometry from the ground state geometry to the excited state geometry. Given in eV.
	temperature : float
		This is the temperature of the crystal. Given in K. 

	Returns
	-------
	constant_M : float
		The M constant.
	"""
	constant_M = ((abs(coupling_energy)**2.0)/h_bar) * ((pi/(classical_reorganisation_energy*kB*temperature))**0.5)
	return constant_M

def get_N_constant(uu,vv,vib_state_uu_occupation_prob,huang_rhys_factor):
	"""
	This method is used to obtain the X constant

	Parameters
	----------
	uu : int
		This is the vibrational mode to examine in molecule 1
	vv : int
		This is the vibrational mode to examine in molecule 2
	vib_state_uu_occupation_prob : float
		Add here.
	huang_rhys_factor : float
		https://second.wiki/wiki/huang-rhys-faktor 
	Returns
	-------
	N_constant : float
		The N constant.
	"""
	franck_condon_overlap = get_franck_condon_overlap(huang_rhys_factor, uu, vv)
	N_constant = vib_state_uu_occupation_prob * (abs(franck_condon_overlap) ** 2.0)
	return N_constant

def get_X_constant(classical_reorganisation_energy,temperature):
	"""
	This method is used to obtain the X constant

	Parameters
	----------
	classical_reorganisation_energy : float
		This is the classical component of the amount of energy that is required for the molecule in the excited state to change geometry from the ground state geometry to the excited state geometry. Given in eV.
	temperature : float
		This is the temperature of the crystal. Given in K. 

	Returns
	-------
	X_constant : float
		The X constant.
	"""
	X_constant = 1/(4.0*classical_reorganisation_energy*kB*temperature)
	return X_constant

def get_Y_constant(uu,vv,WW,classical_reorganisation_energy,temperature):
	"""
	This method is used to obtain the Y constant

	Parameters
	----------
	uu : int
		This is the vibrational mode to examine in molecule 1
	vv : int
		This is the vibrational mode to examine in molecule 2
	WW : float
		This is the wang value without the hbar value. explain this better. Given in eV.
	classical_reorganisation_energy : float
		This is the classical component of the amount of energy that is required for the molecule in the excited state to change geometry from the ground state geometry to the excited state geometry. Given in eV.
	temperature : float
		This is the temperature of the crystal. Given in K. 

	Returns
	-------
	Y_constant : float
		The Y constant.
	"""
	Y_constant = (classical_reorganisation_energy + (uu-vv)*WW)/(2.0*classical_reorganisation_energy*kB*temperature)
	return Y_constant

def get_Z_constant(uu,vv,WW,classical_reorganisation_energy,temperature):
	"""
	This method is used to obtain the Z constant.

	Parameters
	----------
	uu : int
		This is the vibrational mode to examine in molecule 1
	vv : int
		This is the vibrational mode to examine in molecule 2
	WW : float
		This is the wang value without the hbar value. explain this better. Given in eV.
	classical_reorganisation_energy : float
		This is the classical component of the amount of energy that is required for the molecule in the excited state to change geometry from the ground state geometry to the excited state geometry. Given in eV.
	temperature : float
		This is the temperature of the crystal. Given in K. 

	Returns
	-------
	Z_constant : float
		The Z constant.
	"""
	Z_constant = ((classical_reorganisation_energy + (uu-vv)*WW)**2.0)/(4.0*classical_reorganisation_energy*kB*temperature)
	return Z_constant




