"""
get_constant_rate_law_data.py, Geoffrey Weal, 14/12/23

This method is designed obtain components of the rate law that are the products of all constant values for the rate law you want to use. 
"""
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_constant_rate_law_data_methods.Marcus_parameters import get_constant_Marcus_parameters
#from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_rate_law_data_methods.rate_law_methods.MLJ_parameters import get_constant_MLJ_parameters

def get_constant_rate_law_data(kinetic_model, kinetics_details):
	"""
	This method will contain the components of the rate law that are the products of all constant values for the rate law you want to use. 

	Parameters
	----------
	kinetic_model : str.
		This is the kinetic model the user would like to use in their kinetic Monte Carlo algorithm upon the crystal.
	kinetics_details : dict.
		These are the details required for performing the kinetic Monte Carlo algorithm.
		
	Returns
	------- 
	constant_rate_data : dict.
		This dictionary contains all the information about the neighbourhoods that surrounded each molecule in your crystal.
	"""

	# Fourth, get the parameters that are constant for any molecule
	if   kinetic_model.lower() == 'marcus':
		constant_rate_data = get_constant_Marcus_parameters(kinetics_details['temperature'])
	elif kinetic_model.lower() == 'mlj':
		constant_rate_data = get_MLJ_parameters(**kinetics_details) # to do.
	else:
		raise Exception('Error: You must set "kinetic_model" to either "Marcus" or "MLJ". Your kinetic_model is set to: '+str(kinetic_model)+'. Check this out.')
	
	# Seventh, return all_local_neighbourhoods
	return constant_rate_data
