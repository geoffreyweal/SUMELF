"""
get_ATC_coupling_data.py, Geoffrey Weal, 5/12/23

This script will provide the settings required for obtaining atomic transition charges. 
"""

def get_ATC_coupling_data(original_coupling_data):
	"""
	This method will provide the settings required for obtaining atomic transition charges. 

	Parameters
	----------
	original_coupling_data : dict.
		This dictionary includes the information required for obtaining atomic transition charges.

	Returns
	------- 
	provided_coupling_data : dict.
		This dictionary includes the information required for obtaining atomic transition charges.
	"""

	# First, check that the dictionary being given for original_coupling_data provides the information for performing atomic transition charge calculations. 
	if not original_coupling_data['model'] == 'ATC':
		raise Exception('Something is weird, this is not an ATC dictionary? original_coupling_data = '+str(original_coupling_data))

	# Second, obtain the setting for the relative electric permittivity.
	if 'relative_permittivity' in original_coupling_data:
		relative_permittivity = original_coupling_data['relative_permittivity']
	else:
		relative_permittivity = 1.0

	# Third, set up the provided_coupling_data dictionary.
	provided_coupling_data = {'relative_permittivity': relative_permittivity}

	# Fourth, return the provided_coupling_data dictionary.
	return provided_coupling_data




