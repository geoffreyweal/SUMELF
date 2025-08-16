"""
check_ECCP_Information_details.py, Geoffrey Weal, 27/3/24

This method will check if the settings that can not change after running Run_ECCP.py for the first time have not changed when rerunning the ECCP_Information data.
"""

def check_ECCP_Information_details(eccp_information, make_dimer_method, environment_settings):
	"""
	This method will check if the settings that can not change after running Run_ECCP.py for the first time have not changed when rerunning the ECCP_Information data.

	Parameters
	----------
	path_to_ECCP_Information_folder_for_crystal : str.
		This is the path to the ECCP_Information data. 
	make_dimer_method : dict.
		This contains the make_dimer setting from the Run_ECCP.py file
	environment_settings : dict.
		This contains the environment_settings setting from the Run_ECCP.py file
	"""

	# First, obtain make_dimer_method and environment_settings from eccp_information
	make_dimer_method_from_eccp_info    = eccp_information['make_dimer_method']
	environment_settings_from_eccp_info = eccp_information['environment_settings']

	# Second, compare make_dimer_method from the Run_ECCP.py file to that from the ECCP_Information.txt file
	if not make_dimer_method == make_dimer_method_from_eccp_info:
		to_string  = 'Error: The make_dimer_method are different between Run_ECCP.py and ECCP_Information.txt.\n'
		to_string += 'make_dimer_method from Run_ECCP.py: '+str(make_dimer_method)+'\n'
		to_string += 'make_dimer_method from ECCP_Information.txt: '+str(make_dimer_method_from_eccp_info)+'\n'
		to_string += 'Check this'
		raise Exception(to_string)

	# Third, compare environment_settings from the Run_ECCP.py file to that from the ECCP_Information.txt file
	if not environment_settings == environment_settings_from_eccp_info:
		to_string  = 'Error: The environment_settings are different between Run_ECCP.py and ECCP_Information.txt.\n'
		to_string += 'environment_settings from Run_ECCP.py: '+str(environment_settings)+'\n'
		to_string += 'environment_settings from ECCP_Information.txt: '+str(environment_settings_from_eccp_info)+'\n'
		to_string += 'Check this'
		raise Exception(to_string)