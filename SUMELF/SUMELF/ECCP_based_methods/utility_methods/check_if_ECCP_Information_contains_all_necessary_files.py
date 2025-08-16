"""
does_contain_ECCP_Information_data.py, Geoffrey Weal, 25/3/24

This method is designed to check if there already exists ECCP_Information data.
"""
import os
from SUMELF.SUMELF.ECCP_based_methods.utility_methods.get_list_of_files_to_check import get_list_of_files_to_check

# Get the list of files to check for in the ECCP information folder. 
check_for_files = get_list_of_files_to_check()

def check_if_ECCP_Information_contains_all_necessary_files(path_to_ECCP_Information_folder_for_crystal):
	"""
	This method is designed to check that all the necessary files needed exist.

	Parameters
	----------
	path_to_ECCP_Information_folder_for_crystal : str.
		This is the path to the ECCP_Information data. 
	"""

	# First, initalise the list to add files that were not found in the ECCP_Information folder
	can_not_find_files = []

	# Second, go through the files in check_for_files and check if these files exist in path_to_ECCP_Information_folder_for_crystal
	for file_to_check in check_for_files:
		if not os.path.exists(path_to_ECCP_Information_folder_for_crystal+'/'+file_to_check):
			can_not_find_files.append(file_to_check)

	# Third, if there are any files that have not been found in the ECCP_Information folder, report this. 
	if len(can_not_find_files) > 0:
		print('=======================================================================')
		print('Error: Could not find the following files in your '+str(path_to_ECCP_Information_folder_for_crystal)+' folder:')
		for can_not_find_file in can_not_find_files:
			print(can_not_find_file)
		print('Check this. This program will finish without beginning')
		print('=======================================================================')
		exit()

