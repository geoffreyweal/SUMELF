"""
does_contain_ECCP_Information_data.py, Geoffrey Weal, 25/3/24

This method is designed to check if there already exists ECCP_Information data.
"""
import os
from types import SimpleNamespace
from SUMELF.SUMELF.ECCP_based_methods.utility_methods.get_list_of_files_to_check import get_list_of_files_to_check

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Get the list of files to check for in the ECCP information folder. 
check_for_files = get_list_of_files_to_check()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def what_ECCP_Information_files_do_we_have(path_to_ECCP_Information_folder_for_crystal):
	"""
	This method is designed to check if there already exists ECCP_Information data.

	Parameters
	----------
	path_to_ECCP_Information_folder_for_crystal : str.
		This is the path to the ECCP_Information data. 

	Returns
	-------
	This will return a namespace that includes:

	have_crystal_file : bool.
		This indicates if you have the crystal.xyz file.
	have_ECCP_Information_file : bool.
		This indicates if you have the ECCP_Information.txt file.
	have_All_Dimer_Information_file : bool.
		This indicates if you have the All_Dimer_Information.txt file.

	have_Structurally_Equivalent_Molecule_Groups_file : bool.
		This indicates if you have the Structurally_Equivalent_Molecule_Groups.txt file.
	have_Conformationally_Equivalent_Molecule_Groups_file : bool.
		This indicates if you have the Conformationally_Equivalent_Molecule_Groups.txt file.
	have_Structurally_Equivalent_Dimer_Groups_file : bool.
		This indicates if you have the Structurally_Equivalent_Dimer_Groups.txt file.

	have_Structurally_Equivalent_Molecule_Pairs_file : bool.
		This indicates if you have the Structurally_Equivalent_Molecule_Pairs.txt file.
	have_Structurally_Equivalent_Dimer_Pairs_file : bool.
		This indicates if you have the Structurally_Equivalent_Dimer_Pairs.txt file.
	"""

	# First, create a dictionary for holding which ECCP Information files have been located
	contains_ECCP_Information_files = {check_for_file : False for check_for_file in check_for_files}

	# Second, check which files have been located
	for file_to_check in check_for_files:
		if os.path.exists(path_to_ECCP_Information_folder_for_crystal+'/'+file_to_check):
			contains_ECCP_Information_files[file_to_check] = True

	# Third, if at least one of the files in check_for_files has been found
	return convert_contains_ECCP_Information_files_into_variables(contains_ECCP_Information_files)

def convert_contains_ECCP_Information_files_into_variables(contains_ECCP_Information_files):
	"""
	This method is designed to convert the booleans from contains_ECCP_Information_files into local variables

	Parameters
	----------
	contains_ECCP_Information_files : dict.
		This dictionary contains all the information about which files have been found. 

	Returns
	-------
	This will return a namespace that includes:

	have_crystal_file : bool.
		This indicates if you have the crystal.xyz file.
	have_ECCP_Information_file : bool.
		This indicates if you have the ECCP_Information.txt file.
	have_All_Dimer_Information_file : bool.
		This indicates if you have the All_Dimer_Information.txt file.

	have_Structurally_Equivalent_Molecule_Groups_file : bool.
		This indicates if you have the Structurally_Equivalent_Molecule_Groups.txt file.
	have_Conformationally_Equivalent_Molecule_Groups_file : bool.
		This indicates if you have the Conformationally_Equivalent_Molecule_Groups.txt file.
	have_Structurally_Equivalent_Dimer_Groups_file : bool.
		This indicates if you have the Structurally_Equivalent_Dimer_Groups.txt file.

	have_Structurally_Equivalent_Molecule_Pairs_file : bool.
		This indicates if you have the Structurally_Equivalent_Molecule_Pairs.txt file.
	have_Structurally_Equivalent_Dimer_Pairs_file : bool.
		This indicates if you have the Structurally_Equivalent_Dimer_Pairs.txt file.
	"""

	# First, determine if the following files are found in the contains_ECCP_Information_files dictionary.
	have_crystal_file                                     = contains_ECCP_Information_files['crystal.xyz']
	have_ECCP_Information_file                            = contains_ECCP_Information_files['ECCP_Information.txt']
	have_All_Dimer_Information_file                       = contains_ECCP_Information_files['All_Dimer_Information.txt']
	have_Structurally_Equivalent_Molecule_Groups_file     = contains_ECCP_Information_files['Equivalence_Group_Information/Structurally_Equivalent_Molecule_Groups.txt']
	have_Conformationally_Equivalent_Molecule_Groups_file = contains_ECCP_Information_files['Equivalence_Group_Information/Conformationally_Equivalent_Molecule_Groups.txt']
	have_Structurally_Equivalent_Dimer_Groups_file        = contains_ECCP_Information_files['Equivalence_Group_Information/Structurally_Equivalent_Dimer_Groups.txt']
	have_Structurally_Equivalent_Molecule_Pairs_file      = contains_ECCP_Information_files['Equivalence_Group_Information/Structurally_Equivalent_Molecule_Pairs.txt']
	have_Structurally_Equivalent_Dimer_Pairs_file         = contains_ECCP_Information_files['Equivalence_Group_Information/Structurally_Equivalent_Dimer_Pairs.txt']

	# Second, create a dictionary to hold the namespace for these inputs.
	dictionary_for_namespace = {'have_crystal_file': have_crystal_file, 'have_ECCP_Information_file': have_ECCP_Information_file, 'have_All_Dimer_Information_file': have_All_Dimer_Information_file, 'have_Structurally_Equivalent_Molecule_Groups_file': have_Structurally_Equivalent_Molecule_Groups_file, 'have_Conformationally_Equivalent_Molecule_Groups_file': have_Conformationally_Equivalent_Molecule_Groups_file, 'have_Structurally_Equivalent_Dimer_Groups_file': have_Structurally_Equivalent_Dimer_Groups_file, 'have_Structurally_Equivalent_Molecule_Pairs_file': have_Structurally_Equivalent_Molecule_Pairs_file, 'have_Structurally_Equivalent_Dimer_Pairs_file' : have_Structurally_Equivalent_Dimer_Pairs_file}

	# Third, return the booleans that indicates which files have been found. 
	return dictionary_for_namespace












