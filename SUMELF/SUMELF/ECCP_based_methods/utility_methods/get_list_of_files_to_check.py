"""
get_list_of_files_to_check.py, Geoffrey Weal, 25/3/24

This method is designed to check if there already exists ECCP_Information data.
"""

def get_list_of_files_to_check():
	"""
	This method will give a list of all the ECCP Information files that need to be checked.

	Parameters
	----------
	check_for_files : list
		These are all the files that should be checked. 
	"""

	# First, obtain the main files in the ECCP_Information folder
	check_for_files = ['crystal.xyz', 'ECCP_Information.txt', 'All_Dimer_Information.txt']

	# Second, get the files to check in the Equivalence_Group_Information
	equivalency_group_files = ['Structurally_Equivalent_Molecule_Groups.txt', 'Conformationally_Equivalent_Molecule_Groups.txt', 'Structurally_Equivalent_Molecule_Pairs.txt', 'Structurally_Equivalent_Dimer_Groups.txt', 'Structurally_Equivalent_Dimer_Pairs.txt']
	for equivalency_group_file in equivalency_group_files:
		check_for_files.append('Equivalence_Group_Information'+'/'+equivalency_group_file)

	# Third, return check_for_files
	return check_for_files