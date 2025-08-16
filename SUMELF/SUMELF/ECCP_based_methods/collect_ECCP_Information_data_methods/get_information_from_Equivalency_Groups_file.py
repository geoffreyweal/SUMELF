"""
get_information_from_Structurally_Unique_Molecule_InformationTXT_file.py, Geoffrey Weal, 9/3/2024

This method is designed to collect all the information needed from the Structurally_Unique_Molecule_Information.txt file.
"""

def get_information_from_Equivalency_Groups_file(path_to_equivalence_group_file, individual_type, total_number_of_individuals=None):
	"""
	This method is designed to collect all the information needed from the path_to_equivalence_group_file file.

	Parameters
	----------
	path_to_equivalence_group_file : str.
		This is the path to the file containing the equivalence group data.
	individual_type : str
		This indicates if molecules or dimers are being processed by this method. 
	total_number_of_individuals : int
		The total number of individual expected. This is used for debugging.

	Returns
	-------
	equivalence_groups : list
		This is the list of equivalence groups of molecules or dimers. 
	"""

	# First, initialise a list to store equivalence groups from path_to_equivalence_group_file into. 
	equivalence_groups = []

	# Second, collect the data from the Structurally_Unique_Molecule_Information text file. 
	with open(path_to_equivalence_group_file, 'r') as Equivalence_GroupTXT:

		# 2.1: Collect the list of equivalence groups from path_to_equivalence_group_file
		for line in Equivalence_GroupTXT:

			# 2.2: Strip the end line '\n' from line.
			line = line.rstrip()

			# 2.3: Remove the initial [/( from the line.
			if line.startswith('(') or line.startswith('['):
				line = line[1:]

			# 2.4: Remove the final )\] from the line.
			if line.endswith(')') or line.endswith(']'):
				line = line[:-1]

			# 2.5: Obtain all the components of the equivalence group.
			equivalence_group = [int(individual) for individual in line.split(',')]

			# 2.6: Add equivalency_group to equivalency_groups
			equivalence_groups.append(equivalence_group)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Third, check that each equivalence group does not have duplicates
	duplicates_in_equivalence_group_list = []
	for equivalence_group in equivalence_groups:
		if not len(equivalence_group) == len(set(equivalence_group)):
			duplicates_in_equivalence_group = [individual for individual, count in Counter(equivalence_group).items() if (count >= 2)]
			duplicates_in_equivalence_group_list.append((equivalence_group, duplicates_in_equivalence_group))
	if len(duplicates_in_equivalence_group_list) > 0:
		to_string  = f'Error: There are duplicate {individual_type}s in one or more equivalence_group.'+'\n'
		for equivalence_group, duplicates_in_equivalence_group in duplicates_in_equivalence_group_list:
			to_string += f'Equivalence Group: {equivalence_group}; Duplicates: {duplicates_in_equivalence_group}\n'
		to_string += f'All Equivalence Groups: {equivalence_groups}\n'
		to_string += 'Check this'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Fourth, obtain all the individuals from equivalence_groups
	all_individuals_from_equivalence_groups = [j for sub in equivalence_groups for j in sub]

	# Fifth, check that there are no duplicates in the full equivalence_groups
	if not len(all_individuals_from_equivalence_groups) == len(set(all_individuals_from_equivalence_groups)):
		to_string  = f'Error: There are repeated molecules across all {individual_type} equivalence group.'
		duplicates_in_all_individuals_from_equivalence_groups = [individual for individual, count in Counter(all_individuals_from_equivalence_groups).items() if (count >= 2)]
		to_string += f'Duplicate {individual_type}s: {duplicates_in_all_individuals_from_equivalence_groups}\n'
		to_string += f'All Equivalence Groups: {equivalence_groups}\n'
		to_string += 'Check this'
		raise Exception(to_string)

	# Sixth, make sure that all individual expected are in equivalence_groups
	if total_number_of_individuals is not None:
		if not sorted(all_individuals_from_equivalence_groups) == list(range(1,total_number_of_individuals+1)):
			to_string  = f'Error: There are missing {individual_type}s compared to that expected\n'
			to_string += f'All {individual_type}s in Equivalence Group: {sorted(all_individuals_from_equivalence_groups)}\n'
			to_string += f'{individual_type.title()}s Expected: {list(range(1,total_number_of_individuals+1))}\n'

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Seventh, return equivalence_groups
	return equivalence_groups

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 













