"""
get_information_from_Structurally_Equivalent_Pairs.py, Geoffrey Weal, 25/3/24

This method is designed to obtain all the equivalent pairs of molecules/dimers. 
"""

def get_information_from_equivalent_pairs(path_to_equivalent_pairs_file, individual_type, total_number_of_individuals=None):
	"""
	This method is designed to obtain all the equivalent pairs of molecules/dimers. 

	Parameters
	----------
	path_to_equivalent_pairs_file : str.
		This is the path to the file containing the equivalent pairs of molecules/dimers.
	individual_type : str
		This indicates if molecules or dimers are being processed by this method. 
	total_number_of_individuals : int
		The total number of individual expected. This is used for debugging.

	Returns
	-------
	equivalent_pairs : list
		This is the list of equivalent pairs of molecules/dimers.
	"""

	# First, open the text file containing equivalent pairs of molecules/dimers.
	Equivalent_PairsTXT = open(path_to_equivalent_pairs_file, 'r')

	# Second, obtain the list of equivalent pairs of molecules/dimers from the txt file. 
	equivalent_pairs = eval(Equivalent_PairsTXT.readline().rstrip())

	# Third, close the txt file
	Equivalent_PairsTXT.close()

	# Fourth, check that the names of all individual are between 1 and total_number_of_individuals.
	if total_number_of_individuals is not None:
		individual_pair_issues = []
		for indiv1_name, indiv2_name in equivalent_pairs:
			if (not (1 <= indiv1_name <= total_number_of_individuals)) or (not (1 <= indiv2_name <= total_number_of_individuals)):
				individual_pair_issues.append((indiv1_name, indiv2_name))
		if len(individual_pair_issues) > 0:
			to_string  = f'Error: Some of the {individual_type}s in '+str(path_to_equivalent_pairs_file)+' have an unexpected name.\n'
			to_string += f'The names of {individual_type}s should be between 1 and {total_number_of_individuals}.\n'
			to_string += f'{individual_type}_pair_issues: {individual_pair_issues}\n'
			to_string += f'equivalent_pairs: {equivalent_pairs}\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# Fifth, return equivalent_pairs
	return equivalent_pairs