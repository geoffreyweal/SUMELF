"""
get_equivalent_dimer_group_data.py, Geoffrey Weal, 27/3/24

This method is designed to obtain all the information about the structurally and conformationally equivalent dimer groups in the crystal. This is obtained from the ECCP_Information folder. 
"""
from SUMELF.SUMELF.ECCP_based_methods.collect_ECCP_Information_data_methods.get_information_from_Equivalency_Groups_file import get_information_from_Equivalency_Groups_file
#from SUMELF.SUMELF.ECCP_based_methods.utility_methods.convert_equivalence_groups_from_name_to_index                      import convert_equivalence_groups_from_name_to_index

from SUMELF.SUMELF.ECCP_based_methods.collect_ECCP_Information_data_methods.get_information_from_equivalent_pairs        import get_information_from_equivalent_pairs
#from SUMELF.SUMELF.ECCP_based_methods.utility_methods.convert_equivalence_pairs_from_name_to_index                       import convert_equivalence_pairs_from_name_to_index

def get_equivalent_dimer_group_data(path_to_ECCP_Information_folder_for_crystal, eccp_information):
	"""
	This method is designewd to obtain all the information about the structurally and conformationally equivalent molecule groups in the crystal. This is obtained from the ECCP_Information folder. 

	Parameters
	----------
	path_to_ECCP_Information_folder_for_crystal : str.
		This is the path to the ECCP_Information folder for the crystal you want to obtain information about the structurally and conformationally equivalent molecule groups.
	eccp_information : dict.
		This dictionary contains all the data from the previous ECCP run. Included in this is the number of dimers that are expected. Used for double-checking purposes.

	Returns
	-------
	structurally_equivalent_dimer_groups : dict. 
		These are the structurally equivalent dimer groups from your previous ECCP run.
	structurally_equivalent_dimer_pairs : list. 
		These are all the pairs of structurally equivalent dimers from your previous ECCP run.
	"""

	# First, obtain the number of dimers as given in the eccp_information file. 
	if 'number_of_dimers' in eccp_information:
		number_of_dimers = eccp_information['number_of_dimers']
	else:
		number_of_dimers = None

	# Second, collect data from the 'Structurally_Equivalent_Dimer_Groups.txt' file.
	structurally_equivalent_dimer_groups        = get_information_from_Equivalency_Groups_file(path_to_ECCP_Information_folder_for_crystal+'/Equivalence_Group_Information/Structurally_Equivalent_Dimer_Groups.txt', individual_type='dimer', total_number_of_individuals=number_of_dimers)
	#structurally_equivalent_dimer_groups        = convert_equivalence_groups_from_name_to_index(structurally_equivalent_dimer_groups)

	# Third, collect data from the 'Structurally_Equivalent_Dimer_Pairs.txt' file.
	structurally_equivalent_dimer_pairs         = get_information_from_equivalent_pairs(path_to_ECCP_Information_folder_for_crystal+'/Equivalence_Group_Information/Structurally_Equivalent_Dimer_Pairs.txt', individual_type='dimer', total_number_of_individuals=number_of_dimers)
	#structurally_equivalent_dimer_pairs         = convert_equivalence_pairs_from_name_to_index(structurally_equivalent_dimer_pairs)

	# Fourth, return structurally_equivalent_dimer_groups, and structurally_equivalent_dimer_pairs
	return structurally_equivalent_dimer_groups, structurally_equivalent_dimer_pairs