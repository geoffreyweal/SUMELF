"""
get_equivalent_molecule_group_data.py, Geoffrey Weal, 27/3/24

This method is designed to obtain all the information about the structurally and conformationally equivalent molecule groups in the crystal. This is obtained from the ECCP_Information folder. 
"""
from SUMELF.SUMELF.ECCP_based_methods.collect_ECCP_Information_data_methods.get_information_from_Equivalency_Groups_file import get_information_from_Equivalency_Groups_file
#from SUMELF.SUMELF.ECCP_based_methods.utility_methods.convert_equivalence_groups_from_name_to_index                      import convert_equivalence_groups_from_name_to_index

from SUMELF.SUMELF.ECCP_based_methods.collect_ECCP_Information_data_methods.get_information_from_equivalent_pairs        import get_information_from_equivalent_pairs
#from SUMELF.SUMELF.ECCP_based_methods.utility_methods.convert_equivalence_pairs_from_name_to_index                       import convert_equivalence_pairs_from_name_to_index


def get_equivalent_molecule_group_data(path_to_ECCP_Information_folder_for_crystal, eccp_information):
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
	structurally_equivalent_molecule_groups : dict. 
		These are the structurally equivalent molecule groups in the crystal. 
	conformationally_equivalent_molecule_groups : dict. 
		These are the conformationally equivalent molecule groups in the crystal. 
	structurally_equivalent_molecule_pairs : list. 
		These are all the pairs of structurally equivalent molecules in the crystal. 
	"""

	eccp_information['number_of_molecules']

	# First, obtain the number of molecules as given in the eccp_information file. 
	if 'number_of_molecules' in eccp_information:
		number_of_molecules = eccp_information['number_of_molecules']
	else:
		number_of_molecules = None

	# First, collect data from the 'Structurally_Equivalent_Molecule_Groups.txt' file.
	structurally_equivalent_molecule_groups     = get_information_from_Equivalency_Groups_file(path_to_ECCP_Information_folder_for_crystal+'/Equivalence_Group_Information/Structurally_Equivalent_Molecule_Groups.txt', individual_type='molecule', total_number_of_individuals=number_of_molecules)
	#structurally_equivalent_molecule_groups     = convert_equivalence_groups_from_name_to_index(structurally_equivalent_molecule_groups)

	# Second, collect data from the 'Conformationally_Equivalent_Molecule_Groups.txt' file.
	conformationally_equivalent_molecule_groups = get_information_from_Equivalency_Groups_file(path_to_ECCP_Information_folder_for_crystal+'/Equivalence_Group_Information/Conformationally_Equivalent_Molecule_Groups.txt', individual_type='molecule', total_number_of_individuals=number_of_molecules)
	#conformationally_equivalent_molecule_groups = convert_equivalence_groups_from_name_to_index(conformationally_equivalent_molecule_groups)

	# Third, collect data from the 'Structurally_Equivalent_Molecule_Pairs.txt' file.
	structurally_equivalent_molecule_pairs      = get_information_from_equivalent_pairs(path_to_ECCP_Information_folder_for_crystal+'/Equivalence_Group_Information/Structurally_Equivalent_Molecule_Pairs.txt', individual_type='molecule', total_number_of_individuals=number_of_molecules)
	#structurally_equivalent_molecule_pairs      = convert_equivalence_pairs_from_name_to_index(structurally_equivalent_molecule_pairs)

	# Fourth, return structurally_equivalent_molecule_groups, conformationally_equivalent_molecule_groups, and structurally_equivalent_molecule_pairs
	return structurally_equivalent_molecule_groups, conformationally_equivalent_molecule_groups, structurally_equivalent_molecule_pairs