"""
collect_ECCP_Information_data.py, Geoffrey Weal, 9/3/2024

This method is designed to collect all the information from the ECCP_Information for a crystal if it has already completed finished.
"""
import os
from ase.io import read

from SUMELF.SUMELF.ECCP_based_methods.utility_methods.check_if_ECCP_Information_contains_all_necessary_files                   import check_if_ECCP_Information_contains_all_necessary_files
from SUMELF.SUMELF.ECCP_based_methods.collect_ECCP_Information_data_methods.get_information_from_ECCP_InformationTXT_file      import get_information_from_ECCP_InformationTXT_file
from SUMELF.SUMELF.ECCP_based_methods.utility_methods.convert_list_from_names_to_indices                                       import convert_list_from_names_to_indices

from SUMELF.SUMELF.ECCP_based_methods.collect_ECCP_Information_data_methods.get_information_from_All_Dimer_InformationTXT_file import get_information_from_All_Dimer_InformationTXT_file
from SUMELF.SUMELF.ECCP_based_methods.utility_methods.convert_dimers_details_from_name_to_index                                import convert_dimers_details_from_name_to_index

from SUMELF.SUMELF.ECCP_based_methods.collect_ECCP_Information_data_methods.get_information_from_Equivalency_Groups_file       import get_information_from_Equivalency_Groups_file
from SUMELF.SUMELF.ECCP_based_methods.utility_methods.convert_equivalence_groups_from_name_to_index                            import convert_equivalence_groups_from_name_to_index
from SUMELF.SUMELF.ECCP_based_methods.collect_ECCP_Information_data_methods.get_information_from_equivalent_pairs              import get_information_from_equivalent_pairs
from SUMELF.SUMELF.ECCP_based_methods.utility_methods.convert_equivalence_pairs_from_name_to_index                             import convert_equivalence_pairs_from_name_to_index

from SUMELF.SUMELF.ECCP_based_methods.collect_ECCP_Information_data_methods.check_consistancy_between_files                    import check_consistancy_between_files

def collect_ECCP_Information_data(path_to_ECCP_Information_folder_for_crystal):
	"""
	This method is designed to collect all the information from the ECCP_Information for a crystal if it has already completed finished.

	Parameters
	----------
	path_to_ECCP_Information_folder_for_crystal : str.
		This is the path to the ECCP_Information folder for the crystal you want to obtain information about.
	
	Returns
	-------
	crystal : ase.Atoms
		This is the Atoms object for the crystal
	eccp_information : list
		This is all the information from the ECCP_Information.txt file
	symmetric_to_unique_molecules : dict.
		This dictionary indicates which spatially symmetric molecules are equivalent to which unique molecules in the crystal. 
	conformationally_equivalent_to_unique_molecules : dict. 
		This dictionary indicates which conformationally equivalent molecules are equivalent to which unique molecules in the crystal. 
	dimers_details : dict.
		This list contains all the information about all the dimers in the molecules based on the 'nearest_atoms_method' setting of the make_dimer_method dictionary in Run_ECCP.py
	"""

	# First, make sure that this path_to_eccp_folder folder exists
	if not os.path.exists(path_to_ECCP_Information_folder_for_crystal):
		raise Exception('Error: path '+str(path_to_ECCP_Information_folder_for_crystal)+' does not exist. Check this')

	# Second, make sure that the second to bottom folder name is called ECCP_Information
	if not path_to_ECCP_Information_folder_for_crystal.split('/')[-2] == 'ECCP_Information':
		to_string  = 'Error: The second to last directory name expected for using this method is suppose to be ECCP_Information.\n'
		to_string += 'Path given: '+str(path_to_ECCP_Information_folder_for_crystal)+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Third, we want to check if the following files exist. These files are the minimal files needed to completely write all files without needing to perform any processes.
	# * crystal.xyz: Contains all the information about the atoms and molecules in the crystal.
	# * ECCP_Information.txt: Contains some information about --> the molecules, including which molecules are solvents.
	#                                                         --> the dimers that are unique, symmetric, and which dimers are symmetric to which unique molecules.
	# * Structurally_Unique_Molecule_Information.txt: Indicates which molecules are structurally the same to qhich unique molecules
	# * Conformationally_Unique_Molecule_Information.txt: Indicates which molecules are conformationally the same to which unique molecules. 
	# * All_Dimer_Information.txt: Contains information about how to construct all dimers in the molecule. 
	check_if_ECCP_Information_contains_all_necessary_files(path_to_ECCP_Information_folder_for_crystal)

	# Fourth, read in the crystal.xyz file. 
	crystal = read(path_to_ECCP_Information_folder_for_crystal+'/'+'crystal.xyz')

	# Fifth, check that crystal file contains all information about the molecules in the crystal, including
	# * MoleculeList: Indicates which atoms go with which molecule in the All_Molecules folder
	# * NeighboursList: Indicates which atoms are bonded to which atoms in the molecule. 
	check_for_missing_attributes(crystal)

	# Sixth, collect data from the 'ECCP_Information.txt' file.
	eccp_information = get_information_from_ECCP_InformationTXT_file(path_to_ECCP_Information_folder_for_crystal+'/'+'ECCP_Information.txt')
	
	# Seventh, convert the solvent lists for molecules and dimers from names to indices.
	solvent_list               = convert_list_from_names_to_indices(eccp_information[3])
	dimers_containing_solvents = convert_list_from_names_to_indices(eccp_information[9])

	# Eighth, update the eccp_information with solvents given by their indices rather than their names.
	eccp_information = tuple([eccp_information[0], eccp_information[1], eccp_information[2], solvent_list, eccp_information[4], eccp_information[5], eccp_information[6], eccp_information[7], eccp_information[8], dimers_containing_solvents, eccp_information[10], eccp_information[11]])

	# Ninth, obtain the number of molecules and dimers in the crystal and previous ECCP run.
	_, _, no_of_molecules, _, _, _, _, _, no_of_dimers, _,  _, _ = eccp_information

	# Tenth, collect data from the 'All_Dimer_Information.txt' file.
	dimers_details = get_information_from_All_Dimer_InformationTXT_file(path_to_ECCP_Information_folder_for_crystal+'/'+'All_Dimer_Information.txt', total_no_of_dimers=no_of_dimers)
	dimers_details = convert_dimers_details_from_name_to_index(dimers_details, total_no_of_dimers=no_of_dimers)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Eleventh, obtain the molecule and dimer equivalence groups

	# 11.1: Collect data from the 'Structurally_Equivalent_Molecule_Groups.txt' file.
	structurally_equivalent_molecule_groups     = get_information_from_Equivalency_Groups_file(path_to_ECCP_Information_folder_for_crystal+'/Equivalence_Group_Information/Structurally_Equivalent_Molecule_Groups.txt', individual_type='molecule', total_number_of_individuals=no_of_molecules)
	structurally_equivalent_molecule_groups     = convert_equivalence_groups_from_name_to_index(structurally_equivalent_molecule_groups)

	# 1.2: Collect data from the 'Conformationally_Equivalent_Molecule_Groups.txt' file.
	conformationally_equivalent_molecule_groups = get_information_from_Equivalency_Groups_file(path_to_ECCP_Information_folder_for_crystal+'/Equivalence_Group_Information/Conformationally_Equivalent_Molecule_Groups.txt', individual_type='molecule', total_number_of_individuals=no_of_molecules)
	conformationally_equivalent_molecule_groups = convert_equivalence_groups_from_name_to_index(conformationally_equivalent_molecule_groups)

	# 1.3: Collect data from the 'Structurally_Equivalent_Dimer_Groups.txt' file.
	structurally_equivalent_dimer_groups        = get_information_from_Equivalency_Groups_file(path_to_ECCP_Information_folder_for_crystal+'/Equivalence_Group_Information/Structurally_Equivalent_Dimer_Groups.txt', individual_type='dimer', total_number_of_individuals=no_of_dimers)
	structurally_equivalent_dimer_groups        = convert_equivalence_groups_from_name_to_index(structurally_equivalent_dimer_groups)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Twelfth, obtain all the equivalent pairs of molecules and dimers.

	# 12.1: Collect data from the 'Structurally_Equivalent_Molecule_Pairs.txt' file.
	structurally_equivalent_molecule_pairs      = get_information_from_equivalent_pairs(path_to_ECCP_Information_folder_for_crystal+'/Equivalence_Group_Information/Structurally_Equivalent_Molecule_Pairs.txt', individual_type='molecule', total_number_of_individuals=no_of_molecules)
	structurally_equivalent_molecule_pairs      = convert_equivalence_pairs_from_name_to_index(structurally_equivalent_molecule_pairs)

	# 12.1: Collect data from the 'Structurally_Equivalent_Dimer_Pairs.txt' file.
	structurally_equivalent_dimer_pairs         = get_information_from_equivalent_pairs(path_to_ECCP_Information_folder_for_crystal+'/Equivalence_Group_Information/Structurally_Equivalent_Dimer_Pairs.txt', individual_type='dimer', total_number_of_individuals=no_of_dimers)
	structurally_equivalent_dimer_pairs         = convert_equivalence_pairs_from_name_to_index(structurally_equivalent_dimer_pairs)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Thirteenth, check that all the files are consistent with each other. 
	check_consistancy_between_files(crystal, eccp_information, dimers_details, structurally_equivalent_molecule_groups, conformationally_equivalent_molecule_groups, structurally_equivalent_dimer_groups, structurally_equivalent_molecule_pairs, structurally_equivalent_dimer_pairs)

	# Fourteenth, return crystal, eccp_information, dimers_details, structurally_equivalent_molecule_groups, conformationally_equivalent_molecule_groups, structurally_equivalent_dimer_groups, structurally_equivalent_molecule_pairs, and structurally_equivalent_dimer_pairs
	return crystal, eccp_information, dimers_details, structurally_equivalent_molecule_groups, conformationally_equivalent_molecule_groups, structurally_equivalent_dimer_groups, structurally_equivalent_molecule_pairs, structurally_equivalent_dimer_pairs

missing_attributes_to_check = ['MoleculeList', 'NeighboursList']
def check_for_missing_attributes(crystal):
	"""
	This method is designed to check if the crystal contains the MoleculeList and NeighboursList lists
	"""

	# First, initialise a list to add missing attributes to
	missing_attributes = []

	# Second, for eaxch attribute to check.
	for missing_attribute_to_check in missing_attributes_to_check:

		# 2.1: Check if missing_attribute_to_check is in the crystal object
		if missing_attribute_to_check not in crystal.arrays.keys():

			# 2.2: Add missing_attribute_to_check to missing_attributes
			missing_attributes.append(missing_attribute_to_check)

	# Third, raise an error if the crystal object does not contain one or more attributes from missing_attributes_to_check.
	if len(missing_attributes) > 0:
		to_string  = 'Error: You are trying to read in the "crystal.xyz" file, however it does not have the following attributes: '+str(missing_attributes)+'\n'
		to_string += 'Check this. This program will finish without beginning'
		raise Exception(to_string)















