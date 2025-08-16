"""
check_consistancy_between_files.py, Geoffrey Weal, 9/3/2024

This program is designed to check the ECCP files against each other and make sure they are consistent with each other. 
"""
from copy import deepcopy
from networkx import Graph, is_connected

def check_consistancy_between_files(have_ECCP_Information_file, have_ECCP_Information_crystal_file, has_neighbouring_molecules, has_unique_molecules, has_unique_dimers, crystal=None, eccp_information=None, dimer_details=None, structurally_equivalent_molecule_groups=None, conformationally_equivalent_molecule_groups=None, structurally_equivalent_dimer_groups=None, structurally_equivalent_molecule_pairs=None, structurally_equivalent_dimer_pairs=None):
	"""
	This program is designed to check the ECCP files against each other and make sure they are consistent with each other. 

	Parameters
	----------
	have_ECCP_Information_file : bool.
		
	have_ECCP_Information_crystal_file : bool.

	has_neighbouring_molecules : bool.

	has_unique_molecules : bool.

	has_unique_dimers : bool.


	crystal : ase.Atoms
		This is the ase.Atoms object of the crystal
	eccp_information : list 
		This is the information from the ECCP_Information.txt file. 
	dimer_details : dict.
		This dictionary contains all the information about how to construct the dimers from the molecules in this crystal for the 'max_dimer_distance' setting in make_dimer_method

	structurally_equivalent_molecule_groups : list.
		These are the structurally equivalent molecule groups. 
	conformationally_equivalent_molecule_groups : list.
		These are the conformationally equivalent molecule groups. 
	structurally_equivalent_dimer_groups : list.
		These are the structurally equivalent dimer groups. 

	structurally_equivalent_molecule_pairs : list.
		These are all the pairs of structurally equivalent molecules. 
	structurally_equivalent_dimer_pairs : list.
		These are all the pairs of structurally equivalent dimers. 
	"""

	# ==========================================================================================================================================================
	# Part 1: Extract the information we need to check

	# First, extract information from the crystal object
	if have_ECCP_Information_crystal_file:
		MoleculeList   = crystal.get_array('MoleculeList')
		NeighboursList = crystal.get_array('NeighboursList')

	# Second, extract the data from eccp_information. This information includes:
	# * make_dimer_method: This dictionary indicates the settings for how dimers were created by the Run_ECCP.py script
	# * environment_settings: This dictionary indicates the environmental settings that indicates how molecules surrounding moelcules and dimers are included in calculations. 
	# * number_of_molecules: This is the number of molecules in the crystal. 
	# * solvent_names: This is the names of the solvent molecules in the crystal.
	# * number_of_structurally_unique_molecules: This is the number of structurally unique molecules.
	# * number_of_structurally_equivalent_molecules: This is the number of structurally equivalent molecules.
	# * number_of_conformationally_unique_molecules: This is the number of conformationally unique molecules. 
	# * number_of_conformationally_equivalent_molecules: This is the number of conformationally equivalent molecules. 
	# * number_of_dimers: This is the number of dimers found for the given value of 'max_dimer_distance' in the make_dimer_method dictionary. 
	# * dimers_that_contain_solvents: These are the dimers that contain solvent molecules in them. 
	# * number_of_structurally_unique_dimers: This is the number of structurally unique dimers.
	# * number_of_structurally_equivalent_dimers: This is the number of structurally equivalent dimers.
	if have_ECCP_Information_file:
		make_dimer_method                               = eccp_information['make_dimer_method']
		environment_settings                            = eccp_information['environment_settings']
		number_of_molecules                             = eccp_information['number_of_molecules']
		solvent_names                                   = eccp_information['solvent_names']
		number_of_structurally_unique_molecules         = eccp_information['number_of_structurally_unique_molecules']
		number_of_structurally_equivalent_molecules     = eccp_information['number_of_structurally_equivalent_molecules']
		number_of_conformationally_unique_molecules     = eccp_information['number_of_conformationally_unique_molecules']
		number_of_conformationally_equivalent_molecules = eccp_information['number_of_conformationally_equivalent_molecules']
		number_of_dimers                                = eccp_information['number_of_dimers']
		dimers_that_contain_solvents                    = eccp_information['dimers_that_contain_solvents']
		number_of_structurally_unique_dimers            = eccp_information['number_of_structurally_unique_dimers']
		number_of_structurally_equivalent_dimers        = eccp_information['number_of_structurally_equivalent_dimers']

	# ==========================================================================================================================================================
	# Part 2: Compare the data between crystal, eccp_information, structurally_equivalent_molecule_groups, conformationally_equivalent_molecule_groups, 
	#         structurally_equivalent_dimer_groups, structurally_equivalent_molecule_pairs, structurally_equivalent_dimer_pairs, and make sure they are consistent 
	#         between each other. 

	# Section 2.1: Checks involving the ECCP_Information.txt file

	# Third, make sure that the number of molecules in the crystal object (in MoleculeList) is the same as given in the ECCP_Information.txt file. 
	if have_ECCP_Information_crystal_file:
		if not number_of_molecules == len(set(MoleculeList)):
			to_string  = 'Error: The number of molecules in the MoleculeList of the crystal.xyz file is not the same as given in "Number of molecules obtained from crystal" in ECCP_Information.txt\n'
			to_string += 'Number of molecules given in crystal.xyz file: '+str(len(set(MoleculeList)))+'\n'
			to_string += 'Number of molecules given in the ECCP_Information.txt file: '+str(number_of_molecules)+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
	# Fourth, make sure that the structurally equivalent molecule groups are the same between the Structurally_Equivalent_Molecule_Groups.txt and ECCP_Information.txt file.

	if have_ECCP_Information_file and has_unique_molecules:

		# 4.1: Get all the molecules from structurally_equivalent_molecule_groups
		all_molecules_in_structurally_equivalent_molecule_groups = [j for sub in structurally_equivalent_molecule_groups for j in sub]

		# 4.2: Make sure that there are the same number of molecules in structurally_equivalent_molecule_groups (in Structurally_Equivalent_Molecule_Groups.txt) and number_of_molecules (in ECCP_Information.txt).
		if not len(all_molecules_in_structurally_equivalent_molecule_groups) == number_of_molecules:
			to_string  = 'Error: There are different numbers of molecules in Structurally_Equivalent_Molecule_Groups.txt and ECCP_Information.txt'
			to_string += f'number_of_molecules (from ECCP_Information.txt): {number_of_molecules}\n'
			to_string += f'Number of molecules in structurally_equivalent_molecule_groups (from Structurally_Equivalent_Molecule_Groups.txt): {len(all_molecules_in_structurally_equivalent_molecule_groups)}\n'
			to_string += f'structurally_equivalent_molecule_groups (from Structurally_Equivalent_Molecule_Groups.txt): {structurally_equivalent_molecule_groups}\n'
			to_string += 'Check your Structurally_Equivalent_Molecule_Groups.txt and ECCP_Information.txt files\n'
			to_string += 'Check this.'
			raise Exception(to_string)

		# 4.3: Check that all the molecules in Structurally_Equivalent_Molecule_Groups.txt have expected name as based on number_of_molecules
		if not sorted(all_molecules_in_structurally_equivalent_molecule_groups) == list(range(1,number_of_molecules+1)):
			to_string  = 'Error: There are some molecules with unexpected names (or duplicates) in Structurally_Equivalent_Molecule_Groups.txt'
			to_string += f'All molecules in Structurally_Equivalent_Molecule_Groups.txt (given by their names): {sorted(all_molecules_in_structurally_equivalent_molecule_groups)}\n'
			to_string += f'Expected molecule names: {list(range(number_of_molecules))}\n'
			to_string += 'Check your Structurally_Equivalent_Molecule_Groups.txt and ECCP_Information.txt files\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	if has_unique_molecules:

		# 4.4: Compare the number of structral equivalence groups with the number of structurally unique molecules.
		if not len(structurally_equivalent_molecule_groups) == number_of_structurally_unique_molecules:
			to_string  = 'Error: The number of structually equivalent molecule groups in Structurally_Equivalent_Molecule_Groups.txt does not reflect the number of structurally unique molecules in ECCP_Information.txt.\n'
			to_string += f'number_of_structurally_unique_molecules (from ECCP_Information.txt): {number_of_structurally_unique_molecules}\n'
			to_string += f'len(structurally_equivalent_molecule_groups) (from Structurally_Equivalent_Molecule_Groups.txt, given by their names): {len(structurally_equivalent_molecule_groups)}\n'
			to_string += f'structurally_equivalent_molecule_groups (from Structurally_Equivalent_Molecule_Groups.txt, given by their names): {structurally_equivalent_molecule_groups}\n'
			to_string += 'Check your Structurally_Equivalent_Molecule_Groups.txt and ECCP_Information.txt files\n'
			to_string += 'Check this.'
			raise Exception(to_string)

		# 4.5: Compare the number of structurally equivalent molecules between the Structurally_Equivalent_Molecule_Groups.txt and ECCP_Information.txt files.
		all_structurally_equivalent_molecules = []
		for structurally_equivalent_molecule_group in structurally_equivalent_molecule_groups:
			structurally_equivalent_molecules = deepcopy(structurally_equivalent_molecule_group)
			structurally_equivalent_molecules.remove(min(structurally_equivalent_molecule_group))
			all_structurally_equivalent_molecules += list(structurally_equivalent_molecules)
		if not len(all_structurally_equivalent_molecules) == number_of_structurally_equivalent_molecules:
			to_string  = 'Error: The number of structually equivalent molecules in Structurally_Equivalent_Molecule_Groups.txt does not reflect the number of structurally equivalent molecules in ECCP_Information.txt.\n'
			to_string += f'number_of_structurally_equivalent_molecules (from ECCP_Information.txt): {number_of_structurally_equivalent_molecules}\n'
			to_string += f'len(all_structurally_equivalent_molecules) (from Structurally_Equivalent_Molecule_Groups.txt): {len(all_structurally_equivalent_molecules)}\n'
			to_string += f'all_structurally_equivalent_molecules (from Structurally_Equivalent_Molecule_Groups.txt, given by their names): {all_structurally_equivalent_molecules}\n'
			to_string += f'structurally_equivalent_molecule_groups (from Structurally_Equivalent_Molecule_Groups.txt, given by their names): {structurally_equivalent_molecule_groups}\n'
			to_string += 'Check your Structurally_Equivalent_Molecule_Groups.txt and ECCP_Information.txt files\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Fifth, make sure that the conformationally equivalent molecule groups are the same between the Conformational_Equivalent_Molecule_Groups.txt and ECCP_Information.txt file.

	if have_ECCP_Information_file and has_unique_molecules:

		# 5.1: Get all the molecules from structurally_equivalent_molecule_groups
		all_molecules_in_conformationally_equivalent_molecule_groups = [j for sub in conformationally_equivalent_molecule_groups for j in sub]

		# 5.2: Make sure that there are the same number of molecules in conformationally_equivalent_molecule_groups (in Conformationally_Equivalent_Molecule_Groups.txt) and number_of_molecules (in ECCP_Information.txt).
		if not len(all_molecules_in_conformationally_equivalent_molecule_groups) == number_of_molecules:
			to_string  = 'Error: There are different numbers of molecules in Conformationally_Equivalent_Molecule_Groups.txt and ECCP_Information.txt'
			to_string += f'number_of_molecules (from ECCP_Information.txt): {number_of_molecules}\n'
			to_string += f'Number of molecules in conformationally_equivalent_molecule_groups (from Conformationally_Equivalent_Molecule_Groups.txt): {len(all_molecules_in_conformationally_equivalent_molecule_groups)}\n'
			to_string += f'conformationally_equivalent_molecule_groups (from Conformationally_Equivalent_Molecule_Groups.txt, given by their names): {conformationally_equivalent_molecule_groups}\n'
			to_string += 'Check your Conformationally_Equivalent_Molecule_Groups.txt and ECCP_Information.txt files\n'
			to_string += 'Check this.'
			raise Exception(to_string)

		# 5.3: Check that all the molecules in Conformationally_Equivalent_Molecule_Groups.txt have expected name as based on number_of_molecules
		if not sorted(all_molecules_in_conformationally_equivalent_molecule_groups) == list(range(1,number_of_molecules+1)):
			to_string  = 'Error: There are some molecules with unexpected names (or duplicates) in Conformationally_Equivalent_Molecule_Groups.txt'
			to_string += f'All molecules in Conformationally_Equivalent_Molecule_Groups.txt: {sorted(all_molecules_in_conformationally_equivalent_molecule_groups)}\n'
			to_string += f'Expected molecule names: {list(range(number_of_molecules))}\n'
			to_string += 'Check your Conformationally_Equivalent_Molecule_Groups.txt and ECCP_Information.txt files\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	if has_unique_molecules:

		# 5.4: Compare the number of conformational equivalence groups with the number of conformational unique molecules.
		if not len(conformationally_equivalent_molecule_groups) == number_of_conformationally_unique_molecules:
			to_string  = 'Error: The number of conformationally equivalent molecule groups in Conformational_Equivalent_Molecule_Groups.txt does not reflect the number of conformational unique molecules in ECCP_Information.txt.\n'
			to_string += f'number_of_conformationally_unique_molecules (from ECCP_Information.txt): {number_of_conformationally_unique_molecules}\n'
			to_string += f'len(conformationally_equivalent_molecule_groups) (from Conformational_Equivalent_Molecule_Groups.txt): {len(conformationally_equivalent_molecule_groups)}\n'
			to_string += f'conformationally_equivalent_molecule_groups (from Conformational_Equivalent_Molecule_Groups.txt, given by their names): {conformationally_equivalent_molecule_groups}\n'
			to_string += 'Check your Conformational_Equivalent_Molecule_Groups.txt and ECCP_Information.txt files\n'
			to_string += 'Check this.'
			raise Exception(to_string)

		# 5.5: Compare the number of conformationally equivalent molecules between the Conformational_Equivalent_Molecule_Groups.txt and ECCP_Information.txt files.
		all_conformationally_equivalent_molecules = []
		for conformationally_equivalent_molecule_group in conformationally_equivalent_molecule_groups:
			conformationally_equivalent_molecules = deepcopy(conformationally_equivalent_molecule_group)
			conformationally_equivalent_molecules.remove(min(conformationally_equivalent_molecule_group))
			all_conformationally_equivalent_molecules += list(conformationally_equivalent_molecules)
		if not len(all_conformationally_equivalent_molecules) == number_of_conformationally_equivalent_molecules:
			to_string  = 'Error: The number of structually equivalent molecules in Conformational_Equivalent_Molecule_Groups.txt does not reflect the number of structurally equivalent molecules in ECCP_Information.txt.\n'
			to_string += f'number_of_conformationally_equivalent_molecules (from ECCP_Information.txt): {number_of_conformationally_equivalent_molecules}\n'
			to_string += f'len(all_conformationally_equivalent_molecules) (from Conformational_Equivalent_Molecule_Groups.txt): {len(all_conformationally_equivalent_molecules)}\n'
			to_string += f'all_conformationally_equivalent_molecules (from Conformational_Equivalent_Molecule_Groups.txt, given by their names): {all_conformationally_equivalent_molecules}\n'
			to_string += f'conformationally_equivalent_molecule_groups (from Conformational_Equivalent_Molecule_Groups.txt, given by their names): {conformationally_equivalent_molecule_groups}\n'
			to_string += 'Check your Conformational_Equivalent_Molecule_Groups.txt and ECCP_Information.txt files\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \
	# Sixth, make sure that the structurally equivalent dimer groups are the same between the Structurally_Equivalent_Dimer_Groups.txt and ECCP_Information.txt file.
	
	if have_ECCP_Information_file and has_unique_dimers:

		# 6.1: Get all the dimers from structurally_equivalent_dimer_groups
		all_dimers_in_structurally_equivalent_dimer_groups = [j for sub in structurally_equivalent_dimer_groups for j in sub]

		# 6.2: Make sure that there are the same number of dimers in structurally_equivalent_dimer_groups (in Structurally_Equivalent_Dimer_Groups.txt) and number_of_dimers (in ECCP_Information.txt).
		if not len(all_dimers_in_structurally_equivalent_dimer_groups) == number_of_dimers:
			to_string  = 'Error: There are different numbers of dimers in Structurally_Equivalent_Dimer_Groups.txt and ECCP_Information.txt'
			to_string += f'number_of_dimers (from ECCP_Information.txt): {number_of_dimers}\n'
			to_string += f'Number of dimers in structurally_equivalent_dimer_groups (from Structurally_Equivalent_Dimer_Groups.txt): {len(all_dimers_in_structurally_equivalent_dimer_groups)}\n'
			to_string += f'structurally_equivalent_dimer_groups (from Structurally_Equivalent_Dimer_Groups.txt, given by their names): {structurally_equivalent_dimer_groups}\n'
			to_string += 'Check your Structurally_Equivalent_Dimer_Groups.txt and ECCP_Information.txt files\n'
			to_string += 'Check this.'
			raise Exception(to_string)

		# 6.3: Check that all the dimers in Structurally_Equivalent_Dimer_Groups.txt have expected name as based on number_of_dimers
		if not sorted(all_dimers_in_structurally_equivalent_dimer_groups) == list(range(1,number_of_dimers+1)):
			to_string  = 'Error: There are some dimers with unexpected names (or duplicates) in Structurally_Equivalent_Dimer_Groups.txt'
			to_string += f'All dimers in Structurally_Equivalent_Dimer_Groups.txt: {sorted(all_dimers_in_structurally_equivalent_dimer_groups)}\n'
			to_string += f'Expected molecule names: {list(range(number_of_dimers))}\n'
			to_string += 'Check your Structurally_Equivalent_Dimer_Groups.txt and ECCP_Information.txt files\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	if has_unique_dimers:

		# 6.4: Compare the number of structral equivalence groups with the number of structurally unique dimers.
		if not len(structurally_equivalent_dimer_groups) == number_of_structurally_unique_dimers:
			to_string  = 'Error: The number of structually equivalent dimer groups in Structurally_Equivalent_Dimer_Groups.txt does not reflect the number of structurally unique dimers in ECCP_Information.txt.\n'
			to_string += f'number_of_structurally_unique_dimers (from ECCP_Information.txt): {number_of_structurally_unique_dimers}\n'
			to_string += f'len(structurally_equivalent_dimer_groups) (from Structurally_Equivalent_Dimer_Groups.txt): {len(structurally_equivalent_dimer_groups)}\n'
			to_string += f'structurally_equivalent_dimer_groups (from Structurally_Equivalent_Dimer_Groups.txt, given by their names): {structurally_equivalent_dimer_groups}\n'
			to_string += 'Check your Structurally_Equivalent_Dimer_Groups.txt and ECCP_Information.txt files\n'
			to_string += 'Check this.'
			raise Exception(to_string)

		# 6.5: Compare the number of structurally equivalent dimers between the Structurally_Equivalent_Dimer_Groups.txt and ECCP_Information.txt files.
		all_structurally_equivalent_dimers = []
		for structurally_equivalent_dimer_group in structurally_equivalent_dimer_groups:
			structurally_equivalent_dimers = deepcopy(structurally_equivalent_dimer_group)
			structurally_equivalent_dimers.remove(min(structurally_equivalent_dimer_group))
			all_structurally_equivalent_dimers += list(structurally_equivalent_dimers)
		if not len(all_structurally_equivalent_dimers) == number_of_structurally_equivalent_dimers:
			to_string  = 'Error: The number of structually equivalent dimers in Structurally_Equivalent_Dimer_Groups.txt does not reflect the number of structurally equivalent dimers in ECCP_Information.txt.\n'
			to_string += f'number_of_structurally_equivalent_dimers (from ECCP_Information.txt): {number_of_structurally_equivalent_dimers}\n'
			to_string += f'len(all_structurally_equivalent_dimers) (from Structurally_Equivalent_Dimer_Groups.txt): {len(all_structurally_equivalent_dimers)}\n'
			to_string += f'all_structurally_equivalent_dimers (from Structurally_Equivalent_Dimer_Groups.txt, given by their names): {all_structurally_equivalent_dimers}\n'
			to_string += f'structurally_equivalent_dimer_groups (from Structurally_Equivalent_Dimer_Groups.txt, given by their names): {structurally_equivalent_dimer_groups}\n'
			to_string += 'Check your Structurally_Equivalent_Dimer_Groups.txt and ECCP_Information.txt files\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Seventh, make sure that the names of all molecules are between 0 and number_of_molecules-1.
	if has_unique_molecules:
		molecule_pair_issues = []
		for mol1_name, mol2_name in structurally_equivalent_molecule_pairs:
			if (not (1 <= mol1_name <= number_of_molecules)) or (not (1 <= mol2_name <= number_of_molecules)):
				molecule_pair_issues.append((mol1_name, mol2_name))
		if len(molecule_pair_issues) > 0:
			to_string  = f'Error: Some of the molecules in structurally_equivalent_molecule_pairs have an unexpected names.\n'
			to_string += f'The names of molecules should be between 1 and {number_of_molecules}.\n'
			to_string += f'molecule_pair_issues: {molecule_pair_issues}\n'
			to_string += f'structurally_equivalent_molecule_pairs: {structurally_equivalent_molecule_pairs}\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# Eighth, make sure that the names of all dimers are between 0 and number_of_dimers-1.
	if has_unique_dimers:
		dimer_pair_issues = []
		for dimer1_name, dimer2_name in structurally_equivalent_dimer_pairs:
			if (not (1 <= dimer1_name <= number_of_dimers)) or (not (1 <= dimer2_name <= number_of_dimers)):
				dimer_pair_issues.append((dimer1_name, dimer2_name))
		if len(dimer_pair_issues) > 0:
			to_string  = f'Error: Some of the dimers in structurally_equivalent_dimer_pairs have an unexpected names.\n'
			to_string += f'The names of dimers should be between 1 and {number_of_dimers}.\n'
			to_string += f'dimer_pair_issues: {dimer_pair_issues}\n'
			to_string += f'structurally_equivalent_dimer_pairs: {structurally_equivalent_dimer_pairs}\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Ninth, check that dimer_details contains the number of dimers given in the All_Dimer_Information.txt file is the same as that given by the ECCP_Information.txt file.
	if has_neighbouring_molecules:
		if not number_of_dimers == len(dimer_details):
			to_string  = 'Error: The number of dimers given in the All_Dimer_Information.txt file is not the same as that given by the ECCP_Information.txt file.\n'
			to_string += 'Number of dimers given in the All_Dimer_Information.txt file: '+str(len(dimer_details))+'\n'
			to_string += 'Number of dimers given in the ECCP_Information.txt file: '+str(number_of_dimers)+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# ----------------------------------------------------------------------------------------------------------------------------------------------------------
	# Section 2.2: Dimer check on dimer_details

	if has_neighbouring_molecules:

		# Tenth, check all dimers contain two molecules that exist.
		for dimer_name, (mol1_name, mol2_name, UCV_i, UCV_j, UCV_k, DV_x, DV_y, DV_z, move_COM_x, move_COM_y, move_COM_z) in dimer_details.items():
			if not (1 <= mol1_name <= number_of_molecules):
				raise Exception('Error: Molecule '+str(mol1_name)+' does not exist. Should be between 1 and '+str(number_of_molecules))
			if not (1 <= mol2_name <= number_of_molecules):
				raise Exception('Error: Molecule '+str(mol2_name)+' does not exist. Should be between 1 and '+str(number_of_molecules))

		# Eleventh, check that the dimer names in dimer_details are in consecutive order from 0 to len(dimer_details)-1
		all_dimer_names = sorted([dimer_name for dimer_name, (mol1_name, mol2_name, UCV_i, UCV_j, UCV_k, DV_x, DV_y, DV_z, move_COM_x, move_COM_y, move_COM_z) in dimer_details.items()])
		if not all_dimer_names == sorted(range(1,len(dimer_details)+1)):
			to_string  = 'Error: The dimers names in dimer_details are not in consecutive order from 1 to len(dimer_details)\n'
			to_string += 'Dimer names in dimer_details: '+str(all_dimer_names)+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# ----------------------------------------------------------------------------------------------------------------------------------------------------------
	# Section 2.3: structurally_equivalent_molecule_pairs checks

	if has_unique_molecules:

		# Twelfth, check that structurally_equivalent_molecule_pairs is consistent. 

		# 12.1: Initialise a list to record all the sets of symmetric groups to.
		symmetric_groups = []

		# 12.2: For each set of symmetric dimer pairs.
		for dimer1_no, dimer2_no in structurally_equivalent_molecule_pairs:

			# 12.2.1: Determine how to put the dimer pair into the symmetric_groups list.
			found_symmetric_group = False
			for index, symmetric_group in enumerate(symmetric_groups):
				if (dimer1_no in symmetric_group) or (dimer2_no in symmetric_group):
					found_symmetric_group = True
					if (dimer1_no in symmetric_group) and (dimer2_no in symmetric_group):
						continue
					elif (dimer1_no in symmetric_group):
						symmetric_groups[index].append(dimer2_no)
					elif (dimer2_no in symmetric_group):
						symmetric_groups[index].append(dimer2_no)

			# 12.2.2: If symmetric_groups is not found, add it as a new set to symmetric_groups.
			if not found_symmetric_group:
				symmetric_groups.append([dimer1_no, dimer2_no])
		
		# 12.3: Check if all sets do not have repeating numbers.
		symmetric_groups_with_repeats = []
		for symmetric_group in symmetric_groups:
			if not len(symmetric_group) == len(set(symmetric_group)):
				symmetric_groups_with_repeats.append(symmetric_group)
		if len(symmetric_groups_with_repeats) > 0:
			to_string  = 'Error: One or more of the sets of structurally equivalent molecule groups has repeated molecules in them.\n'
			to_string += 'This is fine, but its weird, probably indicates a possible programming error\n'
			to_string += 'Structurally equivalent molecule groups with repeats: '+str(symmetric_groups_with_repeats)+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

		# 12.4: Convert the lists in symmetric_groups into sets:
		for index in range(len(symmetric_groups)):
			symmetric_groups[index] = set(symmetric_groups[index])

		# 12.5: Check if all dimers do not belong to multiple sets. 
		all_molecules_in_multiple_sets = []
		for index1 in range(len(symmetric_groups)):
			symmetric_group_of_focus = symmetric_groups[index1]
			for index2 in range(index1+1, len(symmetric_groups)):
				symmetric_group_to_compare = symmetric_groups[index2]
				# 12.5.1: Check if any dimer in symmetric_group_of_focus is found in symmetric_group_to_compare
				molecules_in_multiple_sets = symmetric_group_to_compare & symmetric_group_of_focus
				if not len(molecules_in_multiple_sets) == 0:
					all_molecules_in_multiple_sets.append((molecules_in_multiple_sets, symmetric_group_to_compare, symmetric_group_of_focus))
		if len(all_molecules_in_multiple_sets) > 0:
			to_string  = 'Error: There are one or more dimers that are found in more than one symmetric set. This should not happen and indicates a programming error in the main ECCP program.\n'
			to_string += 'Dimers in mutliple sets:\n'
			for molecules_in_multiple_sets, symmetric_group_to_compare, symmetric_group_of_focus in all_molecules_in_multiple_sets:
				to_string += str(molecules_in_multiple_sets)+' --> \t'+str(symmetric_group_to_compare)+'\t'+str(symmetric_group_of_focus)+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# ----------------------------------------------------------------------------------------------------------------------------------------------------------
	# Section 2.4: structurally_equivalent_dimer_pairs checks

	if has_unique_dimers:

		# Thirteenth, check that structurally_equivalent_dimer_pairs is consistent. 

		# 13.1: Initialise a list to record all the sets of symmetric groups to.
		symmetric_groups = []

		# 13.2: For each set of symmetric dimer pairs.
		for dimer1_no, dimer2_no in structurally_equivalent_dimer_pairs:

			# 13.2.1: Determine how to put the dimer pair into the symmetric_groups list.
			found_symmetric_group = False
			for index, symmetric_group in enumerate(symmetric_groups):
				if (dimer1_no in symmetric_group) or (dimer2_no in symmetric_group):
					found_symmetric_group = True
					if (dimer1_no in symmetric_group) and (dimer2_no in symmetric_group):
						continue
					elif (dimer1_no in symmetric_group):
						symmetric_groups[index].append(dimer2_no)
					elif (dimer2_no in symmetric_group):
						symmetric_groups[index].append(dimer1_no)

			# 13.2.2: If symmetric_groups is not found, add it as a new set to symmetric_groups.
			if not found_symmetric_group:
				symmetric_groups.append([dimer1_no, dimer2_no])
		
		# 13.3: Check if all sets do not have repeating numbers.
		symmetric_groups_with_repeats = []
		for symmetric_group in symmetric_groups:
			if not len(symmetric_group) == len(set(symmetric_group)):
				symmetric_groups_with_repeats.append(symmetric_group)
		if len(symmetric_groups_with_repeats) > 0:
			to_string  = 'Error: One or more of the sets of symmetric dimers has dimers repeated in them.\n'
			to_string += 'This is fine, but its weird, probably indicates a possible programming error\n'
			to_string += 'Symmetric groups with repeats: '+str(symmetric_groups_with_repeats)+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

		# 13.4: Convert the lists in symmetric_groups into sets:
		for index in range(len(symmetric_groups)):
			symmetric_groups[index] = set(symmetric_groups[index])

		# 13.5: Check if all dimers do not belong to multiple sets. 
		all_dimers_in_multiple_sets = []
		for index1 in range(len(symmetric_groups)):
			symmetric_group_of_focus = symmetric_groups[index1]
			for index2 in range(index1+1, len(symmetric_groups)):
				symmetric_group_to_compare = symmetric_groups[index2]
				# 13.5.1: Check if any dimer in symmetric_group_of_focus is found in symmetric_group_to_compare
				dimers_in_multiple_sets = symmetric_group_to_compare & symmetric_group_of_focus
				if not len(dimers_in_multiple_sets) == 0:
					all_dimers_in_multiple_sets.append((dimers_in_multiple_sets, symmetric_group_to_compare, symmetric_group_of_focus))
		if len(all_dimers_in_multiple_sets) > 0:
			to_string  = 'Error: There are one or more dimers that are found in more than one symmetric set. This should not happen and indicates a programming error in the main ECCP program.\n'
			to_string += 'Dimers in mutliple sets:\n'
			for dimers_in_multiple_sets, symmetric_group_to_compare, symmetric_group_of_focus in all_dimers_in_multiple_sets:
				to_string += str(dimers_in_multiple_sets)+' --> \t'+str(symmetric_group_to_compare)+'\t'+str(symmetric_group_of_focus)+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# ----------------------------------------------------------------------------------------------------------------------------------------------------------
	# Section 2.5: Solvent checks

	# Fourteenth, check that the dimers with solvents in them are consistent with which molecules are solvents. 
	if have_ECCP_Information_crystal_file:
		dimers_that_have_solvent_molecules = []
		for dimer_name, (mol1_name, mol2_name, UCV_i, UCV_j, UCV_k, DV_x, DV_y, DV_z, move_COM_x, move_COM_y, move_COM_z) in dimer_details.items():
			if (mol1_name in solvent_names) or (mol2_name in solvent_names):
				dimers_that_have_solvent_molecules.append(dimer_name) 
		if not sorted(dimers_that_contain_solvents) == sorted(dimers_that_have_solvent_molecules):
			to_string  = 'Error: There is inconsistency between the solvent molecules, and the dimers that have been recorded as having solvents in them.\n'
			to_string += 'Dimers that have solvent molecules as obtained using the "Molecules that are solvents" entry in ECCP_Information.txt: '+str(dimers_that_have_solvent_molecules)+'\n'
			to_string += 'Dimers that have solvent molecules as given by the "Dimers that contain solvents" entry in ECCP_Information.txt: '+str(dimers_that_contain_solvents)+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# ----------------------------------------------------------------------------------------------------------------------------------------------------------
	# Section 2.6: Checks on the crystal object

	if have_ECCP_Information_crystal_file:

		# Fifteenth, make networkx.Graph objects for each molecule using data from MoleculeList and NeighboursList.

		# 15.1: Get the nodes from MoleculeList.
		node_to_molecule = {}
		molecule_nodes = {}
		for atom_index, molecule_name in enumerate(MoleculeList):
			molecule_nodes.setdefault(molecule_name,[]).append(atom_index)
			if atom_index in node_to_molecule:
				raise Exception('Programming error?')
			node_to_molecule[atom_index] = molecule_name

		# 15.2: Get the edges from NeighboursList.
		molecule_edges_unassigned = []
		for atom_index, bonded_atom_indices in enumerate(NeighboursList):
			bonded_atom_indices = eval('['+str(bonded_atom_indices)+']')
			for bonded_atom_index in bonded_atom_indices:

				# 15.2.1: If atom_index < bonded_atom_index, record it.
				if atom_index < bonded_atom_index:
					if (atom_index, bonded_atom_index) in molecule_edges_unassigned:
						raise Exception('Error: '+str((atom_index, bonded_atom_index))+' already in molecule_edges_unassigned')
					molecule_edges_unassigned.append((atom_index,bonded_atom_index))
					continue

				# 15.2.2: If atom_index > bonded_atom_index, make sure (bonded_atom_index, atom_index) is already in molecule_edges_unassigned.
				if atom_index > bonded_atom_index:
					if not (bonded_atom_index, atom_index) in molecule_edges_unassigned:
						raise Exception('Error: '+str((bonded_atom_index, atom_index))+' is not in molecule_edges_unassigned. It should be.')
					continue

				# 15.2.3: If atom_index == bonded_atom_index, raise an error as atom is bonded to itself.
				if atom_index == bonded_atom_index:
					raise Exception('Error: Atom '+str(atom_index)+' is bonded to itself. That is a problem.')
					continue

		# 15.3: Assign edges in molecule_edges_unassigned to moelcules.
		molecule_edges = {}
		for atom1_index, atom2_index in molecule_edges_unassigned:
			mol1_name = node_to_molecule[atom1_index]
			mol2_name = node_to_molecule[atom2_index]
			if not mol1_name == mol2_name:
				to_string  = 'Error: Molecule '+str(mol1_name+1)+' and '+str(mol2_name+1)+' are bonded to eachother through bond '+str((atom1_index, atom2_index))+'. Check this.'
				raise Exception(to_string)
			molecule_edges.setdefault(mol1_name,[]).append((atom1_index, atom2_index))

		# 15.4: Make the dictionary for storing graph objects.
		molecule_graphs = {}
		for index in molecule_nodes.keys():
			molecule_graph = Graph()
			molecule_graph.add_nodes_from(molecule_nodes[index])
			molecule_graph.add_edges_from(molecule_edges[index])
			if not sorted(molecule_graph.nodes()) == sorted(molecule_nodes[index]):
				raise Exception('Error: There are edge atoms in molecule_edges['+str(index)+'] that are not in molecule_nodes['+str(index)+']. Check this')
			molecule_graphs[index] = molecule_graph

		# 15.5: Check that all graphs for each molecule are fully connected.
		for mol_name, molecule_graph in molecule_graphs.items():
			if not is_connected(molecule_graph):
				raise Exception('Error: molecule '+str(mol_name)+' is not fully connected.')

	# ----------------------------------------------------------------------------------------------------------------------------------------------------------



