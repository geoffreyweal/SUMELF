"""
run_which_ECCP_operations.py, Geoffrey Weal, 27/3/24

This method is designed to determine which componets of the ECCP need to be performed from scratch, and which can be obtained from file. 
"""

def run_which_ECCP_operations(have_crystal_file, have_ECCP_Information_file, have_All_Dimer_Information_file, have_Structurally_Equivalent_Molecule_Groups_file, have_Conformationally_Equivalent_Molecule_Groups_file, have_Structurally_Equivalent_Molecule_Pairs_file, have_Structurally_Equivalent_Dimer_Groups_file, have_Structurally_Equivalent_Dimer_Pairs_file):
	"""
	This method is designed to determine which componets of the ECCP need to be performed from scratch, and which can be obtained from file. 

	Parameters
	----------
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

	Returns
	-------
	have_ECCP_Information_file: bool.
		This indicates if you have the ECC_Information.txt file.
	have_crystal_file : bool.
		Indicates if the ECCP program has access to a crystal file of the system
	has_neighbouring_molecules : bool.
		Indicates if the ECCP program has a list of the molecules neighbouring molecules in the crystal on file. 
	has_unique_molecules : bool.
		Indicates if the ECCP program has the structurally and conformationally equivalent molecule groups on file. 
	has_unique_dimers : bool.
		Indicates if the ECCP program has the structurally equivalent dimer groups on file. 
	"""

	# First, initalise the booleans indicating which data the ECCP program has available to it. 
	has_neighbouring_molecules = False
	has_unique_molecules       = False
	has_unique_dimers          = False

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Second, determine if All_Dimer_Information.txt exist for determining the neighbours surrounding each molecule in the crystal. 
	if have_All_Dimer_Information_file:
		has_neighbouring_molecules = True

	# Third, determine if the files containing molecule equivalence groups information exist.
	if have_Structurally_Equivalent_Molecule_Groups_file or have_Conformationally_Equivalent_Molecule_Groups_file or have_Structurally_Equivalent_Molecule_Pairs_file:
		if not (have_Structurally_Equivalent_Molecule_Groups_file and have_Conformationally_Equivalent_Molecule_Groups_file and have_Structurally_Equivalent_Molecule_Pairs_file):
			to_string  = 'Error: Some of the files needed for retrieving molecule equivalence groups are missing:\n'
			to_string += '\n'
			to_string += f'Structurally_Equivalent_Molecule_Groups.txt: {have_Structurally_Equivalent_Molecule_Groups_file}\n'
			to_string += f'Conformationally_Equivalent_Molecule_Groups.txt: {have_Conformationally_Equivalent_Molecule_Groups_file}\n'
			to_string += f'Structurally_Equivalent_Molecule_Pairs.txt: {have_Structurally_Equivalent_Molecule_Pairs_file}\n'
			to_string += '\n'
			to_string += 'Recommendation: Remove all these files and run the ECCP again.\n'
			to_string += '\n'
			to_string += 'Check this.\n'
			raise Exception(to_string)
		has_unique_molecules = True

	# Fourth, determine if the files containing dimer equivalence groups information exist.
	if have_Structurally_Equivalent_Dimer_Groups_file or have_Structurally_Equivalent_Dimer_Pairs_file:
		if not (have_Structurally_Equivalent_Dimer_Groups_file and have_Structurally_Equivalent_Dimer_Pairs_file):
			to_string  = 'Error: Some of the files needed for retrieving dimer equivalence groups are missing:\n'
			to_string += '\n'
			to_string += f'Structurally_Equivalent_Dimer_Groups.txt: {have_Structurally_Equivalent_Dimer_Groups_file}\n'
			to_string += f'Structurally_Equivalent_Dimer_Pairs.txt: {have_Structurally_Equivalent_Dimer_Pairs_file}\n'
			to_string += '\n'
			to_string += 'Recommendation: Remove all these files and run the ECCP again.\n'
			to_string += '\n'
			to_string += 'Check this.\n'
			raise Exception(to_string)
		has_unique_dimers = True

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Fifth, if have_ECCP_Information_file is False, but one or more or the above are True, there is a problem and the previous ECCP run should be reset.
	if not have_ECCP_Information_file:
		if any([has_neighbouring_molecules, has_unique_molecules, has_unique_dimers]):
			to_string  = 'Error: The ECCP program can find some of the files used for inputs, but can not find the ECCP_Information.txt file.\n'
			to_string += 'The ECCP_Information.txt file is required if you want to use the other files in the ECCP program.\n'
			to_string += '\n'
			to_string += f'Found ECCP_Information.txt: {have_ECCP_Information_file}\n'
			to_string += '\n'
			to_string += 'Other input files found:\n'
			to_string += f'crystal.xyz: {have_crystal_file}\n'
			to_string += f'All_Dimer_Information.txt: {have_All_Dimer_Information_file}\n'
			to_string += '\n'
			to_string += f'Structurally_Equivalent_Molecule_Groups.txt: {have_Structurally_Equivalent_Molecule_Groups_file}\n'
			to_string += f'Conformationally_Equivalent_Molecule_Groups.txt: {have_Conformationally_Equivalent_Molecule_Groups_file}\n'
			to_string += f'Structurally_Equivalent_Molecule_Pairs.txt: {have_Structurally_Equivalent_Molecule_Pairs_file}\n'
			to_string += '\n'
			to_string += f'Structurally_Equivalent_Dimer_Groups.txt: {have_Structurally_Equivalent_Dimer_Groups_file}\n'
			to_string += f'Structurally_Equivalent_Dimer_Pairs.txt: {have_Structurally_Equivalent_Dimer_Pairs_file}\n'
			to_string += '\n'
			to_string += 'Recommendation: Reset your previous ECCP run and run again.\n'
			to_string += '\n'
			to_string += 'Check this.\n'
			raise Exception(to_string)

	# Sixth, if have_ECCP_Information_file is False, but one or more or the above are True, there is a problem and the previous ECCP run should be reset.
	if not have_crystal_file:
		if any([have_ECCP_Information_file, has_neighbouring_molecules, has_unique_molecules, has_unique_dimers]):
			to_string  = 'Error: The ECCP program can find some of the files used for inputs, but can not find the ECCP_Information.txt file.\n'
			to_string += 'The crystal.xyz file is required if you want to use the other files in the ECCP program.\n'
			to_string += '\n'
			to_string += f'Found crystal.xyz: {have_crystal_file}\n'
			to_string += '\n'
			to_string += 'Other input files found:\n'
			to_string += f'ECCP_Information.txt: {have_ECCP_Information_file}\n'
			to_string += f'All_Dimer_Information.txt: {have_All_Dimer_Information_file}\n'
			to_string += '\n'
			to_string += f'Structurally_Equivalent_Molecule_Groups.txt: {have_Structurally_Equivalent_Molecule_Groups_file}\n'
			to_string += f'Conformationally_Equivalent_Molecule_Groups.txt: {have_Conformationally_Equivalent_Molecule_Groups_file}\n'
			to_string += f'Structurally_Equivalent_Molecule_Pairs.txt: {have_Structurally_Equivalent_Molecule_Pairs_file}\n'
			to_string += '\n'
			to_string += f'Structurally_Equivalent_Dimer_Groups.txt: {have_Structurally_Equivalent_Dimer_Groups_file}\n'
			to_string += f'Structurally_Equivalent_Dimer_Pairs.txt: {have_Structurally_Equivalent_Dimer_Pairs_file}\n'
			to_string += '\n'
			to_string += 'Recommendation: Reset your previous ECCP run and run again.\n'
			to_string += '\n'
			to_string += 'Check this.\n'
			raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Seventh, make sure that if All_Dimer_Information.txt exists that the crystal.xyz file also exists.
	if has_neighbouring_molecules and not have_crystal_file:
		to_string  = 'Error: The ECCP program found the All_Dimer_Information.txt, but not the crystal.xyz.\n'
		to_string  = 'crystal.xyz is needed if you want to recover dimers from All_Dimer_Information.txt.\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Eighth, if you have found information about the molecule equivalence groups, make sure that All_Dimer_Information.txt and crystal.xyz exists.
	if has_unique_molecules and any([not have_crystal_file, not has_neighbouring_molecules]):
		to_string  = 'Error: The ECCP program found files for obtaining unique molecules, but the following files were not found (which should have been):\n'
		to_string += '\n'
		to_string += f'crystal.xyz: {have_crystal_file}\n'
		to_string += f'All_Dimer_Information.txt: {have_All_Dimer_Information_file}\n'
		to_string += '\n'
		to_string += f'Structurally_Equivalent_Molecule_Groups.txt: {have_Structurally_Equivalent_Molecule_Groups_file}\n'
		to_string += f'Conformationally_Equivalent_Molecule_Groups.txt: {have_Conformationally_Equivalent_Molecule_Groups_file}\n'
		to_string += f'Structurally_Equivalent_Molecule_Pairs.txt: {have_Structurally_Equivalent_Molecule_Pairs_file}\n'
		to_string += '\n'
		to_string += 'Recommendation: Remove all files above (except for crystal.xyz) and run again.\n'
		to_string += '\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Ninth, if you have found information about the dimer equivalence groups, make sure that All_Dimer_Information.txt and crystal.xyz exists.
	if has_unique_dimers and any([not have_crystal_file, not has_neighbouring_molecules]):
		to_string  = 'Error: The ECCP program found files for obtaining unique dimers, but the following files were not found (which should have been):\n'
		to_string += '\n'
		to_string += f'crystal.xyz: {have_crystal_file}\n'
		to_string += f'All_Dimer_Information.txt: {have_All_Dimer_Information_file}\n'
		to_string += '\n'
		to_string += f'Structurally_Equivalent_Dimer_Groups.txt: {have_Structurally_Equivalent_Dimer_Groups_file}\n'
		to_string += f'Structurally_Equivalent_Dimer_Pairs.txt: {have_Structurally_Equivalent_Dimer_Pairs_file}\n'
		to_string += '\n'
		to_string += 'Recommendation: Remove all files above (except for crystal.xyz) and run again.\n'
		to_string += '\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Ninth, return has_neighbouring_molecules, has_unique_molecules, and has_unique_dimers
	return has_neighbouring_molecules, has_unique_molecules, has_unique_dimers






























