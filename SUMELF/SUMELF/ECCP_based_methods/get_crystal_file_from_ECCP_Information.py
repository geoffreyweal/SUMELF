
from SUMELF import read_crystal

def get_crystal_file_from_ECCP_Information(path_to_ECCP_Information_folder):
	"""
	This method is designed to collect all the information needed from the All_Dimer_Information.txt file.

	Parameters
	----------
	path_to_ECCP_Information_folder : str.
		This is the path to the crystal.xyz file in the ECCP_Information folder. 

	Returns
	-------
	crystal : ase.Atoms
		This is the ase.Atoms object for the crystal. 
	"""

	# First, write the full path to the crystal.xyz file in the ECCP_Information folder.
	path_to_crystal_xyz_file_in_ECCP_Information = path_to_ECCP_Information_folder+'/'+'crystal.xyz'

	# Second, obtain the crystal from file. 
	crystal = read_crystal(path_to_crystal_xyz_file_in_ECCP_Information)

	# Third, check that crystal file contains all information about the molecules in the crystal, including
	# * MoleculeList: Indicates which atoms go with which molecule in the All_Molecules folder
	# * NeighboursList: Indicates which atoms are bonded to which atoms in the molecule. 
	check_for_missing_attributes(crystal)

	# Fourth, return crystal
	return crystal


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
