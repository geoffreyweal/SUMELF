"""
get_dimer_details_data.py, Geoffrey Weal, 27/3/2024

This method is designed to collect all the information needed from the All_Dimer_Information.txt file
"""
import numpy as np

check_header_components        = ['Dimer', 'No', 'First', 'Molecule', 'Second', 'Molecule', 'UCV(i)', 'UCV(j)', 'UCV(k)', 'DV(x)', 'DV(y)', 'DV(z)', 'move_dimer_COM(x)', 'move_dimer_COM(y)', 'move_dimer_COM(z)']
check_header_components_string = ['Dimer No',    'First Molecule',    'Second Molecule ',   'UCV(i)', 'UCV(j)', 'UCV(k)', 'DV(x)', 'DV(y)', 'DV(z)', 'move_dimer_COM(x)', 'move_dimer_COM(y)', 'move_dimer_COM(z)']
def get_dimer_details_data(path_to_All_Dimer_Information_folder, eccp_information):
	"""
	This method is designed to collect all the information needed from the All_Dimer_Information.txt file.

	Parameters
	----------
	path_to_All_Dimer_Information_folder : str.
		This is the path to the All_Dimer_Information text file. 
	eccp_information : dict.
		This dictionary contains all the data from the previous ECCP run. Included in this is the number of dimers that are expected. Used for double-checking purposes.

	Returns
	-------
	dimers_details : dict
		This dictionary contains information about the dimers. 
	"""

	# First, obtain the number of dimers as given in the eccp_information file. 
	if 'number_of_dimers' in eccp_information:
		number_of_dimers = eccp_information['number_of_dimers']
	else:
		number_of_dimers = None

	# Second, initialise a dictionary to hold which symmetric molecules go with which unique molecule in the "path_to_Conformationally_Unique_Molecule_Information_file" file.
	dimers_details = {}

	# Third, write the full path to the All_Dimer_Information.txt file. 
	path_to_All_Dimer_Information_file = path_to_All_Dimer_Information_folder+'/'+'All_Dimer_Information.txt'

	# Fourth, collect the data from the Structurally_Unique_Molecule_Information text file. 
	with open(path_to_All_Dimer_Information_file, 'r') as All_Dimer_InformationTXT:
		
		# 4.1: Check that the first line just a title line. 
		line = All_Dimer_InformationTXT.readline()
		line = line.rstrip()
		if not (line.split() == check_header_components):
			to_string  = 'Error: First line should be "'+'\t'.join(check_header_components_string)+'"\n'
			to_string += f'However, it is: {line}\n'
			to_string += 'Check this'
			print(to_string)
			import pdb; pdb.set_trace()
			raise Exception(to_string)

		# 4.2: Collect the list of symmetric molecule: unique molecule from the rest of the "path_to_Conformationally_Unique_Molecule_Information_file" file.
		for line in All_Dimer_InformationTXT:
			line = line.rstrip()

			# 4.2.1: Collect the information about the dimer from this line.
			dimer_no, mol1_name, mol2_name, UCV_i, UCV_j, UCV_k, DV_x, DV_y, DV_z, move_COM_x, move_COM_y, move_COM_z = line.split()

			# 4.2.2: Convert the information from the line into ints and floats.
			dimer_no   =   int(dimer_no)
			mol1_name  =   int(mol1_name)
			mol2_name  =   int(mol2_name)
			UCV_i      =   int(UCV_i)
			UCV_j      =   int(UCV_j)
			UCV_k      =   int(UCV_k)
			DV_x       = float(DV_x)
			DV_y       = float(DV_y)
			DV_z       = float(DV_z)
			move_COM_x = float(move_COM_x)
			move_COM_y = float(move_COM_y)
			move_COM_z = float(move_COM_z)

			# 4.2.3: Make sure that dimer_no has not been entered more than once in the dimers_details dictionary.
			if dimer_no in dimers_details.keys():
				raise Exception('Error: '+str(dimer_no)+' already in dimers_details. Make sure that '+str(dimer_no)+' has not been entered into '+str(path_to_All_Dimer_Information_file)+' more than once.')

			# 4.2.4: Enter information about the dimer into the dimers_details dictionary. 
			dimers_details[dimer_no] = (mol1_name, mol2_name, UCV_i, UCV_j, UCV_k, DV_x, DV_y, DV_z, move_COM_x, move_COM_y, move_COM_z)

	# Fifth, make sure that the dimers entered are consecutive. If not, this indicates there is an issue.
	if number_of_dimers is not None:
		expected_dimers = list(range(1,number_of_dimers+1))
	else:
		expected_dimers = list(range(1,len(dimers_details)+1))
	if not sorted(dimers_details.keys()) == expected_dimers:
		to_string  = 'Issue: The names of dimers entered into '+str(path_to_All_Dimer_Information_file)+' are not in consecutive order. This is indicative of an issue\n'
		to_string += 'Expected Dimers: '+str(expected_dimers)+'\n'
		to_string += 'Dimer names in file: '+str(sorted(dimers_details.keys()))+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Sixth, obtain all the dimer indices in dimers_details
	dimer_names = sorted(dimers_details.keys())

	# Seventh, check that dimer_names does not contain any duplicates
	if not len(dimer_names) == len(set(dimer_names)):
		to_string  = 'Issue: There are duplicate dimer names in dimers_details\n'
		duplicates_in_dimers_details_with_dimer_indices = sorted([dimer_index for dimer_index, count in Counter(dimers_details).items() if (count >= 2)])
		to_string += 'Duplicates in dimers_details: '+str(duplicates_in_dimers_details_with_dimer_indices)+'\n'
		to_string += 'Dimer names recorded in dimers_details: '+str(dimer_names)+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Eighth, return dimers_details
	return dimers_details

