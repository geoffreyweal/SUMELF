"""
get_ECCP_Information_data.py, Geoffrey Weal, 9/3/2024

This method is designed to collect all the information needed from the ECCP_Information.txt file
"""
import os
from types import SimpleNamespace
#from SUMELF.SUMELF.ECCP_based_methods.utility_methods.convert_list_from_names_to_indices import convert_list_from_names_to_indices

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# variables_read_signals is given as {line_prefix: [variable name, prefix for variable, suffix for variable, values that if given convert to another value]}
variables_read_signals = {}
variables_read_signals['make_dimer_method']                                                     = ('make_dimer_method', '', '', [])
variables_read_signals['environment_settings']                                                  = ('environment_settings', '', '', [])

variables_read_signals['Number of molecules obtained from crystal']                             = ('number_of_molecules', '', '', [])
variables_read_signals['Molecules that are solvents']                                           = ('solvent_names', '[', ']', [('None', '')])

variables_read_signals['Number of structurally unique molecules obtained from crystal']         = ('number_of_structurally_unique_molecules', '', '', [])
variables_read_signals['Number of structurally equivalent molecules obtained from crystal']     = ('number_of_structurally_equivalent_molecules', '', '', [])

variables_read_signals['Number of conformationally unique molecules obtained from crystal']     = ('number_of_conformationally_unique_molecules', '', '', [])
variables_read_signals['Number of conformationally equivalent molecules obtained from crystal'] = ('number_of_conformationally_equivalent_molecules', '', '', [])

variables_read_signals['Number of dimers obtained from crystal']                                = ('number_of_dimers', '', '', [])
variables_read_signals['Dimers that contain solvents']                                          = ('dimers_that_contain_solvents', '[', ']', [('None', '')])

variables_read_signals['Number of structurally unique dimers obtained from crystal']            = ('number_of_structurally_unique_dimers', '', '', [])
variables_read_signals['Number of structurally equivalent dimers obtained from crystal']        = ('number_of_structurally_equivalent_dimers', '', '', [])
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def get_ECCP_Information_data(path_to_ECCP_Information_folder):
	"""
	This method is designed to collect all the information needed from the ECCP_Information.txt file.

	Parameters
	----------
	path_to_ECCP_Information_folder : str.
		This is the path to the ECCP_Information folder. 

	Returns
	-------
	eccp_information_data : dict. 
		This contains all the information from the ECCP_Information.txt file. The data in this dictionary are:

			make_dimer_method : dict. 
				This dictionary indices the settings for how dimers were created by the Run_ECCP.py script
			environment_settings : dict.
				This dictionary indicates the environmental settings that indicates how molecules surrounding moelcules and dimers are included in calculations. 

			number_of_molecules : int
				This is the number of molecules in the crystal. 
			solvent_indices : list of ints
				This is the names of the solvent molecules in the crystal.

			number_of_structurally_unique_molecules : int
				This is the number of structurally unique molecules.
			number_of_structurally_equivalent_molecules : int
				This is the number of structurally equivalent molecules.

			number_of_conformationally_unique_molecules : int
				This is the number of conformationally unique molecules. 
			number_of_conformationally_equivalent_molecules : int
				This is the number of conformationally equivalent molecules. 

			number_of_dimers : int
				This is the number of dimers found for the given value of 'max_dimer_distance' in the make_dimer_method dictionary. 
			dimers_that_contain_solvents : list of int
				These are the dimers that contain solvent molecules in them. 

			number_of_structurally_unique_dimers : int
				This is the number of structurally unique dimers.
			number_of_conformationally_unique_dimers : int
				This is the number of structurally equivalent dimers.
	"""

	# First, write the full path to the ECCP_Information.txt file. 
	path_to_ECCP_Information_file = path_to_ECCP_Information_folder +'/'+'ECCP_Information.txt'

	# Second, check if the ECCP_Information.txt file exists. If it does not, there is an issue
	if not os.path.exists(path_to_ECCP_Information_file):
		to_string  = 'Error: The ECCP_Information.txt found was not found in '+str(path_to_ECCP_Information_file)+'. Check this.'
		raise Exception(to_string)

	# Third, read through the 'ECCP_Information.txt' file and obtain all the information from it.
	namespace = {}
	with open(path_to_ECCP_Information_file, 'r') as ECCP_InformationTXT:

		# 3.1: For each line in ECCP_InformationTXT
		for line in ECCP_InformationTXT:

			# 3.2: Remove the endline from line
			line = line.rstrip()

			# 3.3: For each variable to get as indicated by variables_read_signals
			for line_prefix, (variable_name, prefix, suffix, convert_inputs) in variables_read_signals.items():

				# 3.4: If line begins with line_prefix, we want to get the variable from this line
				if line.startswith(line_prefix):

					# 3.5: Get variable from ECCP_InformationTXT
					variable = str(line.replace(line_prefix+':','').lstrip())

					# 3.6: Convert the value to another based on variables_read_signals[-1]
					for input_for_variable, convert_to in convert_inputs:
						if variable == input_for_variable:
							variable = convert_to
							break

					# 3.7: Add variable to namespace
					variable = prefix + variable + suffix
					try:
						namespace[variable_name] = eval(variable)
					except Exception as expection:
						raise Exception('Can not obtain the variable for '+str(variable_name)+'\nvariable: '+str(variable)+'\nCheck this.')

	# Fourth, convert the namespace from a dictionary to SimpleNamespace object.
	namespace = SimpleNamespace(**namespace)

	# Fifth, obtain the compoents from namespace
	make_dimer_method                               = namespace.make_dimer_method
	environment_settings                            = namespace.environment_settings
	number_of_molecules                             = namespace.number_of_molecules
	solvent_names                                   = namespace.solvent_names
	number_of_structurally_unique_molecules         = namespace.number_of_structurally_unique_molecules
	number_of_structurally_equivalent_molecules     = namespace.number_of_structurally_equivalent_molecules
	number_of_conformationally_unique_molecules     = namespace.number_of_conformationally_unique_molecules
	number_of_conformationally_equivalent_molecules = namespace.number_of_conformationally_equivalent_molecules
	number_of_dimers                                = namespace.number_of_dimers
	dimers_that_contain_solvents                    = namespace.dimers_that_contain_solvents
	number_of_structurally_unique_dimers            = namespace.number_of_structurally_unique_dimers
	number_of_structurally_equivalent_dimers        = namespace.number_of_structurally_equivalent_dimers

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Sixth, convert the solvent lists for molecules and dimers from names to indices.
	#solvent_indices              = convert_list_from_names_to_indices(solvent_names)
	#dimers_that_contain_solvents = convert_list_from_names_to_indices(dimers_that_contain_solvents)

	# Seventh, remove solvent_names from namespace and replace it with solvent_indices
	#del namespace.solvent_names
	#namespace.solvent_indices = solvent_indices

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Eighth, perform checks to make sure that the ECCP_Information.txt file is consistent with itself.

	# 8.1: Check that the molecule numbers add up correctly for the structural equivalence molecule groups.
	if not number_of_molecules == number_of_structurally_unique_molecules + number_of_structurally_equivalent_molecules:
		to_string  = 'Error: The total number of molecules does not equal the number of structurally unique and equivalent molecules.\n'
		to_string += 'number_of_molecules: '+str(number_of_molecules)+'\n'
		to_string += 'number_of_structurally_unique_molecules: '+str(number_of_structurally_unique_molecules)+'\n'
		to_string += 'number_of_structurally_equivalent_molecules: '+str(number_of_structurally_equivalent_molecules)+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 8.2: Check that the molecule numbers add up correctly for the conformational equivalence molecule groups.
	if not number_of_molecules == number_of_conformationally_unique_molecules + number_of_conformationally_equivalent_molecules:
		to_string  = 'Error: The total number of molecules does not equal the number of conformationally unique and equivalent molecules.\n'
		to_string += 'number_of_molecules: '+str(number_of_molecules)+'\n'
		to_string += 'number_of_conformationally_unique_molecules: '+str(number_of_conformationally_unique_molecules)+'\n'
		to_string += 'number_of_conformationally_equivalent_molecules: '+str(number_of_conformationally_equivalent_molecules)+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 8.3: Check that the dimer numbers add up correctly for the structural equivalence dimer groups. 
	if not number_of_dimers == number_of_structurally_unique_dimers + number_of_structurally_equivalent_dimers:
		to_string  = 'Error: The total number of dimers does not equal the number of structurally unique and equivalent dimers.\n'
		to_string += 'number_of_dimers: '+str(number_of_dimers)+'\n'
		to_string += 'number_of_structurally_unique_dimers: '+str(number_of_structurally_unique_dimers)+'\n'
		to_string += 'number_of_structurally_equivalent_dimers: '+str(number_of_structurally_equivalent_dimers)+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 8.4: Check that the name of the solvents in solvent_names is not outside of 1 to number_of_molecules
	irregular_solvent_indices = [solvent_name for solvent_name in solvent_names if not (1 <= solvent_name <= number_of_molecules)]
	if len(irregular_solvent_indices) > 0:
		to_string  = 'Error: There are solvents that have indices outside of that expected.\n'
		to_string += 'irregular_solvent_indices: '+str(irregular_solvent_indices)+'\n'
		to_string += f'Solvent indices are expected to be between 1 and {number_of_molecules}.\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 8.5: Check that the name of the solvents in the dimers is not outside of 0 to number_of_molecules
	irregular_dimer_solvent_indices = [dimer_name for dimer_name in dimers_that_contain_solvents if not (1 <= dimer_name <= number_of_dimers)]
	if len(irregular_dimer_solvent_indices) > 0:
		to_string  = 'Error: There are dimers containing solvents that have indices outside of that expected.\n'
		to_string += 'irregular_dimer_solvent_indices: '+str(irregular_dimer_solvent_indices)+'\n'
		to_string += f'Solvent indices are expected to be between 1 and {number_of_dimers}.\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Eighth, convert namespace into dictionary
	eccp_information_data = vars(namespace)

	# Ninth, return information from ECCP_Information.txt
	return eccp_information_data






















