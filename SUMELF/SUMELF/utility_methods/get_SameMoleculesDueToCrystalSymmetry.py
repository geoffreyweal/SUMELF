"""
get_SameMoleculesDueToCrystalSymmetry.py, Geoffrey Weal, 1/4/24

This method is designed to convert the inputfrom "SameMoleculesDueToCrystalSymmetry" from a crystal into dictionary form.
"""
from collections import Counter

def get_SameMoleculesDueToCrystalSymmetry(SameMoleculesDueToCrystalSymmetry_as_string):
	"""
	This method is designed to convert the inputfrom "SameMoleculesDueToCrystalSymmetry" from a crystal into dictionary form.

	Parameters
	----------
	SameMoleculesDueToCrystalSymmetry_as_string : str.
		This string gives the input for SameMoleculesDueToCrystalSymmetry as given directly from the crystal file in string form. 

	SameMoleculesDueToCrystalSymmetry : dict.
		This record all the molecules that are the same due to the symmetry of the crystal. This is given as {unique molecule: kist of structurally equivalent molecules due to the symmetry of the crystal}
	"""

	# First, if SameMoleculesDueToCrystalSymmetry_as_string starts with {, remove it from the string
	if SameMoleculesDueToCrystalSymmetry_as_string.startswith('{'):
		SameMoleculesDueToCrystalSymmetry_as_string = SameMoleculesDueToCrystalSymmetry_as_string[1:]

	# Second, if SameMoleculesDueToCrystalSymmetry_as_string ends with }, remove it from the string
	if SameMoleculesDueToCrystalSymmetry_as_string.endswith('}'):
		SameMoleculesDueToCrystalSymmetry_as_string = SameMoleculesDueToCrystalSymmetry_as_string[:-1]

	# Third, split up the list based on commas.
	SameMoleculesDueToCrystalSymmetry_comma_split = SameMoleculesDueToCrystalSymmetry_as_string.split(',')

	# Fourth, create a list that hold all the components of the SameMoleculesDueToCrystalSymmetry string
	SameMoleculesDueToCrystalSymmetry_components = []
	for variable in SameMoleculesDueToCrystalSymmetry_comma_split:
		if ':' in variable:
			SameMoleculesDueToCrystalSymmetry_components.append([])
		SameMoleculesDueToCrystalSymmetry_components[-1].append(variable)
	SameMoleculesDueToCrystalSymmetry_components = [','.join(component) for component in SameMoleculesDueToCrystalSymmetry_components]

	# Fifth, initialise the SameMoleculesDueToCrystalSymmetry dictionary
	SameMoleculesDueToCrystalSymmetry = {}

	# Sixth, record a list of all the molecules that have been recorded.
	molecules_recorded = []

	# Seventh, for each compoment in SameMoleculesDueToCrystalSymmetry_components:
	for component in SameMoleculesDueToCrystalSymmetry_components:

		# 7.1: Split the string based on the colon symbol, :.
		if component.endswith(':'):
			unique_mol_name              = int(component.replace(':',''))
			list_of_equivalent_mol_names = []
		else:
			unique_mol_name, list_of_equivalent_mol_names = component.split(':')
			unique_mol_name              = int(unique_mol_name)
			list_of_equivalent_mol_names = eval(list_of_equivalent_mol_names)

		# 7.2: Make sure that unique_mol_name is not already in SameMoleculesDueToCrystalSymmetry:
		if unique_mol_name in SameMoleculesDueToCrystalSymmetry:
			to_string += f'Error: Molecule {unique_mol_name} may have been entered twice or more in SameMoleculesDueToCrystalSymmetry\n'
			to_string += f'SameMoleculesDueToCrystalSymmetry_as_string: {SameMoleculesDueToCrystalSymmetry_as_string}\n'
			to_string += 'Check this.'
			raise Exception(to_string)

		# 7.3: Add unique_mol_name to SameMoleculesDueToCrystalSymmetry.
		SameMoleculesDueToCrystalSymmetry[unique_mol_name] = list_of_equivalent_mol_names

		# 7.4: Add the unique_mol_name group to molecules_recorded for double checking later on.
		molecules_recorded += [unique_mol_name] + list(list_of_equivalent_mol_names)

	# Eighth, check that there are no duplicates in molecules_recorded. 
	#          * If there are, there may be a problem
	if not len(molecules_recorded) == len(set(molecules_recorded)):
		to_string  = 'Error: There are duplicates molecule entries in SameMoleculesDueToCrystalSymmetry.\n'
		duplicates_in_molecules_recorded = sorted([mol_name for mol_name, count in Counter(molecules_recorded).items() if (count >= 2)])
		to_string += f'Duplicate Molecule Names: {duplicates_in_molecules_recorded}\n'
		to_string += f'SameMoleculesDueToCrystalSymmetry: {SameMoleculesDueToCrystalSymmetry}\n'
		to_string += 'Check this.'
		print(to_string)
		import pdb; pdb.set_trace()
		raise Exception(to_string)

	# Ninth, return SameMoleculesDueToCrystalSymmetry
	return SameMoleculesDueToCrystalSymmetry