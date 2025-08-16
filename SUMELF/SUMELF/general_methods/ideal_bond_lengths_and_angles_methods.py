"""
ideal_bond_lengths_and_angles_methods.py, Geoffrey Weal, 19/7/22

This script is designed to hold methods useful for determining the ideal bond lengths and angles between atoms in a molecule.
"""
from math import pi

def create_bond_lengths():
	"""
	This method will create a double dictionary that contains all the bond lengths typical between two elements.

	Returns
	-------
	bond_lengths : dict.
		The bond lengths between elements.
	"""

	# First, initialise the overal bond lengths dictionary.
	bond_lengths = {element: {} for element in ['H', 'C', 'O']}

	# Second, create all the bond distances between elements and hydrogen.
	bond_lengths['H']['O'] = {0: None, 1: 0.97, 2: 0.96, 3: 0.96}
	bond_lengths['H']['S'] = {0: None, 1: None, 2: 1.34, 3: 1.41}
	#bond_lengths['H']['N'] = {0: None, 1: None, 2: 1.02, 3: 1.02, 4: 1.03}
	bond_lengths['H']['N'] = {0: None, 1: None, 2: 0.97, 3: 0.97, 4: 0.97}
	#bond_lengths['H']['C'] = {0: None, 1: None, 2: 1.06, 3: 1.08, 4: 1.10}
	bond_lengths['H']['C'] = {0: None, 1: None, 2: 0.97, 3: 0.97, 4: 0.97}
	bond_lengths['H']['B'] = {0: None, 1: None, 2: None, 3: 1.19, 4: 1.17}

	# Third, create all the bond distances between elements and carbon.
	bond_lengths['C']['O'] = {0: None, 1: None, 2: 1.43, 3: None}
	bond_lengths['C']['S'] = {0: None, 1: None, 2: 1.50, 3: None}
	bond_lengths['C']['N'] = {0: None, 1: None, 2: 1.28, 3: 1.47, 4: None}
	bond_lengths['C']['C'] = {0: None, 1: None, 2: 1.20, 3: 1.34, 4: 1.53}
	bond_lengths['C']['B'] = {0: None, 1: None, 2: None, 3: None, 4: None}

	# Fourth, create all the bond distances between elements and oxygen.
	bond_lengths['O']['O'] = {0: None, 1: None, 2: None, 3: None}
	bond_lengths['O']['S'] = {0: None, 1: None, 2: 1.43, 3: None}
	bond_lengths['O']['N'] = {0: None, 1: None, 2: 1.40, 3: None, 4: None}
	bond_lengths['O']['C'] = bond_lengths['C']['O']
	bond_lengths['O']['B'] = {0: None, 1: None, 2: None, 3: None, 4: None}

	# Fifth, return the bond_lengths dictionary
	return bond_lengths

bond_lengths = create_bond_lengths()
def get_bond_lengths(element_to_attach=None):
	"""
	This method is designed to return the bond lengths for atoms between an atom, and the atom to attach.

	Parameters
	----------
	element_to_attach : str.
		This is the element you are wanting to add to your molecule/crystal/system.

	Returns
	-------
	Return a dictionary of bondlengths that include the element_to_attach in it.
	"""

	# First, if element_to_attach is None, this is a problem and we need to sort out what is happening in the code. 
	if element_to_attach is None:
		raise Exception('element_to_attach needs to be some element. element_to_attach = '+str(element_to_attach))

	# Second, if element_to_attach is not included in bond_lengths, raise this as this will need to be programmed into the program.
	if element_to_attach not in bond_lengths.keys():
		raise Exception('element_to_attach needs to be either H, C, or O. element_to_attach = '+str(element_to_attach))

	# Third, return the bond lengths for element_to_attach.
	return bond_lengths[element_to_attach]

# ---------------------------------------------------------------------------------------------------------------------------

def get_bond_angles_to_element():
	"""
	This method is designed to provide all the bond angle information for various elements with various number of neighbours bonded to it.

	The bond_angles dictionary is given by:

	* element_to_attach                       : The atom that will be attached.
	* element                                 : The atom that we are attaching the atom to.
	* total_no_of_neighbours                  : The total number of neighbouring atoms and lone pair of electons that need to be bound to element.
	* no_of_pairs_of_lone_electrons (optional): This is the number of lone pair of electrons contained about element.

	# The angle return is that for A-E-A, where E = element and A = element_to_attach

	Returns
	-------
	bond_angles : dict. 
		This dictionary contains all the bond angle information for various elements with various number of neighbours bonded to it. 
	"""

	# First, initialise the overall bond angles dictionary.
	bond_angles = {element: {} for element in ['H', 'C', 'O']}

	# Second, create all the bond angles between elements and hydrogen.
	bond_angles['H']['O'] = {0: None, 1: None, 2: 104.5, 3: 113.0}
	bond_angles['H']['S'] = {0: None, 1: None, 2: 92.1 , 3: 106.0}
	bond_angles['H']['N'] = {0: None, 1: None, 2: {0: 'flat', 1: 'flat', 2: 110.5}, 3: {0: 'flat', 1: 107.8}, 4: 109.5}
	bond_angles['H']['C'] = {0: None, 1: None, 2: {0: 'flat', 1: 'flat', 2: 109.5}, 3: {0: 120.0, 1: 109.5}, 4: 109.5}
	bond_angles['H']['B'] = {0: None, 1: None, 2: 'flat', 3: 'flat', 4: 109.5}

	# Third, create all the bond angles between elements and carbons.
	bond_angles['C']['O'] = {0: None, 1: None, 2: 113.0, 3: None}
	bond_angles['C']['S'] = {0: None, 1: None, 2: 97.7 , 3: None}
	bond_angles['C']['N'] = {0: None, 1: None, 2: {0: 'flat', 1: 'flat', 2: 110.5}, 3: {0: 'flat', 1: 107.0}, 4: 109.5}
	bond_angles['C']['C'] = {0: None, 1: None, 2: {0: 'flat', 1: 'flat', 2: 109.5}, 3: {0: 120.0, 1: 109.5}, 4: 109.5}
	bond_angles['C']['B'] = {0: None, 1: None, 2: 'flat', 3: 'flat', 4: None}

	# Fourth, create all the bond angles between elements and oxygen.
	bond_angles['O']['O'] = {0: None, 1: None, 2: None , 3: None}
	bond_angles['O']['S'] = {0: None, 1: None, 2: 119.5, 3: None}
	bond_angles['O']['N'] = {0: None, 1: None, 2: {0: 'flat', 1: 'flat', 2: 102.2}, 3: {0: None, 1: None}, 4: None}
	bond_angles['O']['C'] = {0: None, 1: None, 2: {0: 'flat', 1: 'flat', 2: 108.9}, 3: {0: None, 1: None}, 4: None}
	bond_angles['O']['B'] = {0: None, 1: None, 2: 'flat', 3: 'flat', 4: None}

	# Fifth, convert all angles from degrees to radians.
	for element_to_attach in bond_angles.keys():
		for element in bond_angles[element_to_attach].keys():
			for total_number_of_atom_surrounding_central_atom in bond_angles[element_to_attach][element].keys():
				if bond_angles[element_to_attach][element][total_number_of_atom_surrounding_central_atom] is None: 
					continue
				elif isinstance(bond_angles[element_to_attach][element][total_number_of_atom_surrounding_central_atom], float):
					bond_angles[element_to_attach][element][total_number_of_atom_surrounding_central_atom] *= (pi/180.0)
				elif isinstance(bond_angles[element_to_attach][element][total_number_of_atom_surrounding_central_atom], dict):
					for no_of_pairs_of_lone_electrons in bond_angles[element_to_attach][element][total_number_of_atom_surrounding_central_atom].keys():
						if isinstance(bond_angles[element_to_attach][element][total_number_of_atom_surrounding_central_atom][no_of_pairs_of_lone_electrons], float):
							bond_angles[element_to_attach][element][total_number_of_atom_surrounding_central_atom][no_of_pairs_of_lone_electrons] *= (pi/180.0)

	# Sixth, return the bond_angles dictionary
	return bond_angles

bond_angles = get_bond_angles_to_element()
def get_bond_angle(element, total_no_of_neighbours, no_of_pairs_of_lone_electrons=None, element_to_attach=None):
	"""
	This method is designed to return the desired bond angle between atoms surrounding the central atom.

	Parameters
	----------
	element : str.
		This is the element of the central atom
	total_no_of_neighbours : int.
		This is the number of atoms a
	no_of_pairs_of_lone_electrons : int.
		This is the number of lone pair of electrons in the central atom
	element_to_attach : str.
		This is the element you are wanting to add to your molecule/crystal/system.

	Returns
	-------
	The angle expected between bonds attach to the central atom of a specific element.
	"""

	# First, check that element_to_attach is not None, otherwise there is a script that does not provide the element_to_attach
	if element_to_attach is None:
		raise Exception('element_to_attach needs to be some element. element_to_attach = '+str(element_to_attach))

	# Second, if element_to_attach is not included in bond_lengths, raise this as this will need to be programmed into the program.
	if element_to_attach not in bond_angles.keys():
		raise Exception('element_to_attach needs to be either H, C, or O. element_to_attach = '+str(element_to_attach))
	
	# Third, obtain the bond angle between atoms. This is described by
	specific_bond_angles = bond_angles[element_to_attach][element][total_no_of_neighbours]

	# Fourth, do something if no_of_pairs_of_lone_electrons is not None.
	if no_of_pairs_of_lone_electrons is None:
		if not isinstance(specific_bond_angles,float):
			raise Exception('huh')
		return specific_bond_angles

	# Fifth, return the bond angle.
	if isinstance(specific_bond_angles,float):
		return specific_bond_angles
	else:
		return specific_bond_angles[no_of_pairs_of_lone_electrons]

# ---------------------------------------------------------------------------------------------------------------------------




