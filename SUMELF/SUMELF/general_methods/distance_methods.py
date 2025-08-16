"""
distance_methods.py, Geoffrey Weal, 19/7/22

This script is designed to hold methods useful for obtaining distances between points and bonds lengths, as well as determining if two atoms are likely bonded together. 
"""

from ase.data import atomic_numbers, covalent_radii

def get_distance(position1, position2):
	"""
	This method gives the bond length between two atoms in the molecule. Does not use Minimum Image Convention. 

	Parameters
	----------
	position1 : np.array
		A 1x3 row matrix containing the x, y, and z coordinates of the atom1. Given in Å. 
	position2 : np.array
		A 1x3 row matrix containing the x, y, and z coordinates of the atom2. Given in Å. 

	Returns
	-------
	float
		Returns the distance between those atom1 and atom2. Given in Å. 
		
	"""
	#return (sum([(float(p1)-float(p2))**2.0 for p1, p2 in zip(position1,position2)]))**0.5
	return (sum([(p1-p2)**2.0 for p1, p2 in zip(position1,position2)]))**0.5

def get_xyz_distances(position1, position2, give_absolute_distances=False):
	"""
	This method gives the x, y, and z distances between two positions. Does not use Minimum Image Convention. 

	Parameters
	----------
	position1 : np.array
		A 1x3 row matrix containing the x, y, and z coordinates of the atom1. Given in Å. 
	position2 : np.array
		A 1x3 row matrix containing the x, y, and z coordinates of the atom2. Given in Å. 
	give_absolute_distances : bool.
		If True, give the absolute distance between the each position (no negative values given). If false, give the relative distance between the each position (positive and negative values possible).

	Returns
	-------
	list of three floats
		Returns the x, y, and z distances between those atom1 and atom2. Given in Å. 
		
	"""

	# First, check that the lengths of position1 and position2 are the same.
	if not (len(position1) == len(position2)):
		raise Exception('Error, position lengths are not the same.\nposition1 = '+str(position1)+'\nposition1 length = '+str(len(position1))+'\nposition2 = '+str(position2)+'\nposition2 length = '+str(len(position2)))
	
	# Second, obtain the difference in distance for each coordinate
	differences_in_distance = []
	for atom_index in range(len(position1)):

		# 2.1: Get the distance for a coordinate between position1 and position2
		#distance = float(position1[atom_index]) - float(position2[atom_index])
		distance = position1[atom_index] - position2[atom_index]

		# 2.2: Obtain the absolute position if desired
		if give_absolute_distances:
			distance = abs(distance)

		# 2.3: Append the distance for this coordinate to differences_in_distance
		differences_in_distance.append(distance)

	# Third, return differences_in_distance
	return differences_in_distance

def less_than_or_equal_to_max_bondlength(distance_between_atoms, element1, element2, factor=None):
	"""
	This method will determine if the distance between two atoms can be considered a bond or not.

	Parameters
	----------
	distance_between_atoms : float
		The length between two atoms being considered. Given in Å. 
	element1 : str. or int.
		The element or the index of the first atom. element1 and element2 need to be both strings or ints.
	element2 : str.	or int.
		The element or the index of the second atom. element1 and element2 need to be both strings or ints.
	factor : float
		This is the factor to multiply the maximum covalent bond distance by. If you dont want to give this, set to None. Default: None

	Returns
	-------
	Returns True if the distance between two atoms is short enough to be considered bonding. False if not. 
		
	"""
	if factor is None:
		max_bondlength = get_covalent_bond_distance(element1, element2) * get_bond_length_factor_extension(element1, element2)
	else:
		max_bondlength = get_covalent_bond_distance(element1, element2) * factor
	return (distance_between_atoms <= max_bondlength)

element_bonding_distances = {}
def get_covalent_bond_distance(element1, element2):
	"""
	This method will return the maximum covalent bond distance based on the covalent radii of each element in the bond.

	Parameters
	----------
	element1: str.
		The element of the first atom.
	element2: str.
		The element of the second atom.

	Attributes
	----------
	covalent_bond_distance : float
		The maximum covalent bond length.
	"""

	# First, obtain the key for element_bonding_distances, which is the tuple of each element. 
	key = tuple(sorted([element1,element2]))

	# Second, get the bond length between the elements given if they have already been recorded
	max_bondlength = element_bonding_distances.get(key, None)

	# Third, if key is not in element_bonding_distances, get it and record it to the element_bonding_distances dict. 
	if max_bondlength is None:
		max_bondlength = covalent_radii[atomic_numbers[element1]]+covalent_radii[atomic_numbers[element2]]
		element_bonding_distances[key] = max_bondlength

	# Fourth, return max_bondlength
	return max_bondlength

maximum_bondlength_over_cut_factor = 1.2
def get_bond_length_factor_extension(element1, element2):
	"""
	Increase the max bond length to be considered a bond by the amount given by this method

	Attributes
	----------
	element1 : str.
		The element of atom 1.
	element2 : str.
		The element of atom 2.

	Returns
	-------
	The factor of the maximum bond length to multiply by.
	"""
	if element1 == 'O' and element2 == 'O':
		return 1.02
	else:
		return maximum_bondlength_over_cut_factor

# ---------------------------------------------------------------------------------------------------------------------------

def are_two_values_within_eachother(value1, value2, max_diff):
	"""
	This method is designed to determine if two values are within max_diff of each other. 

	Parameters
	----------
	value1 : float
		This is the first value.
	value2 : float
		This is the second value.
	value1 : float
		This is the maximum difference between value1 and value2.

	Returns
	-------
	True if abs(value1 - value2) <= max_diff. False if not.
	"""
	return (abs(value1 - value2) <= max_diff)

def are_two_lists_within_eachother(list1, list2, max_diff):
	"""
	This method is designed to determine if the values in two list (at the same indices) are within max_diff of each other.

	Parameters
	----------
	list1 : list
		This is the first list.
	list2 : list
		This is the second list.
	max_diff : float
		This is the maximum difference between values in the same index in both list1 and list2.

	Returns
	-------
	True if all abs(list1[index] - list2[index]) <= max_diff. False if not.
	"""
	# First, make sure that the length of each list is the same
	if not (len(list1) == len(list2)):
		raise Exception('issue. list are not the same length.')

	# Second, return if all values at the same index in lists list1 and list2 are within max_diff of each other. 
	return all( (abs(float(value1) - float(value2)) <= max_diff) for (value1, value2) in zip(list1, list2) )

# ---------------------------------------------------------------------------------------------------------------------------





