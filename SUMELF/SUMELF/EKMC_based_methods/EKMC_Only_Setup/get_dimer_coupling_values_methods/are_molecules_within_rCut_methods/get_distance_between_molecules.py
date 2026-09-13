"""
get_distance_between_molecules.py, Geoffrey Weal, 14/12/23

This script contains methods for measureing the distance between two molecules in a crystal. 
"""

from SUMELF import get_distance

def diff_in_centre_of_mass(molecule1, molecule2, displacement_vector):
	"""
	This method will return the distance between two molecules centre of masses, relative to the displacement_vector of the unit cells they are in.

	Parameters
	----------
	molecule1 : ase.Atoms
		This is the first molecule.
	molecule1 : ase.Atoms
		This is the second molecule.
	displacement_vector : numpy.array
		This is the displacement vector to move molecule 2 by.

	Returns
	-------
	distance : float
		The short distance between two molecules. 
	"""

	# First, get the centre of masses of each molecule. 
	centre_of_mass_molecule1 = molecule1.get_center_of_mass()
	centre_of_mass_molecule2 = molecule2.get_center_of_mass()

	# Second, get the distance between the centre of masses of molecules 1 and 2.
	distance = round(get_distance(centre_of_mass_molecule1, centre_of_mass_molecule2 + displacement_vector),4)

	# Third, return the distance between the centre of masses of molecules 1 and 2.
	return distance

def diff_in_centre_of_molecule(molecule1, molecule2, displacement_vector):
	"""
	This method will return the distance between two molecules centre of molecules (excluding hydrogen atoms), relative to the displacement_vector of the unit cells they are in.

	Parameters
	----------
	molecule1 : ase.Atoms
		This is the first molecule.
	molecule1 : ase.Atoms
		This is the second molecule.
	displacement_vector : numpy.array
		This is the displacement vector to move molecule 2 by.

	Returns
	-------
	distance : float
		The short distance between two molecules. 
	"""

	# First, get the centre of molecule of each molecule, excluding hydrogens
	centre_of_molecule_molecule1 = centre_of_molecule(molecule1, include_hydrogen=False)
	centre_of_molecule_molecule2 = centre_of_molecule(molecule2, include_hydrogen=False)

	# Second, get the distance between the centre of masses of molecules 1 and 2.
	distance = round(get_distance(centre_of_molecule_molecule1, centre_of_molecule_molecule2 + displacement_vector),4)

	# Third, return the distance between the centre of masses of molecules 1 and 2.
	return distance

def nearest_distance_between_molecules(positions1, positions2, displacement_vector):
	"""
	This method will return the nearest distance between molecules.

	Parameters
	----------
	positions1 : numpy.array
		These are the positions of atoms in molecule 1.
	positions2 : numpy.array
		These are the positions of atoms in molecule 2.
	displacement_vector : numpy.array
		This is the displacement vector to move molecule 2 by.

	Returns
	-------
	shortest_distance : float
		The short distance between two molecules. 
	"""

	# First, determine the shortest distance between the two molecules.
	shortest_distance = float('inf')
	for pos_index1 in range(len(positions1)):
		for pos_index2 in range(len(positions2)):
			distance = round(get_distance(positions1[pos_index1], positions2[pos_index2] + displacement_vector),4)
			# If we have found a distance that is less than max_dimer_distance, then these two molcules are a dimer.
			if (distance <= shortest_distance):
				shortest_distance = distance

	# Second, determine the shortest distance between the two molecules of interest. 
	return shortest_distance

# ============================================================================================================================================

