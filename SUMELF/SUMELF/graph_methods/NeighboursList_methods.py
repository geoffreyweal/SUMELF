"""
neighbours_from_xyz_file.py, Geoffrey Weal, 4/8/22

This script is designed read and write the neighbours from the neighbourList of the xyz file.
"""
import numpy as np
from collections import Counter

def get_neighbours_from_NeighboursList(NeighboursList_string_form):
	"""
	This method is designed to convert the neighbours from the neighbourList of the xyz file.

	Parameters
	----------
	NeighboursList_string_form : list or numpy.array
		This is a list of the neighbours for each atom, given as a string in each list entry.

	Returns
	-------
	NeighboursList : list of tuple of ints
		This list contains all the indices of all the neighbours for each atom in the NeighboursList_string_form list.
	"""

	# First, convert NeighboursList to a list if it is currently an array.
	if not isinstance(NeighboursList_string_form, list):
		NeighboursList_string_form = NeighboursList_string_form.tolist()

	# Second, obtain the NeighboursList
	NeighboursList = []
	for neighbours in NeighboursList_string_form:
		if neighbours == '-':
			NeighboursList.append([])
		else:
			neighbours_as_ints = [int(neighbour) for neighbour in neighbours.split(',')]
			NeighboursList.append(neighbours_as_ints)

	# Third, return the NeighboursList
	return NeighboursList

def convert_NeighboursList_to_bonds(NeighboursList):
	"""
	This method is designed to convert the:
		--> NeighboursList: 2D list; where position in list is the index of atom 1 in the bond, and the list contains the indices of atoms bonded to atom 1
		       to
		--> bonds: list o tuple, tuples give bonds between two atoms
	"""

	# First, initialise forward and reverse bond lists.
	bonds_forward = []
	bonds_reverse = []

	# Second, for each atom in NeighboursList
	for atom1_index, bonded_indices in enumerate(NeighboursList):

		# 2.1: For each atom bonded to index 1
		for atom2_index in bonded_indices:

			# 2.2: Add bond to either bonds_forward or bonds_reverse
			if   atom1_index < atom2_index:
				bonds_forward.append((atom1_index, atom2_index))
			elif atom1_index > atom2_index:
				bonds_reverse.append((atom2_index, atom1_index))
			else:
				raise Exception('Error: Atom '+str(atom1_index)+' seems to be bonded to itself in NeighboursList. Check this')

	# Third, check that the same bond has not been entered twice into bonds_forward or bonds_reverse.
	if not (len(bonds_forward) == len(set(bonds_forward))):
		raise Exception('Error: forward bond added more than once in bonds_forward. Check NeighboursList for: '+str([item for item, count in dict(Counter(bonds_forward)).items() if count > 1]))
	if not (len(bonds_reverse) == len(set(bonds_reverse))):
		raise Exception('Error: reverse bond added more than once in bonds_reverse. Check NeighboursList for: '+str([item for item, count in dict(Counter(bonds_reverse)).items() if count > 1]))
	
	# Fourth, sort bonds_forward and bonds_reverse
	bonds_forward.sort()
	bonds_reverse.sort()

	# Fifth, check that bonds_forward and bonds_reverse are the same (they should be)
	if not bonds_forward == bonds_reverse:
		raise Exception('Error, one or more of the bonds in bonds_forward is not found in bonds_reverse. Check NeighboursList for atoms: '+str(set(bonds_forward)))

	# Sixth, return bonds_forward (both bonds_forward and bonds_reverse are the same).
	return bonds_forward


