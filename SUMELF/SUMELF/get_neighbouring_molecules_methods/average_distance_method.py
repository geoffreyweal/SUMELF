"""
average_distance_method.py, Geoffrey Weal, 23/2/22

This script obtains neighbours by looking at the average distance between non-hydrogen atoms between two molecules. 
"""
import numpy as np

from SUMELF import get_cell_corner_points, remove_hydrogens, get_distance

def get_neighbours_average_distance_method(molecules, max_neighbour_distance):
	"""
	This method will obtain neighbouring pairs of molecules in the crystal based on the average distance between atoms in each molecule. 

	Parameters
	----------
	molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal that you want to determine neighbouring pairs of molecules for.
	max_neighbour_distance : float.
		This is the maximum distance that any of the atoms between the two molecules can be on average to be considered as a neighbouring pair of molecules. Given in Å. 

	Returns
	-------
	neighbour_pairs : list
		This is a list of all the neighbouring pairs of molecules identified.
	"""

	# ------------------------------------------------------------------------------------------------
	# REMINDER: This method is unfinished and needs looking into before it can be used.
	#
	# What still needs to be decided/fixed below:
	#   * The semantics: the code below accepts a pair if ANY single atom's mean distance to all
	#     atoms in the other molecule is within max_neighbour_distance. Confirm this is the
	#     intended definition of "average distance".
	#   * Undefined names: Neighbour_Generator (never imported), molecule_graphs (not a parameter),
	#     dimer_generator (the generator built above is called neighbour_generator), and
	#     max_dimer_distance (the parameter is called max_neighbour_distance).
	#   * "molecules" is now a dict keyed by molecule name, not a list indexed by integer.
	#   * The return shape must match the other neighbour methods: a flat list of
	#     (mol_name1, mol_name2, unit_cell_displacement, displacement, distance) tuples.
	#     This currently returns a 2-tuple of (dimer_pairs, average_distances) with 3-element entries.
	#   * There is a leftover "import pdb; pdb.set_trace()" further down.
	#   * Decide whether this method should support multiprocessing (no_of_cpus), as
	#     nearest_atoms_method, centre_of_mass_method and centre_of_molecule_method now do.
	# ------------------------------------------------------------------------------------------------
	raise NotImplementedError(
		'The "average_distance_method" neighbour/dimer method is NOT IMPLEMENTED and needs looking into. '
		'See the REMINDER comment at the top of this method for the list of things to resolve. '
		'In the meantime use "nearest_atoms_method", "centre_of_mass" or "centre_of_molecule".'
	)



	# First, create the neighbour generator that will create all the neighbouring pair of molecules that could exist in the crystal.
	neighbour_generator_object = Neighbour_Generator(molecules, molecule_graphs, molecules[0].get_cell())
	neighbour_generator = neighbour_generator_object.generator()

	# Second, remove all hydrogen from the molecule, as we dont want to include hydeogens in our analysis
	molecule_with_no_hydrogens = []
	for molecule in molecules:
		copied_molecule = remove_hydrogens(molecule)
		molecule_with_no_hydrogens.append(copied_molecule)

	# Third, identify dimers based on if there are any atoms in each molecule that are within max_dimer_distance distance of each other
	# I DONT THINK this is WORKING PROPERLY, NEED TO CHECK IT OUT.
	dimer_pairs = []
	shortest_average_distances = []
	for index1, index2, positions1, positions2, displacement, unit_cell_displacement in dimer_generator:

		# if the average distance between atoms between each molecule are on average less than max_dimer_distance, the molecules are considered dimers.
		is_dimer = False
		average_distances_from_mol_2_to_mol_1 = {index: [] for index in range(len(positions2))}
		average_distances = []

		# check average distances of atoms in molecule 1
		for pos_index1 in range(len(positions1)):
			distances = []
			for pos_index2 in range(len(positions2)):
				distance = round(get_distance(positions1[pos_index1],positions2[pos_index2] + displacement),4)
				distances.append(distance)
				average_distances_from_mol_2_to_mol_1[pos_index2].append(distance)
			average_distance = get_average_distance(distances)
			if (average_distance <= max_dimer_distance):
				dimer_pairs.append((index1,index2,displacement))
				average_distances.append(average_distance)
				is_dimer = True
				break

		# if dimer not found, check average distances of atoms in molecule 2
		if not is_dimer:
			for index in range(len(average_distances_from_mol_2_to_mol_1)):
				average_distance = get_average_distance(average_distances_from_mol_2_to_mol_1[index])
				if (average_distance <= max_dimer_distance):
					dimer_pairs.append((index1,index2,displacement))
					average_distances.append(average_distance)
					break

		# 2.3: Send the result of if the dimer was accepted or not back to the generator
		end_of_for_loop_check = dimer_generator.send(is_dimer)
		if not (end_of_for_loop_check == 'Go to for loop'):
			raise Exception('Communication error of Dimer_Generator generator with for loop.')



	# Return the list of dimer pairs.
	return dimer_pairs, average_distances

def get_average_distance(distances):
	"""
	This method will give the average distance from a list of distances. 

	Parameters
	----------
	distances : list
		This is a list of all the distances between an atom in one molecule and all atoms in another molecule. All distances are given in Å. 

	Returns
	-------
	float
		The average distance between an atom in one molecule and all atoms in another molecule. Given in Å. 
	"""
	return sum(distances)/float(len(distances))


