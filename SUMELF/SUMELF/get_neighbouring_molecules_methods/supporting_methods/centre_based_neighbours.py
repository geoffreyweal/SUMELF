"""
centre_based_neighbours.py, Geoffrey Weal, 13/9/25

This script holds the machinery shared by the neighbour methods that compare molecules by a single
representative point each (their centre of mass, or the centre of the molecule).

Both of those methods do the same thing once each molecule has been reduced to one point: compare
every pair of points across the surrounding unit cells and record the pairs that lie within
max_neighbour_distance of each other. That shared work lives here so the two methods only need to
say how they turn a molecule into a point.
"""
import sys
import numpy as np
import multiprocessing as mp
from tqdm import tqdm

from SUMELF import get_distance

def get_neighbours_from_centres(mol_names, centres, cell_points, unit_cell_displacements, max_neighbour_distance, no_of_cpus=1):
	"""
	This method will obtain the neighbouring pairs of molecules in the crystal, based on the distance
	between the representative centre point of each molecule.

	Parameters
	----------
	mol_names : list of str.
		These are the names of the molecules in the crystal, in a consistent order.
	centres : dict.
		This is the representative point of each molecule in the crystal, keyed by molecule name.
	cell_points : list of numpy.array
		These are the displacements of the surrounding unit cells, given in Å.
	unit_cell_displacements : list
		These are the ijk values of the surrounding unit cells.
	max_neighbour_distance : float.
		This is the maximum distance between the centres of two molecules for them to be considered neighbours. Given in Å.
	no_of_cpus : int.
		This is the number of cpus available to use. In most cases this should just be set to 1 cpu, however for very large systems you may want to implement multiple cpus.

	Returns
	-------
	neighbourhood_molecules_info : list
		This is a list of all the neighbouring pairs of molecules identified, each given as
		(mol_name1, mol_name2, unit_cell_displacement, displacement, distance).
	"""

	# First, use either the single or multi-cpu method for obtaining neighbouring molecules.
	if no_of_cpus == 1:
		return obtain_neighbours_with_single_cpu(mol_names, centres, cell_points, unit_cell_displacements, max_neighbour_distance)
	return obtain_neighbours_with_multi_cpu(mol_names, centres, cell_points, unit_cell_displacements, max_neighbour_distance, no_of_cpus)

# ===============================================================================================================

def get_inputs(mol_names, centres, cell_points, unit_cell_displacements, max_neighbour_distance, neighbourhood_molecules_info=None):
	"""
	This generator is designed to obtain the inputs for comparing each molecule against the others in the crystal.

	Each yielded entry covers one molecule and every molecule at or after it in mol_names. index2 starts at
	index1 rather than index1+1 so that a molecule can pair with itself in an adjacent unit cell. Pairs where
	index2 < index1 are not included, because a (index1,index2) pair is the same as a (index2,index1) pair,
	even if the two molecules are chemically different.
	"""

	# First, yield the comparison work for each molecule in turn.
	for index1 in range(len(mol_names)):
		mol_names2 = [mol_names[index2] for index2 in range(index1, len(mol_names))]
		yield (mol_names[index1], mol_names2, centres, cell_points, unit_cell_displacements, max_neighbour_distance, neighbourhood_molecules_info)

def compare_molecule_against_others(mol_name1, mol_names2, centres, cell_points, unit_cell_displacements, max_neighbour_distance):
	"""
	This method will obtain all the neighbouring pairs between one molecule and the other molecules given.

	Returns
	-------
	A list of the neighbouring pairs found for this molecule.
	"""

	# First, initialise the list to record neighbouring pairs in, and obtain the centre of molecule 1.
	neighbouring_pairs = []
	origin_cell_point  = np.array((0,0,0))
	centre1            = centres[mol_name1]

	# Second, compare molecule 1 against each of the other molecules given.
	for mol_name2 in mol_names2:
		centre2 = centres[mol_name2]
		for displacement, unit_cell_displacement in zip(cell_points, unit_cell_displacements):

			# Third, do not include a neighbouring pair of a molecule with itself in the same position.
			if (mol_name1 == mol_name2) and (displacement == origin_cell_point).all():
				continue

			# Fourth, if the distance between the centres is within max_neighbour_distance, these two molecules are a neighbouring pair.
			distance = round(get_distance(centre1, centre2 + displacement), 4)
			if distance <= max_neighbour_distance:
				neighbouring_pairs.append((mol_name1, mol_name2, unit_cell_displacement, displacement, distance))

	# Fifth, return the neighbouring pairs found for this molecule.
	return neighbouring_pairs

# ===============================================================================================================

def obtain_neighbours_with_single_cpu(mol_names, centres, cell_points, unit_cell_displacements, max_neighbour_distance):
	"""
	This method will obtain the neighbouring pairs of molecules in the crystal using a single cpu.
	"""

	# First, initialise the list to record the neighbouring pairs found.
	neighbourhood_molecules_info = []

	# Second, compare each molecule against the others.
	inputs = get_inputs(mol_names, centres, cell_points, unit_cell_displacements, max_neighbour_distance)
	for mol_name1, mol_names2, _, _, _, _, _ in tqdm(inputs, total=len(mol_names), desc='Obtaining neighbouring pairs of molecules', unit='molecule'):
		neighbourhood_molecules_info += compare_molecule_against_others(mol_name1, mol_names2, centres, cell_points, unit_cell_displacements, max_neighbour_distance)

	# Third, return the list of neighbouring pairs of molecules.
	return neighbourhood_molecules_info

def obtain_neighbours_method_for_multi_cpu(input_variables):
	"""
	This method will obtain the neighbouring pairs for one molecule against the others, for running across multiple cpus.
	"""

	# First, obtain the variables for processing from the input_variables tuple.
	mol_name1, mol_names2, centres, cell_points, unit_cell_displacements, max_neighbour_distance, neighbourhood_molecules_info = input_variables

	# Second, record the neighbouring pairs found for this molecule in the shared list.
	neighbourhood_molecules_info += compare_molecule_against_others(mol_name1, mol_names2, centres, cell_points, unit_cell_displacements, max_neighbour_distance)

def obtain_neighbours_with_multi_cpu(mol_names, centres, cell_points, unit_cell_displacements, max_neighbour_distance, no_of_cpus=1):
	"""
	This method will obtain the neighbouring pairs of molecules in the crystal across multiple cpus.
	"""

	# First, create the manager to save lists to.
	with mp.Manager() as manager:

		# Second, create the list to collect information on the neighbouring molecules in the crystal.
		neighbourhood_molecules_info = manager.list()

		# Third, write a warning message to the user.
		print('Obtaining neighbourhoods between molecules (Please wait until after 100%, as the process will still be running.)', file=sys.stderr)

		# Fourth, obtain the input value generator.
		input_values = tqdm(get_inputs(mol_names, centres, cell_points, unit_cell_displacements, max_neighbour_distance, neighbourhood_molecules_info), total=len(mol_names), desc='Obtaining neighbouring pairs of molecules', unit='molecule')

		# Fifth, run the multiprocessing jobs. The result of map_async is waited on with .get() so that
		# any exception raised inside a worker is re-raised here rather than being silently discarded.
		pool = mp.Pool(processes=no_of_cpus)
		result = pool.map_async(obtain_neighbours_method_for_multi_cpu, input_values)
		pool.close()
		result.get()
		pool.join()

		# Sixth, save the list from a multiprocessing Manager list to a regular list.
		neighbourhood_molecules_info = list(neighbourhood_molecules_info)

	# Seventh, return the list of neighbouring pairs of molecules.
	return neighbourhood_molecules_info

# ===============================================================================================================
