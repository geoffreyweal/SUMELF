"""
centre_of_mass_method.py, Geoffrey Weal, 17/2/22

This script includes methods for obtaining neighbouring pairs of molecules using the centre of mass.
"""
import numpy as np

from SUMELF import get_cell_corner_points, get_centre_of_mass
from SUMELF.SUMELF.get_neighbouring_molecules_methods.supporting_methods.centre_based_neighbours import get_neighbours_from_centres

def get_neighbours_centre_of_mass_method(molecules, max_neighbour_distance, include_hydrogens_in_neighbour_analysis=False, no_of_cpus=1):
	"""
	This method will obtain neighbouring pairs of molecules in the crystal based on the distances between their centre of masses.

	Parameters
	----------
	molecules : dict. of ase.Atoms objects
		These are all the individual molecules identified in the crystal that you want to determine neighbouring pairs of molecules for, keyed by molecule name.
	max_neighbour_distance : float.
		This is the maximum distance that the centre of mass between two molecules can be to be considered a neighbouring pair of molecules. Given in Å.
	include_hydrogens_in_neighbour_analysis : bool.
		This tag indicates if you want to include hydrogens when determining the centre of mass of each molecule. Default: False
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large systems you may want to implement multiple cpus.

	Returns
	-------
	neighbourhood_molecules_info : list
		This is a list of all the neighbouring pairs of molecules identified, each given as
		(mol_name1, mol_name2, unit_cell_displacement, displacement, distance).
	"""

	# First, obtain the names of the molecules in a consistent order.
	mol_names = sorted(molecules.keys())

	# Second, get the cell points of the super cell around the origin unit cell with reach 1.
	crystal_cell_lattice = molecules[mol_names[0]].get_cell()
	cell_points, unit_cell_displacements = get_cell_corner_points(crystal_cell_lattice, super_cell_reach=1, get_corrspeonding_ijk_values=True)

	# Third, reduce each molecule down to its centre of mass.
	centres = {mol_name: get_centre_of_mass_of_molecule(molecules[mol_name], include_hydrogens_in_neighbour_analysis) for mol_name in mol_names}

	# Fourth, obtain and return the neighbouring pairs of molecules based on those centres of mass.
	return get_neighbours_from_centres(mol_names, centres, cell_points, unit_cell_displacements, max_neighbour_distance, no_of_cpus=no_of_cpus)

# ===============================================================================================================

def get_centre_of_mass_of_molecule(molecule, include_hydrogens_in_neighbour_analysis):
	"""
	This method will obtain the centre of mass of a molecule, optionally ignoring its hydrogens.

	Parameters
	----------
	molecule : ase.Atoms object
		This is the molecule to obtain the centre of mass of.
	include_hydrogens_in_neighbour_analysis : bool.
		This tag indicates if you want to include hydrogens when determining the centre of mass.

	Returns
	-------
	The centre of mass of the molecule, as a numpy.array.
	"""

	# First, obtain the elements and positions of the atoms in this molecule.
	elements  = molecule.get_chemical_symbols()
	positions = molecule.get_positions()

	# Second, remove the hydrogens from consideration if the user does not want to include them.
	if not include_hydrogens_in_neighbour_analysis:
		non_hydrogens = [index for index, element in enumerate(elements) if (element != 'H')]
		# Only exclude the hydrogens if doing so would leave some atoms to work with.
		if len(non_hydrogens) > 0:
			elements  = [elements[index]  for index in non_hydrogens]
			positions = [positions[index] for index in non_hydrogens]

	# Third, return the centre of mass of the atoms being considered.
	return get_centre_of_mass(elements, np.array(positions))

# ===============================================================================================================
