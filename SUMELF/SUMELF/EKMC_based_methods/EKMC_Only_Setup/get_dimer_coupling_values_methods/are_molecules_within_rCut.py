"""
are_molecules_within_rCut.py, Geoffrey Weal, 15/11/23

This script is designed to determine if two molecules that are a relative cell distance away (given by cell_point) are within range of eachother.
"""

import numpy as np
from SUMELF import remove_hydrogens
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_dimer_coupling_values_methods.are_molecules_within_rCut_methods.get_distance_between_molecules import diff_in_centre_of_mass, diff_in_centre_of_molecule, nearest_distance_between_molecules

def are_molecules_within_rCut(molecule1, molecule2, cell_point, crystal_cell_lattice, short_range_rCut, long_range_rCut, rCut_mol_dist_description):
	"""
	This method will determine if two molecules that are a relative cell distance away (given by cell_point) are within range of eachother.

	Parameters
	----------
	molecule1 : ase.Atoms
		This is the first molecule to compare.  This molecule will be found in the (0,0,0) relative unit cell.
	molecule2 : ase.Atoms
		This is the second molecule to compare. This molecule will be found in the (cell_point) relative unit cell.
	cell_point : 3 tuple
		This is the relative cell distance that molecule1 and molecule2 are from each other. 
	crystal_cell_lattice : ase.Cell
		This is the unit cell matrix
	short_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the short-range coupling model.
	long_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the large-range coupling model.
	rCut_mol_dist_description : str.
		This is the description given to how the distance between two molecules is obtained. 
		
	Returns
	------- 
		'Short Range' : dimer_distance <= short_range_rCut
		'Long Range'  : short_range_rCut < dimer_distance <= long_range_rCut
		None          : long_range_rCut < dimer_distance
	"""

	# First, obtain the displacement vector for the difference in cell from molecule 1 to molecule 2.
	displacement_vector = np.matmul(np.array(cell_point), crystal_cell_lattice)

	# Second, check that rCut_mol_dist_description is one of the following, and get the distance between the molecules in the dimer.
	if rCut_mol_dist_description == 'centre_of_mass':
		dimer_distance = diff_in_centre_of_mass(molecule1, molecule2, displacement_vector)
	elif rCut_mol_dist_description == 'centre_of_molecule':
		dimer_distance = diff_in_centre_of_molecule(molecule1, molecule2, displacement_vector)
	elif rCut_mol_dist_description == 'nearest_atoms_method':
		dimer_distance = nearest_distance_between_molecules(remove_hydrogens(molecule1.copy()).get_positions(), remove_hydrogens(molecule2.copy()).get_positions(), displacement_vector)
	else:
		toString  = 'Error: Your input for "rCut_mol_dist_description" must be either:\n'
		toString += '\t* "centre_of_mass":       Take the dimer distance from the centre of masses of each molecule.\n'
		toString += '\t* "centre_of_molecule":   Take the dimer distance from the centre of molecules of each molecule (excluding hydrogen atoms).\n'
		toString += '\t* "nearest_atoms_method": Take the dimer distance as the closest distance between atoms in each molecule (excluding hydrogen atoms).\n'
		toString += 'Check this in your input python script and rerun.'
		raise Exception(toString)

	# Third, obtain the coupling for this dimer set based on dimer_distance vs short_range_rCut/long_range_rCut
	if dimer_distance <= short_range_rCut: 
		return 'Short Range'
	elif short_range_rCut < dimer_distance <= long_range_rCut: 
		return 'Long Range'
	
	# Fourth, if you got here, the dimer pair is outside long_range_rCut, so return None.
	return None


