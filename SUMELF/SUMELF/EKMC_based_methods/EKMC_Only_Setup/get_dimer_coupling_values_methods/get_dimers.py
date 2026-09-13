"""
get_dimers.py, Geoffrey Weal, 14/11/23

This script is designed to obtain all the dimers that are important in your crystal (for your system values of rCut)
"""

from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_dimer_coupling_values_methods.get_dimers_methods.get_expanded_cell_corner_points import get_expanded_cell_corner_points

zero_point = (0, 0, 0)
def get_dimers(molecules, crystal_cell_lattice, long_range_rCut):
	"""
	This method will obtain all the dimers that are in unit cells, where the unit cells are within long_range_rCut of each other. 

	Parameters
	----------
	molecules : list of ase.Atoms
		This is a list of all the molecules in your crystal
	crystal_cell_lattice : ase.Cell
		This is the unit cell matrix
	long_range_rCut : float.
		Any molecule pairs within this distance of each other will have a coupling given by the large-range coupling model.
		
	Returns
	------- 
	A tuple containing the name of the two molecules in the dimer, and their relative cell ijk difference from molecule 1 --> molecule 2
	"""

	# First, get all unit cells that are within the neighbourhood rCut (long_range_rCut) of the origin unit cell.
	cell_points, surrounding_cell_points = get_expanded_cell_corner_points(crystal_cell_lattice, long_range_rCut)
	all_cell_points = cell_points + surrounding_cell_points

	# Second, convert molecules into a list. 
	molecules_list = tuple(molecules.items())

	# Third, obtain all the kinetic information needed to obtain rate constants between molecule in a crystal for a EKMC simulation.
	for index1, (mol_1_name, molecule1) in enumerate(molecules_list):

		# 3.1: Yield the list of dimers for the same molecule in different cells that are within long_range_rCut
		for cell_point in all_cell_points:
			dimer1 = (mol_1_name, mol_1_name, molecule1.copy(), molecule1.copy(), cell_point)
			yield (dimer1,)

		# 3.2: For other different molecules in the dimer:
		for (mol_2_name, molecule2) in molecules_list[index1+1:]:

			# 3.3: Yield the list of dimers for different molecules in the same (0,0,0) and different cells that are within long_range_rCut 
			for cell_point in [zero_point] + all_cell_points: 
				dimer1 = (mol_1_name, mol_2_name, molecule1.copy(), molecule2.copy(), cell_point)
				minus_cell_point = tuple(-uc for uc in cell_point)
				dimer2 = (mol_2_name, mol_1_name, molecule2.copy(), molecule1.copy(), minus_cell_point)
				yield (dimer1, dimer2)
