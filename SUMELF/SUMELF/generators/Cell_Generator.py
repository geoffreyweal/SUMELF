"""
Cell_Generator.py, Geoffrey Weal, 17/2/22

This script is designed to generate the cells to look at.
"""
from itertools import product

class Cell_Generator:

	"""
	This class is designed to generate the cells to look at.
	"""
	def __init__(self):
		self.directions = ['pi', 'ni', 'pj', 'nj', 'pk', 'nk'] # p means the positive direction, n means the negative direction.
		self.ijk_directions = {direction: 1    for direction in self.directions}
		self.ijk_expand     = {direction: True for direction in self.directions}
		self.all_are_molecules_within_max_dimer_distance = []

	def add_result(self, are_molecules_within_max_dimer_distance):
		"""
		This method is designed to add the are_molecules_within_max_dimer_distance boolean to the self.all_are_molecules_within_max_dimer_distance list

		Parameters
		----------
		are_molecules_within_max_dimer_distance : bool.
			Did the dimer created using the yielded cell from the generate_next_cell method generate a dimer that was stored by the ECCP program. 
		"""
		self.all_are_molecules_within_max_dimer_distance.append(are_molecules_within_max_dimer_distance)

	def generate_next_ijk_points(self):
		"""
		This method will generate the next cell to look at.

		This method works by gradually expanding the ijk points given from the origin in the postive and negative i, j, and k values. This is like expanding a box that has 6 sides.

		If no new dimers are found by expanding in one direction, then we dont worry about further expanding in that direction, as no new dimers will be found by expanding the cell.

		This is like expanding the box in one direction until no new dimers are made, then that side of the box is not expanded any further. 

		Returns
		-------
		cell_point : tuple of three ints
			This is the cell ijk point for the Dimer_Generator to use next. 
		"""

		# First, yield the origin cell point
		yield (0,0,0)

		# Second, gradually expand the reach of this method to examine ijk points that could provide a dimer
		while any(self.ijk_expand.values()):

			# Third, get the current postive and negative i, j, and k values. 
			pi, ni, pj, nj, pk, nk = [self.ijk_directions[direction] for direction in self.directions]

			# Fourth, get the range of values between the positive and negative values in the i, j, and k directions. 
			i_range = tuple(range(-ni,pi+1,1))
			j_range = tuple(range(-nj,pj+1,1))
			k_range = tuple(range(-nk,pk+1,1))

			# Fifth, obtain all the cell points for the positive and negative values in the i, j, and k walls.
			pi_wall_cell_points = sorted(product((pi,), j_range, k_range))
			ni_wall_cell_points = sorted(product((-ni,), j_range, k_range))
			pj_wall_cell_points = sorted(product(i_range, (pj,), k_range))
			nj_wall_cell_points = sorted(product(i_range, (-nj,), k_range))
			pk_wall_cell_points = sorted(product(i_range, j_range, (pk, )))
			nk_wall_cell_points = sorted(product(i_range, j_range, (-nk, )))
			all_cell_points = [pi_wall_cell_points, ni_wall_cell_points, pj_wall_cell_points, nj_wall_cell_points, pk_wall_cell_points, nk_wall_cell_points]
			all_cell_points = {ijk_value: cell_points for ijk_value, cell_points in zip(self.directions, all_cell_points)}

			# Sixth, go through all directions. 
			for direction in self.directions:
				
				# 6.1: If you are still wanting to expand in direction, ...
				if self.ijk_expand[direction]:

					# 6.1.1: First, clear the all_are_molecules_within_max_dimer_distance list for saving if a dimer was made from this cell or not by the get_dimer method you use
					self.all_are_molecules_within_max_dimer_distance = []
					
					# 6.1.2: Yield all the cell points that cover the wall in direction (for example, if pi, all the cell points on the wall of the positive i direction).
					# During the for loop of the Dimer_Generator, the result of if a dimer was created should be appended to the all_are_molecules_within_max_dimer_distance 
					# list using the add_result method. 
					for cell_point in all_cell_points[direction]:
						yield cell_point

					# 6.1.3: Determine if this direction need to be expanded anymore or if all the dimers that could be created in this direction has been created.
					if not any(self.all_are_molecules_within_max_dimer_distance):
						self.ijk_expand[direction] = False
					else:
						self.ijk_directions[direction] += 1
