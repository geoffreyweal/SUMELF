"""
get_expanded_cell_corner_points.py, Geoffrey Weal, 12/8/22

This script is will provide all the points around the origin (0,0,0) points that would allow molecules to be within neighbourhood_rCut distance from a molecule in the unit cell (if it's COM is considered at the origin).
"""
import numpy as np
from itertools import product, count
from SUMELF import get_distance

zero_point = (0,0,0)
zero_point_np = np.array((0,0,0))
def get_expanded_cell_corner_points(crystal_cell_lattice, neighbourhood_rCut):
	"""
	This method will provide all the points around the origin (0,0,0) points that would allow molecules to be within neighbourhood_rCut distance from a molecule in the unit cell (if it's COM is considered at the origin).

	Parameters
	----------
	crystal_cell_lattice : ase.Cell
		This is the unit cell matrix
	neighbourhood_rCut : float
		This is the circumference around the COM of each molecule to add to the neighbourhood.

	Returns
	-------
	cell_points : list
		These are the cell points around the origin within neighbourhood_rCut distance of the origin. Does not include the origin in this list (0,0,0). 
	surrounding_cell_points : list
		These are the cell points around the origin that are slightly outside the neighbourhood_rCut radius of the origin, but may still likely to contain molecules that would be within the neighbourhood_rCut radius.
	"""

	# =================================================================================================================================
	# First, obtain all the lattice point around the (0,0,0) point that are within neighbourhood_rCut distance from this origin point.
	cell_points = [(0,0,0)]
	for max_index in count(start=1, step=1):
		continue_to_expand = False
		# Get all the values around the origin in the current max layer of lattice points
		range_to_scan_over = list(range(-max_index,max_index+1))
		ranges_to_scan_over = [range_to_scan_over for _ in range(3)]
		for cell_point in product(*ranges_to_scan_over):
			# Only consider those lattice point that contain max_index or -max_index in it
			# We want to only consider those points in the max_index layer
			# as we have already analysed those points in the lower layers.
			if not ((max_index in cell_point) or (-max_index in cell_point)):
				continue
			displacement_vector = np.matmul(np.array(cell_point),crystal_cell_lattice)
			distance_from_origin = get_distance(displacement_vector,(0,0,0))
			if distance_from_origin <= neighbourhood_rCut:
				continue_to_expand = True
				cell_points.append(cell_point)
		# If none of the point in the max_index layer were added to the cell_points, then no further/higher
		# layers will be within the circumference of neighbourhood_rCut from the (0,0,0) point, so can stop here
		if not continue_to_expand:
			break

	# =================================================================================================================================
	# Second, we want to also consider all lattice point just adjacent to all points surrounding the cell_points points.
	# This is because after some math the molecules in these cells may just be within the neighbourhood_rCut radius. 
	surrounding_cell_points = []
	for cell_point in cell_points:
		range_to_scan_over = tuple(range(-1,1+1,1))
		cell_point_to_add = [range_to_scan_over for _ in range(3)]
		for cell_point_to_add in product(*cell_point_to_add):
			extra_cell_point = tuple(a+b for a, b in zip(cell_point,cell_point_to_add))
			if not extra_cell_point in cell_points+surrounding_cell_points:
				surrounding_cell_points.append(extra_cell_point)
	if zero_point in surrounding_cell_points:
		surrounding_cell_points.remove(zero_point)

	# =================================================================================================================================
	# Third, remove the (0,0,0) cell from the cell_points. We only want to include non-orgin cells in cell_points (and surrounding_cell_points)
	cell_points.remove((0,0,0))

	# =================================================================================================================================
	# Fourth, return cell_points and surrounding_cell_points
	return sorted(cell_points), sorted(surrounding_cell_points)

	# =================================================================================================================================





