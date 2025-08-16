"""
angle_methods.py, Geoffrey Weal, 19/7/22

This script is designed to hold methods useful for obtaining and dealing with bond angles.
"""
import numpy as np

from SUMELF.SUMELF.general_methods.distance_methods import get_distance
from SUMELF.SUMELF.general_methods.geometry_methods import get_unit_vector

distance_tolerance = 0.1
def get_angle(position_left, position_centre, position_right):
	"""
	This method will give the angle between the position_left - position_centre vector and the position_right - position_centre vector.

	Parameters
	----------
	position_left : np.array
		A 1x3 row matrix containing the x, y, and z coordinates of the atom1. Given in Å. 
	position_centre : np.array
		A 1x3 row matrix containing the x, y, and z coordinates of the atom2. Given in Å. 
	position_right : np.array
		A 1x3 row matrix containing the x, y, and z coordinates of the atom3. Given in Å. 

	Returns
	-------
	angle : float
		This is the angle between these positions, in degrees.
	"""

	# First, make sure that the positions are no too close to each other that it would be unrealistic
	if (get_distance(position_left, position_centre) <= distance_tolerance) or (get_distance(position_centre, position_right) <= distance_tolerance) or (get_distance(position_right, position_left) <= distance_tolerance):
		import pdb; pdb.set_trace()
		raise Exception('Error: The position_left and/or position_right and/or position_centre positions are too close to each other. Uncertainties in positions may lead to inaccurate angle output. Check this out.')

	# Second, get the unit vectors for the position_left - position_centre vector and the position_right - position_centre vector.
	unit_vector_1 = get_unit_vector(position_left  - position_centre)
	unit_vector_2 = get_unit_vector(position_right - position_centre)

	# Third, get the angle between unit_vector_1 and unit_vector_2
	angle = np.arccos(np.dot(unit_vector_1, unit_vector_2))

	# Fourth, return angle in degrees
	return float(np.degrees(angle))

def above_angle_tolerance(vector1, vector2, angle_tolerance, upper_angle_limit, lower_angle_limit=0.0):
	"""
	This method is designed to determine the angle between vector1 and vector2, and determine if this greater than the angle_tolerance.

	Parameters
	----------
	vector1 : np.array
		The first vector.
	vector2 : np.array
		The second vector.
	angle_tolerance : float
		This is the angle that the two vectors must be greater than. Must be given in radians. 
	upper_angle_limit : float
		This is the upper angle limit. Must be given in radians. 
	lower_angle_limit : float
		This is the lower angle limit. Must be given in radians. 
	"""
	return ((lower_angle_limit+angle_tolerance) <= np.arccos(np.dot(vector1, vector2)) <= (upper_angle_limit-angle_tolerance))


