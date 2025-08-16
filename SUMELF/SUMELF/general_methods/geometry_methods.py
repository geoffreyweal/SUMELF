"""
geometry_methods.py, Geoffrey Weal, 19/7/22

This script is designed to hold methods useful for obtaining and manipulating vectors.
"""

import numpy as np
from numpy.linalg import svd
from scipy.spatial.transform import Rotation

# ---------------------------------------------------------------------------------------------------------------------------

def get_unit_vector(vector):
	"""
	This method will turn a vector into a unit vector of length 1.0.

	Parameters
	----------
	vector : np.array
		This is the vector you want to obtained the unit vector for.

	Returns
	-------
	vector_hat : np.array
		This is the unit vector.
	"""

	# First, obtain the unit vector.
	vector_hat = vector / np.linalg.norm(vector)

	# Second, raise an exception if any of the entries in the vector_hat is not a number (NaN).
	if np.isnan(vector_hat).any():
		raise Exception('Error: vector hat has come as NaN. vector_hat = '+str(vector))

	# Third, return the unit vector.
	return vector_hat
	
# ---------------------------------------------------------------------------------------------------------------------------

def rotate_vector_around_axis(vector, angle, axis):
	"""
	This method will rotate a vector around an axis by an angle (in radians). 

	Parameters
	----------
	vector : np.array
		This is the vector you want to rotate around the axis.
	angle : float
		This is the angle you want to rotate you vector around the axis by. This must be given in radians. If you are in degrees, use the numpy.radians method to convert your degrees into radians
	axis : np.array
		This is the axis you want to rotate your vector about. 

	Returns
	-------
	rotated_vector : np.array
		This is the vector that you have rotate around the axis by an angle (in radians). 
	"""

	# First, obtain the rotation vector by multiplying the axis by the angle (in radians). 
	rotation_vector = angle * axis

	# Second, obtain the rotation matrix object.
	rotation = Rotation.from_rotvec(rotation_vector)

	# Third, obtain the vector after it has been rotated.
	rotated_vector = rotation.apply(vector)

	# Fourth, return the vector after it has been rotated.
	return rotated_vector

def project_point_onto_line(point, line_direction, origin):
	"""
	This method will project a point onto a line. 

	Parameters
	----------
	point : np.array
		This is the point that you want to project onto the plane.
	line_direction : np.array
		This is the normal of the plane you want to project the point on.
	origin : np.array
		This is the origin point on the plane.

	Returns
	-------
	projected_point : np.array
		This is the projected point on the line.
	"""

	# First, get the vector from the point to the origin.
	vector = point - origin

	# Second, make sure the line direction is a unit vector
	line_direction_unit_vector = get_unit_vector(line_direction)

	# Third, get the distance of the point from the closest point on the line.
	#        * This is equivalent to projecting vector onto line_direction_unit_vector
	dist = project_u_onto_v(vector, line_direction_unit_vector)
	#dist = np.dot(vector,line_direction_unit_vector)

	# Fourth, get the projection of point onto the line with direction line_direction
	projected_point = point - dist*line_direction_unit_vector

	# Fifth, return the projected point.
	return projected_point

def project_u_onto_v(u, v):
	"""
	This method will obtain the result for projecting the u vector onto the v vector

	Parameters
	----------
	u : numpy.array
		This is the u vector
	v : numpy.array
		This is the vector v to project the vector u upon.

	Returns
	-------
	proj_of_u_on_v : float
		This is the value of the projection of u onto v
	"""

	# First, find norm of the vector v
	v_norm = np.sqrt(sum(v**2.0))    
	  
	# Second, apply the formula as mentioned above for projecting a vector onto another vector find dot product using np.dot()
	proj_of_u_on_v = np.dot(u, v)/(v_norm**2.0)

	# Third, return the value of projecting u upon v
	return proj_of_u_on_v

def project_point_onto_plane(point, normal, origin):
	"""
	This method will project a point onto a plane.

	Parameters
	----------
	point : np.array
		This is the point that you want to project onto the plane.
	normal : np.array
		This is the normal of the plane you want to project the point on.
	origin : np.array
		This is the origin point on the plane.

	Returns
	-------
	projected_point : np.array
		This is the projected point on the plane.
	"""

	# First, get the vector from the point to the origin.
	vector = point - origin

	# Second, make sure the normal is a unit vector
	normal_unit_vector = get_unit_vector(normal)

	# Third, get the distance of the point from the closest point on the plane.
	dist = np.dot(vector,normal_unit_vector)

	# Fourth, get the projection of point onto the plane with normal
	projected_point = point - dist*normal

	# Fifth, return the projected point.
	return projected_point

# ---------------------------------------------------------------------------------------------------------------------------

def get_reflection_matrix_from_plane(normal_vector):
    """
    This method will give the reflection matrix for a given normal vector.

    See https://math.stackexchange.com/questions/2015960/finding-the-matrix-of-a-reflection-in-a-plane for more information, where d = n. 

    Parameters
    ----------
    normal_vector : numpy.array
        This is the normal vector of the plane to reflect over.

    Returns
    -------
    reflection_matrix : 2D numpy.array
        This is the reflection matrix to reflect points across this plane.
    """

    # First, get the unit vector of the normal of the plane to reflect across
    normal_unit_vector = get_unit_vector(normal_vector)[np.newaxis]

    # Second, create the reflection matrix to reflect points across the plane.
    reflection_matrix = np.eye(3) - 2*np.matmul(normal_unit_vector.T,normal_unit_vector)

    # Third, return the reflection matrix
    return reflection_matrix

# ---------------------------------------------------------------------------------------------------------------------------

def planeFit(points):
    """
    This method is designed to obtain the best fitting plane for a set number of points.

	Parameters
	----------
	points : list of numpy.arrays
		These are all the positions involved in the plane

	Returns
	-------
	centre : numpy.array 
		This is the centre of the plane
	normal_vector : numpy.array 
		This is the normal vector for the plane
    """

    # First, get the centre of all the points given
    centre = sum(points)/len(points)

    # Second, re-orientate all point so that the origin is about the centre vector
    points_around_centre = (points - centre).T

    # Third, obtain the normal vector for the plane
    svd_result = svd(points_around_centre)
    left = svd_result[0]
    normal_vector = get_unit_vector(left[:,-1])

    # Fourth, return the centre position and the normal vector for the plane
    return centre, normal_vector

def project_point_onto_plane(vector, centre_position, normal_vector, get_angle_between_vector_and_plane=False):
	"""
	This method will return the vector from projecting vector onto the plane with centre_position and normal_unit_vector.

	Parameters
	----------
	vector : numpy.array
		The point to project onto the plane
	centre_position : numpy.array
		This is the central point on the plane
	normal_vector : numpy.array
		This is the normal for the plane. 

	Returns
	-------
	centre : numpy.array 
		This is the centre of the plane
	normal_vector : numpy.array 
		This is the normal vector for the plane
	"""

	# First, determine the position of the vector with respect to the centre of the plane.
	moved_vector = vector - centre_position

	# Second, make sure that the normal vector is a unit vector.
	normal_unit_vector = get_unit_vector(normal_vector)

	# Third, get the scalar distance from moved_vector to plane along the normal
	dist = np.dot(moved_vector, normal_vector)

	# Fourth, multiply the unit normal vector by the distance, and subtract that vector from your point.
	projected_point = moved_vector - (dist * normal_vector)
	
	# Fifth, obtain the angle between the vector and the plane if desired, and return values.
	if get_angle_between_vector_and_plane:
		angle_between_bond_and_plane = np.arccos(np.dot(get_unit_vector(moved_vector), get_unit_vector(projected_point)))
		return (projected_point + centre_position), angle_between_bond_and_plane
	else:
		return (projected_point + centre_position)

# ---------------------------------------------------------------------------------------------------------------------------

def get_rotation_matrix_around_arbitrary_axis(vector, angle_radians):
	"""
	This method will give the rotation matrix for rotating a point about a vector by an angle.

	Parameters
	----------
	vector : numpy.array 
		This is the vector to rotate a point around
	angle_radians : float
		This is the angle to rotation the point around the vector by. This needs to be given in radians

	Returns
	-------
	rotation_matrix : numpy.array 
		This is the rotation matrix for rotating a point about a vector by an angle.
	"""

	# First, make sure that the vector is a unit vector.
	unit_vector = get_unit_vector(vector)

	# Second, get the rotation matrix.
	# See https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
	rotation_matrix = np.cos(angle_radians) * np.eye(3,3) + np.sin(angle_radians) * get_cross_product_matrix(unit_vector) + (1 - np.cos(angle_radians)) * (np.cross(unit_vector, unit_vector))

	# Third, turn the rotation matrix.
	return rotation_matrix

# ---------------------------------------------------------------------------------------------------------------------------

def get_cross_product_matrix(vector):
	"""
	This method will return the cross product matrix for a given vector.

	Parameters
	----------
	vector : numpy.array 
		This is the vector to obtain the cross product matrix for.

	Returns
	-------
	cross_product_matrix : numpy.array 
		This is the cross product matrix for the given vector.
	"""

	# First, check the input vector is a 3 by 1 vector.
	if not ((vector.shape == (3,)) or (vector.shape == (3,1)) or (vector.shape == (1,3))):
		raise Exception('Huh')

	# Second, obtain the cross product matrix for the given vector
	cross_product_matrix = np.array([[0, -vector[3-1], vector[2-1]], [vector[3-1], 0, -vector[1-1]], [-vector[2-1], vector[1-1], 0]])

	# Third, return the cross product matrix.
	return cross_product_matrix

# ---------------------------------------------------------------------------------------------------------------------------






