"""
supporting_methods.py, Geoffrey Weal, 5/7/2022

This script contains methods used by the python scripts in this folder.
"""
from math import pi
import numpy as np

from SUMELF.SUMELF.general_methods.distance_methods import get_distance

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

	raise Exception('HGer this from geometry methods')

	# First, find norm of the vector v
	v_norm = np.sqrt(sum(v**2.0))    
	  
	# Second, apply the formula as mentioned above for projecting a vector onto another vector find dot product using np.dot()
	proj_of_u_on_v = np.dot(u, v)/(v_norm**2.0)

	# Third, return the value of projecting u upon v
	return proj_of_u_on_v

same_atom_distance_tolerance = 0.05
def record_hydrogen_positions(acceptor_info):
	"""
	This method is designed to record all the positiions of hydrogens bound to hydrogen acceptors in your H-bonding system.

	Parameters
	----------
	acceptor_info : list
		This is all the information about the atoms involved in the H-bonding system

	Returns
	-------
	hydrogens_bound_to_acceptors_positions : list of numpy.arrays
		This is a list of all the positions of all the hydrogens involved in your H-bonding system.
	"""

	# First, obtain all the positions of all the hydrogens involved in your H-bonding system.
	hydrogens_bound_to_acceptors_positions = []
	for _, _, hydrogen_positions in acceptor_info:
		for hydrogen_position in hydrogen_positions:
			# Check that this hydrogen atom has not already been recorded.
			for recorded_hydrogen_position in hydrogens_bound_to_acceptors_positions:
				if get_distance(hydrogen_position,recorded_hydrogen_position) <= same_atom_distance_tolerance:
					break
			else:
				hydrogens_bound_to_acceptors_positions.append(hydrogen_position)

	# Second, return hydrogens_bound_to_acceptors_positions
	return hydrogens_bound_to_acceptors_positions

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------







