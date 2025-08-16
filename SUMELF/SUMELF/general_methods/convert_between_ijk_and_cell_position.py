"""
Cell_Generator_associated_methods.py, Geoffrey Weal, 17/2/22

This script contains methods for dealing with ijk and cell points.
"""
import numpy as np
from numpy.linalg import inv

def convert_ijk_to_displacement_vector(ijk_lengths, crystal_cell_lattice):
	"""
	This method is designed to convert the ijk values of a unit cell into its cartesian positions.

	Parameters
	----------
	ijk_lengths: tuple
		This is the ijk lengths. 
	crystal_cell_lattice : np.array
		A 3x3 row matrix giving the cell vectors in the i, j, and k basis set directions (in cartesian xyz coordinates). Given in Å. 

	Returns
	-------
	position : numpy.array
		This is the position you want to convert into ijk lengths. Given in Å. 
	"""
	# First, obtain the position, given the i, j, k values
	position = crystal_cell_lattice.T @ np.array(ijk_lengths)

	# Second, return position
	return position

def convert_position_into_ijk_lengths(position, crystal_cell_lattice, should_be_interger_length=False):
	"""
	This method will convert the (x,y,z) position given into ijk lengths

	Parameters
	----------
	position : numpy.array
		This is the position you want to convert into ijk lengths. Given in Å. 
	crystal_cell_lattice : np.array
		A 3x3 row matrix giving the cell vectors in the i, j, and k basis set directions (in cartesian xyz coordinates). Given in Å. 
	should_be_interger_length : bool.
		This boolean indicates if you expect the output ijk lengths to be interger values. Default: False. 

	Returns
	-------
	ijk_lengths: tuple
		This is the ijk lengths. 
	"""

	# First, obtain the ijk lengths using matrix linear algebra.
	ijk_lengths_array = inv(crystal_cell_lattice.T) @ position

	# Second, convert the lengths from numpy.array to a list of either floats or ints (depending on your input for should_be_interger_length).
	ijk_lengths = []
	for vector_length in ijk_lengths_array:

		# 2.1: Round the vector length to a sensible number of decimal places.
		vector_length = round(vector_length,6)

		# 2.2: Write your vector_length to the ijk_lengths list
		if should_be_interger_length:
			
			# 2.2.1.1: If you expect your ijk lengths to be intergers, perform this check
			if not vector_length % 1.0 == 0.0:
				raise Exception('not integer unit length: '+str(vector_length))

			# 2.2.1.2: Save the ijk length as an int
			ijk_lengths.append(int(vector_length))

		else:

			# 2.2.2: Save the ijk length as a float
			ijk_lengths.append(vector_length)

	# Third, return the ijk_lengths list as a tuple
	return tuple(ijk_lengths)



