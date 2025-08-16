"""
perform_symmetry_operation_upon_molecule.py, Geoffrey Weal, 29.5.22

This script is designed will move the molecule to the new position after the symmetry_operation symmetry operation.
"""
import numpy as np

def perform_symmetry_operation_upon_molecule(molecule, symmetry_operation):
	"""
	This method will move the molecule to the new position after the symmetry_operation symmetry operation.

	Parameters
	----------
	molecule : ase.Atoms
		This is the molecule you want to move with the symmetry operation. The molecule object should have the unit cell given in it.
	symmetry_operation : numpy.array(3*3), numpy.array(3*1)
		These are the rotation/reflection and translation matrices for the symmetry operation you want to perform on the molecule. 
	
	Returns
	-------
	molecule : ase.Atoms
		This is the molecule that has now had the symmetry_operation symmetry operation performed upon it.
	"""

	# Initial step, check that the molecule contains a cell.
	if molecule.get_cell() is None:
		raise Exception('Error. To use the perform_symmetry_operation_upon_molecule method, the molecule must have the unit cell given. See https://wiki.fysik.dtu.dk/ase/ase/atoms.html')

	# First, extract the rotation/reflection matrix and the translation vector from symmetry_operation 
	scaled_rotation_matrix, scaled_translation_vector = symmetry_operation

	# Second, obtain the cartesian translation that has been repeated so it can be added to each column in obtaining fractional_positions_after_operation
	scaled_translation_column_vector = scaled_translation_vector.reshape((3,1))
	scaled_translation_vectors = np.repeat(scaled_translation_column_vector, len(molecule), axis=1)

	# Third, obtain the positions of the molecule after the symmetry operation. 
	fractional_positions_after_operation = ((scaled_rotation_matrix @ molecule.get_scaled_positions(wrap=False).T) + scaled_translation_vectors).T

	# Fourth, get the positions of this new molecule so that the atom have been moved by the current symmetry operation.
	molecule.set_scaled_positions(fractional_positions_after_operation)

	# Fifth, return the molecule that has now been moved by the symmetry operation.
	return molecule
