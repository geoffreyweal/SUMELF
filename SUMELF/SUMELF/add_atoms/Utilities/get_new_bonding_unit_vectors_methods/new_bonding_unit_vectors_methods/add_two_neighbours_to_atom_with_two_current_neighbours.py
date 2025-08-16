"""
add_two_neighbours_to_atom_with_two_current_neighbours.py, Geoffrey Weal, 5/7/2022

This script contains a method for determining the bonding vectors for adding two new hydrogen atoms to the central atom that already contains two neighbours.
"""
import numpy as np
from math import pi

from SUMELF.SUMELF.general_methods.geometry_methods                      import get_unit_vector
from SUMELF.SUMELF.general_methods.geometry_methods                      import planeFit
from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods import get_bond_angle
from SUMELF.SUMELF.general_methods.geometry_methods                      import get_rotation_matrix_around_arbitrary_axis
from SUMELF.SUMELF.general_methods.distance_methods                      import get_distance

def add_two_neighbours_to_atom_with_two_current_neighbours(bonding_unit_vectors, molecule, molecule_graph, index_to_attach_Hs_to, atom_element, no_of_Hs_to_attach, number_of_lone_pairs_of_electrons, new_bonding_unit_vectors_recently_been_created, element_to_attach=None, logger=None):
	"""
	This method will return the bonding vectors for adding two new hydrogens to the atom.

	Note that:

				neighbours_to_add = no_of_Hs_to_attach + number_of_lone_pairs_of_electrons

				total_no_of_neighbours = no_of_Hs_to_attach + no_of_neighbouring_atoms_before_adding_hydrogens + number_of_lone_pairs_of_electrons # (this may exclude the number_of_lone_pairs_of_electrons value)
				

	The types of system that might use this method are:

	no_of_Hs_to_attach | no_of_neighbouring_atoms_before_adding_hydrogens | number_of_lone_pairs_of_electrons
	-------------------|--------------------------------------------------|----------------------------------
			1 		   |                         2                        |                 1                
			2 		   |                         2                        |                 0                

	Parameters
	----------
	bonding_unit_vectors : list of numpy.array 
		These are the bonding unit vectors for the current neighbours that are already bonded to the atom of focus.
	molecule : ase.Atoms
		This is the molecule you want to remove the aliphatic carbons to.
	molecule_graph : networkx.Graph
		This is the graph of this molecule.
	index_to_attach_Hs_to : int.
		This is the index that you want to add hydrogens to.
	atom_element : str. 
		This is the element for the central atom.
	no_of_Hs_to_attach : int.
		This is the number of hydrogens that will be attached to the central atom (number of new bonding unit vectors to obtain).
	number_of_lone_pairs_of_electrons : int.
		These are the number of pairs of lone electrons about the central atom
	element_to_attach : str.
		This is the element you are using this method to attach a new atom to your molecule. If None given, their will be an issue given.
		
	Returns
	-------
	new_bonding_unit_vectors : list of numpy.array 
		These are the bonding unit vectors to add new hydrogens to in the current atom. 
	"""

	# Before beginning, obtain the number of atoms already surrounding the central atom.
	no_of_neighbouring_atoms_before_adding_hydrogens = len(bonding_unit_vectors)
	assert no_of_neighbouring_atoms_before_adding_hydrogens == 2

	# First, get the centre of the bonding_unit_vectors
	central_point = sum(bonding_unit_vectors)/len(bonding_unit_vectors)
	central_point = get_unit_vector(central_point)

	# Second, obtain the plane for the atoms in the double bond and neighbours
	_, normal = planeFit([np.array((0,0,0))] + bonding_unit_vectors)

	# Third, get the perpendicular vector to both central_point and normal with the cross product.
	perpendicular_vector = get_unit_vector(np.cross(central_point, normal))

	# Fourth, obtain the angle that you want between molecules.
	new_bonding_angle = get_bond_angle(atom_element, no_of_Hs_to_attach+no_of_neighbouring_atoms_before_adding_hydrogens, number_of_lone_pairs_of_electrons, element_to_attach=element_to_attach)
	if not isinstance(new_bonding_angle,float):
		import pdb; pdb.set_trace()
		raise Exception('Error: new_bonding_angle need to be a float')

	# Fifth, obtain the matrix for rotating the central_point vector around the perpendicular_vector axis.
	if   no_of_Hs_to_attach == 1: 
		# The new atom will be placed so that it has an angle of new_bonding_angle from each of the existing neighbouring atoms. 
		positive_rotation_matrix_around_perpendicular_vector = get_rotation_matrix_for_two_neighbours_one_atom_one_lone_pair_of_electrons(bonding_unit_vectors, perpendicular_vector, new_bonding_angle)
	elif no_of_Hs_to_attach == 2: 
		# Get the bonding angle for the two bonding vectors to be from each other.
		rotation_angle_for_central_point =  pi - (new_bonding_angle/2.0) # <-- new_bonding_angle given in randians.  180.0-((new_bonding_angle)/2.0)) * (pi/180.0) if new_bonding_angle is in degrees
		positive_rotation_matrix_around_perpendicular_vector = get_rotation_matrix_around_arbitrary_axis(perpendicular_vector, rotation_angle_for_central_point)
		negative_rotation_matrix_around_perpendicular_vector = get_rotation_matrix_around_arbitrary_axis(perpendicular_vector,-rotation_angle_for_central_point)
	else:
		raise Exception('Huh?')

	# Sixth, obtain the positions for bonding the two new hydrogen atoms too, which are obtained by rotating the central_point by 120 and 240 degree around the perpendicular_vector vector. 
	if   no_of_Hs_to_attach == 1: 
		new_hydrogen_1_bond_vector = positive_rotation_matrix_around_perpendicular_vector @ central_point
	elif no_of_Hs_to_attach == 2: 
		new_hydrogen_1_bond_vector = positive_rotation_matrix_around_perpendicular_vector @ central_point
		new_hydrogen_2_bond_vector = negative_rotation_matrix_around_perpendicular_vector @ central_point
	else:
		raise Exception('Huh?')

	# Seventh, return the bonding vectors for the hydrogens, dont worry about returning the bonding vectors for the lone pairs of electrons.
	if   no_of_Hs_to_attach == 1: 
		return [get_unit_vector(new_hydrogen_1_bond_vector)]
	elif no_of_Hs_to_attach == 2: 
		return [get_unit_vector(new_hydrogen_1_bond_vector), get_unit_vector(new_hydrogen_2_bond_vector)]
	else:
		raise Exception('Huh?')

distance_tolerance = 0.001
def get_rotation_matrix_for_two_neighbours_one_atom_one_lone_pair_of_electrons(bonding_unit_vectors, perpendicular_vector_test, theta):
	"""
	This method will return the rotation matrix for an atom that contains two neighbours and has one lone pair of electrons.

	Parameters
	----------
	bonding_unit_vectors : list of numpy.arrays
		These are the bonding unit vectors from the central atom to neighbours that currently exist.
	perpendicular_vector_test : 
		This vector should create a line in the same direction as the perpendicular_vector that will be created in this method. This is for checking.
	theta : float
		This is the angle that we want to rotate the new bonding vector by (given in radians).
	"""

	# Before starting, make sure that this method is accepting only two bonding vectors
	if not (isinstance(bonding_unit_vectors,list) and (len(bonding_unit_vectors) == 2)):
		raise Exception('Can only use this method with two bonding_unit_vectors')

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# First, we will obtain two variables:
	#        * length_of_middle_vector: This should be equal to 1
	#        * perpendicular_vector : This is the perpendicular vector to rotate vectors about.
	#                                 --> This will be checked against perpendicular_vector_test to make sure if it the expected value. 

	# First, get the vector that point up between the vectors in bonding_unit_vectors
	# The length of this vector will be set to 1 (unit vector).
	middle_vector = get_unit_vector(sum(bonding_unit_vectors))
	length_of_middle_vector = np.linalg.norm(middle_vector)

	# Second, name the first vector in bonding_unit_vectors as vector1
	vector1_unit = get_unit_vector(bonding_unit_vectors[0])

	# Third, get the angle between vector1 and middle_vector
	angle_between_central_point = np.arccos(np.dot(vector1_unit, middle_vector))

	# Fourth, based on a right angle triangle, get the length of vector1 so it reaches the same height as central_point
	length_of_vector1 = length_of_middle_vector / np.cos(angle_between_central_point)
	vector1 = length_of_vector1 * vector1_unit

	# Fifth, obtain vector2, which is the top vector for the right angle triangle that contains middle_vector on a side and vector1 on the hypotenuse.
	# The unit vector for vector2 is perpendicular_vector
	length_of_vector2 = np.linalg.norm(middle_vector) * np.tan(angle_between_central_point)
	perpendicular_vector = get_unit_vector((bonding_unit_vectors[0]-bonding_unit_vectors[1])/2)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	vector2 = length_of_vector2 * perpendicular_vector
	if not ((get_distance(perpendicular_vector, perpendicular_vector_test) < distance_tolerance) or (get_distance(-perpendicular_vector, perpendicular_vector_test) < distance_tolerance)):
		print('perpendicular_vector = '+str(perpendicular_vector))
		print('perpendicular_vector_test = '+str(perpendicular_vector_test))
		import pdb; pdb.set_trace()
		raise Exception('perpendular vector should be the same as the test.')

	# Sixth, obtain the length of vector3, which is the hypotenuse of the second triangle that is at right angles to the first triangle along middle_vector
	# Equations about based on trigonomety involving the four trianglar faces of a tetrahedron, where three of the faces are right angle triangles.
	# The equations used ito get this equation are:
	#	1. a^2    = (v1)^2 + (v3)^2 - 2(v1)(v3)cos(theta)
	#	2. a^2    = (v2)^2 + b^2
	#	3. (v3)^2 = m^2 + b^2
	length_of_vector3 = ((length_of_vector1 ** 2.0) + (length_of_middle_vector ** 2.0) - (length_of_vector2 ** 2.0))/(2.0*length_of_vector1*np.cos(theta))

	# Seventh, get the angle between vector3 and middle_vector, which is a right angle triangle. 
	# We will use trigonometry of a right angle triangle with the known lengths of this triangle to get this angle. 
	omega = np.arccos(length_of_middle_vector/length_of_vector3)

	# Eighth, we can now rotate middle_vector about the perpendicular_vector by the angle omega
	rotation_matrix = get_rotation_matrix_around_arbitrary_axis(perpendicular_vector, omega)
	
	# Ninth, we can now rotate middle_vector about the perpendicular_vector by the angle omega
	return rotation_matrix


