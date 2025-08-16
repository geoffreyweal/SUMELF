"""
supporting_methods.py, Geoffrey Weal, 5/7/2022

This script contains methods used by the python scripts in this folder.
"""
import numpy as np
from random import uniform
from math import pi as math_PI

from SUMELF.SUMELF.general_methods.angle_methods                                                                                                               import above_angle_tolerance
from SUMELF.SUMELF.general_methods.geometry_methods                                                                                                            import get_unit_vector, project_u_onto_v
from SUMELF.SUMELF.general_methods.distance_methods                                                                                                            import get_distance
from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods                                                                                       import get_bond_lengths
from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods                                                                                       import get_bond_angle
from SUMELF.SUMELF.general_methods.geometry_methods                                                                                                            import get_rotation_matrix_around_arbitrary_axis
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.get_new_bonding_vector_methods.supporting_methods import record_hydrogen_positions

same_atom_distance_tolerance = 0.05
angle_tolerance = 15.0 * (math_PI/180.0)
def get_best_bonding_vector_with_one_neighbour(molecule, donor_index, neighbour_index, total_no_of_neighbours, no_of_pairs_of_lone_electrons, wrapped_crystal, crystal_graph, acceptor_info, element_to_attach=None):
	"""
	This method will provide a optimal bonding unit vector for hydrogen bonding system, where the H-bond donor is bound to one neighbouring atom.

	Parameters
	----------
	molecule : ase.Atoms
		This is the molecule of interest in the crystal.
	donor_index : int
		This is the index of the H-bond donor
	neighbour_index : int
		This is the index of the atom that neighbours the donor_index atom.
	total_no_of_neighbours : int
		This is the total number of neighbours that we want to surround the central atom
	no_of_pairs_of_lone_electrons : int
		This is the number of lone pairs of electrons surrounding the central atom.
	wrapped_crystal : ase.Atoms
		This is the crystal object.
	crystal_graph : networkx.Graph
		This is the networkx graph for the full crystal.
	acceptor_info : list of ints
		These are the indices of all the H-bond acceptors in the crystal, and their positions.
	element_to_attach : str.
		This is the element you are using this method to attach a new atom to your molecule. If None given, their will be an issue given.

	Returns
	-------
	new_bonding_unit_vector : numpy.array
		This is the newly created bonding unit vector
	"""

	# First, get the bond length betwen atoms and element_to_attach
	bond_lengths = get_bond_lengths(element_to_attach=element_to_attach)

	# Second, obtain the element and position of the H-bond donor
	donor_symbol   = molecule[donor_index].symbol
	donor_position = molecule[donor_index].position

	# Third, obtain the charges and the positions of the nearby hydrogen bonding acceptors
	acceptor_charges = [round(wrapped_crystal[acceptor_index].charge,5) for acceptor_index, _, _ in acceptor_info]
	acceptor_charges = [(acceptor_charge - 0.5) for acceptor_charge in acceptor_charges]
	acceptor_positions = [acceptor_position for _, acceptor_position, _ in acceptor_info]

	# Fourth, obtain the charges and the positions of any hydrogens that are bound to hydrogen bonding acceptors.
	hydrogens_bound_to_acceptors_positions = record_hydrogen_positions(acceptor_info)
	hydrogens_bound_to_acceptors_charges = [+0.5 for hydrogen_position in hydrogens_bound_to_acceptors_positions]		

	# Fifth, obtain the position of the neighbour atom currently bonded to the H-bond donor
	if isinstance(neighbour_index,int):
		neighbour_position = molecule[neighbour_index].position
	elif isinstance(neighbour_index,np.ndarray):
		new_bond_vector_length = 1.0 # This can arbitrary be 1.0, as we only need neighbour_position for obtaining the normal_axis, which is a unit vector.
		neighbour_position = donor_position + (new_bond_vector_length * neighbour_index)
	else:
		raise Exception('Huh')

	# Sixth, get the unit vector along the neighbour_position donor_position axis, point away from neighbour_position
	normal_axis = -get_unit_vector(neighbour_position - donor_position)

	# Seventh, get the distances between donor and acceptor atoms
	force_of_acceptor_upon_donor = []
	for acceptor_position, acceptor_charge in zip(acceptor_positions+hydrogens_bound_to_acceptors_positions, acceptor_charges+hydrogens_bound_to_acceptors_charges):
		if get_distance(acceptor_position, donor_position) <= same_atom_distance_tolerance:
			continue
		force_of_acceptor_upon_donor.append(((0.5*acceptor_charge)/(get_distance(acceptor_position,donor_position)**2.0), acceptor_position))

	# Eighth, sort the acceptor_positions by how close they are to donor_position and their charge
	force_of_acceptor_upon_donor.sort(reverse=True, key=lambda x:x[0])

	# ------------------------------------------------------------------------------------------------
	# Ninth, determine the vector for obtaining a vector that is perpendicular to the normal_axis

	# 9.1: Obtain an acceptor_unit_vector that will produce a good perpendicular vector
	for _, acceptor_position in force_of_acceptor_upon_donor:
		acceptor_unit_vector = get_unit_vector(acceptor_position - donor_position)
		if above_angle_tolerance(acceptor_unit_vector, normal_axis, angle_tolerance, upper_angle_limit=math_PI):
			break
	else:
		for acceptor_unit_vector in [get_unit_vector(np.array((0,0,1))), get_unit_vector(np.array((1,0,0))), get_unit_vector(np.array((0,1,0)))]:
			if above_angle_tolerance(acceptor_unit_vector, normal_axis, angle_tolerance, upper_angle_limit=math_PI):
				break
		else:
			for _ in range(100):
				acceptor_unit_vector = get_unit_vector(np.array((uniform(-1.0,1.0), uniform(-1.0,1.0), uniform(-1.0,1.0))))
				import pdb; pdb.set_trace()
				if above_angle_tolerance(acceptor_unit_vector, normal_axis, angle_tolerance, upper_angle_limit=math_PI):
					break
			else:
				raise Exception('Huh? There has been repetition in a neverending for loop.')

	# 9.2: Project the acceptor_unit_vector onto the normal_axis
	projection_of_acceptor_unit_vector_onto_normal_axis = project_u_onto_v(acceptor_unit_vector, normal_axis)

	# 9.3: Obtain the perpendicular_vector that is perpendicular to the normal_axis, pointing towards the acceptor_position
	perpendicular_vector = acceptor_unit_vector - projection_of_acceptor_unit_vector_onto_normal_axis * normal_axis

	# 9.4: Turn the perpendicular_vector into a unit vector.
	perpendicular_unit_vector = get_unit_vector(perpendicular_vector)
	# ------------------------------------------------------------------------------------------------

	# Tenth, determine the bond length and angle for the new bond vector.
	bond_length = bond_lengths[donor_symbol][total_no_of_neighbours]
	bond_angle  = get_bond_angle(donor_symbol, total_no_of_neighbours, no_of_pairs_of_lone_electrons, element_to_attach=element_to_attach)

	# Eleventth, get the initial bonding vector
	angle = math_PI - bond_angle
	normal_direction = np.cos(angle) * normal_axis
	perp_direction = np.sin(angle) * perpendicular_unit_vector

	# Twelfth, determine the forces upon the hydrogen as the bond_vector is rotated around the normal_axis
	force_vs_angle = {}
	for angle_degrees in range(3600):
		angle_degrees = (angle_degrees/10.0)

		# 12.1: Obtain the rotation matrix for a point about the normal_axis
		rotation_matrix = get_rotation_matrix_around_arbitrary_axis(normal_axis, angle_degrees*(math_PI/180.0))

		# 12.2: Get the potential position for the hydrogen atom.
		bond_vector = normal_direction + (rotation_matrix @ perp_direction)
		potential_H_position = donor_position + bond_length * bond_vector

		# 12.3: Obtain the unit vectors and the distances from the hydrogen to the H-bond acceptor.
		unit_vectors_from_donor_to_acceptors = [get_unit_vector(acceptor_position-potential_H_position) for acceptor_positions in (acceptor_positions+hydrogens_bound_to_acceptors_positions)]
		distances_between_donor_and_acceptor = [get_distance(acceptor_position,potential_H_position) for acceptor_positions in (acceptor_positions+hydrogens_bound_to_acceptors_positions)]

		# 12.4: Get the force vectors upon the hydrogen atom from the H-bond acceptors.
		force_vectors = [((-1*acceptor_charge)/(distance_between_centre_and_acceptor ** 2.0))*unit_vector for acceptor_charge, distance_between_centre_and_acceptor, unit_vector in zip(acceptor_charges+hydrogens_bound_to_acceptors_charges, distances_between_donor_and_acceptor, unit_vectors_from_donor_to_acceptors)]

		# 12.5: Get the overall force vector upon the hydrogen atom.
		force_vector = sum(force_vectors)

		# 12.6 : Get the magnitude of the force.
		force_magnitude = np.linalg.norm(force_vector)

		# 12.7: Record the force vector at each angle.
		force_vs_angle[angle_degrees] = force_magnitude

	# Thirteenth, pick the angle that has the force with the largest magnitude
	max_amplitude_angle, max_amplitude = max(force_vs_angle.items(), key=lambda x:x[1])

	# Fourteenth, obtain the rotation matrix for a point about the normal_axis for the bonding vector to use
	rotation_matrix = get_rotation_matrix_around_arbitrary_axis(normal_axis, max_amplitude_angle*(math_PI/180.0))

	# Fifteenth, get the position for the hydrogen atom
	bond_vector = normal_direction + (rotation_matrix @ perp_direction)
	H_position = donor_position + bond_length * bond_vector

	# Sixteenth, get the bonding unit vector.
	new_bonding_unit_vector = get_unit_vector(H_position - donor_position)

	'''
	import matplotlib.pyplot as plt
	plt.plot(force_vs_angle.keys(), force_vs_angle.values())
	plt.show()
	import pdb; pdb.set_trace()
	raise Exception('Check this method is working.')
	'''

	# Seventeenth, return the newly created unit vector
	return new_bonding_unit_vector

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------



