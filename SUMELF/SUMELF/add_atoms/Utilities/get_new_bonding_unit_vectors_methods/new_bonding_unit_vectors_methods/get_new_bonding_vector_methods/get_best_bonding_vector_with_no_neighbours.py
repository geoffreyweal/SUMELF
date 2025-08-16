"""
supporting_methods.py, Geoffrey Weal, 5/7/2022

This script contains methods used by the python scripts in this folder.
"""
import numpy as np

from SUMELF.SUMELF.general_methods.geometry_methods                                                                                                            import get_unit_vector
from SUMELF.SUMELF.general_methods.distance_methods                                                                                                            import get_distance
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors_methods.new_bonding_unit_vectors_methods.get_new_bonding_vector_methods.supporting_methods import record_hydrogen_positions

same_atom_distance_tolerance = 0.05
def get_best_bonding_vector_with_no_neighbours(molecule, donor_index, wrapped_crystal, acceptor_info, element_to_attach=None):
	"""
	This method will provide a optimal bonding unit vector for hydrogen bonding system, where the H-bond donor is bound to no atoms.

	Parameters
	----------
	molecule : ase.Atoms
		This is the wrapped molecule of interest in the crystal.
	donor_index : int
		This is the index of the H-bond donor.
	wrapped_crystal : ase.Atoms
		This is the crystal object.
	acceptor_info : list of ints
		These are the indices of all the H-bond acceptors in the crystal, and their positions.
	element_to_attach : str.
		This is the element you are using this method to attach a new atom to your molecule. If None given, their will be an issue given.

	Returns
	-------
	new_bonding_unit_vector : numpy.array
		This is the newly created bonding unit vector.
	"""

	# First, obtain the position of the H-bond donor.
	donor_position = molecule[donor_index].position

	# Second, obtain the charges and the positions of the nearby hydrogen bonding acceptors.
	acceptor_charges = [round(wrapped_crystal[acceptor_index].charge,5) for acceptor_index, _, _ in acceptor_info]
	acceptor_charges = [(acceptor_charge - 0.5) for acceptor_charge in acceptor_charges]
	acceptor_positions = [acceptor_position for _, acceptor_position, _ in acceptor_info]

	# Third, obtain the charges and the positions of any hydrogens that are bound to ydrogen bonding acceptors.
	hydrogens_bound_to_acceptors_positions = record_hydrogen_positions(acceptor_info)
	hydrogens_bound_to_acceptors_charges = [+0.5 for hydrogen_position in hydrogens_bound_to_acceptors_positions]
	
	# Fourth, get the force vectors upon the H-bond donor atom from the H-bond acceptors
	force_vector_of_acceptor_upon_donor = np.array((0.,0.,0.))
	for acceptor_charge, acceptor_position in zip(acceptor_charges+hydrogens_bound_to_acceptors_charges, acceptor_positions+hydrogens_bound_to_acceptors_positions):
		if get_distance(acceptor_position, donor_position) <= same_atom_distance_tolerance:
			continue
		unit_vector = get_unit_vector(acceptor_position-donor_position)
		distance = get_distance(acceptor_position,donor_position)
		force_vector_of_acceptor_upon_donor += (-acceptor_charge/(distance**2.0)) * unit_vector

	if not np.any(force_vector_of_acceptor_upon_donor):
		raise Exception('Error, Huh?')

	# Fifth, get the bonding unit vector, which we will assume is close to being the same vector as the force vector upon the donor atom by all other acceptor atoms
	new_bonding_unit_vector = get_unit_vector(force_vector_of_acceptor_upon_donor)

	# Sixth, return the newly created unit vector
	return new_bonding_unit_vector

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------




