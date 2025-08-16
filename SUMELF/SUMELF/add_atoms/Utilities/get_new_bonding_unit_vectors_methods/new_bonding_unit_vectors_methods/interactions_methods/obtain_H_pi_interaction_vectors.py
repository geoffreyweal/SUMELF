"""
obtain_H_pi_interaction_vectors.py, Geoffrey Weal, 6/7/2022

This script is designed to identify new H-pi interactions with aromatic rings.
"""
from math import pi
import numpy as np
import logging

#from shapely.geometry import Point, Polygon

from networkx import cycle_basis

from SUMELF.SUMELF.general_methods.geometry_methods                      import planeFit
from SUMELF.SUMELF.general_methods.geometry_methods                      import get_unit_vector
from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods import get_bond_angle
from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods import get_bond_lengths
from SUMELF.SUMELF.general_methods.general_molecules_methods             import get_atomic_rings_in_ase_object

def obtain_H_pi_interaction_vectors(molecule, molecule_graph, index_to_attach_Hs_to, total_no_of_neighbours, no_of_pairs_of_lone_electrons, crystal, crystal_graph, element_to_attach=None):
	"""
	This method is designed to identify new H-pi interactions with aromatic rings.

	Parameters
	----------
	molecule : ase.Atoms
		This is the molecule you want to check if its hydrogens are involved in hydrogen-bonding.
	molecule_graph : networkx.Graph
		This is the graph of the molecule.
	index_to_attach_Hs_to : int 
		This is the index to check its neighbouring hydrogens for hydrogen-bonding
	total_no_of_neighbours : int
		This is the total number of neighbours that we want to surround the central atom
	no_of_pairs_of_lone_electrons : int
		This is the number of lone pairs of electrons surrounding the central atom.
	crystal : ase.Atoms
		This is the full crystal object.
	crystal_graph : networkx.Graph
		This is the networkx graph for the full crystal.
	element_to_attach : str.
		This is the element you are using this method to attach a new atom to your molecule. If None given, their will be an issue given.

	Returns
	-------
	new_bonding_unit_vectors : list of numpy.array
		These are all the new bonding unit vector for adding missing hydrogens to the system.
	"""

	# Before beginning, get the bond length betwen atoms and element_to_attach
	bond_lengths = get_bond_lengths(element_to_attach=element_to_attach)

	# First, obtain the element for the central atom.
	central_atom_element = molecule[index_to_attach_Hs_to].symbol

	# Second, determine the bond able between atoms about the central atom.
	bond_angle = get_bond_angle(central_atom_element, total_no_of_neighbours, no_of_pairs_of_lone_electrons=no_of_pairs_of_lone_electrons, element_to_attach=element_to_attach)

	# Third, determine the position for index_to_attach_Hs_to
	index_to_attach_Hs_to_position = molecule[index_to_attach_Hs_to].position

	#################################################################################
	# PART 4: Obtain hydrogen bonding vectors for atoms that are adjacent to pi rings
	#################################################################################

	# 4.1: Obtain all the rings inolved in the molecule of interest
	molecule_rings = [ring for ring in cycle_basis(molecule_graph) if (len(ring) <= 7)]

	# 4.2: Determine all the rings that index_to_attach_Hs_to is bonded to (index_to_attach_Hs_to is not a hydrogen).
	# This is currently commented out, as for vcarious reason we want this algorithm to work on aromatic rings as well as aliphatic rings (due to sterics in aliphatic rings).
	# molecule_aromatic_rings = get_atomic_rings_in_ase_object(molecule_rings, hybridisation_of_atoms_in_molecule, molecule)

	# 4.3: Determine all the aromatic rings that neighbour index_to_attach_Hs_to 
	molecule_aromatic_rings_that_neighbour_index_to_attach_Hs_to = []
	for ring in molecule_rings:
		for ring_atom_index in molecule_graph[index_to_attach_Hs_to]:
			if ring_atom_index in ring:
				molecule_aromatic_rings_that_neighbour_index_to_attach_Hs_to.append((ring, ring_atom_index))
				break

	# 4.4: If this neighbouring carbon is bound to more than one aromatic ring, don't worry about using this method, as their are probably pi-pi effects that are more important than this effect
	if len(molecule_aromatic_rings_that_neighbour_index_to_attach_Hs_to) >= 2:
		raise Exception('Need to figure out what to do at this point, using an example.')

	# 4.5: For each ring:
	new_bonding_unit_vectors = []
	for ring, ring_atom_bonded_to_adjacent_atom_index in molecule_aromatic_rings_that_neighbour_index_to_attach_Hs_to:

		# 4.5.1: Determine the positions for the atoms in the ring.
		positions_in_ring = [molecule[ring_index].position for ring_index in ring]

		# 4.5.2: Determine the plane for the ring.
		centre_point, normal_vector = planeFit(positions_in_ring)
		normal_vector = get_unit_vector(normal_vector)

		# 4.5.3: Get the unit vector that point along the bond from the atom in the ring to the index_to_attach_Hs_to atom.
		parallel_vector = get_unit_vector(index_to_attach_Hs_to_position - molecule[ring_atom_bonded_to_adjacent_atom_index].position)

		# 4.5.4: Determine the new bonding vector for this hydrogen.
		angle_from_normal = bond_angle - (pi/2.0)
		new_bonding_vector = np.sin(angle_from_normal)*parallel_vector + np.cos(angle_from_normal)*normal_vector

		# 4.5.5: Get the unit vector for the new_bonding_vector
		new_bonding_unit_vector = get_unit_vector(new_bonding_vector)

		# 4.5.6: Append new_bonding_unit_vector to new_bonding_unit_vectors
		new_bonding_unit_vectors.append(new_bonding_unit_vector)

	#################################################################################

	# Fifth, return new_bonding_unit_vectors
	return new_bonding_unit_vectors

# ---------------------------------------------------------------------------------------------------------------------------

