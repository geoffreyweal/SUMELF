"""
check_for_H_pi_bonding.py, Geoffrey Weal, 6/7/2022

This script is designed to determine which hydrogens are involved in adjacent pi bonding.
"""
from math import pi
import numpy as np

from networkx import cycle_basis

from SUMELF.SUMELF.general_methods.geometry_methods          import planeFit
from SUMELF.SUMELF.general_methods.geometry_methods          import get_unit_vector
from SUMELF.SUMELF.general_methods.geometry_methods          import project_point_onto_plane
from SUMELF.SUMELF.general_methods.general_molecules_methods import get_hybridisation_from_CSD
from SUMELF.SUMELF.general_methods.distance_methods          import get_distance
from SUMELF.SUMELF.general_methods.unit_cell_methods         import get_cell_corner_points
from SUMELF.SUMELF.general_methods.general_molecules_methods import get_atomic_rings_in_ase_object

adjacent_atom_H_ring_normal_angle_tolerance = 45.0 * (pi/180.0)
def check_for_H_pi_bonding(copy_molecule, copy_molecule_graph, index_to_attach_Hs_to, crystal, crystal_graph, hybridisation_of_atoms_in_molecule):
	"""
	This method will check if their are any hydrogens that are near or adjacent to a conjugated ring.

	Parameters
	----------
	copied_molecule : ase.Atoms
		This is the molecule you want to check if its hydrogens are involved in hydrogen-bonding.
	molecule_graph : networkx.Graph
		This is the graph of the copied_molecule.
	index_to_attach_Hs_to : int 
		This is the index to check its neighbouring hydrogens for hydrogen-bonding
	crystal ; ase.Atoms
		This is the full crystal object.
	crystal_graph : networkx.Graph
		This is the networkx graph for the full crystal.
	hybridisation_of_atoms_in_molecule : list of str
		This list provides all the hybridisations given by the user (if available). s

	Returns
	-------
	To do.
	"""

	# First, determine the indices of hydrogens bound to index_to_attach_Hs_to
	bound_hydrogen_indices = [index for index in copy_molecule_graph[index_to_attach_Hs_to] if (copy_molecule[index].symbol in ['H', 'D'])]

	# Second, determine the position for index_to_attach_Hs_to
	index_to_attach_Hs_to_position = copy_molecule[index_to_attach_Hs_to].position

	#######################################################
	# PART 3: Check hydrogens that are adjacent to pi rings
	#######################################################

	# 3.1: Obtain all the rings inolved in the molecule of interest
	molecule_rings = [ring for ring in cycle_basis(copy_molecule_graph) if (len(ring) <= 7)]

	# 3.2: Determine all the rings that index_to_attach_Hs_to is bonded to (index_to_attach_Hs_to is not a hydrogen).
	molecule_aromatic_rings = get_atomic_rings_in_ase_object(molecule_rings, hybridisation_of_atoms_in_molecule)

	# 3.3: Determine all the aromatic rings that neighbour index_to_attach_Hs_to 
	molecule_aromatic_rings_that_neighbour_index_to_attach_Hs_to = []
	for ring in molecule_aromatic_rings:
		for ring_atom_index in copy_molecule_graph[index_to_attach_Hs_to]:
			if ring_atom_index in ring:
				molecule_aromatic_rings_that_neighbour_index_to_attach_Hs_to.append(ring)
				break

	# 3.4: Determine the hydrogens that are likely interacting with the pi cloud of aromatic rings
	hydrogens_that_have_interactions_with_pi_rings = []
	for neighbour_hydrogen_index in bound_hydrogen_indices:

		# 3.4.1: Obtain the position of the neighbouring hydrogen atom
		neighbour_hydrogen_position = copy_molecule[neighbour_hydrogen_index].position

		# 3.4.2: Determine the unit vector from the neihgbouring atom and its hydrogen
		neighbouring_atom_H_bond_unit_vector = get_unit_vector(index_to_attach_Hs_to_position - neighbour_hydrogen_position)

		# For each ring:
		for ring in molecule_aromatic_rings_that_neighbour_index_to_attach_Hs_to:

			# 3.4.1: Determine the positions for the atoms in the ring.
			positions_in_ring = [copy_molecule[ring_index].position for ring_index in ring]

			# 3.4.2: Determine the plane for the ring.
			centre_point, normal_vector = planeFit(positions_in_ring)

			# 3.4.3: We now want to project neighbouring_atom_H_bond_unit_vector onto the plane that is based on normal_vector and cross(normal_vector,neighbouring_atom_H_bond_unit_vector)
			perpendicular_unit_vector = get_unit_vector(np.cross(normal_vector, neighbouring_atom_H_bond_unit_vector))
			perp_centre_vector, perp_normal_vector = planeFit([np.array((0,0,0)), normal_vector, perpendicular_unit_vector])

			# 3.4.5: Obtain the neighbouring_atom_H_bond_unit_vector when projected onto the perpendicular plane.
			projected_vector_on_perp_plane_unit_vector = get_unit_vector(project_point_onto_plane(neighbouring_atom_H_bond_unit_vector, np.array((0,0,0)), perp_normal_vector, get_angle_between_vector_and_plane=False))

			# 3.4.6: Determine the dihedral angle between the H-atom bond and the 
			dihedral_angle = abs(np.arccos(np.dot(normal_vector, neighbouring_atom_H_bond_unit_vector)))
			if dihedral_angle > (pi/2.0):
				dihedral_angle = pi - dihedral_angle

			if dihedral_angle <= adjacent_atom_H_ring_normal_angle_tolerance:

				import pdb; pdb.set_trace()

	########################################################
	# PART 4: Check hydrogens that are nearby aromatic rings
	########################################################

	# 4.1: Wrap the molecule in the origin unit cell
	wrapped_copy_molecule = copy_molecule.copy()
	wrapped_copy_molecule.set_cell(crystal.get_cell())
	wrapped_copy_molecule.set_pbc(True)
	wrapped_copy_molecule.wrap()

	# 4.2: Get the list of hybridisations of atoms in the crystal
	hybridisation_of_atoms_in_crystal = crystal.arrays['Hybridisation']

	# 4.3: Obtain all the rings in the crystal
	crystal_rings = [ring for ring in cycle_basis(crystal_graph) if (len(ring) <= 7)]

	# 4.4: Obtain all the rings in the crystal that are aromatic
	crystal_aromatic_rings = get_atomic_rings_in_ase_object(crystal_rings, hybridisation_of_atoms_in_crystal)

	# 4.5: Determine the hydrogen atoms that are near pi rings in the crystal
	crystal_H_pi_ring_indices = []
	for neighbour_hydrogen_index in bound_hydrogen_indices:

		# 4.5.1: Obtain the position of the neighbouring hydrogen atom
		neighbour_hydrogen_position = wrapped_copy_molecule[neighbour_hydrogen_index].position

		# 4.5.2: Determine if this hydrogen is close to a pi ring
		if is_there_a_H_pi_interaction_between_crystal_and_H(neighbour_hydrogen_index, neighbour_hydrogen_position, crystal_aromatic_rings, crystal):
			crystal_H_pi_ring_indices.append(neighbour_hydrogen_index)

	if not len(hydrogens_that_have_interactions_with_pi_rings + crystal_H_pi_ring_indices) == 0:
		from ase.visualize import view
		view(copy_molecule)
		import pdb; pdb.set_trace()

# ---------------------------------------------------------------------------------------------------------------------------
		
H_pi_ring_distance_tolerance = 3.0 # A
T_H_above_ring_angle_tolerance = 45.0 * (pi/180.0)
def is_there_a_H_pi_interaction_between_crystal_and_H(neighbour_hydrogen_index, neighbour_hydrogen_position, crystal_aromatic_rings, crystal):
	"""
	This method will determine if their are any hydrogens that may be interacting with a pi ring.

	Parameters
	----------
	neighbour_hydrogen_position : np.array
		This is the position of the hydrogen atom of interest
	crystal_aromatic_rings : list of list of ints
		These are the indices of atoms in the crystal that are apart of the aromatic rings to focus on.
	crystal : ase.Atoms
		This is the crystal object.

	Returns
	-------
	aromatic_rings : list of lists of int.
		These are teh indices of the aromatic rings in the ase.atoms object.
	"""

	# First, determine the positions of the the cell corners for the crystal
	cell_corner_points = get_cell_corner_points(crystal_cell_lattice=crystal.get_cell())

	# Second, make a copy of the crystal and wrap it. Record the position of the wrapped crystal
	wrapped_crystal = crystal.copy()
	wrapped_crystal.wrap()
	crystal_wrapped_elements  = wrapped_crystal.get_chemical_symbols()
	crystal_wrapped_positions = wrapped_crystal.get_positions()

	# Third, look through all the aromatic rings in the crystal.
	for ring in crystal_aromatic_rings:

		# 3.1: Obtain the positions of the atoms involved in the aromatic ring
		aromatic_ring_position = [crystal[ring_atom_index].position for ring_atom_index in ring]

		# 3.2: Obtain the plane that best describes the aromatic ring.
		centre_position, normal_vector = planeFit(aromatic_ring_position)

		# 3.3: For each atom in the ring
		for ring_atom_index in ring:

			# 3.3.1: Get the position of the atom in the ring
			ring_atom_position = crystal_wrapped_positions[ring_atom_index]

			# 3.3.2: For each cell movement by one unit cell about the origin unit cell.
			for cell_corner_point in cell_corner_points:

				# 3.2.3: Get the position of the atom in one of the cells
				position_of_atom_in_a_cell = ring_atom_position + cell_corner_point

				# 3.3.4: Get the distance between the atom in the ring and the hydrogen. 
				distance = get_distance(position_of_atom_in_a_cell, neighbour_hydrogen_position)

				# 3.2.5: If the ring atom is within distance of the hydrogen atom
				if distance <= H_pi_ring_distance_tolerance:

					return False #True
					'''
					# 3.2.6: Get the angle between the hydrogen position and the plane of the ring, where the centre position is given as the position_of_atom_in_a_cell
					projected_point, _ = project_point_onto_plane(neighbour_hydrogen_position, centre_position, normal_vector, get_angle_between_vector_and_plane=True)
					vec1 = get_unit_vector(neighbour_hydrogen_position - position_of_atom_in_a_cell)
					vec2 = get_unit_vector(projected_point             - position_of_atom_in_a_cell)
					angle_between_point_and_plane_from_ring_atom = np.arccos(np.dot(vec1,vec2))

					# Get the angle between the ring plane's normal and the hydrogen position, where the centre position is given as the position_of_atom_in_a_cell
					angle_between_point_and_plane_normal_from_ring_atom = (pi/2.0) - abs(angle_between_point_and_plane_from_ring_atom)

					if angle_between_point_and_plane_normal_from_ring_atom > (pi/2.0):
						angle_between_point_and_plane_normal_from_ring_atom = pi - angle_between_point_and_plane_normal_from_ring_atom

					import pdb; pdb.set_trace()

					if angle_between_point_and_plane_normal_from_ring_atom <= T_H_above_ring_angle_tolerance:
						return True
					'''
					
	return False

# ---------------------------------------------------------------------------------------------------------------------------
	




