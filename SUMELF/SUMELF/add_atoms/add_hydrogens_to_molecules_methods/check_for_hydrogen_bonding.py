"""
check_for_hydrogen_bonding.py, Geoffrey Weal, 6/7/2022

This script is designed to determine which hydrogens are involved in hydrogen bonding.
"""
from SUMELF.SUMELF.general_methods.unit_cell_methods import get_cell_corner_points
from SUMELF.SUMELF.general_methods.angle_methods     import get_angle
from SUMELF.SUMELF.general_methods.geometry_methods  import get_unit_vector
from SUMELF.SUMELF.general_methods.distance_methods  import get_distance

hydrogen_bond_donor    = ['N', 'O', 'F'] # 'C',
hydrogen_bond_acceptor = ['N', 'O', 'F']
distance_tolerance = 0.1
def check_for_hydrogen_bonding(copied_molecule, copy_molecule_graph, index_to_attach_Hs_to, crystal):
	"""
	This method is designed to determine hydrogens in the molecule that could be inolved in hydrogen bonding.

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

	Returns
	-------
	To do.
	"""

	# First, only continue with this method if index_to_attach_Hs_to is a hydrogen_bond_donor atom.
	if copied_molecule[index_to_attach_Hs_to].symbol not in hydrogen_bond_donor:
		return 

	# Second, determine the positions of the the cell corners for the crystal
	cell_corner_points = get_cell_corner_points(crystal_cell_lattice=crystal.get_cell())

	# Third, make a copy of the crystal and wrap it.Record the position of the wrapped crystal
	wrapped_crystal = crystal.copy()
	wrapped_crystal.wrap()
	crystal_wrapped_elements  = wrapped_crystal.get_chemical_symbols()
	crystal_wrapped_positions = wrapped_crystal.get_positions()

	# Fourth, get the position for the central atom.
	central_atom_position = copied_molecule[index_to_attach_Hs_to].position

	# Fifth, get the positinos of the neighbouring atoms to the central atom.
	neighbouring_atom_positions = [copied_molecule[index].position for index in copy_molecule_graph[index_to_attach_Hs_to]]

	# Sixth, determining the indices of the hydrogen that are attached to index_to_attach_Hs_to
	hydrogens_bound_to_index_to_attach_Hs_to = [index for index in copy_molecule_graph[index_to_attach_Hs_to] if copied_molecule[index].symbol == 'H']

	# Seventh, determine the hydrogen atoms likely inolved in hydrogen bonding in the molecule. 
	hydrogen_bonding_atoms = {}
	molecule_positions = copied_molecule.get_positions()
	for H_bound_index in hydrogens_bound_to_index_to_attach_Hs_to:

		# 7.1: Obtain the position of the bound hydrogen of interest.
		bound_H_position = molecule_positions[H_bound_index]

		# 7.2: Look through all the atoms in the crystal
		for third_index in range(len(wrapped_crystal)):

			# 7.2.1: If the third atom is a hydrogen, move on to the nex
			if crystal_wrapped_elements[third_index] not in hydrogen_bond_acceptor:
				continue 

			# 7.2.2: obtain the position of this third atom. 
			original_third_atom_position = crystal_wrapped_positions[third_index]

			# 7.2.3: Try out all the positions from trnslations about the original unit cell by one unit cell length in each i,j,k direction. 
			for corner_position in cell_corner_points:

				# 7.2.3.1: Get the third_atom_position after translation.
				third_atom_position = original_third_atom_position + corner_position

				# 7.2.3.2: If the third atom is index_to_attach_Hs_to or a atom attached to index_to_attach_Hs_to, (i.e., atom 1 or atom 2), continue on
				if any([(get_distance(third_atom_position, position) <= distance_tolerance) for position in ([central_atom_position] + neighbouring_atom_positions)]):
					continue

				# 7.2.3.3: If the hydrogen and the third atom have the right distance and angle, this is a potential hydrogen bonding pair
				if (get_distance(bound_H_position, third_atom_position) <= 3.5) and (170.0 <= get_angle(central_atom_position, bound_H_position, third_atom_position) <= 190.0):
					hydrogen_bonding_atoms.setdefault(H_bound_index,[]).append((central_atom_position, third_atom_position))



	if not len(hydrogen_bonding_atoms) == 0:
		from ase.visualize import view
		import pdb; pdb.set_trace()






