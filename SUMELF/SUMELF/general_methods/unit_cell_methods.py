"""
unit_cell_methods.py, Geoffrey Weal, 19/7/22

This script is designed to hold methods useful for determining features of unit cells. 
"""
import numpy as np
from itertools import product
from scipy.spatial import Delaunay

from ase import Atoms

from SUMELF.SUMELF.general_methods.remove_hydrogens import remove_hydrogens
from SUMELF.SUMELF.general_methods.distance_methods import get_distance

# ---------------------------------------------------------------------------------------------------------------------------

def get_cell_corner_points(crystal_cell_lattice, super_cell_reach=1, bottom_of_range=None, get_corrspeonding_ijk_values=False): # get_unit_cell_displacements --> get_corrspeonding_ijk_values
	"""
	This method gives the corner point of super cells surrounding a unit cell. 

	This method helps with moving atoms to positions outside the original unit cell in order to construct a connecting molecule.

	Parameters
	----------
	crystal_cell_lattice : np.array
		A 3x3 row matrix giving the cell vectors in the i, j, and k basis set directions (in cartesian xyz coordinates). Given in Å. 
	super_cell_reach : int.
		This parameter indicates how big you want the super cell to be. Default: 1.
	bottom_of_range : int.
		This parameter allows you to change the bottom of the range of the super cell. If set to None, bottom_of_range = -super_cell_reach. Default: None.
	get_corrspeonding_ijk_values : bool.
		If True, return the ijk cells as well as the cell points in Å. If False, only return the cell points in Å.

	Returns
	-------
	cell_points: list
		This is a list of all the distances between identical points in each ijk section of the supercell from the origin unit cell. Given in Å. 
	cells : list
		This is a list of all the ijk representations for the points given in cell_points
	"""

	# First, determine the ijk coordinates to range over.
	top_of_range = super_cell_reach
	if bottom_of_range is None:
		bottom_of_range = -super_cell_reach
	range_to_scan_over = tuple(range(bottom_of_range,top_of_range+1,1))

	# Second, determine ijk coordinates that you want to include in your super cell.
	cells = [range_to_scan_over for _ in range(3)]
	cells = sorted(list(product(*cells)))

	# 2.1, move the (0,0,0) unit to the front of the list.
	if (0,0,0) in cells:
		cells.remove((0,0,0))
		cells.insert(0, (0,0,0))

	# Third, convert the ijk coordinates to xyz coordinates
	cell_points = []
	for subcell in cells:
		displacement = np.matmul(np.array(subcell),crystal_cell_lattice)
		'''
		if add_slight_shortening:
			for index in range(len(subcell)):
				if subcell[index] == super_cell_reach:
					displacement[index] -= 0.00000001 # This is so that the point is just on the edge of the next cell
		'''
		cell_points.append(displacement)

	# Fourth, return the cell points and corresponding ijk values (if desired).
	if get_corrspeonding_ijk_values:
		return cell_points, cells
	else:
		return cell_points

# ---------------------------------------------------------------------------------------------------------------------------

zero_displacement = np.array((0,0,0))
def centre_molecule_in_cell(connected_molecule, crystal_cell_lattice, move_molecule=True, dimer_index=28):
	"""
	This method will centre the molecule as close to the middle of the unit cell as possible while maintain its pbc consistent coordinates.

	You can also choose not to move the molecule if all you need to know is the displacement that is needed to move the molecule.

	Parameters
	----------
	connected_molecule : ase.Atoms
		This is the molecule you want to centre as close to the middle of the unit cell as possible.
	crystal_cell_lattice : np.array
		The lattice vectors along the i, j, and k directions (in cartesian xyz coordinates). Given in Å. 
	move_molecule: bool.
		This tag indicates if you want to move connected_molecule as close to the middle of the unit cell as possible while maintain its pbc consistent coordinates. Default: True.
	dimer_index : int
		This is the index of the dimer. This was used for debugging a floating point error using ABIBEW (and for MUPMOC) from the CCDC. Can remove when happy. 

	Returns
	-------
	displacement_molecule_to_centre_molecule_close_to_centre_of_cell: np.array
		This is a 1x3 row vector that gives the displacement required to move the molecule as close to the middle of the unit cell as possible while maintain its pbc consistent coordinates. Given in Å. 
	"""

	# First, get the centre position of the origin unit cell.
	origin_cell_cell_points = get_cell_corner_points(crystal_cell_lattice,super_cell_reach=1,bottom_of_range=0)
	origin_cell_centre_of_unit_cell = sum(origin_cell_cell_points)/float(len(origin_cell_cell_points))
	origin_cell_object = Delaunay(origin_cell_cell_points)

	# Second, get all the displacements that surround the cell around the molecule.
	cell_points = get_cell_corner_points(crystal_cell_lattice,super_cell_reach=1)

	# ---------------------------------------------------------------------------------------------------------------------
	# Third, determine the displacement that is required to move the molecule as close to the centre of the origin unit cell as possible.
	# Will center the molecule excluding the hydrogen atoms.

	# 3.1: Obtain a dummy molecule that has had its hydrogens removed.
	dummy_molecule = remove_hydrogens(connected_molecule.copy())

	# 3.2: Note the values of the previous details.
	displacement_molecule_to_centre_molecule_close_to_centre_of_cell = np.array((0.0,0.0,0.0))
	previous_number_of_atoms_inside_unit_cell = 0
	previous_distance_from_centre = float('inf')
	previous_distance_from_zero_point = float('inf')
	previous_com_y_position = float('inf')
	previous_com_x_position = float('inf')
	previous_counter = 0

	# 3.3: Note the values of the most current details.
	current_number_of_atoms_inside_unit_cell = 0
	current_distance_from_centre = float('inf') 
	current_distance_from_zero_point = float('inf') 
	current_com_y_position = float('inf')
	current_com_x_position = float('inf')

	#if dimer_index == 28:
	#	print(f'looking to move com of {dimer_index}')

	while True:

		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		# 3.3: Determine which displacement in cell_points will move the molecule closer to the origin cells centre of cell position
		single_cell_displacement = np.array((0,0,0))
		centre_of_mass = dummy_molecule.get_center_of_mass()

		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		# 3.4: Try all displacement vectors on the dummy molecule and determine
		#      * The number of atoms of that molecule inside the origin unit cell.
		#      * The distance of the molecules centre of mass to the centre of the origin unit cell.
		#      * The distance of the molecules centre of mass to the (0.,0.,0.) point.
		for displacement in cell_points:

			#if dimer_index == 28:
			#	print(f'displacement = {displacement}')

			moved_com_position = centre_of_mass + displacement
			number_of_atoms_inside_unit_cell = number_of_points_in_cell(dummy_molecule, origin_cell_object, displacement)
			distance_to_centre = round(get_distance(origin_cell_centre_of_unit_cell,moved_com_position),4)
			distance_to_zero_point = round(get_distance(zero_displacement,moved_com_position),4)
			com_y_position = moved_com_position[1]
			com_x_position = moved_com_position[0]

			# Determine if their are more atoms in the unit cell
			if number_of_atoms_inside_unit_cell > current_number_of_atoms_inside_unit_cell:

				#if dimer_index == 28:
				#	print('1: number_of_atoms_inside_unit_cell > current_number_of_atoms_inside_unit_cell')
				#	import pdb; pdb.set_trace()

				# If 
				#	1. There are more points inside this cell, 
				#      --------------------------------------
				# then,
				# --> This is where to move the molecule to currently. Record it.

				single_cell_displacement                 = displacement
				current_number_of_atoms_inside_unit_cell = number_of_atoms_inside_unit_cell
				current_distance_from_centre             = distance_to_centre
				current_distance_from_zero_point         = distance_to_zero_point
				current_com_y_position                   = com_y_position
				current_com_x_position                   = com_x_position

			elif number_of_atoms_inside_unit_cell == current_number_of_atoms_inside_unit_cell:

				#if dimer_index == 28:
				#	print('2: number_of_atoms_inside_unit_cell == current_number_of_atoms_inside_unit_cell')

				# If:
				#	1. The number of atoms are the same, 
				#      --------------------------------
				# then,
				# --> Consider which position is closer to the centre of the cell.
				# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
				if distance_to_centre < current_distance_from_centre:

					#if dimer_index == 28:
					#	print('3: distance_to_centre < current_distance_from_centre')
					#	import pdb; pdb.set_trace()

					# If: 
					#   1. The number of atoms are the same, and
					#   2. This position is closer to the centre of the cell, 
					#      -------------------------------------------------
					# then
					# --> This is where to move the molecule to currently. Record it.

					single_cell_displacement                 = displacement
					current_number_of_atoms_inside_unit_cell = number_of_atoms_inside_unit_cell
					current_distance_from_centre             = distance_to_centre
					current_distance_from_zero_point         = distance_to_zero_point
					current_com_y_position                   = com_y_position
					current_com_x_position                   = com_x_position

				elif distance_to_centre == current_distance_from_centre:

					#if dimer_index == 28:
					#	print('4: distance_to_centre == current_distance_from_centre')
					#	import pdb; pdb.set_trace()

					# If: 
					#   1. The number of atoms are the same, and 
					#   2. If the displacement to the centre of the unit cell is the same, 
					#      ---------------------------------------------------------------
					# then,
					# --> Determine if this displacement is closer to the (0.,0.,0.) point.
					# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
					if distance_to_zero_point < current_distance_from_zero_point:

						#if dimer_index == 28:
						#	print('5: distance_to_zero_point < current_distance_from_zero_point')

						# If: 
						#   1. The number of atoms are the same, and
						#   2. This position to the centre of the unit cell is the same, 
						#   3. This position is closer to the (0.,0.,0.) point,
						#      -----------------------------------------------
						# then, 
						# --> This is where to move the molecule to currently. Record it.

						single_cell_displacement = displacement
						current_number_of_atoms_inside_unit_cell = number_of_atoms_inside_unit_cell
						current_distance_from_centre = distance_to_centre
						current_distance_from_zero_point = distance_to_zero_point
						current_com_y_position                   = com_y_position
						current_com_x_position                   = com_x_position

					elif distance_to_zero_point == current_distance_from_zero_point:

						#if dimer_index == 28:
						#	print('6: distance_to_zero_point == current_distance_from_zero_point')

						# If: 
						#   1. The number of atoms are the same, 
						#   2. If the displacement to the centre of the unit cell is the same, and
						#   3. If the displacement to the (0.,0.,0.) point is the same, 
						#      -------------------------------------------------------
						# then,
						# --> Determine if the y position of the centre of mass is closer to -float('inf')
						# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
						if com_y_position < current_com_y_position:

							#if dimer_index == 28:
							#	print('7: com_y_position < current_com_y_position')

							# If: 
							#   1. The number of atoms are the same,  
							#   2. If the displacement to the centre of the unit cell is the same, 
							#   3. If the displacement to the (0.,0.,0.) point is the same, and
							#   4. If this centre of mass is closer to -float('inf') in the y position, 
							#      -------------------------------------------------------------------
							# then,
							# --> This is where to move the molecule to currently. Record it.

							single_cell_displacement = displacement
							current_number_of_atoms_inside_unit_cell = number_of_atoms_inside_unit_cell
							current_distance_from_centre = distance_to_centre
							current_distance_from_zero_point = distance_to_zero_point
							current_com_y_position                   = com_y_position
							current_com_x_position                   = com_x_position

						elif com_y_position == current_com_y_position:

							#if dimer_index == 28:
							#	print('8: com_y_position == current_com_y_position')

							# If: 
							#   1. The number of atoms are the same, 
							#   2. If the displacement to the centre of the unit cell is the same,
							#   3. If the displacement to the (0.,0.,0.) point is the same, and
							#   4. If the centre of mass is the same distance from y = -float('inf'),
							#      -----------------------------------------------------------------
							# then,
							# --> Determine if the x position of the centre of mass is closer to -float('inf')
							# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
							if com_x_position < current_com_x_position:

								#if dimer_index == 28:
								#	print('9: com_x_position < current_com_x_position')

								# If: 
								#   1. The number of atoms are the same, 
								#   2. If the displacement to the centre of the unit cell is the same,
								#   3. If the displacement to the (0.,0.,0.) point is the same, 
								#   4. If the centre of mass is the same distance from y = -float('inf'), and
								#   5. If this centre of mass is closer to -float('inf') in the x position, 
								#      -----------------------------------------------------------------
								# then,
								# --> This is where to move the molecule to currently. Record it.

								single_cell_displacement = displacement
								current_number_of_atoms_inside_unit_cell = number_of_atoms_inside_unit_cell
								current_distance_from_centre = distance_to_centre
								current_distance_from_zero_point = distance_to_zero_point
								current_com_y_position                   = com_y_position
								current_com_x_position                   = com_x_position

							'''
							elif com_x_position == current_com_x_position:

								# If: 
								#   1. The number of atoms are the same, 
								#   2. If the displacement to the centre of the unit cell is the same,
								#   3. If the displacement to the (0.,0.,0.) point is the same, 
								#   4. If the centre of mass is the same distance from y = -float('inf'), and
								#   5. If the centre of mass is the same distance from x = -float('inf'),
								#      -----------------------------------------------------------------
								# then,
								# --> The old position and the new position are the same??? 
								#     * There is a programming issue

								to_string  = 'Error: The old and new movement that are being analysed in the "centre_molecule_in_cell" method are the same.\n'
								to_string += 'This should never happen, and indicate a programming error.\n'
								to_string += 'Check this.'
								print(to_string)
								import pdb; pdb.set_trace()
								raise Exception(to_string)
							'''

		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		# 3.5: Move the dummy molecule and determine if it has been placed as much as possible in the unit cell.
		if not (single_cell_displacement == zero_displacement).all():

			# 3.5.1: Move the dummy molecule
			dummy_molecule.set_positions(dummy_molecule.get_positions() + single_cell_displacement)

			# 3.5.2: Update the displacement to the centre of the cell
			displacement_molecule_to_centre_molecule_close_to_centre_of_cell += single_cell_displacement

			# 3.5.3: Perform boolean checks to figure out if the current position is equivalent to the previous position. 
			are_same_number_of_atoms_inside_unit_cell = current_number_of_atoms_inside_unit_cell == previous_number_of_atoms_inside_unit_cell
			are_same_distance_from_centre             = current_distance_from_centre             == previous_distance_from_centre
			are_same_distance_from_zero_point         = current_distance_from_zero_point         == previous_distance_from_zero_point
			are_same_com_y_position                   = current_com_y_position                   == previous_com_y_position
			#are_same_com_x_position                   = current_com_x_position                   == previous_com_x_position

			#if dimer_index == 28:
			#	print(f'Moving the com of {dimer_index} to {displacement_molecule_to_centre_molecule_close_to_centre_of_cell}')

			# 3.5.4: Determine if this is a new better position for the molecule to move to.
			if are_same_number_of_atoms_inside_unit_cell and are_same_distance_from_centre and are_same_distance_from_zero_point and are_same_com_y_position: # and are_same_com_x_position:
				previous_counter += 1
			else:
				previous_counter = 0
				previous_number_of_atoms_inside_unit_cell = current_number_of_atoms_inside_unit_cell
				previous_distance_from_centre             = current_distance_from_centre
				previous_distance_from_zero_point         = current_distance_from_zero_point
				previous_com_y_position                   = current_com_y_position
				previous_com_x_position                   = current_com_x_position
			
			# 3.5.2: If identically good positions have been found too much, just pick one and move on, because both are as good as each other
			if previous_counter == 6:
				raise Exception('Error: Got here. I want to check why. Geoff 21/5/24')
				break

		else:

			#if dimer_index == 28:
			#	print(f'Final Movement of com of {dimer_index} at {displacement_molecule_to_centre_molecule_close_to_centre_of_cell}')

			break

		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# ---------------------------------------------------------------------------------------------------------------------

	# Fourth, if you want to move connected_molecule close to the centre of the cell, do it.
	if move_molecule:
		connected_molecule.set_positions(connected_molecule.get_positions() + displacement_molecule_to_centre_molecule_close_to_centre_of_cell)

	#if dimer_index == 28:
	#	import pdb; pdb.set_trace()

	# Fifth, return the displacement required to move your molecule as close to the centre of the unit cell as possible.
	return displacement_molecule_to_centre_molecule_close_to_centre_of_cell

def number_of_points_in_cell(dummy_molecule, origin_cell_object, displacement):
	"""
	This method will allow us to determine if the atoms in the dummy molecule are within the unit cell or not.
	
	This method is primarily used by the centre_molecule_in_cell method above.

	This method was written with advice from https://stackoverflow.com/questions/29311682/finding-if-point-is-in-3d-poly-in-python
	
	Parameters
	----------
	dummy_molecule : ase.Atoms
		This is the molecule you want to see the number of atoms that are inside the unit cell.
	origin_cell_object : np.spatial.Delaunay
		This object describes the unit cell volume.
	displacement: np.array
		This is the extra displacement to add to the positions of every atom in the dummy_molecule .

	Returns
	-------
	number_of_atoms_inside_unit_cell: int
		This is the number of atom from the dummy molecule that are inside the unit cell described by origin_cell_object
	"""
	# First, set the number_of_atoms_inside_unit_cell to 0
	number_of_atoms_inside_unit_cell = 0

	# Second, get the positinons in atoms given in dummy_molecule.
	if isinstance(dummy_molecule, Atoms):
		positions = dummy_molecule.get_positions()
	else:
		positions = dummy_molecule

	# Third, determine which atoms are inside and outside of the cell object
	for position in positions:
		# For each atom in the molecule, is in inside or outside of origin_cell_object
		if origin_cell_object.find_simplex(position + displacement) >= 0:
			number_of_atoms_inside_unit_cell += 1
	return number_of_atoms_inside_unit_cell

# ---------------------------------------------------------------------------------------------------------------------------
