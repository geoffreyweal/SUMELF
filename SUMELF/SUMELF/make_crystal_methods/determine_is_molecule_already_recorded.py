"""
determine_is_molecule_already_recorded.py, Geoffrey Weal, 29/5/22

This script methods for determining if a molecule has already been included in a list of molecules. 
"""

from SUMELF.SUMELF.general_methods.distance_methods  import get_distance
from SUMELF.SUMELF.general_methods.unit_cell_methods import get_cell_corner_points
from SUMELF.SUMELF.general_methods.remove_hydrogens  import remove_hydrogens

def determine_is_molecule_already_recorded(original_new_molecule, obtained_molecules, crystal_cell_lattice, super_cell_reach, include_hydrogen, consider_elements, return_same_molecule_details=False):
	"""
	This method will determine if you have already recorded this molecule in the crystal. 

	Parameters
	----------
	original_new_molecule : ase.Atoms
		This is the new molecule that has just been obtained from a crystal symmetry operation.
	obtained_molecules : list of ase.Atoms
		This is a list of all the unique ase.Atoms objects in the crystal that have currently been recorded.
	crystal_cell_lattice : numpy.array
		This is the unit cell vectors.
	super_cell_reach

	include_hydrogen

	consider_elements : bool.
		If True, consider the element of the atom when comparing atoms. If False, don't consider the element of the atom when comparing atoms. Default: True.
	return_same_molecule_details :


	Returns
	-------
	True if new_molecule is already in obtained_molecules, False if not.
	"""

	# First, get all the displacements that surround the cell around the molecule.
	cell_points = get_cell_corner_points(crystal_cell_lattice, super_cell_reach=super_cell_reach)

	# Second, make a copy of new_molecule
	new_molecule = original_new_molecule.copy()

	# Third, get a copy of new_molecule with its hydrogens removed. This is required for examples like CIWBOC.
	if not include_hydrogen:
		new_molecule = remove_hydrogens(new_molecule)

	# Fourth, look through all the obtain molecules to see if any are the same as new_molecule
	for obtained_molecule_name, original_obtained_molecule in sorted(obtained_molecules.items(), reverse=True):

		# 4.1: Make a copy of obtained_molecule
		obtained_molecule = original_obtained_molecule.copy()

		# 4.2: Get a copy of obtained_molecule with its hydrogens removed.
		if not include_hydrogen:
			obtained_molecule = remove_hydrogens(obtained_molecule)

		# 4.3: If new_molecule does not have the same chemical formula as obtained_molecule, move on
		if consider_elements:
			if not (sorted(new_molecule.get_chemical_symbols()) == sorted(obtained_molecule.get_chemical_symbols())):
				continue

		# 4.4: The obtained_molecule can be in any translated position given by a cell_point in cell_points. Need to try them all out.
		for cell_point in cell_points:

			# 4.5: Copy the obtained_molecule
			obtained_molecule_at_cell_point = obtained_molecule.copy()

			# 4.6: Place obtained_molecule_at_cell_point in the translated position
			obtained_molecule_at_cell_point.set_positions(obtained_molecule_at_cell_point.get_positions() + cell_point)

			# 4.7: Are new_molecule and obtained_molecule_at_cell_point the same molecule in the same place.
			if same_molecules(new_molecule, obtained_molecule_at_cell_point, consider_elements=consider_elements):
				
				# 4.8: Found a molecule in obtained_molecules that is the same as new_molecule, so return True
				if return_same_molecule_details:
					return True, obtained_molecule_name, cell_point
				else:
					return True

	# Fifth, if new_molecule is not the same as any obtained_molecule in obtained_molecules, even when translated by any of the cell_point.
	if return_same_molecule_details:
		return False, None, None
	else:
		return False

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

max_diff_distance = 0.0001
def same_molecules(new_molecule, obtained_molecule, consider_elements=True):
	"""
	This method will determine if these two molecules are the same molecules, particularly spatial so that two molecules are not on top of each other.

	Note: If consider_elements, elements lists will still be manipulated by this method, but it wont be of any concerquence and will not change the elements of new_molecule and obtained_molecule

	Parameters
	----------
	new_molecule : ase.atoms
		This is the new molecule that has just been obtained from a crystal symmetry operation.
	obtained_molecule : ase.Atoms
		This is one of the already obtained molecules that is currently in the obtained_molecules list
	consider_elements : bool.
		If True, consider the element of the atom when comparing atoms. If False, don't consider the element of the atom when comparing atoms. Default: True

	Returns
	-------
	True if the molecules are the same, False if not.
	"""

	# First, get the elements of the molecules.
	new_molecule_elements      = new_molecule.get_chemical_symbols()
	obtained_molecule_elements = obtained_molecule.get_chemical_symbols()

	# Second, get the positions of the molecules.
	new_molecule_positions      = new_molecule.get_positions().tolist()
	obtained_molecule_positions = obtained_molecule.get_positions().tolist()

	# Third, check that the two molecules do not overlap in elements and positions.
	for new_element, new_position in zip(new_molecule_elements, new_molecule_positions):

		# 3.1: For each atom index in the obtained molecule.
		for obtained_atom_index, (obtained_element, obtained_position) in enumerate(zip(obtained_molecule_elements, obtained_molecule_positions)):

			# 3.2: Are the two elements the same (or set to True if consider_elements is False)
			are_elements_the_same = (obtained_element == new_element) if consider_elements else True

			# 3.3: If "new" atom is the same the "obtained" atom, we have found the atom in the "obtained"
			#      molecule, so remove it from obtained_molecule_elements and obtained_molecule_positions.
			if are_elements_the_same and (get_distance(obtained_position, new_position) < max_diff_distance):
				del obtained_molecule_elements[obtained_atom_index]
				del obtained_molecule_positions[obtained_atom_index]
				break

		else:

			# 3.4: If any one atom in obtained_molecule does not fit with the new_molecule, then the two molecules are different
			return False

	# Fourth, if you got to hear, you were able to match up every atom in new_molecule with every atom in obtained_molecule. 
	#         * Therefore, new_molecule and obtained_molecule are teh same.
	return True

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


