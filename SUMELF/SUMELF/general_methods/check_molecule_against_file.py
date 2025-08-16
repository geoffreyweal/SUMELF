"""
check_molecule_against_file.py , Geoffrey Weal, 28/3/24

This method is designed to compare the molecule object with the molecule on file to check that they are the same and if not, why. 
"""
import os
from ase.io        import read
from ase.visualize import view

# The following variable privide the thresholds for position and charge
position_max_threshold = 0.0000001
charge_max_threshold   = 0.0000001

def check_molecule_against_file(molecule, path_to_molecule_on_file):
	"""
	This method is designed to compare the molecule object with the molecule on file to check that they are the same and if not, why. 

	Parameters
	----------
	molecule : ase.Atoms
		This is the ase.Atoms object. 
	path_to_molecule_on_file : str.
		This is the path of the object save on file. This would have been created in a previous run of the ECCP. 
	"""

	# First, if the molecule does not exist on file, dont worry about running this method.
	if not os.path.exists(path_to_molecule_on_file):
		return

	# Second, if path_to_molecule_on_file ends with '.gjf' (Gaussian) or '.inp' (ORCA), check total charge rather than individual charges
	if path_to_molecule_on_file.endswith('.gjf') or path_to_molecule_on_file.endswith('.inp'):
		check_atom_charge = False
	else:
		check_atom_charge = True

	# Third, read the molecule on file.
	molecule_on_file = read(path_to_molecule_on_file)

	# Fourth, make sure that the molecule and the molecule on file have the same number of atoms.
	if not len(molecule) == len(molecule_on_file):
		to_string  = 'Error: The ase.Atoms object has a different number of atoms compared to that on file:\n'
		to_string += f'Number of atoms in ase.Atoms object: {len(molecule)}\n'
		to_string += f'Number of atoms of system on file: {len(molecule_on_file)}\n'
		to_string += f'Path to system on file: {path_to_molecule_on_file}\n'
		to_string += 'Check this'
		raise Exception(to_string)

	# Fifth, set up the list to record what atoms have issues.
	atoms_with_issues = []

	# Sixth, compare the molecule object with that on file. 
	for atom_index, (atom, atom_on_file) in enumerate(zip(molecule, molecule_on_file)):

		# 6.1: Check if the atoms are the same element:
		if atom.symbol != atom_on_file.symbol:
			atoms_with_issues.append((atom_index, atom, atom_on_file))

		# 6.2: Make sure that atoms have the same x position (within threshold)
		if abs(atom.x - atom_on_file.x) > position_max_threshold:
			atoms_with_issues.append((atom_index, atom, atom_on_file))

		# 6.3: Make sure that atoms have the same y position (within threshold)
		if abs(atom.y - atom_on_file.y) > position_max_threshold:
			atoms_with_issues.append((atom_index, atom, atom_on_file))

		# 6.4: Make sure that atoms have the same z position (within threshold)
		if abs(atom.z - atom_on_file.z) > position_max_threshold:
			atoms_with_issues.append((atom_index, atom, atom_on_file))

		# 6.5: Make sure that atoms have the charge position (within threshold)
		if check_atom_charge:
			if abs(atom.charge - atom_on_file.charge) > charge_max_threshold:
				atoms_with_issues.append((atom_index, atom, atom_on_file))

	# Seventh, raise an issue if there is a difference between the molecule and the molecule on file.
	if len(atoms_with_issues) > 0:
		to_string  = 'Error: There are one or more differences between the ase.Atoms object and that on file.\n'
		to_string += '\n'
		to_string += f'Path to system on file: {path_to_molecule_on_file}\n'
		to_string += '\n'
		to_string += 'The issues are:\n'
		to_string += '\n'
		to_string += 'atom index: atom in object --> equivalent atom in file\n'
		for atom_index, atom, atom_on_file in atoms_with_issues:
			to_string += f'{atom_index}: {str(atom)} --> {str(atom_on_file)}\n'
		to_string += '\n'
		to_string += 'Check this.'
		view([molecule, molecule_on_file])
		raise Exception(to_string)

	# Eighth, obtain the total charge for the original molecule and the molecule on file.
	total_charge_of_original_molecule = sum(molecule.get_initial_charges())
	total_charge_of_molecule_on_file  = sum(molecule_on_file.get_initial_charges())

	# Ninth, check the total charge of the two systems are the same. 
	if abs(total_charge_of_original_molecule - total_charge_of_molecule_on_file) > charge_max_threshold:
		to_string  = 'Error: The total charge between the two systems are different\n'
		to_string += '\n'
		to_string += f'Path to system on file: {path_to_molecule_on_file}\n'
		to_string += '\n'
		to_string += f'Total charge of original molecule: {total_charge_of_original_molecule}\n'
		to_string += f'Total charge of the molecule on file: {total_charge_of_molecule_on_file}\n'
		to_string += '\n'
		to_string += 'Check this.'
		view([molecule, molecule_on_file])
		raise Exception(to_string)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

