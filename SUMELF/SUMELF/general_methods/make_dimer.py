"""
make_dimer.py, Geoffrey Weal, 19/7/22

This script is designed to hold methods useful for making dimers.
"""
import numpy as np
from SUMELF.SUMELF.general_methods.unit_cell_methods import centre_molecule_in_cell

centre_of_origin = np.array((0.,0.,0.))
def make_dimer(molecules, mol1_name, mol2_name, displacement, move_centre_of_mass_by=centre_of_origin):
	"""
	Make the dimer from two molecules in the molecules list, where molecule 2 has been moved by displacement.

	Parameters
	----------
	molecules : dict. of ase.Atoms
		These are all the molecules in the crystal to make dimers of.
	mol1_name : int
		This is the name of the first molecule in the dimer from the molecules dictionary.
	mol2_name : int
		This is the name of the second molecule in the dimer from the molecules dictionary.
	displacement : numpy.array
		This is the displacement to displace the second molecule by to create the dimer.
	move_centre_of_mass_byt : numpy.array
		This is the displacement required to centre the dimer as much as possible in the origin unit cell, mainly to make viewing the molecule easier.

	Returns
	-------
	dimer : ase.Atoms
		This is the dimer
	molecule1 : ase.Atoms
		This is first molecule in the dimer.
	molecule2 : ase.Atoms
		This is second molecule in the dimer.
	"""

	# 
	if isinstance(move_centre_of_mass_by,str) and (move_centre_of_mass_by == 'to_centre_of_unit_cell'):
		move_to_centre_of_unit_cell = True
		move_centre_of_mass_by = centre_of_origin
	else:
		move_to_centre_of_unit_cell = False

	# First, get both molecules in the dimer
	molecule1 = molecules[mol1_name].copy()
	molecule2 = molecules[mol2_name].copy()

	# Second, update the positions of the molecule so they reflect being in the dimer. 
	molecule1.set_positions(molecule1.get_positions() + move_centre_of_mass_by)
	molecule2.set_positions(molecule2.get_positions() + move_centre_of_mass_by + displacement)

	# Third, get the dimer.
	dimer = (molecule1 + molecule2).copy()

	if move_to_centre_of_unit_cell:
		print('Write here')
		import pdb; pdb.set_trace()
		raise Exception('To write this method here')

	# Fourth, return the dimer and the two molecules used in the dimer. 
	return dimer, molecule1, molecule2