"""
import_CHG_file.py, Geoffrey Weal, 14/3/23

This script is designed to import a chg file as a ASE object
"""

from ase import Atoms
from ase.visualize import view

def import_CHG_file(CHG_filepath):
	"""
	This method is designed to import a chg file as a ASE object

	Parameters
	----------
	CHG_filepath : str.
		This is the path to the .chg file

	Returns
	-------
	atoms : ase.Atoms
		This is the system that is represented by the .chg file, as an ase.Atoms object
	"""

	# First, initialise the lists to return the chemical symbols, positions, and charges of atoms in the system.
	symbols = []; positions = []; charges = []

	# Second, record the chemical symbols, positions, and charges of atoms from the .chg file
	with open(CHG_filepath, 'r') as CHG_FILE:

		# 2.1: For each line in CHG_FILE
		for line in CHG_FILE:

			# 2.2: obtain the data about the atom from the line
			symbol, xx, yy, zz, charge = line.rstrip().split()

			# 2.3: Convert x, y, z and charge variables from string to floats
			position = (float(xx), float(yy), float(zz))
			charge = float(charge) / (2.0 ** 0.5) # This is to manually correct an issue. See the Multiwfn 3.8 dev manual, Section 4.A.9, page 960.
			
			# 2.4: Add information about the atom to the corresponding lists. 
			symbols  .append(symbol)
			positions.append(position)
			charges  .append(charge)

	# Third, create an Atoms object that takes the symbols, positions, and charges for the system
	atoms = Atoms(symbols=symbols, positions=positions, charges=charges)

	# Fourth, return Atoms object
	return atoms







