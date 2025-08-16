"""
utilities.py, Geoffrey Weal, 8/3/2024

This script contains methods to help with processing crystals with the methods in the process_crystal_methods folder. 
"""
import numpy as np

def get_solvent_molecule_names_from_SolventsList(crystal, molecules):
	"""
	This method is designed to get the solvents in the crystal. 

	Parameters
	----------
	crystal : ase.Atoms
		This is your crystal you are analysing with the ECCP method.
	molecules : dict. of {int: ase.Atoms}
		This dictionary contains all the molecules in the crystal, where the keys indicates the name of the crystal (as an int).
	
	Returns
	-------
	solvent_molecule_names : list of ints
		This list contains all the indices of the solvents in your crystal 
	"""

	# First, get the SolventsList from crystal file if it has it.	
	SolventsList = crystal.info['SolventsList'] if ('SolventsList' in crystal.info.keys()) else []
	if   isinstance(SolventsList,list):
		pass
	elif type(SolventsList) == type(np.array([])):
		SolventsList = list(SolventsList)
	else:
		SolventsList = [int(SolventsList)]

	# Second, only include the molecule indices that are in molecules
	solvent_molecule_names = [mol_name for mol_name in SolventsList if (mol_name in molecules.keys())]

	# Third, return solvent_molecule_names
	return solvent_molecule_names