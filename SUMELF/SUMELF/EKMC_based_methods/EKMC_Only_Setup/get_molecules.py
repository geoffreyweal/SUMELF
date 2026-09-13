"""
get_molecules.py, Geoffrey Weal, 17/2/22

This script is able to obtain the molecules from the crystal
"""
import os
from copy import deepcopy
from ase.io import read
from SUMELF import less_than_or_equal_to_max_bondlength

def get_molecules(molecules_path, include_solvents=True):
	"""
	This method will obtain the molecules from the crystal

	Parameters
	----------
	molecules_path : str.
		This is the path to the folders where the molecules are stored. 
	include_solvents : bool.
		This boolean indicates if you want to include solvents in your KMC model. Default: True.

	Returns
	-------
	molecules : list of ase.Atoms
		This is a list of the molecules in the crystal.
	"""

	# First, these are the folder names.
	all_molecule_folder = 'All_Molecules'

	# Second, obtain the list of all the molecules in the "All Molecules" folder. 
	all_molecules_list = []
	for file in os.listdir(molecules_path+'/'+all_molecule_folder):
		if (os.path.isfile(molecules_path+'/'+all_molecule_folder+'/'+file) and file.startswith('molecule_') and file.endswith('.xyz')):
			molecule_name = str(file.replace('molecule_','').replace('.xyz',''))
			all_molecules_list.append((molecule_name, file))

	# Third, sort the molecules by their name, without solvent tag
	all_molecules_list.sort(key=lambda x: int(x[0].replace('S','')))

	# Fourth, this is a check to see if their are no molecules that are obviously missing.
	# Molecules may still be missing however, so check to make sure you have enough molecules
	if not ([int(x[0].replace('S','')) for x in all_molecules_list] == list(range(1,len(all_molecules_list)+1))):
		print('Error in EKMC: You may be missing some molecules in the '+str(all_molecule_folder)+' folder.')
		print('The list of molecules in the '+str(all_molecule_folder)+' should be consecutive. This may be a sign you are missing a molecule')
		print('Molecules in '+all_molecule_folder+' folder: '+str(all_molecules_list))
		print('Check this out.')
		exit('This program will finish without beginning.')

	# Fifth, get the paths to molecules to use
	molecule_names = []
	molecule_filepaths = []
	for molecule_name, molecule_filename in all_molecules_list:
		# 4.1: if you dont want to include solvents, continue on if you find a solvent. 
		if (not include_solvents) and ('S' in molecule_name):
			continue
		molecule_names.append(molecule_name)
		molecule_filepaths.append(molecules_path+'/'+all_molecule_folder+'/'+molecule_filename)
	
	# Sixth, get molecules
	molecules = [read(molecule_filepath) for molecule_filepath in molecule_filepaths]

	# Seventh, return molecules and their names. 
	return molecule_names, {int(molecule_name.replace('S','')): molecule for molecule_name, molecule in zip(molecule_names, molecules)}

	# ------------------------------------------------------------------------------------------------------


