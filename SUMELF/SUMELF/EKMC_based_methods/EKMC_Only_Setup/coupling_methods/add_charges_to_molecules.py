"""
add_charges_to_molecules.py, Geoffrey Weal, 20/3/22

This script will assign charges to the atoms in each molecule.
"""
import os
from numpy import array

from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.coupling_methods.add_charges_to_molecules_methods_v2.import_CHG_file   import import_CHG_file
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.coupling_methods.add_charges_to_molecules_methods_v2.invariance_method import assigned_ATCs_to_molecules_invariance_method

def add_charges_to_molecules(molecules, molecule_graphs, ATC_folder_path, ECCP_Information_folder_path, crystal_name, functional_and_basis_set, include_solvents, max_disparity=None):
	"""
	This method is designed to assign charges to each atom in the molecules in your crystal based on DFT calculations.

	The charges of your molecules are given in the .chg files given in your ATC_folder_path folder

	This method is able to determine the symmetries in a crystal and determine which molecules are symmetrically the same. This means you only have to determine the charges of molecules that are not symmetric in the crystal.

	This method will apply the correct charges to the correct atoms in the symmetrically equivalent molecules.

	Parameters
	----------
	molecules : dict of int : ase.Atoms
		This is a dictionary of all the molecules in your crystal. Keys are the names of the molecules (excluding solvent indicator). 
	molecule_graphs : list of networkx.Graph
		This is the list of the undirected graph representations of each molecule.
	ATC_folder_path : str.
		This is the folder that contains all the chg files to obtain the charges of atoms in each molecule.
	ECCP_Information_folder_path : str.
		This is the path to the ECCP_Information folder. This stores the Structurally_Unique_Molecule_Information.txt file that contains useful information about which molecules are structurally the same as each other. 
	functional_and_basis_set : str.
		This is the folder of information for the given functional and basis set that you would like to use. 
	max_disparity : float
		This is the maximum disparity between two dimers to be considered invariant. 
	"""

	# First, get how the molecule indices are assigned to the ATC. Also get ATC ase objects.

	# 1.1: Gather the ATCs.
	ATC_ase_objects = {}
	unique_mol_nos = []
	ATC_molecule_names = sorted([ATC_molecule_name for ATC_molecule_name in os.listdir(ATC_folder_path+'/'+crystal_name) if (os.path.isdir(ATC_folder_path+'/'+crystal_name+'/'+ATC_molecule_name) and ATC_molecule_name.startswith('molecule_'))])
	for index in range(len(ATC_molecule_names)-1,-1,-1):
		ATC_molecule_name = int(ATC_molecule_names[index].replace('molecule_','').replace('S',''))
		for molecule_name in molecules.keys():
			if molecule_name == ATC_molecule_name:
				break
		else:
			del ATC_molecule_names[index]
	for ATC_molecule_name in ATC_molecule_names:
		ATC_molecule = import_CHG_file(ATC_folder_path+'/'+crystal_name+'/'+ATC_molecule_name+'/'+functional_and_basis_set+'/output.chg')
		unique_mol_name = int(ATC_molecule_name.replace('molecule_','').replace('S',''))
		unique_mol_no   = unique_mol_name
		unique_mol_nos.append(unique_mol_no)
		ATC_ase_objects[unique_mol_no] = ATC_molecule
	
	# 1.2: Gather the information about which molecules are similar to which unique molecule.
	symmetric_to_structurally_unique = obtain_symmetric_to_structurally_unique(ECCP_Information_folder_path, unique_mol_nos, molecules, include_solvents, crystal_name)

	# 1.3: Determine which atoms in which molcules go together between the ATC_ase_objects and molecules. 
	molecule_to_atc = assigned_ATCs_to_molecules_invariance_method(ATC_ase_objects, molecules, molecule_graphs, symmetric_to_structurally_unique)

	# Second, apply the charges from the ATC to the appropriate molecule.
	for ATC_index, mol_index, ATC_idx in molecule_to_atc:
		ATC_charges = ATC_ase_objects[ATC_index].get_initial_charges()
		non_hydrogen_ATC_reordered_charges = [ATC_charges[index] for index in ATC_idx]
		molecules[mol_index].set_initial_charges(non_hydrogen_ATC_reordered_charges)

# ------------------------------------------------------------------------------------------------------------------------

def obtain_symmetric_to_structurally_unique(ECCP_Information_folder_path, unique_mol_nos, molecules, include_solvents, crystal_name):
	"""

	"""

	# First, obtain the path to the structurally unique molecule information that describes which molecules are structurally unique to eachother in the crystal unit cell.
	path_to_Structurally_Unique_Molecule_Information_file = ECCP_Information_folder_path+'/'+'Structurally_Unique_Molecule_Information.txt' # ECCP_Information_folder_path+'/'+crystal_name+'/'+'Structurally_Unique_Molecule_Information.txt'

	# Second, if the Structurally_Unique_Molecule_Information.txt file does not exist, return None to indicate this information was not given. 
	if not (os.path.exists(path_to_Structurally_Unique_Molecule_Information_file) and os.path.isfile(path_to_Structurally_Unique_Molecule_Information_file)):
		return None

	# Third, write the structurally unique molecules to the symmetric_to_structurally_unique dictionary, which are structurally unique to themselves.
	symmetric_to_structurally_unique = {unique_mol_name: unique_mol_name for unique_mol_name in unique_mol_nos}

	# Fourth, obtain the information about what molecules are structurally unique from the Structurally_Unique_Molecule_Information.txt file
	with open(path_to_Structurally_Unique_Molecule_Information_file,'r') as Structurally_Unique_Molecule_InformationFILE:
		Structurally_Unique_Molecule_InformationFILE.readline()
		for line in Structurally_Unique_Molecule_InformationFILE:

			# 4.1: Get the structurally equivalent and corresponding unique molecle from Structurally_Unique_Molecule_Information.txt line.
			line              = line.rstrip().split()
			equivalent_mol_no = int(line[0])
			unique_mol_no     = int(line[2])

			# 4.2: If you do not want to include solvents, exclude molecules not in the molecules dictionary as these will be solvents. 
			if (not include_solvents) and ((equivalent_mol_no not in molecules.keys()) or (unique_mol_no not in molecules.keys())):
				continue

			# 4.3: Check that the same equivalent_mol_no have been entered into Structurally_Unique_Molecule_Information.txt multiple times, this would indicate an error in Structurally_Unique_Molecule_Information.txt
			if equivalent_mol_no in symmetric_to_structurally_unique:
				raise Exception('Error: '+str(equivalent_mol_no)+' is already in symmetric_to_structurally_unique.\nCurrent symmetric_to_structurally_unique: '+str(symmetric_to_structurally_unique)+'\nNew input: '+str(equivalent_mol_no)+': '+str(unique_mol_no)+'\nCheck your Structurally_Unique_Molecule_Information.txt file\nfilepath: '+str(path_to_Structurally_Unique_Molecule_Information_file))

			# 4.4: Add equivalent_mol_no and unique_mol_no to symmetric_to_structurally_unique
			symmetric_to_structurally_unique[equivalent_mol_no] = unique_mol_no

	# Fifth, if you include solvents, check to see that all molecules are entered in symmetric_to_structurally_unique as keys. 
	if include_solvents and (not sorted(symmetric_to_structurally_unique.keys()) == list(range(1,len(symmetric_to_structurally_unique)+1))):
		print(sorted(symmetric_to_structurally_unique.keys()))
		print(list(range(len(symmetric_to_structurally_unique))))
		raise Exception('Error: your molecule names for this molecule may not be consecutively named and start from 1. Check your '+str(crystal_name)+' ECCP molecule files.')

	# Sixth, return symmetric_to_structurally_unique
	return symmetric_to_structurally_unique

# ------------------------------------------------------------------------------------------------------------------------

def read_chg_file(CHG_filename):
	"""
	This method will read in the chg files 

	Parameters
	----------
	CHG_filename : str
		This is the name of the chg file to be read

	Returns
	-------
	atc_data : tuple of tuples
		This is an object that contains the charges of all atoms in your molecule, given as (symbol, position, charge) for each atom.
	"""

	# First, initialise the ATC data list.
	atc_data = []

	# Second, read the ATC charge data from the .chg file.
	with open(CHG_filename,'r') as CHGfile:
		for line in CHGfile:
			symbol, xx, yy, zz, atc = line.rstrip().split()
			xx = float(xx); yy = float(yy); zz = float(zz)
			position = array([xx, yy, zz])
			atc = float(atc)
			atc /= (2.0 ** 0.5) # This is to manually correct an issue. See the Multiwfn 3.8 dev manual, Section 4.A.9, page 960.
			atc_data.append((symbol, position, atc))

	# Third, return the ATC data list
	return tuple(atc_data)

# ------------------------------------------------------------------------------------------------------------------------
