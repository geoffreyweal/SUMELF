"""
invariance_method.py, Geoffrey Weal, 23/2/22

This script is designed to use the procrustes analysis to determine if dimers are rotationally, translationally, and reflectively invarient.
"""
import numpy as np
from tqdm import trange, tqdm

from copy import deepcopy

from itertools import product

from ase import Atoms

from scipy.linalg import orthogonal_procrustes
from scipy.spatial import procrustes

from SUMELF import obtain_graph

from SUMELF import GraphMatcher, remove_hydrogens, get_distance

def assigned_ATCs_to_molecules_invariance_method(ATC_ase_objects, molecules, molecule_graphs, symmetric_to_structurally_unique, max_disparity=None, print_progress=True):
	"""
	This method uses the procrustes analysis to determine if a ATC are rotationally, translationally, and reflectively invarient to a molecule.

	From using the procrustes analysis, this method will determine which ATCs and molecules are equivalent

	This analysis will ignore hydrogen atoms. 

	The procrustes analysis use is from scipy. See the website below for more information: 
		* https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html

	Parameters
	----------
	ATC_ase_objects : dictionary of ase.Atoms
		The ase.Atoms object of your ATC, including info about the chemical symbols, positions, and charges of each atom in the ATC system.
	molecules : list of ase.Atoms
		This is the list of molecules that are in the crystal
	molecule_graphs : list of networkx.Graph
		This is the list of the undirected graph representations of each molecule
	symmetric_to_structurally_unique : dictionary
		This dictionary contains information about what molecules are the same as those structurally unique molecules
	max_disparity : float
		This is the maximum disparity between a ATC and molecule to be considered invariant. 

	Returns
	-------
	symmetric_ATCs : list
		A list of the indices of equivalent ATCs in the ATCs list. 
	"""

	# First, set max_disparity to a default value if it is given initially as None. 
	if max_disparity is None:
		max_disparity = 0.01

	##################################################################################################
	# Convert ATC lists into objects that can be used in this method

	if print_progress:
		print('-------------')
		print(('Determining equivalent ATCs using invarience method').upper())
	ATC_graphs = {}
	for ATC_index, ATC_ase_object in ATC_ase_objects.items():
		# Note, ATC_index should be in molecule_graphs, so we will use it, as the user SHOULD NOT have modified the molecule file created by EKMC
		#unique_mol_index = symmetric_to_structurally_unique[ATC_index]
		#ATC_graph        = deepcopy(molecule_graphs[unique_mol_index])
		ATC_graphs[ATC_index] = deepcopy(molecule_graphs[ATC_index])
		'''
		ATC_graph = obtain_graph(ATC_ase_object)
		ATC_graphs[ATC_index]      = ATC_graph
		ATC_ase_objects[ATC_index] = ATC_ase_object
		'''

	# Obtain all molecules without hydrogen
	non_hydrogen_ATCs = {} # for debugging
	non_hydrogen_ATCs_elements  = {}
	non_hydrogen_ATCs_positions = {}
	non_hydrogen_ATC_graphs     = {}

	for ATC_index, ATC_ase_object in ATC_ase_objects.items():
		non_hydrogen_ATC, non_hydrogen_ATC_graph = remove_hydrogens(ATC_ase_object, graph=ATC_graphs[ATC_index])
		non_hydrogen_ATCs[ATC_index]           = non_hydrogen_ATC.copy()
		non_hydrogen_ATCs_elements[ATC_index]  = non_hydrogen_ATC.get_chemical_symbols()
		non_hydrogen_ATCs_positions[ATC_index] = non_hydrogen_ATC.get_positions()
		non_hydrogen_ATC_graphs[ATC_index]     = non_hydrogen_ATC_graph

	##################################################################################################
	# Convert molecule ase objects into objects that can be used in this method

	# Obtain all molecules without hydrogen
	non_hydrogen_molecules = {} # for debugging
	non_hydrogen_molecules_elements  = {}
	non_hydrogen_molecules_positions = {}
	non_hydrogen_molecules_graphs    = {}
	for mol_index, mol_ase_object in molecules.items():
		non_hydrogen_molecule, non_hydrogen_molecule_graph = remove_hydrogens(mol_ase_object, graph=molecule_graphs[mol_index])
		non_hydrogen_molecules[mol_index]           = non_hydrogen_molecule.copy()
		non_hydrogen_molecules_elements[mol_index]  = non_hydrogen_molecule.get_chemical_symbols()
		non_hydrogen_molecules_positions[mol_index] = non_hydrogen_molecule.get_positions()
		non_hydrogen_molecules_graphs[mol_index]    = non_hydrogen_molecule_graph

	##################################################################################################
	# Now compare the ATCs and molecules and match the ATCs and their indices to molecules and their indices.
	# Will use rotational, translational, and reflective symmetry to determine this.
	# The ATC positions should be exactly the same as their molecule counterparts, 
	# so keep comparisons very tight.

	nn = len(molecules) if (symmetric_to_structurally_unique is not None) else len(ATC_ase_objects)*len(molecules)
	if print_progress:
		toString = 'Comparing translational, rotational, and reflective invarience between ATCs. '
		if symmetric_to_structurally_unique is None:
			toString += str(len(ATC_ase_objects))+' ATCs to be examined.'
		toString += 'This can take a while with large and complex ATCs.'
		print(toString)
		pbar = tqdm(total=nn,unit='task')
	molecule_to_atc = []
	for mol_index in molecules.keys():
		# get the position data for the molecule.
		non_hydrogen_molecule_elements = non_hydrogen_molecules_elements[mol_index]
		non_hydrogen_molecule_positions = non_hydrogen_molecules_positions[mol_index]
		non_hydrogen_molecule_graph = non_hydrogen_molecules_graphs[mol_index]

		# determine which ATC indices to check, either all or just the structurally unique molecule as assigned by 
		if symmetric_to_structurally_unique is not None:
			ATC_indices_to_check = tuple([symmetric_to_structurally_unique[mol_index]])
		else:
			ATC_indices_to_check = tuple(range(len(ATC_ase_objects)))

		# Check all the combinations of ATC to mol indices to look through. 
		is_varient_without_hydrogens = False
		for ATC_index in ATC_indices_to_check:
			# If a varient is already found for this molecule, just increase the update pbar counter and move on
			if is_varient_without_hydrogens:
				if print_progress:
					pbar.set_description('Found Varient for molecule '+str(mol_index+1)+'. Do not need to compare ATC '+str(ATC_index+1)+' and molecule '+str(mol_index+1))
					pbar.update(1)
				continue

			# get the position data for the ATC.
			non_hydrogen_ATC_elements  = non_hydrogen_ATCs_elements[ATC_index]
			non_hydrogen_ATC_positions = non_hydrogen_ATCs_positions[ATC_index]
			non_hydrogen_ATC_graph     = non_hydrogen_ATC_graphs[ATC_index]

			# write description.
			if print_progress:
				pbar.set_description('Comparing ATC '+str(ATC_index)+' and molecule '+str(mol_index))

			# --------------------------------------------------------
			# THis is weird, to fix latter

			print('Warning. Need to fix SUMELF/SUMELF/EKMC_based_methods/EKMC_Only_Setup/coupling_methods/add_charges_to_molecules_methods_v2/invariance_method.py tags problem.')

			for index in non_hydrogen_molecule_graph.nodes:
				if 'tags' in non_hydrogen_molecule_graph.nodes[index]:
					del non_hydrogen_molecule_graph.nodes[index]['tags']

			for index in non_hydrogen_ATC_graph.nodes:
				if 'tags' in non_hydrogen_ATC_graph.nodes[index]:
					del non_hydrogen_ATC_graph.nodes[index]['tags']

			# --------------------------------------------------------

			# Get all the indices that are equivalent to eachother in each molecule in each dimer.
			GM = GraphMatcher(non_hydrogen_molecule_graph, non_hydrogen_ATC_graph)
			emATC_indices_molecule = sorted(GM.get_all_unique_isomorphic_graphs(), key=lambda x: tuple([value for (key, value) in sorted(x.items())]))

			# In this for loop, go though all the realistic ways that the atoms in the molecule 
			# can be permuted to give the ATC.
			for non_hydrogen_ATC_to_mol_comparison in emATC_indices_molecule:
				non_hydrogen_ATC_idx = get_permutated_indices_list(non_hydrogen_ATC_to_mol_comparison)
				non_hydrogen_ATC_reordered_elements = [non_hydrogen_ATC_elements[index] for index in non_hydrogen_ATC_idx]
				non_hydrogen_ATC_reordered_positions = deepcopy(non_hydrogen_ATC_positions)[non_hydrogen_ATC_idx, :]

				# Create variable for if a mol to ATC comparison is found. 
				# This will allow us to break out of this for loop if such a comparison is found. 
				found_mol_to_ATC_comparison = False
				
				# Check if your molecule and ATC image are rotationally invarient with and without reflection applied to the ATC image. 
				for reflection_matrix in [np.array([1., 1., 1.]), np.array([1., 1., -1.])]:

					# Check if the molecule and ATC image are rotationally invarient.
					reflected_non_hydrogen_ATC_reordered_positions = deepcopy(non_hydrogen_ATC_reordered_positions) * reflection_matrix
					is_varient_without_hydrogens, reflected_R_matrix, reflected_mean_mol, reflected_mean_ATC = are_molecule_and_ATC_variant(non_hydrogen_molecule_elements, non_hydrogen_ATC_reordered_elements, non_hydrogen_molecule_positions, reflected_non_hydrogen_ATC_reordered_positions, max_disparity)
					
					if is_varient_without_hydrogens:

						# Determine if the hydrogen atoms allign with each other also.
						reflected_ATC_ase_object = ATC_ase_objects[ATC_index].copy()
						reflected_ATC_ase_object.set_positions(reflected_ATC_ase_object.get_positions() * reflection_matrix)
						associated_molecule = molecules[mol_index].copy()

						# Determine if the hydrogen atoms allign with each other also.
						ATC_to_mol_comparison = assign_hydrogens_to_ATC_to_mol_comparison(reflected_ATC_ase_object, associated_molecule, non_hydrogen_ATC_to_mol_comparison, reflected_R_matrix, reflected_mean_mol, reflected_mean_ATC)
						if hydrogen_atoms_allign(associated_molecule, reflected_ATC_ase_object, ATC_to_mol_comparison, reflected_R_matrix, reflected_mean_mol, reflected_mean_ATC, max_disparity):
							ATC_idx = get_permutated_indices_list(ATC_to_mol_comparison)
							mol_idx = reverse_idx(ATC_idx)
							molecule_to_atc.append((ATC_index, mol_index, mol_idx))
							found_mol_to_ATC_comparison = True
							break

				# If a comparison between mol to ATC is found, we can break out of this for loop. 
				if found_mol_to_ATC_comparison:
					break

			# If you contains symmetric_to_structurally_unique, This will determine if there was a problem in getting a varient from your expected structurally unique molecule.
			if (symmetric_to_structurally_unique is not None) and (not is_varient_without_hydrogens):
				to_string  = 'Error: For some reason, could not find a relationship between molecule '+str(mol_index)+' and the ATC_molecule '+str(ATC_index)+'\n'
				to_string += 'Check the molecules files in your ECCP_Information folder and the molecules in your Unique_ATC_Calc_Jobs/All_ATC_Calc_Jobs folder, as this may indicates these molecule files are different. They should probably contain the same atoms and elements.\n'
				to_string += '. Check this out'
				print(to_string)
				raise Exception(to_string)

			if print_progress:
				pbar.update(1)

	if print_progress:
		pbar.close()

	# Check to make sure that each molecule has been assigned to a ATC.
	molecules_that_have_been_assigned_an_ATC_to = [mol_index for ATC_index, mol_index, ATC_idx in molecule_to_atc]
	if not sorted(molecules_that_have_been_assigned_an_ATC_to) == sorted(molecules.keys()):
		print('Error in def assigned_ATCs_to_molecules_invariance_method, in invariance_method.py')
		print('Not all molecules have been assign an ATC')
		unassigned_molecules = list(set(molecules.keys()) | set(molecules_that_have_been_assigned_an_ATC_to))
		print('Molecules that have not been assigned: '+str(unassigned_molecules))
		print('Molecules that have been assigned: '+str(molecules_that_have_been_assigned_an_ATC_to))
		exit('This program with finished without completing.')

	# hopefully all ATC's have been matched with the appropriate molecule.
	return molecule_to_atc

##############################################################################################################

def get_atom_data_from_ATC(atomic_transition_charge):
	symbols = []
	positions = []
	charges = []
	for symbol, position, atc in atomic_transition_charge:
		symbols.append(symbol)
		positions.append(position)
		charges.append(atc)
	positions = np.array(positions)
	return symbols, positions, charges

##############################################################################################################

def get_permutated_indices_list(comparison): 
	"""
	This method will provide the permutation list of how to reorder the atoms indices in a molecule from a dictionary that tell this program how the indices of one molecule translate into another equivalent molecule. 

	Parameters
	----------
	comparison : dict.
		A dictionary that contains how the indices of one molecule translate into another identical molecule. 

	Returns
	-------
	idx : list
		The permutation list that tells the program how to reorder atoms in a list.

	"""
	if not (sorted(comparison.keys()) == list(range(len(comparison)))):
		raise Exception('huh?')
	if not (sorted(comparison.values()) == list(range(len(comparison)))):
		raise Exception('huh?')
	permutation = [mol_index for ATC_index, mol_index in sorted(comparison.items())]
	return np.array(permutation)

##############################################################################################################

def are_molecule_and_ATC_variant(mol_elements, ATC_elements, mol_distances, ATC_distances, max_distance_disparity): 
	"""
	This method will check that the elements in the ATC and the molecule are translationally, rotationally, and reflectively varient
	
	Parameters
	----------
	ATC_elements : list
		A list of the elements of atoms in the ATC
	mol_elements : list
		A list of the elements of atoms in the molecule
	ATC_distances : np.array
		A numpy array of the positions of atoms in the ATC
	mol_distances : np.array
		A numpy array of the positions of atoms in the molecule
	max_distance_disparity: float
		This is the maximum that any two "could be equivalent" atoms can be between dimer 1 and dimer 2 for dimers 1 and 2 to be considered variant

	Returns
	-------
	idx : list
		The permutation list that tells the program how to reorder atoms in a list.

	"""
	#If the order of atoms are not the same, something actually may have gone wrong, so check this out
	if not (mol_elements == ATC_elements):
		print('Error in invarient_method.py')
		print('The list of elements for mol_elements and ATC_elements are not the same. ')
		print('However, they should be the same as the networkx graph nodes have been given information about the element for each atom.')
		print('Therefore, this should have been picked up by GraphMatcher object.')
		print('Check this out')
		exit('This program will finish without completing.')

	# determine if the dimer are varient given the particular ordering of atoms in dimers 1 and 2.
	positions_are_variant, R_matrix, mean_mol, mean_ATC = determine_if_positions_are_variant(mol_distances, ATC_distances, max_distance_disparity=max_distance_disparity)

	return positions_are_variant, R_matrix, mean_mol, mean_ATC

def determine_if_positions_are_variant(data1, data2, max_distance_disparity=0.01):
	"""
	This method will determine if two dimers/molecules are translationally, rotationally, and reflectively variant.

	This method have been modified from the _procrustes.py from scipy, see below: https://github.com/scipy/scipy/blob/v1.8.0/scipy/spatial/_procrustes.py#L15-L130

	Parameters
	----------
	data1 : numpy.array
		A array of positions from dimer1/molecule1
	data1 : numpy.array
		A array of positions from dimer2/molecule2
	max_distance_disparity: float
		This is the maximum that any two "could be equivalent" atoms can be between dimer 1 and dimer 2 for dimers 1 and 2 to be considered variant

	Returns
	-------
	True if the two dimers are translationally, rotationally, and reflectively variant. False if not. 

	"""
	mtx1 = np.array(data1, dtype=np.double, copy=True)
	mtx2 = np.array(data2, dtype=np.double, copy=True)

	if mtx1.ndim != 2 or mtx2.ndim != 2:
		raise ValueError("Input matrices must be two-dimensional")
	if mtx1.shape != mtx2.shape:
		raise ValueError("Input matrices must be of same shape")
	if mtx1.size == 0:
		raise ValueError("Input matrices must be >0 rows and >0 cols")

	# translate all the data to the origin
	mean_mtx1 = np.mean(mtx1, 0)
	mean_mtx2 = np.mean(mtx2, 0)
	mtx1 -= mean_mtx1
	mtx2 -= mean_mtx2

	'''
	# We dont require scaling variance in our analysis, so ignored this part of the scipy code. 
	norm1 = np.linalg.norm(mtx1)
	norm2 = np.linalg.norm(mtx2)

	if norm1 == 0 or norm2 == 0:
	    raise ValueError("Input matrices must contain >1 unique points")

	# change scaling of data (in rows) such that trace(mtx*mtx') = 1
	mtx1 /= norm1
	mtx2 /= norm2
	'''

	# transform mtx2 to minimize disparity. This method also takes into account reflective variance
	#R_matrix, ss = orthogonal_procrustes(mtx1, mtx2)
	R_matrix, ss = procrustes(mtx1, mtx2, scaling=False, reflection=False) # Use this method because we want to exclude reflections

	# Obtain the rotated version of data2
	mtx2_rotated = np.dot(mtx2, R_matrix.T) ##* ss
	#mtx2_rotated = np.matmul(R_matrix,mtx2.T).T 

	# get the difference between atom positions between each dimer
	difference_in_atom_positions = np.linalg.norm(mtx1 - mtx2_rotated, axis=1)

	# True if the two dimers are translationally, rotationally, and reflectively variant. False if not. 
	are_varient = all([(distance < max_distance_disparity) for distance in difference_in_atom_positions])

	return are_varient, R_matrix, mean_mtx1, mean_mtx2

##############################################################################################################

def assign_hydrogens_to_ATC_to_mol_comparison(ATC_molecule, molecule, non_hydrogen_ATC_to_mol_comparison, R_matrix, mean_mol, mean_ATC):
	"""
	This method will match the charge in the ATC file to the right atom in the molecule. 

	Parameters
	----------
	ATC_molecule : ase.Atoms
		This is the ase object that represents the ATC. This inlcudes the charges
	molecule : ase.Atoms
		This is the molecule that you want to assign charges from the ATC to.
	R_matrix : numpy.array
		This is the rotation matrix which will rotate the ATC_molecule onto the molecule. Relectivity has been included in this matrix if needed.
	mean_mol : numpy.array
		This is the translational vector to move the molecule to the origin
	mean_ATC : numpy.array
		This is the translational vector to move the ATC object to the origin

	Returns
	-------
	molecule_charges_from_ATC : list
		This is a list of the charges from the ATC to assign to the molecule.
	"""

	# First, get the elements from the molecule and ATC ase objects.
	ATC_symbols = ATC_molecule.get_chemical_symbols()
	mol_symbols = molecule.get_chemical_symbols()

	# Second, create dictionaries to convert indices in non-hydrogen and hydrogen systems.
	def get_non_hydrogen_to_hydrogen_conversion_dict(ATC_symbols):
		non_hydrogen_to_hydrogen_index_comparisons_ATC_molecule = {}
		non_hydrogen_index = 0
		for index, ATC_symbol in zip(range(len(ATC_symbols)), ATC_symbols):
			if not ATC_symbols[index] == 'H':
				non_hydrogen_to_hydrogen_index_comparisons_ATC_molecule[non_hydrogen_index] = index
				non_hydrogen_index += 1
		return non_hydrogen_to_hydrogen_index_comparisons_ATC_molecule
	non_hydrogen_to_hydrogen_index_comparisons_ATC_molecule = get_non_hydrogen_to_hydrogen_conversion_dict(ATC_symbols)
	non_hydrogen_to_hydrogen_index_comparisons_molecule     = get_non_hydrogen_to_hydrogen_conversion_dict(mol_symbols)

	# Third, create a comparison between ATC_molecule to molecule in a hydrogen-inclusive system.
	non_hydrogen_ATC_to_mol_comparison = {non_hydrogen_to_hydrogen_index_comparisons_ATC_molecule[non_hydrogen_ATC_index]: non_hydrogen_to_hydrogen_index_comparisons_molecule[non_hydrogen_mol_index] for non_hydrogen_mol_index, non_hydrogen_ATC_index in non_hydrogen_ATC_to_mol_comparison.items()}

	# Fourth, get the positons from the molecule and ATC ase objects.
	mol_positions = molecule.get_positions()
	ATC_positions = ATC_molecule.get_positions()

	# Fifth, move the ATC object on top of the molecule object.
	ATC_positions_moved_on_top_of_mol = np.dot(ATC_positions - mean_ATC, R_matrix.T) + mean_mol

	# Sixth, assign indices of ATC_molecules to molecules after rotation of ATC_molecules.
	ATC_to_mol_comparison_spatial = {}
	non_hydrogen_ATC_to_mol_comparison_spatial = {}
	for mol_index, mol_symbol, mol_position in zip(range(len(molecule)), mol_symbols, mol_positions):

		# 6.1: Determine which atoms the equivalent spatially using a spatial comparison.
		shortest_ATC_index = -1
		shortest_distance  = float('inf')
		for ATC_index, ATC_symbol, ATC_position in zip(range(len(ATC_molecule)), ATC_symbols, ATC_positions_moved_on_top_of_mol):
			if not (ATC_symbol == mol_symbol):
				continue
			distance = get_distance(mol_position, ATC_position)
			if distance < shortest_distance:
				shortest_ATC_index = ATC_index
				shortest_distance  = distance

		# 6.2: Check that the shortest distance is very short, as these two atoms should be on top of each other.
		if shortest_distance > 0.01:
			shortest_ATC_symbol = ATC_symbols[shortest_ATC_index]
			shortest_ATC_position = ATC_positions_moved_on_top_of_mol[shortest_ATC_index]
			molecule_modified = molecule.copy()
			molecule_modified.set_cell(None)
			ATC_molecule_modified = ATC_molecule.copy()
			ATC_molecule_modified.set_positions(ATC_positions_moved_on_top_of_mol)
			print('Huh?, shortest_distance: '+str(shortest_distance))
			print('Showing molecule and ATC_molecule. Please wait for these to load in a viewer.')
			print('Problems:')
			print('* molecule     --> index: '+str(mol_index)+', symbol: '+str(mol_symbol)+', position: '+str(mol_position))
			print('* ATC_molecule --> index: '+str(shortest_ATC_index)+', symbol: '+str(shortest_ATC_symbol)+', position: '+str(shortest_ATC_position))
			raise Exception('Huh?, shortest_distance: '+str(shortest_distance))

		# 6.3: If shortest_ATC_index is -1, we could not match up an atom in AYC_molecule to molecule. 
		if shortest_ATC_index == -1:
			raise Exception('Huh?, shortest_ATC_index = -1')

		# 6.4: Assign index in ATC_molecule to molecule.
		ATC_to_mol_comparison_spatial[shortest_ATC_index] = mol_index

		# 6.4: Assign index in ATC_molecule to molecule for non hydrogen atoms. 
		if not (mol_symbol == 'H'):
			non_hydrogen_ATC_to_mol_comparison_spatial[shortest_ATC_index] = mol_index
				
	# Seventh, check that the non-hydrogen atoms from spatial method match what we expect from graph theory of non-hydrogen systems. 
	if not (non_hydrogen_ATC_to_mol_comparison_spatial == non_hydrogen_ATC_to_mol_comparison):
		raise Exception('huh, dictionaries')

	# Eighth, Check that the keys and values in ATC_to_mol_comparison_spatial are consecutive and start at 0.
	if not (sorted(ATC_to_mol_comparison_spatial.keys()) == list(range(len(ATC_molecule)))):
		raise Exception('huh, dictionaries')
	if not (sorted(ATC_to_mol_comparison_spatial.values()) == list(range(len(molecule)))):
		raise Exception('huh, dictionaries')

	# Ninth, return ATC_to_mol_comparison_spatial
	return ATC_to_mol_comparison_spatial

##############################################################################################################

def hydrogen_atoms_allign(associated_molecule, ATC_ase_object, ATC_to_mol_comparison, R_matrix, mean_mol, mean_ATC, max_distance_disparity):
	"""
	This method is designed to check if the hydrogens that are mapped from the ATC object onto the molecule's hydrogens allign.

	If they do, great! If not, we may not have the correct allignment yet. 

	Parameters
	----------
	associated_molecule : ase.Atoms
		This is the molecule you are wanting to map the ATC object on to
	ATC_ase_object : ase.Atoms
		This Atoms object holds the ATC data for your molecule. 
	R_matrix : 3 x 3 numpy.array
		This is the rotation matrix required to rotate ATC_ase_object onto associated_molecule
	mean_mol : 3 x 1 numpy.array
		This is the centre of position for associated_molecule
	mean_ATC : 3 x 1 numpy.array
		This is the centre of position for ATC_ase_object

	Returns
	-------
	True if the hydrogen atoms in associated_molecule and ATC_ase_object are within tolerance of each other. False if not.
	"""

	# Preamble, initialise a variable to record if the hydrogen atoms are found in the right place. 
	hydrogens_within_distance_criteria = True

	# First, look through the ATC_to_mol_comparison dictionary and compare atoms in molecule and ATC object
	for ATC_index, mol_index in ATC_to_mol_comparison.items():

		# 1.1: Check that the synbols in molecule and ATC object are the same. 
		mol_symbol = associated_molecule[mol_index].symbol
		ATC_symbol = ATC_ase_object[ATC_index].symbol
		if not mol_symbol == ATC_symbol:
			raise Exception('Mol and ATC symbols not the same.')

		# 1.2: Get the position of the current atoms of interest in molecule and ATC object.
		mol_position = associated_molecule[mol_index].position
		ATC_position = np.dot(ATC_ase_object[ATC_index].position - mean_ATC, R_matrix.T) + mean_mol

		# 1.3: Check if the position of the atoms of interest in molecule and ATC object are in the same place or not.
		if get_distance(mol_position, ATC_position) > max_distance_disparity:
			if mol_symbol in ['H', 'D']:
				hydrogens_within_distance_criteria = False
			else:
				raise Exception('Error: '+str(mol_symbol)+' in is not in the same place in molecule and ATC object.')

	return hydrogens_within_distance_criteria

##############################################################################################################

def reverse_idx(ATC_idx):
	"""
	This will reverse ATC_idx to mol_idx

	Parameters
	----------
	ATC_idx : np.array
		This is the index list to turn a ATC into a mol.

	Returns
	-------
	mol_idx : np.array
		This is the index list to turn a mol into a ATC.
	"""

	# First, make sure that all the indices that are expected are in the ATC_idx list.
	if not (sorted(ATC_idx) == list(range(len(ATC_idx)))):
		raise Exception('huh, dictionaries')

	# Second, create a list of ATC indices and associated mol indices.
	ATC_to_mol_list = [(ATC_index, mol_index) for ATC_index, mol_index in zip(ATC_idx, range(len(ATC_idx)))]
	
	# Third, sort the list by the ATC indices.
	ATC_to_mol_list.sort(key=lambda x:x[0])

	# Fourth, create a np.array that will allow mol indices to be converted into ATC indices.
	mol_idx = np.array([mol_index for ATC_index, mol_index in ATC_to_mol_list])

	# Fifth, return mol_idx
	return mol_idx

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

def procrustes(X, Y, scaling=False, reflection='best'):
    """

    https://stackoverflow.com/questions/18925181/procrustes-analysis-with-numpy

    A port of MATLAB's `procrustes` function to Numpy.

    Procrustes analysis determines a linear transformation (translation,
    reflection, orthogonal rotation and scaling) of the points in Y to best
    conform them to the points in matrix X, using the sum of squared errors
    as the goodness of fit criterion.

        d, Z, [tform] = procrustes(X, Y)

    Inputs:
    ------------
    X, Y    
        matrices of target and input coordinates. they must have equal
        numbers of  points (rows), but Y may have fewer dimensions
        (columns) than X.

    scaling 
        if False, the scaling component of the transformation is forced
        to 1

    reflection
        if 'best' (default), the transformation solution may or may not
        include a reflection component, depending on which fits the data
        best. setting reflection to True or False forces a solution with
        reflection or no reflection respectively.

    Outputs
    ------------
    d       
        the residual sum of squared errors, normalized according to a
        measure of the scale of X, ((X - X.mean(0))**2).sum()

    Z
        the matrix of transformed Y-values

    tform   
        a dict specifying the rotation, translation and scaling that
        maps X --> Y

    """

    n,m = X.shape
    ny,my = Y.shape

    muX = X.mean(0)
    muY = Y.mean(0)

    X0 = X - muX
    Y0 = Y - muY

    ssX = (X0**2.).sum()
    ssY = (Y0**2.).sum()

    # centred Frobenius norm
    normX = np.sqrt(ssX)
    normY = np.sqrt(ssY)

    # scale to equal (unit) norm
    X0 /= normX
    Y0 /= normY

    if my < m:
        Y0 = np.concatenate((Y0, np.zeros(n, m-my)),0)

    # optimum rotation matrix of Y
    A = np.dot(X0.T, Y0)
    U,s,Vt = np.linalg.svd(A,full_matrices=False)
    V = Vt.T
    T = np.dot(V, U.T)

    if reflection != 'best':

        # does the current solution use a reflection?
        have_reflection = np.linalg.det(T) < 0

        # if that's not what was specified, force another reflection
        if reflection != have_reflection:
            V[:,-1] *= -1
            s[-1] *= -1
            T = np.dot(V, U.T)

    traceTA = s.sum()

    if scaling:

        # optimum scaling of Y
        b = traceTA * normX / normY

        # standarised distance between X and b*Y*T + c
        d = 1 - traceTA**2

        # transformed coords
        Z = normX*traceTA*np.dot(Y0, T) + muX

    else:
        b = 1
        d = 1 + ssY/ssX - 2 * traceTA * normY / normX
        Z = normY*np.dot(Y0, T) + muX

    # transformation matrix
    if my < m:
        T = T[:my,:]
    c = muX - b*np.dot(muY, T)
    
    #transformation values 
    tform = {'rotation':T, 'scale':b, 'translation':c}
   
    #return d, Z, tform
    if get_distance(muX, muY) > 0.01:
    	raise Exception('issue here')

    return T, b

