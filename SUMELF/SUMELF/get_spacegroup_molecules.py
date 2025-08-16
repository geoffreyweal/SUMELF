"""
get_spacegroup_molecules.py, Geoffrey Weal, 17/2/22

This script will create teh molecules in the crystals that exist due to symmetry operations of the crystal's spacegroup. 
"""

import numpy as np
import networkx as nx

from ase import Atom, Atoms
from ase.cell import Cell
from ase.geometry import cellpar_to_cell

from SUMELF.SUMELF.make_crystal_methods.perform_symmetry_operation_upon_molecule import perform_symmetry_operation_upon_molecule
from SUMELF.SUMELF.make_crystal_methods.determine_is_molecule_already_recorded   import determine_is_molecule_already_recorded
from SUMELF.SUMELF.general_methods.distance_methods                              import get_distance
from SUMELF.SUMELF.general_methods.unit_cell_methods                             import get_cell_corner_points

# Set up the numpy aways for the Identity symmetry operation
Identity_Rotation    = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
Identity_Translation = np.array([0., 0., 0.])

def get_spacegroup_molecules(molecules, molecule_graphs, SolventsList, symmetry_operations, cell):
	"""
	This method is designed to obtain all the molecules that exist in a crystal due to the symmetry operations of the crystal's spacegroup. 

	This method has been designed because some methods give only part of the crystal, where the other molecules are given by the symmetry operations in symmetry_operations.

	Parameters
	----------
	molecules : dict. of ase.Atoms
		These are the molecules that make up the crystal (before any spacegroup symmetry is considered).
	molecule_graphs : dict. of networkx.Graph
		These are the graph that are associated with the molecules in the molecules list. 
	SolventsList : list of ints
		This list indicates which molecules in the molecules list are solvents.
	symmetry_operations : list of numpy.arrays
		These are the symmetry operations that the moelcules obey in the crystal.
	cell : list of floats
		This is the cell lengths and angles of the unit cell.
	
	Returns
	-------
	molecules : list of ase.Atoms
		These are the all molecules that make up the crystal, including those molecules that exist due to symmetry operations of the crystal's spacegroup. 
	molecule_graphs : list of networkx.Graph
		These are the graph that are associated with the molecules in the molecules list. 
	crystal : ase.Atoms
		This is the crystal object with all arrays and info filled in.
	SolventsList : list of ints
		This list indicates which molecules in the molecules list are solvents.
	"""

	# First, prepare the lists and dictionaries for recording data to.
	new_molecules                    = {}
	new_molecule_graphs              = {}
	MoleculeList                     = []
	new_SolventsList                 = []
	SameMoleculeSpacegroupSymOpsList = {}

	# Second, make the crystal cell matrix from the cell parameters.
	if isinstance(cell, Cell):
		crystal_cell_lattice = cell[:]
	elif len(cell) == 6:
		crystal_cell_lattice = cellpar_to_cell(cell)
	elif len(cell) == 9:
		crystal_cell_lattice = cell
	else:
		raise Exception('Error: The cell provided must be either a ase.cell.Cell object, a cellpar with a length of 6, or the cell matrix with a length of 9. cell = '+str(cell))

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Third, perform several checks to make sure that everything is ok

	# 3.1: Check that the first symmetry operation in symmetry_operations is the identity operation
	if (not (symmetry_operations[0][0] == Identity_Rotation).all()) or (not (symmetry_operations[0][1] == Identity_Translation).all()):
		to_string  = 'Error: The first symmetry operations in symmetry_operations is not the identity.\n'
		to_string += f'First rotation matrix in symmetry_operations (below):\n{symmetry_operations[0][0]}\n'
		to_string += f'First translation matrix in symmetry_operations: {symmetry_operations[0][1]}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Fourth, perform each symmetry_operation upon the main component of the crystal to create the full crystal. 

	# 4.1: Set up a counter for recording the name of the current molecule baing recorded to the crystal
	new_molecule_name = 1

	# 4.2: For each crystal symmetry operation in symmetry_operations
	for symmetry_operation_index, (scaled_rotation_matrix, scaled_translation_vector) in enumerate(symmetry_operations):

		# 4.3: For each molecule in molecules
		for original_molecule_name in sorted(molecules.keys()):

			# 4.4: Check that if we are using the identify crystal symmetry operation, that the name of original_molecule_name is the same as new_molecule_name
			if symmetry_operation_index == 0:
				if not original_molecule_name == new_molecule_name:
					to_string  = 'Error: The original molecule (with the identity symmetry operation) does not have the expected name when applied to the full crystal (including all the symmetry operations).\n'
					to_string += 'This may indicate that in the original crystal there is/are molecule(s) that are duplicated on top of each other (have the same positions).\n'
					to_string += f'Original names of molecules: {sorted(molecules.keys())}\n'
					from ase.visualize import view
					sorted_molecules = sorted(molecules.items(), key=lambda x:x[0])
					view([molecule for mol_name, molecule in sorted_molecules])
					to_string += f'Problem discovered for {original_molecule_name} where it wants to be named {new_molecule_name}\n'
					to_string += 'check this'
					raise Exception(to_string)

			# 4.5: Obtain the molecule input needed.
			molecule = molecules[original_molecule_name]
			molecule_graph = molecule_graphs[original_molecule_name]

			# 4.8: Create a new molecule that is a copy of the original molecule, and set it up to be added to the crystal.
			new_molecule = molecule.copy()
			new_molecule.set_cell(crystal_cell_lattice)
			new_molecule.set_pbc(True)

			# 4.9: perform the symmetry operation upon the molecule. 
			new_molecule = perform_symmetry_operation_upon_molecule(new_molecule, (scaled_rotation_matrix, scaled_translation_vector))

			# 4.10: Translate the new_molecule so that the centre of molecule lies inside the unit cell
			add_translation = get_translation_to_move_COM_inside_unit_cell(new_molecule, crystal_cell_lattice)
			new_molecule.set_positions(new_molecule.get_positions() + add_translation)

			# 4.11: Check to see if this molecule hasn't already been obtained from a previous symmetry operation.
			if determine_is_molecule_already_recorded(new_molecule, new_molecules, crystal_cell_lattice, super_cell_reach=1, include_hydrogen=True, consider_elements=True):
				continue

			# 4.12: Add a unique tag to identify this separate molecule in the crystal
			#new_molecule.set_tags(new_molecule_name)
			MoleculeList += [new_molecule_name]*len(new_molecule)

			# 4.13: Add an original_molecule_name that indicates molecules that are the same type of molecule in the crystal based on symmetry operations.
			if new_molecule_name in SameMoleculeSpacegroupSymOpsList:
				raise Exception(f'Error: {new_molecule_name} is already in SameMoleculeSpacegroupSymOpsList: {SameMoleculeSpacegroupSymOpsList}')
			SameMoleculeSpacegroupSymOpsList[new_molecule_name] = original_molecule_name

			# 4.14: If their are solvents in the crystal, record which molecules are solvent molecules
			if len(SolventsList) > 0:
				if original_molecule_name in SolventsList:
					new_SolventsList.append(new_molecule_name)

			# 4.15: Add this new_molecule to the lists of new_molecules.
			if new_molecule_name in new_molecules:
				raise Exception(f'Error: {new_molecule_name} is already in new_molecules: {new_molecules}')
			new_molecules[new_molecule_name] = new_molecule

			# 4.16: Add the graph for this molecule to the new_molecule_graphs list.
			if new_molecule_name in new_molecule_graphs:
				raise Exception(f'Error: {new_molecule_name} is already in new_molecule_graphs: {new_molecule_graphs}')
			new_molecule_graphs[new_molecule_name] = molecule_graph

			# 4.17: Increment new_molecule_name
			new_molecule_name += 1

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Fifth, make the crystal from the new_molecules. 

	# 5.1: Set up the ase.Atoms file that will work as the crystal.
	crystal = Atoms(cell=crystal_cell_lattice, pbc=True)
	for mol_name, molecule in sorted(new_molecules.items(), key=lambda x: x[0]):
		crystal += molecule.copy()

	# 5.2: Wrap the crystal object so that all atoms lie inside the unit cell.
	crystal.wrap()

	# 5.3: If the atoms in the crystal are all the same tag, set Tags to None
	if all(crystal.get_tags() == 0):
		crystal.set_tags(None)

	# 5.4: Add the spacegroup to the object
	# crystal.info['spacegroup'] = spacegroup.symbol

	# 5.5: Add molecule lists to the object which indicate which atoms are the same atoms across molecules, and which atoms form a molecule.
	crystal.arrays['MoleculeList'] = np.array(MoleculeList)

	# 5.6: Add information which molecules are the same due to the spacegroup of the crystal
	crystal.info['SameMoleculeSpacegroupSymOpsList'] = ','.join([(str(k)+':'+str(v)) for k, v in sorted(SameMoleculeSpacegroupSymOpsList.items(), key=lambda x: x[0])])

	# 5.7: If their are solvents in the crystal, add a tag in the info that includes which molecules are solvents
	if len(new_SolventsList) > 0:
		crystal.info['SolventsList'] = ' '.join(str(index) for index in new_SolventsList)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Sixth, create the graph for the crystal.  

	# 6.1: initialise a list to place molecule graphs where the node names have been changed to that for the crystal
	new_molecule_graphs_for_crystal = []

	# 6.2: For each molecule graph in new_molecule_graphs
	index_counter = 0
	for molecule_name in sorted(new_molecules.keys()):

		# 6.2.1: Get the set of new molecules and associated graphs.
		new_molecule       = new_molecules[molecule_name]
		new_molecule_graph = new_molecule_graphs[molecule_name]

		# 6.2.2: Make sure that the number of atoms in new_molecule_graph is the same as new_molecule:
		if not len(new_molecule_graph) == len(new_molecule):
			to_string  = 'Error: The number of atoms in molecule is not the same as its respective molecule graph.\n'
			to_string += 'Number of atoms in new_molecule: '+str(len(new_molecule))+'\n'
			to_string += 'Number of nodes in new_molecule_graph: '+str(len(new_molecule_graph))+'\n'
			to_string += 'Check this'
			raise Exception(to_string)

		# 6.2.3: Check that the indices in the graph are in consecutive order from 0 to len(molecule)-1
		if not sorted(new_molecule_graph.nodes()) == list(range(len(new_molecule_graph))):
			to_string  = 'Error: The new molecule graph does not have nodes named consecutively from 0 to len(new_molecule_graph)-1.\n'
			to_string += 'nodes in new_molecule_graph: '+str(sorted(new_molecule_graph.nodes()))+'\n'
			to_string += 'Number of nodes in new_molecule_graph: '+str(len(new_molecule_graph))+'\n'
			to_string += 'Check this'
			raise Exception(to_string)

		# 6.2.4: Make the dictionary for mapping the atoms in the molecule to the atoms in the crystal. 
		mapping = {index: index+index_counter for index in range(len(new_molecule_graph))}

		# 6.2.5: Relabel the atom in the molecule to match that in the crystal
		new_molecule_graph = nx.relabel_nodes(new_molecule_graph, mapping)

		# 6.2.6: Append the updated graph to new_molecule_graphs_for_crystal
		new_molecule_graphs_for_crystal.append(new_molecule_graph)

		# 6.2.7: Update index_counter for the next molecule in the crystal
		index_counter += len(new_molecule_graph)

	# 6.3: Make the full crystal_graph from all the molecule's graphs. 
	crystal_graph = nx.compose_all(new_molecule_graphs_for_crystal)
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Seventh, return quantities
	return new_molecules, new_molecule_graphs, crystal, crystal_graph, new_SolventsList

# -------------------------------------------------------------------------------------------------------------------------

def get_translation_to_move_COM_inside_unit_cell(new_molecule, crystal_cell_lattice):
	"""
	This method will translate the compoent inside of the unit cell. 

	This will give a ase.Atoms object that will look more like CCDC Mercury displays its crystals

	Parameters
	----------
	new_molecule : ase.Atoms
		This is the new molecule of the crystal to replicate.
	crystal_cell_lattice : np.array
		This is the matrix that includes the vector molecules of the unit cell.

	Returns
	-------
	add_translation : np.array
		This is a column vector that includes the translation vector required to move the molecule inside of the unit cell while maintaining the crystal spacegroup.
	"""

	# First, get the centre position of the origin unit cell.
	origin_cell_cell_points = get_cell_corner_points(crystal_cell_lattice,super_cell_reach=1,bottom_of_range=0)
	origin_cell_centre_of_unit_cell = sum(origin_cell_cell_points)/float(len(origin_cell_cell_points))

	# Second, get all the displacements that surround the cell around the molecule.
	cell_points = get_cell_corner_points(crystal_cell_lattice,super_cell_reach=1)

	# Third, determine the position of the molecule that is closest to the centre of the unit cell while maintaining the crystal spacegroup.
	distance_from_closest_COM_to_MOC = float('inf')
	add_translation = None
	for cell_point in cell_points:

		# 3.1: Copy the new_molecule
		new_molecule_copy = new_molecule.copy()

		# 3.2: Place it in the translated position
		new_molecule_copy.set_positions(new_molecule_copy.get_positions() + cell_point)

		# 3.3: Obtain the centre of molecule (excluding hydrogens)
		if all(symbol in ['H','D'] for symbol in set(new_molecule_copy.get_chemical_symbols())):
			centre_of_molecule = [atom.position for atom in new_molecule_copy]
		else:
			centre_of_molecule = [atom.position for atom in new_molecule_copy if (not atom.symbol == 'H')]
		centre_of_molecule = sum(centre_of_molecule)/float(len(centre_of_molecule))

		# 3.4: Determine the distance from the centre of molecule to the centre of the unit cell
		current_distance_from_COM_to_MOC = get_distance(centre_of_molecule, origin_cell_centre_of_unit_cell)
		
		# 3.5: If you have found a shorter centre of molecule to centre of the unit cell distance, record that translation vector
		if current_distance_from_COM_to_MOC < distance_from_closest_COM_to_MOC:
			distance_from_closest_COM_to_MOC = current_distance_from_COM_to_MOC
			add_translation = cell_point

	# Fourth, return the translation vector that will move the molecule inside of the unit cell while maintaining the crystal spacegroup.
	return add_translation

# -------------------------------------------------------------------------------------------------------------------------



