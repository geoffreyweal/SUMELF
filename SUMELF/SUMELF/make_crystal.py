"""
make_crystal.py, Geoffrey Weal, 29.5.22

This script id designed to construct the crystal from the main component of the crystal. 
"""
import itertools
import numpy as np
from copy import deepcopy
from collections import Counter

from ase import Atoms, Atom
from ase.cell import Cell
from ase.geometry import cellpar_to_cell

from networkx import Graph, relabel_nodes, compose, connected_components

from SUMELF.SUMELF.make_crystal_methods.perform_symmetry_operation_upon_molecule import perform_symmetry_operation_upon_molecule 
from SUMELF.SUMELF.graph_methods.obtain_graph                                    import obtain_graph 
from SUMELF.SUMELF.process_crystal_methods.make_molecule                         import make_molecule 
from SUMELF.SUMELF.general_methods.general_molecules_methods                     import get_translation_to_move_COM_inside_unit_cell
from SUMELF.SUMELF.make_crystal_methods.determine_is_molecule_already_recorded   import determine_is_molecule_already_recorded
from SUMELF.SUMELF.graph_methods.get_properties_from_graph                       import add_node_properties_to_graph
from SUMELF.SUMELF.general_methods.get_symmetry_operations                       import get_symmetry_operators_string
#from SUMELF.SUMELF.general_methods.assign_bonding_system_methods                 import assign_bonding_system

def make_crystal(molecules, symmetry_operations, cell, wrap=False, solvent_components=[], remove_solvent=False, molecule_graphs=None, return_all_molecules=False):
	"""
	This method is designed to construct the crystal from the main component of the crystal. 

	Parameters
	----------
	molecules : dict. of ase.Atoms
		This dictionary contains all the molecules to place into a crystal. 
	symmetry_operations : list of (numpy.array, numpy.array)
		These are all the symmetry operations as numpy array objects
	cell : list of floats
		This is the cell lengths and angles of the unit cell.
	wrap : bool
		If True, wrap the crystal. If False, do not. Default: False
	solvent_components : list of ints
		This is the list of molecules that have been recorded as solvents in the crystal. Default: []
	remove_solvent : bool.
		This tag indicates if you want to remove the solvents from the crystal, as marked by the solvent_components list. Default: False
	molecule_graphs : dict. of networkx.Graph or None
		These are the graphs that go with each molecule in the molecules list. If None given, the molecule_graphs will not be used in this method. Default: None.
	return_all_molecules : bool.
		If True, also return all the molecules in the crystal and their graphs.

	Return
	------
	crystal : ase.Atoms
		This is the full crystal object
	"""

	'''
	# Preliminary Step 1: Check to make sure that the molecules in the molecules dictionary are named consecutively from 1 to len(molecules)
	if not (sorted(molecules.keys()) == list(range(1,len(molecules)+1))):
		to_string  = f'Error: The molecules in the molecules dictionary are not consecutively ordered from 1 to len(molecules) (which is {len(molecules)}).\n'
		to_string += f'Molecule names in the molecules dictionary: {sorted(molecules.keys())}\n'
		to_string += f'Expected molecule names: {list(range(1,len(molecules)+1))}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Preliminary Step 2: Check to make sure that the graphs in molecule_graphs are named consecutively from 1 to len(molecule_graphs)
	if molecule_graphs is not None:
		if not (sorted(molecule_graphs.keys()) == list(range(1,len(molecule_graphs)+1))):
			to_string  = f'Error: The molecule graphs in molecule_graphs are not consecutively ordered from 1 to len(molecule_graphs) (which is {len(molecule_graphs)}).\n'
			to_string += f'Molecule names in molecule_graphs: {sorted(molecule_graphs.keys())}\n'
			to_string += f'Expected molecule names: {list(range(1,len(molecule_graphs)+1))}\n'
			to_string += 'Check this.'
			raise Exception(to_string)
	'''

	# Preliminary Step 3: Check to make sure all the molecules have an associated graph if molecule_graphs not None. 
	if molecule_graphs is not None:
		if not (sorted(molecules.keys()) == sorted(molecule_graphs.keys())):
			to_string  = f'Error: Some of the molecules and/or their associated graphs are missing\n'
			to_string += f'Molecule names in molecules: {sorted(molecules.keys())}\n'
			to_string += f'Molecule names in molecule_graphs: {sorted(molecule_graphs.keys())}\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# Preliminary Step 4: Modifications are made to the molecules list during this method, so make a copy of the molecules list.
	copied_molecules = {mol_name: molecules[mol_name].copy() for mol_name in molecules.keys()}

	# Preliminary Step 5: If you do not want to have your molecules wrapped, make sure that there are molecule_graphs
	if (not wrap) and (molecule_graphs is None):
		to_string  = 'Error: You are wanting to have attached, non wrapped molecules. However, you do not have the graphs of your molecules.\n'
		to_string += 'Molecule graphs are needed for reconstructing the molecules after being wrapped.\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# ==================================================================================================================================

	# First, copy the molecule graphs if given, and check them out.
	if molecule_graphs is not None:
		copied_molecule_graphs = deepcopy(molecule_graphs)
		check_copied_graphs(copied_molecules, copied_molecule_graphs)

	# Second, make the crystal cell matrix from the cell parameters
	if isinstance(cell, Cell):
		crystal_cell_lattice = cell[:]
	elif len(cell) == 6:
		crystal_cell_lattice = cellpar_to_cell(cell)
	elif len(cell) == 9:
		crystal_cell_lattice = cell
	else:
		raise Exception('Error: The cell provided must be either a ase.cell.Cell object, a cellpar with a length of 6, or the cell matrix with a length of 9. cell = '+str(cell))

	# Third, remove solvents if you do not want to include solvents in your crystal.
	if remove_solvent:
		remove_solvent_from_molecule_dicts(solvent_components, copied_molecules, copied_molecule_graphs)

	# ==================================================================================================================================

	# Fourth, initalise a dictinoary to hold the molecules that will be placed in the crystal object.
	all_molecules_in_crystal = {}

	# Fifth, initalise a list to store which atoms belong with which molecules.
	MoleculeList = []

	# Sixth, if solvents have been given, initalise a list to store which molecules are solvent. 
	if len(solvent_components) > 0:
		SolventsList = []

	# Seventh, if the graphs of the molecules have been given, set up the objects for storing the crystal graph. 
	if copied_molecule_graphs is not None:
		crystal_graph = Graph()
		all_molecule_graphs_for_individual_molecules = {}

	# Eighth, make a counter for recording which molecule we are recording in the crystal.
	molecule_name_in_crystal = 1

	# Ninth, make a counter for recording the current number of atoms recorded in the crystal.
	current_total_number_of_atoms = 0

	identity_molecules = []

	# Tenth, perform each symmetry_operation upon the main component of the crystal to create the full crystal. 
	for symmetry_index, symmetry_operation in enumerate(symmetry_operations):

		# 10.1: For each molecule in copied_molecules.
		for molecule_name, copied_molecule in sorted(copied_molecules.items()):

			# 10.2: Create a new molecule that is a copy of the original molecule, and set it up to be added to the crystal.
			new_molecule = copied_molecule.copy()
			new_molecule.set_cell(crystal_cell_lattice)
			new_molecule.set_pbc(True)

			# 10.3: Check that the molecule contains atoms
			if len(new_molecule) == 0:
				raise Exception('new_molecule has no atoms in it?')

			# 10.4: Perform the symmetry operation upon the molecule (if it is not the identity).
			if symmetry_index != 0:

				# 10.4.1: Wrap the molecule within the crystal unit cell.
				new_molecule.wrap()

				# 10.4.2: Perform the symmetry operation given by symmetry_operation upon new_molecule.
				new_molecule = perform_symmetry_operation_upon_molecule(new_molecule, symmetry_operation)

				# 10.4.3: Make sure the operation still has a wrapped molecule in the crystal after performing the crystal operation. 
				new_molecule.wrap()

			# 10.5: Reconnect the molecule after performing the symmetry operation.
			if not wrap:
				new_molecule, molecule_graph = unwrap_molecule(new_molecule, copied_molecule_graphs, molecule_name)

			# 10.6: Translate the new_molecule so that the centre of molecule lies inside the unit cell
			add_translation = get_translation_to_move_COM_inside_unit_cell(new_molecule, crystal_cell_lattice)
			new_molecule.set_positions(new_molecule.get_positions() + add_translation)

			# 10.7: Check to see if this molecule hasn't already been obtained from a previous symmetry operation.
			if determine_is_molecule_already_recorded(new_molecule, all_molecules_in_crystal, crystal_cell_lattice, super_cell_reach=2, include_hydrogen=False, consider_elements=True, return_same_molecule_details=False):
				#if (symmetry_index != 0) and (molecule_name not in identity_molecules):
				#	import pdb; pdb.set_trace()
				continue

			#if symmetry_index == 0:
			#	identity_molecules.append(molecule_name)

			# 10.8: Add this new_molecule to the lists of all_molecules_in_crystal.
			all_molecules_in_crystal[molecule_name_in_crystal] = new_molecule

			# 10.9: Add a unique molecule_name_in_crystal to the MoleculeList to identify this separate molecule in the crystal.
			MoleculeList += [molecule_name_in_crystal]*len(new_molecule)

			# 10.10: If the atoms in the crystal are all the same tag, set Tags to None
			if all(new_molecule.get_tags() == 0):
				new_molecule.set_tags(None)

			# 10.11: If their are solvents in the crystal, record which molecules are solvent molecules.
			if len(solvent_components) > 0:
				if molecule_name in solvent_components:
					SolventsList.append(molecule_name_in_crystal)

			# 10.12: Add the graph for this molecule to the crystal_graph.
			if copied_molecule_graphs is not None:

				# 10.12.1: Make a copy of the molecule graph
				molecule_graph_to_add = deepcopy(copied_molecule_graphs[molecule_name])

				# 10.12.2: Obain the number of the nodes in the molecule of interest (Nodes are atoms here).
				nodes_list = tuple(molecule_graph_to_add.nodes)

				# 10.12.3: Obtain how to convert the indices of the molecule in the molecule to that in the crystal
				mapping = {nodes_list[atom_index]: current_total_number_of_atoms+atom_index for atom_index in range(len(nodes_list))}

				# 10.12.4: Update the names of the nodes in the molecule graph so the indices of atoms in the molecule are those of the same atoms indices in the crystal
				molecule_graph_to_add = relabel_nodes(molecule_graph_to_add, mapping)

				# 10.12.5: Add the molecule_graph_to_add to the currently being built crystal_graph
				crystal_graph = compose(crystal_graph, molecule_graph_to_add)

				# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

				# 10.12.6: Make sure that molecule_name_in_crystal is not already in all_molecule_graphs_for_individual_molecules
				if molecule_name_in_crystal in all_molecule_graphs_for_individual_molecules:
					raise Exception(f'Error: Molecule name {molecule_name_in_crystal} already in all_molecule_graphs_for_individual_molecules. This is probably a programming error. Check this.')

				# 10.12.7: Add to the molecule graph to the all_molecule_graphs_for_individual_molecules list.
				all_molecule_graphs_for_individual_molecules[molecule_name_in_crystal] = deepcopy(copied_molecule_graphs[molecule_name])

				# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

			# 10.13: Increment the molecule_name number by 1.
			molecule_name_in_crystal += 1

			# 10.14: Increase the current_total_number_of_atoms counter by the number of atoms in this molecule.
			current_total_number_of_atoms += len(copied_molecule)

	# ==================================================================================================================================

	# Eleventh, initialise the dictionary for recording the symmetric molecules due to the crystal symmetry.
	equivalent_molecule_groups_due_to_crystal_symmetry           = {} #molecule_name : [] for molecule_name, copied_molecule in sorted(copied_molecules.items())}
	equivalent_molecule_groups_due_to_crystal_symmetry_molecules = {}

	# Twelfth, check if any of the molecules are the same through crystal symmetries.
	for molecule_name, molecule in all_molecules_in_crystal.items():

		# 12.1: Create a tag to note if molecule_name is equivalent to a molecule in equivalent_molecule_groups_due_to_crystal_symmetry after a crystal symmetry operation.
		found_crystal_symmetric_molecule = False

		# 12.2: For ->
		#       * each unique molecule in the equivalent_molecule_groups_due_to_crystal_symmetry dictionary, and
		#       * each symmetry operation in symmetry_operations.
		for cry_sym_unique_mol_name, (symmetry_index, symmetry_operation) in itertools.product(equivalent_molecule_groups_due_to_crystal_symmetry.keys(), enumerate(symmetry_operations)):

			# 12.3: Make a copy of the molecule.
			molecule_copy = molecule.copy()
			molecule_copy.set_cell(crystal_cell_lattice)
			molecule_copy.set_pbc(True)

			# 12.4: Perform the symmetry operation upon the molecule (if it is not the identity).
			if symmetry_index != 0:

				# 12.4.1 Wrap the molecule within the crystal unit cell.
				molecule_copy.wrap()

				# 12.4.2: Perform the symmetry operation given by symmetry_operation upon molecule_copy.
				molecule_copy = perform_symmetry_operation_upon_molecule(molecule_copy, symmetry_operation)

				# 12.4.3: Make sure the operation still has a wrapped molecule in the crystal after performing the crystal operation. 
				molecule_copy.wrap()

			# 12.5: Reconnect the molecule after performing the symmetry operation.
			if not wrap:
				molecule_copy, molecule_graph_copy = unwrap_molecule(molecule_copy, all_molecule_graphs_for_individual_molecules, molecule_name)

			# 12.6: Translate molecule_copy so that the centre of molecule lies inside the unit cell.
			add_translation = get_translation_to_move_COM_inside_unit_cell(molecule_copy, crystal_cell_lattice)
			molecule_copy.set_positions(molecule_copy.get_positions() + add_translation)

			# 12.7: Determine if molecule_name is the as any of the molecules in equivalent_molecule_groups_due_to_crystal_symmetry_molecules.
			is_molecule_already_recorded, same_molecule_name, same_molecule_cell_point = determine_is_molecule_already_recorded(molecule_copy, equivalent_molecule_groups_due_to_crystal_symmetry_molecules, crystal_cell_lattice, super_cell_reach=2, include_hydrogen=False, consider_elements=True, return_same_molecule_details=True)

			# 12.8: Check to see if this molecule hasn't already been obtained from a previous symmetry operation.
			if is_molecule_already_recorded:

				# 12.8.1: If here, molecule_name can be obtained via a crystal symmetry operation upon a molecule in the crystal with the name same_molecule_name.
				equivalent_molecule_groups_due_to_crystal_symmetry[same_molecule_name].append(molecule_name)

				# 12.8.2: Set found_crystal_symmetric_molecule to True, as we have found a molecule in equivalent_molecule_groups_due_to_crystal_symmetry that is the same as molecule_name after a crystal symmetry operation. 
				found_crystal_symmetric_molecule = True

				# 12.8.2: Break out of this for loop, as we have been able to identify that molecule_name can be obtained by performing a crystal operation upon same_molecule_name.
				break

		else:

			# 12.9.1: If we get to this point, then molecule_name can not be obtained from another molecule in the crystal using a crystal symmetry operation. So it is unique. 
			equivalent_molecule_groups_due_to_crystal_symmetry[molecule_name] = []

			# 12.9.2: Make a copy of this molecule and add it to the equivalent_molecule_groups_due_to_crystal_symmetry_molecules dictionary. 
			equivalent_molecule_groups_due_to_crystal_symmetry_molecules[molecule_name] = molecule.copy()

	# ==================================================================================================================================
	# Ninth, perform several check upon the molecules in the crystal.

	'''
	# 9.1: Check to make sure that the molecules in all_molecules_in_crystal are named consecutively from 1 to len(all_molecules_in_crystal)
	if not (sorted(all_molecules_in_crystal.keys()) == list(range(1,len(all_molecules_in_crystal)+1))):
		to_string  = f'Error: The molecules in all_molecules_in_crystal are not consecutively ordered from 1 to len(all_molecules_in_crystal) (which is {len(all_molecules_in_crystal)}).\n'
		to_string += f'Molecule names in all_molecules_in_crystal: {sorted(all_molecules_in_crystal.keys())}\n'
		to_string += f'Expected molecule names: {list(range(1,len(all_molecules_in_crystal)+1))}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# 9.2: Check to make sure that the graphs in all_molecule_graphs_for_individual_molecules are named consecutively from 1 to len(all_molecule_graphs_for_individual_molecules)
	if molecule_graphs is not None:
		if not (sorted(all_molecule_graphs_for_individual_molecules.keys()) == list(range(1,len(all_molecule_graphs_for_individual_molecules)+1))):
			to_string  = f'Error: The molecule graphs in all_molecule_graphs_for_individual_molecules are not consecutively ordered from 1 to len(all_molecule_graphs_for_individual_molecules) (which is {len(all_molecule_graphs_for_individual_molecules)}).\n'
			to_string += f'Molecule names in all_molecule_graphs_for_individual_molecules: {sorted(all_molecule_graphs_for_individual_molecules.keys())}\n'
			to_string += f'Expected molecule names: {list(range(1,len(all_molecule_graphs_for_individual_molecules)+1))}\n'
			to_string += 'Check this.'
			raise Exception(to_string)
	'''

	# 9.3: Check to make sure all the molecules have an associated graph if molecule_graphs not None. 
	if molecule_graphs is not None:
		if not (sorted(all_molecules_in_crystal.keys()) == sorted(all_molecule_graphs_for_individual_molecules.keys())):
			to_string  = f'Error: Some of the molecules and/or their associated graphs are missing in the crystal\n'
			to_string += f'Molecule names in molecules: {sorted(all_molecules_in_crystal.keys())}\n'
			to_string += f'Molecule names in molecule_graphs: {sorted(all_molecule_graphs_for_individual_molecules.keys())}\n'
			to_string += 'Check this.'
			raise Exception(to_string)

	# ==================================================================================================================================

	# Eleventh, set up the ase.Atoms file that will work as the crystal.
	crystal = Atoms(cell=crystal_cell_lattice, pbc=True)
	for molecule_name in sorted(all_molecules_in_crystal.keys(),reverse=False):
		crystal += all_molecules_in_crystal[molecule_name].copy()

	# Twelfth, add MoleculeList to crystal graph nodes
	if copied_molecule_graphs is not None:
		add_node_properties_to_graph(crystal_graph, 'MoleculeList', MoleculeList)
	else:
		crystal.set_array('MoleculeList', np.array(MoleculeList))

	# Thirteenth, wrap the crystal object so that all atoms lie inside the unit cell.
	if wrap:
		crystal.wrap()

	# Fourteenth, add information about which molecules are the same due to symmetry in the crystal. 
	crystal.info['CrystalSymmetryOperations'] = ' '.join(get_symmetry_operators_string(symmetry_operations))
	crystal.info['SameMoleculesDueToCrystalSymmetry'] = str(equivalent_molecule_groups_due_to_crystal_symmetry).replace(' ','').replace('[],',',').replace('{','').replace('}','')

	# Fiftheenth, if their are solvents in the crystal, add a tag in the info that includes which molecules are solvents
	if len(solvent_components) > 0:
		crystal.info['SolventsList'] = ' '.join(str(index) for index in SolventsList)

	# Sixteenth, remove momenta, as atoms in crystal structures don't have momentum
	if 'masses' in crystal.arrays:
		del crystal.arrays['masses']

	# Seventeenth, remove momenta, as atoms in crystal structures don't have momentum
	if 'momenta' in crystal.arrays:
		del crystal.arrays['momenta']

	# Eighteenth, remove initial_magmoms is the magmom for each atom is 0
	if ('initial_magmoms' in crystal.arrays) and all([(round(magmom,8) == 0.0) for magmom in crystal.arrays['initial_magmoms'].tolist()]):
		del crystal.arrays['initial_magmoms']

	# Sixteenth, return crystal and crystal_graph
	return crystal, crystal_graph

# ======================================================================================================================================
# ======================================================================================================================================

def check_copied_graphs(copied_molecules, copied_molecule_graphs):
	"""
	This method is designed to check that the copied molecules and associated graphs are consistent with each other

	Parameters
	----------
	copied_molecules : list of ase.Atoms
		This is a list of molecules as ase.Atoms objects
	copied_molecule_graphs : networkx.Graph or None
		These are the graphs that go with each molecule in the molecules list. If None given, the molecule_graphs will not be used in this method. Default: None.
	"""

	# First, check that the number of graphs given for each molecule is the name as the number of molecules.
	if not len(copied_molecules) == len(copied_molecule_graphs):
		to_string  = 'Error: The number of molecules given is not the same as the number of molecule graphs given.\n'
		to_string += 'No of Molecules: '+str(len(copied_molecules))+'\n'
		to_string += 'No of Molecule Graphs: '+str(len(copied_molecule_graphs))+'\n'
		to_string += 'Check this.'
		raise Exception(to_string)
	
	# Second, make sure that the name of the molcules go with that of ther graphs in copied_molecule_graphs
	if not sorted(copied_molecules.keys()) == sorted(copied_molecule_graphs.keys()):
		to_string  = 'Error: There are missing molecules and/or molecule graphs.\n'
		to_string += f'Molecule names: {sorted(copied_molecules.keys())}\n'
		to_string += f'Molecule graph names: {sorted(copied_molecule_graphs.keys())}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Third, Check the number of atoms in each molecule is the same for each graph
	for name in copied_molecules.keys():
		if not len(copied_molecules[name]) == len(copied_molecule_graphs[name]):
			to_string  = 'Error: The number of atom in the molecule is not the same as in the molecule_graph.\n'
			to_string += 'Number of atoms in molecules['+str(name)+']: '+str(len(copied_molecules[name]))+'\n'
			to_string += 'Number of atoms in molecule_graphs['+str(name)+']: '+str(len(copied_molecule_graphs[name]))+'\n'
			to_string += 'Check this.'
			raise Exception(to_string)

def remove_solvent_from_molecule_dicts(solvent_components, copied_molecules, copied_molecule_graphs=None):
	"""
	This method is designed to remove the solvent molecules from solvent_components from copied_molecules and copied_molecule_graphs.

	This method will make modifications to copied_molecules and copied_molecule_graphs.

	Parameters
	----------
	solvent_components : list of ints
		This is the list of molecule names in the crystal that are solvents that we dont want to include in the crystal structure.
	copied_molecules : dict. of ase.Atoms
		These are the molecules as ase.Atoms objects
	copied_molecule_graphs : networkx.Graph or None
		These are the graphs that go with each molecule in the molecules list. If None given, the molecule_graphs will not be used in this method. Default: None.
	"""

	# First, for each molecule index in solvent_components (given in reverse order).
	for solvent_name in sorted(solvent_components.keys(), reverse=True):

		# 2.1: Delete the molecule from the copied_molecules list
		del copied_molecules[solvent_name]

		# 2.2: If copied_molecule_graphs exists, delete the molecule's graph from the copied_molecule_graphs list
		if copied_molecule_graphs is not None:
			del copied_molecule_graphs[solvent_name]

	# Second, gather the mapping for changing the names of the molecules after removing the solvents from the crystal.
	mapping = {old_name: new_name for old_name, new_name in zip(sorted(solvent_components.keys(),reverse=False), range(1,len(solvent_components)+1))}

	# Third, convert the name of the molecules in copied_molecules from old_name --> new_name
	for old_name, new_name in sorted(mapping.items(), reverse=True): 
		molecule = copied_molecules.pop(old_name)
		if new_name in solvent_components.keys():
			raise Exception(f'Error: {new_name} is already in solvent_components. There is a programming issue. Check this.')
		copied_molecules[new_name] = molecule

	# Fourth, convert the name of the moelcule graphs in copied_molecule_graphs from old_name --> new_name
	if copied_molecule_graphs is not None:
		for old_name, new_name in sorted(mapping.items(), reverse=True): 
			molecule_graph = copied_molecule_graphs.pop(old_name)
			if new_name in copied_molecule_graphs.keys():
				raise Exception(f'Error: {new_name} is already in copied_molecule_graphs. There is a programming issue. Check this.')
			copied_molecule_graphs[new_name] = molecule_graph

	'''
	# Fifth, check to make sure that the molecules in copied_molecules are named consecutively from 1 to len(copied_molecules)
	if sorted(copied_molecules.keys()) == list(range(1,len(copied_molecules)+1)):
		to_string  = f'Error: The molecules in copied_molecules are not consecutively ordered from 1 to len(copied_molecules) (which is {len(copied_molecules)}).\n'
		to_string += f'Molecule names in copied_molecules: {sorted(copied_molecules.keys())}\n'
		to_string += f'Expected molecule names: {list(range(1,len(copied_molecules)+1))}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Sixth, check to make sure that the graphs in copied_molecule_graphs are named consecutively from 1 to len(copied_molecule_graphs)
	if sorted(copied_molecule_graphs.keys()) == list(range(1,len(copied_molecule_graphs)+1)):
		to_string  = f'Error: The molecule graphs in copied_molecule_graphs are not consecutively ordered from 1 to len(copied_molecule_graphs) (which is {len(copied_molecule_graphs)}).\n'
		to_string += f'Molecule names in copied_molecule_graphs: {sorted(copied_molecule_graphs.keys())}\n'
		to_string += f'Expected molecule names: {list(range(1,len(copied_molecule_graphs)+1))}\n'
		to_string += 'Check this.'
		raise Exception(to_string)
	'''

	import pdb; pdb.set_trace()
	raise Exception('Check that everything is working as expected')

def unwrap_molecule(new_molecule, copied_molecule_graphs, name):
	"""
	This will unwrap the molecule, meaning atoms that are bonded to eachother will be put back together.

	This method will allow atoms to be found outside the origin unit cell.

	Parameters
	----------
	new_molecule : ase.Atoms
		This is the molecule of interest.
	copied_molecule_graphs : list of networkx.Graph or None
		These are the graphs for all the molecules in the crystal. If None given, a custom graph will be obtained for new_molecule. Default: None.
	name : int
		This is the name of new_molecule in copied_molecule_graphs. This is used to get the graph of new_molecule from copied_molecule_graphs.

	Return
	------
	new_molecule : ase.Atoms
		This is the unwrapped, reconstructed molecule of interest. Atoms have been moved so the bonded atoms are all neighbouring each other. A molecule may have part of it's structure located outside of the origin unit cell.
	molecule_graph : networkx.Graph
		This is the graph for new_molecule
	"""

	# First, get the graph for the molecule.
	if copied_molecule_graphs is not None:
		molecule_graph = copied_molecule_graphs[name]
	else:
		molecule_graph = obtain_graph(new_molecule)
		import pdb; pdb.set_trace()
		raise Exception('Check this out.')

	# Second, reconnect the molecule after performing the symmetry operation.
	atom_indices = list(range(len(new_molecule)))
	new_molecule, molecule_graph = make_molecule(atom_indices, molecule_graph, new_molecule, take_shortest_distance=True, make_molecule_method='component_assembly_approach', logger=None)
	
	# Third, return new_molecule and molecule_graph
	return new_molecule, molecule_graph

# ======================================================================================================================================
