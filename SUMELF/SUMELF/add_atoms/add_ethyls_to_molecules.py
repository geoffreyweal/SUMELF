"""
add_ethyls_to_crystal.py, Geoffrey Weal, 3/7/2022

This script is designed to automatically add ethyls to missing sites if possible. 
"""
import os, warnings
from ase.io import write
from copy   import deepcopy

from SUMELF.SUMELF.make_crystal                                                     import make_crystal
from SUMELF.SUMELF.graph_methods.get_properties_from_graph                          import remove_node_properties_from_graph
from SUMELF.SUMELF.add_atoms.add_ethyls_to_molecules_methods.add_ethyls_to_molecule import add_ethyls_to_molecule

def add_ethyls_to_molecules(no_of_ethyls_to_add_to_mols, molecules, molecule_graphs, original_crystal, original_crystal_graph, symmetry_operations, cell, add_hydrogens=True, logger=None, identifier=None):
	"""
	This method is designed to add missing ethyls to crystals. 

	These are ethyls that have been recorded in the crystal but do not have coordinates. These ethyls have been bonded to atoms in the crystal. 

	Parameters
	----------
	no_of_ethyls_to_add_to_mols : dict. of dict. of {int: int}
		This dictionary contains the number of ethyls you want to add to each atom in each molecule in the crystal. 
	molecules : list of ase.Atoms
		These are the molecules in the crystal
	molecule_graphs : list of networkx.graph
		These are the graphs associated with the molecules in the crystal. 
	original_crystal : ase.Atoms
		This is the original crystal
	original_crystal_graph : networkx.graph
		This is the graph associated with the original crystal.
	symmetry_operations : list of (3-by-3 numpy.array, 3-by-1 numpy.array)
		This list contains the rotation/reflection and translation matrices that describe the symmetries present in a crystal structure. 
	cell : numpy.array
		This is the unit cell lattice vectors.
	add_hydrogens : bool.
		this indicates if you want to add hydrogens in this method, or if you will do this later on.
	logger : logging
		This is the log object to report information about the running of the program to.
	identifier : str.
		This is the identifier name of the crystal in the CCDC database

	Returns
	-------
	updated_molecules : list of ase.Atoms
		These are the molecules in the crystal
	updated_molecule_graphs : list of networkx.graph
		These are the graphs associated with the molecules in the crystal. 
	were_ethyls_added : bool.
		This boolean indicates if ethyls were added to this system.
	"""

	# First, check to make sure that the molecules in the molecules dictionary are named consecutively from 1 to len(molecules)
	if not (sorted(molecules.keys()) == list(range(1,len(molecules)+1))):
		to_string  = f'Error: The molecules in the molecules dictionary are not consecutively ordered from 1 to len(molecules) (which is {len(molecules)}).\n'
		to_string += f'Molecule names in molecules: {sorted(molecules.keys())}\n'
		to_string += f'Expected molecule names: {list(range(1,len(molecules)+1))}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Second, check to make sure that the graphs in molecule_graphs are named consecutively from 1 to len(molecule_graphs)
	if not (sorted(molecule_graphs.keys()) == list(range(1,len(molecule_graphs)+1))):
		to_string  = f'Error: The molecule graphs in molecule_graphs are not consecutively ordered from 1 to len(molecule_graphs) (which is {len(molecule_graphs)}).\n'
		to_string += f'Molecule names in molecule_graphs: {sorted(molecule_graphs.keys())}\n'
		to_string += f'Expected molecule names: {list(range(1,len(molecule_graphs)+1))}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Third, make sure that you can all the same molecules and associated graphs
	if not (sorted(molecules.keys()) == sorted(molecule_graphs.keys())):
		to_string  = f'Error: There are some missing molecules and/or molecule graphs.\n'
		to_string += f'Molecule names: {sorted(molecules.keys())}\n'
		to_string += f'Associated molecule graphs: {sorted(molecule_graphs.keys())}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Fourth, check that the atom indices and molecules given in no_of_ethyls_to_add_to_mols are consistent with molecules
	
	# 4.1: Initalise a dictionary that contains which atoms in no_of_ethyls_to_add_to_mols are missing in molecules
	inputs_in__no_of_ethyls_to_add__that_are_missing_in_molecules = {}

	# 4.2: For each molecule given in inputs_in__no_of_ethyls_to_add__that_are_missing_in_molecules:
	for mol_name, no_of_ethyls_to_add_in_mol in no_of_ethyls_to_add_to_mols.items():

		# 4.2.1: Check that the molecule given by mol_name exists in molecules
		if not mol_name in molecules.keys():
			inputs_in__no_of_ethyls_to_add__that_are_missing_in_molecules[mol_name] = 'No molecule'
			continue

		# 4.2.2: For each atom index in no_of_ethyls_to_add_in_mol:
		for atom_index, no_of_ethyls_to_add_to_atom in no_of_ethyls_to_add_in_mol.items():

			#	# 4.2.2.1: Check that atom_index can be found in the molecule
			if   not (0 <= atom_index < len(molecules[mol_name])):
				inputs_in__no_of_ethyls_to_add__that_are_missing_in_molecules.setdefault(mol_name,{})[atom_index] = 'Does not exist'

				# 4.2.2.2: Make sure the user has entered an integer for number of ethyls to add.
			elif not isinstance(no_of_ethyls_to_add_to_atom,int):
				inputs_in__no_of_ethyls_to_add__that_are_missing_in_molecules.setdefault(mol_name,{})[atom_index] = f'Input given not an integer => {no_of_ethyls_to_add_to_atom}.'

				# 4.2.2.3: Make sure the user has not entered in 0 for number of ethyl to add to atom, as this does not make any sense.
			elif no_of_ethyls_to_add_to_atom == 0:
				inputs_in__no_of_ethyls_to_add__that_are_missing_in_molecules.setdefault(mol_name,{})[atom_index] = 'Input given is 0.'

				# 4.2.2.4: Make sure the user has not entered a negative number for number of ethyl to add to atom, as this does not make any sense.
			elif no_of_ethyls_to_add_to_atom < 0:
				inputs_in__no_of_ethyls_to_add__that_are_missing_in_molecules.setdefault(mol_name,{})[atom_index] = f'Input given is an integer number => {no_of_ethyls_to_add_to_atom}.'

	# 4.3: If any problems are given in inputs_in__no_of_ethyls_to_add__that_are_missing_in_molecules, raise an Exception
	if len(inputs_in__no_of_ethyls_to_add__that_are_missing_in_molecules) > 0:
		to_string  = 'Error: Some of the molecules and/or atoms given in the "no_of_ethyls_to_add_to_mols" dictionary have problems in this method.\n'
		to_string += 'These are:\n'
		to_string += '\n'
		to_string += 'molecule name\n'
		to_string += 'Atom index: Reason \n'
		to_string += '------------------\n'
		for mol_name, atom_indices_with_problems in inputs_in__no_of_ethyls_to_add__that_are_missing_in_molecules.items():
			to_string += f'{mol_name}\n'
			for atom_index, reason_for_issue in atom_indices_with_problems.items():
				to_string += f'\t{atom_index}: {reason_for_issue}\n'
			to_string += '\n'
		to_string += 'Check this.\n'
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	# Fifth, make sure that all the atoms in all molecules have the same node and edge properties.

	# 5.1: Get all the names of the nodes found across the graphs in molecule_graphs. 
	node_property_names = set()
	for mol_name, molecule_graph in molecule_graphs.items():
		for node_name, node_property in molecule_graph.nodes.items():
			node_property_names.update(tuple(node_property.keys()))

	# 5.2: Initialise a dictionary to record which atoms do not contain property names.
	node_property_names_in_these_atoms     = {}
	node_property_names_not_in_these_atoms = {}

	# 5.3: Determine which atoms have what properties.
	for mol_name, molecule_graph in molecule_graphs.items():
		for node_name, node_property in molecule_graph.nodes.items():
			for node_property_name in node_property.keys():
				if node_property_name in node_property_names:
					node_property_names_in_these_atoms.setdefault(node_property_name,{}).setdefault(mol_name,[]).append(node_name)
				else:
					node_property_names_not_in_these_atoms.setdefault(node_property_name,{}).setdefault(mol_name,[]).append(node_name)

	# 5.4: Get all the names of the edge found across the graphs in molecule_graphs. 
	edge_property_names = set()
	for mol_name, molecule_graph in molecule_graphs.items():
		for edge_name, edge_property in molecule_graph.edges.items():
			edge_property_names.update(tuple(edge_property.keys()))

	# 5.5: Initialise a dictionary to record which atoms do not contain property names.
	edge_property_names_in_these_atoms     = {}
	edge_property_names_not_in_these_atoms = {}

	# 5.6: Determine which atoms have what properties.
	for mol_name, molecule_graph in molecule_graphs.items():
		for edge_name, edge_property in molecule_graph.edges.items():
			for edge_property_name in edge_property.keys():
				if edge_property_name in edge_property_names:
					edge_property_names_in_these_atoms.setdefault(edge_property_name,{}).setdefault(mol_name,[]).append(edge_name)
				else:
					edge_property_names_not_in_these_atoms.setdefault(edge_property_name,{}).setdefault(mol_name,[]).append(edge_name)

	# 5.7: If there are any atoms or bonds with missing nodes and/or edges, report it here:
	if (len(node_property_names_not_in_these_atoms) > 0) or (len(edge_property_names_not_in_these_atoms) > 0):
		to_string  = 'Error: There are missing node and/or edge properties in some of the atoms in the crystal.\n'
		to_string += 'These are:\n'
		to_string += '\n'
		if len(node_property_names_not_in_these_atoms) > 0:
			to_string += 'Node Properties found in atoms (molecule name: list of atom indices)\n'
			to_string += '--------------------------------------------------------------------\n'
			for mol_name, list_of_atom_indices in node_property_names_in_these_atoms.items():
				to_string += f'{mol_name}: {list_of_atom_indices}\n'
			to_string += '\n'
			to_string += 'Node Properties NOT found in atoms (molecule name: list of atom indices)\n'
			to_string += '------------------------------------------------------------------------\n'
			for mol_name, list_of_atom_indices in node_property_names_not_in_these_atoms.items():
				to_string += f'{mol_name}: {list_of_atom_indices}\n'
			to_string += '\n'
		if len(edge_property_names_not_in_these_atoms) > 0:
			to_string += 'Edge Properties found in atoms (molecule name: list of bonds between atom indices)\n'
			to_string += '----------------------------------------------------------------------------------\n'
			for mol_name, list_of_atom_indices_between_bonds in edge_property_names_in_these_atoms.items():
				to_string += f'{mol_name}: {list_of_atom_indices_between_bonds}\n'
			to_string += '\n'
			to_string += 'Edge Properties NOT found in atoms (molecule name: list of bonds between atom indices)\n'
			to_string += '--------------------------------------------------------------------------------------\n'
			for mol_name, list_of_atom_indices_between_bonds in edge_property_names_not_in_these_atoms.items():
				to_string += f'{mol_name}: {list_of_atom_indices_between_bonds}\n'
			to_string += '\n'
		to_string += 'Check this.\n'
		print(to_string)
		import pdb; pdb.set_trace()
		raise Exception(to_string)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	# Sixth, make a copy of the original crystal and it's associated graph.
	crystal       = original_crystal.copy()
	crystal_graph = deepcopy(original_crystal_graph)

	# Seventh, initalise dictionaries for holding updated molecules and their associated graphs. 
	updated_molecules       = {}
	updated_molecule_graphs = {}

	# Eighth, initialise a boolean for noting if the crystal has been updated with extra ethyl groups
	were_ethyls_added = False

	# Ninth, initalise a dictionary to record all the atoms added to this molecule along with the hydrogens attached to them.
	all_hydrogens_to_add_to_atoms = {}

	# Tenth, for each molecule in the molecules dictionary. 
	for molecule_name in sorted(molecules.keys()):

		# 10.1: Obtain the molecule to focus on at this point.
		molecule                     = molecules[molecule_name]

		# 10.2: Obtain the associated molecule_graph
		molecule_graph               = molecule_graphs[molecule_name]

		# 10.3: Check that there are any ethyl groups to add to molecule_name.
		#      * If there are no atoms to add to molecule_name, move on to the next molecule. 
		if molecule_name in no_of_ethyls_to_add_to_mols.keys():

			# 10.3.1: Add the number of ethyl groups you want to add to the molecule
			no_of_ethyls_to_add_to_mol = no_of_ethyls_to_add_to_mols[molecule_name]

			# 10.3.2: Determine if ethyls need to be added to the molecule and if so, add ethyl groups to the molecule.
			updated_molecule, updated_molecule_graph, were_ethyls_added_to_this_molecule, hydrogens_to_add_to_atoms = add_ethyls_to_molecule(no_of_ethyls_to_add_to_mol, molecule, molecule_graph, crystal, crystal_graph, molecule_name, add_hydrogens=add_hydrogens, logger=logger)

			# 10.3.3: Add either the original molecule or its updated form depending on if ethyl groups were added to the molecule. 
			if not were_ethyls_added_to_this_molecule: # If the molecule has been updated, do the following:
				to_string  = f'Error: No ethyl groups were added to molecule {molecule_name} '
				if identifier is not None:
					to_string += f'(in crystal: {identifier})'
				to_string += '\n'
				to_string += f'This was unexpected, as the following ethyl groups were to be added to this molecule: {no_of_ethyl_to_add_to_mol}\n'
				to_string += 'Check this.\n'
				raise Exception(to_string)

			# Testing
			'''
			warnings.warn('adding molecules to METHOD_add_no_coord_ethyls_to_molecules_test_folder for checking')
			test_folder = 'METHOD_add_no_coord_ethyls_to_molecules_test_folder'
			if not os.path.exists(test_folder):
				os.makedirs(test_folder)
			write(test_folder+'/'+str(identifier)+'_mol_'+str(molecule_name)+'.xyz', [molecule, updated_molecule])
			'''

			# 10.3.4: Add the updated molecule to updated_molecules.
			updated_molecules[molecule_name]       = updated_molecule.copy()

			# 10.3.5: Add the updated molecule graph to updated_molecule_graphs.
			updated_molecule_graphs[molecule_name] = deepcopy(updated_molecule_graph)

			# 10.3.6: Update the crystal.
			crystal, crystal_graph = make_crystal(updated_molecules, symmetry_operations=symmetry_operations, cell=cell, wrap=False, solvent_components=[], remove_solvent=False, molecule_graphs=updated_molecule_graphs)

			# 10.3.7: Add message to logger that non-coordinated were found in the molecule
			if (logger is not None) and (not were_ethyls_added):
				logger.info('Added ethyl group(s) to crystal structure of '+str(identifier))

			# 10.3.8: Set were_ethyls_added to True, as one or more ethyl groups have been added to the crystal. 
			were_ethyls_added = True

			# 10.3.9: Add the atoms that were added to the molecule along with it's hydrogen atom indices for all molecules
			for atom_index, list_of_hydrogen_indices in hydrogens_to_add_to_atoms.items():

				# 10.3.9.1: If molecule_name is not in all_hydrogens_to_add_to_atoms, add it.
				if molecule_name not in all_hydrogens_to_add_to_atoms.keys():
					all_hydrogens_to_add_to_atoms[molecule_name] = {}

				# 10.3.9.2: Check that atom_index is not already in atoms_with_added_hydrogens_indices_for_this_molecule
				if atom_index in all_hydrogens_to_add_to_atoms[molecule_name].keys():
					raise Exception('atom_index in all_hydrogens_to_add_to_atoms[molecule_name]')

				# 10.3.9.3: Add (atom_index, list_of_hydrogen_indices) to atoms_with_added_hydrogens_indices_for_this_molecule
				all_hydrogens_to_add_to_atoms[molecule_name][atom_index] = list_of_hydrogen_indices

		else: # Add the original molecule to updated_molecules

			# 10.3.10: Add the original molecule to updated_molecules.
			updated_molecules[molecule_name]       = molecule.copy()

			# 10.3.11: Add the original molecule graph to updated_molecule_graphs.
			updated_molecule_graphs[molecule_name] = deepcopy(molecule_graph)

	# Eleventh, return updated_molecules, updated_molecule_graphs, and were_ethyls_added.
	return updated_molecules, updated_molecule_graphs, were_ethyls_added, all_hydrogens_to_add_to_atoms

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 





