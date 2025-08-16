"""
process_crystal_with_moleculeList.py, Geoffrey Weal, 2/6/2022

This script is designed to process the crystal given into the molecules it is made up. 

This script takes advantage of the MoleculeList attached to the crystal object.
"""
import sys
import numpy as np
from copy import deepcopy
from tqdm import trange, tqdm

import multiprocessing as mp
from tqdm.contrib.concurrent import process_map

from networkx import Graph
from ase import Atoms, Atom
from ase.spacegroup import Spacegroup

from SUMELF.SUMELF.graph_methods.get_properties_from_graph               import get_node_properties_from_graph
from SUMELF.SUMELF.utility_methods.get_SameMoleculesDueToCrystalSymmetry import get_SameMoleculesDueToCrystalSymmetry
from SUMELF.SUMELF.general_methods.get_symmetry_operations               import get_symmetry_operations
from SUMELF.SUMELF.process_crystal_methods.unwrap_molecules              import unwrap_molecules
from SUMELF.SUMELF.process_crystal_methods.utilities                     import get_solvent_molecule_names_from_SolventsList

def process_crystal_with_moleculeList(crystal, crystal_graph, bonds_to_ignore=None, produce_all_molecules=False, take_shortest_distance=False, make_molecule_method='component_assembly_approach', no_of_cpus=1, logger=None, print_progress=True):
	"""
	This method is designed to process the crystal given into the molecules it is made up. 

	This method takes advantage of the MoleculeList attached to the crystal object.

	Parameters
	----------
	crystal : ase.Atoms
		This is ase object of the crystal. 
	crystal_graph : networkx.graph
		This is the crystal_graph that has already been processed. If you want this method to process the graph for this crystal, set crystal_graph to None.
	bonds_to_ignore : list
		This is a list of the atom pairs (given as indices) to ignore bonding between. This is given as a list of index pairs. Default: None
	produce_all_molecules : bool.
		This boolean indicates if you want to return all the molecules in the crystal, or just the ones that can be replicated using the symmetry operations. Default: False
	take_shortest_distance : bool
		If true, take the shortest distance that was found if a distance shorter than the maximum distance between elements if found. If false, give an error if the distance between two atoms was greater than the max distance expected for those two elements to be bonded to each other. Default: False
	make_molecule_method : str.
		This is the name of the method you want to use to create the molecule. See https://github.com/geoffreyweal/ECCP for more information. Default: 'component_assembly_approach'. 
	no_of_cpus : int.
		This is the number of cpus available to use on this program. In most cases this should just be set to 1 cpu, however for very large system you may want to implement multiple cpus.
	logger : logging/None
		This contains the object for reporting information and warning messages to.
	print_progress : bool.
		Indicates if you want to write what the method is currently doing to the terminal. Default: True

	Returns
	-------
	molecules : list of ase.Atoms
		These are the molecules that make up the crystal (before any spacegroup symmetry is considered).
	molecule_graphs : list of networkx.Graph
		These are the graph that are associated with the molecules in the molecules list. 
	SolventsList : list of ints
		This list indicates which molecules in the molecules list are solvents.
	symmetry_operations : list of (numpy.array(3*3), numpy.array(3*1))
		These are the rotation/reflection and translation matrices for each symmetry operation the crystal contans. 
	"""

	# First, obtain the MoleculeList from crystal_graph
	MoleculeDict = get_node_properties_from_graph(crystal_graph, 'MoleculeList')

	# Second, obtain the indices of atoms in each molecule from the crystal.
	molecule_names = {}
	for atom_index, mol_name in MoleculeDict.items():
		molecule_names.setdefault(mol_name,[]).append(atom_index)

	# ================================================================================================================

	# Third, make sure that all the molecules names are greater than or equal to 1
	if any([(mol_name < 1) for mol_name in molecule_names.keys()]):
		to_string  = 'Error: A molecule given in the crystal object has been labelled a value less than 0.\n'
		to_string += 'All molecules should have a molecule name that is an integer than is greater than 1.\n'
		to_string += f'Molecule names: {sorted(set(molecule_names.keys()))}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# ================================================================================================================
	# Fourth, create all the ASE Atoms objects for molecules in the crystal.

	# 4.1: Initialse the dictionary to record which atom indices in the crystal go with each atom in each molecule. 
	molecule_to_crystal_mappings = {mol_name: {} for mol_name in sorted(molecule_names.keys())}
	
	# 4.2: Obtain the molecules from the crystal.
	for atom_index_in_crystal in range(len(crystal_graph)):

		# 4.2.1: Obtain the molecule index for this atom
		mol_name = MoleculeDict[atom_index_in_crystal]

		# 4.2.2: Obtain the index for this atom in the molecule
		atom_index_in_molecule = len(molecule_to_crystal_mappings[mol_name])

		# 4.2.4: Map each atom in each molecule 
		molecule_to_crystal_mappings[mol_name][atom_index_in_molecule] = atom_index_in_crystal

	# ================================================================================================================
	# Fifth, if the "SameMoleculesDueToCrystalSymmetry" list is given in the crystal file, 
	#        we can use this information to remove some of these molecules from the molecules dictionary. 
	#
	# Note: If a molecule is not in the keys in "SameMoleculesDueToCrystalSymmetry" (and produce_all_molecules is set to True), it is an equivalent molecule.

	# 5.1: If we only want to look at unique molecules, do the following
	if not produce_all_molecules:

		# 5.2: Determine if "SameMoleculesDueToCrystalSymmetry" exists in the crystal file.
		has_SameMoleculesDueToCrystalSymmetry = ('SameMoleculesDueToCrystalSymmetry' in crystal.info)

		# 5.3: Remove all symmetrically similar structures from the molecule.
		if has_SameMoleculesDueToCrystalSymmetry:

			# 5.4: Obtain the dictionary of equivalent molecules and their unique counterparts
			SameMoleculesDueToCrystalSymmetry = get_SameMoleculesDueToCrystalSymmetry(crystal.info['SameMoleculesDueToCrystalSymmetry'])

			# 5.5: For each entry in SameMoleculesDueToCrystalSymmetry
			for unique_mol_name, list_of_equivalent_mol_names in sorted(SameMoleculesDueToCrystalSymmetry.items(), key=lambda x:x[0], reverse=True):

				# 5.6: For each molecule in list_of_equivalent_mol_names:
				for equivalent_mol_name in sorted(list_of_equivalent_mol_names, reverse=True):

					# 5.7: If same_mol_name > mol_name, their is a programming error.
					if unique_mol_name >= equivalent_mol_name:
						to_string  = 'Error: Unique molecule has higher or equal index than equivalent.\n'
						to_string += 'This program is designed so the unique molecule has the lowest index.\n'
						to_string += 'Check this: Equivalent molecule index (should be higher): '+str(equivalent_mol_name)+'\n'
						to_string += 'Unique molecule index (should be lower): '+str(unique_mol_name)+'\n'
						raise Exception(to_string)
					
					# 5.8: Remove the equivalent molecule from molecule_to_crystal_mappings, as we dont want to record it.
					del molecule_to_crystal_mappings[equivalent_mol_name]

	# ================================================================================================================
	# Sixth, make a reverse dictionary for looking up molecule index and atom index in that molecule from the atom index in crystal.

	# 6.1: Initialise a dictionary that will allow the user to make an atom in the crystal to the atom in the molecule it has been mapped to.
	crystal_mappings = {}

	# 6.2: Obtain the mapping from crystal atom index --> (molecule number, molecule atom index)
	for mol_name, atom_mapping in molecule_to_crystal_mappings.items(): 
		for atom_index_in_mol, atom_index_in_crystal in sorted(atom_mapping.items(), key=lambda x:x[0]):
			crystal_mappings[atom_index_in_crystal] = (mol_name, atom_index_in_mol)

	# ================================================================================================================
	# Seventh, obtain the ASE Atoms objects for each molecule in the crystal
	
	# 7.1: Set up each of the ase objects for each molecule in the crystal
	molecules = {mol_name: Atoms(cell=crystal.get_cell(),pbc=True) for mol_name in molecule_to_crystal_mappings.keys()}

	# 7.2: Obtain the molecules from the crystal.
	for atom_index, crystal_atom in enumerate(crystal):

		# 7.2.1: Obtain the molecule index for this atom
		mol_name = MoleculeDict[atom_index]

		# 7.2.2: If mol_name is not in molecule_to_crystal_mappings.keys() (and produce_all_molecules is set to True), 
		#        We dont want to process the atoms in this equivalent molecule. 
		if mol_name not in molecule_to_crystal_mappings.keys():
			continue

		# 7.2.3: Copy the atom from Crystal for use in the molecule.
		mol_atom = Atom(symbol=crystal_atom.symbol, position=crystal_atom.position, tag=crystal_atom.tag, momentum=crystal_atom.momentum, mass=crystal_atom.mass, magmom=crystal_atom.magmom, charge=crystal_atom.charge)

		# 7.2.4: Append the atom to the molecule.
		molecules[mol_name].append(mol_atom)

	# ================================================================================================================
	# Eighth, obtain the graph for each molecule in the crystal
	
	# 8.1: Initialise a dictionary to record all the node information for each atom in each molecule. 
	molecule_graphs_nodes = {mol_name: [] for mol_name in molecule_to_crystal_mappings.keys()}

	# 8.2: Obtain the node information for each atom in each molecule.
	for atom_index_crystal, atom_properties in sorted(crystal_graph.nodes.items(), key=lambda x:x[0]):

		# 8.2.1: Check if atom_index_crystal is in crystal_mappings.
		#        If it is not, it is assumed to be in an equivalent molecule, and we dont want to process it.
		if atom_index_crystal not in crystal_mappings:
			continue

		# 8.2.2: Get the molecule name and atom index in the molecule for the atom in the crystal.
		mol_name, atom_index_in_molecule = crystal_mappings[atom_index_crystal]

		# 8.2.3: Make a copy of the atom_properties
		copy_of_atom_properties = dict(atom_properties)

		# 8.2.4: Remove MoleculeList from copy_of_atom_properties
		if not 'MoleculeList' in copy_of_atom_properties:
			to_string  = 'Error: MoleculeList may not be in crystal_graph.nodes['+str(atom_index_crystal)+'].'
			to_string += 'This should not be the case, and if it is there is a programming issue.\n'
			to_string += 'crystal_graph.nodes['+str(atom_index_crystal)+'] = '+str(crystal_graph.nodes[atom_index_crystal])
			raise Exception(to_string)
		del copy_of_atom_properties['MoleculeList']

		# 8.2.5: Assign the atom properties to the correct atom in the molecule in molecule_graphs_nodes
		molecule_graphs_nodes[mol_name].append((atom_index_in_molecule, copy_of_atom_properties))

	# 8.3: Initialise a dictionary to record all the edge information for each bond in each molecule. 
	molecule_graphs_edges = {mol_name: [] for mol_name in molecule_to_crystal_mappings.keys()}

	# 8.4: Obtain the edge information for each bond in each molecule.
	for (atom_index1_crystal, atom_index2_crystal), bond_properties in sorted(crystal_graph.edges.items(), key=lambda x:x[0]):

		# 8.4.1: Check if atom_index1_crystal and atom_index2_crystal are not in crystal_mappings.
		#        If they are not, it is assumed they are located in an equivalent molecule, and we dont want to process it.
		if (atom_index1_crystal not in crystal_mappings.keys()) and (atom_index2_crystal not in crystal_mappings.keys()):
			continue

		# 8.4.2: If only one of atom_index1_crystal or atom_index2_crystal is found in crystal_mappings, there is a programming error
		if (atom_index1_crystal not in crystal_mappings) or (atom_index2_crystal not in crystal_mappings):
			toString  = 'Error: An atom found in crystal_mappings is bonded to an atom not found in crystal_mappings.\n;'
			toString += 'This indicates a programming error, as one atom is in a unique molecule while the rother atom in the bond is in a different equivalent molecule.\n'
			toString += 'Main big problem is that this indicate two atoms in different molecules are bonded together.\n'
			toString += 'This should not be the case, every atom in a bond should be in the same molecule. Check this.\n\n'
			toString += 'index1 in the bond: '+str(atom_index1_crystal)+'. Found in crystal_mappings: '+str(atom_index1_crystal not in crystal_mappings)+'\n'
			toString += 'index2 in the bond: '+str(atom_index2_crystal)+'. Found in crystal_mappings: '+str(atom_index2_crystal not in crystal_mappings)
			raise Exception(toString)

		# 8.4.3: Get the molecule names and atom indices in the molecule for the atoms in the crystal for this bond.
		mol1_name, atom_index1_in_molecule = crystal_mappings[atom_index1_crystal]
		mol2_name, atom_index2_in_molecule = crystal_mappings[atom_index2_crystal]

		# 8.4.4: Check that mol1_name and mol2_name are the same. They should be as a bond should be contained in the same molecule. 
		if not mol1_name == mol2_name:
			to_string  = 'Error: The atoms in a bond in a crystal has been assigned to different molecules.'+'\n'
			to_string += ' Check:\t'+'\n'
			to_string += '\t'+str(atom_index1_in_molecule)+'-->'+str(mol1_name)+'\n'
			to_string += '\t'+str(atom_index2_in_molecule)+'-->'+str(mol2_name)
			raise Exception(to_string)

		# 8.4.5: Check that atom_index1_in_molecule is not the same as atom_index2_in_molecule.
		if atom_index1_in_molecule == atom_index2_in_molecule:
			raise Exception('Error: atom '+str(atom_index1_in_molecule)+' in molecule '+str(mol1_name)+' seems to be bonded to itself?')

		# 8.4.6: Create the new bond with updated indices for this molecule, and append this bond to the molecule of interest.
		molecule_graphs_edges[mol1_name].append((atom_index1_in_molecule, atom_index2_in_molecule, dict(bond_properties)))

	# 8.5: Create each graph for each molecule
	molecule_graphs = {}
	for mol_name in molecule_to_crystal_mappings.keys():

		# 8.5.1: Initialise the graph.
		molecule_graph = Graph(name=mol_name)

		# 8.5.2: Add atoms (nodes) to the graph.
		molecule_graph.add_nodes_from(molecule_graphs_nodes[mol_name])

		# 8.5.3: Add bonds (edges) to the graph.
		molecule_graph.add_edges_from(molecule_graphs_edges[mol_name])

		# 8.5.4: Add this graph to molecule_graphs.
		molecule_graphs[mol_name] = molecule_graph

	# ================================================================================================================

	# Ninth, get the symmetry operations of the spacegroup for the crystal. 
	#        If the spacegroup is not given, it is assumed to be 'P 1' (i.e. all atoms in the crystal has been given in the unit cell.)
	if not produce_all_molecules:
		has_SymOpts    = 'CrystalSymmetryOperations' in crystal.info.keys()
		has_spacegroup = 'spacegroup'                in crystal.info.keys()
		if   has_SymOpts:
			symmetry_operations = get_symmetry_operations(crystal.info['CrystalSymmetryOperations'].split())
		elif has_spacegroup:
			spacegroup = Spacegroup(spacegroup=crystal.info['spacegroup'])
			symmetry_operations = spacegroup.get_symop()
		else:
			# Assume P 1 symmetry 
			spacegroup = Spacegroup(spacegroup='P 1')
			symmetry_operations = spacegroup.get_symop()
	else:
		# Will return all molecules in the crystal regardless of symmetry, so symmetry operation given is the identity.
		spacegroup = Spacegroup(spacegroup='P 1')
		symmetry_operations = spacegroup.get_symop()

	# ================================================================================================================
	# Tenth, It is very likely that the molecule is wrapped inside of the unit cell. 
	#        Here, we unwrap the molecule and allow parts of it to stick out of the cell to make it easier to view.
	molecules = unwrap_molecules(molecules, molecule_graphs, take_shortest_distance, make_molecule_method, no_of_cpus, logger, print_progress=print_progress)

	# Eleventh, obtain the list of molecules that are solvent in molecules
	solvent_molecule_names = get_solvent_molecule_names_from_SolventsList(crystal, molecules)

	# ================================================================================================================

	# Twelfth, make sure that the molecules have an associated graph
	if not sorted(molecules.keys()) == sorted(molecule_graphs.keys()):
		import pdb; pdb.set_trace()
		to_string  = 'Error: There are missing molecules and/or molecule graphs.\n'
		to_string += f'Molecule names:       {sorted(molecule_names.keys())}\n'
		to_string += f'Molecule Graph names: {sorted(molecule_graphs.keys())}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# Thirteenth, make sure that solvents in solvent_molecule_names can be found in molecules:
	missing_solvent_mol_names = [solvent_mol_name for solvent_mol_name in solvent_molecule_names if (solvent_mol_name not in molecules.keys())]
	if len(missing_solvent_mol_names) > 0:
		to_string  = 'Error, there are some solvents in the solvent list that are not found in the molecules dictionary\n'
		to_string += f'Missing solvent molecules: {sorted(missing_solvent_mol_names)}\n'
		to_string += f'Recorded solvent molecules: {sorted(solvent_molecule_names)}\n'
		to_string += f'All molecules: {sorted(molecules.keys())}\n'
		to_string += 'Check this.'
		raise Exception(to_string)

	# ================================================================================================================

	# Fourteenth, return the molecules and their molecule_graphs
	return molecules, molecule_graphs, solvent_molecule_names, symmetry_operations

# ====================================================================================================================

