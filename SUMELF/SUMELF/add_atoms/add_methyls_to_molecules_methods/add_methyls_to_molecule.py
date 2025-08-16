"""
add_methyls_to_molecule.py, Geoffrey Weal, 25/6/2022

This program will allow the user to manually add methyl groupss to atoms in molecule mol_index
"""
import numpy as np
from copy import deepcopy
from ase import Atom, Atoms
from collections import Counter
from networkx import relabel_nodes

from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods                                          import get_bond_lengths
from SUMELF.SUMELF.general_methods.general_molecules_methods                                                      import get_number_of_lone_pairs_of_electron_pairs
from SUMELF.SUMELF.add_atoms.Utilities.get_bonding_unit_vectors                                                   import get_bonding_unit_vectors
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors                                               import get_new_bonding_unit_vectors
from SUMELF.SUMELF.add_atoms.add_hydrogens_to_molecules_methods.add_hydrogens_to_molecule                         import add_hydrogens_to_molecule

from SUMELF.SUMELF.add_atoms.Utilities.free_rotation_methods.allow_hydrogens_to_freely_rotate                     import allow_hydrogens_to_freely_rotate

carbon_bond_lengths   = get_bond_lengths('C')
hydrogen_bond_lengths = get_bond_lengths('H')
def add_methyls_to_molecule(no_of_methyl_to_add_to_mol, molecule, molecule_graph, crystal, crystal_graph, molecule_name, add_hydrogens=True, logger=None):
	"""
	This program will allow the user to manually add methyl groups to atoms in your crystals.

	Parameters
	----------
	no_of_methyl_to_add_to_mol : dict of {int: int}
		This dictionary indicates which atoms you want to add methyl groups to, and how many methyl groups you want to add to the atom. 
	molecule : ase.Atoms
		This is the molecule you want to add new methyl groups to
	molecule_graph : networkx.Graph
		This is the graph of this molecule.
	crystal : ase.Atoms
		This is the crystal as an ASE Atoms object
	crystal_graph : networkx.graph
		This is the graph associated with the crystal.
	molecule_name : int
		This is the name of the molecule that methyl groups are being added to.
	add_hydrogens : bool.
		this indicates if you want to add hydrogens in this method, or if you will do this later on.
	logger : logging
		This is the logger that info and warning message are sent to during the program

	Returns
	-------
	copied_molecule : ase.Atoms
		This is the molecule with aliphatic carbons removed
	copied_molecule_graph : networkx.Graph
		This is the modified graph of this molecule.
	were_methyls_added_to_this_molecule : bool.
		This boolean indicates if any methyl groups were added to this molecule.
	"""

	# First get a copy of the molecule.
	copied_molecule = molecule.copy()

	# Second, get a copy of the molecule graph.
	copied_molecule_graph = deepcopy(molecule_graph)

	# Third, initialise a boolean to determine if one or more methyl groups were added to this molecule.
	were_methyls_added_to_this_molecule = False

	# Fourth, set the element to attach as to the atom as C (the central atom of the methyl group).
	element_to_attach = 'C'

	# Fifth, initialise a list for recording the indices connected/apart of the newly created methyl group
	hydrogens_to_add_to_atoms = {}

	# Fifth, for each atom you want to attach atoms to in your molecule
	for node_name, no_of_methyl_to_add_to_atom in sorted(no_of_methyl_to_add_to_mol.items()):

		# 5.1: These are the properties of the atom/node you want to add methyl groups to.
		node_properties = copied_molecule_graph.nodes[node_name]

		# 5.2: Obtain the hybridisation for this atom
		hybridisation = node_properties['hybridisation']

		# 5.3: Obtain the formal charge value for this atom
		formal_charge = int(round(molecule[node_name].charge))

		# 5.4: Get the unit vectors for all bonds surrounding node_name, including attached methyls (with positions).
		all_bonding_unit_vectors = get_bonding_unit_vectors(copied_molecule, copied_molecule_graph, node_name)

		# 5.5: Obtain the total number of neighbours about node_name after adding methyl groups.
		total_number_of_neighbours_after_methyls_added = len(all_bonding_unit_vectors) + no_of_methyl_to_add_to_atom

		# 5.6: Check that the molecule does not contain more than 4 neighbours. We are only deal with molecules that obey the octet rule. 
		if total_number_of_neighbours_after_methyls_added > 4:
			from ase.visualize import view
			view(molecule)
			print('Can not add methyl groups to an atoms that will have more than 4 atoms around it in total. index: '+str(node_name))
			print(node_name, total_number_of_neighbours_after_methyls_added, hybridisation, formal_charge)
			import pdb; pdb.set_trace()
			exit()

		# 5.7: Obtain all the bonding unit vectors for atoms that already exist abount atom node_name in the molecule. 
		bonding_unit_vectors = [position for index, position in all_bonding_unit_vectors]

		# 5.8: Obtain the number of long pairs of electrons for this atom.
		number_of_lone_pairs_of_electrons = get_number_of_lone_pairs_of_electron_pairs(hybridisation, total_number_of_neighbours_after_methyls_added, formal_charge, logger=logger)

		# 5.9: Get the total number of new bonding unit vectors to get. These vectors will be used to:
		#      * Update the position of existing methyl groups attached to atom node_name
		#      * Add new methyl groups to atom node_name
		total_number_of_new_bonding_unit_vectors_to_get = no_of_methyl_to_add_to_atom + 0

		# 5.10: Obtain the bonding unit vectors to atteach methyl groups to the atom.
		new_bonding_unit_vectors = get_new_bonding_unit_vectors(total_number_of_new_bonding_unit_vectors_to_get, bonding_unit_vectors, copied_molecule, copied_molecule_graph, node_name, number_of_lone_pairs_of_electrons, crystal, crystal_graph, element_to_attach=element_to_attach, logger=logger)

		# 5.11: Check that there are not more methyl group C indices in indices_of_neighbouring_methyls_to_update_positions_of than the number of new bonding vectors in new_bonding_unit_vectors.
		#       * If this is the case, there is a programming problem.
		if not (total_number_of_new_bonding_unit_vectors_to_get == len(new_bonding_unit_vectors)):
			to_string  = 'Error: There are unequal numbers of methyls to add/update and new bonding unit vectors.\n'
			to_string += 'len(indices_of_neighbouring_methyls_to_update_positions_of) = '+str(0)+'\n'
			to_string += 'len(new_bonding_unit_vectors) = '+str(len(new_bonding_unit_vectors))+'\n'
			to_string += 'This indicates there is a programming error. Check this out.'
			print(to_string)
			import pdb; pdb.set_trace()
			raise Exception(to_string)

		# --------------------------------------------------------------------------------------------------------

		# 5.12: Obtain the element of the atom to attach the methyl group(s) to.
		origin_element  = copied_molecule[node_name].symbol

		# 5.13: Obtain the position of the atom to attach the methyl group(s) to.
		origin_position = copied_molecule[node_name].position

		# 5.14: Initialise a dictionary for recording newly added methyl group. 
		#       * Here, the index of the newly added C in the methyl group has been added. 
		#added_methyls_C_indices = []

		# 5.15: Obtain the positions for planing atoms in.
		for new_bonding_unit_vector in new_bonding_unit_vectors:

			# 5.16.1: record the index of the newly added C atom for the methyl group.
			were_methyls_added_to_this_molecule = True

			# 5.16.2: Determine the average bond length for a bond between a C atom (methyl group) and this atoms element, 
			#         where the atom contain total_no_of_neighbouring_atoms number of neighbouring atoms bonded to it.
			element_to_C_bond_length = carbon_bond_lengths[origin_element][total_number_of_neighbours_after_methyls_added]

			# 5.16.3: Determine the position to place the new C atom (methyl group) in the molecule, that is bonded to atom with index: node_name
			new_atom_position = origin_position + (element_to_C_bond_length * new_bonding_unit_vector)

			# ADD NEW C ATOM (METHYL GROUP) TO MOLECULE #

			# 5.16.4: Add this new C atom to the molecule.
			copied_molecule.append(Atom(element_to_attach,position=new_atom_position))

			# 5.16.5: Get the index of this newly appended C atom.
			new_carbon_atom_index = len(copied_molecule) - 1

			# 5.16.6: Check that new_carbon_atom_index does not already exist in copied_molecule_graph.
			#         * If it does, there is a programming error somewhere.
			if new_carbon_atom_index in copied_molecule_graph.nodes:
				raise Exception('Error: index '+str(new_carbon_atom_index)+' already exists in the graph for this molecule. Check this, as this indicates there is a programming error.')

			# 5.16.7: Add this index as a node to the molecule's graph
			carbon_node_information = {'E': element_to_attach, 'is_H_acceptor': None, 'involved_in_no_of_rings': False, 'is_H_donor': None, 'is_spiro_atom': False, 'hybridisation': 'sp3', 'added_or_modified': True}
			copied_molecule_graph.add_node(new_carbon_atom_index, **carbon_node_information)

			# 5.16.8: Add this bond to the copied_molecule_graph
			carbon_edge_information = {'bond_type': 'Single', 'is_conjugated': False, 'is_cyclic': False, 'involved_in_no_of_rings': None, 'bond_type_from_sybyl_type': 1}
			copied_molecule_graph.add_edge(node_name, new_carbon_atom_index, **carbon_edge_information)

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			# 5.16.9: Add hydrogens to the methyl group.
			if add_hydrogens:

				# 5.16.9.1: Create a dictionary to indicate we want to add hydrogens to 
				no_of_hydrogen_to_add_to_mol = {new_carbon_atom_index: 3}

				# 5.16.9.2: Get the number of elements in the molecule before adding hydrogens
				element_counter_before_H_addition = Counter(copied_molecule.get_chemical_symbols())

				# 5.16.9.3: Add hydrogens to this newly added carbon atom to form a methyl group in the molecule.
				copied_molecule, copied_molecule_graph, were_hydrogens_added_to_this_molecule, added_H_indices = add_hydrogens_to_molecule(no_of_hydrogen_to_add_to_mol, copied_molecule, copied_molecule_graph, crystal, crystal_graph, molecule_name=molecule_name, allow_hydrogen_free_rotation=False, logger=None)

				# 5.16.9.4: Make sure that the were_hydrogens_added_to_this_molecule tag was True
				#          * If this is False, it indicates that hydrogens were not added, which doesn't m ake sense and indicates a programming issue. 
				if not were_hydrogens_added_to_this_molecule:
					to_string  = f'Error: Hydrogens were not added to atom index {new_carbon_atom_index} in molecule {node_name}.\n'
					to_string += 'Check this.'
					raise Exception(to_string)

				# 5.16.9.5: Obtain the change in elements before and after adding hydrogens
				element_counter_after_H_addition = Counter(copied_molecule.get_chemical_symbols())
				element_counter_diff = {symbol: count for symbol, count in (element_counter_after_H_addition - element_counter_before_H_addition).items() if (count != 0)}

				# 5.16.9.6: Make suer that only three hydrogens have been added to the updated system
				if not element_counter_diff == {'H': 3}:
					to_string  = f'Error: The correct number of elements afetererfedfcedvc.......\n'
					to_string += 'Check this.'
					import pdb; pdb.set_trace()
					raise Exception(to_string)

				# 5.16.9.7: Add the hydrogen atoms connected to this newly added ethyl group to added_methyls_H_indices. 
				#added_methyls_H_indices += added_H_indices # list(range(len(copied_molecule)-3, len(copied_molecule)))

			else:

				# 5.16.9.8: Check that new_carbon_atom_index is not already in hydrogens_to_add_to_atoms
				if new_carbon_atom_index in hydrogens_to_add_to_atoms.keys():
					raise Exception('new_carbon_atom_index already in hydrogens_to_add_to_atoms')

				# 5.16.9.8: Add the index of this newly added carbon atom to added_methyls_C_indices
				hydrogens_to_add_to_atoms[new_carbon_atom_index] = 3

		# --------------------------------------------------------------------------------------------------------

	debugging_molecule = copied_molecule.copy()

	# Sixth, allow newly added and modified hydrogens to freely rotate. 
	#if were_hydrogens_added_to_this_molecule:
	#	copied_molecule = allow_hydrogens_to_freely_rotate(copied_molecule, copied_molecule_graph, added_methyls_H_indices, molecule_name=molecule_name, method_type='add_methyls_to_molecule')

	'''
	from ase.visualize import view
	view([debugging_molecule, copied_molecule])
	import pdb; pdb.set_trace()
	raise Exception('Check this method is working.')
	'''

	# Seventh, return the copied_molecule and copied_molecule_graph with the newly added methyl groups, as well as were_methyls_added_to_this_molecule and hydrogens_to_add_to_atoms
	return copied_molecule, copied_molecule_graph, were_methyls_added_to_this_molecule, hydrogens_to_add_to_atoms



