"""
add_hydrogens_to_molecule.py, Geoffrey Weal, 25/6/2022

This program will allow the user to manually add hydrogens to atoms in molecule mol_index
"""
import numpy as np
from copy import deepcopy
from ase import Atom, Atoms
from networkx import relabel_nodes

from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods                      import get_bond_lengths
from SUMELF.SUMELF.general_methods.general_molecules_methods                                  import get_number_of_lone_pairs_of_electron_pairs
#from SUMELF.SUMELF.add_atoms.add_hydrogens_to_molecules_methods.check_for_hydrogen_bonding    import check_for_hydrogen_bonding
#from SUMELF.SUMELF.add_atoms.add_hydrogens_to_molecules_methods.check_for_H_pi_bonding        import check_for_H_pi_bonding
from SUMELF.SUMELF.add_atoms.Utilities.get_bonding_unit_vectors                               import get_bonding_unit_vectors
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors                           import get_new_bonding_unit_vectors

bond_lengths = get_bond_lengths('H')
def add_hydrogens_to_molecule(no_of_hydrogen_to_add_to_mol, molecule, molecule_graph, crystal, crystal_graph, molecule_name=None, allow_hydrogen_free_rotation=None, logger=None):
	"""
	This program will allow the user to manually add hydrogens to atoms in your crystals.

	Parameters
	----------
	no_of_hydrogen_to_add_to_mol : dict of {int: int}
		This dictionary indicates which atoms you want to add hydrogens to and how many hydrogen you want to add to the atom. 
	molecule : ase.Atoms
		This is the molecule you want to add missing hydrogens to
	molecule_graph : networkx.Graph
		This is the graph of this molecule.
	crystal : ase.Atoms
		This is the crystal as an ASE Atoms object
	crystal_graph : networkx.graph
		This is the graph associated with the crystal.
	molecule_name : int
		This is the name of the molecule that ethyl groups are being added to.
	allow_hydrogen_free_rotation : bool.
		If True, this method will allow the newly added hydrogens to rotate about the atom they are bound to. 
	logger : logging
		This is the logger that info and warning message are sent to during the program

	Returns
	-------
	copied_molecule : ase.Atoms
		This is the molecule with aliphatic carbons removed
	copied_molecule_graph : networkx.Graph
		This is the modified graph of this molecule.
	were_hydrogens_added_to_this_molecule : bool.
		This boolean indicates if any hydrogens were added to this molecule.
	"""

	# First get a copy of the molecule.
	copied_molecule = molecule.copy()

	# Second, get a copy of the molecule graph.
	copied_molecule_graph = deepcopy(molecule_graph)

	# Third, initialise a boolean to determine if one or more hydrogens were added to this molecule.
	were_hydrogens_added_to_this_molecule = False

	# Fourth, set the element to attach as hydrogen
	element_to_attach = 'H'

	# Fifth, initialise a list of the indices of all the newly added hydrogens. 
	atoms_with_added_hydrogens_indices = {}

	# Sixth, for each atom you want to attach atoms to in your molecule
	for node_index, no_of_hydrogen_to_add_to_atom in sorted(no_of_hydrogen_to_add_to_mol.items(), key=lambda x: (x[1], x[0])):

		# 6.1: These are the properties of the atom/node you want to add hydrogens to.
		node_properties = copied_molecule_graph.nodes[node_index]

		# 6.2: Obtain the hybridisation for this atom
		hybridisation = node_properties['hybridisation']

		# 6.3: Obtain the formal charge value for this atom
		formal_charge = int(round(molecule[node_index].charge))

		# 6.4: Get the unit vectors for all bonds surrounding node_index, including attached hydrogens (with positions).
		all_bonding_unit_vectors = get_bonding_unit_vectors(copied_molecule, copied_molecule_graph, node_index)

		# 6.5: Obtain the total number of neighbours about node_index after adding hydrogens.
		total_number_of_neighbours_after_hydrogens_added = len(all_bonding_unit_vectors) + no_of_hydrogen_to_add_to_atom

		# 6.6: Check that the molecule does not contain more than 4 neighbours. We are only deal with molecules that obey the octet rule. 
		if total_number_of_neighbours_after_hydrogens_added > 4:
			from ase.visualize import view
			view(molecule)
			print('Can not add hydrogens to an atoms that will have more than 4 atoms around it in total. index: '+str(node_index))
			print(node_index, total_number_of_neighbours_after_hydrogens_added, hybridisation, formal_charge)
			import pdb; pdb.set_trace()
			exit()

		# 6.7: Only include hydrogens the unit vectors for non-hydrogen bonds. 
		indices_of_neighbouring_hydrogens_to_update_positions_of = []
		bonding_unit_vectors = []
		for index, position in all_bonding_unit_vectors:
			if molecule[index].symbol == element_to_attach:
				indices_of_neighbouring_hydrogens_to_update_positions_of.append(index)
			else:
				bonding_unit_vectors.append(position)

		# 6.8: Obtain the number of long pairs of electrons for this atom.
		number_of_lone_pairs_of_electrons = get_number_of_lone_pairs_of_electron_pairs(hybridisation, total_number_of_neighbours_after_hydrogens_added, formal_charge, logger=logger)

		# 6.9: Get the total number of new bonding unit vectors to get. These vectors will be used to:
		#      * Update the position of existing hydrogens attached to atom node_index
		#      * Add new hydrogens to atom node_index
		total_number_of_new_bonding_unit_vectors_to_get = no_of_hydrogen_to_add_to_atom + len(indices_of_neighbouring_hydrogens_to_update_positions_of)

		# 6.10: Obtain the bonding unit vectors to attach hydrogens to the atom.
		new_bonding_unit_vectors = get_new_bonding_unit_vectors(total_number_of_new_bonding_unit_vectors_to_get, bonding_unit_vectors, copied_molecule, copied_molecule_graph, node_index, number_of_lone_pairs_of_electrons, crystal, crystal_graph, element_to_attach=element_to_attach, logger=logger)

		# 6.11: Check that there are not more hydrogen indices in indices_of_neighbouring_hydrogens_to_update_positions_of than the number of new bonding vectors in new_bonding_unit_vectors.
		#       * If this is the case, there is a programming problem.
		if not (total_number_of_new_bonding_unit_vectors_to_get == len(new_bonding_unit_vectors)):
			toString  = 'Error: There are unequal numbers of hydrogens to add/update and new bonding unit vectors.\n'
			toString += 'len(indices_of_neighbouring_hydrogens_to_update_positions_of) = '+str(len(indices_of_neighbouring_hydrogens_to_update_positions_of))+'\n'
			toString += 'len(new_bonding_unit_vectors) = '+str(len(new_bonding_unit_vectors))+'\n'
			toString += 'This indicates there is a programming error. Check this out.'
			print(toString)
			import pdb; pdb.set_trace()
			raise Exception(toString)

		# --------------------------------------------------------------------------------------------------------

		# 6.12: Obtain the element of the atom to attach the hydrogen(s) to.
		origin_element  = copied_molecule[node_index].symbol

		# 6.13: Obtain the position of the atom to attach the hydrogen(s) to.
		origin_position = copied_molecule[node_index].position

		# 6.14: Initialise a dictionary for recording newly added or updated hydrogens. 
		modified_and_added_hydrogens = {}

		# 6.15: Obtain the positions for planing atoms in.
		for new_bonding_unit_vector in new_bonding_unit_vectors:

			# 6.15.1: Record the index of the newly added hydrogen
			were_hydrogens_added_to_this_molecule = True

			# 6.15.2: Determine the average bond length for a bond between a hydrogen and this atoms element, 
			#         where the atom contain total_no_of_neighbouring_atoms number of neighbouring atoms bonded to it.
			element_to_H_bond_length = bond_lengths[origin_element][total_number_of_neighbours_after_hydrogens_added]

			# 6.15.3: Determine the position to place the new hydrogen atom in the molecule, that is bonded to atom with index: node_index
			new_atom_position = origin_position + (element_to_H_bond_length * new_bonding_unit_vector)

			# 6.15.4: If their exist hydrogen indices in indices_of_neighbouring_hydrogens_to_update_positions_of, update its position with new_atom_position 
			if len(indices_of_neighbouring_hydrogens_to_update_positions_of) > 0:

				# 6.15.4.1: UPDATE AN EXISTING HYDROGEN

				# 6.15.4.1.1: Pop a index from indices_of_neighbouring_hydrogens_to_update_positions_of
				hydrogen_to_update_index = indices_of_neighbouring_hydrogens_to_update_positions_of.pop()

				# 6.15.4.1.2: Check that the element of the atom in th emolecule to update is the same as element_to_attach
				#             * If this is not the case, there is a programming error.
				if not (copied_molecule[hydrogen_to_update_index].symbol == element_to_attach):
					toString  = 'Error: The atom you are updating (index '+str(hydrogen_to_update_index)+') does not have the expected element.\n'
					toString += 'copied_molecule['+str(hydrogen_to_update_index)+'].symbol = '+str(copied_molecule[hydrogen_to_update_index].symbol)+'\n'
					toString += 'element_to_attach = '+str(element_to_attach)+'\n'
					toString += 'These element enteries should be the same. There is a programming error. Check this.'
					raise Exception(toString)

				# 6.15.4.1.3: Update the position of atom hydrogen_to_update_index to new_atom_position
				copied_molecule[hydrogen_to_update_index].position = new_atom_position

				# 6.15.4.1.4: Update that this hydrogen atom has been modified to the molecule's graph
				copied_molecule_graph.nodes[hydrogen_to_update_index]['added_or_modified'] = True

				# 6.15.4.1.5: Record the index of the updated hydrogen.
				modified_and_added_hydrogens[hydrogen_to_update_index] = 'updated'

			else:

				# 6.15.4.2: ADD NEW HYDROGEN TO MOLECULE

				# 6.15.4.2.1: Add this new hydrogen atom to the molecule.
				copied_molecule.append(Atom(element_to_attach,position=new_atom_position))

				# 6.15.4.2.2: Get the index of this newly appended hydrogen atom.
				new_hydrogen_atom_index = len(copied_molecule) - 1

				# 6.15.4.2.3: Check that new_hydrogen_atom_index does not already exist in copied_molecule_graph.
				#             * If it does, there is a programming error somewhere.
				if new_hydrogen_atom_index in copied_molecule_graph.nodes:
					raise Exception('Error: index '+str(new_hydrogen_atom_index)+' already exists in the graph for this molecule. Check this, as this indicates there is a programming error.')

				# 6.15.4.2.4: Add this index as a node to the molecule's graph
				hydrogen_node_information = {'E': element_to_attach, 'is_H_acceptor': None, 'involved_in_no_of_rings': False, 'is_H_donor': None, 'is_spiro_atom': False, 'hybridisation': '-', 'added_or_modified': True}
				copied_molecule_graph.add_node(new_hydrogen_atom_index, **hydrogen_node_information)

				# 6.15.4.2.5: Add this bond to the copied_molecule_graph
				hydrogen_edge_information = {'bond_type': 'Single', 'is_conjugated': False, 'is_cyclic': False, 'involved_in_no_of_rings': None, 'bond_type_from_sybyl_type': 1}
				copied_molecule_graph.add_edge(node_index, new_hydrogen_atom_index, **hydrogen_edge_information)

				# 6.15.4.1.5: Record the index of the newly added hydrogen.
				modified_and_added_hydrogens[new_hydrogen_atom_index] = 'new'

		# 6.16: Check that there are no indices in indices_of_neighbouring_hydrogens_to_update_positions_of
		#       All of the atoms in indices_of_neighbouring_hydrogens_to_update_positions_of should have have been updated
		#       and therefore remove from this list. If this is not the case, there is a programming error.
		#       Note: This is a double check, as their should have been caught by the error check at 5.10
		if len(indices_of_neighbouring_hydrogens_to_update_positions_of) > 0:
			toString  = 'Error: There are still indices in the indices_of_neighbouring_hydrogens_to_update_positions_of list.\n'
			toString += 'indices_of_neighbouring_hydrogens_to_update_positions_of = '+str(indices_of_neighbouring_hydrogens_to_update_positions_of)+'\n'
			toString += 'This indicates a programming error. Check this.'
			raise Exception(toString)

		# 6.17: Check that node_index is not already in atoms_with_added_hydrogens_indices
		if node_index in atoms_with_added_hydrogens_indices.keys():
			raise Exception('Error: node_index already in atoms_with_added_hydrogens_indices')

		# 6.17: Get the indices of the hydrogens atoms added to node_index
		atoms_with_added_hydrogens_indices[node_index] = list(modified_and_added_hydrogens.keys())

	# --------------------------------------------------------------------------------------------------------

	# Seventh, return the copied_molecule and copied_molecule_graph with the newly added hydrogen atoms, and were_hydrogens_added_to_this_molecule
	return copied_molecule, copied_molecule_graph, were_hydrogens_added_to_this_molecule, atoms_with_added_hydrogens_indices





