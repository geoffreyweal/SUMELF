"""
add_ethyls_to_molecule.py, Geoffrey Weal, 15/4/2024

This program will allow the user to manually add ethyl groupss to atoms in molecule mol_index
"""
from math import pi as PI
import numpy as np
from copy import deepcopy
from ase import Atom, Atoms
from collections import Counter
from networkx import relabel_nodes

from SUMELF.SUMELF.general_methods.geometry_methods                                       import get_unit_vector, project_u_onto_v, rotate_vector_around_axis
from SUMELF.SUMELF.general_methods.distance_methods                                       import get_distance
from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods                  import get_bond_lengths
from SUMELF.SUMELF.general_methods.general_molecules_methods                              import get_number_of_lone_pairs_of_electron_pairs
from SUMELF.SUMELF.add_atoms.Utilities.get_bonding_unit_vectors                           import get_bonding_unit_vectors
from SUMELF.SUMELF.add_atoms.Utilities.get_new_bonding_unit_vectors                       import get_new_bonding_unit_vectors
from SUMELF.SUMELF.add_atoms.add_hydrogens_to_molecules_methods.add_hydrogens_to_molecule import add_hydrogens_to_molecule

carbon_bond_lengths   = get_bond_lengths('C')
hydrogen_bond_lengths = get_bond_lengths('H')
def add_ethyls_to_molecule(no_of_ethyl_to_add_to_mol, molecule, molecule_graph, crystal, crystal_graph, molecule_name, add_hydrogens=True, logger=None):
	"""
	This program will allow the user to manually add ethyl groups to atoms in your crystals.

	Parameters
	----------
	no_of_ethyl_to_add_to_mol : dict of {int: int}
		This dictionary indicates which atoms you want to add ethyl groups to, and how many ethyl groups you want to add to the atom. 
	molecule : ase.Atoms
		This is the molecule you want to add new ethyl groups to
	molecule_graph : networkx.Graph
		This is the graph of this molecule.
	crystal : ase.Atoms
		This is the crystal as an ASE Atoms object
	crystal_graph : networkx.graph
		This is the graph associated with the crystal.
	molecule_name : int
		This is the name of the molecule that ethyl groups are being added to.
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
	were_ethyls_added_to_this_molecule : bool.
		This boolean indicates if any ethyl groups were added to this molecule.
	"""

	# First get a copy of the molecule.
	copied_molecule = molecule.copy()

	# Second, get a copy of the molecule graph.
	copied_molecule_graph = deepcopy(molecule_graph)

	# Third, initialise a boolean to determine if one or more ethyl groups were added to this molecule.
	were_ethyls_added_to_this_molecule = False

	# Fourth, set the element to attach as to the atom as C (the central atom of the ethyl group).
	element_to_attach = 'C'

	# Fifth, initialise a list for recording the which hydrogens were added to which atoms in the molecule.
	hydrogens_to_add_to_atoms = {}

	# Fifth, for each atom you want to attach atoms to in your molecule
	for node_name, no_of_ethyl_to_add_to_atom in sorted(no_of_ethyl_to_add_to_mol.items()):

		# 5.1: These are the properties of the atom/node you want to add ethyl groups to.
		node_properties = copied_molecule_graph.nodes[node_name]

		# 5.2: Obtain the hybridisation for this atom
		hybridisation = node_properties['hybridisation']

		# 5.3: Obtain the formal charge value for this atom
		formal_charge = int(round(molecule[node_name].charge))

		# 5.4: Get the unit vectors for all bonds surrounding node_name, including attached ethyls (with positions).
		all_bonding_unit_vectors = get_bonding_unit_vectors(copied_molecule, copied_molecule_graph, node_name)

		# 5.5: Obtain the total number of neighbours about node_name after adding ethyl groups.
		total_number_of_neighbours_after_ethyls_added = len(all_bonding_unit_vectors) + no_of_ethyl_to_add_to_atom

		# 5.6: Check that the molecule does not contain more than 4 neighbours. We are only deal with molecules that obey the octet rule. 
		if total_number_of_neighbours_after_ethyls_added > 4:
			from ase.visualize import view
			view(molecule)
			print('Can not add ethyl groups to an atoms that will have more than 4 atoms around it in total. index: '+str(node_name))
			print(node_name, total_number_of_neighbours_after_ethyls_added, hybridisation, formal_charge)
			import pdb; pdb.set_trace()
			exit()

		# 5.7: Obtain all the bonding unit vectors for atoms that already exist abount atom node_name in the molecule. 
		bonding_unit_vectors = [position for index, position in all_bonding_unit_vectors]

		# 5.8: Obtain the number of long pairs of electrons for this atom.
		number_of_lone_pairs_of_electrons = get_number_of_lone_pairs_of_electron_pairs(hybridisation, total_number_of_neighbours_after_ethyls_added, formal_charge, logger=logger)

		# 5.9: Get the total number of new bonding unit vectors to get. These vectors will be used to:
		#      * Update the position of existing ethyl groups attached to atom node_name
		#      * Add new ethyl groups to atom node_name
		total_number_of_new_bonding_unit_vectors_to_get = no_of_ethyl_to_add_to_atom + 0

		# 5.10: Obtain the bonding unit vectors to atteach ethyl groups to the atom.
		new_bonding_unit_vectors = get_new_bonding_unit_vectors(total_number_of_new_bonding_unit_vectors_to_get, bonding_unit_vectors, copied_molecule, copied_molecule_graph, node_name, number_of_lone_pairs_of_electrons, crystal, crystal_graph, element_to_attach=element_to_attach, logger=logger)

		# 5.11: Check that there are not more ethyl group C indices in indices_of_neighbouring_ethyls_to_update_positions_of than the number of new bonding vectors in new_bonding_unit_vectors.
		#       * If this is the case, there is a programming problem.
		if not (total_number_of_new_bonding_unit_vectors_to_get == len(new_bonding_unit_vectors)):
			to_string  = 'Error: There are unequal numbers of ethyls to add/update and new bonding unit vectors.\n'
			to_string += 'len(indices_of_neighbouring_ethyls_to_update_positions_of) = '+str(0)+'\n'
			to_string += 'len(new_bonding_unit_vectors) = '+str(len(new_bonding_unit_vectors))+'\n'
			to_string += 'This indicates there is a programming error. Check this out.'
			print(to_string)
			import pdb; pdb.set_trace()
			raise Exception(to_string)

		# --------------------------------------------------------------------------------------------------------

		# 5.12: Obtain the element of the atom to attach the ethyl group(s) to.
		origin_element  = copied_molecule[node_name].symbol

		# 5.13: Obtain the position of the atom to attach the ethyl group(s) to.
		origin_position = copied_molecule[node_name].position

		# 5.14: Obtain the positions for adding atoms in.
		for new_bonding_unit_vector in new_bonding_unit_vectors:

			# 5.14.1: record the index of the newly added C atom for the ethyl group.
			were_ethyls_added_to_this_molecule = True

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			# GET POSITION OF NEW alpha-C ATOM

			# 5.14.2: Determine the average bond length for a bond between a C atom (ethyl group) and this atoms element, 
			#         where the atom contain total_no_of_neighbouring_atoms number of neighbouring atoms bonded to it.
			element_to_C_bond_length = carbon_bond_lengths[origin_element][total_number_of_neighbours_after_ethyls_added]

			# 5.14.3: Determine the position to place the new C atom (ethyl group) in the molecule, that is bonded to atom with index: node_name
			alpha_carbon_atom_position = origin_position + (element_to_C_bond_length * new_bonding_unit_vector)

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			# ADD NEW alpha-C ATOM (ETHYL GROUP) TO MOLECULE 

			# 5.14.4: Add this new C atom to the molecule.
			copied_molecule.append(Atom(element_to_attach,position=alpha_carbon_atom_position))

			# 5.14.5: Get the index of this newly appended C atom.
			alpha_carbon_atom_index = len(copied_molecule) - 1

			# 5.14.6: Check that alpha_carbon_atom_index does not already exist in copied_molecule_graph.
			#         * If it does, there is a programming error somewhere.
			if alpha_carbon_atom_index in copied_molecule_graph.nodes:
				raise Exception('Error: index '+str(alpha_carbon_atom_index)+' already exists in the graph for this molecule. Check this, as this indicates there is a programming error.')

			# 5.14.7: Add this index as a node to the molecule's graph
			carbon_node_information = {'E': element_to_attach, 'is_H_acceptor': None, 'involved_in_no_of_rings': False, 'is_H_donor': None, 'is_spiro_atom': False, 'hybridisation': 'sp3', 'added_or_modified': True}
			copied_molecule_graph.add_node(alpha_carbon_atom_index, **carbon_node_information)

			# 5.14.8: Add this bond to the copied_molecule_graph
			carbon_edge_information = {'bond_type': 'Single', 'is_conjugated': False, 'is_cyclic': False, 'involved_in_no_of_rings': None, 'bond_type_from_sybyl_type': 1}
			copied_molecule_graph.add_edge(node_name, alpha_carbon_atom_index, **carbon_edge_information)

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			# GET POSITION OF NEW beta-C ATOM

			# 5.14.10: Get the bond length between two sp3 carbons. 
			CC_bond_length = carbon_bond_lengths['C'][3]

			# 5.14.11: Get any unit vector that is in the plane with normal given as new_bonding_unit_vector
			perp_unit_vector = get_unit_vector(vector_in_plane(new_bonding_unit_vector))

			# 5.14.12: Rotate new_bonding_unit_vector about perp_beta_bond_unit_vector by 109.5 degrees to given a beta unit bonding vector for 
			starting_beta_bond_unit_vector = get_unit_vector(rotate_vector_around_axis(-new_bonding_unit_vector, np.radians(109.5), perp_unit_vector))

			# 5.14.13: Initalise a dictionary to record the coulomb force as starting_beta_bond_unit_vector is rotated about new_bonding_unit_vector.
			force_vs_angle = {}

			# 5.14.14: Record the coulomb force upon the beta-C as it is rotated about the new_bonding_unit_vector bond
			for angle_degrees in range(3600):

				# 5.14.14.1: Divide by 10 
				angle_degrees = float(angle_degrees)/10.0

				# 5.14.14.2: Get the bonding vector after rotating starting_beta_bond_unit_vector about new_bonding_unit_vector by angle_degrees
				rotated_starting_beta_bond_unit_vector = get_unit_vector(rotate_vector_around_axis(deepcopy(starting_beta_bond_unit_vector), np.radians(angle_degrees), new_bonding_unit_vector))

				# 5.14.14.3: Get the position of the beta carbon after this rotation. 
				beta_carbon_atom_position = alpha_carbon_atom_position + (CC_bond_length * rotated_starting_beta_bond_unit_vector)

				# 5.14.14.4: Obtain the unit vectors and the distances from the beta-C to each atom in copied_molecule
				unit_vectors_from_atoms_to_beta_C = [get_unit_vector(beta_carbon_atom_position-position) for position in copied_molecule.get_positions()]
				distances_from_atoms_to_beta_C    = [get_distance(beta_carbon_atom_position,position)    for position in copied_molecule.get_positions()]

				# 5.14.14.5: Get the force vectors upon the beta-C atom from the atoms in copied_molecule 
				force_vectors = [(1/(distance_from_atom_to_beta_C ** 2.0))*unit_vector_from_atom_to_beta_C for distance_from_atom_to_beta_C, unit_vector_from_atom_to_beta_C in zip(distances_from_atoms_to_beta_C, unit_vectors_from_atoms_to_beta_C)]

				# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
				# This is an old method for obtaining the total magnitude based on sum of vectors.

				# Get the overall force vector upon the beta-C atom.
				#force_vector = sum(force_vectors)

				# Get the magnitude of the force.
				#force_magnitude = np.linalg.norm(force_vector)

				# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
				# This is the new method for obtaining the total magnitude based on sum of vectors magnitudes.

				# 5.14.14.6: Obtain the magnitude of each force 
				force_magnitudes = [np.linalg.norm(force_vector) for force_vector in force_vectors]

				# 5.14.14.7: Obtain the total sum of the force magnitude. 
				force_magnitude = sum(force_magnitudes)

				# 5.16.14.8: Record the force vector at each angle.
				force_vs_angle[angle_degrees] = force_magnitude

			# 5.16.15: Get the angle that has the force with the largest magnitude.
			max_amplitude_angle, max_amplitude = min(force_vs_angle.items(), key=lambda x:x[1])

			# 5.16.16: Get the vector that minimises the force upon the beta carbon and the other atoms in the molecule. 
			rotated_starting_beta_bond_unit_vector = get_unit_vector(rotate_vector_around_axis(starting_beta_bond_unit_vector, np.radians(max_amplitude_angle), new_bonding_unit_vector))

			# 5.16.17: Get the position of the beta carbon. 
			beta_carbon_atom_position = alpha_carbon_atom_position + (CC_bond_length * rotated_starting_beta_bond_unit_vector)

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			# ADD NEW beta-C ATOM (ETHYL GROUP) TO MOLECULE 

			# 5.16.18: Add this new C atom to the molecule.
			copied_molecule.append(Atom(element_to_attach,position=beta_carbon_atom_position))

			# 5.16.19: Get the index of this newly appended C atom.
			beta_carbon_atom_index = len(copied_molecule) - 1

			# 5.16.20: Check that beta_carbon_atom_index does not already exist in copied_molecule_graph.
			#         * If it does, there is a programming error somewhere.
			if beta_carbon_atom_index in copied_molecule_graph.nodes:
				raise Exception('Error: index '+str(beta_carbon_atom_index)+' already exists in the graph for this molecule. Check this, as this indicates there is a programming error.')

			# 5.16.21: Add this index as a node to the molecule's graph.
			carbon_node_information = {'E': element_to_attach, 'is_H_acceptor': None, 'involved_in_no_of_rings': False, 'is_H_donor': None, 'is_spiro_atom': False, 'hybridisation': 'sp3', 'added_or_modified': True}
			copied_molecule_graph.add_node(beta_carbon_atom_index, **carbon_node_information)

			# 5.16.22: Add this bond to the copied_molecule_graph.
			carbon_edge_information = {'bond_type': 'Single', 'is_conjugated': False, 'is_cyclic': False, 'involved_in_no_of_rings': None, 'bond_type_from_sybyl_type': 1}
			copied_molecule_graph.add_edge(alpha_carbon_atom_index, beta_carbon_atom_index, **carbon_edge_information)

			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			# 5.16.23: Add hydrogens to the methyl group.
			if add_hydrogens:

				# 5.16.23.1: Create a dictionary to indicate we want to add hydrogens to.
				no_of_hydrogen_to_add_to_mol = {alpha_carbon_atom_index: 2, beta_carbon_atom_index: 3}

				# 5.16.23.2: Get the number of elements in the molecule before adding hydrogens.
				element_counter_before_H_addition = Counter(copied_molecule.get_chemical_symbols())

				# 5.16.23.3: add hydrogens to this newly added carbon atom to form a ethyl group in the molecule.
				copied_molecule, copied_molecule_graph, were_hydrogens_added_to_this_molecule, atoms_with_added_hydrogens_indices = add_hydrogens_to_molecule(no_of_hydrogen_to_add_to_mol, copied_molecule, copied_molecule_graph, crystal, crystal_graph, molecule_name=molecule_name, allow_hydrogen_free_rotation=False, logger=None)

				# 5.16.23.4: Make sure that the were_hydrogens_added_to_this_molecule tag was True
				#          * If this is False, it indicates that hydrogens were not added, which doesn't m ake sense and indicates a programming issue. 
				if not were_hydrogens_added_to_this_molecule:
					to_string  = f'Error: Hydrogens were not added to atom index {new_carbon_atom_index} in molecule {node_name}.\n'
					to_string += 'Check this.'
					raise Exception(to_string)

				# 5.16.23.5: Obtain the change in elements before and after adding hydrogens.
				element_counter_after_H_addition = Counter(copied_molecule.get_chemical_symbols())
				element_counter_diff = {symbol: count for symbol, count in (element_counter_after_H_addition - element_counter_before_H_addition).items() if (count != 0)}

				# 5.16.23.6: Make sure that only three hydrogens have been added to the updated system.
				if not element_counter_diff == {'H': 5}:
					to_string  = f'Error: A different number of elements have been added to the molecule after attempting to add 5 hydrogens to the molecule.\n'
					to_string += f'Element change after adding hydrogens to the system: {element_counter_diff}.\n'
					to_string += 'What this should be: '+str({'H': 5})+'.\n'
					to_string += 'Check this.'
					import pdb; pdb.set_trace()
					raise Exception(to_string)

			else:

				# 5.16.23.7: Check that alpha_carbon_atom_index is not already in hydrogens_to_add_to_atoms.
				if alpha_carbon_atom_index in hydrogens_to_add_to_atoms.keys():
					raise Exception('atom_index in hydrogens_to_add_to_atoms')

				# 5.16.23.8: Add 2 hydrogen atoms to alpha_carbon_atom_index in hydrogens_to_add_to_atoms.
				hydrogens_to_add_to_atoms[alpha_carbon_atom_index] = 2

				# 5.16.23.9: Check that beta_carbon_atom_index is not already in hydrogens_to_add_to_atoms.
				if beta_carbon_atom_index in hydrogens_to_add_to_atoms.keys():
					raise Exception('atom_index in hydrogens_to_add_to_atoms')

				# 5.16.23.10: Add 3 hydrogen atoms to beta_carbon_atom_index in hydrogens_to_add_to_atoms.
				hydrogens_to_add_to_atoms[beta_carbon_atom_index] = 3

	# ------------------------------------------------------------------------------------------------------------

	# Sixth, return the copied_molecule and copied_molecule_graph with the newly added ethyl groups, and as well as were_ethyls_added_to_this_molecule, hydrogens_to_add_to_atoms
	return copied_molecule, copied_molecule_graph, were_ethyls_added_to_this_molecule, hydrogens_to_add_to_atoms

# ----------------------------------------------------------------------------------------------------------------

def vector_in_plane(normal):
	# Generate a random vector not parallel to the normal vector
	if np.allclose(normal, [0, 0, 0]):
		raise ValueError("The normal vector cannot be the zero vector.")

	# Choose a random vector
	v_vector = np.random.rand(3)  # Generate a random vector
	while np.dot(v_vector, normal) == 0:
		v_vector = np.random.rand(3)

	# Project v onto the plane by subtracting off the component parallel to n
	vector_in_plane = v_vector - (project_u_onto_v(v_vector, normal) * normal)

	# Return vector_in_plane
	return vector_in_plane

# ----------------------------------------------------------------------------------------------------------------

