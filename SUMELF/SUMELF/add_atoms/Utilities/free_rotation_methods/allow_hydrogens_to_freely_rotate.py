"""
allow_hydrogens_to_freely_rotate.py, Geoffrey Weal, 15/4/2024

This method is designed to allow newly added and update hydrogen atoms to freely rotate.
"""
import sys
from ase.optimize                      import FIRE, LBFGS
from ase.constraints                   import FixAtoms, FixBondLengths, FixInternals
from SUMELF.SUMELF.calculators.Coulomb import Coulomb

def allow_hydrogens_to_freely_rotate(molecule, molecule_graph, hydrogen_indices_to_allow_free_rotation, molecule_name=None, method_type=None):
	"""
	This method is designed to allow newly added and update hydrogen atoms to freely rotate.

	Parameters
	----------
	molecule : ase.Atoms
		This is the molecule you want to allow hydrogens to freely rotate.
	molecule_graph : networkx.Graph
		This is the graph of this molecule.
	hydrogen_indices_to_allow_free_rotation : list of ints
		These are the hydrogens you would like to allow to freely rotate.
	molecule_name : int
		This is the name of the molecule that methyl groups are being added to.
	method_type : str
		This is the method that is asking for atoms to freely rotate. 

	Returns
	-------
	updated_molecule : ase.Atoms
		This is the molecule of interest, where hydrogens in modified_and_added_hydrogens have been allowed to freely rotate. 
	"""

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# First, obtain all the central atoms to allow the hydrogens to freely rotate about.

	raise Exception('Broken')

	# 1.1: Initalise a dictionary to hold which atoms to rotate, and around which atoms to rotate.
	central_atoms = {}

	# 1.2: For each index in the hydrogen_indices_to_allow_free_rotation list.
	for hydrogen_index in hydrogen_indices_to_allow_free_rotation:

		# 1.2.1: Check that hydrogen_index is infact a hydrogen atom.
		if molecule[hydrogen_index].symbol not in ['H', 'D', 'T']:
			to_string  = f'Error: atom index {hydrogen_index} is not a hydrogen.\n'
			to_string += f'Element of atom {hydrogen_index}: {molecule[hydrogen_index].symbol}\n'
			to_string += 'Check this'
			raise Exception(to_string)

		# 1.2.2: Get the indices of the neighbours about this hydrogen atom.
		neighbours_to_hydrogen_indices = tuple(molecule_graph[hydrogen_index].keys())

		# 1.2.3: Check that there are neighbours about this hydrogen, otherwise there is an issue. 
		if len(neighbours_to_hydrogen_indices) == 0:
			to_string  = 'Error: This hydrogen does not neighbour any atoms.\n'
			to_string += 'This means that this hydrogen can not rotate around an atom\n'
			to_string += f'Problematic hydrogen index: {hydrogen_index}\n'
			to_string += 'Check this'
			raise Exception(to_string)

		# 1.2.4: Make that this hydrogen only has one neighbour. 
		if len(neighbours_to_hydrogen_indices) >= 2:
			to_string  = 'Error: This hydrogen contain more than one neighbour.\n'
			to_string += f'Problematic hydrogen index: {hydrogen_index}\n'
			to_string += f'Number of neighbouring atoms: {len(neighbours_to_hydrogen_indices)}\n'
			to_string += f'indices of neighbouring atoms: {neighbours_to_hydrogen_indices}\n'
			to_string += 'Check this'
			raise Exception(to_string)

		# 1.2.5: Get the index of the neighbouring atom to this hydrogen atom
		neighbour_index = neighbours_to_hydrogen_indices[0]

		# 1.2.6: Make sure that hydrogen_index is not already in central_atoms[neighbour_index]
		if (neighbour_index in central_atoms) and (hydrogen_index in central_atoms[neighbour_index]):
			to_string  = f'Error: atom index {hydrogen_index} is already in central_atoms[{neighbour_index}]\n'
			to_string += f'central_atoms[{neighbour_index}] = {central_atoms[{neighbour_index}]}\n'
			to_string += 'Check this'
			raise Exception(to_string)

		# 1.2.7: Make sure that hydrogen_index is not in any of the values in central_atoms
		if any([(hydrogen_index in central_atoms[neighbour_index]) for neighbour_index in central_atoms.keys()]):
			to_string  = f'Error: atom index {hydrogen_index} is already in central_atoms.values()\n'
			to_string += f'central_atoms.values() = {central_atoms.values()}\n'
			to_string += f'central_atoms = {central_atoms}\n'
			to_string += 'Check this'
			raise Exception(to_string)

		# 1.2.8: Place hydrogen_index in central_atoms[neighbour_index]
		central_atoms.setdefault(neighbour_index,[]).append(hydrogen_index)

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Third, make a copy of the molecule that is allowed to be optimised
	copied_molecule = molecule.copy()

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Second, get the atoms, bonds, and angles to fix. 

	# 2.1: Obtain all the atoms that are bonded to each atom in central_atoms.keys() from molecule_graph
	atom_neighbours = {atom_index: tuple(molecule_graph[atom_index].keys()) for atom_index in central_atoms.keys()}

	# 2.2: Atoms involved in fixed bonds
	atoms_indices_involve_in_fixed_bonds = set()

	# 2.3: Initialise the list to add the bond lengths to fix.
	bond_lengths_to_fix = []

	# 2.3: Initialise the list to add the bond angle to fix.
	#bond_angles_to_fix = []

	# 2.4: For each atom in atom_neighbours:
	for central_atom_index, neighbouring_atom_indices in atom_neighbours.items(): 

		# 2.4.1: Get the atom indices of bond lengths to fix.
		for neighbouring_atom_index in neighbouring_atom_indices:
			bond_lengths_to_fix.append(tuple(sorted([central_atom_index, neighbouring_atom_index])))

		# 2.4.2: Fix the bond angles between atoms about central_atom_index:
		for i1 in range(len(neighbouring_atom_indices)):
			for i2 in range(i1+1,len(neighbouring_atom_indices)):
				atom_index_1, atom_index_2 = sorted([neighbouring_atom_indices[i1], neighbouring_atom_indices[i2]])

				# Add indices to atoms_indices_involve_in_fixed_bonds
				atoms_indices_involve_in_fixed_bonds.add(atom_index_1)
				atoms_indices_involve_in_fixed_bonds.add(atom_index_2)

				if (copied_molecule[atom_index_1].symbol not in ('H', 'D', 'T')) or (copied_molecule[atom_index_2].symbol not in ('H', 'D', 'T')):
					continue
				#bond_angles_to_fix.append((atom_index_1, central_atom_index, atom_index_2))
				bond_lengths_to_fix.append((atom_index_1, atom_index_2))



	# 2.5: Remove repeated bond lengths from bond_lengths_to_fix
	#bond_lengths_to_fix = tuple([(molecule.get_distance(index1, index2),     (index1, index2))         for index1, index2         in sorted(set(bond_lengths_to_fix))])
	bond_lengths_to_fix = tuple(sorted(set(bond_lengths_to_fix)))

	# 2.6: Remove repeated bond angles from bond_angles_to_fix
	#bond_angles_to_fix = tuple([(molecule.get_angle(index1, index2, index3), (index1, index2, index3)) for index1, index2, index3 in sorted(set(bond_angles_to_fix)) ])

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



	# Fourth, remove atoms that are far enough away from important atoms that we can ignore them.
	atoms_to_remove_from_model = []
	for index in range(len(copied_molecule)):
		if index in atoms_indices_involve_in_fixed_bonds:
			continue
		for atom_involved_in_rotation_index in atoms_indices_involve_in_fixed_bonds:
			distance = copied_molecule.get_distance(index, atom_involved_in_rotation_index)
			if distance <= 2.0:
				break
		else:
			atoms_to_remove_from_model.append(index)

	# 
	mapping = list(range(len(copied_molecule)))
	for index in sorted(atoms_to_remove_from_model,reverse=True):
		del copied_molecule[index]
		mapping.remove(index)

	#
	mapping         = {new_index: old_index for new_index, old_index in enumerate(mapping)}
	mapping_reverse = {old_index: new_index for new_index, old_index in mapping.items()}



	new_hydrogen_indices_to_allow_free_rotation = [mapping_reverse[index] for index in hydrogen_indices_to_allow_free_rotation]

	new_bond_lengths_to_fix = [(mapping_reverse[index1], mapping_reverse[index2]) for index1, index2 in bond_lengths_to_fix]

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


	# Third, constrain all the atoms that have not been added or modified in the molecule. 
	fixed_atoms = FixAtoms(indices=sorted([atom.index for atom in copied_molecule if atom.index not in new_hydrogen_indices_to_allow_free_rotation]))

	# 4.7: Fix the bond lengths and angles between atoms in bonds_to_fix and bond_angles_to_fix
	#fixed_internals = FixInternals(bonds=bond_lengths_to_fix, angles_deg=bond_angles_to_fix)

	fixed_internals = FixBondLengths(pairs=new_bond_lengths_to_fix)

	# 4.8: Contrain the required atoms, bond lengths, and bond angles. 
	copied_molecule.set_constraint([fixed_internals, fixed_atoms])

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Fifth, apply Thole polarisation values to atoms in the molecule
	molecule_thole_polarisations = [+10.0 for atom in copied_molecule] # [element_thole_polarisations[atom.symbol] for atom in copied_molecule]
	copied_molecule.set_initial_charges(molecule_thole_polarisations)

	# Sixth, set the coulombic calculator to the H-bonding system
	copied_molecule.calc = Coulomb()

	# Seventh, perform an optimisation upon the 
	optimisation = FIRE(copied_molecule, logfile='optimisation.log', trajectory='optimisation.traj')
	to_string = 'Optimising hydrogen positions'
	if molecule_name is not None:
		to_string += ' in molecule '+str(molecule_name)
	if method_type is not None:
		to_string += ' using the '+str(method_type)+' method'
	print(to_string, file=sys.stderr)
	try:
		optimisation.run(fmax=0.001, steps=100)
	except:
		pass

	debugging_molecule = molecule.copy()

	# Eighth, write the new positions of the hydrogens in hydrogen_indices_to_allow_free_rotation into the original molecules object
	molecule_positions = molecule.get_positions()
	for hydrogen_index in hydrogen_indices_to_allow_free_rotation:
		molecule_positions[hydrogen_index] = copied_molecule[mapping_reverse[hydrogen_index]].position
	molecule.set_positions(molecule_positions)


	from ase.visualize import view
	view([debugging_molecule, copied_molecule])
	import pdb; pdb.set_trace()
	raise Exception('Check this method is working.')


	# Ninth, return molecule
	return molecule















