"""
allow_atoms_to_freely_rotate.py, Geoffrey Weal, 15/4/2024

This method is designed to allow selected atoms to freely rotate about the molecules. 
"""
from ase.optimize                      import FIRE
from ase.constraints                   import FixAtoms, FixInternals
from SUMELF.SUMELF.calculators.Coulomb import Coulomb

def allow_atoms_to_freely_rotate(molecule, molecule_graph, atom_indices_to_allow_free_rotation, connecting_indices=[], molecule_name=None, method_type=None):
	"""
	This method is designed to allow selected atoms to freely rotate about the molecules. 

	Parameters
	----------
	molecule : ase.Atoms
		This is the molecule you want to allow atoms to freely rotate.
	molecule_graph : networkx.Graph
		This is the graph of this molecule.
	atom_indices_to_allow_free_rotation : list of ints
		These are the atoms you would like to allow to freely rotate.
	molecule_name : int
		This is the name of the molecule that methyl groups are being added to.
	method_type : str
		This is the method that is asking for atoms to freely rotate. 

	Returns
	-------
	updated_molecule : ase.Atoms
		This is the molecule of interest, where atoms in atom_indices_to_allow_free_rotation have been allowed to freely rotate. 
	"""

	# First, make a copy of the molecule that is allowed to be optimised
	copied_molecule = molecule.copy()

	# Second, constrain all the atoms that have not been added or modified in the molecule. 
	fixed_atoms = FixAtoms(indices=sorted([atom.index for atom in copied_molecule if (atom.index not in atom_indices_to_allow_free_rotation)]))

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Third, get the atoms, bonds, and angles to fix. 

	# 3.1: Obtain all the atoms that are bonded to each atom in central_atoms.keys() from molecule_graph
	atom_neighbours = {atom_index: tuple(molecule_graph[atom_index].keys()) for atom_index in (atom_indices_to_allow_free_rotation + connecting_indices)}

	# 3.2: Initialise the list to add the bond lengths to fix.
	bond_lengths_to_fix = []

	# 3.3: Initialise the list to add the bond angle to fix.
	bond_angles_to_fix = []

	# 3.4: For each atom in atom_neighbours:
	for central_atom_index, neighbouring_atom_indices in atom_neighbours.items(): 

		# 3.4.1: Get the atom indices of bond lengths to fix.
		for neighbouring_atom_index in neighbouring_atom_indices:
			bond_lengths_to_fix.append(tuple(sorted([central_atom_index, neighbouring_atom_index])))

		# 3.4.2: Fix the bond angles between atoms about central_atom_index:
		for i1 in range(len(neighbouring_atom_indices)):
			for i2 in range(i1+1,len(neighbouring_atom_indices)):
				atom_index_1, atom_index_2 = sorted([neighbouring_atom_indices[i1], neighbouring_atom_indices[i2]])
				bond_angles_to_fix.append((atom_index_1, central_atom_index, atom_index_2))

	# 3.5: Remove repeated bond lengths from bond_lengths_to_fix
	bond_lengths_to_fix = tuple([(molecule.get_distance(index1, index2),     (index1, index2))         for index1, index2         in sorted(set(bond_lengths_to_fix))])

	# 3.6: Remove repeated bond angles from bond_angles_to_fix
	bond_angles_to_fix = tuple([(molecule.get_angle(index1, index2, index3), (index1, index2, index3)) for index1, index2, index3 in sorted(set(bond_angles_to_fix)) ])

	# 3.7: Fix the bond lengths and angles between atoms in bonds_to_fix and bond_angles_to_fix
	fixed_internals = FixInternals(bonds=bond_lengths_to_fix, angles_deg=bond_angles_to_fix)

	# 3.8: Contrain the required atoms, bond lengths, and bond angles. 
	copied_molecule.set_constraint([fixed_internals, fixed_atoms])

	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	# Fourth, apply Thole polarisation values to atoms in the molecule
	molecule_thole_polarisations = [+10 for atom in copied_molecule] # [element_thole_polarisations[atom.symbol] for atom in copied_molecule]
	copied_molecule.set_initial_charges(molecule_thole_polarisations)

	# Fifth, set the coulombic calculator to the H-bonding system
	copied_molecule.calc = Coulomb()

	# Sixth, perform an optimisation upon the 
	optimisation = FIRE(copied_molecule, logfile=None) #trajectory='test.traj', logfile=None)
	to_string = 'Optimising atom positions'
	if molecule_name is not None:
		to_string += ' in molecule '+str(molecule_name)
	if method_type is not None:
		to_string += ' using the '+str(method_type)+' method'
	try:
		optimisation.run(fmax=0.001, steps=100)
	except:
		pass

	debugging_molecule = molecule.copy()

	# Seventh, write the new positions of the hydrogens in hydrogen_indices_to_allow_free_rotation into the original molecules object
	molecule_positions = molecule.get_positions()
	for atom_index in atom_indices_to_allow_free_rotation:
		molecule_positions[atom_index] = copied_molecule[atom_index].position
	molecule.set_positions(molecule_positions)

	#'''
	from ase.visualize import view
	view([debugging_molecule, molecule])
	import pdb; pdb.set_trace()
	raise Exception('Check this method is working Angles are probably not conserved.')
	#'''

	# Eighth, return molecule
	return molecule



