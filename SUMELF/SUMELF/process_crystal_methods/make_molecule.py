"""
make_molecule.py, Geoffrey Weal, 17/2/22

This script will create the individual molecules that are found within a crystal structure.
"""
from ase import Atoms

from SUMELF.SUMELF.process_crystal_methods.make_molecule_methods.super_lattice_approach      import create_molecule_using_super_lattice_approach
from SUMELF.SUMELF.process_crystal_methods.make_molecule_methods.component_assembly_approach import create_molecule_using_component_assembly_approach

def make_molecule(atom_indices_in_crystal_to_extract_as_molecule, crystal_graph, crystal, take_shortest_distance=False, make_molecule_method='component_assembly_approach', logger=None):
	"""
	This method will create the individual molecules that are found within a crystal structure.

	This method does this using the results from the results given when the crystals graph was put through the connected_components method.

	This method also give the option of removing the aliphatic side chains from the molecules.

	Parameters
	----------
	atom_indices_in_crystal_to_extract_as_molecule : list
		A list of indices of the atoms in each molecule in the crystal. The indices given are the indices of atoms in the crystal
	crystal_graph : networkx.Graph
		The graph for the crystal
	crystal : ase.Atoms
		The crystal in Atomic Simulation Environment format
	take_shortest_distance : bool
		If true, take the shortest distance that was found if a distance shorter than the maximum distance between elements if found. If false, give an error if the distance between two atoms was greater than the max distance expected for those two elements to be bonded to each other. Default: False
	make_molecule_method : str.
		This is the name of the method you want to use to create the molecule. Their are two options for this: 'super_lattice_approach' and 'component_assembly_approach'. See https://github.com/geoffreyweal/SUMELF for more information. Default: 'component_assembly_approach'. 
	logger : logging/None
		This contains the object for reporting information and warning messages to.

	Returns
	-------
	complete_molecule : ase.Atoms
		This is the molecule in bonded form, with or without aliphatic side chains based on your choice
	molecule_graph : networkx.Graph
		The updated graph for the molecule in the crystal
	"""

	# First, get information about crystal.
	crystal_cell_lattice = crystal.get_cell()

	# Second, get the tags from the crystal (if they were given in the crystal).
	if 'tags' in crystal.arrays.keys():
		crystal_original_tags = crystal.get_tags()

	# Third, set up the molecule in the current crystal.
	molecule_in_crystal = Atoms()
	molecule_in_crystal.set_cell(crystal_cell_lattice)

	# Fourth, create a mapping to map the atom in the crystal to the atom in the molecule
	mapping_atom_in_crystal_TO_atom_in_molecule_indices = {}
	for atom_index_in_molecule in range(len(atom_indices_in_crystal_to_extract_as_molecule)):

		# 4.1: Obtain the index of the atom in the crystal.
		atom_in_crystal_index = atom_indices_in_crystal_to_extract_as_molecule[atom_index_in_molecule]

		# 4.2: Obtain the atom from the crystal and add it to the molecule.
		molecule_in_crystal.append(crystal[atom_in_crystal_index])

		# 4.3: Map atom_in_crystal_index --> atom_index_in_molecule in the mapping dictionary.
		mapping_atom_in_crystal_TO_atom_in_molecule_indices[atom_in_crystal_index] = atom_index_in_molecule

	# Fifth, the molecule may be disconnected as it is given in reference to the crystal. 
	#        * Here, we will recreate the molecule, keeping the atoms in the same positions in the crystal. 
	#        * This will retain the molecules positional information in the crystal, while making the molecule easier to process and look at.
	if make_molecule_method == 'super_lattice_approach':
		complete_molecule, molecule_graph = create_molecule_using_super_lattice_approach(molecule_in_crystal, crystal_graph, crystal, mapping_atom_in_crystal_TO_atom_in_molecule_indices, take_shortest_distance=take_shortest_distance, to_print=True)
	elif make_molecule_method == 'component_assembly_approach':
		complete_molecule, molecule_graph = create_molecule_using_component_assembly_approach(molecule_in_crystal, crystal_graph, crystal, mapping_atom_in_crystal_TO_atom_in_molecule_indices, take_shortest_distance=take_shortest_distance, to_print=True, logger=logger)
	else:
		print('Error: The input method for the make molecule method can be either:')
		print('\t* super_lattice_approach: The Super Lattice Method')
		print('\t* component_assembly_approach: The Component Assembly Method')
		print('See https://github.com/geoffreyweal/ECCP for more information')
		exit('This program will finish without completing')

	# Sixth, place the original tags back on each atom in this molecule (if they were given in the crystal).
	if 'tags' in crystal.arrays.keys():
		tags = [crystal_original_tags[index] for index in atom_indices_in_crystal_to_extract_as_molecule]
		complete_molecule.set_tags(tags)

	# Seventh, return the connected molecule and the molecule with all atoms contained within the origin crystal, and the associated graph for the molecule.
	return complete_molecule, molecule_graph

# ---------------------------------------------------------------------------------------------------------------------------


