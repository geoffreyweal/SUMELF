"""
add_charges_to_molecules.py, Geoffrey Weal, 20/3/22

This script will assign charges to the atoms in each molecule.
"""
#from numpy import array

def add_zero_charges_to_molecules(molecules, molecule_graphs, max_disparity=None):
	"""
	This method is designed to assign 0 charges to each atom in the molecules in your crystal.

	Parameters
	----------
	molecules : list of ase.Atoms
		This is a list of all the molecules in your crystal.
	molecule_graphs : list of networkx.Graph
		This is the list of the undirected graph representations of each molecule.
	max_disparity : float
		This is the maximum disparity between two dimers to be considered invariant. 
	"""

	# First, apply zero charges from the ATC to the appropriate molecule.
	for mol_index in range(len(molecules)):
		molecules[mol_index].set_initial_charges([0]*len(molecules[mol_index]))

# ------------------------------------------------------------------------------------------------------------------------
