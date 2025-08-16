"""
remove_solvents_from_molecules_list.py, Geoffrey Weal, 6/6/22

This script is designed to remove the solvents from the molecules list.
"""

def remove_solvents_from_molecules_list(molecules, molecule_graphs, SolventsList):
	"""
	This method is designed to remove the solvents from the molecules list.

	Parameters
	----------
	molecules : list of ase.Atoms
		These are the molecules that make up the crystal (before any spacegroup symmetry is considered).
	molecule_graphs : list of networkx.Graph
		These are the graph that are associated with the molecules in the molecules list. 
	SolventsList : list of ints
		This list indicates which molecules in the molecules list are solvents.
	
	Returns
	-------
	molecules : list of ase.Atoms
		These are the all molecules that make up the crystal, including those molecules that exist due to symmetry operations of the crystal's spacegroup. 
	molecule_graphs : list of networkx.Graph
		These are the graph that are associated with the molecules in the molecules list. 
	"""

	raise Exception('Check if this method is being used')

	# First, remove the molecules and graphs based on the molecule indices in the SolventsList list. 
	for index in sorted(SolventsList,reverse=True):
		del molecules[index]
		del molecule_graphs[index]

	# Second, return the lists.
	return molecules, molecule_graphs