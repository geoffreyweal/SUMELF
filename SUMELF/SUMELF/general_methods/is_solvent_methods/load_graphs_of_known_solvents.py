"""
load_graphs_of_known_solvents.py, Geoffrey Weal, 15/2/24

This script is designed to load molecules from the graph txt files into memory. 
"""
import os
from SUMELF.SUMELF.general_methods.is_solvent_methods.solvent_graph_methods import read_graph, create_graph_from_mol_files

# These variables indicate the folders where solvent files are placed.
general_path_to_solvent_mol_files   = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'solvent_files', 'solvent_mol_files'))
general_path_to_solvent_graph_files = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'solvent_files', 'solvent_graph_files'))

def load_graphs_of_known_solvents():
	"""
	This method is designed to load molecules from the graph txt files into memory. 

	Returns
	-------
	solvent_graphs : list of networkx.Graph
		This is a list of all the graphs of the solvents in the graph folder. 
	"""
	
	# First, make sure that the graph folder exists.
	if not os.path.exists(general_path_to_solvent_graph_files):

		# 1.1: If the graph folder does not exist, create it and convert all mol files to graphs and save it into here. 
		create_graph_from_mol_files()

	# Second, initialise the solvents list/
	solvent_graphs = []

	# Third, for each file in the graphs folder. 
	for file in os.listdir(general_path_to_solvent_graph_files):

		# Fourth, reach the graph file of the solvent. 
		solvent_graph = read_graph(general_path_to_solvent_graph_files+'/'+file.replace('.mol','.txt'))

		# Fifth, append the solvent graph to the solvent_graphs list. 
		solvent_graphs.append(solvent_graph)

	# Sixth, return solvent_graphs
	return solvent_graphs