"""
get_molecules.py, Geoffrey Weal, 1/4/24

This subsidiary program is designed to obtain the structurally unique molecules from the crystal files obtained from the ACSD program

This is designed to allow the user to check the quality of the molecules given in the crystal. This helps to determine how to fix the atoms in the crystal. 
"""
import os, shutil
from tqdm   import tqdm
from ase.io import read, write
from SUMELF import obtain_graph, process_crystal, add_graph_to_ASE_Atoms_object
#from ECCP.ECCP.get_simple_molecule_graphs import get_simple_molecule_graphs
#from SUMELF import obtain_unique_molecules
#from ECCP.ECCP.get_unique_molecules import get_unique_molecules

# This is the default name of the folder containing the crystal structures.
crystal_database_foldername_DEFAULT = ['crystal_database']
output_molecules_crystal_database_foldername_DEFAULT = [crystal_database_foldername_DEFAULT[0]+'_molecules']

class CLICommand:
	"""This method is designed to repair crystals that you have indicated need repairing, based on how you have indicated your crystal needs to be repaired in the "repair_crystals.py" file.
	"""

	@staticmethod
	def add_arguments(parser):
		parser.add_argument('--crystal_database_foldername',                  nargs=1, help='This is the folder that contains the crystal xyz files you want to split into individual molecules.', default=crystal_database_foldername_DEFAULT)
		parser.add_argument('--output_molecules_crystal_database_foldername', nargs=1, help='This is the folder to save the xyz files of individual molecules into.', default=output_molecules_crystal_database_foldername_DEFAULT)

	@staticmethod
	def run(args):
		
		# First, obtain the variables for crystal_database_foldername and output_molecules_crystal_database_foldername from args.
		crystal_database_foldername                  = args.crystal_database_foldername[0]
		output_molecules_crystal_database_foldername = args.output_molecules_crystal_database_foldername[0]

		# Second, run the get_molecules method.
		Run_method(crystal_database_foldername=crystal_database_foldername, output_molecules_crystal_database_foldername=output_molecules_crystal_database_foldername)

def Run_method(crystal_database_foldername=crystal_database_foldername_DEFAULT, output_molecules_crystal_database_foldername=output_molecules_crystal_database_foldername_DEFAULT, molecule_equivalence_method={'method': 'invariance_method', 'type': 'combination'}, no_of_cpus=1):
	"""
	This method is designed to obtain the structurally unique molecules from the crystal files obtained from the ACSD program

	This is designed to allow the user to check the quality of the molecules given in the crystal. This helps to determine how to fix the atoms in the crystal. 

	Parameters
	----------
	crystal_database_foldername : str.
		This is the path to the folder that contains all the crystal files created by the ACSD program.
	make_molecule_method : str.
		This is the name of the method you want to use to create the molecule. See https://github.com/geoffreyweal/ECCP for more information. Default: 'component_assembly_approach'. 
	"""

	# First, obtain the list of crystal files to examine.
	crystal_files = sorted([file for file in os.listdir(crystal_database_foldername) if (os.path.isfile(crystal_database_foldername+'/'+file) and file.endswith('.xyz'))])

	# Second, create a folder that contains all the molecules files for the crystals in crystal_database_foldername
	if os.path.exists(output_molecules_crystal_database_foldername):
		shutil.rmtree(output_molecules_crystal_database_foldername)
	os.makedirs(output_molecules_crystal_database_foldername)

	# Third, print to the user what you are doing.
	print(f'Extracting the molecules from the crystal files in {crystal_database_foldername}')
	print(f'Crystals to extract: {len(crystal_files)}')

	# Fourth, make the progress bar for processing crystals.
	pbar = tqdm(crystal_files, desc='Processing Molecules from Crystal Files', unit=' crystals')

	# Fifth, for each crystal in crystal_files
	for crystal_file in pbar:

		# 5.1: Get the name of the crystal file.
		crystal_name = crystal_file.replace('.xyz','')

		# 5.3: Indicate what crystal is being processed
		pbar.set_description(crystal_name)

		# 5.3: Get the crystal file.
		crystal = read(crystal_database_foldername+'/'+crystal_file)

		# 5.4: Get the graph of the crystal.
		crystal, crystal_graph = obtain_graph(crystal,name=crystal_name)

		# 5.5: Get the molecules and the graphs associated with each molecule in the crystal.
		molecules, molecule_graphs, SolventsList, symmetry_operations, cell = process_crystal(crystal,crystal_graph=crystal_graph,take_shortest_distance=True,return_list=False,logger=None,print_progress=False)

		'''
		# 5.6: Obtain the simplified graph for the crystals
		simple_molecule_graphs = get_simple_molecule_graphs(molecule_graphs)

		# 5.7: Only include structurally unique molecules on the molecules recorded to disk.
		structurally_unique_molecules_indices, structurally_equivalent_molecule_groups, _, _, _ = get_unique_molecules(molecules, simple_molecule_graphs, crystal, molecule_equivalence_method=molecule_equivalence_method, neighbouring_molecules_about_molecules=None, include_hydrogens_in_uniqueness_analysis=False, no_of_cpus=no_of_cpus)
		'''

		# 5.6: Add the node and edge information from the molecules graph back to the molecule
		for molecule_name in molecules.keys():
			add_graph_to_ASE_Atoms_object(molecules[molecule_name], molecule_graphs[molecule_name])

		# 5.7: Create the folder to store molecule xyz data to.
		os.makedirs(output_molecules_crystal_database_foldername+'/'+crystal_name)

		# 5.8: Save each molecule from the crystal to disk
		for molecule_name, molecule in molecules.items():
			write(output_molecules_crystal_database_foldername+'/'+crystal_name+'/'+str(molecule_name)+'.xyz', molecule)