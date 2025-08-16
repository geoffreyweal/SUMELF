# The information about the Supporting Methods for Electronic Functions (SUMELF) program

__name__    = 'SUMELF (Supporting Methods for Electronic Functions) Program'
__version__ = '0.84'
__author__  = 'Dr. Geoffrey Weal, Dr. Chayanit Wechwithayakhlung, Dr. Josh Sutton, Dr. Daniel Packwood, Dr. Paul Hume, Prof. Justin Hodgkiss'

import sys, importlib

if sys.version_info[0] == 2:
	toString = ''
	toString += '\n'
	toString += '================================================'+'\n'
	toString += 'This is the Supporting Methods for Electronic Functions (SUMELF) program'+'\n'
	toString += 'Version: '+str(__version__)+'\n'
	toString += '\n'
	toString += 'The Supporting Methods for Electronic Functions (SUMELF) program requires Python3. You are attempting to execute this program in Python2.'+'\n'
	toString += 'Make sure you are running the Supporting Methods for Electronic Functions (SUMELF) program in Python3 and try again'+'\n'
	toString += 'This program will exit before beginning'+'\n'
	toString += '================================================'+'\n'
	raise ImportError(toString)
if sys.version_info[1] < 4:
	toString = ''
	toString += '\n'
	toString += '================================================'+'\n'
	toString += 'This is the Supporting Methods for Electronic Functions (SUMELF) program'+'\n'
	toString += 'Version: '+str(__version__)+'\n'
	toString += '\n'
	toString += 'The Supporting Methods for Electronic Functions program requires Python 3.4 or greater.'+'\n'
	toString += 'You are using Python '+str('.'.join(sys.version_info))
	toString += '\n'
	toString += 'Use a version of Python 3 that is greater or equal to Python 3.4.\n'
	toString += 'This program will exit before beginning'+'\n'
	toString += '================================================'+'\n'
	raise ImportError(toString)

# ------------------------------------------------------------------------------------------------------------------------
# A check for ASE

ase_spec = importlib.util.find_spec("ase")
ase_found = (ase_spec is not None)
if not ase_found:
	toString = ''
	toString += '\n'
	toString += '================================================'+'\n'
	toString += 'This is the Supporting Methods for Electronic Functions (SUMELF) program'+'\n'
	toString += 'Version: '+str(__version__)+'\n'
	toString += '\n'
	toString += 'The Supporting Methods for Electronic Functions program requires ASE.'+'\n'
	toString += '\n'
	toString += 'Install ASE through pip by following the instruction in https://github.com/GardenGroupUO/SUMELF'+'\n'
	toString += 'These instructions will ask you to install ase by typing the following into your terminal\n'
	toString += 'pip3 install --user --upgrade ase\n'
	toString += '\n'
	toString += 'This program will exit before beginning'+'\n'
	toString += '================================================'+'\n'
	raise ImportError(toString)	

import ase
ase_version_minimum = '3.19.0'
from packaging import version
if version.parse(ase.__version__) < version.parse(ase_version_minimum):
	toString = ''
	toString += '\n'
	toString += '================================================'+'\n'
	toString += 'This is the Supporting Methods for Electronic Functions (SUMELF) program'+'\n'
	toString += 'Version: '+str(__version__)+'\n'
	toString += '\n'
	toString += 'The SUMELF program requires ASE greater than or equal to '+str(ase_version_minimum)+'.'+'\n'
	toString += 'The current version of ASE you are using is '+str(ase.__version__)+'.'+'\n'
	toString += '\n'
	toString += 'Install ASE through pip by following the instruction in https://github.com/GardenGroupUO/SUMELF'+'\n'
	toString += 'These instructions will ask you to install ase by typing the following into your terminal\n'
	toString += 'pip3 install --user --upgrade ase\n'
	toString += '\n'
	toString += 'This program will exit before beginning'+'\n'
	toString += '================================================'+'\n'
	raise ImportError(toString)

# ------------------------------------------------------------------------------------------------------------------------

scipy_spec = importlib.util.find_spec("scipy")
scipy_found = (scipy_spec is not None)
if not scipy_found:
	toString = ''
	toString += '\n'
	toString += '================================================'+'\n'
	toString += 'This is the Supporting Methods for Electronic Functions (SUMELF) program'+'\n'
	toString += 'Version: '+str(__version__)+'\n'
	toString += '\n'
	toString += 'The SUMELF program requires the "scipy" program.'+'\n'
	toString += '\n'
	toString += 'Install scipy by typing the following into your terminal\n'
	toString += '\n'
	toString += 'pip3 install --user --upgrade scipy\n'
	toString += '\n'
	toString += 'This program will exit before beginning'+'\n'
	toString += '================================================'+'\n'
	raise ImportError(toString)	

# ------------------------------------------------------------------------------------------------------------------------

networkx_spec = importlib.util.find_spec("networkx")
networkx_found = (networkx_spec is not None)
if not networkx_found:
	toString = ''
	toString += '\n'
	toString += '================================================'+'\n'
	toString += 'This is the Supporting Methods for Electronic Functions (SUMELF) program'+'\n'
	toString += 'Version: '+str(__version__)+'\n'
	toString += '\n'
	toString += 'The SUMELF program requires the "networkx" program.'+'\n'
	toString += '\n'
	toString += 'Install networkx through pip by following the instruction in https://github.com/GardenGroupUO/SUMELF'+'\n'
	toString += 'These instructions will ask you to install networkx by typing the following into your terminal\n'
	toString += '\n'
	toString += 'pip3 install --user --upgrade networkx\n'
	toString += '\n'
	toString += 'This program will exit before beginning'+'\n'
	toString += '================================================'+'\n'
	raise ImportError(toString)	

# ------------------------------------------------------------------------------------------------------------------------

pymatgen_spec = importlib.util.find_spec("pymatgen")
pymatgen_found = (pymatgen_spec is not None)
if not pymatgen_found:
	toString = ''
	toString += '\n'
	toString += '================================================'+'\n'
	toString += 'This is the Supporting Methods for Electronic Functions (SUMELF) program'+'\n'
	toString += 'Version: '+str(__version__)+'\n'
	toString += '\n'
	toString += 'The SUMELF program requires the "pymatgen" program.'+'\n'
	toString += '\n'
	toString += 'Install pymatgen through pip by following the instruction in https://github.com/GardenGroupUO/SUMELF'+'\n'
	toString += 'These instructions will ask you to install pymatgen by typing the following into your terminal\n'
	toString += '\n'
	toString += 'pip3 install --user --upgrade pymatgen\n'
	toString += '\n'
	toString += 'This program will exit before beginning'+'\n'
	toString += '================================================'+'\n'
	raise ImportError(toString)	

# ------------------------------------------------------------------------------------------------------------------------

packaging_spec = importlib.util.find_spec("packaging")
packaging_found = (packaging_spec is not None)
if not packaging_found:
	toString = ''
	toString += '\n'
	toString += '================================================'+'\n'
	toString += 'This is the Supporting Methods for Electronic Functions (SUMELF) program'+'\n'
	toString += 'Version: '+str(__version__)+'\n'
	toString += '\n'
	toString += 'The SUMELF program requires the "packaging" program.'+'\n'
	toString += '\n'
	toString += 'Install packaging through pip by following the instruction in https://github.com/GardenGroupUO/SUMELF'+'\n'
	toString += 'These instructions will ask you to install packaging by typing the following into your terminal\n'
	toString += '\n'
	toString += 'pip3 install --user --upgrade packaging\n'
	toString += '\n'
	toString += 'This program will exit before beginning'+'\n'
	toString += '================================================'+'\n'
	raise ImportError(toString)	

# ------------------------------------------------------------------------------------------------------------------------

tqdm_spec = importlib.util.find_spec("tqdm")
tqdm_found = (tqdm_spec is not None)
if not tqdm_found:
	toString = ''
	toString += '\n'
	toString += '================================================'+'\n'
	toString += 'This is the Supporting Methods for Electronic Functions (SUMELF) program'+'\n'
	toString += 'Version: '+str(__version__)+'\n'
	toString += '\n'
	toString += 'The SUMELF program requires the "tqdm" program.'+'\n'
	toString += '\n'
	toString += 'Install tqdm through pip by following the instruction in https://github.com/GardenGroupUO/SUMELF'+'\n'
	toString += 'These instructions will ask you to install tqdm by typing the following into your terminal\n'
	toString += '\n'
	toString += 'pip3 install --user --upgrade tqdm\n'
	toString += '\n'
	toString += 'This program will exit before beginning'+'\n'
	toString += '================================================'+'\n'

# ------------------------------------------------------------------------------------------------------------------------

tqdm_spec = importlib.util.find_spec("xlsxwriter")
tqdm_found = (tqdm_spec is not None)
if not tqdm_found:
	toString = ''
	toString += '\n'
	toString += '================================================'+'\n'
	toString += 'This is the Supporting Methods for Electronic Functions (SUMELF) program'+'\n'
	toString += 'Version: '+str(__version__)+'\n'
	toString += '\n'
	toString += 'The SUMELF program requires the "xlsxwriter" program.'+'\n'
	toString += '\n'
	toString += 'Install xlsxwriter through pip by following the instruction in https://github.com/GardenGroupUO/SUMELF'+'\n'
	toString += 'These instructions will ask you to install xlsxwriter by typing the following into your terminal\n'
	toString += '\n'
	toString += 'pip3 install --user --upgrade xlsxwriter\n'
	toString += '\n'
	toString += 'This program will exit before beginning'+'\n'
	toString += '================================================'+'\n'

# ------------------------------------------------------------------------------------------------------------------------

mprof_spec = importlib.util.find_spec("memory_profiler")
mprof_found = (tqdm_spec is not None)
if not mprof_found:
	toString = ''
	toString += '\n'
	toString += '================================================'+'\n'
	toString += 'This is the Supporting Methods for Electronic Functions (SUMELF) program'+'\n'
	toString += 'Version: '+str(__version__)+'\n'
	toString += '\n'
	toString += 'The SUMELF program requires the "memory_profiler" program.'+'\n'
	toString += '\n'
	toString += 'Install memory_profiler through pip by following the instruction in https://github.com/GardenGroupUO/SUMELF'+'\n'
	toString += 'These instructions will ask you to install memory_profiler by typing the following into your terminal\n'
	toString += '\n'
	toString += 'pip3 install --user --upgrade memory_profiler\n'
	toString += '\n'
	toString += 'This program will exit before beginning'+'\n'
	toString += '================================================'+'\n'

# ------------------------------------------------------------------------------------------------------------------------

__author_email__ = 'geoffrey.weal@vuw.ac.nz'
__license__ = 'GNU AFFERO GENERAL PUBLIC LICENSE'
__url__ = 'https://github.com/GardenGroupUO/SUMELF'
__doc__ = 'See https://github.com/GardenGroupUO/SUMELF for the documentation on this program'

# ================================================================================================
# Making and Processing Crystal Methods
from SUMELF.SUMELF.process_crystal                                               import process_crystal
from SUMELF.SUMELF.make_crystal                                                  import make_crystal
from SUMELF.SUMELF.get_spacegroup_molecules                                      import get_spacegroup_molecules
# ================================================================================================
# Are molecules inbetween dimers methods
from SUMELF.SUMELF.molecules_inbetween_dimers                                    import molecules_inbetween_dimers
# ================================================================================================
# General methods
from SUMELF.SUMELF.general_methods.angle_methods                                 import get_angle, above_angle_tolerance
from SUMELF.SUMELF.general_methods.check_molecule_against_file                   import check_molecule_against_file
from SUMELF.SUMELF.general_methods.convert_between_ijk_and_cell_position         import convert_ijk_to_displacement_vector, convert_position_into_ijk_lengths
from SUMELF.SUMELF.general_methods.distance_methods                              import get_distance, get_xyz_distances, less_than_or_equal_to_max_bondlength, are_two_values_within_eachother, are_two_lists_within_eachother
from SUMELF.SUMELF.general_methods.general_molecules_methods                     import read_crystal, get_centre_of_mass, get_centre_of_molecule, get_number_of_lone_pairs_of_electron_pairs, get_hybridisation_from_ASE, get_hybridisation_from_CSD, get_bond_type_from_CSD, get_atomic_rings_in_ase_object, get_translation_to_move_COM_inside_unit_cell, is_the_same_molecule_exact
from SUMELF.SUMELF.general_methods.geometry_methods                              import get_unit_vector, rotate_vector_around_axis, project_point_onto_line, project_u_onto_v, project_point_onto_plane, get_reflection_matrix_from_plane, planeFit, project_point_onto_plane, get_rotation_matrix_around_arbitrary_axis, get_cross_product_matrix
from SUMELF.SUMELF.general_methods.get_symmetry_operations                       import get_symmetry_operations
from SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods         import get_bond_lengths, get_bond_angle
from SUMELF.SUMELF.general_methods.is_solvent                                    import is_solvent
from SUMELF.SUMELF.general_methods.make_dimer                                    import make_dimer
from SUMELF.SUMELF.general_methods.remove_hydrogens                              import remove_hydrogens
from SUMELF.SUMELF.general_methods.unique_system_methods                         import obtain_unique_molecules, obtain_unique_dimers
from SUMELF.SUMELF.general_methods.unit_cell_methods                             import get_cell_corner_points, centre_molecule_in_cell, number_of_points_in_cell
# ================================================================================================
# Methods relating to graphs associated with ASE Atoms objects
from SUMELF.SUMELF.graph_methods.add_graph_to_ASE_Atoms_object                   import add_graph_to_ASE_Atoms_object
from SUMELF.SUMELF.graph_methods.get_properties_from_graph                       import get_node_properties_from_graph, remove_node_properties_from_graph
from SUMELF.SUMELF.graph_methods.obtain_graph                                    import obtain_graph
# ================================================================================================
# Method for graph matching between molecules
from SUMELF.SUMELF.isomorphvf2_GRW                                               import GraphMatcher
# ================================================================================================
# Methods from make_crystal_methods folder
from SUMELF.SUMELF.make_crystal_methods.determine_is_molecule_already_recorded   import determine_is_molecule_already_recorded
from SUMELF.SUMELF.make_crystal_methods.perform_symmetry_operation_upon_molecule import perform_symmetry_operation_upon_molecule
# ================================================================================================
# Methods for adding Hydrogens and Methyl groups to a molecule.
from SUMELF.SUMELF.add_atoms.add_hydrogens_to_molecules                          import add_hydrogens_to_molecules
from SUMELF.SUMELF.add_atoms.add_methyls_to_molecules                            import add_methyls_to_molecules
from SUMELF.SUMELF.add_atoms.add_ethyls_to_molecules                             import add_ethyls_to_molecules
# ================================================================================================
# Mehthods for obtaining ASE objects from chg files, and adding charges to molecules in other systems.
from SUMELF.SUMELF.add_charges_to_molecules_methods.import_CHG_file              import import_CHG_file
from SUMELF.SUMELF.add_charges_to_molecules_methods.invariance_method            import assigned_ATCs_to_molecules_invariance_method
# ================================================================================================
# Cell Generator methods
from SUMELF.SUMELF.generators.Cell_Generator                                     import Cell_Generator
# ================================================================================================
# Calculators
from SUMELF.SUMELF.calculators.Coulomb                                           import Coulomb
# ================================================================================================
# Utility methods
from SUMELF.SUMELF.utility_methods.folder_methods                                import make_folder, remove_folder, move_folder
from SUMELF.SUMELF.utility_methods.get_SameMoleculesDueToCrystalSymmetry         import get_SameMoleculesDueToCrystalSymmetry
# ================================================================================================
# Methods for reading Electronic_Crystal_Calculation_Prep files 
from SUMELF.SUMELF.ECCP_based_methods.what_ECCP_Information_files_do_we_have     import what_ECCP_Information_files_do_we_have
from SUMELF.SUMELF.ECCP_based_methods.run_which_ECCP_operations                  import run_which_ECCP_operations
from SUMELF.SUMELF.ECCP_based_methods.get_ECCP_Information_data                  import get_ECCP_Information_data
from SUMELF.SUMELF.ECCP_based_methods.check_ECCP_Information_details             import check_ECCP_Information_details
from SUMELF.SUMELF.ECCP_based_methods.get_crystal_file_from_ECCP_Information     import get_crystal_file_from_ECCP_Information
from SUMELF.SUMELF.ECCP_based_methods.get_dimer_details_data                     import get_dimer_details_data
from SUMELF.SUMELF.ECCP_based_methods.get_equivalent_molecule_group_data         import get_equivalent_molecule_group_data
from SUMELF.SUMELF.ECCP_based_methods.get_equivalent_dimer_group_data            import get_equivalent_dimer_group_data
from SUMELF.SUMELF.ECCP_based_methods.check_consistancy_between_files            import check_consistancy_between_files
# ================================================================================================

# Making and Processing Crystal Methods
process_crystal_methods                           = ['process_crystal']
make_crystal_methods                              = ['make_crystal']
get_spacegroup_molecules_methods                  = ['get_spacegroup_molecules']
making_and_processing_crystal_methods             = process_crystal_methods + make_crystal_methods + get_spacegroup_molecules_methods

# Are molecules inbetween dimers methods
molecules_inbetween_dimers_methods                = ['molecules_inbetween_dimers']

# General methods
angle_methods                                     = ['get_angle', 'above_angle_tolerance']
distance_methods                                  = ['get_distance', 'get_xyz_distances', 'less_than_or_equal_to_max_bondlength', 'are_two_values_within_eachother', 'are_two_lists_within_eachother']
general_molecules_methods                         = ['read_crystal', 'get_centre_of_mass', 'get_centre_of_molecule', 'get_number_of_lone_pairs_of_electron_pairs',  'get_hybridisation_from_ASE', 'get_hybridisation_from_CSD', 'get_bond_type_from_CSD', 'get_atomic_rings_in_ase_object', 'get_translation_to_move_COM_inside_unit_cell', 'is_the_same_molecule_exact']
geometry_methods                                  = ['get_unit_vector', 'rotate_vector_around_axis', 'project_point_onto_line', 'project_u_onto_v', 'project_point_onto_plane', 'get_reflection_matrix_from_plane', 'planeFit', 'project_point_onto_plane', 'get_rotation_matrix_around_arbitrary_axis', 'get_cross_product_matrix']
get_symmetry_operations_methods                   = ['get_symmetry_operations']
ideal_bond_lengths_and_angles_methods             = ['get_bond_lengths', 'get_bond_angle']
is_solvent_methods                                = ['is_solvent']
make_dimer_methods                                = ['make_dimer']
remove_hydrogens_methods                          = ['remove_hydrogens']
unique_system_methods                             = ['obtain_unique_molecules', 'obtain_unique_dimers']
unit_cell_methods                                 = ['get_cell_corner_points', 'centre_molecule_in_cell', 'number_of_points_in_cell']
general_methods                                   = angle_methods + distance_methods + general_molecules_methods + geometry_methods + get_symmetry_operations_methods + ideal_bond_lengths_and_angles_methods + is_solvent_methods + make_dimer_methods + remove_hydrogens_methods + unique_system_methods + unit_cell_methods

# Methods relating to graphs associated with ASE Atoms objects
add_graph_to_ASE_Atoms_object_methods             = ['add_graph_to_ASE_Atoms_object']
get_properties_from_graph_methods                 = ['get_node_properties_from_graph', 'remove_node_properties_from_graph']
obtain_graph_methods                              = ['obtain_graph']
ase_associated_graph_methods                      = add_graph_to_ASE_Atoms_object_methods + get_properties_from_graph_methods + obtain_graph_methods

# Method for graph matching between molecules
graph_matching_methods                            = ['GraphMatcher']

# Methods from make_crystal_methods folder
determine_is_molecule_already_recorded_methods    = ['determine_is_molecule_already_recorded']
perform_symmetry_operation_upon_molecule_methods  = ['perform_symmetry_operation_upon_molecule']
make_crystal_associated_methods                   = determine_is_molecule_already_recorded_methods + perform_symmetry_operation_upon_molecule_methods

# Methods for adding Hydrogens and Methyl groups to a molecule.
add_hydrogens_to_molecules_methods                = ['add_hydrogens_to_molecules']
add_methyls_to_molecules_methods                  = ['add_methyls_to_molecules']
add_ethyls_to_molecules_methods                   = ['add_ethyls_to_molecules']
add_atoms_methods                                 = add_hydrogens_to_molecules_methods + add_methyls_to_molecules_methods + add_ethyls_to_molecules_methods

# Methods for obtaining ASE objects from chg files, and adding charges to molecules in other systems.
import_CHG_file_methods                           = ['import_CHG_file']
assigned_ATCs_to_molecules_invariance_methods     = ['assigned_ATCs_to_molecules_invariance_method']
ATC_and_charge_methods                            = import_CHG_file_methods + assigned_ATCs_to_molecules_invariance_methods

# Cell Generator methods
cell_generator_methods                            = ['Cell_Generator', 'convert_ijk_to_displacement_vector', 'convert_position_into_ijk_lengths']

# Calculators
calculator_methods                                = ['Coulomb']

# Methods related to the Electronic_Crystal_Calculation_Prep program that are used in other programs
ECCP_read_files_methods                            = ['what_ECCP_Information_files_do_we_have', 'run_which_ECCP_operations', 'get_ECCP_Information_data', 'check_ECCP_Information_details', 'get_crystal_file_from_ECCP_Information', 'get_dimer_details_data', 'get_equivalent_molecule_group_data', 'get_equivalent_dimer_group_data', 'check_consistancy_between_files']

# Utility methods
folder_methods                                    = ['make_folder', 'remove_folder', 'move_folder']
get_SameMoleculesDueToCrystalSymmetry_methods     = ['get_SameMoleculesDueToCrystalSymmetry']
utility_methods                                   = folder_methods + get_SameMoleculesDueToCrystalSymmetry_methods

# Methods to access
__all__ = making_and_processing_crystal_methods + molecules_inbetween_dimers_methods + general_methods + ase_associated_graph_methods + graph_matching_methods + make_crystal_associated_methods + add_atoms_methods + ATC_and_charge_methods + cell_generator_methods + calculator_methods + ECCP_read_files_methods + utility_methods

# ------------------------------------------------------------------------------------------------------------------------




