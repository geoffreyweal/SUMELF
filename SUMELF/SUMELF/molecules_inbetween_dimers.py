"""
molecules_inbetween_dimers.py, Geoffrey Weal, 19/7/22

This script is designed to determine if there are any external molecules found inbetween an existing dimer.
"""
import os
import numpy as np
import multiprocessing as mp

from tqdm import tqdm

from SUMELF.SUMELF.general_methods.remove_hydrogens import remove_hydrogens
from SUMELF.SUMELF.molecules_inbetween_dimers_method.get_molecules_in_vicinity_of_dimer import get_molecules_in_vicinity_of_dimer
from SUMELF.SUMELF.general_methods.unit_cell_methods import number_of_points_in_cell
from scipy.spatial import Delaunay

from ase import Atoms
from ase.visualize import view

zero_point = np.array((0.,0.,0.))
def molecules_inbetween_dimers(all_dimers_info, molecules, molecule_graphs, crystal_cell_lattice, no_of_cpus=1, write_to_txt_filename='molecules_inbetween_dimers_information.txt'):
	"""
	This method is designed to determine if there are any external molecules found inbetween an existing dimer.

	Parameters
	----------
	all_dimers_info : list
		A list of dimers, consisting of the two molecules involved in the dimer.
	molecules : list of ase.Atoms
		These are the molecules that can be used to make dimers from.
	molecule_graphs : list of network.Graph
		These are the graphs associated with the molecules in the molecules list that can be used to make dimers from.
	crystal_cell_lattice : numpy.array
		This is the unit cell matrix, containing the unit cell vectors that describe the unit cell.

	Returns
	-------

	"""

	raise Exception('Need to check that this method works.')

	molecules_inbetween_dimers_information = os.getcwd()+'/'+write_to_txt_filename
	if os.path.exists(molecules_inbetween_dimers_information):
		os.remove(molecules_inbetween_dimers_information)
	with open(molecules_inbetween_dimers_information,'w') as informationTXTfile:
		informationTXTfile.write('')

	print('.---------------------------------------------------------.')
	print('Checking for molecules that are inbetween dimers')

	# First, obtain all molecules without hydrogen atoms.
	non_hydrogen_molecules           = {}
	non_hydrogen_molecules_elements  = {}
	non_hydrogen_molecules_positions = {}
	non_hydrogen_graphs              = {}
	for mol_name in range(1,len(molecules)+1):
		non_hydrogen_molecule, non_hydrogen_graph = remove_hydrogens(molecules[mol_name], graph=molecule_graphs[mol_name])
		non_hydrogen_molecules[mol_name]           = non_hydrogen_molecule
		non_hydrogen_molecules_elements[mol_name]  = non_hydrogen_molecule.get_chemical_symbols()
		non_hydrogen_molecules_positions[mol_name] = non_hydrogen_molecule.get_positions()
		non_hydrogen_graphs[mol_name]              = non_hydrogen_graph

	number_of_atoms_inbetween_dimers = []

	# Second, for each dimer obtained by the ECCP method
	pbar = tqdm(all_dimers_info)
	for dimer_no, d_m1_name, d_m2_name, cell_point, displacement_vector, move_com_by_1 in pbar:

		# Third, obtain the elements of the molecules involved in the dimer.
		d_m1_elements  = non_hydrogen_molecules_elements[d_m1_name]
		d_m2_elements  = non_hydrogen_molecules_elements[d_m2_name]

		# Fourth, obtain the positions of the molecules involved in the dimer.
		d_m1_positions = non_hydrogen_molecules_positions[d_m1_name] + move_com_by_1
		d_m2_positions = non_hydrogen_molecules_positions[d_m2_name] + move_com_by_1 + displacement_vector

		#dimer = Atoms(d_m1_elements+d_m2_elements,d_m1_positions.tolist() + d_m2_positions.tolist())

		# Fifth, obtain all the external molecules that lies within the vicinity of this dimer.
		molecules_in_vicinity_of_dimer = get_molecules_in_vicinity_of_dimer(non_hydrogen_molecules, d_m1_elements, d_m2_elements, d_m1_positions, d_m2_positions, crystal_cell_lattice, no_of_cpus=no_of_cpus)

		# Sixth, write what you are doing to pbar.
		pbar.set_description('Determine number of atoms inbetween dimer')

		# Seventh, obtain the Delaunay object, containing 
		dimer_object = Delaunay(np.concatenate((d_m1_positions,d_m2_positions)).tolist())

		# Eighth, obtain the number of atoms of each external molecule that lie inbetween this dimer.
		pool = mp.Pool(no_of_cpus)
		input_variables_generator = get_inputs(molecules_in_vicinity_of_dimer, cell_point, dimer_object)
		number_of_atoms_inbetween_dimer = pool.map(get_molecules_inbetween_dimers, input_variables_generator)
		pool.close()
		pool.join()

		for index in range(len(number_of_atoms_inbetween_dimer)-1,-1,-1):
			if number_of_atoms_inbetween_dimer[index] is None:
				del number_of_atoms_inbetween_dimer[index]

		with open(molecules_inbetween_dimers_information,'a') as informationTXTfile:
			informationTXTfile.write(str(dimer_no)+': '+str(number_of_atoms_inbetween_dimer)+'\n')

		# Ninth, record the number_of_atoms_inbetween_dimer to number_of_atoms_inbetween_dimers
		#number_of_atoms_inbetween_dimers.append(number_of_atoms_inbetween_dimer)

	# Tenth, return lists
	#return number_of_atoms_inbetween_dimers


def get_inputs(molecules_in_vicinity_of_dimer, unit_cell_displacement, dimer_object):
	for new_non_hydrogen_molecule, index, unit_cell_displacement in molecules_in_vicinity_of_dimer:
		yield (new_non_hydrogen_molecule, index, unit_cell_displacement, dimer_object)

def get_molecules_inbetween_dimers(input_variables):

	new_non_hydrogen_molecule, index, unit_cell_displacement, dimer_object = input_variables

	# 8.1: Obtain the positions of the external molecule.
	position = new_non_hydrogen_molecule.get_positions().tolist()

	# 8.2: Determine if their are any atoms in the Delaunay object.
	number_of_atoms_inbetween_dimer_for_external_molecule = number_of_points_in_cell(position, dimer_object, zero_point)

	# 8.3: If there are atoms in the Delaunay object, record these
	if number_of_atoms_inbetween_dimer_for_external_molecule > 0:
		return (number_of_atoms_inbetween_dimer_for_external_molecule, index, unit_cell_displacement)
	else:
		return None















