import multiprocessing as mp

from tqdm import tqdm

from SUMELF.SUMELF.general_methods.distance_methods                      import get_distance
from SUMELF.SUMELF.generators.Cell_Generator                             import Cell_Generator
from SUMELF.SUMELF.general_methods.convert_between_ijk_and_cell_position import convert_ijk_to_displacement_vector
from SUMELF.SUMELF.general_methods.general_molecules_methods             import is_the_same_molecule_exact

def get_molecules_in_vicinity_of_dimer(non_hydrogen_molecules, d_m1_elements, d_m2_elements, d_m1_positions, d_m2_positions, crystal_cell_lattice, no_of_cpus):
	"""
	This method is designed to obtain the external molecules that are invicity of the dimer.

	Parameters
	----------
	non_hydrogen_molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal. These should not contain hydrogens.
	d_m1_elements : list
		These are the elements of the atoms in the first molecule in the dimer. These do not contain hydrogen atoms.
	d_m2_elements : list
		These are the elements of the atoms in the second molecule in the dimer. These do not contain hydrogen atoms.
	d_m1_positions : np.array
		These are the positions of the atoms in the first molecule in the dimer. These should not contain the position of hydrogens.
	d_m2_positions : np.array
		These are the positions of the atoms in the second molecule in the dimer. These should not contain the position of hydrogens.
	crystal_cell_lattice : numpy.array
		This is the matrix of the unit cell for the crystal.

	Returns
	-------

	"""

	# First, obtain the furthest distance between atoms of the two molecules involved in the dimer.
	farthest_distance_between_molecules_in_dimer = get_farthest_distance_between_two_molecules(d_m1_positions, d_m2_positions)

	# Second, initialise the list to hold the external molecules that are in the vicinity of the dimer.
	external_molecules_within_vicinity = []

	# Third, for each molecule in the non_hydrogen_molecules list.
	#with mp.Manager() as manager: 
	#external_molecules_within_vicinity = manager.list()
	pool = mp.Pool(no_of_cpus)
	input_variables_generator = get_inputs(non_hydrogen_molecules, d_m1_elements, d_m2_elements, d_m1_positions, d_m2_positions, crystal_cell_lattice, farthest_distance_between_molecules_in_dimer, external_molecules_within_vicinity)
	#[is_molecule_in_the_vicinity_of_dimer(value) for value in input_variables_generator]
	external_molecules_within_vicinity = pool.map(is_molecule_in_the_vicinity_of_dimer, input_variables_generator)
	#pool.map_async(is_molecule_in_the_vicinity_of_dimer, tqdm(input_variables_generator, total=len(non_hydrogen_molecules), desc='determining molecules in vicinity of current dimer', unit='calc'))
	pool.close()
	pool.join()
	#external_molecules_within_vicinity = list(external_molecules_within_vicinity)

	raise Exception('Check that instances of index here is correct, or if it should be name. (name = index+1 for both molecules and dimers)')

	for index in range(len(external_molecules_within_vicinity)-1,-1,-1):
		if external_molecules_within_vicinity[index] is None:
			del external_molecules_within_vicinity[index]

	# Fourth, return external_molecules_within_vicinity
	return external_molecules_within_vicinity

def get_inputs(non_hydrogen_molecules, d_m1_elements, d_m2_elements, d_m1_positions, d_m2_positions, crystal_cell_lattice, farthest_distance_between_molecules_in_dimer, external_molecules_within_vicinity):
	"""
	This is a generator designed to provide the variables for input parameters.
	"""
	for index in range(1,len(non_hydrogen_molecules)+1):
		yield (index, non_hydrogen_molecules, d_m1_elements, d_m2_elements, d_m1_positions, d_m2_positions, crystal_cell_lattice, farthest_distance_between_molecules_in_dimer, external_molecules_within_vicinity)

def is_molecule_in_the_vicinity_of_dimer(input_variables):
	"""
	This method is designed to determine if an external molecule is in the invicity of a dimer.

	Parameters
	----------
	index : int
		This is the index of the molecule to investigate
	non_hydrogen_molecules : list of ase.Atoms objects.
		These are all the individual molecules identified in the crystal. These should not contain hydrogens.
	d_m1_elements : list
		These are the elements of the atoms in the first molecule in the dimer. These do not contain hydrogen atoms.
	d_m2_elements : list
		These are the elements of the atoms in the second molecule in the dimer. These do not contain hydrogen atoms.
	d_m1_positions : np.array
		These are the positions of the atoms in the first molecule in the dimer. These should not contain the position of hydrogens.
	d_m2_positions : np.array
		These are the positions of the atoms in the second molecule in the dimer. These should not contain the position of hydrogens.
	crystal_cell_lattice : numpy.array
		This is the matrix of the unit cell for the crystal.
	"""

	# First, obtain the variables from input_variables
	index, non_hydrogen_molecules, d_m1_elements, d_m2_elements, d_m1_positions, d_m2_positions, crystal_cell_lattice, farthest_distance_between_molecules_in_dimer, external_molecules_within_vicinity = input_variables

	# Second, obtain the molecule from the non_hydrogen_molecules list.
	non_hydrogen_molecule = non_hydrogen_molecules[index]

	# Third, create the cell generator object
	cell_generator = Cell_Generator()

	# Fourth, get the displacement of the molecule index2 from its original unit cell (being from indexB-indexA). 
	for unit_cell_displacement in cell_generator.generate_next_ijk_points(): 

		# Fifth, get the relative displacement of molecules index1 and index2.
		relative_displacement = convert_ijk_to_displacement_vector(unit_cell_displacement, crystal_cell_lattice)

		# Sixth, obtain the elements and the positions of the atoms in the external molecule. Here, the external molecule has been moved to the relative_displacement based on the cell to translate it into.
		external_molecule_elements  = non_hydrogen_molecule.get_chemical_symbols()
		external_molecule_positions = non_hydrogen_molecule.get_positions() + relative_displacement

		# Seventh, make sure that the external molecule is not the same as one of the molecules in the dimer.
		if is_the_same_molecule_exact(external_molecule_elements, d_m1_elements, external_molecule_positions, d_m1_positions) or is_the_same_molecule_exact(external_molecule_elements, d_m2_elements, external_molecule_positions, d_m2_positions):
			continue

		# Eighth, obtain the distance between the external molecule and each molecule in the dimer. 
		distance_between_molecule_1_and_external_molecule = get_shortest_distance_between_two_molecules(external_molecule_positions, d_m1_positions)
		distance_between_molecule_2_and_external_molecule = get_shortest_distance_between_two_molecules(external_molecule_positions, d_m2_positions)

		# Tenth, determine if this external molecule is in the vacinity of the dimer
		is_molecule_within_max_dimer_distance = (distance_between_molecule_1_and_external_molecule <= farthest_distance_between_molecules_in_dimer) and (distance_between_molecule_2_and_external_molecule <= farthest_distance_between_molecules_in_dimer)

		# Eleventh, send the result of if a dimer was generated from unit_cell_displacement to the cell_generator
		cell_generator.add_result(is_molecule_within_max_dimer_distance)

		# Tenth, if this external molecule is in the vacinity of the dimer, record it 
		if is_molecule_within_max_dimer_distance:
			new_non_hydrogen_molecule = non_hydrogen_molecule.copy()
			new_non_hydrogen_molecule.set_positions(external_molecule_positions)
			return (new_non_hydrogen_molecule, index, unit_cell_displacement)
		else:
			return None

# ====================================================================================================================================================================================

def get_farthest_distance_between_two_molecules(m1_positions, m2_positions):
	"""
	This method is designed to determine the furtherest distance of atoms between two molecules.

	Parameters
	----------
	m1_positions : np.array
		These are the positions of the atoms in the first molecule. 
	m2_positions : np.array
		These are the positions of the atoms in the second molecule. 

	Returns
	-------
	farthest_distance : float
		This is the farthest distance of atoms between two molecules.
	"""

	# First, initialise the farthest_distance variable
	farthest_distance = 0.0

	# Second, for each atomic position in molecule 1
	for m1_position in m1_positions:

		# Third, for each atomic position in molecule 2
		for m2_position in m2_positions:

			# Fourth, get the distance between the two selected atoms in each molecule
			distance = get_distance(m1_position, m2_position)

			# Fifth, record this new distance if it is greater than farthest_distance, otherwise move on.
			if distance > farthest_distance:
				farthest_distance = distance

	# Sixth, return farthest_distance
	return farthest_distance


def get_shortest_distance_between_two_molecules(m1_positions, m2_positions):
	"""
	This method is designed to determine the shortest distance of atoms between two molecules.

	Parameters
	----------
	m1_positions : np.array
		These are the positions of the atoms in the first molecule. 
	m2_positions : np.array
		These are the positions of the atoms in the second molecule. 

	Returns
	-------
	farthest_distance : float
		This is the shortest distance of atoms between two molecules.
	"""

	# First, initialise the shortest_distance variable
	shortest_distance = float('inf')

	# Second, for each atomic position in molecule 1
	for m1_position in m1_positions:

		# Third, for each atomic position in molecule 2
		for m2_position in m2_positions:

			# Fourth, get the distance between the two selected atoms in each molecule
			distance = get_distance(m1_position, m2_position)

			# Fifth, record this new distance if it is less than shortest_distance, otherwise move on.
			if distance < shortest_distance:
				shortest_distance = distance

	# Sixth, return shortest_distance
	return shortest_distance

# ====================================================================================================================================================================================



