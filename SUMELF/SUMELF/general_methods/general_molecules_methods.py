"""
general_molecules_methods.py, Geoffrey Weal, 19/7/22

This script is designed to hold methods useful for obtaining information about molecules in ASE.
"""
import numpy as np

from ase.io import read

from SUMELF.SUMELF.general_methods.general_data      import atomic_masses, atomic_numbers

from SUMELF.SUMELF.general_methods.unit_cell_methods import get_cell_corner_points
from SUMELF.SUMELF.general_methods.distance_methods  import get_distance
from SUMELF.SUMELF.general_methods.geometry_methods  import planeFit, project_point_onto_plane

# ---------------------------------------------------------------------------------------------------------------------------

def read_crystal(filepath):
    """
    This method is designed to read in crystal files.

    Parameters
    ----------
    filepath : str.
        This is the path to the crystal file.

    Returns
    -------
    crystal : ase.Atoms
        This is the crystal object.
    """

    # First, read the file from disk and import it as a ASE Atoms object
    if filepath.endswith('.cif'):
        #crystal = read(filepath,disorder_groups='remove_disorder')
        crystal = read(filepath)
    else:
        crystal = read(filepath)

    # Second, make sure that the periodic boundary condition setting is set to True for the crystal.
    crystal.set_pbc(True)

    # Third, return the crystal object
    return crystal

# ---------------------------------------------------------------------------------------------------------------------------

def get_centre_of_mass(elements, positions):
    """
    This method will give the centre of mass of the dimer/system.

    Parameters
    ----------
    elements : list of str.
        These are the elements of the atoms in the dimer/system.
    positions : 2D numpy.array
        These are the positions of the atoms in the dimer/system.

    Returns
    -------
    com_position : 1D numpy.array
        This is the centre of mass of the dimer/system.
    """

    # First, make sure that the number of elements is the same as the number of positions given.
    if not len(positions) == len(elements):
        raise Exception('huh')

    # Second, get the atomic masses of each atom in the molecule.
    atomic_masses_of_system = np.array([atomic_masses[atomic_numbers[element]] for element in elements])

    # Third, get the centre of mass position for the molecule.
    com_position = (atomic_masses_of_system @ positions) / float(sum(atomic_masses_of_system))

    # Fourth, return the centre of mass position for the molecule.
    return com_position

# ---------------------------------------------------------------------------------------------------------------------------

def get_centre_of_molecule(molecule, include_hydrogens=True):
    """
    Get the centre of the molecule.

    Parameters
    ----------
    molecule : ase.Atoms object
        This is the molecule you want to get the centre of molecule for.

    Returns
    -------
    centre_of_molecule : numpy.array
        This is the centre of the molecule (in xyz coordinates). 
        
    """

    # First, get the positions of atoms in the molecule.
    atom_positions = molecule.get_positions()

    # Second, remove any atoms that are hydrogens from the atom_positions list (if you dont want to include hydrogens)
    if not include_hydrogens:
        elements = atom_positions.get_chemical_symbols()
        rows_to_remove = [atom_index for atom_index in range(len(elements)-1,-1,-1) if (elements[atom_index] in ['H', 'D'])]
        atom_positions = np.delete(atom_positions, rows_to_remove, axis=1)

    # Third, get the centre of the molecule.
    centre_of_molecule = np.sum(atom_positions, axis=0)/float(len(atom_positions))

    # Fourth, return the centre of the molecule position.
    return centre_of_molecule

# ---------------------------------------------------------------------------------------------------------------------------

def get_number_of_lone_pairs_of_electron_pairs(hybridisation, total_no_of_neighbours, formal_charge, logger=None):
    """
    This method will return the bonding vectors for adding new hydrogens to the atom.

    Parameters
    ----------
    hybridisation : str.
        This is the hybridisation of the atom.
    total_no_of_neighbours : int
        This is the total number of atoms surrounding the central atom. 
    formal_charge : int
        This is the formal charge for this atom.

    Returns
    -------
    number_of_lone_pairs_of_electrons : list of numpy.array 
        This is the number number of lone pairs of electrons
    """

    # First, determine the number of sigma bonds based on the atoms hybridisation
    if   hybridisation == 'sp3':
        no_of_sigma_bonds = 4
    elif hybridisation == 'sp2':
        no_of_sigma_bonds = 3
    elif hybridisation == 'sp':
        no_of_sigma_bonds = 2
    else:
        import pdb; pdb.set_trace()
        raise Exception('Error: Only sp3, sp2, and sp hybridisations have been coded for currently. hybridisation given = '+str(hybridisation))

    # Second, determine the number of lone pairs of electrons required based on the hybridisation and the total number of atoms surrounding the molecule.
    number_of_lone_pairs_of_electrons = no_of_sigma_bonds - total_no_of_neighbours

    # Third, if the formal charge of the atom is not zero, change the number of number_of_lone_pairs_of_electrons
    if not formal_charge == 0:
        print('Write this part of code.')
        import pdb; pdb.set_trace()
        raise Exception('Dont contribute until this part of the code is written')

    # Fourth, if number_of_lone_pairs_of_electrons, the crystallographer may have placed to many neighbours around an atom in 
    # If so, give a warning to the user and set number_of_lone_pairs_of_electrons to 0. 
    if number_of_lone_pairs_of_electrons < 0:
        warning_message = 'Warning: The crystallographer has placed '+str(-number_of_lone_pairs_of_electrons)+' extra neighbours about an atom in this crystal. Check this crystal out to make sure everything is all good.'
        if logger is not None:
            logger.warning(warning_message)
        else:
            print(warning_message)
        number_of_lone_pairs_of_electrons = 0

    # Fifth, return number_of_lone_pairs_of_electron_pairs
    return number_of_lone_pairs_of_electrons

# ---------------------------------------------------------------------------------------------------------------------------

def get_hybridisation_from_ASE():
    """
    This method will return the hybridisation for several atoms based on their hybridisation.

    To do if needed.
    """
    raise Exception('Method to write.')

def get_hybridisation_from_CSD(atom_CSD):
    """
    This method will return the hybridisation for several atoms based on their hybridisation.

    This information is obtained using the Sybyl atom type.

    Sybyl atom type information is given at https://tccc.iesl.forth.gr/education/local/quantum/molecular_modeling/guide_documents/SYBYL_data_document.html

    Parameters
    ----------
    atom_CSD : ccdc.molecule.Atom
        This is the atom that we want to obtain the hybridsation of using the Sybyl atom type

    Return
    ------
    hybridisation :str.
        This is the hybridisation for this atom.
    """

    # First, obtain the Sybyl type for this atom as given by the CCDC. 
    sybyl_type = atom_CSD.sybyl_type

    # Second, obtain the hybridisation from the sybyl_type
    hybridisation = get_hybridisation(sybyl_type)

    # Third, return the hybridisation
    return hybridisation

def get_hybridisation(sybyl_type):
    """
    This method will return the hybridisation for several atoms based on their hybridisation.

    This information is obtained using the Sybyl atom type.

    Sybyl atom type information is given at https://tccc.iesl.forth.gr/education/local/quantum/molecular_modeling/guide_documents/SYBYL_data_document.html

    Parameters
    ----------
    sybyl_type : str
        This is the Sybyl type for the atom.

    Return
    ------
    hybridisation :str.
        This is the hybridisation for this atom.
    """

    # First, if the atom is a carbon, returns its hybridisation.
    if sybyl_type == 'C.1':
        return 'sp'
    if sybyl_type == 'C.2':
        return 'sp2'
    if sybyl_type == 'C.3':
        return 'sp3'
    if sybyl_type == 'C.ar': # Aromatic
        return 'sp2'
    if sybyl_type == 'C.cat': # Cation
        return 'sp2'

    # Second, if the atom is a nitrogen, returns its hybridisation.
    if sybyl_type == 'N.1':
        return 'sp'
    if sybyl_type == 'N.2':
        return 'sp2'
    if sybyl_type == 'N.3':
        return 'sp3'
    if sybyl_type == 'N.4':
        return 'sp3'
    if sybyl_type == 'N.am': # amide nitrogen
        return 'sp2'
    if sybyl_type == 'N.ar': # aromatic nitrogen
        return 'sp2'
    if sybyl_type == 'N.pl3': # trigonal nitrogen
        return 'sp2'

    # Third, if the atom is a phosphorus, returns its hybridisation.
    if sybyl_type == 'P.3':
        return 'sp3'

    # Fourth, if the atom is a oxygen, returns its hybridisation.
    if sybyl_type == 'O.2':
        return 'sp2'
    if sybyl_type == 'O.3':
        return 'sp3'
    if sybyl_type == 'O.co2':
        return 'sp3'
    if sybyl_type == 'O.spc':
        return 'sp3'
    if sybyl_type == 'O.t3p':
        raise Exception('Check this')
        return 'sp3'

    # Fifth, if the atom is a sulphur, returns its hybridisation.
    if sybyl_type == 'S.2':
        return 'sp2'
    if sybyl_type == 'S.3':
        return 'sp3'
    if sybyl_type == 'S.o':
        return 'sp3'
    if sybyl_type == 'S.o2':
        return 'sp3'

    # Sixth, if none of these options have been found, return dash line to indicate no hybridisation can be obtained
    return '-'

# ---------------------------------------------------------------------------------------------------------------------------

def get_bond_type_from_CSD(bond_CSD):
    """
    This method will return the bond type for several atoms based on their bond sybyl_type.

    This information is obtained using the Sybyl bond type.

    Sybyl bond type information is given at https://tccc.iesl.forth.gr/education/local/quantum/molecular_modeling/guide_documents/SYBYL_data_document.html

    Parameters
    ----------
    bond_CSD : ccdc.molecule.Bond
        This is the bond that we want to obtain the bond type of using the Sybyl bond type

    Return
    ------
    bond_type :str.
        This is the bond type for this bond.
    """

    # First, obtain the Sybyl type for this atom as given by the CCDC. 
    sybyl_type = bond_CSD.sybyl_type

    # Second, obtain the hybridisation from the sybyl_type
    bond_type = get_bond_type(sybyl_type)

    # Third, return the bond_type
    return bond_type

def get_bond_type(sybyl_type):
    """
    This method will return the bond type for several bond types

    This information is obtained using the Sybyl bond type.

    Sybyl bond type information is given at https://tccc.iesl.forth.gr/education/local/quantum/molecular_modeling/guide_documents/SYBYL_data_document.html

    Parameters
    ----------
    sybyl_type : str
        This is the Sybyl type for the bond.

    Return
    ------
    bond_type :str.
        This is the bond type for this bond.
    """

    # First, return sybyl_type into a more user friendly int/string 
    if sybyl_type == '1':
        return 1
    if sybyl_type == '2':
        return 2
    if sybyl_type == '3':
        return 3
    if sybyl_type == 'am': # Amide
        return 'amide'
    if sybyl_type == 'ar': # Aromatic
        return 'aromatic'
    if sybyl_type == 'du': # Dummy
        return 'dummy'
    if sybyl_type == 'nc': # Dummy
        return 'not_connected'
    if sybyl_type == 'un': # Dummy
        return 'unknown'

    # Seventh, if none of these options have been found, return dash line to indicate no bond type can be obtained
    return '-'

# ----------------------------------------------------------------------------------------------------------------
def get_atomic_rings_in_ase_object(rings, hybridisation_of_atoms, molecule):
    """
    Determine all the aromatic rings in the molecule.

    Requires the hybridisation to be correctly given in the crystal file.

    Parameters
    ----------
    rings : list of lists of int.
        These are all the indices of the rings in the ase.atoms object.
    hybridisation_of_atoms : list of str.
        These are the hybridisations of the atoms in the ase.atoms object.
    molecule : ase.Atoms
        This is the molecule of interest.

    Returns
    -------
    aromatic_rings : list of lists of int.
        These are the indices of the aromatic rings in the ase.atoms object.
    """

    # First, for each ring identified in the molecule
    aromatic_rings = []
    for ring in rings:

        # Second, obtain the hybridisations of atoms in rings
        hybridisation_of_atoms_in_ring = [hybridisation_of_atoms[ring_atom_index] for ring_atom_index in ring]

        # Third, determine if all the atoms in the ring are sp2 hybridised
        if all([(atom_hybridisation == 'sp2') for atom_hybridisation in hybridisation_of_atoms_in_ring]):
            aromatic_rings.append(ring)
            continue
    
    # Fourth, return the aromatic rings. 
    return aromatic_rings

# ----------------------------------------------------------------------------------------------------------------

def get_translation_to_move_COM_inside_unit_cell(new_molecule, crystal_cell_lattice):
    """
    This method will translate the compoent inside of the unit cell. 

    This will give a ase.Atoms object that will look more like CCDC Mercury displays its crystals

    Parameters
    ----------
    new_molecule : ase.Atoms
        This is the new molecule of the crystal to replicate.
    crystal_cell_lattice : np.array
        This is the matrix that includes the vector molecules of the unit cell.

    Returns
    -------
    add_translation : np.array
        This is a column vector that includes the translation vector required to move the molecule inside of the unit cell while maintaining the crystal periodicity.
    """

    # First, get the centre position of the origin unit cell.
    origin_cell_centre_of_unit_cell = get_the_centre_of_the_unit_cell(crystal_cell_lattice)

    # Second, get all the displacements that surround the cell around the molecule.
    cell_points = get_cell_corner_points(crystal_cell_lattice,super_cell_reach=1)

    # Third, determine the position of the molecule that is closest to the centre of the unit cell while maintaining the crystal periodicity.
    distance_from_closest_COM_to_MOC = float('inf')
    add_translation = None
    for cell_point in cell_points:

        # 3.1: Copy the new_molecule
        new_molecule_copy = new_molecule.copy()

        # 3.2: Place it in the translated position
        new_molecule_copy.set_positions(new_molecule_copy.get_positions() + cell_point)

        # 3.3: Obtain the centre of molecule (excluding hydrogens)
        if all(symbol in ['H','D'] for symbol in set(new_molecule_copy.get_chemical_symbols())):
            centre_of_molecule = [atom.position for atom in new_molecule_copy]
        else:
            centre_of_molecule = [atom.position for atom in new_molecule_copy if (not atom.symbol == 'H')]
        centre_of_molecule = sum(centre_of_molecule)/len(centre_of_molecule)

        # 3.4: Determine the distance from the centre of molecule to the centre of the unit cell
        current_distance_from_COM_to_MOC = get_distance(centre_of_molecule, origin_cell_centre_of_unit_cell)

        # 3.5: If you have found a shorter centre of molecule to centre of the unit cell distance, record that translation vector
        if current_distance_from_COM_to_MOC < distance_from_closest_COM_to_MOC:
            distance_from_closest_COM_to_MOC = current_distance_from_COM_to_MOC
            add_translation = cell_point

    # Fourth, return the translation vector that will move the molecule inside of the unit cell while maintaining the crystal periodicity.
    return add_translation

def get_the_centre_of_the_unit_cell(crystal_cell_lattice):
    """
    This method is designed to obtain the centre (middle) of the unit cell.

    Parameters
    ----------
    crystal_cell_lattice : np.array
        This is the matrix that includes the vector molecules of the unit cell.

    Returns
    -------
    origin_cell_centre_of_unit_cell : np.array
        This is the centre (middle) of the origin unit cell.
    """

    # First, get the 8 corner points of the origin unit cell
    origin_cell_cell_points = get_cell_corner_points(crystal_cell_lattice,super_cell_reach=1,bottom_of_range=0)

    # Second, get the middle of the origin unit cell.
    origin_cell_centre_of_unit_cell = sum(origin_cell_cell_points)/float(len(origin_cell_cell_points))

    # Third, return the centre (middle) of the origin unit cell.
    return origin_cell_centre_of_unit_cell

# -------------------------------------------------------------------------------------------------------------------------

distance_tolerance = 0.01
def is_the_same_molecule_exact(molecule1_elements, molecule2_elements, molecule1_positions, molecule2_positions):
    """
    This method is designed to indicate if two molecules are pretty much exactly the same, as in same atoms and overlapping in space.

    Parameters
    ----------

    """

    print('Finish off the documentation here.')
    import pdb; pdb.set_trace()
    raise Exception('This same method can be found in make_crystal_methods/same_molecules. Consider updating so there are no duplicate methods.')

    indices_of_molecule_2_to_investigate = list(range(len(molecule2_elements)))

    for molecule1_element, molecule1_position in zip(molecule1_elements, molecule1_positions):
        for index in indices_of_molecule_2_to_investigate:
            molecule2_element  = molecule2_elements[index]
            molecule2_position =  molecule2_positions[index]
            if (molecule1_element == molecule2_element) and (get_distance(molecule1_position, molecule2_position) <= distance_tolerance):
                indices_of_molecule_2_to_investigate.remove(index)
                break
        else:
            return False

    if not len(indices_of_molecule_2_to_investigate) == 0:
        raise Exception('huh?')

    return True

# ---------------------------------------------------------------------------------------------------------------------------

