"""
get_symmetry_operations.py, Geoffrey Weal, 18/6/22

This script is designed to return the matrices for performing symmetry operations
"""

import numpy as np

def parse_sitesym_element(element):
    """Parses one element from a single site symmetry in the form used
    by the International Tables.

    From https://gitlab.com/ase/ase/-/blob/master/ase/spacegroup/spacegroup.py
    
    Examples:
    
    >>> parse_sitesym_element("x")
    ([(0, 1)], 0.0)
    >>> parse_sitesym_element("-1/2-y")
    ([(1, -1)], -0.5)
    >>> parse_sitesym_element("z+0.25")
    ([(2, 1)], 0.25)
    >>> parse_sitesym_element("x-z+0.5")
    ([(0, 1), (2, -1)], 0.5)
    
    Parameters
    ----------
    element: str
        Site symmetry like "x" or "-y+1/4" or "0.5+z".
      
    Returns
    -------
    list[tuple[int, int]]
        Rotation information in the form '(index, sign)' where index is
        0 for "x", 1 for "y" and 2 for "z" and sign is '1' for a positive
        entry and '-1' for a negative entry. E.g. "x" is '(0, 1)' and
        "-z" is (2, -1).
      
    float
        Translation information in fractional space. E.g. "-1/4" is
        '-0.25' and "1/2" is '0.5' and "0.75" is '0.75'.
    """
    element = element.lower()
    is_positive = True
    is_frac = False
    sng_trans = None
    fst_trans = []
    snd_trans = []
    rot = []
    
    for char in element:
        if char == "+":
            is_positive = True
        elif char == "-":
            is_positive = False
        elif char == "/":
            is_frac = True
        elif char in "xyz":
            rot.append((ord(char)-ord("x"), 1 if is_positive else -1))
        elif char.isdigit() or char == ".":
            if sng_trans is None:
                sng_trans = 1.0 if is_positive else -1.0
            if is_frac:
                snd_trans.append(char)
            else:
                fst_trans.append(char)
    
    trans = 0.0 if not fst_trans else (sng_trans * float("".join(fst_trans)))
    if is_frac:
        trans /= float("".join(snd_trans))
        
    return rot, trans

def parse_sitesym_single(sym, sep=",", force_positive_translation=False):
    """Parses a single site symmetry in the form used by International 
    Tables and overwrites 'out_rot' and 'out_trans' with data.

    From https://gitlab.com/ase/ase/-/blob/master/ase/spacegroup/spacegroup.py
    
    Parameters
    ----------  
    sym: str
        Site symmetry in the form used by International Tables (e.g. "x,y,z", "y-1/2,x,-z").
    sep: str
        String separator ("," in "x,y,z").
    force_positive_translation: bool
        Forces fractional translations to be between 0 and 1 (otherwise negative values might be accepted).
        Defaults to 'False'.
      
    Returns
    -------
    Nothing is returned: 'out_rot' and 'out_trans' are changed inplace.
    """

    out_rot = np.zeros((3, 3), dtype='int')
    out_trans = np.zeros(3)

    out_rot[:] = 0.0
    out_trans[:] = 0.0
    
    for i, element in enumerate(sym.split(sep)):
        e_rot_list, e_trans = parse_sitesym_element(element)
        for rot_idx, rot_sgn in e_rot_list:
            out_rot[i][rot_idx] = rot_sgn
        out_trans[i] = (e_trans % 1.0) if force_positive_translation else e_trans

    return out_rot, out_trans
        
def get_symmetry_operations(symmetry_operators):
    """
    Obtain the symmetry operations from the given symmetry_operator string.

    Parameters
    ----------
    symmetry_operators : list of str
        These are the symmetry operation obtained from the crystal object from the CSD Python API.

    Returns
    -------
    symmetry_operations : list of (numpy.array, numpy.array)
        These are all the symmetry operations as numpy array objects
    """
    return [parse_sitesym_single(symmetry_operator) for symmetry_operator in symmetry_operators]

# --------------------------------------------------------------------------------------------------------------

def get_symmetry_operator(symmetry_operation):
    """
    This method is designed to convert the rotation and translation numpy.arrays of a symmetry_operation into a string that symbolises this operation.

    Parameters
    ----------
    symmetry_operation : (numpy.array, numpy.array)
        This is the rotation and translation numpy.arrays foir a symmetry operation

    Returns
    -------
    symmetry_operators : str.
        This is the symmetry operation in string form.
    """

    # First, obtain the rotation and translation numpy.array for this operation
    rotation, translation = symmetry_operation

    # Second, convert the operations upon the x, y and z axis into a string 
    axis_strings = []
    for axis_index in range(3):

        # 2.1: Get the row from the rotation and tralation matrix for the axis of interest.
        line_rotation = rotation[axis_index,:]
        translation_value = translation[axis_index]

        # 2.2: Convert the rotation conponent of the symmetry operation into a string
        symmetry_string = ''
        for value, axis in zip(line_rotation, ['x', 'y', 'z']):
            if value == 0:
                continue
            elif value == 1:
                value = '+'
            elif value == -1:
                value = '-'
            else:
                value = str(value)
            symmetry_string += value + axis

        # 2.3: Convert the translation conponent of the symmetry operation into a string
        if translation_value > 0.0:
            symmetry_string += '+'+str(translation_value)
        elif translation_value < 0.0:
            symmetry_string += str(translation_value)

        # 2.4: If nothing was added to the string, set it to 0
        if symmetry_string == '':
            symmetry_string = '0'

        # 2.5: If symmetry_string starts with +, remove the +
        if symmetry_string.startswith('+'):
            symmetry_string = symmetry_string[1:]
        
        # 2.5: Return the symmetry operation for this axis as a string
        axis_strings.append(symmetry_string)

    # Third, convert the axis_strings into a single full string.
    symmetry_operator = ','.join(axis_strings)

    # Fourth, return the symmetry operator as a string.
    return symmetry_operator

def get_symmetry_operators_string(symmetry_operations):
    """
    This method is designed to convert the rotation and translation numpy.arrays in the symmetry_operations list into strings.

    Parameters
    ----------
    symmetry_operations : list of (numpy.array, numpy.array)
        These are all the symmetry operations as numpy array objects

    Returns
    -------
    symmetry_operators : list of str.
        These are the symmetry operation in string form.
    """
    return tuple([get_symmetry_operator(symmetry_operation) for symmetry_operation in symmetry_operations])

# --------------------------------------------------------------------------------------------------------------



























