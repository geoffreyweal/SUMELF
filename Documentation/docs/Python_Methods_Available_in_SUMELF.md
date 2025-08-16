# Python Methods Available in SUMELF

The SUMELF module contains variable methods that are used by programs, including the:

* the Access Cambridge Structural Database (ACSD) program, 
* the Remove Sidechain Groups from Crystals (RSGC) program, 
* the Electronic Crystal Calculation Prep (ECCP) program, and
* the Excitonic Kinetic Monte Carlo (EKMC) program.

This page include information about each of the methods included in this module.

If you would like to use any of the methods given in this page, you can use them by adding at the top of your python script:

```python
from SUMELF import XYZ
```

where ``XYZ`` is the method you would like to use. For example, if you would like to use the  ``process_crystal`` method, write into your python script:

```python
from SUMELF import process_crystal
```

If you would like to import more than one method into your program, for example, you want to import the ``process_crystal``, ``above_angle_tolerance``, and ``centre_molecule_in_cell`` methods, write at the top of your python script:

```python
from SUMELF import process_crystal, above_angle_tolerance, centre_molecule_in_cell
```


## process_crystal

.. autoclass:: SUMELF.SUMELF.process_crystal.process_crystal
   :members:


## make_crystal

.. autoclass:: SUMELF.SUMELF.make_crystal.make_crystal
   :members:


## make_crystal_from_molecules

.. autoclass:: SUMELF.SUMELF.make_crystal_from_molecules.make_crystal_from_molecules
   :members:


## general_methods

There are numerous methods that cover many general methods that may be useful.


### angle_methods

.. automodule:: SUMELF.SUMELF.general_methods.angle_methods
   :members:


### general_molecules_methods

.. automodule:: SUMELF.SUMELF.general_methods.general_molecules_methods
   :members:


### get_symmetry_operations

.. automodule:: SUMELF.SUMELF.general_methods.get_symmetry_operations
   :members:


### ideal_bond_lengths_and_angles_methods

.. automodule:: SUMELF.SUMELF.general_methods.ideal_bond_lengths_and_angles_methods
   :members:


### make_dimer

.. automodule:: SUMELF.SUMELF.general_methods.make_dimer
   :members:


### NeighboursList_methods

.. automodule:: SUMELF.SUMELF.general_methods.NeighboursList_methods
   :members:


### obtain_graph

.. automodule:: SUMELF.SUMELF.general_methods.obtain_graph
   :members:


### process_method_in_parallel

.. automodule:: SUMELF.SUMELF.general_methods.process_method_in_parallel
   :members:


### remove_hydrogens

.. automodule:: SUMELF.SUMELF.general_methods.remove_hydrogens
   :members:


### unique_system_methods

.. automodule:: SUMELF.SUMELF.general_methods.unique_system_methods
   :members:


### unit_cell_methods

.. automodule:: SUMELF.SUMELF.general_methods.unit_cell_methods
   :members:





