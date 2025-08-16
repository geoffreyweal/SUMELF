# Programs Available in SUMELF

The SUMELF program is designed to hold methods and modules that are used by multiple programs used in the grand scheme for simulating/calculating electronic properties as well as charge and exciton diffusion in molecular crystals. 

## The ``get_molecules`` module

It is possible to separate a crystal into its separate molecules using the ``get_molecules`` module. This module will extract the individual molecules in a crystal structure, and save them into separate ``.xyz`` files. 

To use the ``get_molecules`` module, type the following into your terminal:

```bash
# First, change directory into the directory contain the folder with your crystal files
cd path_to_crystal_database

# Second, run the get_molecules program:
SUMELF get_molecules crystal_database
```

This will separate the crystal into its individual molecules, and save them as individual ``xyz`` file in folder called ``crystal_database_molecules``. You can view these molecules using your favourite GUI viewer, mine is ``ase gui`` ([click here for more about ``ase gui``](https://wiki.fysik.dtu.dk/ase/ase/gui/basics.html)).

NEED TO ADD OPTIONS TO THIS.

!!! info

	You can also use [Mercury (a CCDC program)](https://www.ccdc.cam.ac.uk/solutions/software/free-mercury/) to view the crystal. However, you will need to use the ``get_molecules`` to get the names of the molecules and the atom indices to make modifications to for the  ``Repair_Crystals`` module in the ReCrystals program

