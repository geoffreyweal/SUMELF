# The Supporting Methods for Electronic Functions (SUMELF) Program

Authors: Dr. Geoffrey Weal<sup>\*,†</sup>, Dr. Chayanit Wechwithayakhlung<sup>†</sup>, Dr. Josh Sutton<sup>\*</sup>, Dr. Daniel Packwood<sup>†</sup>, Dr. Paul Hume<sup>\*</sup>, Prof. Justin Hodgkiss<sup>\*</sup>

<sup>\*</sup> Victoria University of Wellington, Wellington, New Zealand; The MacDiarmid Institute for Advanced Materials and Nanotechnology, Wellington, New Zealand. 

<sup>†</sup> Institute for Integrated Cell-Material Sciences (iCeMS), Kyoto University, Kyoto, Japan.

Group pages: https://people.wgtn.ac.nz/paul.hume/grants, https://www.packwood.icems.kyoto-u.ac.jp/, https://people.wgtn.ac.nz/justin.hodgkiss/grants


## What is the Supporting Methods for Electronic Functions (SUMELF) Program

The Supporting Methods for Electronic Functions (SUMELF) program contains various supporting program that are used by various programs. This includes:


**ACSD: The Access Cambridge Structural Database Program**

* Used to collect crystals directly from the CCDC. 


**ReCrystals: The Repair Crystals Program**

* Used to repair and modify the molecules in the crystal.


**RSGC: The Remove SideGroups from Crystals Program**

* Used to remove sidegroups from the crystals, in particular aliphatic sidechains that are not significant to semiconductor properties and can make electronic calculations difficult. 


**ECCP: The Electronic Crystal Calculation Prep Program**

* Used to extract the molecules from the crystal and collect electronic calculation data for understanding coupling, charge-transfer, and energy-trasnfer for exciton and charge diffusion.


**EKMC: The Exciton Kinetic Monte Carlo Program**

* Used to simulate exciton diffusion using kinetic Monte Carlo.


**SORE: The Sum-Over-Rates Equation Program**

* Used to simulate exciton diffusion using a simple sum-over-rate equation model. 


## Installation

It is recommended to read the installation page before using the SUMELF program. See [Installation: Setting Up SUMELF and Pre-Requisites Packages](Installation.md#Installation: Setting Up SUMELF and Pre-Requisites Packages) for more information. Note that you can install SUMELF through ``pip3`` and ``conda``. 

## Guide To Using SUMELF

The SUMELF program contains methods used for multiple programs, so you don't use SUMELF directly. However, SUMELF does contain several modules that can be helpful across these programs. See [How To Use The SUMELF Program](Using_The_SUMELF_Program.md) to learn about how to use these modules. 

## The Grand Scheme

The SUMELF program is used as part of a grand scheme for calculating the excited-state electronic properties of molecules in a crystal. This includes simulations of exciton and charge diffusion through crystal structures, in particular for organic molecules (but not limited to them). This scheme is shown below, along with where the SUMELF program is used in this scheme. 

<img alt="Schematic of Grand Scheme" src="Shared_Images/Grand_Scheme/Grand_Scheme.png?raw=true#only-light" />
<img alt="Schematic of Grand Scheme" src="Shared_Images/Grand_Scheme/Grand_Scheme_Dark.png?raw=true#only-dark" />
