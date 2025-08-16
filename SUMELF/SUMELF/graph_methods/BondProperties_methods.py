"""
neighbours_from_xyz_file.py, Geoffrey Weal, 4/8/22

This script is designed read and write the neighbours from the BondProperties of the xyz file.
"""
import numpy as np

def get_neighbours_from_BondProperties(BondProperties):
	"""
	This method is designed to convert the neighbours from the neighbourList of the xyz file.

	Parameters
	----------
	BondProperties : dict.
		This is a dict of the properties for each bond in the system. Format: bond_property_name -> bond: bond_property

	Returns
	-------
	bonds : list of tuple of ints
		This is the list of bonds
	bonds_as_key : dict.
		This contains all the propeties for each bond in BondProperties. Format: bond -> bond_property_name: bond_property
	"""

	# First, initialise the bonds list
	bonds        = None
	bonds_as_key = None

	# Second, for each property in BondProperties.
	for bond_property_name, bond_property_details in BondProperties.items():

		# 2.1: Get the bonds between atoms in the system from BondProperties.
		bonds_from_property = sorted([bond for (bond, bond_information) in bond_property_details])

		# 2.2: If bonds is None, record bonds as bonds_from_property
		if bonds is None:
			bonds        = bonds_from_property
			bonds_as_key = {bond: {} for bond in bonds}
		
		# 2.3: Check if the bonds recorded in bonds_from_property are the same as bonds
		#      * They should be. If not, there is a problem with the system ASE atoms object. 
		if not (bonds_from_property == bonds):
			raise Exception('Error here')

		# 2.4: Get the bond information, where the key is given as the bond indices
		for (bond, bond_information) in bond_property_details:
			bonds_as_key[bond][bond_property_name] = bond_information

	# Third, return bonds
	return bonds, bonds_as_key

def convert_BondProperties_for_ASE_Atoms_Object(bonds_as_key):
	"""
	This method is designed to convert the BondProperties format from:
		* bond -> bond_property_name: bond_property
			to 
		* bond_property_name -> bond: bond_property

	Parameters
	----------
	bonds_as_key : dict.
		This contains all the propeties for each bond in BondProperties. Format: bond -> bond_property_name: bond_property

	Returns
	-------
	BondProperties : dict.
		This is a dict of the properties for each bond in the system. Format: bond_property_name -> bond: bond_property
	"""

	# First, get all the properties given in the BondProperties
	all_bond_properties = set()
	for bonds, bond_properties in bonds_as_key.items():
		all_bond_properties.update(bond_properties.keys())

	# Second, make sure each bond contains all the bond properties obtained
	bond_property_check = {}
	for bond, bond_properties in bonds_as_key.items():
		differences = set(bond_properties.keys()) ^ set(all_bond_properties)
		if len(differences) > 0:
			bond_property_check[bond] = sorted(differences)
	if len(bond_property_check) > 0:
		raise Exception('Error: Bond bonds in bonds_as_key do not contain or contain extra bond properties: '+str(bond_property_check))

	# Third, initialise the bond_properties dictionary.
	BondProperties = {property_name: [] for property_name in all_bond_properties}

	# Fourth, convert bonds_as_key so that it goes bond_property_name -> bond: bond_property
	for bond, bond_properties in bonds_as_key.items(): 
		for bond_property_name, bond_property in bond_properties.items(): 
			BondProperties[bond_property_name].append((bond, bond_property))

	# Fifth, sort all the lists for each bond property name
	for property_name, property_info in BondProperties.items():
		BondProperties[property_name] = sorted(property_info)

	# Sixth, return BondProperties
	return BondProperties
