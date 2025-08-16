"""
unique_system_methods.py, Geoffrey Weal, 19/7/22

This script is designed to hold methods for determining if two system are unique or not.
"""

def obtain_unique_molecules(molecule_equivalence_method):
	"""
	Refer to obtain_unique_systems definition.
	"""
	return obtain_unique_systems(molecule_equivalence_method, 'molecule_equivalence_method')

def obtain_unique_dimers(dimer_equivalence_method):
	"""
	Refer to obtain_unique_systems definition.
	"""
	return obtain_unique_systems(dimer_equivalence_method, 'dimer_equivalence_method')

def obtain_unique_systems(equivalence_method, equivalence_method_name):
	"""
	This method will determine if a list of systems will be analysed to determine which are equivalent and which one are not

	Parameters
	----------
	equivalence_method : dict.
		This dictionary include information about the equivalence method that one wants to use
	equivalence_method_name : str.
		This is the name of the equivalent method that was used. This is used for debugging/programming breaking responces.

	Returns
	-------
	True if you want to locate unique systems, False if not.
	"""

	# First, if equivalence_method has not been given, return False
	if equivalence_method is None:
		return False

	# Second, if 'method' is not in the equivalence_method string, their may be an issue
	if not 'method' in equivalence_method:
		print('ECCP Issue: You have not given the method you want to use in the '+str(equivalence_method_name)+' dictionary.')
		print('Look at your '+str(equivalence_method_name)+' dictionary and run this program again.')
		exit('This program will finish without completing.')

	# Third, if equivalence_method is a dictinoary and is given as None, return False
	if equivalence_method['method'] in [None, 'None', 'none', 'NONE']:
		return False

	# Fourth, if nothing else, return True.
	return True

