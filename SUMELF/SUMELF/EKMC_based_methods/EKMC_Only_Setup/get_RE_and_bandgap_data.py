"""
get_EET_coupling_data.py, Geoffrey Weal, 19/4/22

This script will obtain the electronic coupling energies as calculated using the eet function in Gaussian between pairs of neighbouring molecules (dimers).
"""
import os
from ase.io import read
from SUMELF.SUMELF.EKMC_based_methods.EKMC_Only_Setup.get_RE_and_bandgap_data_methods.helpful_functions import get_reorganisation_energy, convert_hartree_to_eV

def get_RE_and_bandgap_data(reorganisation_and_bandgap_energy_details, crystal_name, functional_and_basis_set, molecules_path, original_molecule_names):
	"""
	This method will obtain the ground and excited sstructure energy gaps requierd for obtaining reorganisation energies between molecules in crystal.

	Parameters
	----------
	reorganisation_and_bandgap_energy_details : str.
		These are the details about the reorganisation energies between dimers and the bandgap energies of molecules data.
	crystal_name : str.
		This is the name of the crystal.
	functional_and_basis_set : str.
		This is the name of the functional and basis set used in calculations to obtain EET and ATC information. 
	molecules_path : str.
		This is the path to all the molecules in the crystal.
	original_molecule_names : list of ints
		These are the names of the molecules in the crystal. 

	Returns
	------- 
	reorganisation_energy_data : dict.
		These are all the reorgansation energies for the dimers in the crystal.
	bandgap_energy_data : dict.
		These are all the bandgap energies for the molecules in the crystal.
	conformationally_equivalent_molecules : dict. 
		This dictionary contains informatino about which conformationally equivalent molecules are assigned to which conformatinoally unique molecules (and theirfore have the same reorganisation energy data assigned to them)
	"""

	# Preamble, remove solvent data from molecule names in the original_molecule_names list.
	molecule_names = [int(str(molecule_name).replace('S','')) for molecule_name in original_molecule_names]

	# First, collect the information about which conformationally equivalent molecules are the same as the conformationally unique molecules.
	conformationally_equivalent_molecules = {}
	with open(molecules_path+'/Conformationally_Unique_Molecule_Information.txt','r') as Conformationally_Unique_Molecule_InformationTXT:
		Conformationally_Unique_Molecule_InformationTXT.readline()
		for line in Conformationally_Unique_Molecule_InformationTXT:
			line = line.rstrip().split()
			symmetric_molecule_no = int(line[0].replace('S',''))
			unique_molecule_no    = int(line[2].replace('S',''))
			if (symmetric_molecule_no in molecule_names) and (unique_molecule_no in molecule_names):
				conformationally_equivalent_molecules[symmetric_molecule_no] = unique_molecule_no

	# Second, quick check to make sure Conformationally_Unique_Molecule_Information.txt is all good.
	# Values are unique molecules, Keys are equivalent molecules
	if len(set(conformationally_equivalent_molecules.keys()) & set(conformationally_equivalent_molecules.values())) > 0:
		raise Exception('Error: there are equivalent molecules that have also been assigned as unique in Conformationally_Unique_Molecule_Information.txt.\nCheck this out\nconformationally_equivalent_molecules = '+str(conformationally_equivalent_molecules))

	# Third, set up the reorganisation_energy_data and the bandgap_energy_data dictionaries. 
	reorganisation_energy_data = {}
	bandgap_energy_data        = {}

	# ==================================================================================================================================

	# Fourth, if reorganisation_and_bandgap_energy_details['path_to_RE_folder'] exists, gather reorganisation energy and bandgap energy from reorganisation_and_bandgap_energy_details
	if 'path_to_RE_folder' in reorganisation_and_bandgap_energy_details:

		# 4.1: Get the path to the reorganisation energy folder.
		reorganisation_energy_path = reorganisation_and_bandgap_energy_details['path_to_RE_folder']

		# 4.2: Check if the relavant Egaps.txt file exists in the Individual_RE_Data folder of either the Unique_RE_Gaussian_Jobs or All_RE_Gaussian_Jobs folder. 
		if not os.path.exists(reorganisation_energy_path):
			raise Exception('Error: Can not find the reorganisation energy data folder: '+str(reorganisation_energy_path))
		if not os.path.exists(reorganisation_energy_path+'/Individual_RE_Data'):
			raise Exception('Error: Can not find the reorganisation energy data folder: '+reorganisation_energy_path+'/Individual_RE_Data')
		energies_TXT_filename = crystal_name+'_'+functional_and_basis_set+'_energies.txt'
		if not os.path.exists(reorganisation_energy_path+'/Individual_RE_Data/'+energies_TXT_filename):
			raise Exception('Error: Can not find the energies.txt file (contains energy gap data for molecule in crystal): '+reorganisation_energy_path+'/Individual_RE_Data/'+energies_TXT_filename)
		if not os.path.isfile(reorganisation_energy_path+'/Individual_RE_Data/'+energies_TXT_filename):
			raise Exception('Error: energies.txt (contains energy gap data for molecule in crystal) is not a file (maybe a folder?): '+reorganisation_energy_path+'/Individual_RE_Data/'+energies_TXT_filename)

		# 4.3: Record the energy for ground and excited structures for each molecule in the crystal.
		reorganisation_energy_molecule_data = {}
		with open(reorganisation_energy_path+'/Individual_RE_Data/'+energies_TXT_filename,'r') as EgapTXT:
			EgapTXT.readline()
			for line in EgapTXT:
				molecule_name, eGS_gGS_energy, eGS_gES_energy, eES_gGS_energy, eES_gES_energy, band_gap_energy = line.rstrip().split()
				molecule_no = int(molecule_name.replace('molecule_','').replace('S',''))
				if molecule_no in reorganisation_energy_molecule_data:
					raise Exception('Error: Molecule'+str(molecule_no)+' seems to appear twice in '+reorganisation_energy_path+'/Individual_RE_Data/'+energies_TXT_filename)
				reorganisation_energy_molecule_data[molecule_no] = (float(eGS_gGS_energy), float(eGS_gES_energy), float(eES_gGS_energy), float(eES_gES_energy))

		# 4.4: Obtain all the molecules that we have reorganisation energy data for on disk.
		molecules_in_reorganisation_energy_molecule_data = sorted(reorganisation_energy_molecule_data.keys())

		# 4.5: Record the reorganisation energies of the dimers in the crystal.
		for index1 in range(len(molecules_in_reorganisation_energy_molecule_data)):

			# 4.5.1: Get the name of molecule 1.
			mol1_name = molecules_in_reorganisation_energy_molecule_data[index1]

			# 4.5.2: Check if we want to obtain the reorganisation energy of mol1_name (we may not want to if this is a solvent).
			if mol1_name not in molecule_names:
				continue
			
			# 4.5.3: Check if mol1_name is not an equivalent molecule.
			if mol1_name in conformationally_equivalent_molecules.keys():
				raise Exception('Error: molecule '+str(mol1_name)+' is an equivalent molecule, so we dont want to include its reorganisation energy data from the disk.')
			
			# 4.5.4: Get the energies of the excited and ground state of mol1_name 
			mol1_eGS_gGS, mol1_eGS_gES, mol1_eES_gGS, mol1_eES_gES = reorganisation_energy_molecule_data[mol1_name]
			
			for index2 in range(index1, len(molecules_in_reorganisation_energy_molecule_data)):

				# 4.5.5: Get the name of molecule 2.
				mol2_name = molecules_in_reorganisation_energy_molecule_data[index2]

				# 4.5.6: Check if we want to obtain the reorganisation energy of mol6_name (we may not want to if this is a solvent).
				if mol2_name not in molecule_names:
					continue

				# 4.5.7: Check if mol2_name is not an equivalent molecule.
				if mol2_name in conformationally_equivalent_molecules.keys():
					raise Exception('Error: molecule '+str(mol2_name)+' is an equivalent molecule, so we dont want to include its reorganisation energy data from the disk.')

				# 4.5.8: Get the energies of the excited and ground state of mol2_name 
				mol2_eGS_gGS, mol2_eGS_gES, mol2_eES_gGS, mol2_eES_gES = reorganisation_energy_molecule_data[mol2_name]

				# 4.5.9: Get the dimer with arrangement (mol1_name, mol2_name), check it is not already in reorganisation_energy_data, and record the reorganisation energy for a molecule moving from mol1 --> mol2.
				dimer_pair1 = (mol1_name, mol2_name)
				if dimer_pair1 in reorganisation_energy_data:
					raise Exception('Error: dimer_pair1 '+str(dimer_pair1)+' is being entered into reorganisation_energy_data twice?\nmolecules_in_reorganisation_energy_molecule_data_dict = '+str(molecules_in_reorganisation_energy_molecule_data))
				reorganisation_energy_data[dimer_pair1] = convert_hartree_to_eV(get_reorganisation_energy(mol1_eGS_gGS, mol1_eGS_gES, mol2_eES_gGS, mol2_eES_gES))

				# 4.5.10: If mol1_name == mol2_name, dont need to do the next part of this for loop.
				if mol1_name == mol2_name:
					continue

				# 4.5.11: Get the dimer with arrangement (mol2_name, mol1_name), check it is not already in reorganisation_energy_data, and record the reorganisation energy for a molecule moving from mol2 --> mol1.
				dimer_pair2 = (mol2_name, mol1_name)
				if dimer_pair2 in reorganisation_energy_data:
					raise Exception('Error: dimer_pair2 '+str(dimer_pair2)+' is being entered into reorganisation_energy_data twice?\nmolecules_in_reorganisation_energy_molecule_data_dict = '+str(molecules_in_reorganisation_energy_molecule_data))
				reorganisation_energy_data[dimer_pair2] = convert_hartree_to_eV(get_reorganisation_energy(mol2_eGS_gGS, mol2_eGS_gES, mol1_eES_gGS, mol1_eES_gES))

		# 4.6: Record the bandgap energies of the molecules in the crystal.
		for mol_name in molecules_in_reorganisation_energy_molecule_data:
			eGS_gGS, eGS_gES, eES_gGS, eES_gES = reorganisation_energy_molecule_data[mol_name]
			if mol_name in bandgap_energy_data:
				raise Exception('Error: mol_name '+str(mol_name)+' is being entered into bandgap_energy_data twice?\bandgap_energy_data = '+str(bandgap_energy_data))
			if mol_name in molecule_names:
				bandgap_energy_data[mol_name] = convert_hartree_to_eV(eES_gES-eGS_gGS)

	# ==================================================================================================================================

	# Fifth, if reorganisation or bandgap energy data is given manually, record unique/equivalent molecules data
	if ('manual_RE_data' in reorganisation_and_bandgap_energy_details) or ('manual_bandgap_data' in reorganisation_and_bandgap_energy_details):

		# 5.1: Obtain the list of equivalent and unique molecules in the crystal. 
		equivalent_molecules = sorted(set(conformationally_equivalent_molecules.keys()))
		unique_molecules = []
		for molecule_name in molecule_names:
			if molecule_name not in equivalent_molecules:
				unique_molecules.append(molecule_name)
		unique_molecules = sorted(set(unique_molecules))

	# ==================================================================================================================================

	# Sixth, if reorganisation energy data is given manually, record this and do consistancy checks. 
	if 'manual_RE_data' in reorganisation_and_bandgap_energy_details:

		# 6.1: Record the manually entered reorganisation energy data.
		manual_RE_data = {}
		for key, value in reorganisation_and_bandgap_energy_details['manual_RE_data'].items():
			if isinstance(key,tuple):
				molecule1_name = int(str(key[0]).replace('S',''))
				molecule2_name = int(str(key[1]).replace('S',''))
				dimer_pair = (molecule1_name, molecule2_name)
				manual_RE_data[dimer_pair] = value
		for key, value in reorganisation_and_bandgap_energy_details['manual_RE_data'].items():
			if not isinstance(key,tuple):
				molecule_name = int(str(key).replace('S',''))
				dimer_pair = (molecule_name, molecule_name)
				if dimer_pair in manual_RE_data:
					raise Exception('Error: Reorganisation energy for molecule '+str(molecule_name)+' is already entered in reorganisation_and_bandgap_energy_details["manual_RE_data"] as '+str(dimer_pair)+'.\nCheck this out.\nreorganisation_and_bandgap_energy_details["manual_RE_data"] = '+str(reorganisation_and_bandgap_energy_details["manual_RE_data"]))
				manual_RE_data[dimer_pair] = value

		# 6.2: Check to make sure that a molecule entered into reorganisation_and_bandgap_energy_details["manual_RE_data"]
		#      actually exists in the crystal.
		inconsistant_molecules_REs = []
		for (molecule1_name, molecule2_name) in manual_RE_data.keys():
			if molecule1_name not in unique_molecules+equivalent_molecules:
				inconsistant_molecules_REs.append(molecule1_name)
			if molecule2_name not in unique_molecules+equivalent_molecules:
				inconsistant_molecules_REs.append(molecule2_name)
		if len(inconsistant_molecules_REs) > 0:
			inconsistant_molecules_REs = sorted(set(inconsistant_molecules_REs))
			toString  = 'Some of the molecules entered into reorganisation_and_bandgap_energy_details["manual_RE_data"] were not found in the crystal.'+'\n'
			toString += 'Inconsistent molecules in "manual_RE_data": '+str(inconsistant_molecules_REs)+'\n'
			toString += 'Possible dimers between molecules '+str(sorted(unique_molecules+equivalent_molecules))+'\n'
			toString += 'Those unique molecules: '+str(unique_molecules)
			raise Exception(toString)

		# 6.3: Add manually entered reorganisation energies of unique molecules to reorganisation_energy_data
		inconsistant_dimers_REs = []
		for mol1_name in sorted(molecule_names):
			for mol2_name in sorted(molecule_names):

				# 6.3.1: Obtain the reorganisation energies for the unique molecules
				if (mol1_name in unique_molecules) and (mol2_name in unique_molecules):

					# 6.3.1.1: If both mol1_name and mol2_name are unique, record the reorganisation energy
					dimer = (mol1_name, mol2_name)
					reorganisation_energy_data[dimer] = manual_RE_data[dimer]

				else:

					# 6.3.1.2: If either mol1_name or mol2_name are equivalent, check that it is not inconsistent with the unique molecules. 
					unique_mol1_unique = conformationally_equivalent_molecules[mol1_name] if (mol1_name in conformationally_equivalent_molecules.keys()) else mol1_name
					unique_mol2_unique = conformationally_equivalent_molecules[mol2_name] if (mol2_name in conformationally_equivalent_molecules.keys()) else mol2_name

					# 6.3.1.2: Get the equivalent dimer and its associated unique dimer.
					equivalent_dimer = (mol1_name, mol2_name)
					unique_dimer     = (unique_mol1_unique, unique_mol2_unique)

					# 6.3.1.3: If the user doesnt give a manual reorganisation energy for equivalent_dimer, give the computational version for manual_RE_data.
					if equivalent_dimer not in manual_RE_data:
						continue

					# 6.3.1.3: If this issue arises, then there is a missing unique dimer in reorganisation_and_bandgap_energy_details['manual_RE_data'].
					if unique_dimer not in manual_RE_data:
						toString  = 'Error: You have not included some unique dimers in your reorganisation_and_bandgap_energy_details["manual_RE_data"] dictionary.\n'
						toString += 'Check this.\n'
						toString += 'reorganisation_and_bandgap_energy_details["manual_RE_data"] = '+str(reorganisation_and_bandgap_energy_details["manual_RE_data"])+'\n'
						toString += 'Problematic dimer: '+str(unique_dimer)+'\n'
						toString += 'unique_molecules = '+str(unique_molecules)+'\n'
						toString += 'equivalent_molecules = '+str(equivalent_molecules)
						raise Exception(toString)

					# 6.3.1.4: Obtain the reorganisation energy for the equivalent and unique dimer
					reorganisation_energy_equivalent = manual_RE_data[equivalent_dimer]
					reorganisation_energy_unique     = manual_RE_data[unique_dimer]

					# 6.3.1.5: If the reorganisation energies for the equivalent and unique dimer are not the same, generate error
					if not (reorganisation_energy_equivalent == reorganisation_energy_unique):
						inconsistant_dimers_REs.append((equivalent_dimer, unique_dimer))

		# 6.4: Indicate if there were any inconsistancies between manually entered reorganisation energies for equivalent molecules.
		if len(inconsistant_dimers_REs) > 0:
			toString  = 'Error: There seems to be inconsistancies between reorganisation energies of equivalent dimers that should be the same.\n'
			toString += 'inconsistant_dimers_REs (equivalent_dimer: unique_dimer): \n'
			for equivalent_dimer, unique_dimer in inconsistant_dimers_REs:
				toString += str(equivalent_dimer)+': '+str(unique_dimer)+'\n'
			raise Exception(toString)

	# ==================================================================================================================================

	# Seventh, if bandgap energy data is given manually, record this and do consistancy checks. 
	if 'manual_bandgap_data' in reorganisation_and_bandgap_energy_details:

		# 7.1: Record the manually entered bandgap energy data.
		manual_bandgap_data = {}
		for key, value in reorganisation_and_bandgap_energy_details['manual_bandgap_data'].items():
			molecule_name = int(str(key).replace('S',''))
			if molecule_name in manual_bandgap_data:
				raise Exception('Error: Bandgap energy for molecule '+str(molecule_name)+' is already entered in reorganisation_and_bandgap_energy_details["manual_bandgap_data"].\nCheck this out.\nreorganisation_and_bandgap_energy_details["manual_bandgap_data"] = '+str(reorganisation_and_bandgap_energy_details["manual_bandgap_data"]))
			manual_bandgap_data[molecule_name] = value

		# 4.2: Check to make sure that a molecule entered into reorganisation_and_bandgap_energy_details["manual_bandgap_data"]
		#      actually exists in the crystal.
		inconsistant_molecules_bandgaps = []
		for (molecule1_name, molecule2_name) in manual_RE_data.keys():
			if molecule1_name not in unique_molecules+equivalent_molecules:
				inconsistant_molecules_bandgaps.append(molecule1_name)
			if molecule2_name not in unique_molecules+equivalent_molecules:
				inconsistant_molecules_bandgaps.append(molecule2_name)
		if len(inconsistant_molecules_bandgaps) > 0:
			inconsistant_molecules_bandgaps = sorted(set(inconsistant_molecules_bandgaps))
			toString  = 'Some of the molecules entered into reorganisation_and_bandgap_energy_details["manual_bandgap_data"] were not found in the crystal.'+'\n'
			toString += 'Inconsistent molecules in "manual_bandgap_data": '+str(inconsistant_molecules_bandgaps)+'\n'
			toString += 'Possible dimers between molecules '+str(sorted(unique_molecules+equivalent_molecules))+'\n'
			toString += 'Those unique molecules: '+str(unique_molecules)
			raise Exception(toString)

		# 7.3: Add manually entered bandgap energies of unique molecules to reorganisation_energy_data
		inconsistant_molecules_bandgaps = []
		for mol_name in sorted(molecule_names):

			# 7.3.1: Obtain the reorganisation energies for the unique molecules
			if (mol_name in unique_molecules):

				# 7.3.1.1: If mol_name is unique, record the reorganisation energy
				bandgap_energy_data[mol_name] = manual_bandgap_data[mol_name]

			else:

				# 7.3.1.2: If mol_name is equivalent, check that it is not inconsistent with the unique molecules. 
				unique_mol_unique = conformationally_equivalent_molecules[mol_name] if (mol_name in conformationally_equivalent_molecules.keys()) else mol_name

				# 7.3.1.3: If the user did not give a manual bandgap energy, give this for further checking reasons.
				if mol_name not in manual_bandgap_data:
					continue
					
				# 7.3.1.4: If this issue arises, then there is a missing unique molecule in reorganisation_and_bandgap_energy_details["manual_bandgap_data"].
				if unique_mol_unique not in manual_bandgap_data:
					toString  = 'Error: You have not included some unique molecules in your reorganisation_and_bandgap_energy_details["manual_bandgap_data"] dictionary.\n'
					toString += 'Check this.\n'
					toString += 'reorganisation_and_bandgap_energy_details["manual_bandgap_data"] = '+str(reorganisation_and_bandgap_energy_details["manual_bandgap_data"])+'\n'
					toString += 'Problematic molecule: '+str(unique_mol_unique)+'\n'
					toString += 'unique_molecules = '+str(unique_molecules)+'\n'
					toString += 'equivalent_molecules = '+str(equivalent_molecules)
					raise Exception(toString)

				# 7.3.1.5: Obtain the reorganisation energy for the equivalent and unique dimer
				reorganisation_energy_equivalent = manual_bandgap_data[mol_name]
				reorganisation_energy_unique     = manual_bandgap_data[unique_mol_unique]

				# 7.3.1.6: If the reorganisation energies for the equivalent and unique dimer are not the same, generate error
				if not (reorganisation_energy_equivalent == reorganisation_energy_unique):
					inconsistant_molecules_bandgaps.append((mol_name, unique_mol_unique))

		# 7.4: Indicate if there were any inconsistancies between manually entered bandgap energies for equivalent molecules.
		if len(inconsistant_molecules_bandgaps) > 0:
			toString  = 'Error: There seems to be inconsistancies between bandgap energies of equivalent molecules that should be the same.\n'
			toString += 'inconsistant_molecules_bandgaps (equivalent_molecule: unique_molecule): \n'
			for equivalent_molecule, unique_molecule in inconsistant_molecules_bandgaps:
				toString += str(equivalent_molecule)+': '+str(unique_molecule)+'\n'
			raise Exception(toString)

	# ==================================================================================================================================

	# Eighth, return dictionaries of the energy gaps for ground and excited structures for each molecule in the crystal. 
	return bandgap_energy_data, reorganisation_energy_data, conformationally_equivalent_molecules

	# ==================================================================================================================================




