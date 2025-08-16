'''
get_no_of_electrons.py, Geoffrey Weal, 20/6/2024

This program will count the total number of electrons in your system
'''
import os, csv
from tqdm import tqdm
from ase.io import read
from collections import Counter

class CLICommand:
    """Will determine which geometric optimisation jobs have completed and which ones have not.
    """

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('--print_progress', nargs=1, help='Determines if you want to print the process.', default=['True'])

    @staticmethod
    def run(args):

        # First, determine if the user want to print the information about the process to the screen.
        print_progress = args.print_progress
        if len(print_progress) == 1:
            print_progress = print_progress[0]
        else:
            raise Exception('Error: Only one input can be given for print_progress')
        if print_progress.lower() in ['true', 't']:
            print_progress = True
        elif print_progress.lower() in ['false', 'f']:
            print_progress = False
        else:
            raise Exception('Error: print_progress must be either True or False')

        # Second, run the program
        Run_method(print_progress=print_progress)

def Run_method(print_progress=True):
    """
    This method will count the number of electrons in the system. 

    Parameters
    ----------
    print_progress : bool.
        This indicates if you want to print the information about the progress to the terminal.
    """

    # First, initialise lists to record results to.
    electron_info_details = open('electron_information.csv', 'w')

    # Second, set up the file to be used as a csv file
    fieldnames = ['filepath', 'No_of_electrons', 'No_of_base_electrons', 'No_of_excess_electrons (from charge)', 'Lowest_Multiplicity']
    csv_writer = csv.writer(electron_info_details)
    csv_writer.writerow(fieldnames)

    # Third, obtain the current working directory. 
    original_path = os.getcwd()

    # Fourth, obtain the os.walk object.
    pbar = os.walk(original_path)

    # Fifth, if you want to print the progress of this script to a progress bar, do it here.
    if print_progress:
        pbar = tqdm(pbar, unit=' paths examined')

    # Sixth, determine all the Gaussian jobs to check. 
    for root, dirs, files in pbar:

        # 6.1: Sort dirs and files so that things come out in alphabetical order.
        dirs.sort()
        files.sort()

        # 6.2: For each file in the current directory.
        for file in files:

            # 6.2.1: Check if this is file is appropriate to process.
            is_xyz_file          = file.endswith('.xyz')
            #is_gaussian_gjf_file = file.endswith('.gjf')
            #is_orca_inp_file     = file.endswith('.inp')

            # 6.2.2: Check this file is not appropriate, move on to the next file. 
            if not (is_xyz_file): # or is_gaussian_gjf_file or is_orca_inp_file):
                continue

            # 6.2.3: Open the file in ASE
            system = read(root+'/'+file)

            # 6.2.4: Get the elements in this molecule.
            element_information = Counter(system.get_chemical_symbols())

            # 6.2.5: This is the atomic numbers of the atoms in the system. This is the same as the number of electrons the atoms have, 
            atom_number_information = system.get_atomic_numbers()

            # 6.2.6: Obtain any excess charges on the atoms in the system
            get_atom_charges = system.get_initial_charges()

            # 6.2.7: Check that all atoms have an interger charge.
            if any([not (charge % 1.0 == 0.0) for charge in get_atom_charges]):
                raise Exception('Error: Some of the atoms in molecule '+str(file)+' do not have an integer charge. This is a problem. See: '+root+'/'+file)

            # 6.2.8: Get the number of electrons in the system. The sum(get_atom_charges) is negative because electrons has a negative charge, so 
            #        * negative charge means extra electrons.
            #        * positive charge means missing electrons.

            number_of_base_electrons   = sum(atom_number_information) 

            number_of_excess_electrons = -int(sum(get_atom_charges))


            total_number_of_electrons  = number_of_base_electrons + number_of_excess_electrons

            #import pdb; pdb.set_trace()

            # 6.2.9: Double check again that there is an integer number of electrons in the system.
            if not (total_number_of_electrons % 1.0 == 0.0):
                raise Exception('Error: Some of the atoms in molecule '+str(file)+' do not have an integer charge. This is a problem. See: '+root+'/'+file)

            # 6.2.10: Determine what the predicted lowest multiplicity is for the molecule.
            if total_number_of_electrons % 2 == 0.0:
                lowest_multiplicity = 1
            else:
                lowest_multiplicity = 2

            # 6.2.11: Save the electron information to file.
            csv_writer.writerow([root+'/'+file, total_number_of_electrons, number_of_base_electrons, number_of_excess_electrons, lowest_multiplicity])

    # Seventh: Close the electron information text file. 
    electron_info_details.close()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 



