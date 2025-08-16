'''
SUMELF_Did_Complete.py, Geoffrey Weal, 11/04/2022

This program will determine which of your dimers have been successfully calculated in Gaussian.
'''
import os

from ECCP.ECCP_Programs.Did_Complete_Main import determine_if_all_gaussian_job_completed

class CLICommand:
    """Will determine which ATC and EET jobs have completed and which ones have not.
    """

    @staticmethod
    def add_arguments(parser):
        pass

    @staticmethod
    def run(args):
        Run_method()

def Run_method():
    """
    This method will determine which of your dimers have been successfully calculated in Gaussian.
    """
    general_path = os.getcwd()

    print('########################################################################')
    print('########################################################################')
    print('Checking if OPV_Dimer_Pairer job have completed in Gaussian')
    print('-----------------------------------------------------------')

    atc_jobs_finished_successfully, atc_jobs_finished_unsuccessfully, eet_jobs_finished_successfully, eet_jobs_finished_unsuccessfully, re_jobs_finished_successfully, re_jobs_finished_unsuccessfully, unalligned_jobs = determine_if_all_gaussian_job_completed(general_path, print_progress=True)

    print('########################################################################')
    print('########################################################################')
    print('-------------------')
    print('Jobs that completed'.upper())
    print('-------------------')
    print('------------------')
    print('     ATC Jobs     ')
    print('------------------')
    for dirpath in atc_jobs_finished_successfully:
        print(dirpath)
    print('------------------')
    print('     EET Jobs     ')
    print('------------------')
    for dirpath in eet_jobs_finished_successfully:
        print(dirpath)
    print('------------------')
    print('      RE Jobs     ')
    print('------------------')
    for dirpath in re_jobs_finished_successfully:
        print(dirpath)
    if (not len(atc_jobs_finished_unsuccessfully+re_jobs_finished_unsuccessfully+eet_jobs_finished_unsuccessfully) == 0):
        print('########################################################################')
        print('########################################################################')
        print('--------------------------')
        print('Jobs that did not complete'.upper())
        print('--------------------------')
        if (not len(atc_jobs_finished_unsuccessfully) == 0):
            print('------------------')
            print('     ATC Jobs     ')
            print('------------------')
            for dirpath in atc_jobs_finished_unsuccessfully:
                print(dirpath)
        if (not len(eet_jobs_finished_unsuccessfully) == 0):
            print('------------------')
            print('     EET Jobs     ')
            print('------------------')
            for dirpath in eet_jobs_finished_unsuccessfully:
                print(dirpath)
        if (not len(re_jobs_finished_unsuccessfully) == 0):
            print('------------------')
            print('      RE Jobs     ')
            print('------------------')
            for dirpath, results in re_jobs_finished_unsuccessfully:
                if len(results) == 2:
                    did_opt_complete, did_sp_calc_complete = results
                    did_opt_complete = 'YES' if did_opt_complete else 'NO'
                    did_sp_calc_complete = 'YES' if did_sp_calc_complete else 'NO'
                    print(dirpath+' (OPT: '+str(did_opt_complete)+', SP: '+str(did_sp_calc_complete)+', FREQ: '+str(did_sp_calc_complete)+')')
                else:
                    did_opt_complete, did_freq_complete, did_sp_calc_complete = results
                    did_opt_complete = 'YES' if did_opt_complete else 'NO'
                    did_freq_complete = 'YES' if did_freq_complete else 'NO'
                    did_sp_calc_complete = 'YES' if did_sp_calc_complete else 'NO'
                    print(dirpath+' (OPT: '+str(did_opt_complete)+', FREQ: '+str(did_freq_complete)+', SP: '+str(did_sp_calc_complete)+')')
    if (not len(unalligned_jobs) == 0):
        print('########################################################################')
        print('########################################################################')
        print('--------------------------------------')
        print('Jobs that were neither ATC, RE, or EET jobs, or have not begun'.upper())
        print('--------------------------------------')
        for dirpath in unalligned_jobs:
            print(dirpath)
    print('########################################################################')
    print('########################################################################')





