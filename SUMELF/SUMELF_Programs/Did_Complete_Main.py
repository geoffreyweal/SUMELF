'''
Did_Complete_Main.py, Geoffrey Weal, 08/03/2019

This program will determine which of your dimers have been successfully calculated in Gaussian.
'''

import os, sys

from SUMELF.SUMELF_Programs.Did_Complete_Main_methods.shared_general_methods import reverse_readline

def determine_what_the_job_is(path_to_file, break_string='Input orientation:'):
    """

    """
    looking_at_ATC_output = False
    looking_at_EET_output = False
    looking_at_ICT_output = False

    with open(path_to_file,'r') as outputLOG:
        path_including_hash = False
        for line in outputLOG:
            # Get input of line from output file
            if   path_including_hash:
                if '----------------------' in line:
                    path_including_hash = False
                else:
                    current_line += line
                    continue
            elif '#' in line:
                current_line  = line
                path_including_hash = True
                continue
            else:
                current_line = line
            # Perform task with input given
            if 'density=(transition=1)' in current_line:
                looking_at_ATC_output = True
                break
            if 'eet' in current_line:
                looking_at_EET_output = True
                break

            if break_string == 'empty':
                if line.strip():
                    break
            else:
                if (break_string in current_line):
                    break

    return looking_at_ATC_output, looking_at_EET_output, looking_at_ICT_output


def did_finish_calc_on_system(path, gjffile_name, logfile_name):
    """
    This method will go through the output.log file from a electronic coupling calculation performed by Gaussian on a dimer and will determine if it completed successfully or not.
 
    Parameters
    ----------
    path : str.
        This is the path to the output.log file.

    Returns
    -------
    True if finished successfully. False if not.
    """

    path_to_gjffile   = path+'/'+gjffile_name
    path_to_outputLOG = path+'/'+logfile_name
    path_to_outputCHG = path+'/output.chg'

    # -------------------------------------------------------------------------------

    # Determine what calculations we are looking at.
    looking_at_RE_output  = False
    for file in os.listdir(path):
        if file in ['GS_GS.gjf', 'ES_ES.gjf']:
            looking_at_RE_output = True
            GS_or_ES_type = file.split('_')[0]
            break
    else:
        if os.path.exists(path_to_outputLOG):
            looking_at_ATC_output, looking_at_EET_output, looking_at_ICT_output = determine_what_the_job_is(path_to_outputLOG,break_string='Input orientation:')
        else:
            looking_at_ATC_output, looking_at_EET_output, looking_at_ICT_output = False, False, False

    no_of_true_statements = sum([looking_at_ATC_output, looking_at_EET_output, looking_at_RE_output])

    # Program likely has not begun, figure out what the job is and then return None
    if no_of_true_statements == 0:
        '''
        looking_at_ATC_output, looking_at_EET_output, looking_at_ICT_output = determine_what_the_job_is(path_to_gjffile,break_string='empty')    
        if   looking_at_ATC_output:
            file_type = "ATC"
        elif looking_at_EET_output:
            file_type = 'EET'
        elif looking_at_ICT_output:
            file_type = 'ICT'
        elif looking_at_RE_output:
            file_type = 'RE'
        else:
            file_type = None
        return file_type, None
        '''
        return None, None

    if no_of_true_statements >= 2:
        raise Exception('Huh?')

    # -------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------

    if looking_at_ATC_output:
        # Check to see if this ATC calculation has finished successfully.
        did_terminate_normally = False
        have_ATC_details = False
        counter = 0
        for line in reverse_readline(path_to_outputLOG):
            # Check if the gaussian file has terminated normally
            if 'Normal termination of Gaussian' in line:
                did_terminate_normally = True
                break
            if counter >= 20 and not did_terminate_normally:
                break

        if not did_terminate_normally:
            return 'ATC', False

        have_ATC_details = os.path.exists(path_to_outputCHG)
        return 'ATC', (did_terminate_normally and have_ATC_details)

        # ---------------------------------------------------------------------------

    elif looking_at_EET_output:
        # Check to see if this EET calculation has finished successfully.
        did_terminate_normally = False
        have_EET_details = False
        counter = 0
        for line in reverse_readline(path_to_outputLOG):
            # Check if the gaussian file has terminated normally
            if 'Normal termination of Gaussian' in line:
                did_terminate_normally = True
            # If not found the termination signal after 20 lines from the end of file, 
            # The job probably did not terminate properly
            if counter >= 20 and not did_terminate_normally:
                break
            # Indicate if EET header is found, indicate true
            if 'Electronic Coupling for Excitation Energy Tranfer' in line:
                have_EET_details = True
            # If this is found, prbably should have found the EET by then, so dont bother anymore
            if 'initial guesses have been made' in line:
                break
            # If these tags are true, finish
            if did_terminate_normally and have_EET_details:
                break
            counter += 1

        return 'EET', (did_terminate_normally and have_EET_details)

        # ---------------------------------------------------------------------------

    elif looking_at_RE_output:
        # Check to see if this RE calculation has finished successfully.
        optimisation_file = GS_or_ES_type+'_'+GS_or_ES_type+'.log'
        if GS_or_ES_type == 'ES':
            freq_file = GS_or_ES_type+'_'+GS_or_ES_type+'_freq.log'
        single_point_file = GS_or_ES_type+'_'+convert_GS_and_ES(GS_or_ES_type)+'.log'
        # --------------
        # First, check if the optimisation finished
        did_RE_opt_terminate_normally = False
        have_RE_opt_details = False
        counter = 0
        if os.path.exists(path+'/'+optimisation_file):
            for line in reverse_readline(path+'/'+optimisation_file):
                # Check if the gaussian file has terminated normally
                if 'Normal termination of Gaussian' in line:
                    did_RE_opt_terminate_normally = True
                # If not found the termination signal after 20 lines from the end of file, 
                # The job probably did not terminate properly
                if counter >= 20 and not did_RE_opt_terminate_normally:
                    break
                # Indicate if EET header is found, indicate true
                if 'Stationary point found' in line:
                    have_RE_opt_details = True
                # If this is found, prbably should have found the EET by then, so dont bother anymore
                if 'Predicted change in Energy' in line:
                    break
                # If these tags are true, finish
                if did_RE_opt_terminate_normally and have_RE_opt_details:
                    break
                counter += 1
        # --------------
        # Second, check if the freq calc finished if you performed opt on a ES_ES.gjf file
        if GS_or_ES_type == 'ES':
            did_RE_freq_terminate_normally = False
            counter = 0
            if os.path.exists(path+'/'+freq_file):
                for line in reverse_readline(path+'/'+freq_file):
                    # Check if the gaussian file has terminated normally
                    if 'Normal termination of Gaussian' in line:
                        did_RE_freq_terminate_normally = True
                    # If not found the termination signal after 20 lines from the end of file, 
                    # The job probably did not terminate properly
                    if counter >= 20 and not did_RE_freq_terminate_normally:
                        break
                    # If these tags are true, finish
                    if did_RE_freq_terminate_normally:
                        break
                    counter += 1
        # --------------
        # Third, check if the single point calc finished if you performed opt on a ES_ES.gjf file
        did_RE_single_point_terminate_normally = False
        counter = 0
        if os.path.exists(path+'/'+single_point_file):
            for line in reverse_readline(path+'/'+single_point_file):
                # Check if the gaussian file has terminated normally
                if 'Normal termination of Gaussian' in line:
                    did_RE_single_point_terminate_normally = True
                # If not found the termination signal after 20 lines from the end of file, 
                # The job probably did not terminate properly
                if counter >= 20 and not did_RE_single_point_terminate_normally:
                    break
                # If these tags are true, finish
                if did_RE_single_point_terminate_normally:
                    break
                counter += 1
        # -------------- 
        # Return results of calculations
        if GS_or_ES_type == 'ES':
            return 'RE', ((did_RE_opt_terminate_normally and have_RE_opt_details), did_RE_freq_terminate_normally, did_RE_single_point_terminate_normally)
        else:
            return 'RE', ((did_RE_opt_terminate_normally and have_RE_opt_details), did_RE_single_point_terminate_normally)
        # -------------- 

    # -------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------

def convert_GS_and_ES(input_type):
    if   input_type == 'GS':
        return 'ES'
    elif input_type == 'ES':
        return 'GS'
    else:
        raise Exception('Error in def convert_GS_and_ES, in Did_Complete_Main.py. Can only be ES or GS input_type.')

def determine_if_all_gaussian_job_completed(general_path, print_progress=False):
    """
    This method will go through folders in search of output.log files, and will determine from those output.log files if they had finished successfully or not.

    Parameters
    ----------
    general_path : str.
        This is the overall directory to search through for output.log files.

    Returns
    -------
    
    """

    atc_jobs_finished_successfully = []
    atc_jobs_finished_unsuccessfully = []

    eet_jobs_finished_successfully = []
    eet_jobs_finished_unsuccessfully = []

    re_jobs_finished_successfully = []
    re_jobs_finished_unsuccessfully = []

    unalligned_jobs = []

    for root, dirs, files in os.walk(general_path):

        # If output.log not found, move on to the next entry in the for loop
        if   'GS_GS.gjf' in files:
            logfile_name = 'GS_GS.log'
        elif 'ES_ES.gjf' in files:
            logfile_name = 'ES_ES.log'
        else:
            for file in files:
                if '.gjf' in file:
                    gjffile_name = file
                    logfile_name = 'output.log' #file.replace('.log', '.gjf')
                    break
            else:
                continue

        # Print details of where the program is up to:
        if print_progress:
            sys.stdout.write("\r                                                                  ")
            sys.stdout.flush()
            sys.stdout.write("\r"+str(root))
            sys.stdout.flush()

        # Go through the output.log file to see if the job finished successfully or not.
        atc_eet_or_re, did_finish = did_finish_calc_on_system(root, gjffile_name, logfile_name)

        # Add job path to appropriate list
        if atc_eet_or_re == 'ATC':
            atc_jobs_finished_successfully.append(root) if did_finish else atc_jobs_finished_unsuccessfully.append(root)
        elif atc_eet_or_re == 'EET':
            eet_jobs_finished_successfully.append(root) if did_finish else eet_jobs_finished_unsuccessfully.append(root)
        elif atc_eet_or_re == 'RE':
            re_jobs_finished_successfully.append(root) if all(did_finish) else re_jobs_finished_unsuccessfully.append((root, did_finish))
        else:
            unalligned_jobs.append(root)
        # This will prevent the program looking further down the directories.
        dirs[:] = []
        files[:] = []

    atc_jobs_finished_successfully.sort()
    atc_jobs_finished_unsuccessfully.sort()

    eet_jobs_finished_successfully.sort()
    eet_jobs_finished_unsuccessfully.sort()

    re_jobs_finished_successfully.sort()
    re_jobs_finished_unsuccessfully.sort()

    unalligned_jobs.sort()

    return atc_jobs_finished_successfully, atc_jobs_finished_unsuccessfully, eet_jobs_finished_successfully, eet_jobs_finished_unsuccessfully, re_jobs_finished_successfully, re_jobs_finished_unsuccessfully, unalligned_jobs




