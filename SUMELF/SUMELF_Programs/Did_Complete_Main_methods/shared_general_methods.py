'''
Geoffrey Weal, processing_OPV_Dimer_data.py, 9/3/22

This script contains methods for processing_OPV_Dimer_data.py

'''
import os

# -----------------------------------------------------------------

def found_a_gaussian_job_that_has_run(files, root):
    '''
    This method is designed to determine if a Gaussian job is running or is currently running.

    Parameters
    ----------
    files : list
        This is a list of all the filenames to look through to see it contains Gaussian input and output files.
    root : str.
        This is the path to where the calculation files are located

    Returns
    -------
    True if both the input .gjf file and the output .log files are found. 
    '''
    does_gjf_exist = False
    does_log_exist = False

    for file in files:
        if file.endswith('.gjf'):
            does_gjf_exist = True
            continue
        elif file.endswith('.log'):
            does_log_exist = True
            continue

    return (does_gjf_exist and does_log_exist)

# -----------------------------------------------------------------

def did_gaussian_job_complete(log_filepath):
    '''
    This method is designed to quickly determine if a Gaussian job has completed

    Parameters
    ----------
    log_filepath : str
        This is the path to the log file.

    Returns
    -------
    True if the log file indicates the Gaussian log file indicates that their has been normal termination, otherwise return False
    '''
    last_line = get_lastline(log_filepath)
    return ('Normal termination of Gaussian' in last_line)

def get_lastline(log_filepath):
    '''
    This method is designed to obtain the last line in a file efficiently. 

    Parameters
    ----------
    log_filepath : str
        This is the path to the log file.

    Returns
    -------
    The last line in a file.
    '''
    with open(log_filepath, 'rb') as f:
        try:  # catch OSError in case of a one line file 
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        last_line = f.readline().decode()
    return last_line

# -----------------------------------------------------------------

gaussian_temp_types = ['.chk','.d2e','.int','.rwf','.skr']
other_files_to_remove = ['core.']
def gaussian_temp_files_to_remove(root, files, remove_fort7_file=False, print_to_display=True):
    '''
    This method is designed to remove temporary gaussian files. 

    These are files that will not need to be used again, and many are very large so removing is very advantagous.

    Parameters
    ----------
    root : str
        This is the path to the gaussian temp files
    files : list
        This is a list of all the filenames to look through to see it contains Gaussian temp files.
    remove_fort7_file : bool.
        If True, remove the fort.7 file if found. If False, do not remove this file.
    print_to_display : bool.
        Print what is being remove to screen. Default: True
    '''

    # For determining if to remove wfn file as multiwfn file has already been processed for ATC calculation
    wfn_file = None
    have_chg_file = False

    # Get files to remove
    temp_files_to_remove = []
    for file in files:
        for gaussian_temp_type in gaussian_temp_types:
            if file.endswith(gaussian_temp_type):
                temp_files_to_remove.append(file)
        for other_file_to_remove in other_files_to_remove:
            if other_file_to_remove in file:
                temp_files_to_remove.append(file)
        if file.endswith('.wfn'): 
            wfn_file = file
        if file.endswith('.chg'): 
            have_chg_file = True
        if (file == 'fort.7') and remove_fort7_file:
            temp_files_to_remove.append(file)

    # Determine if to remove wfn file
    if have_chg_file and (wfn_file is not None):
        temp_files_to_remove.append(wfn_file)

    # Remove files if they exist
    if len(temp_files_to_remove) > 0:
        if print_to_display:
            print('Tidying files in: '+str(root))
            print('Removing Gaussian temp files and other unnecessary files: '+str(temp_files_to_remove))
        for temp_file_to_remove in temp_files_to_remove:
            path_to_temp_file_to_remove = root+'/'+temp_file_to_remove
            if os.path.exists(path_to_temp_file_to_remove):
                os.remove(path_to_temp_file_to_remove)
    else:
        if print_to_display:
            print('No files need to be removed to tidy: '+str(root))

# -----------------------------------------------------------------

def reverse_readline(filename, buf_size=8192):
    """
    A generator that returns the lines of a file in reverse order

    Parameters
    ----------
    filename : str
        This is the path to the file you want to read.
    buf_size : int
        This is the buffer size to read in.

    Returns
    -------
    Returns each line in the file in reverse order.
    """
    
    with open(filename) as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first 
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment

# -----------------------------------------------------------------




