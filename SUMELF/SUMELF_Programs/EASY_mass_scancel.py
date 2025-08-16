'''
SUMELF_mass_scancel.py, Geoffrey Weal, 14/8/22

This program is designed to allow the user to cancel a number of jobs between two slurm ids.
'''
from subprocess import Popen #call, 

class CLICommand:
    """Will determine which ATC and EET jobs have completed and which ones have not.
    """

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('lower_id',  nargs='*', help='This is the Slurm ID to begin from.')
        parser.add_argument('higher_id', nargs='*', help='This is the Slurm ID to end on.')
        parser.add_argument('time_between_250_cancels', nargs='*', help='This is the wait time you want to give per 250 cancellations.')

    @staticmethod
    def run(args):
        id_low  = args_submit.lower_id
        id_high = args_submit.higher_id

        time_between_250_cancels = args_submit.time_between_250_cancels
        if len(time_between_250_cancels) > 0:
            time_between_250_cancels = int(time_between_250_cancels[0])
        else:
            time_between_250_cancels = 60

        Run_method(id_low, id_high, time_between_250_cancels)

def Run_method(id_low, id_high, time_between_250_cancels):
    """
    This program is designed to allow the user to cancel a number of jobs in slurm between two slurm ids.
    """
    for scancel_id in range(id_low,id_high+1):
        args = ["scancel",str(scancel_id)]
        print('Cancelling Job id: '+str(scancel_id))
        Popen(args)
        process.wait()
        if counter >= 250:
            pbar = trange(time_between_250_cancels)
            pbar.set_description('Have cancelled 250 jobs consecutively. ')
            for _ in pbar:
                time.sleep(1)
            counter = 0
        else:
            counter += 1