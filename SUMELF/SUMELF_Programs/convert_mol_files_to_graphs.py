'''
convert_mol_files_to_graphs.py, Geoffrey Weal, 15/2/24

This script is designed to convert mol files to graphs.
'''
from SUMELF.SUMELF.general_methods.is_solvent_methods.solvent_graph_methods import create_graph_from_mol_files

class CLICommand:
    """Will mconvert mol files to graphs.
    """

    @staticmethod
    def add_arguments(parser):
        pass

    @staticmethod
    def run(args):
        Run_method()

def Run_method():
    """
    This method is designed to convert mol files to graphs.
    """
    create_graph_from_mol_files()
