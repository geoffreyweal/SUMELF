
def get_reorganisation_energy(donor_eGS_gGS_energy, donor_eGS_gES_energy, acceptor_eES_gGS_energy, acceptor_eES_gES_energy):
    '''
    This algorithm is designed to return the reorganisation for an exciton to jump from the donor molecule to the acceptor molecule.
    
    See https://pubs.rsc.org/en/content/articlepdf/2021/tc/d0tc05697a for more details on getting the reorganisation energy.
    
    Parameters
    ----------
    donor_eGS_gGS_energy : float
        This is the ground state energy of the donor molecule in the ground state geometry.
    donor_eGS_gES_energy : float
        This is the ground state energy of the donor molecule in the excited state geometry.
    acceptor_eES_gGS_energy : float
        This is the excited state energy of the acceptor molecule in the ground state geometry.
    acceptor_eES_gES_energy : float
        This is the excited state energy of the acceptor molecule in the excited state geometry.

    Returns
    -------
    The reorganisation energy for the exciton to jump from the donor molecule to the acceptor molecule (in eV). 
    '''
    return (donor_eGS_gES_energy - donor_eGS_gGS_energy) + (acceptor_eES_gGS_energy - acceptor_eES_gES_energy)

hartree_to_eV = 27.211386245988
def convert_hartree_to_eV(hartree_energy):
    return hartree_energy * hartree_to_eV

def get_hartree_to_eV_conversion_value():
    return hartree_to_eV