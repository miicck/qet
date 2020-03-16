import os
from   qet.test_systems import lah10_fm3m
from   qet.calculations import relax, phonon_grid

def run_test_lah10_fm3m():
    
    # Create and move to a directory to
    # run the test from
    old_path = os.getcwd()
    os.system("mkdir lah10_fm3m_test")
    os.chdir("lah10_fm3m_test")

    # Relax the structure
    rlx = relax(lah10_fm3m)
    rlx.run()

    # Calculate phonons on a grid
    ph  = phonon_grid(lah10_fm3m)
    ph.run()

    # Restore old directory
    os.chdir(old_path)
