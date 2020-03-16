import os
from   qet.test_systems import lah10_fm3m
from   qet.calculations import relax, phonon_grid, dos, proj_dos
from   qet.logs         import log

def run_test_lah10_fm3m():
    
    # Create and move to a directory to
    # run the test from
    old_path = os.getcwd()
    os.system("mkdir lah10_fm3m_test")
    os.chdir("lah10_fm3m_test")

    # Relax the structure
    rlx = relax(lah10_fm3m)
    log(rlx.run())

    # Calculate the density of states
    ds  = dos(lah10_fm3m)
    log(ds.run())

    # Calculate the projected density of states
    pdos = proj_dos(lah10_fm3m)
    log(pdos.run())

    # Calculate phonons on a grid
    ph  = phonon_grid(lah10_fm3m)
    log(ph.run())

    # Restore old directory
    os.chdir(old_path)
