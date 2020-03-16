import os
from   qet.test_systems import lah10_fm3m
from   qet.calculations import relax, phonon_grid

def run_test_lah10_fm3m():
    
    old_path = os.getcwd()
    os.system("mkdir lah10_fm3m_test")
    os.chdir("lah10_fm3m_test")

    rlx = relax(lah10_fm3m)
    rlx.run()
