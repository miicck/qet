import os
from   qet.test.test_systems import lah10_fm3m
from   qet.calculations      import relax, phonon_grid, dos, proj_dos, calculate_tc
from   qet.logs              import log

def test_tc_calculator():
    
    # Create and move to a directory to
    # run the test from
    old_path = os.getcwd()
    os.system("mkdir lah10_fm3m_tc_test")
    os.chdir("lah10_fm3m_tc_test")

    # Use cheapish parameters
    lah10_fm3m["kpts_per_qpt"] = 3
    lah10_fm3m["ecutwfc"]      = 40

    try: calculate_tc(lah10_fm3m)
    except: pass

    # Restore old directory
    os.chdir(old_path)

def run_test_lah10_fm3m():
    
    # Create and move to a directory to
    # run the test from
    old_path = os.getcwd()
    os.system("mkdir lah10_fm3m_test")
    os.chdir("lah10_fm3m_test")
    try:

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

    except: pass

    # Restore old directory
    os.chdir(old_path)

def print_test_lah10_fm3m():

    # Print the parameters
    print(lah10_fm3m)

    # Print relax.in file
    rlx = relax(lah10_fm3m)
    print("\n"+str(rlx))

    # Print dos.in file
    ds  = dos(lah10_fm3m)
    print("\n"+str(ds))

    # Print proj_dos.in file
    pdos = proj_dos(lah10_fm3m)
    print("\n"+str(pdos))

    # Print ph.in file
    ph  = phonon_grid(lah10_fm3m)
    print("\n"+str(ph))

def test_alchemy():
    
    from qet.alchemy           import mutations
    from qet.alchemy.optimize  import optimize
    from qet.test.test_systems import lah10_fm3m

    # Start with Fm3m phase of LaH10
    start = lah10_fm3m

    # Run optimizer
    path = optimize(
        start,                                  # Start structure
        lambda s : len(s["atoms"]),             # Minimize number of atoms
        [
            mutations.remove_random_atom,       # Allowed mutations
            mutations.duplicate_random_atom,
        ],
        lambda s : s["atom_counts"]["H"] > 0,   # Must have at least 1 H

        # We consider all structures with the same atoms to be the same
        structure_compare = lambda l1, a1, l2, a2 : True

        )
