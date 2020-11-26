import multiprocessing as mp
from qet.plots import plot_tc_vs_smearing
from qet.plots import plot_phonon_mode_atoms
import os
import time
import sys

dir_data = []
for f in os.listdir("."):
    if not os.path.isdir(f): continue
    dir_data.append(f)

starting_from = None
for a in sys.argv[1:]:
    if a.startswith("starting_from="):
        starting_from = a.split("=")[-1]
        sys.argv.remove(a)

selected_dirs = []
for f in sys.argv[1:]:
    found = None
    for d in dir_data:
        if d.startswith(f):
            found = d
            break
    if found is None:
        print(f+" not found!")
    else:
        print(f+" -> "+found)
        selected_dirs.append(found)

if len(selected_dirs) > 0: dir_data = selected_dirs

dir_data.sort()
if not starting_from is None:
    dir_data = dir_data[dir_data.index(starting_from):]

cpus    = 2
running = []
for f in dir_data:

    while len(running) >= cpus:

        for p in running:
            if not p.is_alive():
                running.remove(p)

        time.sleep(0.2)

    p = mp.Process(target=plot_tc_vs_smearing, args=([f],))
    p.start()
    running.append(p)

    for k in os.listdir(f):
        if not k.startswith("kpq"): continue
        m = f+"/"+k+"/ph_interp.modes"
        if not os.path.isfile(m): continue
        p2 = mp.Process(target=plot_phonon_mode_atoms, args=(m,))
        p2.start()
        running.append(p2)
        break
