import multiprocessing as mp
from qet.plots import plot_tc_vs_smearing
import os
import time

dir_data = []
for f in os.listdir("."):
    if not os.path.isdir(f): continue
    dir_data.append(f)
dir_data.sort()

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
