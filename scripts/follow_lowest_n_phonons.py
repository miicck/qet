from qet.params import parameters
import sys
import os

try:
    param_file = sys.argv[1]
    modes_file = sys.argv[2]
    n = int(sys.argv[3])
except:
    print("Usage: python ~/../follow_phonon.py param.in matdyn.modes n")
    quit()

p = parameters(param_file)

hsp = {}
hsp_coords = p["high_symmetry_bz_points"]
for i in hsp_coords:
    name = hsp_coords[i]
    if name in hsp: continue
    hsp[name] = p["bz_path"][i]

for mev in [1.0, 10.0, 100.0]:
    for name in hsp:
        for i in range(n):
            q, f, vec, disp = p.apply_phonon(hsp[name], i, modes_file=modes_file, 
                use_closest=True, return_freq=True, return_evec=True, return_disp=True, target_energy=mev)

            max_d = 0.0
            for j in range(0, len(q["atoms"])):
                if q["atoms"][j][0] == "La":
                    d = sum(x*x for x in disp[j])
                    d = d**0.5
                    if d > max_d:
                        max_d = d

            d = "{0}meV".format(int(mev))
            os.system("mkdir "+d+" 2>/dev/null")
            fname = "{0}/{1}_phonon_mode_{2}.cell".format(d, name, i)
            q.save_to_cell(fname)

            with open(fname.replace(".cell", ".info"), "w") as inf:
                inf.write("Q-point {0} {1} mode {2}\n".format(name, hsp[name], i))
                inf.write("frequency = {0} cm^-1\n".format(f))
                inf.write("Maxmum La displacement = {0} angstrom\n".format(max_d))

            os.system("c2x --cif "+fname+" "+fname.replace(".cell", ".cif"))
