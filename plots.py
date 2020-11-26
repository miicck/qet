from qet.params import parameters as params
import qet.parser    as parser
import qet.constants as constants
import os
import math

# Get a series of plottable colors
# indexed by n, cycling
def color_cycle(n):
    arr = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 1.0, 0.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 0.5, 0.0],
        [1.0, 0.0, 0.5],
        [0.5, 1.0, 0.0],
        [0.0, 1.0, 0.5],
        [0.5, 0.0, 1.0],
        [0.0, 0.5, 1.0]
    ]

    return arr[n % len(arr)]

# Plot the density of states calculated
# with output file given by filename
def plot_dos(filename="./dos.out"):
    import matplotlib.pyplot as plt
    
    # Parse the DOS
    out = parser.dos_out(filename)

    plt.plot(out["DOS energies"], out["DOS (energy)"], color="blue", label="DOS ($E$)")
    plt.axvline(out["fermi energy"], color="red", label="Fermi energy")
    plt.axhline(out["DOS (E_F)"], color="green", label="DOS ($E_f$)")
    plt.legend()
    plt.show()

# Plot the phonon density of states in 
# the given .dos file
def plot_pdos(filename="./ph_interp.dos"):
    import matplotlib.pyplot as plt
    
    out = parser.phonon_interp_dos_out(filename)
    plt.plot(out["frequencies"], out["dos"])
    plt.xlabel("Frequency (cm$^{-1}$)")
    plt.ylabel("Phonon DOS")
    plt.show()

# Plot the projected density of states calculated
# with output file given by filename
def plot_proj_dos(filename="./proj_dos.out"):
    import matplotlib.pyplot as plt
    
    # Parse the projected density of states
    out = parser.proj_dos_out(filename)

    dos_plot     = plt.subplot(311)
    df_plot      = plt.subplot(312)
    tot_dos_plot = plt.subplot(313)
    
    labels = []
    dfs    = []

    e_fermi   = out["fermi energy"]
    total_dos = None

    # Loop over atoms/wavefunctions
    for atom_num in sorted(out["PDOS energies"]):
        for wfn_num in sorted(out["PDOS energies"][atom_num]):

            # Get energies, projected dos, atom name 
            # and wavefunction name
            es = out["PDOS energies"][atom_num][wfn_num]
            ps = out["PDOS (energy)"][atom_num][wfn_num]
            an = out["PDOS atom names"][atom_num][wfn_num]
            wf = out["PDOS wfc names"][atom_num][wfn_num]
            df = out["PDOS (fermi energy)"][atom_num][wfn_num]

            # Plot in eV relative to fermi level
            es = [(x - e_fermi)*constants.RY_TO_EV for x in es]

            if total_dos is None: total_dos = ps
            else: total_dos = [t+p for t,p in zip(total_dos, ps)]

            # Plot this projected DOS
            label = "Atom {0} ({1}) wfn {2} ({3})"
            label = label.format(atom_num, an, wfn_num, wf)
            dos_plot.plot(es,ps,label=label)

            labels.append(label)
            dfs.append(df)

    dos_plot.axvline(0.0, linestyle=":", label="Fermi energy")
    dos_plot.set_xlabel("Energy (eV)")
    dos_plot.set_ylabel("PDOS")
    dos_plot.set_title("PDOS")
    dos_plot.legend(ncol=int(len(labels)**0.5))

    tot_dos_plot.plot(es, total_dos)
    tot_dos_plot.set_xlabel("Energy (eV)")
    tot_dos_plot.set_ylabel("Total DOS")
    tot_dos_plot.axvline(0.0, linestyle=":", label="Fermi energy")
    tot_dos_plot.legend()

    # Plot the densities of states at the fermi level
    x = range(0, len(labels))
    df_plot.bar(x, dfs)
    df_plot.set_xticks(x)
    df_plot.set_xticklabels([])
    df_plot.set_ylabel("PDOS")
    df_plot.set_title("Density of states $(E_{fermi})$")
    
    for i in range(0, len(labels)):
        df_plot.text(x[i], 0, "  " + labels[i], rotation=90)

    # Show the plot
    plt.show()

def plot_h_derived_dos(filename="./proj_dos.out", plot=True):
    if plot:
        import matplotlib.pyplot as plt
    
    # Parse the projected density of states
    out = parser.proj_dos_out(filename)

    # Get the fermi energy
    e_fermi = out["fermi energy"]

    h_dos = None
    non_h_dos = None
    
    # Loop over atoms/wavefunctions
    for atom_num in sorted(out["PDOS energies"]):
        for wfn_num in sorted(out["PDOS energies"][atom_num]):

            # Get energies, projected dos, atom name 
            # and wavefunction name
            es = out["PDOS energies"][atom_num][wfn_num]
            ps = out["PDOS (energy)"][atom_num][wfn_num]
            an = out["PDOS atom names"][atom_num][wfn_num]
            wf = out["PDOS wfc names"][atom_num][wfn_num]
            df = out["PDOS (fermi energy)"][atom_num][wfn_num]

            # Plot in eV relative to fermi level
            es = [(x - e_fermi)*constants.RY_TO_EV for x in es]

            # Accumulate hydrogen, or non-hydrogen dos
            if an.lower() == "h":
                if h_dos is None: h_dos = ps
                else: h_dos = [h + p for h,p in zip(h_dos, ps)]
            else:
                if non_h_dos is None: non_h_dos = ps
                else: non_h_dos = [n + p for n,p in zip(non_h_dos, ps)]

    if (h_dos is None):
        print(filename+" has no hydrogen!")
        return [0, False]

    if (non_h_dos is None):
        print(filename+" is purely hydrogen!")
        return [0, False]

    # Work out the dos ratio
    ratio = [0 if abs(h) < 10e-10 else h/(h+n) for h,n in zip(h_dos, non_h_dos)]

    # Work out the index of the fermi energy
    i_fermi = 0
    min_abs = float("inf") 
    for i in range(0, len(es)):
        if abs(es[i]) < min_abs:
            min_abs = abs(es[i])
            i_fermi = i

    val_band_max = i_fermi
    con_band_min = i_fermi

    dos_thresh = 10e-3
    if h_dos[i_fermi] + non_h_dos[i_fermi] < dos_thresh:

        # This is an insulator, identify the band edges
        while con_band_min < len(es):
            if h_dos[con_band_min] + non_h_dos[con_band_min] > dos_thresh:
                break
            con_band_min = con_band_min + 1

        while val_band_max >= 0:
            if h_dos[val_band_max] + non_h_dos[val_band_max] > dos_thresh:
                break
            val_band_max = val_band_max - 1

    # Will contain doping values and coresponding hydrogen-dos ratio changes
    doping_res = []

    # Work out the result of positive doping
    for i in range(con_band_min, len(es)):
        doping = es[i] - es[con_band_min]
        if abs(doping) > 0.2: break # Too much doping
        doping_res.append([doping, ratio[i] - ratio[con_band_min], h_dos[i]])
    
    # Work out the result of negative doping
    for i in range(val_band_max, -1, -1):
        doping = es[i] - es[val_band_max]
        if abs(doping) > 0.2: break # Too much doping
        doping_res.append([doping, ratio[i] - ratio[val_band_max], h_dos[i]])

    # Work out if this is a metal, or an insulator
    is_metal = val_band_max == con_band_min

    doping_res.sort()
    doping, dosinc, doping_hdos = zip(*doping_res)

    # Work out a score based on the maximum increase in the DOS
    i_max = dosinc.index(max(dosinc))
    score = dosinc[i_max]/abs(doping[i_max])

    if plot:

        # Plot stuff
        plt.subplot(311)
        plt.xlabel("Energy (eV)")
        plt.ylabel("DOS")
        plt.plot(es, h_dos, label="Hydrogen derived DOS")
        plt.plot(es, non_h_dos, label="Non-hydrogen DOS")
        plt.axvline(0, color="black", linestyle=":", label="Fermi energy")
        plt.axvline(es[val_band_max], color="blue")
        plt.axvline(es[con_band_min], color="blue")
        plt.legend()

        plt.subplot(312)
        plt.xlabel("Energy (eV)")
        plt.ylabel("Hydrogen-derived DOS ratio")
        plt.plot(es, ratio)
        plt.axvline(0, color="black", linestyle=":", label="Fermi energy")
        plt.axvline(es[val_band_max], color="blue")
        plt.axvline(es[con_band_min], color="blue")
        plt.legend()

        plt.subplot(313)
        plt.xlabel("Doping (eV)")
        plt.ylabel("Change in Hydrogen DOS ratio")
        plt.plot(doping, dosinc)
        plt.axhline(0, color="black", label="No DOS ratio change")
        plt.axvline(doping[i_max], color="green", label="Max DOS increase, score = {0}".format(score))
        plt.legend()

        plt.show()

    return [score, is_metal]

def rank_doping(dirs):

    score_name = []
    for d in dirs:
        s, m = plot_h_derived_dos(filename=d+"/proj_dos.out", plot=False)
        score_name.append([s,d, m])

    score_name.sort()
    for s, d, m in score_name:
        print(s,d,"metal" if m else "insulator")

def plot_ebands(filename="./e_bands.in"):
    import matplotlib.pyplot as plt

    start     = False
    exp_count = None
    kpoints   = []
    labels    = {}
    with open(filename) as f:
        for line in f:
            if not start:
                if "K_POINTS" in line:
                    start = True
                continue

            if exp_count is None:
                exp_count = int(line)
                continue

            if "!" in line:
                labels[len(kpoints)] = line.split("!")[-1].strip()

            k = [float(x) for x in line.split()[0:3]]
            kpoints.append(k)

    e_fermi = None
    outf = filename.replace(".in", ".out")
    if os.path.isfile(outf):
        with open(outf) as f:
            for line in f:
                if "the Fermi energy is" in line:
                    e_fermi = float(line.split()[-2])

    if e_fermi is None:
        print("Could not identify the fermi energy!")
        e_fermi = 0.0
    else:
        plt.axhline(0.0, color="black")

    out = parser.extract_bands_out(filename)
    bands = zip(*out["bands"])
    for b in bands:
        plt.plot([x - e_fermi for x in b], marker="+", linestyle="none")

    for i in labels:
        plt.axvline(i, color="black", linestyle=":")
    plt.xticks([i for i in labels], [labels[i] for i in labels])
    plt.ylabel("Energy (eV)")

    plt.show()

def plot_a2f(filename="./a2F.dos1"):
    import matplotlib.pyplot as plt
    
    out = parser.a2f_dos_out(filename)

    # Convert frequencies to Cm^-1
    ws = out["frequencies"] 
    ws = [w*constants.RY_TO_CMM1 for w in ws]

    # Run over mode-resolved eliashberg functions
    total  = None
    i_mode = 1
    while "a2f_mode_{0}".format(i_mode) in out:

        # Get this constribution and the cumulative total
        a2f_mode = out["a2f_mode_{0}".format(i_mode)]
        if total is None: total = [0.0 for a in a2f_mode]
        new_total = [t+a for t,a in zip(total, a2f_mode)]

        # Fill in this contribution
        plt.fill_between(ws, total, new_total)
    
        # Move to next node
        total = new_total
        i_mode += 1

    plt.xlabel("Frequency $\\omega$ (cm$^{-1}$)")
    plt.ylabel("$\\alpha^2F(\\omega)$\n(colored by mode)")
    plt.show()

# Returns the squared distnace between two points
# in reciprocal space, including folding
def delta_q(q, p):
    q = [x - math.floor(x) for x in q]
    p = [x - math.floor(x) for x in p]
    d = [abs(q[i] -p[i]) for i in range(3)]
    return sum(x*x for x in d)

def plot_phonon_dispersion(filename="./ph_interp.freq"):
    import matplotlib.pyplot as plt        

    # Load the structure from an scf.in file
    # in the current directory
    p = params()
    p.load_from_qe_input("scf.in")

    # Parse filename for the phonon modes
    q_point = None
    modes = {}
    with open(filename) as f:
        for line in f:
            try:
                splt = [float(x) for x in line.split()]
                if len(splt) == 3 and all(x <= 1.0 and x >= 0.0 for x in splt):
                    q_point = tuple(splt)
                else:
                    if not q_point in modes:
                        modes[q_point] = []
                    modes[q_point].extend(splt)
            except:
                pass

    # Perform nearest neighbour interpolation
    # along the BZ path, to obtain a phonon dispersion
    dispersion = []
    for q in p["bz_path"]:
        q = tuple(q)

        found = None
        min_dis = float("inf")
        for q2 in modes:
            dis = delta_q(q,q2)
            if dis < min_dis:
                min_dis = dis
                found = q2
        
        dispersion.append(modes[found])

    # Plot the phonon dispersion
    plt.ylabel("Frequency (cm^-1)")
    plt.axhline(0, color="black")
    for mode in zip(*dispersion):
        plt.plot(mode)

    # Label high-symmetry points
    hsp = p["high_symmetry_bz_points"]
    ticks = []
    labels = []
    for pt in hsp:
        plt.axvline(pt, color=[0.5,0.5,0.5])
        ticks.append(pt)
        labels.append(hsp[pt])


    plt.xticks(ticks=ticks, labels=labels, rotation=90)
    plt.show()

def plot_phonon_mode_atoms(filename="./ph_interp.modes"):
    import matplotlib.pyplot as plt
    import numpy as np
    
    freqs = []
    evecs = []
    with open(filename) as f:
        for line in f:
            line = line.strip()

            if line.startswith("("):
                # parse eigenvector for this atom
                line = line.replace("(", "")
                line = line.replace(")", "")
                x = [float(w) for w in line.split()]

                # Record it's magnitude
                x = (x[0]*x[0] + x[2]*x[2] + x[4]*x[4])**0.5
                evecs[len(freqs)-1].append(x)
                continue

            if line.startswith("freq"):
                # parse mode frequency in cm^-1
                w = float(line.split("=")[-1].split("[")[0])
                freqs.append(w)
                evecs.append([])
                continue

    BINS = 1000

    minf  = min(freqs)
    maxf  = max(freqs)
    atoms = len(evecs[0])
    bins  = np.zeros((atoms, BINS+1))

    # Loop over modes
    for i in range(0, len(freqs)):
        f = freqs[i]

        b = (f-minf)/(maxf-minf)
        b = int(b * BINS)

        for j in range(0, atoms):
            bins[j][b] += evecs[i][j]

    # Attempt to parse atom names from scf.in file
    atom_names = []
    started = False
    scf_file = os.path.dirname(filename)+"/scf.in"
    if os.path.isfile(scf_file):
        with open(scf_file) as f:
            for line in f:
                if len(atom_names) >= atoms: break
                if "ATOMIC_POSITIONS" in line:
                    started = True
                    continue
                if started: atom_names.append(line.split()[0])

    # Failed, just use numbers
    if len(atom_names) < atoms:
        atom_names = range(atoms)

    # Colors for the various atoms
    atom_colors = {}
    for n in atom_names:
        if n in atom_colors: continue
        atom_colors[n] = color_cycle(len(atom_colors)-1)

    # Plot atom-resolved eigenvectors
    total = np.zeros(BINS+1)
    freqs = np.linspace(minf, maxf, BINS+1)
    for j in range(0, atoms):

        color = atom_colors[atom_names[j]]
        new_total = total + bins[j]

        plt.fill_between(freqs, total, new_total, 
            label="atom: {0}".format(atom_names[j]),
            color=color)

        atom_colors[atom_names[j]] = [c * 0.8 for c in color]
        total = new_total

    plt.legend(ncol=int(atoms**0.5))
    plt.show()

def plot_tc_vs_smearing(directories=["./"], 
    force_allen_dynes=False, ask=False, plot_relative=False):
    from qet.calculations import tc_from_a2f_allen_dynes

    # If only one directory is given, 
    # include also subdirectories
    if len(directories) == 1:
        for d in os.listdir(directories[0]):
            d = directories[0]+"/"+d
            if not os.path.isdir(d): continue
            directories.append(d)

    if len(directories) == 0:
        print("No directories to plot in tc_vs_smearing!")
        return
    
    if ask:
        txt = input("Would you like to plot {0} (n to skip)?".format(directories[0]))
        if txt == "n":
            return

    import matplotlib.pyplot as plt

    print("first directory: "+directories[0])
    title = directories[0]
    if os.path.isfile(directories[0]+".in"):
        p = params(directories[0]+".in")
        title = p["space_group_name"]+" "+p["stoichiometry_string"]
    plt.suptitle(title)

    if plot_relative:
        print("Plotting relative")
    tc_rel = None

    # Will contain the method used to evaluate Tc
    method = "None"

    for directory in directories:
        if not os.path.isdir(directory):
            print("Directory "+directory+" not found, skipping...")
            continue

        tcs1 = []
        tcs2 = []

        method = "Allen-Dynes"

        # Loop over a2f.dos files
        for f in os.listdir(directory):
            if not "a2F.dos" in f: continue
            f = directory+"/"+f 
            n = int(f.split("a2F.dos")[-1])

            # Attempt to read Tc from elk directory
            elk_dir = directory + "/tc_dos_{0}".format(n)
            if os.path.isdir(elk_dir) and not force_allen_dynes:

                method = "Eliashberg"
                for mu_dir in os.listdir(elk_dir):
                    mu_dir = elk_dir + "/" + mu_dir
                    tc_f   = mu_dir + "/tc.out"
                    if not os.path.isfile(tc_f): continue
                    
                    with open(tc_f) as f:
                        for line in f:
                            if "liashberg" in line:
                                tc = float(line.split()[0])
                                if len(tcs1) == len(tcs2): tcs1.append([n, tc])
                                else: tcs2.append([n, tc])

            # No elk directory, calculate using allen-dynes
            else:

                # Some of the a2F.dos files had Eliashberg, some didn't
                if method != "Allen-Dynes": method = "Mixed"

                try:
                    # Get tc for two different mu* values
                    tcs = tc_from_a2f_allen_dynes(f, mu_stars=[0.1, 0.15])
                    tcs1.append([n, tcs[0.1]])
                    tcs2.append([n, tcs[0.15]])
                except Exception as e:
                    print(e)
                    continue

        # No data => skip
        if len(tcs1) == 0 or len(tcs2) == 0:
            print("No data in "+directory)
            continue

        # Sort by acending smearing
        tcs1.sort()
        tcs2.sort()

        ns, tc1 = zip(*tcs1)
        ns, tc2 = zip(*tcs2)

        # Attempt to find el_ph_sigma in .in files
        el_ph_sigma = None
        for f in os.listdir(directory):
            if not f.endswith(".in"): continue
            f = directory+"/"+f

            with open(f) as of:
                for line in of:
                    if "el_ph_sigma" in line:
                        el_ph_sigma = float(line.split("=")[-1].replace(",",""))
                        break

        # Label with what we're plotting on the x-axis
        if el_ph_sigma is None: 
            plt.xlabel("Smearing number")
        else:
            ns = [el_ph_sigma*n for n in ns]
            plt.xlabel("Smearing width $\\sigma$ (Ry)")

        if plot_relative and tc_rel is None:
            tc_rel = tc1[-1]
        elif not (tc_rel is None):
            diff = tc1[-1] - tc_rel
            tc1 = [t1 - diff for t1 in tc1]
            tc2 = [t2 - diff for t2 in tc2]

        plt.fill_between(ns, tc1, tc2, alpha=0.25, label=directory) 
        plt.plot(ns, tc2, linestyle="none", marker="+")
        plt.ylabel("$T_C$ ({0} with $\\mu^* \\in \; [0.1, 0.15]$)".format(method))

    if plt.ylim()[1] > 1000.0:
        print("Found T_C > 1000 K, rescaling axis")
        print(tc1)
        plt.ylim([-10,1000])

    plt.legend()
    plt.show()

def plot_alch_network(directory=None):
    from qet.alchemy.network import plot_alch_network
    plot_alch_network(directory=directory, pickle=False)

# Run like a program
def main():
    import sys

    ask = "ask" in sys.argv
    rel = "rel" in sys.argv

    # The possible tasks to invoke
    invoke_list = {
        "tc_vs_smearing"    : lambda : plot_tc_vs_smearing(sys.argv[2:], ask=ask, plot_relative=rel),
        "tc_vs_smearing_ad" : lambda : plot_tc_vs_smearing(sys.argv[2:], force_allen_dynes=True, ask=ask, plot_relative=rel),
        "a2f"               : lambda : plot_a2f(sys.argv[2]),
        "alch_network"      : lambda : plot_alch_network(sys.argv[2]),
        "proj_dos"          : lambda : plot_proj_dos(),
        "proj_dos_h"        : lambda : plot_h_derived_dos(),
        "rank_doping"       : lambda : rank_doping(sys.argv[2:]),
        "pdos"              : lambda : plot_pdos(),
        "ebands"            : lambda : plot_ebands(),
        "phonon_atoms"      : lambda : plot_phonon_mode_atoms(),
        "phonon_dispersion" : lambda : plot_phonon_dispersion()
    }

    # Check arguments
    if len(sys.argv) < 2 or not sys.argv[1] in invoke_list:
        print("The first argument to plots.py should be one of:")
        for i in invoke_list:
            print("    "+i)
        return

    invoke_list[sys.argv[1]]()

# Check if we are directly invoking this script, if so
# run plots.py as a program
if __name__ == "__main__": main()
