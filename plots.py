import qet.parser        as parser
import qet.constants     as constants
import os

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
        plt.plot([x - e_fermi for x in b])

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

    # Plot atom-resolved eigenvectors
    total = np.zeros(BINS+1)
    freqs = np.linspace(minf, maxf, BINS+1)
    for j in range(0, atoms):
        new_total = total + bins[j]
        plt.fill_between(freqs, total, new_total, label="atom: {0}".format(atom_names[j]))
        total = new_total

    plt.legend()
    plt.show()

def plot_tc_vs_smearing(directories=["./"], force_allen_dynes=False, ask=False):
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
    plt.suptitle(directories[0])

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

    # The possible tasks to invoke
    invoke_list = {
        "tc_vs_smearing"    : lambda : plot_tc_vs_smearing(sys.argv[2:], ask="ask" in sys.argv),
        "tc_vs_smearing_ad" : lambda : plot_tc_vs_smearing(sys.argv[2:], force_allen_dynes=True),
        "a2f"               : lambda : plot_a2f(sys.argv[2]),
        "alch_network"      : lambda : plot_alch_network(sys.argv[2]),
        "proj_dos"          : lambda : plot_proj_dos(),
        "pdos"              : lambda : plot_pdos(),
        "ebands"            : lambda : plot_ebands(),
        "phonon_atoms"      : lambda : plot_phonon_mode_atoms()
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
