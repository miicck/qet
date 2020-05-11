import qet.parser        as parser
import qet.constants     as constants
import matplotlib.pyplot as plt
import os

# Plot the density of states calculated
# with output file given by filename
def plot_dos(filename="./dos.out"):
    
    # Parse the DOS
    out = parser.dos_out(filename)

    plt.plot(out["DOS energies"], out["DOS (energy)"], color="blue", label="DOS ($E$)")
    plt.axvline(out["fermi energy"], color="red", label="Fermi energy")
    plt.axhline(out["DOS (E_F)"], color="green", label="DOS ($E_f$)")
    plt.legend()
    plt.show()

# Plot the projected density of states calculated
# with output file given by filename
def plot_proj_dos(filename="./proj_dos.out"):
    
    # Parse the projected density of states
    out = parser.proj_dos_out(filename)

    dos_plot = plt.subplot(211)
    df_plot  = plt.subplot(212)
    
    labels = []
    dfs    = []

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

            # Plot this projected DOS
            label = "Atom {0} ({1}) wfn {2} ({3})"
            label = label.format(atom_num, an, wfn_num, wf)
            dos_plot.plot(es,ps,label=label)

            labels.append(label)
            dfs.append(df)

    dos_plot.axvline(out["fermi energy"], linestyle=":", label="Fermi energy")
    dos_plot.set_xlabel("Energy (Ry)")
    dos_plot.set_ylabel("PDOS")
    dos_plot.set_title("PDOS")
    dos_plot.legend(ncol=int(len(labels)**0.5))

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

def plot_a2f(filename="./a2F.dos1"):
    
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

def plot_tc_vs_smearing(directories=["./"], force_allen_dynes=False):
    from qet.calculations import tc_from_a2f_allen_dynes

    # Ensure primrary grids are treated first
    def sort_key(d):
        if "primary" in d: return "0"
        return d
    directories.sort(key=sort_key)

    tc_offset = None
    ax_normal = plt.subplot(221)
    ax_offset = plt.subplot(222)

    for directory in directories:

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
                if method != "Allen-Dynes": method = "Mixed"
                try:
                    # Get tc for two different mu* values
                    tcs = tc_from_a2f_allen_dynes(f, mu_stars=[0.1, 0.15])
                    tcs1.append([n, tcs[0.1]])
                    tcs2.append([n, tcs[0.15]])
                except: continue

        if len(tcs1) == 0 or len(tcs2) == 0:
            continue

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

        if el_ph_sigma is None: 
            plt.xlabel("Smearing number")
        else:
            ns = [el_ph_sigma*n for n in ns]
            plt.xlabel("Smearing width $\\sigma$ (Ry)")

        ax_normal.fill_between(ns, tc1, tc2, alpha=0.5, label=directory) 

        if tc_offset is None: 
            tc_offset = tc1[-1]
        else:
            dt = tc1[-1] - tc_offset
            tc1 = [t-dt for t in tc1]
            tc2 = [t-dt for t in tc2]

        ax_offset.fill_between(ns, tc1, tc2, alpha=0.5, label=directory)


    for ax in [ax_offset, ax_normal]:
        ax.set_ylabel("$T_C$ ({0} with $\\mu^* \\in \; [0.1, 0.15]$)".format(method))
        ax.legend()

        if ax.get_ylim()[1] > 1000.0:
            print("Found T_C > 1000 K, rescaling axis")
            ax.set_ylim([0, 1000.0])

    plt.show()

def plot_alch_network(directory=None):
    from qet.alchemy.network import plot_alch_network
    plot_alch_network(directory=directory, pickle=False)

# Run like a program
def main():
    import sys

    # The possible tasks to invoke
    invoke_list = {
        "tc_vs_smearing"    : lambda : plot_tc_vs_smearing(sys.argv[2:]),
        "tc_vs_smearing_ad" : lambda : plot_tc_vs_smearing(sys.argv[2:], force_allen_dynes=True),
        "a2f"               : lambda : plot_a2f(sys.argv[2]),
        "alch_network"      : lambda : plot_alch_network(sys.argv[2]),
        "proj_dos"          : lambda : plot_proj_dos(),
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
