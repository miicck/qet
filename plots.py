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

def plot_tc_vs_smearing(directory="./"):
    from qet.calculations import tc_from_a2f_allen_dynes

    tcs1 = []
    tcs2 = []

    # Loop over a2f.dos files
    for f in os.listdir(directory):
        if not "a2F.dos" in f: continue
        f = directory+"/"+f 

        try:
            # Get tc for two different mu* values
            tcs        = tc_from_a2f_allen_dynes(f, mu_stars=[0.1, 0.15])
            n          = int(f.split("a2F.dos")[-1])
            tcs1.append([n, tcs[0.1]])
            tcs2.append([n, tcs[0.15]])

        except: continue

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

    plt.fill_between(ns, tc1, tc2, alpha=0.5)
    plt.ylabel("$T_C$ (allen-dynes with $\\mu^* \\in [0.1, 0.15]$)")
    plt.show()
