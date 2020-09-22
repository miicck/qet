# QET
Welcome to the QET (quantum espresso tools) library, for carrying out high-throughput calculations using Quantum Espresso.
QET automates the process of setting up, parallelising and running jobs in Quantum Espresso and streamlines the set up of
calculations in series to build up high-throughput workflows.

# Overview
In order to get an idea of how the various parts of QET work together, we're going to take a look at how the 
`examples/ternary_optimizer.py` example works. This example takes a set of ternary-hydride superconductors and applies
alchemical mutations to construct a network of new materials. This network will grow in such a way as to maximize the fraction
of electronic density of states at the fermi level that originates from the hydrogen atoms. This is closely related to
how good the resulting superconductor is. First, we import the various tools we'll need:

```
from qet.calculations               import relax, proj_dos       # Basic QE calculation types that we will use
from qet.alchemy.network            import alch_network          # The alchemical optimization method that we will use
from qet.alchemy                    import mutations             # The alchemical moves that we will employ
from qet.examples.initial_ternaries import get_initial_ternaries # Initial seed structures for this example
```

The modules in QET are organised in a hierarchical structure, with more specific functionality further down the hierarchy.
For example, `qet.alchemy` contains the basic mutations that all alchemical methods can use, whereas `qet.alchemy.network`
contains tools specific to an alchemical network. Next, we define a our objective function, which uses the basic calculation
types contained within `qet.calculations` to evaluate the hydrogen-derived density of states at the fermi level:

```
def h_dos(structure):

    # Relax the structure
    res = relax(structure).run()
    
    # Update the structure to the relaxed version
    structure["atoms"]   = res["relaxed atoms"]
    structure["lattice"] = res["relaxed lattice"]

    # Calculate the DOS
    res = proj_dos(structure).run()

    # Get the fraction of the DOS due to hydrogens
    total_dos    = 0.0
    hydrogen_dos = 0.0

    # Loop over atoms and atomic wavefunctions
    for atom_number in res["PDOS (fermi energy)"]:
        for wfc_number in res["PDOS (fermi energy)"][atom_number]:

            # Get the type of atom and the contribution to the DOS at E_F
            atom_name = res["PDOS atom names"][atom_number][wfc_number]
            dos_ef    = res["PDOS (fermi energy)"][atom_number][wfc_number]

            # Accumulate the total and hydrogen DOS
            total_dos += dos_ef
            if atom_name == "H":
                hydrogen_dos += dos_ef

    return hydrogen_dos/total_dos
```

This function takes an argument `strcuture` which will be of the `parameters` type defined in `qet.params`. This
`parameters` object defines all of the various settings used by Quantum Espresso, including the crystal structure
of the material in question. It's state is defined by a minimal structure and parameter set and it exposes many useful
methods for accessing that information and derived information. Using these methods, we setup a constaint on the 
types of structures that we are willing to consider. In this case, this consists of only ternary structures with 
up to 6 non-hydrogen elements and between 1 and 8 hydrogens per non-hydrogen:
```
def is_valid(structure):
    
    MAX_NON_H = 6
    MIN_H_PER_NON_H = 1
    MAX_H_PER_NON_H = 8

    species = structure["species"]
    if len(species) != 3: return False # Only ternaries

    counts = structure["atom_counts"]
    non_hs = sum([counts[e] for e in counts if not e == "H"])

    # Limit on non-hydrogens
    if non_hs > MAX_NON_H : return False

    # Limits on number of hydrogens
    if counts["H"] > MAX_H_PER_NON_H * non_hs: return False
    if counts["H"] < MIN_H_PER_NON_H * non_hs: return False

    return True
```

Finally, we construct an alchemical network and expand it in order to maximize the hydrogen-derived density of
states at the fermi level:
```
# Load/create the network
nw = alch_network("network")

# Add the seed structures
for structure in get_initial_ternaries():
    nw.create_vertex(structure)

# Setup the function to minimize
to_min = lambda s : -h_dos(s)

# Optimize the network
for n in range(0, 100):
    nw.expand_to_minimize(to_min, mutations, is_valid=is_valid)
```

This will create a set of subdirectories, each continaing a new structure and the corresponding files from the calculation
of the hydrogen-derived density of states. As it progresses, the structures with the highest hydrogen-DOS will be mutated
to create new structures. Many calculations can be ran simultaneously in the parent directory and will each pick a different
structure to expand, this allows many nodes of a supercomputer to used on the same alchemical network.



