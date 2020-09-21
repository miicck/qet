from qet.calculations               import relax, proj_dos       # Basic QE calculation types that we will use
from qet.alchemy.network            import alch_network          # The alchemical optimization method that we will use
from qet.alchemy                    import mutations             # The alchemical moves that we will employ
from qet.examples.initial_ternaries import get_initial_ternaries # Initial seed structures for this example

def minus_h_dos(structure):

    # Relax the structure
    res = relax(structure).run()
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
            if atom_name.strip().lower() == "h":
                hydrogen_dos += dos_ef

    return -hydrogen_dos/total_dos

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

def maximize_hdos(muts, init_structures=None):

    # Load/create the network
    nw = alch_network("network")

    # Add the seed structures
    if init_structures is None: init_structures = get_initial_ternaries()
    for structure in init_structures:
        nw.create_vertex(structure)

    # Optimize the network
    for n in range(0, 100):
        nw.expand_to_minimize(minus_h_dos, muts, is_valid=is_valid)

def fd3m_search():
    
    # General search around the Fd3m structure
    init = get_initial_ternaries()
    maximize_hdos(mutations.all_mutations, init_structures = [init[0]])

def general_search():

    # Search using all mutations available
    maximize_hdos(mutations.all_mutations)

def stoichiometry_search():

    # Allow stoichiometry modifications only
    maximize_hdos([mutations.substitute_random_species,
                   mutations.remove_random_atom,
                   mutations.duplicate_random_atom])

def substitution_search():

    # Allow substitutions only, dont modify the number of atoms
    maximize_hdos([mutations.substitute_random_species])
