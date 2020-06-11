from qet.calculations               import relax, proj_dos
from qet.alchemy.network            import alch_network
from qet.alchemy                    import mutations
from qet.examples.initial_ternaries import get_initial_ternaries
from qet.test.test_systems          import lah10_fm3m

def is_valid(structure):
    
    MAX_NON_H = 6
    MIN_H_PER_NON_H = 1
    MAX_H_PER_NON_H = 16

    species = structure["species"]
    if len(species) > 5: return False # Maximum of 4 non-hydrogen species
    if len(species) < 2: return False # At least binary

    counts = structure["atom_counts"]
    non_hs = sum([counts[e] for e in counts if not e == "H"])

    # Limit on non-hydrogens
    if non_hs > MAX_NON_H : return False

    # Limits on number of hydrogens
    if counts["H"] > MAX_H_PER_NON_H * non_hs: return False
    if counts["H"] < MIN_H_PER_NON_H * non_hs: return False

    return True

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

def substitute_species(s):
    # Restricted set of elements to search
    restricted_elms = ["B","N","H","La","C","Ga","Pt","Au"]
    return mutations.substitute_random_species(s, restricted_elms)

def substitute_atom(s):
    # Restricted set of elements to search
    restricted_elms = ["B","N","H","La","C","Ga","Pt","Au"]
    return mutations.substitute_random_atom(s, restricted_elms)

def search():
    
    # The allowed mutations
    muts = [mutations.remove_random_atom, 
            mutations.duplicate_random_atom, 
            substitute_species,
            substitute_atom]

    lah10_fm3m["pseudo_dirs"] = [
         "/home/mjh261/rds/rds-t2-cs084/pseudopotentials/gbrv", 
         "/home/mjh261/rds/rds-t2-cs084/pseudopotentials/ssp_efficiency", 
         "/home/mjh261/rds/rds-t2-cs084/pseudopotentials/qe_website"]

    # Maximize the hydrogen density of states, starting with LaH10 Fm-3m
    maximize_hdos(muts, init_structures = [lah10_fm3m])
