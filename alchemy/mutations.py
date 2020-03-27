from   qet.alchemy.elements import propose_substitution
from   qet.elements         import elements
from   qet.logs             import log
import random, copy

def substitute_random_species(structure):

    # Choose the type to be replaced
    to_replace = structure["species"]
    to_replace = to_replace[random.randrange(len(to_replace))][0]

    # Choose the type to replace it with
    sub = propose_substitution(to_replace)

    fs = "Replacing {0} with {1} in {2}"
    fs = fs.format(to_replace, sub, structure["stochiometry_string"])
    log(fs, "alchemy.log")

    # Create a new set of atoms with the replacement
    new_structure = structure.copy()
    for i, a in enumerate(new_structure["atoms"]):
        if a[0] == to_replace: new_structure["atoms"][i][0] = sub

    # Try to maintain a sensible volume
    volume_factor = new_structure["covalent_volume"]/structure["covalent_volume"]
    new_structure["lattice"] *= volume_factor ** (1.0/3.0)

    return new_structure

def remove_random_atom(structure):
    
    # Dont remove the last atom
    if len(structure["atoms"]) == 1:
        return None

    # Copy the atom list and remove a random atom from the result
    new_structure = structure.copy()
    i_rem         = random.randrange(len(structure["atoms"]))
    atom_removed  = new_structure["atoms"][i_rem]

    fs = "Removing atom {0} in {1}\n    Removed {2} @ {3:8.6f} {4:8.6f} {5:8.6f}"
    fs = fs.format(i_rem, structure["stochiometry_string"], atom_removed[0], *atom_removed[1])
    log(fs, "alchemy.log")

    del new_structure["atoms"][i_rem]
    return new_structure


def duplicate_random_atom(structure):
    
    # Copy a random atom
    new_structure = structure.copy()
    new_atoms     = new_structure["atoms"]
    i_dupe        = random.randrange(len(new_atoms))
    new_atom      = copy.deepcopy(new_atoms[i_dupe])

    fs  = "Duplicating atom {0} in {1}".format(i_dupe, structure["stochiometry_string"])
    fs += "\n    Duplicated {0} @ {1:8.6f} {2:8.6f} {3:8.6f}".format(new_atom[0], *new_atom[1])

    # Displace it by a gaussian
    for j in range(3):
        new_atom[1][j] += random.gauss(0, 0.1)
    new_atoms.append(new_atom)

    fs += "\n    New atom   {0} @ {1:8.6f} {2:8.6f} {3:8.6f}".format(new_atom[0], *new_atom[1])
    log(fs, "alchemy.log")

    return new_structure

def shuffle_atoms(structure):
    
    # No point shuffling a single atom
    if len(structure["atoms"]) == 1:
        return None

    # Shuffle the atoms into each others locations
    new_structure = structure.copy()
    new_atoms     = new_structure["atoms"]
    new_pos       = [a[1] for a in new_atoms]
    random.shuffle(new_pos)
    for i in range(len(new_pos)): new_atoms[i][1] = new_pos[i]

    log("Shuffled atoms in {0}".format(structure["stochiometry_string"]), "alchemy.log")

    return new_structure

