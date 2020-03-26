from qet.alchemy.elements       import propose_substitution
from qet.logs                   import log
import copy, random

def substitute_random_atom_type(stucture):

    # Choose the type to be replaced
    to_replace = structure["species"]
    to_replace = to_replace[random.randrange(len(to_replace))][0]

    # Choose the type to replace it with
    sub = propose_substitution(to_replace)

    # Create a new set of atoms with the replacement
    new_structure = copy.deepcopy(structure)
    for i, a in enumerate(new_structure["atoms"]):
        if a[0] == to_replace: new_structure["atoms"][i][0] = sub

    return new_structure

def remove_random_atom(structure):
    
    # Dont remove the last atom
    if len(structure["atoms"]) == 1:
        return None

    # Copy the atom list and remove a random atom from the result
    new_structure = copy.deepcopy(structure)
    i_rem         = random.randrange(len(structure["atoms"]))

    del new_structure["atoms"][i_rem]
    return new_structure


def duplicate_random_atom(structure):
    
    # Copy a random atom
    new_structure = copy.deepcopy(structure)
    new_atoms     = new_structure["atoms"]
    i_dupe        = random.randrange(len(new_atoms))
    new_atom      = copy.deepcopy(new_atoms[i_dupe])

    # Displace it by a gaussian
    for j in range(3):
        new_atom[1][j] += random.gauss(0, 0.1)
    new_atoms.append(new_atom)

    return new_structure

def shuffle_atoms(structure):
    
    # No point shuffling a single atom
    if len(structure["atoms"]) == 1:
        return None

    # Shuffle the atoms into each others locations
    new_structure = copy.deepcopy(structure)
    new_atoms     = new_structure["atoms"]
    new_pos       = [a[1] for a in new_atoms]
    random.shuffle(new_pos)
    for i in range(len(new_pos)): new_atoms[i][1] = new_pos[i]

    return new_structure

