from qet.elements     import atom_substitutions, elements
from qet.params       import parameters
from qet.calculations import relax
from qet.logs         import log
import random
import numpy as np

MAX_HYDROGENS    = 32
MIN_HYDROGENS    = 4
MAX_NON_HYDROGEN = 2
MIN_NON_HYDROGEN = 1

class ternary:

    # Create a structure with the given 
    # crystal lattice and atoms given
    # in fractional coordinates
    def __init__(self, lattice, atoms):
        
        self.params            = parameters()
        self.params["lattice"] = lattice
        self.params["atoms"]   = atoms
        self.shake_count       = 0

        # Count the numbers of each atom
        self.atom_counts = {}
        for a in self.params["atoms"]:
            if not a[0] in self.atom_counts: 
                self.atom_counts[a[0]]  = 1
            else: 
                self.atom_counts[a[0]] += 1
        if len(self.atom_counts) > 3:
            raise Exception("Not a ternary!")

        # Count the elements in e1,n1,e2,n2,nh form
        self.n1 = None
        self.n2 = None
        self.nh = self.atom_counts["H"]
        for a in self.atom_counts:
            if a == "H": continue
            if self.n1 is None:
                self.e1 = a
                self.n1 = self.atom_counts[a]
                continue
            if self.n2 is None:
                self.e2 = a
                self.n2 = self.atom_counts[a]
                continue
            raise Exception("Types of atoms > 3!")

        # Make sure e1 is always the heaviest atom
        if elements[self.e1]["atomic number"] < elements[self.e2]["atomic number"]:

            # Wrong way around, swap
            tmpe = self.e1
            tmpn = self.n1
            self.e1 = self.e2
            self.n1 = self.n2
            self.e2 = tmpe
            self.n2 = tmpn

        # Get my unique name
        self.name = "{0}_{1}_{2}_{3}_H_{4}_shake_{5}".format(
            self.e1, self.n1, self.e2, self.n2, self.nh, self.shake_count)

        log("Created "+self.name, "alchemy.log")

    # Propose a substitution for atom e
    def get_substitute(self, e):

        # Get the options according to
        # the atom_substitution matrix
        options = dict(atom_substitutions[e])
        if "H" in options:
            del options["H"] # Don't sub hydrogen in

        # Normalize them to a probability
        tot = 0.0
        for o in options: tot += options[o]
        for o in options: options[o] /= tot

        # Pick according to that probability
        rnd = np.random.random()
        tot = 0.0
        for o in options:
            tot += options[o]
            if tot > rnd:
                return o

        raise Exception("No substitute atom found!")

    # Get the probabilities of various move types
    def move_probs(self):

        # Get something proportional to the probabilities
        # of making various kinds of moves
        move_type_probs = {

            # Change e1, e2 with reasonable probability
            "e1_replace"   : 1.0,
            "e2_replace"   : 1.0,

            # Change nh with reasonable probability
            "nh_increase"  : 1.0 if self.nh < MAX_HYDROGENS else 0.0,
            "nh_decrease"  : 1.0 if self.nh > MIN_HYDROGENS else 0.0,

            # Change n1, n2 with low probability
            "n1_increase"  : 0.1 if self.n1 < MAX_NON_HYDROGEN else 0.0,
            "n2_increase"  : 0.1 if self.n2 < MAX_NON_HYDROGEN else 0.0,

            # Vaguley prefer decreasing n1, n2 for speed
            "n1_decrease"  : 0.2 if self.n1 > MIN_NON_HYDROGEN else 0.0,
            "n2_decrease"  : 0.2 if self.n2 > MIN_NON_HYDROGEN else 0.0,
        }

        # Normalize the above
        tot = 0.0
        for mt in move_type_probs: tot += move_type_probs[mt]
        for mt in move_type_probs: move_type_probs[mt] /= tot

        return move_type_probs

    # Suggest a perturbed ternary structure
    # including alchemical pertubations
    def propose_move(self):
    
        # Pick a move type according to
        # probabilities self.move_probs 
        probs = self.move_probs()
        rnd   = np.random.random()
        tot   = 0.0
        for mt in probs:
            tot += probs[mt]
            if tot > rnd:
                move = mt
                break

        if move == "e1_replace" or move == "e2_replace":
            
            # Replace a non-hydrogen atom
            to_replace = self.e1 if move == "e1_replace" else self.e2
            sub = self.get_substitute(to_replace)
            log("replacing "+to_replace+" with "+sub, "alchemy.log")

            # Create the new list of atoms with to_replace replaced
            new_atoms = list(self.params["atoms"])
            for i, a in enumerate(new_atoms):
                if a[0] == to_replace: new_atoms[i][0] = sub

            return ternary(self.params["lattice"], new_atoms)

        if move == "nh_decrease" or move == "n1_decrease" or move == "n2_decrease":

            # Remove a random atom of the given type
            if   move == "nh_decrease": to_rem = "H"
            elif move == "n1_decrease": to_rem = self.e1
            elif move == "n2_decrease": to_rem = self.e2
            else: raise Exception("Unkown atom-removal move!")
            
            new_atoms = list(self.params["atoms"])
            ais = [i for i, a in enumerate(new_atoms) if a[0] == to_rem]
            irem = ais[random.randrange(0, len(ais))]
            log("removing atom {0} {1}".format(irem+1, new_atoms[irem]), "alchemy.log")
            del new_atoms[irem]

            return ternary(self.params["lattice"], new_atoms)

        if 

    def __str__(self):
        return self.params.__str__()


start_lattice = [
    [4.0, 0, 0],
    [0, 4.0, 0],
    [0, 0, 4.0],
]

start_atoms = [
    ["Li", [0,    0,    0   ]],
    ["Mg", [0.5,  0.5,  0.5 ]],
    ["H" , [0.75, 0.5,  0.5 ]],
    ["H" , [0.25, 0.5,  0.5 ]],
    ["H" , [0.5,  0.75, 0.5 ]],
    ["H" , [0.5,  0.25, 0.5 ]],
    ["H" , [0.5,  0.5,  0.75]],
    ["H" , [0.5,  0.5,  0.25]],
]

ter = ternary(start_lattice, start_atoms)
print(ter.propose_move())
