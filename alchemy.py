from qet.elements import atom_substitutions, elements
from qet.params   import parameters
from qet.logs     import log
import random, os

class alch_structure:

    # Create a structure with the given crystal
    # lattice and atoms (in fractional coordinates)
    def __init__(self, lattice, atoms, 
        check_valid=None,
        mutation_type_prob=None,
        mutating=None):
        
        self.params = parameters()
        self.params["lattice"] = lattice
        self.params["atoms"]   = atoms

        # We're mutating a structure, copy various
        # parameters from it
        if not mutating is None:
            check_valid        = mutating.check_valid
            mutation_type_prob = mutating.mutation_type_prob

        # If no validation function is supplied
        # assume all structures are valid
        if check_valid is None:
            self.check_valid = lambda s : True
        else:
            self.check_valid = check_valid

        # If no function specifying mutation type
        # probabilites is supplied, assume equal
        # probability mutation types
        if mutation_type_prob is None:
            self.mutation_type_prob = lambda mt : 1.0
        else:
            self.mutation_type_prob = mutation_type_prob

        # Generate a unique name for this alch_structure
        self.name    = self.params["stochiometry_string"]
        self.version = 1
        for f in os.listdir("."):
            if f.startswith(self.name):
                self.version += 1

        self.name += "_v_{0}".format(self.version)


    # Propose a substitution for atom e
    def get_substitute(self, e):

        # Choose the atom based on the atom substitution matrix
        options = atom_substitutions[e]
        rnd = random.randint(0, sum([options[o] for o in options]))
        tot = 0
        for o in options:
            tot += options[o]
            if tot > rnd: return o

        raise Exception("No substitute atom found!")


    # Propose a mutation type
    def propose_mutation(self):
        
        mutations = [
            "replace_random_atom_type",  # Substitues a random atom type (e.g H2O -> F2O)
            "remove_random_atom",        # Remove a random atom (e.g H2O -> HO)
            "dupe_random_atom",          # Duplicate a random atom (e.g H2O -> H2O2)
            "shuffle_atoms",             # Shuffle the locations of the atoms
        ]

        # Choose a mutation type based on mutation_type_prob
        # (automatically noramlizes the probabilities)
        probs = [self.mutation_type_prob(mt) for mt in mutations]
        rnd   = random.uniform(0, sum(probs))
        tot   = 0.0
        for i, p in enumerate(probs):
            tot += p
            if tot > rnd: return mutations[i]
        
        raise Exception("Proposing mutation type failed!")

    # Return a mutated alch_structure, that doesn't
    # necassarily obey self.check_valid
    def random_mutation(self):

        mutation = self.propose_mutation()

        if mutation == "replace_random_atom_type":

            # Choose the type to be replaced
            to_replace = self.params["species"]
            to_replace = to_replace[random.randrange(len(to_replace))][0]

            # Choose the type to replace it with
            sub = self.get_substitute(to_replace)

            # Create a new set of atoms with the replacement
            new_atoms = self.params["atoms"].copy()
            for i, a in enumerate(new_atoms):
                if a[0] == to_replace: new_atoms[i][0] = sub

            log("replacing {0} with {1} in {2}".format(to_replace, sub, self.name), "alchemy.log")
            return alch_structure(self.params["lattice"], new_atoms, mutating=self)

        if mutation == "remove_random_atom":

            # Dont remove the last atom
            if len(self.params["atoms"]) == 1:
                return None

            # Copy the atom list and remove a random atom from the result
            new_atoms = self.params["atoms"].copy()
            i_rem     = random.randrange(len(new_atoms))
            log("removing atom {0} {1} in {2}".format(i_rem, new_atoms[i_rem], self.name), "alchemy.log")
            del new_atoms[i_rem]
            return alch_structure(self.params["lattice"], new_atoms, mutating=self) 

        if mutation == "dupe_random_atom":
            
            # Copy a random atom
            new_atoms = self.params["atoms"].copy()
            i_dupe    = random.randrange(len(new_atoms))
            new_atom  = new_atoms[i_dupe].copy()

            # Displace it by a gaussian
            for j in range(3):
                new_atom[1][j] += random.gauss(0, 0.1)
            new_atoms.append(new_atom)

            fs = "duplicating atom {0} in {3}\n    Original  {1}\n    Duplicate {2}"
            log(fs.format(i_dupe, new_atoms[i_dupe], new_atom, self.name), "alchemy.log")
            return alch_structure(self.params["lattice"], new_atoms, mutating=self) 

        if mutation == "shuffle_atoms":

            # No point shuffling a single atom
            if len(self.params["atoms"]) == 1:
                return None
            
            # Shuffle the atoms into each others locations
            new_atoms = self.params["atoms"].copy()
            new_pos   = [a[1] for a in new_atoms]
            random.shuffle(new_pos)
            for i in range(len(new_pos)): new_atoms[i][1] = new_pos[i]

            log("shuffling atoms in {0}".format(self.name), "alchemy.log")
            return alch_structure(self.params["lattice"], new_atoms, mutating=self) 

        raise Exception("Unknown mutation: "+mutation)

    # Return a mutation that obeys self.check_valid
    def mutate(self):
        
        while True:
            prop = self.random_mutation()
            if prop is None: continue
            if self.check_valid(prop): return prop

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

alc = alch_structure(start_lattice, start_atoms)
for n in range(0,100):
    alc = alc.mutate()
