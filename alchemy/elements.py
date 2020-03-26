import os
from qet.elements import data

# Atom-substitution counts based on the ICSD
# (useful for alchemical optimization)
#
# from
#      https://tddft.org/bmg/files/data/pettifor/raw_data/substitution.dat
#
# based on the paper
#      https://doi.org/10.1088%2F1367-2630%2F18%2F9%2F093011
#

# Parse from substitution.dat
atom_substitutions = {}
base_dir = os.path.dirname(os.path.abspath(__file__))
with open(base_dir+"/substitution.dat") as f:
    i = -1
    for line in f:
        i += 1
        if i < 1: continue # First row/colum is padded with zeros
        ints = [int(w) for w in line.split()]

        # Get substitute atoms, sorted by frequency
        subs = {data[j-1][2] : c for j, c in enumerate(ints) if c > 0 and j > 0}
        subs = {k: v for k, v in sorted(subs.items(), key=lambda item: -item[1])}

        atom_substitutions[data[i-1][2]] = subs

# Propose a substitution for atom e
def propose_substitution(e):

    # Choose the atom based on the atom substitution matrix
    options = atom_substitutions[e]
    rnd = random.randint(0, sum([options[o] for o in options]))
    tot = 0
    for o in options:
        tot += options[o]
        if tot >= rnd: return o

    print(tot, rnd, sum([options[o] for o in options]))
    raise Exception("No substitute atom found for {0} in {1}!".format(e, options))

