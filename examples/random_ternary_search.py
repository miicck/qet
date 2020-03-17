# 
#   This file builds random cage-like ternary hydride (potential)
#   superconductors, performs a relaxation on each and calculates
#   the projected density of states.
#
import os
from qet.params       import parameters
from qet.logs         import log
from qet.elements     import elements
from qet.calculations import relax
from scipy.optimize   import minimize

PSEUDO_DIR = "/home/mjh261/rds/rds-t2-cs084/pseudopotentials/gbrv"

def get_platonic_hydrogens(nh, scale=0.3333333333):

    # Scale of the displacement from the centre
    a = scale

    # Golden ratio
    p = (1.0 + 5.0**0.5)/2.0

    # Normalize the scale by the length of
    # the pertubation for the given platonic type
    if   nh == 4 : a /= 3.0 ** 0.5
    elif nh == 6 : a /= 1.0
    elif nh == 8 : a /= 3.0 ** 0.5
    elif nh == 12: a /= (1 + p**2)**0.5
    elif nh == 20: a /= 3.0 ** 0.5

    ap  = a * p
    aop = a / p

    if nh == 4:
        # Tetrahedron
        return [
            ["H", [0.5+a, 0.5+a, 0.5+a]],
            ["H", [0.5+a, 0.5-a, 0.5-a]],
            ["H", [0.5-a, 0.5+a, 0.5-a]],
            ["H", [0.5-a, 0.5-a, 0.5+a]]
        ]

    elif nh == 6:
        # Octahedron
        return [
            ["H", [0.5+a, 0.5, 0.5]],
            ["H", [0.5-a, 0.5, 0.5]],

            ["H", [0.5, 0.5+a, 0.5]],
            ["H", [0.5, 0.5-a, 0.5]],

            ["H", [0.5, 0.5, 0.5+a]],
            ["H", [0.5, 0.5, 0.5-a]]
        ]

    elif nh == 8:
        # Cube
        return [
            ["H", [0.5+a, 0.5+a, 0.5+a]],

            ["H", [0.5-a, 0.5+a, 0.5+a]],
            ["H", [0.5+a, 0.5-a, 0.5+a]],
            ["H", [0.5+a, 0.5+a, 0.5-a]],

            ["H", [0.5+a, 0.5-a, 0.5-a]],
            ["H", [0.5-a, 0.5+a, 0.5-a]],
            ["H", [0.5-a, 0.5-a, 0.5+a]],

            ["H", [0.5-a, 0.5-a, 0.5-a]],
        ]

    elif nh == 12:
        # Icosahedron
        return [
            ["H", [0.5, 0.5+a, 0.5+ap]],
            ["H", [0.5, 0.5-a, 0.5+ap]],
            ["H", [0.5, 0.5+a, 0.5-ap]],
            ["H", [0.5, 0.5-a, 0.5-ap]],

            ["H", [0.5+ap, 0.5, 0.5+a]],
            ["H", [0.5-ap, 0.5, 0.5+a]],
            ["H", [0.5+ap, 0.5, 0.5-a]],
            ["H", [0.5-ap, 0.5, 0.5-a]],

            ["H", [0.5+a, 0.5+ap, 0.5]],
            ["H", [0.5-a, 0.5+ap, 0.5]],
            ["H", [0.5+a, 0.5-ap, 0.5]],
            ["H", [0.5-a, 0.5-ap, 0.5]]
        ]

    elif nh == 20:
        # Dodechahedron
        return [
            ["H", [0.5+a, 0.5+a, 0.5+a]],

            ["H", [0.5-a, 0.5+a, 0.5+a]],
            ["H", [0.5+a, 0.5-a, 0.5+a]],
            ["H", [0.5+a, 0.5+a, 0.5-a]],

            ["H", [0.5+a, 0.5-a, 0.5-a]],
            ["H", [0.5-a, 0.5+a, 0.5-a]],
            ["H", [0.5-a, 0.5-a, 0.5+a]],

            ["H", [0.5-a, 0.5-a, 0.5-a]],

            ["H", [0.5, 0.5+aop, 0.5+ap]],
            ["H", [0.5, 0.5-aop, 0.5+ap]],
            ["H", [0.5, 0.5+aop, 0.5-ap]],
            ["H", [0.5, 0.5-aop, 0.5-ap]],

            ["H", [0.5+aop, 0.5, 0.5+ap]],
            ["H", [0.5-aop, 0.5, 0.5+ap]],
            ["H", [0.5+aop, 0.5, 0.5-ap]],
            ["H", [0.5-aop, 0.5, 0.5-ap]],

            ["H", [0.5+aop, 0.5+ap, 0.5]],
            ["H", [0.5-aop, 0.5+ap, 0.5]],
            ["H", [0.5+aop, 0.5-ap, 0.5]],
            ["H", [0.5-aop, 0.5-ap, 0.5]]
        ]

    # Cant make a platonic solid from this numbner
    msg = "Tried to create a platonic solid with {0} verticies!"
    raise ValueError(msg.format(nh))


# Fill in hydrogens using platonic solids.
def generate_structure_platonic(e1, e2, n1, n2, nh):

    # Equal n1, n2 case
    if n1 == n2:
        
        # Put e1 at corner, e2 at centre of cell
        atoms   = [
            [e1, [0,   0,   0  ]],
            [e2, [0.5, 0.5, 0.5]]
        ]

        # If neccassary, combine several platonic solids
        # to obtain a particular hydrogen count.
        if nh in [4, 6, 8, 12, 20]:
            atoms.extend(get_platonic_hydrogens(nh))

        elif nh == 10:
            atoms.extend(get_platonic_hydrogens(4))
            atoms.extend(get_platonic_hydrogens(6))

        elif nh == 14:
            atoms.extend(get_platonic_hydrogens(6))
            atoms.extend(get_platonic_hydrogens(8))

        elif nh == 16:
            atoms.extend(get_platonic_hydrogens(4))
            atoms.extend(get_platonic_hydrogens(12))

        elif nh == 18:
            atoms.extend(get_platonic_hydrogens(6))
            atoms.extend(get_platonic_hydrogens(12))

        else:
            return None

        par = parameters()
        par["lattice"] = [[10,0,0],[0,10,0],[0,0,10]]
        par["atoms"]   = atoms
        return par
        
    else:
        return None
    
# Returns a parameter set with a structure
# with the given stochiometry
def generate_structure(e1, e2, n1, n2, nh):
    return generate_structure_platonic(e1, e2, n1, n2, nh)

# Run a particular ternary hydride
# with stoichiometry e1_n1 e2_n2 H_nh
count = 0
def run(e1, e2, n1, n2, nh):

    # Check we have the pseudopotentials
    if not os.path.isfile(PSEUDO_DIR+"/"+e1+".UPF"): return
    if not os.path.isfile(PSEUDO_DIR+"/"+e2+".UPF"): return

    # Generate the structure, skip if we couldn't
    par = generate_structure(e1, e2, n1, n2, nh)
    if par is None: return

    # Count the structures
    global count
    count += 1
    
    # Get the title of this structure
    fs    = "{0}_{1}_{2}_{3}_H_{4}"
    title = fs.format(e1, n1, e2, n2, nh)
    print(count, title)

    par["pseudo_dir"] = PSEUDO_DIR
    rlx = relax(par)
    with open(title+".relax.in","w") as f:
        f.write(rlx.gen_input_file())

# Choose the elements that can be used to
# build the ternary hydride
non_hydrogen_choices = []
for e in elements:

    an = elements[e]["atomic number"]
    ra = elements[e]["radioactive"]

    if an < 2:   continue # Not hydrogen
    if an > 5: continue # We don't consider atomic numbers > 103
    if ra:       continue # No radioactive elements

    non_hydrogen_choices.append(e)


# Construct hydrides with the given
# ratios of the two non-hydrogen elements
allowed_ratios = [[1,1], [2,1], [1,2]]

# The minimum and maximum allowed numbers
# of hydrogens per non-hydrogen
MIN_H_CONTENT = 1
MAX_H_CONTENT = 8

# Loop over pairs of elements
for i in range(0, len(non_hydrogen_choices)):
    e1 = non_hydrogen_choices[i]

    for j in range(0, i):
        e2 = non_hydrogen_choices[j]

        # Loop over allowed ratios
        for rat in allowed_ratios:
            num_non_h = sum(rat)

            # Loop over allowed hydrogen content
            for nh in range(
                num_non_h * MIN_H_CONTENT, 
                num_non_h * MAX_H_CONTENT + 1):

                # Run this stochiometry
                run(e1, e2, rat[0], rat[1], nh)
