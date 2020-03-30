# 
#   This file builds platonic cage-like ternary hydride (potential)
#   superconductors, performs a relaxation on each and calculates
#   the projected density of states.
#
import os, math
import numpy          as     np
from qet.params       import parameters
from qet.logs         import log
from qet.elements     import elements
from qet.calculations import relax

PSEUDO_DIR = "/home/mjh261/rds/rds-t2-cs084/pseudopotentials/gbrv"

def get_platonic_hydrogens(nh, lattice, centre=[0.5,0.5,0.5], scale=1.0):

    # Scale of the displacement from the centre
    a = 0.33333333333333 * scale * np.linalg.det(lattice)**(1.0/3.0)

    # The centre of the platonic solid in cart coords
    x,y,z = np.matmul(np.array(lattice).T, centre)

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

    # Generate the hydrogens in cartesian coordinates

    if nh == 4:
        # Tetrahedron
        atoms =  [
            ["H", [x+a, y+a, z+a]],
            ["H", [x+a, y-a, z-a]],
            ["H", [x-a, y+a, z-a]],
            ["H", [x-a, y-a, z+a]]
        ]

    elif nh == 6:
        # Octahedron
        atoms =  [
            ["H", [x+a, y,   z]],
            ["H", [x-a, y,   z]],

            ["H", [x,   y+a, z]],
            ["H", [x,   y-a, z]],

            ["H", [x,   y,   z+a]],
            ["H", [x,   y,   z-a]]
        ]

    elif nh == 8:
        # Cube
        atoms =  [
            ["H", [x+a, y+a, z+a]],

            ["H", [x-a, y+a, z+a]],
            ["H", [x+a, y-a, z+a]],
            ["H", [x+a, y+a, z-a]],

            ["H", [x+a, y-a, z-a]],
            ["H", [x-a, y+a, z-a]],
            ["H", [x-a, y-a, z+a]],

            ["H", [x-a, y-a, z-a]],
        ]

    elif nh == 12:
        # Icosahedron
        atoms =  [
            ["H", [x,    y+a,  z+ap]],
            ["H", [x,    y-a,  z+ap]],
            ["H", [x,    y+a,  z-ap]],
            ["H", [x,    y-a,  z-ap]],

            ["H", [x+ap, y,    z+a]],
            ["H", [x-ap, y,    z+a]],
            ["H", [x+ap, y,    z-a]],
            ["H", [x-ap, y,    z-a]],

            ["H", [x+a,  y+ap, z]],
            ["H", [x-a,  y+ap, z]],
            ["H", [x+a,  y-ap, z]],
            ["H", [x-a,  y-ap, z]]
        ]

    elif nh == 20:
        # Dodechahedron
        atoms =  [
            ["H", [x+a,   y+a,   z+a]],

            ["H", [x-a,   y+a,   z+a]],
            ["H", [x+a,   y-a,   z+a]],
            ["H", [x+a,   y+a,   z-a]],

            ["H", [x+a,   y-a,   z-a]],
            ["H", [x-a,   y+a,   z-a]],
            ["H", [x-a,   y-a,   z+a]],

            ["H", [x-a,   y-a,   z-a]],

            ["H", [x,     y+aop, z+ap]],
            ["H", [x,     y-aop, z+ap]],
            ["H", [x,     y+aop, z-ap]],
            ["H", [x,     y-aop, z-ap]],

            ["H", [x+aop, y,     z+ap]],
            ["H", [x-aop, y,     z+ap]],
            ["H", [x+aop, y,     z-ap]],
            ["H", [x-aop, y,     z-ap]],

            ["H", [x+aop, y+ap,  z]],
            ["H", [x-aop, y+ap,  z]],
            ["H", [x+aop, y-ap,  z]],
            ["H", [x-aop, y-ap,  z]]
        ]

    else:
        return None

    # Convert back to fractional coordinates
    linv = np.linalg.inv(np.array(lattice).T)
    for i in range(len(atoms)):
        atoms[i][1] = np.matmul(linv, atoms[i][1])

    return atoms


def get_extended_platonic(n, lattice, centre=[0.5,0.5,0.5]):
    
    atoms = []

    # If neccassary, combine several platonic solids
    # to obtain a particular hydrogen count.
    if nh in [4, 6, 8, 12, 20]:
        atoms.extend(get_platonic_hydrogens(nh, lattice, centre))

    elif nh == 10:
        atoms.extend(get_platonic_hydrogens(4, lattice, centre))
        atoms.extend(get_platonic_hydrogens(6, lattice, centre))

    elif nh == 14:
        atoms.extend(get_platonic_hydrogens(6, lattice, centre))
        atoms.extend(get_platonic_hydrogens(8, lattice, centre))

    elif nh == 16:
        atoms.extend(get_platonic_hydrogens(4,  lattice, centre))
        atoms.extend(get_platonic_hydrogens(12, lattice, centre))

    elif nh == 18:
        atoms.extend(get_platonic_hydrogens(6,  lattice, centre))
        atoms.extend(get_platonic_hydrogens(12, lattice, centre))

    elif nh == 22:
        atoms.extend(get_platonic_hydrogens(12, lattice, centre))
        atoms.extend(get_platonic_hydrogens(6,  lattice, centre))
        atoms.extend(get_platonic_hydrogens(4,  lattice, centre))

    elif nh == 24:
        atoms.extend(get_platonic_hydrogens(20, lattice, centre))
        atoms.extend(get_platonic_hydrogens(4,  lattice, centre, scale=1.5))

    elif nh == 26:
        atoms.extend(get_platonic_hydrogens(20, lattice, centre))
        atoms.extend(get_platonic_hydrogens(6,  lattice, centre))

    elif nh == 28:
        atoms.extend(get_platonic_hydrogens(20, lattice, centre))
        atoms.extend(get_platonic_hydrogens(8,  lattice, centre, scale=1.5))

    elif nh == 30:
        atoms.extend(get_platonic_hydrogens(20, lattice, centre))
        atoms.extend(get_platonic_hydrogens(6,  lattice, centre))
        atoms.extend(get_platonic_hydrogens(4,  lattice, centre, scale=1.5))
        
    else:
        return None

    return atoms


# Fill in hydrogens using platonic solids.
def generate_structure_platonic(e1, e2, n1, n2, nh, volume_guess):

    # Equal n1, n2 case
    if n1 == n2:
        
        # Put e1 at corner, e2 at centre of cell
        atoms   = [
            [e1, [0,   0,   0  ]],
            [e2, [0.5, 0.5, 0.5]]
        ]

        # Use a cubic lattice
        lattice = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
        
    elif [n1,n2] in [[1,2], [2,1]]:
        
        # Work out which we have more of, e1 or e2.
        emost  = e1
        eleast = e2
        if n2 > n1:
            emost  = e2
            eleast = e1

        # Arrange the cell so we have the
        # correct ratio of e1 : e2
        atoms = [
            [eleast, [0.5,  0.5,  0.5 ]],
            [emost,  [0.75, 0.75, 0.75]],
            [emost,  [0.25, 0.25, 0.25]]
        ]

        # Use a lattice with a=b=c alpha=beta=gamma=60
        lattice = [
            [1.0, 0.0,       0.0      ], 
            [0.5, 0.8660254, 0.0      ],
            [0.5, 0.2886751, 0.8164966]]

    else:
        fs = "Don't know how to deal with the n1 = {0}, n2 = {1} case."
        raise Exception(fs.format(n1,n2))

    # Set the lattice volume to the guessed volume
    lattice  = np.array(lattice)
    lat_vol  = np.linalg.det(lattice)
    lattice *= (volume_guess/lat_vol)**(1.0/3.0)

    # Create the hydrogen cage around the central atom
    ph = get_extended_platonic(nh, lattice, centre=[0.5, 0.5, 0.5])
    if ph is None: return None

    atoms.extend(ph)

    par = parameters()
    par["lattice"] = lattice
    par["atoms"]   = atoms
    return par
    
# Returns a parameter set with a structure
# with the given stochiometry
def generate_structure(e1, e2, n1, n2, nh):

    # Estimate the volume based on the covalent radius
    v1 =  elements[e1]["covalent radius"] ** 3
    v2 =  elements[e2]["covalent radius"] ** 3
    vh = elements["H"]["covalent radius"] ** 3
    volume_guess  = n1*v1 + n2*v2 + nh*vh
    volume_guess *= 4.0*3.14159/3.0 # 4/3 pi
    volume_guess *= 5.0 # This factor is from tests

    #return generate_fibonacci_cage_structure(e1, e2, n1, n2, nh, volume_guess)
    return generate_structure_platonic(e1, e2, n1, n2, nh, volume_guess)

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
    title = par["stoichiometry_string"]
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
    if an > 4: continue # We don't consider atomic numbers > 103
    if ra:       continue # No radioactive elements

    non_hydrogen_choices.append(e)


# Construct hydrides with the given
# ratios of the two non-hydrogen elements
allowed_ratios = [[1,1], [2,1], [1,2]]

# The minimum and maximum allowed numbers
# of hydrogens per non-hydrogen
MIN_H_CONTENT = 1
MAX_H_CONTENT = 10

# Loop over pairs of elements
for i in range(0, len(non_hydrogen_choices)):
    e1 = non_hydrogen_choices[i]

    for j in range(0, i):
        e2 = non_hydrogen_choices[j]
    
        # Ensure e1 is the heavier element
        if elements[e2]["atomic number"] > elements[e1]["atomic number"]:
            tmp = e2
            e2  = e1
            e1  = tmp

        # Loop over allowed ratios
        for rat in allowed_ratios:
            num_non_h = sum(rat)

            # Loop over allowed hydrogen content
            for nh in range(
                num_non_h * MIN_H_CONTENT, 
                num_non_h * MAX_H_CONTENT + 1):

                # Run this stochiometry
                run(e1, e2, rat[0], rat[1], nh)
