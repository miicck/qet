# 
#   This file builds random cage-like ternary hydride (potential)
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

# Generate points on a sphere following the fibonacci method
def fibonacci_sphere(samples):

    points = []
    offset = 2./samples
    increment = math.pi * (3. - math.sqrt(5.));

    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2);
        r = math.sqrt(1 - pow(y,2))

        phi = (i % samples) * increment

        x = math.cos(phi) * r
        z = math.sin(phi) * r

        points.append([x,y,z])

    return points

def get_atom_hydrogen_cluster(e, nh, centre, radius):
    
    # Place the atom at the centre of
    # the cluster
    centre = np.array(centre)
    cluster = [[e, centre]]

    # Place hydrogens in a rough sphere around the cluster
    sphere_pts = fibonacci_sphere(nh)
    for n in range(nh):
        hpos = centre + radius*np.array(sphere_pts[n])
        cluster.append(["H", hpos])

    return cluster

# Generate fibonacci hydrogen cages around e1 and e2 with radii equal to the
# covalent radius of e1 and e2 respectively. The ratio of the  number 
# of hydrogens on e1 to the number of hydrogens on e2 is equal to the
# ratio of covalent radii squared.
def generate_fibonacci_cage_structure(e1, e2, n1, n2, nh, volume_guess):

    # Create a list of the non-hydrogen elements
    es = []
    for i in range(0, n1): es.append(e1)
    for i in range(0, n2): es.append(e2)

    # Set the number of hydrogens for each element proportional
    # to the covalent radius of the element
    cs   = [elements[e]["covalent radius"] for e in es]
    norm = sum([c**2.0 for c in cs])
    nhs  = [float(nh)*(c**2.0)/norm for c in cs]
    nhs  = [int(round(x)) for x in nhs]

    # If we have too many hydrogens, remove
    # one from the atom with the largest covalent radius
    if sum(nhs) == nh + 1:
        imax       = cs.index(max(cs))
        nhs[imax] -= 1
    
    # If we have too few hydrogens, remove one
    # from the atom with the smallest covalent radius
    if sum(nhs) == nh - 1:
        imin       = cs.index(min(cs))
        nhs[imin] += 1

    if sum(nhs) != nh:
        raise Exception("Incorrect number of hydrogens produed!")

    if n1 + n2 == 2:
        # Place atoms at corner and BCC position
        centres        = [[0,0,0],[0.5,0.5,0.5]]
        cluster_radius = 0.5
    elif n1 + n2 == 3:
        # Place atoms along diagonal of unit cell
        t = 1.0/3.0
        centres        = [[0,0,0],[t,t,t],[2*t,2*t,2*t]]
        cluster_radius = t
    else:
        raise Exception("Unsupported number of non-hydrogen atoms!")

    # Place each element in the cell, accompanied by it's
    # cluster of hydrogens
    atoms = []
    for i, e in enumerate(es):
        cluster = get_atom_hydrogen_cluster(
            e, nhs[i], centres[i], cluster_radius)
        atoms.extend(cluster)
    
    # Put in a cubic cell
    par            = parameters()
    alat           = volume_guess ** (1.0/3.0)
    par["atoms"]   = atoms
    par["lattice"] = [[alat,0,0],[0,alat,0],[0,0,alat]]

    return par

# Generate a ternary structure such that each of the
# non-hydrogen elements sees the same hydrogen environmnet
def generate_structure_unbiased(e1, e2, n1, n2, nh, volume_guess):

    raise Exception("Implementation unfinished!")

    # These cases are not implemented
    if n1 != n2:    return None 
    if nh % 2 != 0: return None

    # Place the two non-hydrogen atoms on
    # two interpenetrating simple-cubic lattices
    atoms = [
        [e1,[0.0,0.0,0.0]],
        [e2,[0.5,0.5,0.5]]
    ]

    # We can only do even numbers of hydrogens
    # because the environment around the first
    # atom must be the same as the environment
    # around the second atom

    if   nh == 2:
        # Hydrogens either side of central atom
        atoms.extend([
            # Around first atom
            ["H",[0.25,0.25,0.25]],
            # Around second atom
            ["H",[0.75,0.75,0.75]]
        ])

    elif nh == 4:
        # The hydrogens are in a diagonal 
        # square around the central atom
        atoms.extend([
            # Around first atom
            ["H",[0.25,0.25,0.25]],
            ["H",[0.25,0.25,0.75]],
            # Around second atom
            ["H",[0.75,0.75,0.25]],
            ["H",[0.75,0.75,0.75]],
        ])

    elif nh == 6:
        # The hydrogens are in an octahedron
        # around the central atom
        atoms.extend([
            # Around the first atom
            ["H",[0.5, 0.0, 0.0]],
            ["H",[0.0, 0.5, 0.0]],
            ["H",[0.0, 0.0, 0.5]],
        ])
    
    par            = parameters()
    alat           = volume_guess ** (1.0/3.0)
    par["atoms"]   = atoms
    par["lattice"] = [[alat,0,0],[0,alat,0],[0,0,alat]]

def get_platonic_hydrogens(nh, scale=0.3333333333, centre=[0.5,0.5,0.5]):

    # Scale of the displacement from the centre
    a = scale

    # The centre of the platonic solid
    x,y,z = centre

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
def generate_structure_platonic(e1, e2, n1, n2, nh, volume_guess):

    # Estimate the lattice length based on the estimated volume
    alat = volume_guess**(1.0/3.0)

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
        par["lattice"] = [[alat,0,0],[0,alat,0],[0,0,alat]]
        par["atoms"]   = atoms
        return par
        
    else:
        return None
    
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
MAX_H_CONTENT = 10

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
