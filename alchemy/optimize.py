from qet.params   import parameters
from qet.logs     import log
from datetime     import datetime
import os, random

# Generate a random valid mutation of the given structure
def randomly_mutate(structure, mutations, check_valid):

    # Nice header to mutation
    line = "{:^64}".format("Mutating "+structure["stochiometry_string"])
    div  = "".join(["-" for c in line])
    log("\n"+div, "alchemy.log")
    log(line, "alchemy.log")
    log(div, "alchemy.log")

    mutation = None
    while (mutation is None) or (not check_valid(mutation)):
        mutation = mutations[random.randrange(len(mutations))](structure)
    return mutation

# Evaluate the objective function on the given
# structure, creating a directory for the
# objective evaluation. Will check existing 
# directories to see if this objective has 
# already been evaluated.
def eval_objective(structure, objective, structure_compare):
    
    name    = structure["stochiometry_string"]
    version = 1
    log("Evaluating objective for {0}".format(name), "alchemy.log")

    # Find previous versions
    for d in os.listdir("."):
        if not os.path.isdir(d): continue
        if not d.startswith(name): continue
        version += 1

        # Parse the structure to see if it's
        # the same (according to structure_compare), 
        # and we need not re-do the optimization
        lattice = []
        atoms   = []
        with open(d+"/objective.log") as f:
            for line in f:
                splt = line.split()

                if splt[0] == "lattice": 
                    lattice.append([float(w) for w in splt[1:]])
                    continue

                if splt[0] == "atom":
                    atoms.append([splt[1], [float(w) for w in splt[2:]]])
                    continue

                if splt[0] == "objective":
                    obj = float(splt[1])
                    continue

        if structure_compare(
            lattice, atoms, 
            structure["lattice"], structure["atoms"]):

            # Structre is the same as previously evaluated
            log("    From equivalent structure in "+d, "alchemy.log")
            return obj
    
    # Create the directory to evaluate this
    # objective in
    obj_dir = "{0}_v{1}".format(name, version)
    log("    From new directory "+obj_dir, "alchemy.log")
    os.system("mkdir "+obj_dir)
    os.chdir(obj_dir)

    with open("objective.log","w") as f:

        # Note the structure, in case we arrive at the
        # same structure later
        for l in structure["lattice"]:
            f.write("lattice {0} {1} {2}\n".format(*l))
        for a in structure["atoms"]:
            f.write("atom {0} {1} {2} {3}\n".format(a[0],*a[1]))

        # Evaluate and record the objective
        obj = objective(structure)
        f.write("objective {0}".format(obj))

    # Move back to previous directory
    os.chdir("..")

    return obj

def optimize(
    start_structure,                # Structure to begin from
    objective,                      # Objective function to minimize
    mutations,                      # Allowed mutations of the structure (see alchemy.mutations)
    check_valid = lambda s : True,  # Check if a given structure is allowable
    max_iter = 100,                 # Maximum optimization steps

    # Function that determines if two structures are simmilar enough
    # to not need recalculating, takes (lattice_1, atoms_1, lattice_2, atoms_2)
    structure_compare = lambda l1, a1, l2, a2 : False
    ):

    # Create the optimization directory
    opt_dir = "alchemical_optimization"
    n = 0
    while os.path.isdir(opt_dir):
        n += 1
        opt_dir = "alchemical_optimization_{0}".format(n)
    os.system("mkdir "+opt_dir)

    # Set the opt_dir as our working directory
    os.chdir(opt_dir)
    log("Starting optimization...", "alchemy.log")
    
    # Initilize the structure and the
    # value of the objective
    structure = start_structure
    last_obj  = eval_objective(structure, objective, structure_compare)

    # Initilize the path
    path = [{
        "structure" : structure,
        "objective" : last_obj,
        "proposal"  : False
    }]

    # Run optimization iterations
    for iteration in range(max_iter):
        
        # Generate a new structure
        new_structure = randomly_mutate(structure, mutations, check_valid)
        new_obj       = eval_objective(new_structure, objective, structure_compare)

        # Add it as a "proposal" point
        # along the path
        path.append({
            "structure" : new_structure,
            "objective" : new_obj,
            "proposal"  : True
        })

        log("Old objective value: {0}".format(last_obj), "alchemy.log")
        log("New objective value: {0}".format(new_obj),  "alchemy.log")

        if new_obj < last_obj:

            # Accept new structure
            structure = new_structure
            last_obj  = new_obj
            log("Mutation accepted", "alchemy.log")
        
        else:
            log("Mutation rejected", "alchemy.log")

        # Record accepted (or reverted)
        # structure/objective
        path.append({
            "structure" : structure,
            "objective" : last_obj,
            "proposal"  : False
        })

    return path

# Plot information about an optimization path
def plot_path(path):
    import matplotlib.pyplot as plt
    
    # Plot objective function evolution
    p1 = plt.subplot(221)
    p2 = plt.subplot(222)

    p1.plot([p["objective"] for p in path if not p["proposal"]])
    p1.set_ylabel("Objective function (best)")
    p1.set_xlabel("Iteration")

    p2.plot([p["objective"] for p in path])
    p2.set_ylabel("Objective function (current)")
    p2.set_xlabel("Iteration")

    plt.show()

def test():

    from qet.alchemy           import mutations
    from qet.test.test_systems import lah10_fm3m

    # Start with Fm3m phase of LaH10
    start = lah10_fm3m

    # Run optimizer
    path = optimize(
        start,                                  # Start structure
        lambda s : len(s["atoms"]),             # Minimize number of atoms
        [                                            
            mutations.remove_random_atom,       # Allowed mutations
            mutations.duplicate_random_atom,
        ],
        lambda s : s["atom_counts"]["H"] > 0,   # Must have at least 1 H

        # We consider structures with the same atoms to be the same
        structure_compare = lambda l1, a1, l2, a2 : True            

        )

    plot_path(path)
