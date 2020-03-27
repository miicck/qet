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

    # Mutate the structure until a valid
    # mutation is found
    mutation = None
    while True:

        mutation = mutations[random.randrange(len(mutations))](structure)

        # Mutations return None if there was no sensible
        # way to apply them
        if mutation is None:
            log("    Mutation returned None, trying again...", "alchemy.log")
            continue

        # Check validty of the mutation
        if not check_valid(mutation):
            log("    Mutation invalid, trying again...", "alchemy.log")
            continue

        # Mutation valid, return it
        return mutation

    raise Exception("Failed to generate a random mutation.")

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
            log("    Objective = {0}".format(obj), "alchemy.log")
            return obj
           
    # Count previous objective evaluations
    # so we know which evaluation this is
    objective_number = 1
    for d in os.listdir("."):
        if not os.path.isdir(d): continue
        if not os.path.isfile(d+"/objective.log"): continue
        objective_number += 1

    # Create the directory to evaluate this
    # objective in
    obj_dir = "{0}_v{1}".format(name, version)
    log("    From new directory "+obj_dir, "alchemy.log")
    os.system("mkdir "+obj_dir)
    os.chdir(obj_dir)

    # Try to calculate the objective
    # set to inf if failed
    try:
        # Evaluate and record the objective
        obj = objective(structure)

    except Exception as e:

        # Flag failed calculation, but don't stop
        log("    Objective evaluation failed with below error\n    {0}".format(e), "alchemy.log")
        obj = float("inf")

    # Write the results of the objective eval
    with open("objective.log","w") as f:

        # Write the objective/evaluation number to file
        f.write("n {0}\n".format(objective_number))
        f.write("objective {0}\n".format(obj))

        # Note the structure, in case we arrive at the
        # same structure later
        for l in structure["lattice"]:
            f.write("lattice {0} {1} {2}\n".format(*l))
        for a in structure["atoms"]:
            f.write("atom {0} {1} {2} {3}\n".format(a[0],*a[1]))

    # Move back to previous directory
    os.chdir("..")

    log("    Objective = {0}".format(obj), "alchemy.log")
    return obj

# Performs an alchemical optimization starting with the given structure,
# and applying the given mutations in order to minimize the given objective. 
# It will not consider structures for which check_valid(structure) is False
# and will stop after max_iter mutations.
#
# For now, it works as follows:
#   1. Randomly mutate structure
#   2. Calculate the objective using the mutated structure
#   3. Accept mutation if new objective is the lowest we've seen
#   4. Repeat
#
def optimize(
    start_structure,                # Structure to begin from
    objective,                      # Objective function to minimize; called as objective(structure)
    mutations,                      # Allowed mutations of the structure (see alchemy.mutations)
    check_valid = lambda s : True,  # Check if a given structure is allowable
    max_iter = 100,                 # Maximum optimization steps
    # Function that determines if two structures of the same stochiometry are simmilar enough
    # to not need recalculating, takes (lattice_1, atoms_1, lattice_2, atoms_2)
    structure_compare = lambda l1, a1, l2, a2 : False
    ):

    # Check arguments are sensible
    if not isinstance(start_structure, parameters):
        raise TypeError("The input structure should be a 'parameters' object!")

    optimize_start_time = datetime.now()

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
    iteration = 1
    while True:
        
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

        iteration += 1
        if iteration > max_iter:
            log("\nMax iter = {0} hit, stopping.".format(max_iter), "alchemy.log")
            break

    # Output time taken
    optimize_end_time = datetime.now()
    fs = "Optimization compete, total time {0}"
    fs = fs.format(optimize_end_time - optimize_start_time)
    log(fs, "alchemy.log")

    return path

# Plot information about an optimization
def plot():
    import matplotlib.pyplot as plt

    # Parse objective.log files
    data = []

    # Loop over objective.log files
    for d in os.listdir("."):
        if not os.path.isdir(d): continue
        if not os.path.isfile(d+"/objective.log"): continue

        with open(d+"/objective.log") as f:
            for line in f:
                splt = line.split()

                # Parse structure number
                if splt[0] == "n":
                    n = int(splt[1])
                    continue
                
                # Parse objective function
                if splt[0] == "objective":
                    obj = float(splt[1])
                    continue

        # Get the Latex structure name from
        # the directory name
        name = ""
        splt = d.split("_")
        for i in range(0, len(splt)-1, 2):
            name += splt[i]
            if int(splt[i+1]) > 1:
                name += "$_{0}$".format("{"+splt[i+1]+"}")

        # Ignore structures that could not be calculated
        if obj == float("inf"):
            continue

        data.append([n, obj, name])

    # Sort data by number and unzip into columns
    data.sort(key=lambda d:d[0])
    data = list(zip(*data))
    
    p1 = plt.subplot(221)

    # Create some extra space at the
    # top of the plot
    max_y   = max(data[1])
    min_y   = min(data[1])
    dm      = (max_y-min_y)/4.0
    max_y  += 2*dm
    p1.set_ylim([min_y, max_y])

    # Plot the objective function
    p1.plot(data[0], data[1], marker="+")

    # Label the structures
    ap = { "arrowstyle" : "->" }
    for i, (x,y) in enumerate(zip(data[0], data[1])):
        p1.annotate(data[2][i], (x,y), xytext=(x,y+dm), arrowprops=ap, rotation=90)

    # Title the axes
    p1.set_ylabel("Objective function")
    p1.set_xlabel("Structure number")

    plt.show()
