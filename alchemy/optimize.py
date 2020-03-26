from qet.params   import parameters
from qet.logs     import log
import random

# Generate a valid mutation of the given structure
def randomly_mutate(structure, mutations, check_valid):

    mutation = None
    while (mutation is None) or (not check_valid(mutation)):
        mutation = mutations[random.randrange(len(mutations))](structure)
    return mutation
        
def optimize(
    start_structure,                # Structure to begin from
    objective,                      # Objective function to minimize
    mutations,                      # Allowed mutations of the structure (see alchemy.mutations)
    check_valid = lambda s : True,  # Check if a given structure is allowable
    max_iter = 100                  # Maximum optimization steps
    ):
    
    # Initilize the structure and the
    # value of the optimization
    structure = start_structure
    last_obj  = objective(structure)

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
        new_obj       = objective(new_structure)

        # Add it as a "proposal" point
        # along the path
        path.append({
            "structure" : new_structure,
            "objective" : new_obj,
            "proposal"  : True
        })

        if new_obj < last_obj:

            # Accept new structure
            structure = new_structure
            last_obj  = new_obj

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

    from qet.alchemy import mutations

    # Start with LiMgH6
    start = parameters()

    start["lattice"] = [
        [4.0, 0, 0],
        [0, 4.0, 0],
        [0, 0, 4.0],
    ]

    start["atoms"] = [
        ["Li", [0,    0,    0   ]],
        ["Mg", [0.5,  0.5,  0.5 ]],
        ["H" , [0.75, 0.5,  0.5 ]],
        ["H" , [0.25, 0.5,  0.5 ]],
        ["H" , [0.5,  0.75, 0.5 ]],
        ["H" , [0.5,  0.25, 0.5 ]],
        ["H" , [0.5,  0.5,  0.75]],
        ["H" , [0.5,  0.5,  0.25]],
    ]

    # Run optimizer
    path = optimize(
        start,                                  # Start structure
        lambda s : len(s["atoms"]),             # Minimize number of atoms
        [                                            
            mutations.remove_random_atom,       # Allowed mutations
            mutations.duplicate_random_atom,
        ],
        lambda s : s["atom_counts"]["H"] > 0    # Must have at least 1 H
        )

    plot_path(path)
