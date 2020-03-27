from qet.alchemy.optimize  import optimize, plot_path
from qet.alchemy           import mutations
from qet.params            import parameters
from qet.calculations      import relax, dos

start_structure = parameters()

# Start with a poor metal, BCC Lithium
start_structure["atoms"] = [
    ["Li", [0.0, 0.0, 0.0]]
]

start_structure["lattice"] = [
    [-1.8,  1.8,  1.8],
    [ 1.8, -1.8,  1.8],
    [ 1.8,  1.8, -1.8],
]

# Use cheap parameters for this test
start_structure["ecutwfc"]        = 20
start_structure["ecutrho"]        = 200
start_structure["kpoint_spacing"] = 0.05
start_structure["pseudo_dir"]     = "/home/mjh261/rds/rds-t2-cs084/pseudopotentials/gbrv"

# The mutations we're allowed to make to
# the structure during the optimization
allowed_mutations = [
    mutations.substitute_random_species,
    #mutations.remove_random_atom,    # There would be no point doing these mutations
    #mutations.duplicate_random_atom, # as we restrict the serach to 1 atom per cell
    #mutations.shuffle_atoms,         # anyway.
]

# Only consider structures with a single atom per cell
# (not that we can get more with the above mutations anyway)
def is_valid(structure):
    if len(structure["atoms"]) > 1: return False
    return True

# This is the thing we're trying to minimize.
# In our case, relax the structure and calculate
# the density of states at the fermi level.
# - DOS(E_F) is then the objective because,
# for a good metal, we want to maximize DOS(E_F).
def objective(structure):

    # Relax the structure
    rlx        = relax(structure)
    rlx_result = rlx.run()

    # Get the relaxed geometry
    rlx_structure = structure.copy()
    rlx_structure["atoms"]   = rlx_result["relaxed atoms"]
    rlx_structure["lattice"] = rlx_result["relaxed lattice"]

    # Calcualte the density of states
    # of the relaxed geometry
    ds         = dos(rlx_structure)
    ds_result  = ds.run()

    return -ds_result["DOS (E_F)"]

# We consider materials with the same stochiometry equivelant
def compare_structures(lattice_1, atoms_1, lattice_2, atoms_2):
    return True

optimize(
    start_structure,
    objective,
    allowed_mutations,
    check_valid = is_valid,
    structure_compare = compare_structures
)
