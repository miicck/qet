from   qet.params import parameters
from   qet.logs   import log, logging_enabled
from   time       import sleep
from   filelock   import FileLock, Timeout
import numpy      as     np
import os, random, math

# A vertex of an alch_network, designed to deal with multiple processes
# trying to access / modify the vertex
class alch_vertex:

    # Load a vertex from a directory
    def __init__(self, vertex_dir):

        if not os.path.isdir(vertex_dir):
            fs = "Tried to load a vertex from the non-existant directory: {0}"
            raise Exception(fs.format(vertex_dir))

        # Record my vertex directory
        vertex_dir = os.path.abspath(vertex_dir)
        self.dir   = vertex_dir
        self.name  = self.dir.split("/")[-1]
        self.locks = {}

    # Get a lock by name (creates it if it doesn't exist)
    def lock(self, name):
        if name in self.locks: return self.locks[name]
        lock = FileLock(self.dir+"/"+name)
        self.locks[name] = lock
        return lock

    # Convert to vertex
    def __str__(self):
        return "alch_vertex at {0}".format(self.dir)

    # Get the evaluated objectives in dictionary form
    def get_evaluated_objectives(self):

        with self.lock("obj_file_lock"):
            ret = {}
            if os.path.isfile(self.dir+"/objectives"):
                with open(self.dir+"/objectives") as f:
                    for line in f:
                        splt = line.split(":")
                        ret[splt[0]] = float(splt[1])
            return ret

    # Evaluate the given objective function
    # on the parameter set for this vertex
    def objective(self, objective_name, objective_function):

        # Return stored objective value
        stored = self.get_evaluated_objectives()
        if objective_name in stored: 
            return stored[objective_name]

        # Attempt to evaluate the objective. 
        #
        # If evaluataion fails (i.e the underlying objective_function
        # could not be evaluated), the objective is set to inf.
        #
        # If evaluation is already underway on another process (i.e we cannot 
        # acquire objective_lock within some timeout), return inf this time,
        # but do not store inf for future retrieval.
        #
        try:
            with self.lock("objective_lock").acquire(timeout=1.0):

                # Set the working directory to the vertex directory
                old_wd = os.getcwd()
                os.chdir(self.dir)

                try:
                    # Attempt to evaluate the objective
                    log("Evaluating objective in {0}".format(self.name), "alchemy.log")
                    obj = objective_function(self.params)

                    # If the objective has length 2, we assume it has the form
                    # [objective value, new structure]. If it has length 1 assume 
                    # it is of the form [objective value], but issue a warning,
                    # as the function should probably have not returned a list if
                    # there was only one thing to return. All other lengths are errors.
                    if hasattr(obj, "__len__"):
                        if   len(obj) >  2: raise Exception("Objective returned too many results!")
                        elif len(obj) == 0: raise Exception("Objective returned too few results!")
                        elif len(obj) == 1:
                            fs = "objective assumed to be of the form [objective value] in {0}"
                            log(fs.format(self.name), "alchemy.warnings")
                            obj = obj[0]
                        else:
                            fs = "Objective returned [objective value, new structure] in {0}"
                            log(fs.format(self.name), "alchemy.log")
                            self.params = obj[1]
                            obj = obj[0]

                    fs = "evaluated objective in {2}\n    {0} = {1}"
                    log(fs.format(objective_name, obj, self.name), "alchemy.log")

                except Exception as e:
                    # Failed to evaluate the objective, set it to inf
                    fs = "failed to evaluate the objective\n    {0} in {1}"
                    log(fs.format(objective_name, self.name), "alchemy.log")
                    log("error: {0}".format(e), "alchemy.log")
                    obj = float("inf")

                # Write the values of all the objectives to file
                with self.lock("obj_file_lock"):
                    stored = self.get_evaluated_objectives()
                    stored[objective_name] = obj
                    with open("objectives", "w") as f:
                        for o in stored:
                            f.write("{0}:{1}\n".format(o, stored[o]))

                # Restore working directory
                os.chdir(old_wd)
                return obj

        except Timeout:
            fs = "Evaluation already underway in {0}".format(self.name)
            log(fs, "alchemy.log")
            return float("inf")
        
    # Return my parameters, as recorded in the /parameters file
    @property
    def params(self):
        with self.lock("parameters_lock"):
            filename = self.dir + "/parameters"
            if not os.path.isfile(filename):
                raise Exception("Parameters for vertex {0} not found!".format(self.dir))
            return parameters(filename=filename)

    # Set my parameters (by saving them in the /parameters file)
    @params.setter
    def params(self, params):
        with self.lock("parameters_lock"):
            filename = self.dir + "/parameters"
            params.save(filename)

    # The verticies that were mutated to
    # obtain this vertex
    @property
    def parents(self):
        with self.lock("parents_lock"):
            par = []
            if not os.path.isfile(self.dir + "/parents"): return par
            with open(self.dir + "/parents") as f:
                for line in f:
                    line = line.strip()
                    if len(line) == 0: continue
                    if line in par: continue
                    par.append(line)
            return par

    # Add a parent vertex
    def add_parent(self, parent):
        with self.lock("parents_lock"):
            pts = self.parents
            pdir = parent.dir
            pdir = pdir.split("/")[-1]
            if pdir in pts: return
            pts.append(pdir)
            with open(self.dir + "/parents", "w") as f:
                for p in pts:
                    f.write(p+"\n")

    # Get my name, in latex format
    @property
    def latex_name(self):
        splt    = self.name.split("_")
        version = splt[-1]
        splt    = splt[0:-1]
        ltx     = ""
        for i in range(0, len(splt), 2):
            if int(splt[i+1]) > 1:
                ltx += "{0}$_{1}$".format(splt[i], "{"+splt[i+1]+"}")
            else:
                ltx += splt[i]
        return ltx + " ({0})".format(version)
            

# An alch_network is a graph of alch_vertexs, where each vertex corresponds
# to a particular parameter set.
class alch_network:

    def __init__(self, 

        # The directory representing this network, verticies will be
        # represented by subdirectories within base_dir.
        # If None, a new directory of the form network_i will be created
        # If base_dir exists, we will load the network in base_dir
        # Otherwise a new network will be created with name base_dir
        base_dir = None,                            

        # Function used to compare two parameters sets of the same stochiometry
        # returns true if they should be considered as the same vertex.
        parameters_comp = lambda s1, s2 : True,  # Default: same stoichiometry => same vertex

        ):

        # Generate new network name
        if base_dir is None:
            i = 1
            while True:
                if not os.path.isdir("network_{0}".format(i)):
                    base_dir = "network_{0}".format(i)
                    break
                i += 1

        # Ensure the network directory exists
        created = False
        if not os.path.isdir(base_dir):
            os.system("mkdir "+base_dir)
            created = True

        # Use the absolute path from here on
        print(base_dir)
        base_dir  = os.path.abspath(base_dir)
        self.dir  = base_dir
        self.name = base_dir.split("/")[-1]

        # Record the parameters comparison function
        self.parameters_comp = parameters_comp

        if created: log("Created network "+self.dir, "alchemy.log")
        else: 
            fs = "Loaded network {0} with {1} verticies"
            fs = fs.format(self.dir, len(self.verticies))
            log(fs, "alchemy.log")


    # My verticies
    @property
    def verticies(self):

        # Load the verticies from disk every time,
        # in case there are multiple processes working
        # on the network, creating verticies.
        verts = []
        
        # Loop over subdirectories
        for d in os.listdir(self.dir):
            d = self.dir + "/" + d
            if not os.path.isdir(d): continue

            # Load this vertex
            verts.append(alch_vertex(d))

        return verts

    # Return a random vertex
    def random_vertex(self):
        vs = self.verticies
        return vs[random.randrange(len(vs))]

    # Attempt to create a new vertex with the given parameters
    # if the vertex already exists, it will return that
    # if it fails it will return None
    # otherwise returns the newly created vertex
    def create_vertex(self, params):

        if params is None:
            raise Exception("Tried to create a vertex with params = None!")

        # Label directories with stoichiometry
        base_vert_dir = self.dir + "/" + params["stoichiometry_string"] 

        # Loop over verticies with the same stochiometry
        version = 1
        while True:

            vertex_dir = base_vert_dir + "_v{0}".format(version)

            # We've found a free vertex name
            if not os.path.isdir(vertex_dir): break

            # Vertex name exists, increment version
            version += 1
            
            # If these verticies are the same (as far as
            # the network is concerned), return the one that
            # already exists
            vertex = alch_vertex(vertex_dir)
            if self.parameters_comp(params, vertex.params):
                log("Vertex {0} already exists".format(vertex.name), "alchemy.log")
                return vertex

        # Initialize the vertex directory
        os.system("mkdir "+vertex_dir)
                
        new_vertex        = alch_vertex(vertex_dir)
        new_vertex.params = params
        log("Created vertex {0}".format(new_vertex.name), "alchemy.log")
        return new_vertex

    # Expand upon a given vertex by applying
    # the given parameter mutation. Returns 
    # the new vertex if successful, otherwise None
    def expand_vertex(self, vertex, param_mutation, is_valid):

        if not os.path.isdir(vertex.dir):
            raise Exception("Tried to expand a vertex that wasn't in the network!")

        # Nice formatting
        expand_text = "{:^64}".format("Expanding vertex "+vertex.name)
        underline   = "".join("-" for c in expand_text)
        log("\n"+expand_text, "alchemy.log")
        log(underline, "alchemy.log")

        # Apply mutation, return None on failure
        mutation = param_mutation(vertex.params)
        if mutation is None: 
            log("Mutation {0} returned None".format(param_mutation.__name__), "alchemy.log")
            log(underline, "alchemy.log")
            return None

        # Check structure produced is valid
        if not is_valid(mutation):
            log("Mutation {0} produced invalid structure".format(param_mutation.__name__), "alchemy.log")
            log(underline, "alchemy.log")
            return None

        # Adjust the volume of the new structure in an
        # attempt to accomodate the mutation
        cov_volume_factor = mutation["covalent_volume"]/vertex.params["covalent_volume"]
        lat_volume_factor = np.linalg.det(mutation["lattice"])/np.linalg.det(vertex.params["lattice"])
        volume_boost = cov_volume_factor / lat_volume_factor
        mutation["lattice"] *= volume_boost ** (1.0/3.0)

        # Create new vertex
        new_vertex = self.create_vertex(mutation)
        new_vertex.add_parent(vertex)
        log(underline, "alchemy.log")
        return new_vertex

    # Choose a random vertex according to weights generated
    # by the given weight function
    def random_weighted_vertex(self, weight=lambda v : 1.0):  

        verts = self.verticies
        if len(verts) == 0: return None
        weights = [weight(v) for v in verts]
        return random.choices(verts, weights=weights)[0]

    # By default choose vertex by exp(-objective) weighting
    def default_vertex_chooser(self, objective_name, objective_function):
        
        wf = lambda v : np.exp(-v.objective(objective_name, objective_function))
        return self.random_weighted_vertex(weight = wf)

    # Expand the network with an aim to minimize the given objective
    # will choose a vertex with low objective value and expand it using
    # one of the given mutations
    def expand_to_minimize(self,
        objective_name,     # Name of objective to minimize
        objective_function, # Objective to minimize
        mutations,          # List of allowed mutations to vertex parameters
        # Function used to determine if a mutated structure is valid
        is_valid = lambda s : True,
        # The function to choose a vertex to expand from a network
        vertex_chooser = lambda n, on, of : n.default_vertex_chooser(on, of)
        ):

        if len(mutations) == 0:
            raise Exception("No mutations specified in expand_to_minimize")

        # For now, choose a random mutation
        mut = mutations[random.randrange(len(mutations))]

        log("Choosing vertex to expand (minimizing {0})...".format(objective_name), "alchemy.log")
        vert = vertex_chooser(self, objective_name, objective_function)
        log("Vertex chosen: "+vert.name, "alchemy.log")
        self.expand_vertex(vert, mut, is_valid)

    # Plot this network
    def plot(self, objective=None):
        import matplotlib.pyplot as plt
        import networkx          as nx
            
        # Get the verticies
        verts = self.verticies
        names = [v.latex_name for v in verts]
        objs  = [v.get_evaluated_objectives() for v in verts]

        # If objective=None, color according to first objective found
        if objective is None:
            for o in objs:
                if not objective is None:
                    break
                for k in o:
                    objective = k
                    break

        # Obtain colors from objective function values
        # green = lowest value of objective
        # blue  = highest value of objective
        # red   = infinite/not calculated objective
        # white = underway
        vals = []
        for o in objs:
            if objective in o: vals.append(o[objective])
            else: vals.append(float("inf"))
        max_v = max(v for v in vals if math.isfinite(v))
        min_v = min(v for v in vals if math.isfinite(v))
        vals = [(v-min_v)/(max_v-min_v) for v in vals]
        colors = [[0, 1.0-v, v, 0.4] if math.isfinite(v) else [1.0,0,0,0.4] for v in vals]
        for i, o in enumerate(objs):
            if not objective in o:
                colors[i] = [1.0,1.0,1.0,0.4]

        # Construct the graph
        g = nx.DiGraph()
        for v in verts:
            for p in v.parents:
                pname = alch_vertex(self.dir+"/"+p).latex_name
                g.add_edge(pname, v.latex_name)
        g.add_nodes_from(names)

        # Get the node positions
        pos = nx.drawing.nx_agraph.graphviz_layout(g)

        # Setup the plot area
        max_x = max(pos[n][0] for n in names)
        min_x = min(pos[n][0] for n in names)
        max_y = max(pos[n][1] for n in names)
        min_y = min(pos[n][1] for n in names)
        plt.xlim(min_x, max_x)
        plt.ylim(min_y, max_y)
        plt.axis("off")

        # Draw nodes
        for n,c in zip(names,colors):
            args = dict(
                verticalalignment="center", 
                horizontalalignment="center",
            )
            plt.annotate(n, pos[n], bbox=dict(facecolor="white"), **args)
            plt.annotate(n, pos[n], bbox=dict(facecolor=c), **args)

        # Get the edges
        edges = [[pos[e[0]], pos[e[1]]] for e in g.edges()]

        # Draw edges
        for (xy1, xy2) in edges:
            dxy = [xy2[0]-xy1[0],xy2[1]-xy1[1]]
            plt.arrow(xy1[0], xy1[1], dxy[0], dxy[1], head_width=0)

        # Draw arrows
        for (xy1, xy2) in edges:
            dxy    = [xy2[0]-xy1[0],xy2[1]-xy1[1]]
            centre = [(xy1[0]+xy2[0])/2.0, (xy1[1]+xy2[1])/2.0]
            plt.arrow(centre[0]-dxy[0]/128, centre[1]-dxy[1]/128, dxy[0]/128, dxy[1]/128,
                      head_width=10, fc="black", ec=None)


        plt.show()

def plot_alch_network(directory=None):
    if directory is None: directory = os.getcwd()
    logging_enabled(False)
    alch_network(directory).plot()
    logging_enabled(True)

##########################
# TESTS FOR ALCH_NETWORK #
##########################

def min_atoms_test():
    from qet.test.test_systems import lah10_fm3m
    from qet.alchemy           import mutations

    nw = alch_network("test_network")
    v  = nw.create_vertex(lah10_fm3m)

    obj_name = "atom count"
    obj_func = lambda s : len(s["atoms"])
    muts     = [mutations.substitute_random_species,
                mutations.remove_random_atom,
                mutations.duplicate_random_atom]

    for n in range(100):
        nw.expand_to_minimize(obj_name, obj_func, muts)

def get_minus_dos_ef(structure):
    from qet.calculations import relax, dos

    # Relax the structure 
    rlx = relax(structure)
    res = rlx.run()

    # Calculate DOS of the relaxed structure
    structure["atoms"]   = res["relaxed atoms"]
    structure["lattice"] = res["relaxed lattice"]
    dos_calc = dos(structure)
    res      = dos_calc.run()

    # Maximize DOS => minimize -DOS
    # we also return the relaxed parameters, which will
    # replace the parameters for this vertex
    return [-res["DOS (E_F)"], structure]

def metal_test():
    from qet.test.test_systems import bcc_lithium
    from qet.alchemy           import mutations

    #bcc_lithium["cores_per_node"] = 1
    bcc_lithium["pseudo_dir"]     = "/home/mjh261/rds/rds-t2-cs084/pseudopotentials/gbrv"
    nw = alch_network("metals_network")
    v  = nw.create_vertex(bcc_lithium)

    muts = [mutations.substitute_random_species,
            mutations.remove_random_atom,
            mutations.duplicate_random_atom]

    for n in range(100):
        nw.expand_to_minimize("-DOS (E_F)", get_minus_dos_ef, muts)