from qet.params import parameters
from qet.logs   import log
import os, random

# A vertex of an alch_network
class alch_vertex:

    # Load a vertex from a directory
    def __init__(self, vertex_dir):

        if not os.path.isdir(vertex_dir):
            raise Exception("Tried to load a vertex from a non-existant directory!")

        # Record my vertex directory
        vertex_dir = os.path.abspath(vertex_dir)
        self.dir   = vertex_dir

    # Convert to vertex
    def __str__(self):
        return "alch_vertex at {0}".format(self.dir)

    # Get the evaluated objectives in dictionary form
    def get_evaluated_objectives(self):
        
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
            fs = "using previously evaluated objective {0} = {1} in {2}"
            log(fs.format(objective_name, stored[objective_name], self.dir), "alchemy.log")
            return stored[objective_name]

        # Set the working directory to the vertex directory
        old_wd = os.getcwd()
        os.chdir(self.dir)

        try:
            # Attempt to evaluate the objective
            obj = objective_function(self.params)
            fs = "evaluated objective {0} = {1} in {2}"
            log(fs.format(objective_name, obj, self.dir), "alchemy.log")
        except:
            # Failed to evaluate the objective, set it to inf
            fs = "failed to evaluate the objective {0} in {1}"
            log(fs.format(objective_name, self.dir), "alchemy.log")
            obj = float("inf")

        # Write the values of all the objectives to file
        stored[objective_name] = obj
        with open("objectives", "w") as f:
            for o in stored:
                f.write("{0}:{1}\n".format(o, stored[o]))

        # Restore working directory
        os.chdir(old_wd)
        
    # Return my parameters, as recorded in the
    # /parameters file
    @property
    def params(self):
        filename = self.dir + "/parameters"
        if not os.path.isfile(filename):
            raise Exception("Parameters for vertex {0} not found!".format(self.dir))
        return parameters(filename=filename)


# An alch_network is a graph, where each vertex corresponds
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
        if not os.path.isdir(base_dir):
            os.system("mkdir "+base_dir)
            log("Created network in "+base_dir, "alchemy.log")
        else:
            log("Loaded network from "+base_dir,"alchemy.log")

        # Use the absolute path from here on
        base_dir  = os.path.abspath(base_dir)
        self.dir  = base_dir

        # Record the parameters comparison function
        self.parameters_comp = parameters_comp

    # My verticies
    @property
    def verticies(self):

        # Load the verticies from disk every time,
        # in case there are multiple processes working
        # on the network, creating verticies.
        verts = []
        
        # Loop over subdirectories
        for d in os.listdir(self.dir):
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
                return vertex

        # Initialize the vertex directory
        os.system("mkdir "+vertex_dir)
        params.save(vertex_dir+"/parameters")
                
        return alch_vertex(vertex_dir)

    # Expand upon a given vertex by applying
    # the given parameter mutation. Returns 
    # the new vertex if successful, otherwise None
    def expand_vertex(self, vertex, param_mutation):

        if not os.path.isdir(vertex.dir):
            raise Exception("Tried to expand a vertex that wasn't in the network!")

        # Apply mutation, return None on failure
        mutation = param_mutation(vertex.params)
        if mutation is None: return None

        # Create new vertex, return None on failure
        new_vertex = self.create_vertex(mutation)
        if new_vertex is None: return None

        # Add the new vertex to the network and return it
        return new_vertex

def test():
    from qet.test.test_systems import lah10_fm3m
    from qet.alchemy           import mutations
    nw = alch_network("test_network")
    v = nw.create_vertex(lah10_fm3m)

    for n in range(100):
        v = nw.expand_vertex(v, mutations.substitute_random_species)
        v.objective("atom count", lambda s : len(s["atoms"]))
