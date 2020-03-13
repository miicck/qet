import os
import numpy as np
from constants import BOHR_TO_ANGSTROM
from logs import log

# A parameters object contains the specification
# for the current calculation, including both the
# parameters for the DFT code as well as the 
# specification for the actual calculation to run
class parameters:
    
    # Default constructor
    def __init__(self):

        # Set the default parameters
        # any unspecified parameters will be left 
        # unspecified in input files and therefore
        # result in QE default values being used
        self.par = {}

        self.par["outdir"]         = "./"          # outdir = working dir
        self.par["ibrav"]          = 0             # no bravis-lattice index
        self.par["ecutwfc"]        = 60            # plane-wave cutoff (Ry)
        self.par["ecutrho"]        = 600           # density cutoff (Ry)
        self.par["occupations"]    = "smearing"    # treat as metallic
        self.par["deguass"]        = 0.02          # metal smearing width (Ry)
        self.par["kpoint_spacing"] = 0.02          # kpoint spacing 2pi*0.02A^-1

        # By default, assume cores_per_node is
        # equal to the number of cores where the
        # script is running and nodes = 1
        self.par["nodes"] = 1
        try:
            import multiprocessing
            self.par["cores_per_node"] = multiprocessing.cpu_count()
        except ImportError:
            self.par["cores_per_node"] = 1

        # Default the pseudopotential directory to
        # home/pseudopotentials if $HOME is defined
        # otherwise set to "./"
        self.par["pseudo_dir"] = "./"
        if "HOME" in os.environ:
            self.par["pseudo_dir"] = os.environ["HOME"]+"/pseudopotentials"
            

    # Get parameter values with []
    def __getitem__(self, key):

        # Raise an error if the key isn't found
        if not key in self.par:

            # Try to generate the parameter from
            # the parameters we have
            
            # Count the atoms
            if key == "nat" : 
                self["nat"] = len(self["atoms"])

            # Count the atom types
            elif key == "ntyp": 
                self["ntyp"] = len(self["species"])

            # Get a list of the species
            elif key == "species":
                spec = []
                for a in self["atoms"]: 
                    if a[0].capitalize() in spec: continue
                    spec.append(a[0].capitalize())
                self["species"] = spec

            # Generate the kpoint grid
            elif key == "kpoint_grid":
                rlat = np.linalg.inv(self["lattice"]).T
                b2k  = lambda b : int(np.linalg.norm(b)/self["kpoint_spacing"])
                self["kpoint_grid"] = [b2k(b) for b in rlat]

            else:
                # Could not generate, error out
                exept = "Key \"{0}\" not found in parameters object."
                raise ValueError(exept.format(key))

        return self.par[key]


    # Set parameter values with []
    def __setitem__(self, key, value):
        
        # Set the given parameter
        self.par[key] = value


    # Convert a parameter to an input line in a QE file
    def to_input_line(self, key):

        # Return nothing if the key is absent
        # so the QE default value will be used
        try:
            val = self[key]
        except:
            return ""

        # Try to write an integer parameter
        try:
            fval = float(val)
            ival = int(val)
            # If the floating point value differs
            # this was actually a float that has been
            # truncated
            if fval != float(ival): raise ValueError
            return "    {0} = {1},\n".format(key, ival)
        except: pass

        # Try to write a float parameter
        try:
            val = float(val)
            return "    {0} = {1},\n".format(key, val)
        except: pass

        # Write bool parameters
        if val.lower().strip() == "true":
            return "    {0} = .true.,\n".format(key)

        if val.lower().strip() == "false":
            return "    {0} = .false.,\n".format(key)

        # Write string parameters
        return "    {0} = '{1}',\n".format(key, val)


    # Convert parameters object to string
    def __str__(self):

        # Pad format string nicely
        maxl = max([len(p) for p in self.par])
        fs   = "\n    {0:" + str(maxl) + "} : {1}"
        
        s  = "Parameters:"
        s += "\n(parameters not in this list will be left as the QE defaults.)"
        for p in self.par:
            
            # Custom formatting for atoms
            if p == "atoms":
                atoms = self.par["atoms"]
                s += fs.format("atoms", len(atoms))
                afs = "\n        {0:3} {1:10.7f} {2:10.7f} {3:10.7f}"
                for a in atoms:
                    s += afs.format(a[0].capitalize(), *a[1])        
                continue

            # Custom formatting for lattice
            if p == "lattice":
                lattice = self.par["lattice"]
                s += fs.format("lattice", "")
                for l in lattice:
                    s += "\n        {0:10.7f} {1:10.7f} {2:10.7f}".format(*l)
                continue

            # Output key : value
            s += fs.format(p, self.par[p])
        return s

    # Parse the crystal lattice from input file lines
    def parse_lattice(self, lines):

        for i, line in enumerate(lines):

            # Find the start of the lattice block
            if not line.startswith("lattice"):  
                continue

            # Parse the lattice units
            if len(line.split()) < 2:
                raise Exception("Lattice input must specify units!")
            units = line.split()[1]
            if   units == "angstrom": factor = 1.0
            elif units == "bohr"    : factor = BOHR_TO_ANGSTROM

            # Parse the lattice
            lat = []
            for j in range(i+1, i+4):
                lat.append([factor*float(w) for w in lines[j].split()])
            self["lattice"] = np.array(lat)
            break

    # Parse atom coordinates from input file lines
    def parse_atoms(self, lines):

        for i, line in enumerate(lines):

            # Find start of atom block 
            if not line.startswith("atoms"):
                continue
            
            # Parse the units the atoms are in
            if len(line.split()) < 2:
                raise Exception("Atom positions input must specify units!")

            # Get the inverse lattice in the same
            # units, to convert to fractional coordinates
            units = line.split()[1]
            if units == "crystal":
                linv = np.identity(3)
            elif units == "fractional":
                linv = np.identity(3)
            elif units == "angstrom":
                linv = np.linalg.inv(self["lattice"].T)
            elif units == "bohr":
                linv = la.inv(ANGSTROM_TO_BOHR*self["lattice"].T)
            else:
                raise ValueError("Unknown atom coordinate units: "+units)

            # Read the atoms one by one, until end of file or blank line
            atoms = []
            for j in range(i+1, len(lines)):
                line = lines[j]

                # Stop if we reach an empty line
                if len(line) == 0: break
                
                name = line.split()[0]
                pos  = [float(x) for x in line.split()[1:]]
                pos  = np.dot(linv, pos) # Convert to fractional coords
                atoms.append([name,pos])

            self["atoms"] = atoms

def read_parameters(filename):

    # Create the default parameters object
    params = parameters()

    # Read the file into lines
    with open(filename) as f:
        lines = f.read().split("\n")

    # Parse the file, line by line
    log("Reading parameters from {0}:".format(filename))

    # Loop over lines, use a while loop so
    # i can be incremented by the parsing logic
    for i in range(0, len(lines)):
    
        # Clean the line, remove comments/whitespace
        line = lines[i]
        line = line.lower()
        for cc in ["#", "//", "!"]:
            if line.find(cc) >= 0:
                line = line[0:line.find(cc)]
        line = line.strip()
        lines[i] = line

    params.parse_lattice(lines)
    params.parse_atoms(lines)

    return params
