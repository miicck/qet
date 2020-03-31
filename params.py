import os, numbers, copy
import numpy          as     np
from   qet.constants  import BOHR_TO_ANGSTROM
from   qet.type_tools import str_to_type 
from   qet.elements   import elements
from   collections    import defaultdict

# A parameters object contains the specification
# for the current calculation, including both the
# parameters for the DFT code as well as the 
# specification for the actual calculation to run
class parameters:
    
    # Default constructor
    def __init__(self, filename=None):

        # Set the default parameters
        # any unspecified parameters will be left 
        # unspecified in input files and therefore
        # result in QE default values being used
        self.par = {}

        self["outdir"]         = "./"          # outdir = working dir
        self["ibrav"]          = 0             # no bravis-lattice index
        self["ecutwfc"]        = 60            # plane-wave cutoff (Ry)
        self["ecutrho"]        = 600           # density cutoff (Ry)
        self["occupations"]    = "smearing"    # treat as metallic
        self["degauss"]        = 0.02          # metal smearing width (Ry)
        self["qpoint_spacing"] = 0.1           # qpoint spacing (2pi A^-1)
        self["kpts_per_qpt"]   = 6             # ratio of kpt to qpt grid

        # By default, assume cores_per_node is
        # equal to the number of cores where the
        # script is running and nodes = 1
        self["nodes"] = 1
        try:
            import multiprocessing
            self["cores_per_node"] = multiprocessing.cpu_count()
        except ImportError:
            self["cores_per_node"] = 1

        # Default the pseudopotential directory to
        # home/pseudopotentials if $HOME is defined
        # otherwise set to "./"
        self["pseudo_dir"] = "./"
        if "HOME" in os.environ:
            self["pseudo_dir"] = os.environ["HOME"]+"/pseudopotentials"

        if not filename is None:
            self.load(filename)

    # Return a deep copy of myself
    def copy(self):
        return copy.deepcopy(self)

    # Try to generate the parameter from
    # the parameters we have
    def gen_param(self, key):
        
        # Count the atoms
        if key == "nat" : 
            return len(self["atoms"])

        # Count the atom types
        if key == "ntyp": 
            return len(self["species"])

        # Get a list of the species
        # with masses and pseudo names
        if key == "species":
            spec = []
            for a in self["atoms"]: 
                if a[0] in spec: continue
                spec.append(a[0])

            for i, s in enumerate(spec):
                spec[i] = [s, elements[s]["mass number"], s+".UPF"]

            return spec

        # Get a dictionary of the form atom name : count
        if key == "atom_counts":
            atom_counts = defaultdict(lambda: 0)
            for a in self["atoms"]:
                if a[0] in atom_counts: atom_counts[a[0]] += 1
                else: atom_counts[a[0]] = 1
            return atom_counts

        # Reurn the stochiometry 
        # of the given cell as a string
        if key == "stoichiometry_string":
            atom_counts = self["atom_counts"]
            ss = ""
            for a in atom_counts:
                ss += a + "_{0}_".format(atom_counts[a])
            ss = ss[0:-1]
            return ss

        # Get an estimate for the volume of the
        # cell by approximating each atom as 
        # a covalent-radius sphere
        if key == "covalent_volume":
            vol = 0.0
            for a in self["atoms"]:
                vol += elements[a[0]]["covalent radius"]**3.0
            return np.pi * vol * 4.0/3.0

        # Generate the qpoint grid
        if key == "qpoint_grid":
            
            # Generate qpoint grid from spacing
            rlat = np.linalg.inv(self["lattice"]).T
            qps  = float(self["qpoint_spacing"])
            b2q  = lambda b : int(np.linalg.norm(b)/qps)
            return [b2q(b) for b in rlat]

        # Generate the kpoint grid
        if key == "kpoint_grid":

            if "kpoint_spacing" in self.par:

                # Generate kpoint grid from spacing
                rlat = np.linalg.inv(self["lattice"]).T
                kps  = float(self["kpoint_spacing"])
                b2k  = lambda b : int(np.linalg.norm(b)/kps)
                return [b2k(b) for b in rlat]

            elif "kpts_per_qpt" in self.par:

                # Generate kpoint grid from qpoint grid
                kpq = self["kpts_per_qpt"] 
                qpg = self["qpoint_grid"]
                return [kpq * q for q in qpg]

            else:

                msg = "Could not generate k-point grid from parameter set."
                raise RuntimeError(msg)

        # Could not generate, error out
        exept = "Key \"{0}\" cold not be generated in parameters object."
        raise ValueError(exept.format(key))
            

    # Get parameter values with []
    def __getitem__(self, key):

        # Attempt to generate the
        # parameter if key if not found
        if not key in self.par:
            return self.gen_param(key)

        return self.par[key]

    # Validate and standardise a given parameter
    # based on it's key and value
    def validate_and_standardise_param(self, key, value):

        if key == "atoms":
            # Sort atoms in decreasing atomic number
            value.sort(key = lambda a : -elements[a[0]]["atomic number"])

        elif key == "lattice":
            # Lattice should be stored as a np.array
            value = np.array(value)

        return value

    # Set parameter values with []
    def __setitem__(self, key, value):

        # Try to convert a string value to the correct type
        if isinstance(value, str):
            value = str_to_type(value)

        self.par[key] = self.validate_and_standardise_param(key, value)


    # Convert a parameter to an input line in a QE file
    def to_input_line(self, key):

        # Return nothing if the key is absent
        # so the QE default value will be used
        try:
            val = self[key]
        except:
            return ""

        # Bool type
        if isinstance(val, bool):
            if val: return "    {0} = .true.,\n".format(key)
            else  : return "    {0} = .false.,\n".format(key)

        # Numeric types
        if isinstance(val, numbers.Number):
            return "    {0} = {1},\n".format(key, val)

        # String parameters
        return "    {0} = '{1}',\n".format(key, val)


    # Convert parameters object to string
    def __str__(self):

        # Pad format string nicely
        maxl = max([len(p) for p in self.par])
        fs   = "\n    {0:" + str(maxl) + "} : {1} ({2})"
        
        s  = "Parameters:"
        s += "\n(parameters not in this list will be left as the QE defaults.)"
        for p in self.par:
            
            # Custom formatting for atoms
            if p == "atoms":
                atoms = self.par["atoms"]
                s += fs.format("atoms", len(atoms), "atoms")
                afs = "\n        {0:3} {1:10.7f} {2:10.7f} {3:10.7f}"
                for a in atoms:
                    s += afs.format(a[0], *a[1])        
                continue

            # Custom formatting for lattice
            if p == "lattice":
                lattice = self.par["lattice"]
                s += fs.format("lattice", "", "matrix")
                for l in lattice:
                    s += "\n        {0:10.7f} {1:10.7f} {2:10.7f}".format(*l)
                continue

            # Output key : value
            s += fs.format(p, self.par[p], type(self.par[p]))
        return s

    # Parse the crystal lattice from input file lines
    def parse_lattice(self, lines):

        i_dealt_with = []
        for i, line in enumerate(lines):

            # Find the start of the lattice block
            if not line.startswith("lattice"):  
                continue

            # Log the lines that have been dealt with during this parsing
            i_dealt_with.append(i)

            # Parse the lattice units
            if len(line.split()) < 2:
                raise Exception("Lattice input must specify units!")
            units = line.split()[1]
            if   units == "angstrom": factor = 1.0
            elif units == "bohr"    : factor = BOHR_TO_ANGSTROM

            # Parse the lattice
            lat = []
            for j in range(i+1, i+4):
                i_dealt_with.append(j)
                lat.append([factor*float(w) for w in lines[j].split()])
            self["lattice"] = np.array(lat)
            break

        return i_dealt_with

    # Parse atom coordinates from input file lines
    def parse_atoms(self, lines):

        i_dealt_with = []
        for i, line in enumerate(lines):
            
            # Find start of atom block 
            if not line.startswith("atoms"):
                continue

            # Log the lines that have been dealt with during this parsing
            i_dealt_with.append(i)
            
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
                
                i_dealt_with.append(j)
                name = line.split()[0].capitalize()
                pos  = [float(x) for x in line.split()[1:]]
                pos  = np.dot(linv, pos) # Convert to fractional coords
                atoms.append([name,pos])

            self["atoms"] = atoms

        return i_dealt_with

    # Save parameter set to file
    def save(self, filename):

        with open(filename, "w") as f:

            # Write lattice
            lattice = self["lattice"]
            f.write("lattice angstrom\n")
            f.write("{0} {1} {2}\n".format(*lattice[0]))
            f.write("{0} {1} {2}\n".format(*lattice[1]))
            f.write("{0} {1} {2}\n".format(*lattice[2]))
            f.write("\n")

            # Write atoms
            f.write("atoms fractional\n")
            for a in self["atoms"]:
                f.write("{0} {1} {2} {3}\n".format(a[0], *a[1]))
            f.write("\n")

            # Write other parameters
            for p in self.par:
                if p == "atoms":   continue
                if p == "lattice": continue

                f.write("{0} {1}\n".format(p, self.par[p]))

    # Load parameter set from file
    def load(self, filename):

        # Read the file into lines
        with open(filename) as f:
            lines = f.read().split("\n")

        # Tidy lines
        for i in range(0, len(lines)):
        
            # Clean the line, remove comments/whitespace
            line = lines[i]
            line = line.lower()
            for cc in ["#", "//", "!"]:
                if line.find(cc) >= 0:
                    line = line[0:line.find(cc)]
            line = line.strip()
            lines[i] = line

        # Parse special sections
        i_dealt_with = []
        i_dealt_with.extend(self.parse_lattice(lines))
        i_dealt_with.extend(self.parse_atoms(lines))

        # Assume the rest is simple key value form
        for i in range(0, len(lines)):

            # Skip lines already dealt with
            # and blank lines
            if i in i_dealt_with:  continue
            if len(lines[i]) == 0: continue

            # Parse key-value pair
            spl = lines[i].split()
            if len(spl) != 2:
                exept = "Could not parse the line {0} as a key-value pair!"
                raise ValueError(exept.format(lines[i]))

            self[spl[0]] = spl[1]
