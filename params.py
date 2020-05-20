import os, numbers, copy
import numpy          as     np
from   qet.constants  import BOHR_TO_ANGSTROM, ANGSTROM_TO_BOHR
from   qet.type_tools import str_to_type 
from   qet.elements   import elements
from   collections    import defaultdict

# This is thrown when a REQUIRED parameter
# could not be found or generated.
class ParamNotFound(Exception):
    pass

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

        self["outdir"]           = "./"              # outdir = working dir
        self["ibrav"]            = 0                 # no bravis-lattice index
        self["ecutwfc"]          = 60                # plane-wave cutoff (Ry)
        self["occupations"]      = "smearing"        # treat as metallic
        self["degauss"]          = 0.02              # metal smearing width (Ry)
        self["qpoint_spacing"]   = 0.15              # qpoint spacing (2pi A^-1)
        self["force_cube_grids"] = False             # true if grids must be of form NxNxN
        self["kpts_per_qpt"]     = 6                 # ratio of kpt to qpt grid
        self["ldisp"]            = True              # use a grid of q-points
        self["reduce_io"]        = True              # reduce io to a strict minimum
        self["fildvscf"]         = "dvscf"           # potential variation file
        self["electron_phonon"]  = "interpolated"    # electron-phonon method
        self["el_ph_sigma"]      = 0.005             # smearing spacing
        self["el_ph_nsigma"]     = 50                # smearing points
        self["fildyn"]           = "matdyn"          # dynamical matrix prefix
        self["flfrc"]            = "force_constants" # force constants filename
        self["zasr"]             = "simple"          # acoustic sum rule to apply
        self["ph_interp_amt"]    = 8                 # phonon interpolation grid size (as multiple of qpoint_grid)
        self["ndos"]             = 200               # number of energy steps to use when interpolating DOS
        self["ph_interp_prefix"] = "ph_interp"       # the prefix to give to files produced by phonon interpolations
        self["pseudo_dirs"]      = []                # directories to search for pseudopotentials

        # By default, assume cores_per_node is
        # equal to the number of cores where the
        # script is running and nodes = 1
        self["nodes"] = 1
        try:
            import multiprocessing
            self["cores_per_node"] = multiprocessing.cpu_count()
        except ImportError:
            self["cores_per_node"] = 1

        # Add HOME/pseudopotentials and PSEUDO_DIR to the
        # pseudopotential directories to search, if they exist
        if "HOME" in os.environ:
            pd = os.environ["HOME"]+"/pseudopotentials"
            if not pd in self["pseudo_dirs"]: self["pseudo_dirs"].append(pd)
        if "PSEUDO_DIR" in os.environ:
            pd = os.environ["PSEUDO_DIR"]
            if not pd in self["pseudo_dirs"]: self["pseudo_dirs"].append(pd)

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

        # Find the pseudo_dir that contains
        # all of the needed pseudopotentials
        if key == "pseudo_dir":

            # Work out which pseudo_dirs contain
            # which pseudopotentials
            found_in = {}
            for s, m, p in self["species"]:
                found_in[p] = []
                for pd in self["pseudo_dirs"]:
                    if not os.path.isdir(pd): continue
                    if p in os.listdir(pd):
                        found_in[p].append(pd)

            # See if any one pseudo_dir contains
            # all of the needed pseudods
            for pd in self["pseudo_dirs"]:

                has_all = True
                for p in found_in:
                    if not pd in found_in[p]:
                        has_all = False
                        break

                # This pseudo_dir contains all the
                # needed pseudos, go ahead and use it
                if has_all: return pd

            # See if we can combine pseudos from
            # multiple directories
            for p in found_in:

                if len(found_in[p]) == 0:
                    err = "Could not find the pseudopotentail "
                    err += p + " in any of:"
                    for pd in self["pseudo_dirs"]:
                        err += "\n"+pd
                    raise ParamNotFound(err)

            # Create a file with the pseudopotential origin locations
            pof = open("pseudopotential_origin", "w")
             
            # We have found all the pseudos, collect
            # them into the working directory
            for p in found_in:

                # Copy the fist found pseudo to
                # working directory
                os.system("cp "+found_in[p][0]+"/"+p+" .")
                pof.write(p+" from "+found_in[p][0]+"\n")

            pof.close()
            return "./"

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

        # Default to ecutrho = 10*ecutwfc
        if key == "ecutrho": return 10*self["ecutwfc"]

        # Generate the qpoint grid
        if key == "qpoint_grid":
            
            # Generate qpoint grid from spacing
            rlat = np.linalg.inv(self["lattice"]).T
            qps  = float(self["qpoint_spacing"])
            b2q  = lambda b : int(np.linalg.norm(b)/qps)
            grid = [max(1,b2q(b)) for b in rlat]
            if self["force_cube_grids"]: grid = [max(grid) for g in grid]
            return grid

        # Get individual components of qpoint grid
        if key == "nq1": return self["qpoint_grid"][0]
        if key == "nq2": return self["qpoint_grid"][1]
        if key == "nq3": return self["qpoint_grid"][2]

        # Get individial components of interpolated qpoint grid
        if key == "ph_interp_nq1": return self["qpoint_grid"][0]*self["ph_interp_amt"]
        if key == "ph_interp_nq2": return self["qpoint_grid"][0]*self["ph_interp_amt"]
        if key == "ph_interp_nq3": return self["qpoint_grid"][0]*self["ph_interp_amt"]

        # Phonon interpolation output files from prefix
        if key == "ph_interp_dos_file":   return self["ph_interp_prefix"] + ".dos"
        if key == "ph_interp_freq_file":  return self["ph_interp_prefix"] + ".freq"
        if key == "ph_interp_modes_file": return self["ph_interp_prefix"] + ".modes"
        if key == "ph_interp_eig_file":   return self["ph_interp_prefix"] + ".eig"

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
                raise ParamNotFound(msg)

        # This wasn't one of the generatable objects, treat
        # this as a KeyError, so we use the QE default value
        # (if there is one)
        exept = "Key \"{0}\" cold not be generated in parameters object."
        raise KeyError(exept.format(key))
            

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

            # Ensure atom symbols start with a capital letter and don't
            # contain whitespace
            for i in range(0, len(value)):
                value[i][0] = value[i][0].strip().capitalize()

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
    def to_input_line(self, key, name=None):

        # Return nothing if the key is absent
        # so the QE default value will be used
        try:
            val = self[key]
        except ParamNotFound as pnf:
            raise pnf # Something actually went wrong
        except KeyError as e:
            return "" # The key wasn't present, use the default

        # Allow name different from key
        if not name is None: key = name

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
                linv = np.linalg.inv(ANGSTROM_TO_BOHR*self["lattice"].T)
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

    # Parse pseudo
    def parse_pseudo_dirs(self, lines):
        
        i_dealt_with = []
        for i, line in enumerate(lines):
            
            # Only parse pseudo_dir lines
            if not line.startswith("pseudo_dir"):
                continue

            pd = line.split()[1]
            if not pd in self["pseudo_dirs"]:
                self["pseudo_dirs"].append(pd)

            i_dealt_with.append(i)

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

            # Write pseudo_dirs
            pseudo_dirs = self["pseudo_dirs"]
            for pd in pseudo_dirs:
                f.write("pseudo_dir {0}\n".format(pd))
            f.write("\n")

            # Write other parameters
            for p in self.par:
                if p == "atoms":       continue
                if p == "lattice":     continue
                if p == "pseudo_dirs": continue

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
        i_dealt_with.extend(self.parse_pseudo_dirs(lines))

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

    # Load from a CASTEP .cell file
    def load_from_cell(self, filename):
        
        lattice = []
        atoms   = []

        i_done = []
        with open(filename) as f:
            lines = [l.strip().lower() for l in f.read().split("\n")]
            for i, line in enumerate(lines):
                if i in i_done: continue

                if "endblock" in line: 
                    continue

                if "lattice_cart" in line:
                    offset = 1 if "ang" in lines[i+1] else 0
                    for j in range(i+1+offset, i+4+offset):
                        lattice.append(float(w) for w in lines[j].split())
                        i_done.append(j)

                if "positions_frac" in line:
                    for j in range(i+1, len(lines)):
                        i_done.append(j)
                        jline = lines[j]
                        if "endblock" in jline: break
                        atoms.append([jline.split()[0], [float(w) for w in jline.split()[1:4]]])

        self["lattice"] = lattice
        self["atoms"]   = atoms
