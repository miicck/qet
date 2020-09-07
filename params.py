import os, numbers, copy
import numpy          as     np
import time
from   qet.constants  import BOHR_TO_ANGSTROM, ANGSTROM_TO_BOHR
from   qet.type_tools import str_to_type 
from   qet.elements   import elements
from   qet.logs       import log
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
    
    # The first time.time() that __init__ was called
    first_init_time = None

    # Default constructor
    def __init__(self, filename=None):

        # Record the first initialization time
        if parameters.first_init_time is None:
            parameters.first_init_time = time.time()
            log("First init time set to {0}".format(parameters.first_init_time))

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
        self["min_q_grid_size"]  = 2                 # The minimum q-points to a side for a q-point grid
        self["force_cube_grids"] = False             # true if grids must be of form NxNxN
        self["kpts_per_qpt"]     = 10                # ratio of kpt to qpt grid
        self["ldisp"]            = True              # use a grid of q-points
        self["reduce_io"]        = True              # reduce io to a strict minimum
        self["fildvscf"]         = "dvscf"           # potential variation file
        self["electron_phonon"]  = "interpolated"    # electron-phonon method
        self["el_ph_sigma"]      = 0.0025            # smearing spacing
        self["el_ph_nsigma"]     = 50                # smearing points
        self["fildyn"]           = "matdyn"          # dynamical matrix prefix
        self["flfrc"]            = "force_constants" # force constants filename
        self["zasr"]             = "simple"          # acoustic sum rule to apply
        self["ph_interp_amt"]    = 8                 # phonon interpolation grid size (as multiple of qpoint_grid)
        self["ndos"]             = 200               # number of energy steps to use when interpolating DOS
        self["ph_interp_prefix"] = "ph_interp"       # the prefix to give to files produced by phonon interpolations
        self["pseudo_dirs"]      = []                # directories to search for pseudopotentials
        self["bz_path_points"]   = 100               # the approximate number of points along a BZ path
        self["total_walltime"]   = 129600            # allocated walltime in seconds (negative => treat as infinite)
        self["tidy_up_time"]     = 1800              # time reserved for tidying up at end of total_walltime

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

        # Output timing information
        log("Initialized parameters object at time {0}".format(time.time()))
        log("    first_init_time : {0}".format(parameters.first_init_time))
        log("    total_walltime  : {0}".format(self["total_walltime"]))
        log("    end_time        : {0}".format(self["end_time"]))
        log("    max_seconds     : {0}".format(self["max_seconds"]))

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

        # Work out the space group
        if key == "space_group":
            self.eval_symmetry()
            return self["space_group"]

        # Work out the number of symmetry operations
        if key == "sym_ops":
            self.eval_symmetry()
            return self["sym_ops"]

        # Work out a good BZ path
        if key == "bz_path" or key == "high_symmetry_bz_points":
            try: import seekpath
            except ImportError:
                raise ImportError("Could not import SeeKpath!")

            # Convert the structure into a form that SeeKpath can digest
            frac_coords  = []
            atom_numbers = []
            unique_names = []

            for a in self["atoms"]:
                if not a[0] in unique_names:
                    unique_names.append(a[0])

            for a in self["atoms"]:
                frac_coords.append(a[1])
                atom_numbers.append(unique_names.index(a[0]))
            

            # Call SeeKpath to get the BZ path
            structure = (self["lattice"], frac_coords, atom_numbers)
            path = seekpath.get_path(
                structure, 
                with_time_reversal=True,
                symprec=0.001,
                angle_tolerance=0.5,
                threshold=0)

            # Work out how many points we have along each segment of the BZ path
            pc          = path["point_coords"]
            segs        = [[np.array(pc[p[0]]), np.array(pc[p[1]])] for p in path["path"]]
            seg_names   = [[p[0], p[1]] for p in path["path"]] 
            seg_lengths = [np.linalg.norm(s[1]-s[0]) for s in segs]
            tot_length  = sum(seg_lengths)
            seg_counts  = [max(int(self["bz_path_points"]*l/tot_length),2) for l in seg_lengths]

            kpoints = [] # Will contain the k-points in the path
            high_symm_points = {} # Will contain the names of high-symmetry points along the path

            for i, c in enumerate(seg_counts):
                pts = np.linspace(0.0, 1.0, c)
                pts = [segs[i][0] + (segs[i][1]-segs[i][0])*p for p in pts]
                high_symm_points[len(kpoints)] = seg_names[i][0]
                kpoints.extend(pts)

                if i + 1 < len(seg_names) and seg_names[i][1] == seg_names[i+1][0]:
                    kpoints.pop() # Remove repeated high symmetry point
                else:
                    high_symm_points[len(kpoints)-1] = seg_names[i][1]

            self["high_symmetry_bz_points"] = high_symm_points
            self["bz_path"] = kpoints

            if key == "bz_path": return kpoints
            else: return high_symm_points

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
            
            # Generate qpoint grid from spacing/min_q_grid_size
            rlat = np.linalg.inv(self["lattice"]).T
            qps  = float(self["qpoint_spacing"])
            b2q  = lambda b : int(np.linalg.norm(b)/qps)
            grid = [max(self["min_q_grid_size"],b2q(b)) for b in rlat]
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

        # Default to k-point parallelism
        if key == "pools" : return self["cores_per_node"]*self["nodes"]
        if key == "images": return 1

        # Get the default k-point grids for calculating Tc
        if key == "tc_kpqs" : return [self["kpts_per_qpt"]-1, self["kpts_per_qpt"]]

        # Default q-e bin/ path = environment path
        if key == "path_override" : return ""

        # If total_walltime is > 0, this will be the time returned by
        # time.time() when we will run out of walltime. Otherwise is -1.
        if key == "end_time":
            if self["total_walltime"] < 0: 
                return -1
            else: 
                return parameters.first_init_time + self["total_walltime"] 

        # Returns the maximum seconds we let a calculation run for
        # from this moment in time. Is equal to the end_time minus
        # the current time minus the tidy_up_time.
        if key == "max_seconds":
            if self["end_time"] < 0:
                return 10e7 # No end time => essentially infinite time
            return max(self["end_time"] - time.time() - self["tidy_up_time"], 0)

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

    # Check if a key is contained within this parameters object
    # optionally checking to see if it can be generated
    def contains_key(self, key, generate=True):
        if key in self.par:
            return True
        if not generate: return False
        try:
            self.gen_param(key)
            return True
        except:
            return False

    # Validate and standardise a given parameter
    # based on it's key and value
    def validate_and_standardise_param(self, key, value):
        
        # Try to convert a string value to the correct type
        if isinstance(value, str):
            value = str_to_type(value)

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
        self.par[key] = self.validate_and_standardise_param(key, value)

    # Construct special input lines that do not fit the type-based
    # construction used in self.to_input_line
    def special_input_line(self, key, value):
        
        # Alpha_mix is an array that defaults to the first entry
        if key == "alpha_mix":
            return "    alpha_mix(1) = {0},\n".format(value)

        return None

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

        # Check for special input line
        sf = self.special_input_line(key, val)
        if not sf is None: return sf

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

    # Parse pseudopotential directories from input file
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

    # Parse k-point grid sizes to use for Tc calculation
    def parse_tc_kpqs(self, lines):
        
        i_dealt_with = []
        parsed_kpqs  = []
        for i, line in enumerate(lines):
            
            if not line.startswith("tc_kpqs"):
                continue

            parsed_kpqs.extend([int(w) for w in line.split()[1:]])
            i_dealt_with.append(i)

        if len(parsed_kpqs) > 0:
            self["tc_kpqs"] = parsed_kpqs

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
        i_dealt_with.extend(self.parse_tc_kpqs(lines))

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

    # Use c2x to evaluate the symmetry of the crystal
    def eval_symmetry(self):
        
        # Create a temporary .cell file for c2x to parse
        TMP_CELL = "tmp_for_symm.cell"
        TMP_SYMM = "tmp_for_symm.symm"
        self.save_to_cell(TMP_CELL)

        # Work out the number of symmetry ops
        os.system("c2x --symmetry "+TMP_CELL+" > /dev/null 2>"+TMP_SYMM)
        with open(TMP_SYMM) as f:
            self["sym_ops"] = int(f.read().split()[1])

        # Work out the space group
        os.system("c2x --int "+TMP_CELL+" 2>"+TMP_SYMM)
        with open(TMP_SYMM) as f:
            self["space_group"] = f.read().split()[-1]

        # Remove temporary files
        os.system("rm "+TMP_CELL+" "+TMP_SYMM)
        time.sleep(0.25) # Wait for filesystem to catch up

    # Save the geometry to a CASTEP .cell file
    def save_to_cell(self, filename):
        with open(filename,"w") as f:
            f.write("%BLOCK LATTICE_CART\nANG\n")
            for l in self["lattice"]:
                f.write(" ".join([str(w) for w in l])+"\n")
            f.write("%ENDBLOCK LATTICE_CART\n\n")
            f.write("%BLOCK POSITIONS_FRAC\n")
            for a in self["atoms"]:
                f.write("{0} {1} {2} {3}\n".format(a[0], *a[1]))
            f.write("%ENDBLOCK POSITIONS_FRAC\n")

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
                        lattice.append([float(w) for w in lines[j].split()])
                        i_done.append(j)

                if "positions_frac" in line:
                    for j in range(i+1, len(lines)):
                        i_done.append(j)
                        jline = lines[j]
                        if "endblock" in jline: break
                        atoms.append([jline.split()[0], [float(w) for w in jline.split()[1:4]]])

        self["lattice"] = lattice
        self["atoms"]   = atoms
