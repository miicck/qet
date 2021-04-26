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
        self["mpirun"]           = "mpirun"          # The mpi runner that we want to use

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


    ####################
    # PARAMETER ACCESS #
    ####################


    # Get parameter values with []
    def __getitem__(self, key):

        # Attempt to generate the
        # parameter if key if not found
        if not key in self.par:
            return self.gen_param(key)

        return self.par[key]

    # Set parameter values with []
    def __setitem__(self, key, value):
        self.par[key] = self.validate_and_standardise_param(key, value)

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

        if key == "a": return np.linalg.norm(self["lattice"][0])
        if key == "b": return np.linalg.norm(self["lattice"][1])
        if key == "c": return np.linalg.norm(self["lattice"][2])

        if key == "alpha":
            lat = self["lattice"]
            ret = np.dot(lat[1], lat[2])
            ret = np.arccos(ret)*180/np.pi
            return ret

        if key == "beta":
            lat = self["lattice"]
            ret = np.dot(lat[0], lat[2])
            ret = np.arccos(ret)*180/np.pi
            return ret

        if key == "gamma":
            lat = self["lattice"]
            ret = np.dot(lat[0], lat[1])
            ret = np.arccos(ret)*180/np.pi
            return ret

        # Work out the space group
        if key == "space_group_name":
            self.eval_symmetry()
            return self["space_group_name"]

        # Work out the number of symmetry operations
        if key == "sym_ops":
            self.eval_symmetry()
            return self["sym_ops"]

        # Work out a good BZ path
        if key == "bz_path" or key == "high_symmetry_bz_points":
            try: import seekpath
            except ImportError:
                log("Could not import seekpath!")
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
            

    #########################
    # INPUT FILE GENERATION #
    #########################


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

    ############
    # SYMMETRY #
    ############

    # Use c2x to evaluate the symmetry of the crystal
    # (c2x is available at https://www.c2x.org.uk/)
    def eval_symmetry(self):
        
        # Create a temporary .cell file for c2x to parse
        TMP_CELL = "tmp_for_symm.cell"
        TMP_SYMM = "tmp_for_symm.symm"
        self.save_to_cell(TMP_CELL)

        # Work out the number of symmetry ops
        os.system("c2x --symmetry "+TMP_CELL+" > /dev/null 2>"+TMP_SYMM)
        with open(TMP_SYMM) as f:
            read = f.read().strip()
            if "\n" in read: read = read.split("\n")[-1]
            self["sym_ops"] = int(read.split()[1])

        # Work out the space group
        os.system("c2x --int "+TMP_CELL+" 2>"+TMP_SYMM)
        with open(TMP_SYMM) as f:
            self["space_group_name"] = f.read().split()[-1]

        # Remove temporary files
        os.system("rm "+TMP_CELL+" "+TMP_SYMM)
        time.sleep(0.25) # Wait for filesystem to catch up


    # For a wavevector k_p in the B.Z of this crystal, returns 
    # a supercell of this crystal for which the wavevector
    # folds back to the gamma point (i.e returns a supercell
    # that can support explicit pertubations with wavevector k_p).
    # If apply_eigenvector is set, and has dimensions 
    # #atoms x 3, then a phonon mode of the form eigenvector * cos(k_p dot R)
    # will be applied to the resulting structure. Apply_eigenvector is assumeed 
    # to be in units of angstrom.
    def generate_commensurate_supercell(self, k_p, apply_eigenvector=None, return_disp=False):
        import math
        import numpy as np
        from fractions import Fraction

        if len(k_p) != 3:
            raise Exception("Expected k-vector of length 3!")

        # The fundamental logic here is described in
        # 10.1103/PhysRevB.92.184301
        # We essentially brute-force search the 
        # Hermite-normal form (HNF) supercell matricies S 
        # until we find k_s = S k_p such that k_s has 
        # integer entries
        is_int = lambda v : abs(v - round(v)) < 10e-4
        possibilities = []

        # Find the lowest common multiple of the
        # denominators of the fractional entries in k_p
        # (this is the determinant of the supercell matrix needed)
        lcm = None
        for i in range(0, 3):
            f = Fraction(k_p[i]).limit_denominator(100)
            d = f.denominator
            if lcm is None: lcm = d
            else: lcm = lcm * d // math.gcd(lcm, d)

        # Loop over all diagonals of S
        # which give the desired determinant (lcm)
        for s11 in range(1, lcm+1):
            for s22 in range(1, lcm//s11+1):
                s33 = lcm//s11//s22
                if s11*s22*s33 != lcm: continue

                # Loop over off-diagonals that satisfy HNF
                for s12 in range(0, s22):
                    for s13 in range(0, s33):
                        for s23 in range(0, s33):

                            # Evaluate k_s and check if it is integer
                            ks1 = s11*k_p[0] + s12*k_p[1] + s13*k_p[2]
                            if not is_int(ks1): continue
                            ks2 = s22*k_p[1] + s23*k_p[2]
                            if not is_int(ks2): continue
                            ks3 = s33*k_p[2]
                            if not is_int(ks3): continue

                            # Record this supercell
                            possibilities.append([s11, s12, s13, s22, s23, s33])

        # Algorithm failed
        if len(possibilities) == 0:
            fs = "Could not find supercell of size {0} comensurate with k = {1}!"
            raise Exception(fs.format(lcm, k_p))

        norm_sq      = lambda v : v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
        lowest_score = float("inf")
        lat          = self["lattice"]
        ss           = None
        ss_lat       = None

        # Find the possibility with the smallest 
        # average basis vector size
        for p in possibilities:

            # Work out supercell basis vectors
            v1 = p[0]*lat[0] + p[1]*lat[1] + p[2]*lat[2]
            v2 = p[3]*lat[1] + p[4]*lat[2]
            v3 = p[5]*lat[2]

            # Record if score is lowest so far
            score = norm_sq(v1) + norm_sq(v2) + norm_sq(v3)
            if score < lowest_score:
                lowest_score = score
                ss = p
                ss_lat = [v1, v2, v3]

        superlat_inv = np.linalg.inv(np.array(ss_lat).T)
        max_delta = max(ss)*4
        frac_eps  = 1e-4
        ss_atoms  = []
        ss_disps  = []

        # Work out where the atoms go in the supercell
        for ai, a in enumerate(self["atoms"]):

            # The position of the atom in the first unit cell in cartesians
            cart_pos = a[1][0]*lat[0]+a[1][1]*lat[1]+a[1][2]*lat[2]

            # Loop over images of the pimitive cell
            for dx in range(-max_delta, max_delta+1):
                for dy in range(-max_delta, max_delta+1):
                    for dz in range(-max_delta, max_delta+1):

                        # The position of the image of the atom in cartesians
                        # and in fractional supercell coordinates
                        image_pos = cart_pos + dx*lat[0] + dy*lat[1] + dz*lat[2]

                        # Apply phonon pertubation if eigenvector is provided
                        if not apply_eigenvector is None:
                            phase = np.cos(np.dot(k_p, [dx+a[1][0], dy+a[1][1], dz+a[1][2]])*np.pi*2)
                            delta = phase * np.array(apply_eigenvector[ai]) 
                            image_pos += delta

                        # Get fraction coordinates in the supercell
                        frac_pos_sup = np.dot(superlat_inv, image_pos)

                        # If x is in the first supercell, add it to the supercell atoms
                        # (exclude atoms at 1.0 because we're including atoms at 0.0)
                        if all([x > -frac_eps and x < 1-frac_eps for x in frac_pos_sup]):
                            ss_atoms.append([a[0], frac_pos_sup])
                            ss_disps.append(delta)

        # Build the supercell
        sup = self.copy()
        sup["lattice"] = ss_lat
        sup["atoms"]   = ss_atoms

        # Check that the supercell contains the correct number of atoms
        if len(self["atoms"])*lcm != len(sup["atoms"]):
            print("Length of supercell = {0} (target = {1})".format(len(sup["atoms"]), len(self["atoms"])*lcm))
            raise Exception("Supercell does not contain the correct number of atoms!")

        self.eval_symmetry()
        sup.eval_symmetry()
        print("Symmetry before: {0} Symmetry after: {1}".format(self["space_group_name"], sup["space_group_name"]))

        ret = [sup]
        if return_disp: ret.append(ss_disps)
        if len(ret) == 1: return ret[0]
        return ret

    # Apply a phonon of the given wavevector (q) and amplitude (amp).
    # Returns a parameters object describing the perturbed cell
    # the target energy is in meV/atom
    def apply_phonon(self, q, mode, target_energy=1.0, modes_file="ph_interp.modes", use_closest=False, return_freq=False, return_evec=False, return_disp=False):

        # Parse the modes file for the eigenvectors
        from qet.parser import phonon_interp_modes
        modes = phonon_interp_modes(modes_file)

        qpts  = modes["q-points"]
        evecs = modes["eigenvectors"]
        freqs = modes["frequencies"]

        # Find the closest q-point to the requested one
        disp = lambda q1, q2 : sum((q1[i]-q2[i])*(q1[i]-q2[i]) for i in range(3))

        min_disp = float("inf")
        min_i    = None
        for i in range(len(qpts)):
            
            d = disp(qpts[i], q)
            if d < min_disp:
                min_disp = d
                min_i = i

        # Error out if we are not using the closest q-point to that specified
        if not use_closest:
            if min_disp > 10e-6:
                fs  = "Could not find modes for the q-point {0}, the closest was "
                fs += "{1}. Set use_closest=True to approximate the eigenvector at {0} "
                fs += "with the eigenvector at {1}."
                raise Exception(fs.format(q, qpts[min_i]))
        
        # Find the eigenvector corresponding to the 
        # requested (q-point, mode) combination
        modes_per_q = int(len(evecs)/len(qpts))
        if mode < 0 or mode >= modes_per_q:
            raise Exception("Mode index {0} is out of range!".format(mode))

        evec = evecs[min_i*modes_per_q+mode]
        if modes_per_q != len(evec)*3:
            raise Exception("#modes != #atoms * 3")
            
        # Work out the displacement that corresponds to
        # the requested energy

        # Get the frequency of this phonon
        freq = freqs[min_i*modes_per_q+mode]

        # Calculate d^T M d where d is the displacement
        # and m is the mass tensor
        dmd = 0
        masses = {}
        for s in self["species"]:
            masses[s[0]] = s[1]

        for i, v in enumerate(evec):
            m = masses[self["atoms"][i][0]] 
            dmd += m * sum([x*x for x in evec[i]])

        # Work out the amplitude needed to obtain the target energy
        # (I'm afraid of units, so work in S.I)
        from qet.constants import EV_TO_J, AMU_TO_KG, CMM1_TO_HZ, ANGSTROM_TO_M

        # Convert energy to joules per supercell (from meV/atom)
        target_energy = target_energy * len(evec) * EV_TO_J / 1000.0

        # Convert mass to KG and length from Angstrom to M
        # (using angstrom here so that the displacement we pass to 
        #  generate_commensurate_supercell will be in angstrom)
        dmd = dmd * AMU_TO_KG * ANGSTROM_TO_M * ANGSTROM_TO_M

        # Convert freq to HZ
        freq_cmm = freq
        freq = freq * CMM1_TO_HZ 

        # Work out amplitude by inverting E = 0.5 w^2 amp^2 dmd
        amp  = 2 * target_energy / (dmd * freq * freq)
        amp  = amp ** 0.5

        if amp > 10:
            print("Phonon amplitude is large ({0})".format(amp))
            print("Perhaps the phonon frequency ({0} cm^-1) is too small".format(freq_cmm))
            print("Setting the amplitude to zero")
            amp = 0.0

        evec = [[x*amp for x in e] for e in evec]

        ss, disp = self.generate_commensurate_supercell(q, apply_eigenvector=evec, return_disp=True)
        ret = [ss]
        if return_freq: ret.append(freq_cmm)
        if return_evec: ret.append(evec)
        if return_disp: ret.append(disp)
        if len(ret) == 1: return ret[0]
        return ret


    ####################
    # SAVING / LOADING #
    ####################


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

                # Already written
                if p == "atoms":       continue 
                if p == "lattice":     continue
                if p == "pseudo_dirs": continue

                # Don't write
                if p == "bz_path":                 continue
                if p == "high_symmetry_bz_points": continue

                val = self.par[p]

                # Space-seperated list formatting
                if isinstance(val, list):
                    val = " ".join(str(x) for x in val)

                f.write("{0} {1}\n".format(p, val))

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

    # Load from a quantum espresso .in file
    def load_from_qe_input(self, filename):

        lattice = []
        atoms   = []
        
        with open(filename) as f:
            lines = [l.strip().lower() for l in f.read().split("\n")]
            for i, line in enumerate(lines):
                if "cell_parameters" in line:
                    for j in range(i+1, i+4):
                        lattice.append([float(w) for w in lines[j].split()])

                if "atomic_positions" in line:
                    for j in range(i+1, len(lines)):
                        try:
                            name, x, y, z = lines[j].split()
                            atoms.append([name, [float(x), float(y), float(z)]])
                        except:
                            break

        self["lattice"] = lattice
        self["atoms"]   = atoms
