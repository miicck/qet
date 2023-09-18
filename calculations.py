import os
import time
import subprocess
import psutil
import threading
import qet.parser as parser
from   qet.logs   import log

# Pad the equals signs in the input file
# so they line up nicely :)
def pad_input_file(s):
    
    lines = s.split("\n")
    maxi = max([l.index("=") for l in lines if "=" in l])
    snew = ""
    for l in lines:
        if "=" in l:
            i = l.index("=")
            pad = " "*(maxi-i)
            l = l[0:i] + pad + l[i:]
        snew += l + "\n"

    # Ensure there is a newline at the end, otherwise QE might crash 
    return snew.strip() + "\n" 

# Get the geompetry part of an input file
def input_geometry(params, explicit_kpoints=None, kpoint_labels={}):

    # S will eventually contain the text for
    # the geometry part of the input
    s = ""

    # Cell parameters card
    s += "CELL_PARAMETERS (angstrom)\n"

    exp = params["expand"]
    s += "{0} {1} {2}\n".format(*[exp*x for x in params["lattice"][0]])
    s += "{0} {1} {2}\n".format(*[exp*x for x in params["lattice"][1]])
    s += "{0} {1} {2}\n".format(*[exp*x for x in params["lattice"][2]])

    # Atomic species card
    s += "\nATOMIC_SPECIES\n"
    for sp in params["species"]:
        s += "{0} {1} {2}\n".format(*sp)

    # Atomic positions card
    if params.contains_key("space_group"):
        s += "\nATOMIC_POSITIONS (crystal_sg)\n"
    else:
        s += "\nATOMIC_POSITIONS (crystal)\n"
    for a in params["atoms"]:
        s += "{0} {1} {2} {3}\n".format(a[0], *a[1])

    if explicit_kpoints is None:

        # automatic K-points card
        s += "\nK_POINTS automatic\n"
        s += "{0} {1} {2} 0 0 0\n".format(*params["kpoint_grid"])

    else:
        
        # explicit K-points card
        s += "\nK_POINTS (crystal)\n"
        w  = 1.0/float(len(explicit_kpoints))
        s += "{0}\n".format(len(explicit_kpoints))
        for i, k in enumerate(explicit_kpoints):
            s += "{0} {1} {2} {3}".format(k[0],k[1],k[2],w)
            if i in kpoint_labels: s += " ! "+kpoint_labels[i]
            s += "\n"

    return s

# Get the control part of a pw.x input file
def pw_control_input(params, calculation="scf", recover=False):

    # Control namelist
    s  = "&CONTROL\n"
    s += "    calculation = '{0}',\n".format(calculation)
    if recover: s += "    restart_mode = 'restart',\n"
    s += params.to_input_line("outdir")
    s += params.to_input_line("pseudo_dir")
    s += params.to_input_line("forc_conv_thr")
    s += params.to_input_line("max_seconds")
    s += "/\n\n"

    # System namelist
    s += "&SYSTEM\n"
    s += params.to_input_line("ntyp")
    s += params.to_input_line("nat")

    explicit_symm = False
    if params.contains_key("space_group"):
        explicit_symm = True
        s += params.to_input_line("space_group")
    elif params.contains_key("ibrav"):
        if params["ibrav"] != 0:
            explicit_symm = True
        s += params.to_input_line("ibrav")

    if explicit_symm:
        s += params.to_input_line("a")
        s += params.to_input_line("b")
        s += params.to_input_line("c")
        s += params.to_input_line("cosab")
        s += params.to_input_line("cosbc")
        s += params.to_input_line("cosac")

    s += params.to_input_line("ecutwfc")
    s += params.to_input_line("ecutrho")
    s += params.to_input_line("occupations")
    s += params.to_input_line("smearing")
    s += params.to_input_line("degauss")
    s += params.to_input_line("la2F")
    s += params.to_input_line("nr1")
    s += params.to_input_line("nr2")
    s += params.to_input_line("nr3")
    s += params.to_input_line("nr1s")
    s += params.to_input_line("nr2s")
    s += params.to_input_line("nr3s")
    s += params.to_input_line("nosym")
    s += "/\n\n"

    # Electrons namelist
    s += "&ELECTRONS\n"
    s += params.to_input_line("mixing_beta")
    s += params.to_input_line("conv_thr")
    s += "/\n\n"


    if "relax" in calculation:

        # Ions namelist
        s += "&IONS\n"
        s += params.to_input_line("ion_dynamics")
        s += "/\n\n"

    if calculation == "vc-relax":

        # Cell namelist
        s += "&CELL\n"
        s += params.to_input_line("cell_factor")
        s += params.to_input_line("cell_dynamics")
        s += params.to_input_line("press")
        s += params.to_input_line("press_conv_thr")
        s += "/\n\n"

    return s

# Returns True if the given file is a completed
# quantum espresso output
def is_complete(filename):

    # Doesn't exist => not complete
    if not os.path.isfile(filename):
        return False

    # JOB DONE => complete
    with open(filename) as f:
        for line in f:

            # If this appears before JOB DONE, not done
            if "Maximum CPU time exceeded" in line:
                return False

            # No convergence => not done
            if "No convergence has been achieved" in line:
                return False

            # JOB DONE, and no obvious reason to disbelieve
            if "JOB DONE" in line:
                return True

    return False
     

# A thread for tracking cpu/memory
# usage, which can be stopped
class cpu_tracking_thread(threading.Thread):

    def __init__(self, filename="cpu_ram_usage",  *args, **kwargs):
        super(cpu_tracking_thread, self).__init__(*args, **kwargs)
        self.stop_event = threading.Event()
        self.filename = filename
        self.start_time = time.time()

    def stop(self):
        self.stop_event.set()

    def stopped(self):
        return self.stop_event.is_set()

    def record(self, file):
        fs = "{0}, {1}, {2}\n"
        mem = psutil.virtual_memory().percent
        cpu = psutil.cpu_percent()
        tim = round(time.time()-self.start_time, 1)
        file.write(fs.format(tim, mem, cpu))
        file.flush()

    def run(self):
        
        # Time between records
        REPEAT_TIME = 10

        # Open the logging file
        with open(self.filename, "w") as f:
            f.write("Time (s) memory (%) cpu (%)\n")

            while True:

                # Record every REPEAT_TIME seconds
                time.sleep(REPEAT_TIME)
                self.record(f)
                if self.stopped(): return

##################
#  CALCULATIONS  #
##################

# Calculation base class, calclation types derive from this
class calculation:
    
    # Contructor providing parameters
    def __init__(self, in_params):
        
        # Save the input parameters
        self.in_params = in_params

    # Get the string reperesentation of this
    # calculation <=> input file
    def __str__(self):
        s  = self.default_filename()+".in:\n"
        s += self.gen_input_file()
        return s

    # Should return true if image-based parallelization is
    # allowed for this calculation type (q.e will error out
    # if this isn't correct for the exe)
    def image_parallelization_allowed(self):
        return False
    
    # Run the calculation in the given directory
    # with the given name, will check if the
    # calculation exists and, if so, will attempt
    # to continue it if it is incomplete
    def run(self, filename=None, path="./", required=True):

        if filename is None:
            filename = self.default_filename()

        log("Starting {0} {1} calculation ({2}) at {3}:".format(
            "required" if required else "optional", self.exe(), filename, path))
        
        # Remove trailing /'s from path
        path = path.strip()
        path = path.rstrip("/")

        inf  = path+"/"+filename+".in"
        outf = path+"/"+filename+".out"
        
        # Test to see if the calculation is complete
        if is_complete(outf):

            # Log that we are skipping this complete calculation
            msg = "Calculation \"{0}\" is complete, skipping..."
            log(msg.format(outf))
            return self.parse_output(outf)

        else: # Calculation not complete

            # Create input file, run calculation
            recover = os.path.isfile(outf)

            log("Generating input file...")
            with open(inf, "w") as f:
                f.write(self.gen_input_file(recover=recover))

            log("Setting up command to run...")

            # Get number of processes
            np  = self.in_params["cores_per_node"]*self.in_params["nodes"]
            ppn = self.in_params["cores_per_node"]

            # Setup parallelism scheme
            if self.image_parallelization_allowed():
                # Use k-point/image parallelism
                pools  = self.in_params["pools"]
                images = self.in_params["images"]
                qe_flags = "-nk {0} -ni {1}".format(pools, images) 
            else:
                # Just use k-point parallelism
                qe_flags = "-nk {0}".format(np)

            # Apply overriden q-e location
            bin = self.in_params["path_override"]
            if len(bin) > 0: bin += "/"

            # Get the thing that we will run mpi with
            mpirun = self.in_params["mpirun"]

            try:
                # Check if mpirun accepts -ppn flag
                subprocess.check_output([mpirun, "-ppn", "1", "ls"])
                cmd = "cd {0}; {1} -ppn {2} -np {3} {4} {5} -i {6} > {7}"
                cmd = cmd.format(path, mpirun, ppn, np, bin + self.exe(), qe_flags, inf, outf)
            except:
                # Doesn't accept -ppn flag
                cmd = "cd {0}; {1} -np {2} {3} {4} -i {5} > {6}"
                cmd = cmd.format(path, mpirun, np, bin + self.exe(), qe_flags, inf, outf)

            log("Running:\n"+cmd)

            # Start tracking thread
            tracking = cpu_tracking_thread(filename=filename+".cpu")
            tracking.start()

            try:
                # Run calculation, log stdout/stderr
                stdout = subprocess.check_output([cmd], stderr=subprocess.STDOUT, shell=True)
                log(stdout, filename="qet.stdout")
            except subprocess.CalledProcessError as e:
                # Log subprocess errror
                log(e)

            tracking.stop()

            # Check for success
            if not is_complete(outf):

                if required:
                    msg = "Calculation {0} did not complete, stopping!"
                    msg = msg.format(outf)
                    log(msg)
                    raise RuntimeError(msg)
                else:
                    msg = "Calculation {0} did not complete, but isn't required, continuing..."
                    log(msg.format(outf))
                    return None

            else:

                # Parse the output
                return self.parse_output(outf)

# A simple SCF calculation
class scf(calculation):

    # The executable that carries out this calculation
    def exe(self):
        return "pw.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "scf"

    # Parse calculation output
    def parse_output(self, outf):
        return parser.scf_out(outf)

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        # Control
        s = pw_control_input(self.in_params, 
            calculation="scf", recover=recover)

        # Geometry
        s += input_geometry(self.in_params)

        return pad_input_file(s)


class bands(calculation):
    
    # The executable that carries out this calculation
    def exe(self):
        return "pw.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "e_bands"

    # Parse calculation output
    def parse_output(self, outf):
        return None

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):
        
        # Control
        s = pw_control_input(self.in_params, 
            calculation="scf", recover=recover)

        # Geometry
        s += input_geometry(self.in_params, 
             explicit_kpoints=self.in_params["bz_path"], 
             kpoint_labels=self.in_params["high_symmetry_bz_points"])

        return pad_input_file(s)


class extract_bands(calculation):
    
    # The executable that carries out this calculation
    def exe(self):
        return "bands.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "e_bands_extract"

    # Parse calculation output
    def parse_output(self, outf):
        return None

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):
        
        s  = "&BANDS\n"
        s += self.in_params.to_input_line("outdir")
        s += "/\n"
        return pad_input_file(s)


# Calculate the electronic DOS
class dos(calculation):
    
    # The executable that carries out this calculation
    def exe(self):
        return "dos.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "dos"

    # Parse calculation output
    def parse_output(self, outf):
        return parser.dos_out(outf)

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        s  = "&DOS\n"
        s += self.in_params.to_input_line("outdir")
        s += self.in_params.to_input_line("degauss")
        s += self.in_params.to_input_line("DeltaE")
        s += "/\n"

        return pad_input_file(s)

# Calculate the atom-projected electronic DOS
class proj_dos(calculation):
    
    # The executable that carries out this calculation
    def exe(self):
        return "projwfc.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "proj_dos"

    # Parse calculation output
    def parse_output(self, outf):
        return parser.proj_dos_out(outf)

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        s  = "&PROJWFC\n"
        s += self.in_params.to_input_line("outdir")
        s += self.in_params.to_input_line("degauss")
        s += self.in_params.to_input_line("DeltaE")
        s += "/\n"

        return pad_input_file(s)

# Carry out a geometry optimization (vc-relax)
class relax(calculation):

    # The executable that carries out this calculation
    def exe(self):
        return "pw.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "relax"

    # Parse calculation output
    def parse_output(self, outf):
        return parser.relax_out(outf)

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        # Control
        s  = pw_control_input(self.in_params, 
            calculation=self.relax_calculation(), 
            recover=recover)

        # Geometry
        s += input_geometry(self.in_params)

        return pad_input_file(s)

    def relax_calculation(self):
        return "vc-relax"

# Overloaded version of relax, to carry out relaxation
# of just the atomic posiitons
class fixed_cell_relax(relax):

    def relax_calculation(self):
        return "relax"

# Calculate phonons on a grid
class phonon_grid(calculation):

    # The executable that carries out this calculation
    def exe(self):
        return "ph.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "phonons"

    # Images are allowed for phonon calculations
    def image_parallelization_allowed(self):
        return True

    # Parse calculation output
    def parse_output(self, outf):
        return parser.phonon_grid_out(outf)

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        qpg = self.in_params["qpoint_grid"]

        s  = "Calculate phonons on a coarse grid\n"
        s += "&INPUTPH\n"
        s += self.in_params.to_input_line("outdir")
        s += self.in_params.to_input_line("ldisp")
        s += self.in_params.to_input_line("reduce_io")
        s += self.in_params.to_input_line("max_seconds")
        s += self.in_params.to_input_line("nq1")
        s += self.in_params.to_input_line("nq2")
        s += self.in_params.to_input_line("nq3")
        s += self.in_params.to_input_line("search_sym")
        if recover: s += "    recover = .true.,\n"
        s += "/\n"

        # I've found ph.x sometimes crashes if 
        # the input file doesn't end in a blank line
        return pad_input_file(s) + "\n"

# Calculate electron-phonon coupling on a grid
class electron_phonon_grid(calculation):
    
    # The executable that carries out this calculation
    def exe(self):
        return "ph.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "elph"

    # Images are allowed for phonon calculations
    def image_parallelization_allowed(self):
        return True

    # Parse calculation output
    def parse_output(self, outf):
        return None

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        qpg = self.in_params["qpoint_grid"]

        s  = "Calculate electron-phonon coupling on a coarse grid\n"
        s += "&INPUTPH\n"
        s += self.in_params.to_input_line("outdir")
        s += self.in_params.to_input_line("ldisp")
        s += self.in_params.to_input_line("reduce_io")
        s += self.in_params.to_input_line("fildvscf")
        s += self.in_params.to_input_line("electron_phonon")
        s += self.in_params.to_input_line("max_seconds")
        s += self.in_params.to_input_line("nq1")
        s += self.in_params.to_input_line("nq2")
        s += self.in_params.to_input_line("nq3")
        s += self.in_params.to_input_line("el_ph_sigma")
        s += self.in_params.to_input_line("el_ph_nsigma")
        s += self.in_params.to_input_line("search_sym")
        s += self.in_params.to_input_line("alpha_mix")
        if recover: s += "    recover = .true.,\n"
        s += "/\n"

        # I've found ph.x sometimes crashes if 
        # the input file doesn't end in a blank line
        return pad_input_file(s) + "\n"

# Convert reciprocal space quantities to real-space quantitites
class q2r(calculation):
    
    # The executable that carries out this calculation
    def exe(self):
        return "q2r.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "q2r"

    # Images are allowed for phonon calculations
    def image_parallelization_allowed(self):
        return True

    # Parse calculation output
    def parse_output(self, outf):
        return None

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        s  = "&INPUT\n"
        s += self.in_params.to_input_line("zasr")
        s += self.in_params.to_input_line("fildyn")
        s += self.in_params.to_input_line("flfrc")
        s += self.in_params.to_input_line("la2F")
        s += "/\n"

        return pad_input_file(s)

# Interpolate phonon quantities onto a fine grid 
class interpolate_phonon(calculation):
    
    # The executable that carries out this calculation
    def exe(self):
        return "matdyn.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return self.in_params["ph_interp_prefix"]

    # Images are allowed for phonon calculations
    def image_parallelization_allowed(self):
        return True

    # Parse calculation output
    def parse_output(self, outf):
        return None

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        s  = "&INPUT\n"
        s += "    dos = .true.,\n"
        s += self.in_params.to_input_line("flfrc")
        s += self.in_params.to_input_line("zasr",          name="asr")
        s += self.in_params.to_input_line("ph_interp_nq1", name="nk1")
        s += self.in_params.to_input_line("ph_interp_nq2", name="nk2")
        s += self.in_params.to_input_line("ph_interp_nq3", name="nk3")
        s += self.in_params.to_input_line("ndos")
        s += self.in_params.to_input_line("ph_interp_dos_file",   name="fldos")
        s += self.in_params.to_input_line("ph_interp_freq_file",  name="flfrq")
        s += self.in_params.to_input_line("ph_interp_modes_file", name="flvec")
        s += self.in_params.to_input_line("la2F")
        s += self.in_params.to_input_line("el_ph_nsigma")
        s += "/\n"

        return pad_input_file(s)


#########################
# END CALCULATION TYPES #
#########################


# Calculate Tc using the allen-dynes equation from
# a given a2f(frequency)
def tc_allen_dynes(frequencies, a2f, mu_stars=[0.1, 0.15],
    attenuation_freq=None, get_lambda=False):
    import numpy         as     np
    from   qet.constants import RY_TO_K

    tc_ad = {}
    for mu in mu_stars:
        tc_ad[mu] = 0.0

    if attenuation_freq is None:
        wa  = [[w,max(a,0)] for w,a in zip(frequencies, a2f) if w > 0]
        ws  = [w for w,a in wa]
    else:
        ata = lambda a, w : a * (1.0 - np.exp(-w/attenuation_freq))
        wa  = [[w,max(ata(a,w),0)] for w,a in zip(frequencies, a2f) if w > 0]
        ws  = [w for w,a in wa]

    # Use the allen-dynes equation to estimate Tc
    lam  = np.trapz([2*a/w for w, a in wa], x=ws)
    if lam < 10e-10: 
        if get_lambda: 
            return tc_ad, lam
        return tc_ad

    wlog = np.exp((2/lam)*np.trapz([np.log(w)*a/w for w, a in wa], x=ws))
    wrms = ((2/lam)*np.trapz([a*w for w, a in wa], x=ws))**0.5
   
    for mu in mu_stars:
        g1 = 2.46*(1+3.8*mu)
        g2 = 1.82*(1+6.3*mu)*(wrms/wlog)
        f1 = (1+(lam/g1)**(3.0/2.0))**(1.0/3.0)
        f2 = 1 + (wrms/wlog - 1) * (lam**2) / (lam**2 + g2**2)
        tc_ad[mu] = RY_TO_K*f1*f2*(wlog/1.20)*np.exp(-1.04*(1+lam)/(lam-mu-0.62*lam*mu))

    if get_lambda: 
        return tc_ad, lam
    return tc_ad
    
# For a given a2f.dos file, calculate Tc for
# each of the given mu* values, using the allen-dynes
# equation
def tc_from_a2f_allen_dynes(filename, mu_stars=[0.1, 0.15],
    attenuation_freq=None, get_lambda=False):

    # Parse the a2f file, ignore negative frequencies/elishberg function
    out = parser.a2f_dos_out(filename)

    return tc_allen_dynes(out["frequencies"], out["a2f"], 
        mu_stars=mu_stars, attenuation_freq=attenuation_freq, get_lambda=get_lambda)

def tc_from_gap_function(filename, plot=False): 
    from   scipy.optimize import curve_fit
    import numpy          as     np

    # Parse the superconducting gap vs temperature
    ts = []
    gs = []
    with open(filename) as f:
        for l in f:
            vals = [float(w) for w in l.split()]
            if len(vals) != 3: continue
            ts.append(vals[0])
            gs.append(vals[1]/vals[2])

    # Interpolate g(t) data to a finer grid
    #ts_interp = np.linspace(min(ts), max(ts), len(ts)*5)
    #gs_interp = [np.interp(t, ts, gs) for t in ts_interp]
    #ts = ts_interp
    #gs = gs_interp

    # Guess tc from just where gap first goes above
    # 5% of the maximum on decreasing Tc
    tc_guess = 0
    for i in range(len(ts)-1, -1, -1):
        tc_guess = ts[i]
        if gs[i] > 0.05 * max(gs): break

    def gap_model(t, tc, gmax):
        t = np.array([min(ti, tc) for ti in t])
        return gmax * np.tanh(1.74*np.sqrt(tc/t - 1))

    try:
        p0 = [tc_guess, max(gs)*2]
        par, cov = curve_fit(gap_model, ts, gs, p0)
    
        if plot:
            import matplotlib.pyplot as plt
            plt.plot(ts, gs, marker="+", label="data")
            ts_interp = np.linspace(min(ts), max(ts), 200)
            plt.plot(ts_interp, gap_model(ts_interp, *par), label="fit")
            plt.axvline(par[0], color="black", linestyle=":", label="Tc = "+str(par[0]))
            plt.legend()
            plt.show()

    except Warning as warn:
        log("Gap function fit failed with warning:", "tc.log")
        log(str(warn), "tc.log")
        par = [0, p0[1]]
        cov = float("inf")
    except Exception as err:
        log("Gap function fit failed with error:", "tc.log")
        log(str(err), "tc.log")
        raise err

    if np.isfinite(cov).all():
        tc  = par[0]
        err = cov[0][0]**0.5
    else:
        log("Tc covariance is infinite!", "tc.log")
        tc  = tc_guess
        err = float("inf")

    return tc, err

# Calculate TC by solving the eliashberg equations (requires elk)
def tc_from_a2f_file_eliashberg(filename, mu_stars=[0.1, 0.15], force=False, species_dir=None):
    out = parser.a2f_dos_out(filename)
    wa  = [[w,max(a,0)] for w,a in zip(out["frequencies"], out["a2f"]) if w > 0]
    ws  = [w for w,a in wa]
    a2f = [a for w,a in wa]
    dosn = filename.split("dos")[-1]
    basedir = "/".join(filename.split("/")[0:-1])
    tc_dir = basedir + "/tc_dos_" + dosn
    return tc_from_a2f_eliashberg(tc_dir, ws, a2f, mu_stars=mu_stars, force=force, species_dir=species_dir)

# Calculate TC by solving the eliashberg equations (requires elk)
def tc_from_a2f_eliashberg(tc_dir, ws, a2f, mu_stars=[0.1, 0.15], force=False, species_dir=None):
    import warnings
    from   subprocess import check_output

    # Get the smearing number, make a corresponding tc directory
    if os.path.isdir(tc_dir):
        if force: os.system("rm -r "+tc_dir)
        else:
            # Attempt to parse tc from file
            tc_eliashberg = {}

            # Loop over calculated mu* values
            for mu_dir in os.listdir(tc_dir):
                mu_dir = tc_dir + "/" + mu_dir
                if not os.path.isdir(mu_dir): continue
                mu = float(mu_dir.split("/")[-1].split("_")[1])

                # Parse the tc.out file for the eliashberg Tc
                if not os.path.isfile(mu_dir + "/tc.out"): continue
                with open(mu_dir + "/tc.out") as f:
                    for line in f:
                        if "liashberg" in line:
                            tc_eliashberg[mu] = float(line.split()[0])
                            break

            return tc_eliashberg

    os.system("mkdir "+tc_dir)

    # First, estimate Tc using the allen-dynes equation
    tc_ad = tc_allen_dynes(ws,a2f,mu_stars=mu_stars)
    
    if species_dir is None:
        # Find the elk species directory
        elk_base_dir = check_output(["which", "elk"])
        elk_base_dir = elk_base_dir.decode("utf-8").replace("src/elk\n","")
        species_dir  = elk_base_dir + "/species/"

    tc_eliashberg = {}
    for mu in mu_stars:

        # Create the directory for this value of mu
        mu_dir = tc_dir + "/mu_{0}".format(mu)
        os.system("mkdir "+mu_dir)

        # Create elk a2F input file
        with open(mu_dir+"/ALPHA2F.OUT", "w") as f:
            for i in range(len(ws)):
                # Convert Ry to Ha
                f.write("{0} {1}\n".format(ws[i]*0.5,a2f[i]))

        # Create elk input file
        # A lot of the parameters are not actually needed
        # but are neccassary to allow elk to start up
        with open(mu_dir+"/elk.in", "w") as f:
            f.write("tasks\n260\n\nntemp\n20\nmustar\n{0}\n".format(mu))
            f.write("\nwplot\n{0} {1} {2}\n-0.5 0.5\n".format(len(ws), 1, 1))
            f.write("sppath\n'{0}'\n".format(species_dir))
            f.write("atoms\n1\n'La.in'\n1\n0 0 0 0 0 0\n")
            f.write("avec\n1 0 0\n0 1 0\n0 0 1")

        # Run elk
        os.system("cd {0}; elk > /dev/null".format(mu_dir))

        tc, err = tc_from_gap_function(mu_dir+"/ELIASHBERG_GAP_T.OUT")
        tc_eliashberg[mu] = tc

        # Save the result
        with open(mu_dir+"/tc.out", "w") as f:
            f.write("{0} +/- {1} K (Eliashberg)\n".format(tc, err))
            f.write("{0} K (Allen-Dynes)\n".format(tc_ad[mu]))

        # Tidy some of the larger files up
        os.system("rm "+mu_dir+"/ELIASHBERG_GAP_RA.OUT")
        os.system("rm "+mu_dir+"/ELIASHBERG_Z_RA.OUT")

    return tc_eliashberg

# List all files in the given folder or subdirectories of it
def listfiles(folder):
    for root, folders, files in os.walk(folder):
        for filename in folders + files:
            yield os.path.join(root, filename)

# Traverse all subdirectories, calculating Tc for every
# a2F.dos file we find
def tc_from_a2f_eliashberg_recursive(root_dir, species_dir=None):

    # Run over all a2F.dos files and calculate Tc
    for f in listfiles(root_dir):
        if not "a2F.dos" in f: continue
        log("Calculating Tc for "+f, "tc.log")
        try: 
            tc_res = tc_from_a2f_file_eliashberg(f, species_dir=species_dir)
            log("Success", "tc.log")
            for mu in tc_res:
                log("    Mu = {0}, Tc = {1}".format(mu, tc_res[mu]), "tc.log")
        except Exception as e: 
            log("Failed with excpetion:\n"+str(e), "tc.log")
            pass

# Returns true if the TC calculation in the given directory
# has completed or not
def tc_calculation_complete(dirname):

    elph_in = dirname+"/elph.in"
    if not os.path.isfile(elph_in):
        return False

    # Parse the number of sigma values
    # fromthe .in file
    n_sig = None
    with open(elph_in) as f:
        for line in f:
            if "el_ph_nsigma" in line:
                n_sig = int(line.split("=")[-1].replace(",",""))

    if n_sig is None:
        log("Could not parse el_ph_nsigma from "+elph_in)
        return False

    # If there are less than that many a2F files, the calculation
    # has not completed
    all_exist = True
    for i in range(1, n_sig+1):
        if not os.path.isfile(dirname + "/a2F.dos{0}".format(i)):
            all_exist = False
            break

    return all_exist

# Calculate the conventional superconducting critical temeprature
# for a given parameter set
def calculate_tc(parameters, 
    primary_only=False, 
    skip_elph=False, 
    phonons_only=False,
    phonons_from_elph=False,
    tidy_after=True):

    log("Tc calculation using parameters:")
    log(str(parameters))

    # Get the k-point grid sizes
    # needed to determine the most sensible
    # double-delta smearing parameter.
    kpq = parameters["tc_kpqs"]

    kpqs = {}
    for k in kpq:
        kpqs["kpq_{0}".format(k)] = k

    # Save working directory
    base_wd = os.getcwd()

    # Run an electron-phonon coulping calculation
    # for each of the kpoint grid sizes
    for dirname in kpqs:

        try:
            if tc_calculation_complete(dirname):
                log("Calculation " + base_wd + "/" + dirname + " complete, skipping...")
                continue

            # Go into a directory for this kpoint grid
            os.system("mkdir "+dirname)
            os.chdir(dirname)
            log("Switched to directory "+dirname)

            # Setup the kpoint grid for this directory
            parameters["kpts_per_qpt"] = kpqs[dirname]

            # relax the structure
            parameters["la2F"] = False # We don't need the a2F flag for the relaxation
            if parameters["variable_cell"]:
                res = relax(parameters).run()
            else:
                res = fixed_cell_relax(parameters).run()
            parameters["atoms"]   = res["relaxed atoms"]
            if "relaxed lattice" in res:
                parameters["lattice"] = res["relaxed lattice"]

            # Calculate the projected density of states/bandstructure
            try:
                proj_dos(parameters).run(required=False)
                bands(parameters).run(required=False)
                extract_bands(parameters).run(required=False)
            except Exception as ignored_ex:
                log("Encountered an expection when carrying out non-essential task.")
                log(str(ignored_ex))

            if not skip_elph:

                # We're gonna need the Eliashberg function from now on
                # (unless we're just doing the phonons)
                parameters["la2F"] = not phonons_only

                # Run the succesion of neccasary calculations
                scf(parameters).run()
                if phonons_only: phonon_grid(parameters).run()
                else: electron_phonon_grid(parameters).run()

                # We're recovering just the phonon information
                # from an electron-phonon calculation
                if phonons_from_elph:
                    parameters["la2F"] = False

                q2r(parameters).run()
                interpolate_phonon(parameters).run()

                # El-Ph complete, tidy up the huge files created
                if tidy_after:
                    tidy_tc_calculations()

        except Exception as e:

            # Return to the base directory
            os.chdir(base_wd)
            raise e
    
        # Go back to the base directory
        os.chdir(base_wd)

# Recursively searches for files no longer needed
# by Tc calculations from calculate_tc() and deletes them
def tidy_tc_calculations(base_dir=".", remove_incomplete=False, skip_dirs=[], older_than=-1):

    # Find all subdirectories with an elph.in file
    for elph_in in listfiles(base_dir):
        if not elph_in.endswith("elph.in"): continue
        tc_dir = os.path.dirname(elph_in)

        # Check if we should skip this directory
        skip = False
        for to_skip in skip_dirs:
            if to_skip in tc_dir:
                skip = True
                break
        if skip: continue

        # Check if this file is old enough to be removed (in days)
        if older_than > 0:
            age = time.time() - os.path.getmtime(tc_dir)
            age = age / (3600.0*24.0)
            if age < older_than:
                continue

        # Check calculations have completed
        if not remove_incomplete:
        
            # Parse the number of sigma values
            # fromthe .in file
            n_sig = None
            with open(elph_in) as f:
                for line in f:
                    if "el_ph_nsigma" in line:
                        n_sig = int(line.split("=")[-1].replace(",",""))

            if n_sig is None:
                print("Could not parse el_ph_nsigma from "+elph_in)
                continue

            # If there are less than that many a2F files, it's not safe to delete stuff
            all_exist = True
            for i in range(1, n_sig+1):
                if not os.path.isfile(tc_dir + "/a2F.dos{0}".format(i)):
                    all_exist = False
                    break
            if not all_exist:
                print("Not removing unfinshed calculations in "+tc_dir)
                continue

        # Find big files
        files_to_remove = [
            tc_dir+"/dyna2F",
            tc_dir+"/tc_dos_*/mu_*/ELIASHBERG_IA.OUT"
        ]

        for f in os.listdir(tc_dir):
            if f.startswith("pwscf.wfc"):
                files_to_remove.append(tc_dir+"/"+f)

        dirs_to_remove = [
                tc_dir+"/pwscf.save",
                tc_dir+"/_ph0",
                tc_dir+"/elph_dir"
        ]

        # Remove the big stuff
        for f in files_to_remove:
            if os.path.isfile(f):
                print("removing file "+f)
                os.system("rm "+f)

        # Remove big directories
        for d in dirs_to_remove:
            if os.path.isdir(d):
                print("removing directory "+d)
                os.system("rm -r "+d)
