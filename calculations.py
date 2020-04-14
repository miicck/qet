import os
import subprocess
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
    return snew.strip()

# Get the geompetry part of an input file
def input_geometry(params):

    # Cell parameters card
    s  = "CELL_PARAMETERS (angstrom)\n"
    s += "{0} {1} {2}\n".format(*params["lattice"][0])
    s += "{0} {1} {2}\n".format(*params["lattice"][1])
    s += "{0} {1} {2}\n".format(*params["lattice"][2])

    # Atomic species card
    s += "\nATOMIC_SPECIES\n"
    for sp in params["species"]:
        s += "{0} {1} {2}\n".format(*sp)

    # Atomic positions card
    s += "\nATOMIC_POSITIONS (crystal)\n"
    for a in params["atoms"]:
        s += "{0} {1} {2} {3}\n".format(a[0], *a[1])

    # K-points card
    s += "\nK_POINTS automatic\n"
    s += "{0} {1} {2} 0 0 0\n".format(*params["kpoint_grid"])

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
    s += "/\n\n"

    # System namelist
    s += "&SYSTEM\n"
    s += params.to_input_line("ntyp")
    s += params.to_input_line("nat")
    s += params.to_input_line("ibrav")
    s += params.to_input_line("ecutwfc")
    s += params.to_input_line("ecutrho")
    s += params.to_input_line("occupations")
    s += params.to_input_line("smearing")
    s += params.to_input_line("degauss")
    s += params.to_input_line("la2F")
    s += "/\n\n"

    # Electrons namelist
    s += "&ELECTRONS\n"
    s += params.to_input_line("mixing_beta")
    s += params.to_input_line("conv_thr")
    s += "/\n\n"

    # Ions namelist
    s += "&IONS\n"
    s += params.to_input_line("ion_dynamics")
    s += "/\n\n"

    # Cell namelist
    s += "&CELL\n"
    s += params.to_input_line("cell_dynamics")
    s += params.to_input_line("press")
    s += params.to_input_line("press_conv_thr")
    s += "/\n\n"

    return s

# Returns True if the given file is a completed
# quantum espresso output
def is_complete(filename):
    if not os.path.isfile(filename):
        return False
    with open(filename) as f:
        return "JOB DONE" in f.read()

##################
#  CALCULATIONS  #
##################

class calculation:
    
    # Create an scf calculation from the
    # given input parameters
    def __init__(self, in_params):
        
        # Save the input parameters
        self.in_params = in_params

    # Get the string reperesentation of this
    # calculation <=> input file
    def __str__(self):
        s  = self.default_filename()+".in:\n"
        s += self.gen_input_file()
        return s
    
    # Run the calculation in the given directory
    # with the given name, will check if the
    # calculation exists and, if so, will attempt
    # to continue it if it is incomplete
    def run(self, filename=None, path="./"):

        if filename is None:
            filename = self.default_filename()
        
        # Remove trailing /'s from path
        path = path.strip()
        while path.endswith("/"): path = path[0:-1]

        inf  = path+"/"+filename+".in"
        outf = path+"/"+filename+".out"
        
        # Test to see if the calculation is complete
        if is_complete(outf):

            # Log that we are skipping this complete calculation
            msg = "Calculation \"{0}\" is complete, skipping..."
            log(msg.format(outf))

        else: # Calculation not complete

            # Create input file, run calculation
            recover  = os.path.isfile(outf)
            with open(inf, "w") as f:
                f.write(self.gen_input_file(recover=recover))

            # Use k-point parallelism
            qe_flags = "-nk {0}".format(self.in_params["cores_per_node"])

            cmd = "cd {0}; mpirun -np {1} {2} {3} -i {4} > {5}"
            cmd = cmd.format(
                path, self.in_params["cores_per_node"], 
                self.exe(), qe_flags, inf, outf)

            log("Running:\n"+cmd)
            try:
                # Run calculation, log stdout/stderr
                stdout = subprocess.check_output([cmd], stderr=subprocess.STDOUT, shell=True)
                log(stdout, filename="qet.stdout")
            except subprocess.CalledProcessError as e:
                # Log subprocess errror
                log(e)

            if not is_complete(outf):
                msg = "Calculation {0} did not complete, stopping!"
                msg = msg.format(outf)
                log(msg)
                raise RuntimeError(msg)

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
        s  = pw_control_input(self.in_params, 
            calculation="scf", recover=recover)

        # Geometry
        s += input_geometry(self.in_params)

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
            calculation="vc-relax", recover=recover)

        # Geometry
        s += input_geometry(self.in_params)

        return pad_input_file(s)

# Calculate phonons on a grid
class phonon_grid(calculation):

    # The executable that carries out this calculation
    def exe(self):
        return "ph.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "phonons"

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
        s += self.in_params.to_input_line("nq1")
        s += self.in_params.to_input_line("nq2")
        s += self.in_params.to_input_line("nq3")
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
        s += self.in_params.to_input_line("nq1")
        s += self.in_params.to_input_line("nq2")
        s += self.in_params.to_input_line("nq3")
        s += self.in_params.to_input_line("el_ph_sigma")
        s += self.in_params.to_input_line("el_ph_nsigma")
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

    # Parse calculation output
    def parse_output(self, outf):
        return None

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        s  = "&INPUT\n"
        s += self.in_params.to_input_line("fildyn")
        s += self.in_params.to_input_line("flfrc")
        s += self.in_params.to_input_line("zasr")
        s += self.in_params.to_input_line("la2F")
        s += self.in_params.to_input_line("el_ph_nsigma")
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

    # Parse calculation output
    def parse_output(self, outf):
        return None

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        s  = "&INPUT\n"
        s += "dos = .true.\n"
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
    
# For a given a2f.dos file, calculate Tc for
# each of the given mu* values, using the allen-dynes
# equation
def tc_from_a2f_allen_dynes(filename, mu_stars=[0.1, 0.15]):
    import numpy         as     np
    from   qet.constants import RY_TO_K

    # Parse the a2f file, ignore negative frequencies
    out = parser.a2f_dos_out(filename)
    wa  = [[w,a] for w,a in zip(out["frequencies"], out["a2f"]) if w > 0]
    ws  = [w for w,a in wa]

    # Use the allen-dynes equation to estimate Tc
    lam  = np.trapz([2*a/w for w, a in wa], x=ws)
    wlog = np.exp((2/lam)*np.trapz([np.log(w)*a/w for w, a in wa], x=ws))
    wrms = ((2/lam)*np.trapz([a*w for w, a in wa], x=ws))**0.5

    tc_ad = {}
    for mu in mu_stars:
        g1 = 2.46*(1+3.8*mu)
        g2 = 1.82*(1+6.3*mu)*(wrms/wlog)
        f1 = (1+(lam/g1)**(3.0/2.0))**(1.0/3.0)
        f2 = 1 + (wrms/wlog - 1) * (lam**2) / (lam**2 + g2**2)
        tc_ad[mu] = RY_TO_K*f1*f2*(wlog/1.20)*np.exp(-1.04*(1+lam)/(lam-mu-0.62*lam*mu))

    return tc_ad

# Calculate the conventional superconducting critical temeprature
# for a given parameter set
def calculate_tc(parameters):

    # Work out the two k-point grid sizes
    # needed to determine the most sensible
    # double-delta smearing parameter.
    kpq = parameters["kpts_per_qpt"]
    if kpq < 2: kpq = 2
    
    kpqs = {
        "aux_kpts"     : kpq-1, # Do the smaller kpoint grid first
        "primary_kpts" : kpq,   # Then the normal kpoint grid
    }

    # Save working directory
    base_wd = os.getcwd()

    # Run an electron-phonon coulping calculation
    # for each of the kpoint grid sizes
    for dirname in kpqs:

        try:
            # Go into a directory for this kpoint grid
            os.system("mkdir "+dirname)
            os.chdir(dirname)

            # Setup the kpoint grid for this directory
            parameters["kpts_per_qpt"] = kpqs[dirname]

            # relax the structure
            parameters["la2F"] = False # We don't need the a2F flag for the relaxation
            res = relax(parameters).run()
            parameters["atoms"]   = res["relaxed atoms"]
            parameters["lattice"] = res["relaxed lattice"]

            # We're gonna need the Eliashberg function from now on
            parameters["la2F"] = True

            # Run the succesion of neccasary calculations
            scf(parameters).run()
            electron_phonon_grid(parameters).run()
            q2r(parameters).run()
            interpolate_phonon(parameters).run()

        except Exception as e:

            # Return to the base directory
            os.chdir(base_wd)
            raise e
    
        # Go back to the base directory
        os.chdir(base_wd)
