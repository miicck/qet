import os
import subprocess
from   qet.logs   import log
from   qet.parser import scf_out, relax_out, phonon_grid_out, dos_out, proj_dos_out

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

            cmd = "cd {0}; mpirun -np {1} {2} < {3} > {4}"
            cmd = cmd.format(
                path, self.in_params["cores_per_node"], 
                self.exe(), inf, outf)

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
        return scf_out(outf)

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
        return dos_out(outf)

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
        return proj_dos_out(outf)

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
        return relax_out(outf)

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
        return phonon_grid_out(outf)

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        qpg = self.in_params["qpoint_grid"]

        s  = "Calculate phonons on a coarse grid\n"
        s += "&INPUTPH\n"
        s += self.in_params.to_input_line("outdir")
        if recover:
            s += "    recover = .true.,\n"
        s += "    ldisp = .true.,\n"
        s += "    nq1 = {0},\n".format(qpg[0])
        s += "    nq2 = {0},\n".format(qpg[1])
        s += "    nq3 = {0},\n".format(qpg[1])
        s += "/\n"

        # I've found ph.x sometimes crashes if 
        # the input file doesn't end in a blank line
        return pad_input_file(s) + "\n"
