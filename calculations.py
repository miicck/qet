import os
from qet.logs import log

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
        
        recover = os.path.isfile(outf)
        if recover:
            # Test to see if the calculation is complete
            with open(outf) as f:
                if "JOB DONE" in f.read():
                    msg = "Calculation \"{0}\" is complete, skipping..."
                    log(msg.format(outf))
                    return

        with open(inf, "w") as f:
            f.write(self.gen_input_file(recover=recover))

        cmd = "cd {0}; mpirun -np {1} {2} < {3} > {4}"
        cmd = cmd.format(
            path, self.in_params["cores_per_node"], 
            self.exe(), inf, outf)

        log("Running:")
        log(cmd)
        os.system(cmd)

class scf(calculation):

    # The executable that carries out this calculation
    def exe(self):
        return "pw.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "scf"

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        # Control
        s  = pw_control_input(self.in_params, 
            calculation="scf", recover=recover)

        # Geometry
        s += input_geometry(self.in_params)

        return pad_input_file(s)

class relax(calculation):

    # The executable that carries out this calculation
    def exe(self):
        return "pw.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "relax"

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        # Control
        s  = pw_control_input(self.in_params, 
            calculation="vc-relax", recover=recover)

        # Geometry
        s += input_geometry(self.in_params)

        return pad_input_file(s)

class phonon_grid(calculation):

    # The executable that carries out this calculation
    def exe(self):
        return "ph.x"

    # The default filename for calculations of this type
    def default_filename(self):
        return "phonons"

    # Generate the input file for this calculation
    def gen_input_file(self, recover=False):

        qpg = self.in_params["qpoint_grid"]

        s  = "Calculate phonons on a coarse grid\n"
        s += "&INPUTPH\n"
        s += self.in_params.to_input_line("outdir")
        s += "    ldisp=.true.,\n"
        s += "    nq1= {0},\n".format(qpg[0])
        s += "    nq2= {0},\n".format(qpg[1])
        s += "    nq3= {0},\n".format(qpg[1])
        s += "/\n"

        # I've found ph.x sometimes crashes if 
        # the input file doesn't end in a blank line
        return pad_input_file(s) + "\n"
    
