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
        s += "{0}\n".format(sp)


    # Atomic positions card
    s += "\nATOMIC_POSITIONS (crystal)\n"
    for a in params["atoms"]:
        s += "{0} {1} {2} {3}\n".format(a[0], *a[1])

    # K-points card
    s += "\nK_POINTS automatic\n"
    s += "{0} {1} {2} 0 0 0\n".format(*params["kpoint_grid"])

    return s


class relax:
    
    # Create a relax calculation from the
    # given input parameters
    def __init__(self, in_params):

        # Save the input parameters
        self.in_params = in_params


    def gen_input_file(self, filename="relax.in"):

        # Generate the input file for this calculation
        # Control namelist
        s  = "&CONTROL\n"
        s += "    calculation = 'vc-relax',\n"
        s += self.in_params.to_input_line("outdir")
        s += self.in_params.to_input_line("pseudo_dir")
        s += self.in_params.to_input_line("forc_conv_thr")
        s += "/\n\n"

        # System namelist
        s += "&SYSTEM\n"
        s += self.in_params.to_input_line("ntyp")
        s += self.in_params.to_input_line("nat")
        s += self.in_params.to_input_line("ibrav")
        s += self.in_params.to_input_line("ecutwfc")
        s += self.in_params.to_input_line("ecutrho")
        s += self.in_params.to_input_line("occupations")
        s += self.in_params.to_input_line("smearing")
        s += self.in_params.to_input_line("deguass")
        s += "/\n\n"

        # Electrons namelist
        s += "&ELECTRONS\n"
        s += self.in_params.to_input_line("mixing_beta")
        s += self.in_params.to_input_line("conv_thr")
        s += "/\n\n"

        # Ions namelist
        s += "&IONS\n"
        s += self.in_params.to_input_line("ion_dynamics")
        s += "/\n\n"

        # Cell namelist
        s += "&CELL\n"
        s += self.in_params.to_input_line("cell_dynamics")
        s += self.in_params.to_input_line("press")
        s += self.in_params.to_input_line("press_conv_thr")
        s += "/\n\n"

        # Geometry
        s += input_geometry(self.in_params)

        return pad_input_file(s)
