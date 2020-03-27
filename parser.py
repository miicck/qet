from   qet.constants import EV_TO_RY, KBAR_AU3_TO_RY
from   qet.logs      import log
import os

# Base class for output file types
class output_file:
    
    # Initialized with the filename to parse
    def __init__(self, filename):

        # Initialize dictionary of parsed things
        self.dict = {}
        self.filename = filename
        self.parse(filename)

    # Convert this object to a string
    def __str__(self):
        
        max_l = max([len(k) for k in self.dict])
        fs = "    {0:" + str(max_l) + "} : {1}\n"
    
        string = "Output file values ({0}):\n".format(self.filename)
        for k in self.dict:
            
            # Custom formatting for atoms
            if k == "relaxed atoms":
                string += fs.format(k, len(self.dict[k]))
                afs = "        {0}  {1}  {2}  {3}\n"
                for a in self.dict[k]:
                    string += afs.format(a[0],*a[1])
                continue

            # Custom formatting for lattice
            if k == "relaxed lattice":
                string += fs.format(k, "(matrix)")
                lfs = "        {0}  {1}  {2}\n"
                for l in self.dict[k]:
                    string += lfs.format(*l)
                continue

            # Custom formatting for DOS things
            if k.startswith("DOS") or k.startswith("PDOS"):
                dfs = "length = {0}"
                dfs = dfs.format(len(self.dict[k]))
                string += fs.format(k, dfs)
                continue

            string += fs.format(k, self.dict[k])

        return string.strip()

    def __getitem__(self, key):

        if not key in self.dict:
            # Raise error if key not found
            msg = "Key {0} not found in output file {1}."
            msg = msg.format(key, self)
            raise ValueError(msg)

        return self.dict[key]

    def __setitem__(self, key, value):
        self.dict[key] = value

    # Parse things common to all scf type files
    def parse_scf_things(self, filename):   
        
        # Parse scf file line by line
        with open(filename) as f:
            for line in f:

                # Parse the total energy
                if "total energy" in line and "!" in line:
                    tot_e = float(line.split("=")[-1].replace("Ry",""))
                    self["total energy"] = tot_e
                    continue

                # Parse the fermi energy
                if "Fermi energy is" in line:
                    fermi_e = float(line.split("is")[-1].replace("ev",""))
                    self["fermi energy"] = EV_TO_RY * fermi_e
                    continue

                # Parse the volume
                if "unit-cell volume" in line:
                    vol = float(line.split("=")[-1].split()[0])
                    self["unit cell volume"] = vol
                    continue

                # Parse the number of atoms/cell
                if "number of atoms/cell" in line:
                    self["atoms/cell"] = int(line.split("=")[-1])
                    continue

                # Parse alat
                if "lattice parameter (alat)" in line:
                    alat = float(line.split("=")[-1].split("a.u")[0])
                    self["alat"] = alat
                    continue

                # Parse number of symm ops
                if "Sym. Ops." in line:
                    self["symmetry operations"] = int(line.split()[0])
                    continue

                # Parse CPU time/WALL time
                if "PWSCF" in line and "CPU" in line:
                    cpu_time          = line.split()[2]
                    wall_time         = line.split()[4]
                    self["cpu time"]  = cpu_time
                    self["wall time"] = wall_time
                    continue

                # Parse the pressure
                if "P=" in line:
                    self["pressure"] = float(line.split("P=")[-1])
                    continue

        # Compute derived quantities

        calc_ent =              "pressure"         in self.dict
        calc_ent = calc_ent and "unit cell volume" in self.dict
        calc_ent = calc_ent and "total energy"     in self.dict

        if calc_ent:
            
            # Calculate the enthalpy
            pv = self["unit cell volume"] * self["pressure"] * KBAR_AU3_TO_RY
            self["pv"]       = pv
            self["enthalpy"] = self["total energy"] + pv 


        if "atoms/cell" in self.dict:
            
            # Work out total energy per atom
            if "total energy" in self.dict:
                self["energy/atom"] = self["total energy"]/self["atoms/cell"]

            # Work out volume per atom
            if "unit cell volume" in self.dict:
                self["volume/atom"] = self["unit cell volume"]/self["atoms/cell"]

            # Work out enthalpy per atom
            if "enthalpy" in self.dict:
                self["enthalpy/atom"] = self["enthalpy"]/self["atoms/cell"]

class scf_out(output_file):

    # Parse an scf.out file
    def parse(self, filename):
        self.parse_scf_things(filename)

class relax_out(output_file):
    
    # Parse a relax.out file
    def parse(self, filename):
        self.parse_scf_things(filename)

        with open(filename) as f:
            lines = f.read().split("\n")
            for i, line in enumerate(lines):
                
                # Parsed the relaxed cell vectors
                if "CELL_PARAMETERS" in line:
                    lat = []
                    for j in range(i+1,i+4):
                        vec = [float(w) for w in lines[j].split()]
                        lat.append(vec)
                    self["relaxed lattice"] = lat

                # Parse the relaxed atom positions
                if "ATOMIC_POSITIONS" in line:
                    atoms = []
                    for j in range(i+1, i+1+self["atoms/cell"]):
                        name = lines[j].split()[0]
                        vec  = [float(w) for w in lines[j].split()[1:4]]
                        atoms.append([name, vec])
                    self["relaxed atoms"] = atoms

class phonon_grid_out(output_file):
    
    def parse(self, filename):
        pass

class dos_out(output_file):
    
    def parse(self, filename):

        # Get the directory of this calculation
        dos_dir  = os.path.dirname(filename)
        dos_file = dos_dir + "/pwscf.dos" 

        # Throw an error if the .dos file doesn't exist
        if not os.path.isfile(dos_file):
            msg = "Could not find the DOS file {0}"
            msg = msg.format(dos_file)
            log(msg)
            raise RuntimeError(msg)
        
        # Parse the .dos file
        energies  = []
        densities = []
        first_line = True
        with open(dos_file) as f:
            for line in f:

                # Skip the first line
                if first_line:
                    fermi_energy  = float(line.split()[-2])
                    fermi_energy *= EV_TO_RY 
                    first_line    = False
                    continue

                # Parse the energy and density
                e, d = [float(w) for w in line.split()[0:2]]
                energies.append(e * EV_TO_RY)
                densities.append(d)

        # Store the results
        self["DOS energies"] = energies
        self["DOS (energy)"] = densities
        self["fermi energy"] = fermi_energy

        # Linearly interpolate to get the DOS at
        # the fermi level
        for i in range(len(energies)):
            if energies[i] > fermi_energy:
                r = (energies[i] - fermi_energy)/(energies[i] - energies[i-1])
                self["DOS (E_F)"] = densities[i]*(1-r) + densities[i-1]*r
                break

class proj_dos_out(output_file):
    
    def parse(self, filename):

        # Get the directory of this calculation
        dos_dir  = os.path.dirname(filename)

        # Attempt to find either an scf, or relax calculation
        # to parse the fermi energy from
        if os.path.isfile(dos_dir+"/scf.out"):
            scf = scf_out(dos_dir+"/scf.out")
            self["fermi energy"] = scf["fermi energy"]

        elif os.path.isfile(dos_dir+"/relax.out"):
            rel =  relax_out(dos_dir+"/relax.out")
            self["fermi energy"] = rel["fermi energy"]

        # Parse all of the pdos files
        self["PDOS (energy)"]       = {}
        self["PDOS energies"]       = {}
        self["PDOS wfc names"]      = {}
        self["PDOS atom names"]     = {}
        self["PDOS (fermi energy)"] = {}

        for f in os.listdir(dos_dir):
            if not f.startswith("pwscf.pdos_atm"): continue

            # Parse projection details from filename
            atom_num  = int(f.split("#")[1].split("(")[0])
            atom_name = f.split("(")[1].split(")")[0]
            wfc_num   = int(f.split("#")[2].split("(")[0])
            wfc_name  = f.split("(")[2].split(")")[0]

            energies  = []
            ldos      = []

            # Parse the actual PDOS from the file
            with open(dos_dir + "/" + f) as pdos:
                first_line = True
                for line in pdos:
                    
                    # Skip the first line
                    if first_line:
                        first_line = False
                        continue

                    e,l = [float(w) for w in line.split()[0:2]]
                    energies.append(e * EV_TO_RY)
                    ldos.append(l)

            # Initialize dictionaries
            if not atom_num in self["PDOS energies"]:
                self["PDOS energies"][atom_num]       = {}
                self["PDOS (energy)"][atom_num]       = {}
                self["PDOS wfc names"][atom_num]      = {}
                self["PDOS atom names"][atom_num]     = {}
                self["PDOS (fermi energy)"][atom_num] = {}

            # Fill dictionaries
            self["PDOS energies"][atom_num][wfc_num]   = energies
            self["PDOS (energy)"][atom_num][wfc_num]   = ldos
            self["PDOS wfc names"][atom_num][wfc_num]  = wfc_name
            self["PDOS atom names"][atom_num][wfc_num] = atom_name

            # Calculate the DOS at the fermi level for
            # each atom and wavefunction
            if "fermi energy" in self.dict:
                fermi_e = self["fermi energy"]
                for i, e in enumerate(energies):
                    if e > fermi_e:
                        # Linearly interpolate to get the dos 
                        # at the fermi energy
                        k = "PDOS (fermi energy)"
                        f = (fermi_e - energies[i-1])/(e - energies[i-1])
                        d = ldos[i] * f + ldos[i-1] * (1-f)
                        self[k][atom_num][wfc_num] = d
                        break
