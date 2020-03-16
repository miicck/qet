from   qet.params import parameters
import numpy      as     np

# A cubic phase of the hydride superconductor La H_10
# relaxed at 2500 KBar
lah10_fm3m = parameters()
lah10_fm3m["press"] = 2500

lat = []
lat.append([0.000000000,   2.422122552,   2.422122552])
lat.append([2.422122552,  -0.000000000,   2.422122552])
lat.append([2.422122552,   2.422122552,   0.000000000])

lah10_fm3m["lattice"] = np.array(lat)

atoms = []
atoms.append(["La",[-0.5000000000, 0.5000000000, 0.5000000000]])
atoms.append(["H" ,[ 1.6382108382, 0.1205963873, 0.1205963873]])
atoms.append(["H" ,[ 1.1205963873, 0.6382108382, 0.1205963873]])
atoms.append(["H" ,[ 1.1205963873, 0.1205963873, 0.6382108382]])
atoms.append(["H" ,[ 1.1205963873, 0.1205963873, 0.1205963873]])
atoms.append(["H" ,[ 0.3617891618, 0.8794036127,-0.1205963873]])
atoms.append(["H" ,[ 0.8794036127, 0.3617891618,-0.1205963873]])
atoms.append(["H" ,[ 0.8794036127, 0.8794036127,-0.1205963873]])
atoms.append(["H" ,[ 0.8794036127, 0.8794036127,-0.6382108382]])
atoms.append(["H" ,[ 1.2500000000, 0.2500000000, 0.2500000000]])
atoms.append(["H" ,[ 0.7500000000, 0.7500000000,-0.2500000000]])

lah10_fm3m["atoms"] = atoms
