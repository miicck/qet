from qet.params import parameters
import sys

p = parameters(sys.argv[1])
s = p.generate_commensurate_supercell([float(x) for x in sys.argv[2:]])
s.save(sys.argv[1].replace(".in", ".supercell.in"))
