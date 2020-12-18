from qet.params import parameters
import sys

param_file = sys.argv[1]
modes_file = sys.argv[2]
qpt = [float(x) for x in sys.argv[3:6]]
output_file = sys.argv[6]

p = parameters(param_file)
q = p.apply_phonon(qpt, 0, 1.0, modes_file=modes_file, use_closest=True)
q.save(output_file)
