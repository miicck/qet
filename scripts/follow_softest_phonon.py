from qet import parser
from qet.params import parameters
import sys

param_file = sys.argv[1]
modes_file = sys.argv[2]

modes = parser.phonon_interp_modes(modes_file)
qpts  = modes["q-points"]
freqs = modes["frequencies"]
modes_per_q = int(len(freqs)/len(qpts))
p = parameters(param_file)

if len(p["atoms"])*3 != modes_per_q:
    raise Exception("#modes in modes file != 3 x #atoms in params file!")

min_f = float("inf")
min_i = None
for i in range(len(freqs)):
    if freqs[i] < min_f:
        min_f = freqs[i]
        min_i = i

qi = min_i//modes_per_q
qpt = qpts[qi]
mode_index = min_i - modes_per_q*qi

print("Minimum frequency {0} cm^-1 found at q-point {1}, mode {2}".format(min_f, qpt, mode_index))
q = p.apply_phonon(qpt, mode_index, 1.0, modes_file=modes_file)
