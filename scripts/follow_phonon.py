from qet.params import parameters
import sys

try:
    param_file = sys.argv[1]
    modes_file = sys.argv[2]
    qpt = [float(x) for x in sys.argv[3:6]]
    output_file = sys.argv[6]
except:
    print("Usage: python ~/../follow_phonon.py param.in matdyn.modes 0 0 0 param.out")
    quit()

p = parameters(param_file)
q = p.apply_phonon(qpt, 0, modes_file=modes_file, use_closest=True)
q.save(output_file)
