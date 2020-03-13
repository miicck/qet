from qet.params       import read_parameters
from qet.logs         import log
from qet.calculations import relax

params = read_parameters("test.in")
relax  = relax(params)
print(relax.gen_input_file())
