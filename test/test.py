from qet.params       import parameters
from qet.logs         import log
from qet.calculations import scf, relax

params = parameters(filename="test.in")
print(params)
relax  = relax(params)
print(relax.gen_input_file())
#relax.run()
