from qet.params       import parameters
from qet.logs         import log
from qet.calculations import scf, relax, phonon_grid

params = parameters(filename="test.in")

relax  = relax(params)
relax.run()

ph = phonon_grid(params)
ph.run()
