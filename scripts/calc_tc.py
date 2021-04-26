#!/home/mjh261/bin/anaconda3/bin/python
from qet.params       import parameters
from qet.calculations import calculate_tc
from qet.logs         import log
import os
import sys

# Make, and move to the calculation directory
f = sys.argv[1]
dname = f.replace(".in","")
os.system("mkdir "+dname)
os.chdir(dname)

# Start the calculation
log("Starting TC calculation in "+dname)
p = parameters()
p.load("../"+f)
calculate_tc(p)
