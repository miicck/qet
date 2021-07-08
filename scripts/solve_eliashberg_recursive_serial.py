from qet.calculations import tc_from_a2f_eliashberg_recursive
import sys

species_dir = None
for arg in sys.argv:
    if arg.startswith("elk_species="):
        species_dir = arg.split("=")[-1]

tc_from_a2f_eliashberg_recursive(".", species_dir=species_dir)
