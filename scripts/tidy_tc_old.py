from qet.calculations import tidy_tc_calculations
import sys
if len(sys.argv) < 2:
    print("Please enter the age (in days) that a file must be to be removed.")
    quit()
tidy_tc_calculations(older_than=int(sys.argv[1]), remove_incomplete=True)
