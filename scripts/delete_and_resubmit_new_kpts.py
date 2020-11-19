import os, sys

base_dir = sys.argv[1]
run_dir  = sys.argv[1] + "/" + sys.argv[2]
new_kpq  = int(sys.argv[3])
in_file  = run_dir + ".in"
sub_file = run_dir + "_submit"

if not os.path.isfile(in_file):
    raise Exception(in_file + " does not exist!")

if not os.path.isdir(run_dir):
    raise Exception(run_dir + " does not exist!")

if not os.path.isfile(sub_file):
    raise Exception(sub_file + " does not exist!")

# Replace input file with different kpts_per_qpt
with open(in_file) as f:
    lines = f.read().split("\n")

with open(in_file, "w") as f:
    for l in lines:
        if "kpts_per_qpt" in l:
            l = "kpts_per_qpt {0}".format(new_kpq)
        f.write(l+"\n")

# Remove run_dir and resubmit
os.system("rm -r "+run_dir)
os.system("sbatch "+run_dir+"_submit")
