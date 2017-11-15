import datetime
import itertools
import os
import subprocess

# Modify parameters here
benchmark_path = "benchmark"
out_directory = datetime.datetime.now().strftime('benchmark_%Y-%m-%d_%H-%M-%S')
dimension = 3
size = 50
ppc = 100
temperature = 0.0
iterations = 100
representations = ["SoA", "AoS"]
storages = ["unordered", "ordered"]
# add other combinations here
combination_keys = ["-r", "-e"]
combination_values = list(itertools.product(representations, storages))
# End of parameters

# Enumerate all combinations of parameters and run
if not os.path.exists(out_directory):
    os.makedirs(out_directory)
args_base = (benchmark_path, "-d", str(dimension), "-g", str(size), "-p", str(ppc), "-t", str(temperature), "-i", str(iterations))
for i in range(0, len(combination_values)):
    file_name = ""
    args_combination = ()
    for j in range(0, len(combination_values[i])):
        args_combination += (combination_keys[j], combination_values[i][j])
        file_name += combination_values[i][j] + "_"
    args = args_base + args_combination
    popen = subprocess.Popen(args, stdout=subprocess.PIPE, universal_newlines=True)
    popen.wait()
    file_name = file_name[:-1] + ".txt"
    f = open(os.path.join(out_directory, file_name), "w")
    f.write(str(popen.stdout.read()))
