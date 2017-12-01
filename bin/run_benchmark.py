import datetime
import itertools
import os
import subprocess

# Modify parameters here
benchmark_path = "benchmark"
out_directory = datetime.datetime.now().strftime('benchmark_%Y-%m-%d_%H-%M-%S')
num_repetitions = 3
dimension = 3
ncellsx = 10
ncellsy = 10
ncellsz = 10
nppc = 1
temperature = 0.0
niterations = 1
layouts = ["SoA", "AoS"]
sorting_periods = [10, 20, 30]
ordered_combinations = ["ordered --sortingperiod " + str(period) for period in sorting_periods]
supercell_sizes = list(itertools.product([1, 2, 4], repeat=3))
supercell_combinations = ["supercells --ncellssupercellx " + str(size[0]) + " --ncellssupercelly " + str(size[1]) + " --ncellssupercellz " + str(size[2]) for size in supercell_sizes]
orderings = ["unordered"] + ordered_combinations + supercell_combinations
combination_keys = ["--layout", "--ordering"]
combinations = list(itertools.product(layouts, orderings))
combination_values = []
for combination in combinations:
    new_value = (combination[0], )
    new_value += tuple(combination[1].split(" "))
    combination_values.append(new_value)
# End of parameters

# Enumerate all combinations of parameters and run
if not os.path.exists(out_directory):
    os.makedirs(out_directory)
args_base = (benchmark_path, )
args_base += ("--dimension", str(dimension))
args_base += ("--ncellsx", str(ncellsx))
args_base += ("--ncellsy", str(ncellsy))
args_base += ("--ncellsz", str(ncellsz))
args_base += ("--nppc", str(nppc))
args_base += ("--temperature", str(temperature))
args_base += ("--niterations", str(niterations))
for i in range(0, len(combination_values)):
    file_name = ""
    args_combination = ()
    for j in range(0, len(combination_values[i])):
        if j < len(combination_keys):
            args_combination += (combination_keys[j], combination_values[i][j])
        else:
            args_combination += (combination_values[i][j], )
        if "--" not in combination_values[i][j]:
            file_name += combination_values[i][j] + "_"
    args = args_base + args_combination
    file_name = file_name[:-1] + ".txt"
    f = open(os.path.join(out_directory, file_name), "w")
    for repeatition in range(0, num_repetitions):
        popen = subprocess.Popen(args, stdout=subprocess.PIPE, universal_newlines=True)
        popen.wait()
        f.write(str(popen.stdout.read()))
