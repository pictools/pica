import datetime
import itertools
import os
import subprocess

prefix = ""
if os.name != 'nt':
    prefix = "./"

# Modify parameters here
benchmark_path = prefix + "benchmark"
out_directory = datetime.datetime.now().strftime('benchmark_%Y-%m-%d_%H-%M-%S')
num_repetitions = 3
dimension = 3
ncellsx = 40
ncellsy = 40
ncellsz = 40
nppc = 50
temperature = 0.0
niterations = 100
layouts = ["SoA", "AoS"]
sorting_periods = [10, 20, 30]
ordered_combinations = ["ordered --sortingperiod " + str(period) for period in sorting_periods]
supercell_sizes = list(itertools.product([1, 2, 4], repeat=3))
supercell_combinations = ["supercells --preloading --ncellssupercellx " + str(size[0]) + " --ncellssupercelly " + str(size[1]) + " --ncellssupercellz " + str(size[2]) for size in supercell_sizes]
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
    field_solver_time = []
    particle_loop_time = []
    for repetition in range(0, num_repetitions):
        popen = subprocess.Popen(args, stdout=subprocess.PIPE, universal_newlines=True)
        popen.wait()
        output = str(popen.stdout.read())
        # cut off header
        parameters_output = output[output.find("Parameters"):output.find("Performance")]
        if repetition == 0:
            f.write(parameters_output)
        f.write("Run " + str(repetition) + ":\n")
        performance_output = output[output.find("   Field"):]
        f.write(performance_output + "\n")
        field_solver_time += [float(performance_output[len("   Field solver: "):performance_output.find(" sec")])]
        performance_output = performance_output[performance_output.find("   Particle loop"):]
        particle_loop_time += [float(performance_output[len("   Particle loop: "):performance_output.find(" sec")])]      

    total_time = []
    for idx in range(0, len(field_solver_time)):
        total_time.append(field_solver_time[idx] + particle_loop_time[idx])
    from operator import itemgetter
    best_run = min(enumerate(total_time), key=itemgetter(1))[0] 
    f.write("Best run:\n")
    f.write("   Field solver: " + str(field_solver_time[best_run]) + " sec.\n")
    f.write("   Particle loop: " + str(particle_loop_time[best_run]) + " sec.\n")
    f.write("   Overall: " + str(total_time[best_run]) + " sec.\n")