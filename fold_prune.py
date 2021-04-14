# --------------------------------------------------------------------------------
# Version: 2021-03-04
# Author: Christian Steinmetzger, Petzold group
#
# This script reads in trace files from mc-fold 2.32 or mcff and sorts out folding
# states that contain 1-bp stems or differ from one another just
# by an opened/closed base pair but have otherwise identical base pairing patterns
# --------------------------------------------------------------------------------

# Batch jobs can be submitted by using equal-length lists of directories and file names
dir_names = ['/Users/chstei/Postdoc/E. coli ribosome/test']
file_names = ['test_complete.out.txt']

# Only keep structures that include or not include specific base pairs or bulged nucleotides?
# Default: None, specify a base pair between nucleotides i and j as [(i, j)] or a bulge at nucleotide k as [(k, 0)].
# Multiple criteria can be given as, e.g., [(i, j), (k, 0)]
filter_bp = [(14, 30), (21, 24)]
# Default: 'include' to keep only structures that contain all of the requested base pair patterns. Set to 'exclude' to
# keep only structures that contain none of the requested base pair patterns
filter_mode = 'include'

init = parse = make_table = bp_filter = single_bp = breathing = discard = write = 0

# -------
# Imports
# -------
import forgi.graph.bulge_graph as fgb
import time


# -------
# Classes
# -------
class TraceFile:
    def __init__(self, current_dir_name, current_file_name):
        self.dir_name = current_dir_name
        self.file_name = current_file_name
        self.raw_list = []
        self.filter_list = []
        self.bridging_list = []
        self.breathing_list = []
        self.discard_list = []

    def parse_input(self):
        trace_data = self.dir_name + '/' + self.file_name
        with open(trace_data, 'r') as trace_file:
            self.raw_list = [[index,    # List with 1-indexed number [0] to match mcfold online interface convention,
                              item[0],  # dot-bracket structure [1] and corresponding energy [2]
                              item[1]]
                             for index, item in enumerate((line.split(' ')
                                                           for line in trace_file), start=1)]

    def make_bp_table(self):
        for current_structure in self.raw_list:
            bg = fgb.BulgeGraph.from_dotbracket(current_structure[1])
            current_structure.append(bg.to_pair_table()[1:])                # Add bp table [3] and boolean flag for
            current_structure.append(bg.length_one_stem_basepairs() != [])  # single-bp-bridges [4] to the list

    def filter_specific_bp(self):
        for current_structure in self.raw_list:
            print('***************')
            print('Index: ' + str(current_structure[0]))
            print(current_structure[1])
            for entry in filter_bp:
                if (filter_mode == 'include' and current_structure[3][entry[0] - 1] != entry[1]) or (filter_mode == 'exclude' and current_structure[3][entry[0] - 1] == entry[1]):
                    self.filter_list.append(current_structure[0])
                    print('Requested base pair not present')
                    break

    def remove_single_bp_bridges(self):
        for current_structure in self.raw_list:
            print('***************')
            print('Index: ' + str(current_structure[0]))
            print(current_structure[1])
            if current_structure[0] not in self.filter_list and current_structure[4]:        # The boolean flag is used to sort out any structure with a single-bp-stem
                self.bridging_list.append(current_structure[0])
                print('1-bp bridge')

    def remove_bp_breathing(self):
        for current_structure in self.raw_list:
            print('***************')
            print('Index: ' + str(current_structure[0]))
            if current_structure[0] in self.filter_list or current_structure[0] in self.bridging_list or current_structure[0] in self.breathing_list:
                print('Skipped')
            else:
                print(current_structure[1])
                for next_structure in self.raw_list[current_structure[0]:]:  # The index stored in line[0] for each line of raw_list starts at 1, so starting from slice [0] here automatically counts up by 1
                    print('---------------')
                    print('Next: ' + str(next_structure[0]))
                    if next_structure[0] in self.filter_list:
                        print('Skipped: Requested base pair not present')
                    elif next_structure[0] in self.bridging_list:
                        print('Skipped: 1-bp bridge')
                    elif next_structure[0] in self.breathing_list:
                        print('Skipped: Base pair breathing')
                    else:
                        print(next_structure[1])
                        for partner1, partner2 in zip(current_structure[3], next_structure[3]):
                            difference = abs(partner1 - partner2)
                            print(partner1, partner2, difference)
                            if difference != 0 and difference != max(partner1, partner2):
                                print('Different fold')
                                break
                        else:
                            self.breathing_list.append(next_structure[0])
                            print('Base pair breathing')

    def discard_structures(self):
        for line in sorted(self.filter_list + self.bridging_list + self.breathing_list, reverse=True):
            if line in self.filter_list or line in self.bridging_list:
                self.raw_list.pop(line - 1)
            else:
                self.discard_list.append(self.raw_list.pop(line - 1))

    def write_output(self):
        with open(self.dir_name+'/'+self.file_name.replace('complete', 'incl_14-30_21-24'), 'w') as pruned_file:
            for line in self.raw_list:
                pruned_file.write(f'{line[1]} {line[2]}\n')


if filter_mode not in ['include', 'exclude']:
    raise ValueError('Check that the filter mode has been set correctly')

start = time.time()
print(f'Started at {time.ctime(start)}\n')

init = time.time()
for dir_name, file_name in zip(dir_names, file_names):
    input_file = TraceFile(dir_name, file_name)

# -----------------------
# Parse mcfold trace file
# -----------------------
    input_file.parse_input()
    parse = time.time()
    input_file.make_bp_table()
    make_table = time.time()
# -----------------------
# Apply filter algorithms
# -----------------------
    if filter_bp is not None:
        input_file.filter_specific_bp()
        bp_filter = time.time()
    input_file.remove_single_bp_bridges()
    single_bp = time.time()
    input_file.remove_bp_breathing()
    breathing = time.time()
    input_file.discard_structures()
    discard = time.time()
# -----------------
# Write output file
# -----------------
    input_file.write_output()
    write = time.time()

finish = time.time()
print(f'\nFinished at {time.ctime(finish)}')

print(f'Init time: {round(init - start, 3)} s')
print(f'Parse time: {round(parse - init, 3)} s')
print(f'Table time: {round(make_table - parse, 3)} s')
print(f'Filter time: {round(bp_filter - make_table, 3)} s')
print(f'Single time: {round(single_bp - bp_filter, 3)} s')
print(f'Breathing time: {round(breathing - single_bp, 3)} s')
print(f'Discard time: {round(discard - breathing, 3)} s')
print(f'Write time: {round(write - discard, 3)} s')

print(f'Elapsed time: {round(finish - start, 3)} s')
