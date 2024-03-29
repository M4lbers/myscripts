# --------------------------------------------------------------------------------
# Version: 2021-03-04
# Author: Christian Steinmetzger, Petzold group
#
# This script reads in trace files from mc-fold 2.32 or mcff and sorts out folding
# states that contain single-basepair stems or differ from one another just
# by an opened/closed base pair but have otherwise identical base pairing patterns
# --------------------------------------------------------------------------------

dir_names = ['/Users/chstei/Postdoc/E. coli ribosome/h23']
file_names = ['h23-top_complete.out.txt']

init = parse = make_table = single_bp = breathing = discard = write = 0

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

    def remove_single_bp_bridges(self):
        for current_structure in self.raw_list:
            print('***************')
            print('Index: ' + str(current_structure[0]))
            print(current_structure[1])
            if current_structure[4]:        # The boolean flag is used to sort out any structure with a single-bp-stem
                self.bridging_list.append(current_structure[0])
                print('Single-basepair bridge')

    def remove_bp_breathing(self):
        for current_structure in self.raw_list:
            print('***************')
            print('Index: ' + str(current_structure[0]))
            if current_structure[0] not in self.bridging_list and current_structure[0] not in self.breathing_list:
                print(current_structure[1])
                for next_structure in self.raw_list[current_structure[0]:]:  # The index stored in line[0] for each line of raw_list starts at 1, so starting from slice [0] here automatically counts up by 1
                    print('---------------')
                    print('Next: ' + str(next_structure[0]))
                    if next_structure[0] not in self.bridging_list and next_structure[0] not in self.breathing_list:
                        print(next_structure[1])
                        for partner1, partner2 in zip(current_structure[3], next_structure[3]):
                            difference = abs(partner1 - partner2)
                            print(partner1, partner2, difference)
                            if difference != 0 and difference != max(partner1, partner2):
                                print('Different fold')
                                break
                        else:
                            self.breathing_list.append(next_structure[0])
                            print('Basepair breathing')
                    elif next_structure[0] in self.bridging_list:
                        print('Skipped: Single-basepair bridge')
                    else:
                        print('Skipped: Basepair breathing')
            else:
                print('Skipped')

    def var_remove_bp_breathing(self):  # This implementation should do the same and might be slightly faster
        for current_structure in self.raw_list:
            print('***************')
            print('Index: ' + str(current_structure[0]))
            if current_structure[0] in self.bridging_list or current_structure[0] in self.breathing_list:
                print('Skipped')
            else:
                print(current_structure[1])
                for next_structure in self.raw_list[current_structure[0]:]:  # The index stored in line[0] for each line of raw_list starts at 1, so starting from slice [0] here automatically counts up by 1
                    print('---------------')
                    print('Next: ' + str(next_structure[0]))
                    if next_structure[0] in self.bridging_list:
                        print('Skipped: Single-basepair bridge')
                    elif next_structure[0] in self.breathing_list:
                        print('Skipped: Basepair breathing')
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
                            print('Basepair breathing')

    def discard_structures(self):
        for line in sorted(self.bridging_list + self.breathing_list, reverse=True):
            if line in self.bridging_list:
                self.raw_list.pop(line - 1)
            else:
                self.discard_list.append(self.raw_list.pop(line - 1))

    def write_output(self):
        with open(self.dir_name+'/'+self.file_name.replace('complete', 'pruned'), 'w') as pruned_file:
            for line in self.raw_list:
                pruned_file.write(f'{line[1]} {line[2]}\n')


start = time.time()
print(f'Started at {time.ctime(start)}\n')

for dir_name, file_name in zip(dir_names, file_names):
    input_file = TraceFile(dir_name, file_name)
    init = time.time()
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
    input_file.remove_single_bp_bridges()
    single_bp = time.time()
    input_file.var_remove_bp_breathing()
    breathing = time.time()
    input_file.discard_structures()
    discard = time.time()
# -----------------
# Write output file
# -----------------
    input_file.write_output()
    write = time.time()

# -----------------------
# Parse mcfold trace file
# -----------------------
# trace_data = dir_name+'/'+file_name
# with open(trace_data, 'r') as trace_file:
#     raw_list = [[index,
#                  item[0],
#                  item[1]]
#                 for index, item in enumerate((line.split(' ')
#                                               for line in trace_file), start=1)]
# for line in raw_list:
#     bg = fgb.BulgeGraph.from_dotbracket(line[1])
#     line.append(bg.to_pair_table()[1:])
#     line.append(bg.length_one_stem_basepairs() != [])

# -----------------------
# Apply filter algorithms
# -----------------------
# print('# ------------------------------------------')
# print('# Remove bulges with single-basepair bridges')
# print('# ------------------------------------------')
# bridging_list = []
# for current_structure in raw_list:
#     print('***************')
#     print('Index: ' + str(current_structure[0]))
#     print(current_structure[1])
#     if current_structure[4]:
#         bridging_list.append(current_structure[0])
#         print('Single-basepair bridge')

# print('# -----------------------------------------------')
# print('# Remove bulges resulting from basepair breathing')
# print('# -----------------------------------------------')
# breathing_list = []
# for current_structure in raw_list:
#     print('***************')
#     print('Index: ' + str(current_structure[0]))
#     if current_structure[0] not in bridging_list and current_structure[0] not in breathing_list:
#         print(current_structure[1])
#         for next_structure in raw_list[current_structure[0]:]:  # The index stored in line[0] for each line of raw_list starts at 1, so starting from slice [0] here automatically counts up by 1
#             print('---------------')
#             print('Next: ' + str(next_structure[0]))
#             if next_structure[0] not in bridging_list and next_structure[0] not in breathing_list:
#                 print(next_structure[1])
#                 # npdifference = np.absolute(current_structure[3] - next_structure[3])  # TODO the numpy code is actually slower than the for loop implementation below
#                 # npceiling = np.maximum(current_structure[3], next_structure[3])
#                 # npcomp = np.logical_or(npdifference == 0, npdifference == npceiling)
#                 # if not np.all(npcomp):
#                 #     print('Different fold')
#                 for partner1, partner2 in zip(current_structure[3], next_structure[3]):
#                     difference = abs(partner1 - partner2)
#                     print(partner1, partner2, difference)
#                     if difference != 0 and difference != max(partner1, partner2):
#                         print('Different fold')
#                         break
#                 else:
#                     breathing_list.append(next_structure[0])
#                     print('Basepair breathing')
#             elif next_structure[0] in bridging_list:
#                 print('Skipped: Single-basepair bridge')
#             else:
#                 print('Skipped: Basepair breathing')
#     else:
#         print('Skipped')

# discard_list = []
# for line in sorted(bridging_list+breathing_list, reverse=True):
#     if line in bridging_list:
#         raw_list.pop(line-1)
#     else:
#         discard_list.append(raw_list.pop(line-1))

# with open(dir_name+'/'+file_name.replace('complete', 'class'), 'w') as pruned_file:
#     for line in raw_list:
#         pruned_file.write('{0} {1}\n'.format(line[1], line[2]))

finish = time.time()
print(f'\nFinished at {time.ctime(finish)}')

print(f'Init time: {round(init - start, 3)} s')
print(f'Parse time: {round(parse - init, 3)} s')
print(f'Table time: {round(make_table - parse, 3)} s')
print(f'Single time: {round(single_bp - make_table, 3)} s')
print(f'Breathing time: {round(breathing - single_bp, 3)} s')
print(f'Discard time: {round(discard - breathing, 3)} s')
print(f'Write time: {round(write - discard, 3)} s')

print(f'Elapsed time: {round(finish - start, 3)} s')
