# -------------------------------------------------------------------------------
# Version: 2021-03-03
# Author: Christian Steinmetzger, Petzold group
#
# This script submits an RNA sequence to the mcff command line tool to generate a
# specified number or energy range of secondary structures for further analysis
# -------------------------------------------------------------------------------

# Batch jobs can be submitted by using equal-length lists of directories, file names and sequences
dir_name = ['/Users/chstei/Postdoc/E. coli ribosome/h23',
            '/Users/chstei/Postdoc/E. coli ribosome/h23',
            '/Users/chstei/Postdoc/E. coli ribosome/h23']
file_name = ['h23-top_complete.out.txt',
             'h23-top-GC_complete.out.txt',
             'h23-top-GC-GC_complete.out.txt']
nt_seq = ['GGUGUAGCGGUGAAAUGCGUAGAGACC',
          'GUGUAGCGGUGAAAUGCGUAGAGAC',
          'UGUAGCGGUGAAAUGCGUAGAGAC']

# Specify one of these and set the other to None
energy_cutoff = None
predict_folds = 1000

# -------
# Imports
# -------
import subprocess
import sys

if not len(dir_name) == len(file_name) == len(nt_seq):
    print('# --------------------------------------------------------------------')
    print('# Please provide individual directory and file names for each sequence')
    print('# --------------------------------------------------------------------')
    sys.exit()


for dir_name, file_name, nt_seq in zip(dir_name, file_name, nt_seq):
    try:
        with open(dir_name+'/'+file_name, 'w') as output_file:
            if predict_folds is None:
                mcff_output = subprocess.check_output(f'mcff -s {nt_seq} -t {energy_cutoff}', shell=True, text=True)
            elif energy_cutoff is None:
                # Requesting a defined number of output structures adds four lines to the top of the output that need to
                # be removed before the file is written
                mcff_output_raw = subprocess.check_output(f'mcff -s {nt_seq} -ft {predict_folds} -v2', shell=True, text=True)
                mcff_output_cleaned = mcff_output_raw.splitlines(keepends=True)[4:]
                mcff_output = ''.join(mcff_output_cleaned)
            else:
                print('# ------------------------------------------------------------------------------------')
                print('# Please request either an energy range or a number of suboptimal structures, not both')
                print('# ------------------------------------------------------------------------------------')
                break
            output_file.write(mcff_output)
    except FileNotFoundError:
        print('# ----------------------------------------------------' + '-' * len(dir_name))
        print(f'# Please create the folder {dir_name} before running this script')
        print('# ----------------------------------------------------' + '-' * len(dir_name))
        break
