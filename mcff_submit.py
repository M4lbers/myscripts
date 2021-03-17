# -------------------------------------------------------------------------------
# Version: 2021-03-03
# Author: Christian Steinmetzger, Petzold group
#
# This script submits an RNA sequence to the mcff command line tool to generate a
# specified number or energy range of secondary structures for further analysis
# -------------------------------------------------------------------------------

# Batch jobs can be submitted by using equal-length lists of directories, file names and sequences
dir_name = ['/Users/chstei/Postdoc/E. coli ribosome/h29']
file_name = ['h29-ext8_complete.out.txt']
nt_seq = ['GCGGUGGGCUUCGGCUGAAUCGC']

# Specify one of these and set the other to None
energy_cutoff = 8
predict_folds = None

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
                report_string = f' structures within {energy_cutoff} kcal/mol'
            elif energy_cutoff is None:
                # Requesting a defined number of output structures adds four lines to the top of the output that need to
                # be removed before the file is written
                mcff_output_raw = subprocess.check_output(f'mcff -s {nt_seq} -ft {predict_folds} -v2', shell=True, text=True)
                mcff_output_cleaned = mcff_output_raw.splitlines(keepends=True)[4:]
                mcff_output = ''.join(mcff_output_cleaned)
                report_string = f'/{predict_folds} structures'
            else:
                print('# ------------------------------------------------------------------------------------')
                print('# Please request either an energy range or a number of suboptimal structures, not both')
                print('# ------------------------------------------------------------------------------------')
                break
            output_file.write(mcff_output)
        output_size = str(mcff_output.count('\n'))
        print('# ----------------------------' + '-' * (len(output_size) + len(report_string) + len(nt_seq)))
        print(f'# MC-Flashfold predicted {output_size}{report_string} for {nt_seq}')
        print('# ----------------------------' + '-' * (len(output_size) + len(report_string) + len(nt_seq)))
    except FileNotFoundError:
        print('# ----------------------------------------------------' + '-' * len(dir_name))
        print(f'# Please create the folder {dir_name} before running this script')
        print('# ----------------------------------------------------' + '-' * len(dir_name))
        break
