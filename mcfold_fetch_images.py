# -----------------------------------------------------------------------------------------
# Version: 2021-03-04
# Author: Christian Steinmetzger, Petzold group
#
# This script fetches secondary structure images as pdf or ps files from the mcfold website
# -----------------------------------------------------------------------------------------

dir_name = '/Users/chstei/Postdoc/E. coli ribosome/test'
file_name = 'test_pruned.out.txt'

# Location for saving the downloaded images. Make sure this directory exists before starting the script
# Default: '/Downloads'
target_name = dir_name + '/Downloads_test'

nt_seq = 'AGAAUUCCAGGUGUAGCGGUGAAAUGCGUAGAGAUCUGGAGGAAU'
# Download the first n secondary structure images. Don't go higher than a few thousand to avoid flooding the mcfold
# website with too many requests
# Default: 1000
download_folds = 1500
# File extension for figure
# Default: 'pdf', supports export to multiple formats at the same time with, e.g., ['pdf', 'ps']
# Allowed formats: pdf, ps, jpg, ct
ext_figure = ['pdf']
# URL for generating the secondary structure figure. This should not need to be changed
# Default: 'https://major.iric.ca/cgi-bin/2DRender/render.cgi'
url = 'https://major.iric.ca/cgi-bin/2DRender/render.cgi'

# -------
# Imports
# -------
import requests
from bs4 import BeautifulSoup

dotbracket_data = dir_name + '/' + file_name
with open(dotbracket_data, 'r') as dotbracket_file:
    raw_list = [[index,     # List with 1-indexed number [0] to match mcfold online interface convention,
                 item[0],   # dot-bracket structure [1] and corresponding energy [2]
                 item[1].rstrip()]
                for index, item in enumerate((line.split(' ')
                                              for line in dotbracket_file), start=1)]
cut_list = raw_list[:download_folds]

for current_structure in cut_list:
    query = '?structure=>' + file_name.split('.')[0] + '_' + str(current_structure[0]) \
            + '|' + nt_seq \
            + '|' + current_structure[1] \
            + '%20' + current_structure[2] \
            + '&structno=' + str(current_structure[0])
    submit = requests.get(url + query)
    result = BeautifulSoup(submit.text, features='html.parser')
    for ext in ext_figure:
        link = result.find('a', string=ext.upper())
        download = requests.get(link.attrs['href'])
        output_file = target_name \
                      + '/' + str(current_structure[0]) \
                      + '_' + current_structure[2] \
                      + '.' + ext
        try:
            with open(output_file, 'wb') as downloaded:
                downloaded.write(download.content)
        except FileNotFoundError:
            print('# ----------------------------------------------------' + '-' * len(target_name))
            print(f'# Please create the folder {target_name} before running this script')
            print('# ----------------------------------------------------' + '-' * len(target_name))
            break
    print('Downloaded ' + str(current_structure[0]) + '/' + str(cut_list[-1][0]) + ' structures')
