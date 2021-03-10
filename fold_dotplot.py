# ---------------------------------------------------------------------------
# Version: 2021-03-04
# Author: Christian Steinmetzger, Petzold group
#
# This script reads in trace files from mc-fold 2.32 and generates dot plots
# ---------------------------------------------------------------------------

dir_name = '/Users/chstei/Postdoc/E. coli ribosome/h23'
file_name = 'h23-top_pruned.out.txt'

nt_seq = 'GGUGUAGCGGUGAAAUGCGUAGAGACC'
retain_folds = 1000
from_nt1 = 3
to_nt1 = 13
from_nt2 = 14
to_nt2 = 25

plot_limits = ['dot', 'bulge'] # TODO implement switching graphical limits on/off separately

# Show figure?
# Default: True
show_figure = True # TODO make this separate for all three plots
# Save figure?
# Default: True
save_figure = True
# Title for figure
# Default: file_name
title_figure = file_name
# File extension for figure
# Default: 'pdf', supports export to multiple formats at the same time with e.g. ['pdf', 'png']
ext_figure = ['pdf']

# -------
# Imports
# -------
import forgi.graph.bulge_graph as fgb
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


# ---------
# Functions
# ---------
def layer_setup(x_scale=1.0, y_scale=1.0):
    figure, layer = plt.subplots()
    layer.set_title(title_figure + '\n' + str(len(cut_list)) + ' out of ' + str(len(raw_list)) + ' structures',
                    loc='center')
    figure.set_size_inches(figure.get_size_inches()[0] * x_scale,
                           figure.get_size_inches()[0] * y_scale)  # Scale up from default size
    return figure, layer


def make_figure():
    if save_figure:
        for ext in ext_figure:
            plt.savefig(dir_name + '/' + file_name + '_dot.' + ext, bbox_inches='tight')
    if show_figure:
        plt.show()


#nt_nt = [nt for nt in nt_seq]   # Get individual nucleotides from sequence
nt_nt = [nt+str(i+1) for i, nt in enumerate(nt_seq)]   # Get individual nucleotides from sequence

dotbracket_data = dir_name+'/'+file_name
with open(dotbracket_data, 'r') as dotbracket_file:
    raw_list = [line.split(' ')[:2] for line in dotbracket_file]
cut_list = raw_list[:retain_folds]  # Remove the n folds with the highest energy

bp_freq = np.zeros([len(nt_seq), len(nt_seq)])
bp_hidden = np.ones([len(nt_seq), len(nt_seq)])         # Create alpha value matrix for bp limits
bp_background = np.zeros([len(nt_seq), len(nt_seq)])    # Create uniformly colored background for alpha channel
bulge_freq = []
energy_list = []
for line in range(len(cut_list)):
    bg = fgb.BulgeGraph.from_dotbracket(cut_list[line][0])
    bp = bg.to_pair_table()
    for nt, partner in enumerate(bp[1:]):
        if partner == 0:
            bulge_freq.append(nt)
    for nt, partner in enumerate(bp[1:(1+bp[0]) // 2]): # TODO Try to query the bg function for loops instead, or the partner > nt + 1 test
        if partner != 0:            # ignore unpaired nucleotides now, these could otherwise show up at the edges depending on the image drawing method
            if line == 0:           # Lowest-energy structure
                bp_freq[partner - 1, nt] = len(cut_list)
            if partner > nt + 1:    # ignore duplicates from going over the halfway point of the sequence
                bp_freq[nt, partner - 1] += 1
    energy_list.append(float(cut_list[line][1]))
    energy_array = np.array(energy_list)
bp_freq /= len(cut_list)            # Normalize to [0,1] base pairing probability range

for nt in range(len(nt_seq)):
    if nt in range(from_nt1 - 1) or nt in range(to_nt1, from_nt2 - 1) or nt in range(to_nt2, len(nt_seq)):
        bp_hidden[nt, :] = 0.25
        bp_hidden[:, nt] = 0.25
        bp_background[nt, :] = 0.25
        bp_background[:, nt] = 0.25

# --------------------------
# Set up layers for dot plot
# --------------------------
figure1, layer1 = layer_setup(x_scale=1.5, y_scale=1.5)

# -----------------
# Generate dot plot
# -----------------
#dotplot = layer.pcolormesh(np.ma.masked_where(bp_freq == 0, bp_freq), vmin=0.001, vmax=np.max(bp_freq))
background = layer1.imshow(bp_background, vmin=0, vmax=1, cmap='Greys')
dotplot = layer1.imshow(np.ma.masked_where(bp_freq == 0, bp_freq), vmin=0.001, vmax=np.max(bp_freq))
layer1.tick_params(which='major', top=True, labeltop=True, right=True, labelright=True)     # Dual x and y axes
layer1.tick_params(axis='x', labelrotation=45)
# layer.set(xlabel='5\'–3\'',
#           xticks=(np.arange(len(nt_seq)) + 0.5),    # Shift major ticks to the center of each cell
#           xticklabels=(nt_nt),
#           yticks=(np.arange(len(nt_seq)) + 0.5),
#           yticklabels=(nt_nt),
#           aspect='equal')
# layer.invert_yaxis()
layer1.set(xlabel='5\'–3\'',
           xticks=range(len(nt_seq)),  # Shift major ticks to the center of each cell
           xticklabels=nt_nt,
           #xlim=(from_nt1-1.5, to_nt1-0.5),  # Subtract 1 because of indexing and 0.5 because of centered ticks
           yticks=range(len(nt_seq)),
           yticklabels=nt_nt)#,
            # #ylim=(to_nt2-0.5, from_nt2-1.5))

# --------------
# Add grid lines
# --------------
# layer.set_xticks(np.arange(len(nt_seq)), minor=True)   # Add minor ticks for grid generation
# layer.set_yticks(np.arange(len(nt_seq)), minor=True)

layer1.set_xticks(np.arange(len(nt_seq))-0.5, minor=True)
#layer.set_xlim(from_nt1-1.5, to_nt1-0.5)
layer1.set_yticks(np.arange(len(nt_seq))-0.5, minor=True)
#layer.set_ylim(to_nt2-0.5, from_nt2-1.5)
layer1.tick_params(which='minor', length=0)              # Don't show the actual minor ticks
layer1.grid(which='minor', color='k', linewidth=0.5)
layer1.plot([np.arange(len(nt_seq))-0.5, np.arange(len(nt_seq))+0.5], [np.arange(len(nt_seq))-0.5, np.arange(len(nt_seq))+0.5], color='k', linewidth=0.5)    # Add diagonal

# ----------------
# Adjust color bar
# ----------------
divider = make_axes_locatable(layer1)
legend = divider.append_axes('right', size=figure1.get_size_inches()[0]*0.025, pad=figure1.get_size_inches()[0]*0.05)
figure1.colorbar(cm.ScalarMappable(norm=Normalize(vmin=0.001, vmax=1)), label='Paired frequency', ticks=np.linspace(0, 1, 5), cax=legend)

# ---------------
# Generate figure
# ---------------
make_figure()
    
# ---------------------------------
# Set up layers for bulge histogram
# ---------------------------------
figure2, layer2 = layer_setup(x_scale=1.5, y_scale=0.5)

# ------------------------
# Generate bulge histogram
# ------------------------
layer2.hist(bulge_freq, bins=range((len(nt_seq))))    # Explicit bin creating to ensure centering of the histogram bars
layer2.set_xlabel('5\'–3\'')
layer2.set_xticks(np.arange(len(nt_seq)) + 0.5)   # Shift major ticks to the center of each cell
layer2.set_xticklabels(nt_nt)
layer2.set_ylim(0, len(cut_list))
layer2.set_ylabel('Unpaired frequency')
layer2.set_yticks(np.linspace(0, len(cut_list), 5))
layer2.set_yticklabels(np.linspace(0, 1, 5))

# --------------
# Add grid lines
# --------------
layer2.set_xticks(np.arange(len(nt_seq)), minor=True)
layer2.grid(which='minor', color='k', linewidth=0.5)
layer2.tick_params(which='minor', length=1)
layer2.set_xlim(from_nt1-1, to_nt1)

# ---------------
# Generate figure
# ---------------
make_figure()

# -------------------------------------
# Set up layers for 2nd bulge histogram
# -------------------------------------
figure3, layer3 = layer_setup(x_scale=1.5, y_scale=0.5)

# ----------------------------
# Generate 2nd bulge histogram
# ----------------------------
layer3.hist(bulge_freq, bins=range((len(nt_seq))))    # Explicit bin creating to ensure centering of the histogram bars
layer3.set_xlabel('5\'–3\'')
layer3.set_xticks(np.arange(len(nt_seq)) + 0.5)   # Shift major ticks to the center of each cell
layer3.set_xticklabels(nt_nt)
layer3.set_ylim(0, len(cut_list))
layer3.set_ylabel('Unpaired frequency')
layer3.set_yticks(np.linspace(0, len(cut_list), 5))
layer3.set_yticklabels(np.linspace(0, 1, 5))

# --------------
# Add grid lines
# --------------
layer3.set_xticks(np.arange(len(nt_seq)), minor=True)
layer3.grid(which='minor', color='k', linewidth=0.5)
layer3.tick_params(which='minor', length=1)
layer3.set_xlim(from_nt2-1, to_nt2)

# ---------------
# Generate figure
# ---------------
make_figure()

# --------------------------------
# Set up layers for energy diagram
# --------------------------------
figure4, layer4 = layer_setup(x_scale=1.5, y_scale=0.5)

# -----------------------
# Generate energy diagram
# -----------------------
layer4.bar(range(1, np.size(energy_array) + 1), energy_array-energy_array.min())
layer4.set_xlim(0.5, np.size(energy_array) + 0.5)
layer4.set_xlabel('Structure')
layer4.set_xticks(range(1, np.size(energy_array) + 1))
layer4.set_xticklabels(range(1, np.size(energy_array) + 1))
layer4.set_ylabel(r'$\Delta$$G$ (kcal mol$^{-1}$)')     # matplotlib works best with raw strings, $$ takes care of the proper minus sign without using \u2212

# ---------------
# Generate figure
# ---------------
make_figure()
