###Create Matrix graph to show overlapping genes

import networkx as nx
import matplotlib.pyplot as plt
from decimal import Decimal
import numpy as np
import os
import random
import pandas as pd
import math
import sys
import pylab
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib_venn import venn3
import graph
from matplotlib.font_manager import FontProperties

fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
G = nx.read_edgelist(fh, delimiter='\t')
fh.close()

file = 'gwas_Alzheimer.tsv'
g_alz = graph.readGwasAssociationFile(file)
file = 'gwas_autism.tsv'
g_asd = graph.readGwasAssociationFile(file)
file = 'gwas_bipolar_disorder.tsv'
g_bp = graph.readGwasAssociationFile(file)
file = 'gwas_depression.tsv'
g_mdd = graph.readGwasAssociationFile(file)
file = 'gwas_dementia.tsv'
g_dem = graph.readGwasAssociationFile(file)
file = 'gwas_schizophrenia.tsv'
g_scz = graph.readGwasAssociationFile(file)

xx = open('confused.txt', 'r')
genes = xx.read().splitlines()
#print genes
#genes = g_alz + g_asd + g_bp + g_mdd + g_dem + g_scz
H = G.subgraph(genes)
nodes = list(H.nodes())
num_nodes = len(nodes)
#print num_nodes

matrix = np.zeros((num_nodes, 6)) #number of nodes = rows
for i in range (0, num_nodes):
    curr_node = nodes[i]
    if curr_node in g_alz:
        matrix[i][0] = 1
    else:
        matrix[i][0] = 0
    
    if curr_node in g_asd:
        matrix[i][5] = 1
    else:
        matrix[i][5] = 0
    
    if curr_node in g_bp:
        matrix[i][3] = 1
    else:
        matrix[i][3] = 0

    if curr_node in g_dem:
        matrix[i][1] = 1
    else:
        matrix[i][1] = 0
    if curr_node in g_mdd:
        matrix[i][2] = 1
    else:
        matrix[i][2] = 0
    if curr_node in g_scz:
        matrix[i][4] = 1
    else:
        matrix[i][4] = 0
#print matrix
disorders = ["AD", "BPSD", "GD",
              "BP", "SCZ", "ASD"]

#nodes = nodes
fig, ax = plt.subplots()
#im = ax.imshow(matrix)
#matrix = [[0, 0, 0, 0, 1, 0],
# [0, 0, 0, 0, 1, 0],
# [0, 0, 0, 0, 1, 0]]
#cmap = plt.cm.binary
ax.matshow(matrix, interpolation='nearest',aspect='equal',extent=[0,6,0,10])
ax.set_xticks(np.arange(len(disorders)))
#ax.set_yticks(np.arange(len(nodes)))
ax.set_xticklabels(disorders, rotation=45,ha='left', minor=False)
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off')

plt.ylabel("Genes", rotation=90)
#plt.legend()

#ax.set_yticklabels(nodes)
#plt.setp(ax.get_xticklabels(), ha='center', rotation=45)
plt.tight_layout()
plt.subplots_adjust(bottom=0.1, left=0, right=0.9, top=0.9,wspace=5, hspace=5)

#pylab.plot(range(10), label="Plot 1")
#pylab.plot(range(10, 0, -1), label="Plot 2")

p_patch = mpatches.Patch(color='purple', label='Gene not associated with disorder')
y_patch = mpatches.Patch(color='yellow', label='Gene associated with disorder')
#plt.legend(handles=[p_patch, y_patch])
#pylab.legend(loc=9, handles=[p_patch, y_patch], bbox_to_anchor=(0.5, -0.1))
plt.show()
fontP = FontProperties()
fontP.set_size('small')
art = []
lgd = pylab.legend(loc=9, prop=fontP, handles=[p_patch, y_patch], ncol=2)
art.append(lgd)
pylab.savefig(
    "fig.png", additional_artists=art,
    bbox_inches="tight")

plt.show()


