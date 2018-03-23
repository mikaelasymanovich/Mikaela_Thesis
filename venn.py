import networkx as nx
import matplotlib.pyplot as plt
from decimal import Decimal
import numpy as np
import os
import random
import pylab as P
import pandas as pd
import math
import sys
import matplotlib.patches as mpatches
from matplotlib_venn import venn3
import graph

fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
G = nx.read_edgelist(fh, delimiter='\t')
fh.close()

grp_all = grp_12 = grp_13 = grp_23 = grp_1 = grp_2 = grp_3 = alls = []
group_1 = graph.readGwasAssociationFile('gwas_depression.tsv')
group_2 = graph.readGwasAssociationFile('gwas_bipolar_disorder.tsv')
group_3 = graph.readGwasAssociationFile('gwas_schizophrenia.tsv')
all_genes = group_1 + group_2 + group_3


H = G.subgraph(all_genes)
nodes = H.nodes()

for n in nodes:
    if ((n in group_1) and (n in group_2) and (n in group_3)):
        grp_all.append(n)
    elif (n in group_1 and n in group_2):
        grp_12.append(n)
    elif (n in group_1 and n in group_3):
        grp_13.append(n)
    elif (n in group_2 and n in group_3):
        grp_23.append(n)
    elif n in group_3:
        grp_3.append(n)
    elif n in group_1:
        grp_1.append(n)
    elif n in group_2:
        grp_2.append(n)

venn3(subsets = (len(grp_1), len(grp_2), len(grp_12), len(grp_3), len(grp_13), len(grp_23), len(grp_all)),
set_labels = ('Depression', 'Bipolar Disorder', 'Schizophrenia'))
plt.show()
