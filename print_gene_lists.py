## Write the gene lists to files

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
	
fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')

G = nx.read_edgelist(fh, delimiter='\t')

fh.close()

nodes_all = G.nodes()
k_alz = 97 #bipolar
k_aut = 112 #depression
k_dem = 191
k_bp = 153
k_dep = 190
k_scz = 157

g_alz = random.sample(nodes_all, k_alz)
g_aut = random.sample(nodes_all, k_aut)
g_dem = random.sample(nodes_all, k_dem)
g_bp = random.sample(nodes_all, k_bp)
g_dep = random.sample(nodes_all, k_dep)
g_scz = random.sample(nodes_all, k_scz)

f= open("g_alz", 'w')
for gene in g_alz:
	f.write(gene + "\n")

f= open("g_aut", 'w')
for gene in g_aut:
	f.write(gene + "\n")

f= open("g_dem", 'w')
for gene in g_dem:
	f.write(gene + "\n")

f= open("g_bp", 'w')
for gene in g_bp:
	f.write(gene + "\n")

f= open("g_dep", 'w')
for gene in g_dep:
	f.write(gene + "\n")

f= open("g_scz", 'w')
for gene in g_scz:
	f.write(gene + "\n")
