#Probabilities for overlap
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

N = 1000
nodes_all = G.nodes()
k_1 = 153 #bipolar
k_2 = 190 #depression
k_3 = 157 #schizo
grp_all = grp_12 = grp_13 = grp_23 = grp_1 = grp_2 = grp_3 = alls = []
gt_all = gt_12 = gt_13 = gt_23 = gt_1 = gt_2 = gt_3 = pr_largest_component = 0

for i in range(0,N):
    group_1 = random.sample(nodes_all, k_1)
    group_2 = random.sample(nodes_all, k_2)
    #group_3 = random.sample(nodes_all, k_3)
    genes = group_1 + group_2 #+ group_3
    
    # Take subgraph of these combined subsets of genes
    H = G.subgraph(genes)
    # Get first two components
    Gcc=sorted(nx.connected_component_subgraphs(H), key = len, reverse=True)
    G0 = Gcc[0]
    G1 = Gcc[1]
    
    # Graph of only the first two largest components
    F = nx.compose(G0, G1)
    nodes = F.nodes()
    
    if nx.number_of_nodes(Gcc[0]) > 165:
        pr_largest_component += 1
    
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
    
    if len(grp_all) > 11:
        gt_all += 1
    if len(grp_12) > 55:
        gt_12 += 1
    if len(grp_13) > 2:
        gt_13 += 1
    if len(grp_23) > 1:
        gt_23 += 1
    if len(grp_1) > 14:
        gt_1 += 1
    if len(grp_2) > 33:
        gt_2 += 1
    if len(grp_3) > 53:
        gt_3 += 1
    alls.append(len(grp_all))
        

print("Largest Connected Component")
print("P(X>x) = ", float(pr_largest_component)/float(N), "\n")

print("# Nodes in all 3 samples")
print("P(X>x) = ", float(gt_all)/float(N), "\n")

print("# Nodes in group 1 and 2")
print("P(X>x) = ", float(gt_12)/float(N))

print("# Nodes in group 1 and 3")
print("P(X>x) = ", float(gt_13)/float(N))

print("# Nodes in group 2 and 3")
print("P(X>x) = ", float(gt_23)/float(N))

print("# Nodes in group 1")
print("P(X>x) = ", float(gt_1)/float(N))

print("# Nodes in group 2")
print("P(X>x) = ", float(gt_2)/float(N))

print("# Nodes in group 3")
print("P(X>x) = ", float(gt_3)/float(N))

plt.hist(alls)
plt.xlabel('Distribution of Number of Nodes in All 3 Sets')
plt.show()