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

gh = open('xx.txt', 'r')
genes = gh.read().splitlines()
H = G.subgraph(genes)
pos=nx.spring_layout(H)
nx.draw(H, pos=pos)
nx.draw_networkx_labels(H, pos=pos)

plt.show()



print genes
add_genes = []
ggg = ['IL7R', 'NKD1', 'IFIH1', 'HDAC7', 'HDAC9', 'TNXB', 'SMAD7', 'ITGA4', 'SMAD3', 'PLCG2', 'CD226', 'IL23R', 'IL27', 'PTPN2', 'TYK2', 'IKZF3', 'MAP3K8', u'STAT4', u'STAT3', u'IL18RAP', u'CIITA', u'RORC', u'PPBP', u'CXCR2', u'CXCR1', u'OSMR', u'JAK2', u'IL1R1', u'PTPRC', u'NFKB1', u'CDH3', u'IGF2', u'PDGFB', u'CD28', u'RPS6KA4', u'CARD9', u'PRKCB', u'CELSR3', u'CEBPA', u'CAMK2A', u'GPR65', u'CCL20', u'MAML2', u'FCGR2A', u'GNA12', u'IL2RA', u'RGS14', u'ADCY7', u'CXCL5', u'IRF5', 'ADCY3', 'IL10', 'ABI1', 'HNF4A', 'NOTCH1', 'NOTCH4', 'KEAP1', 'COL13A1', 'ITGAL', 'GRB7', 'PDE4A', 'ATF6B']
for file in os.listdir('./Diseases'):
	if (file != '.DS_Store') and (file != 'vogelstein.txt') and (file != 'gwas_breast_cancer.tsv') and (file != 'gwas_crohn.tsv') and (file != 'gwas_waist_hip.tsv') and (file != 'ulcerative_colitis.tsv') and (file != 'gwas_multiple_sclerosis.tsv'):
	 	g_ = graph.readGwasAssociationFile(file)
	 	H = G.subgraph(g_)
	 	nodes = list(H.nodes())
	 	for u in nodes: #genes in disease
	 		for v in genes: #genes in all 5 diseases
	 			if G.has_edge(u, v):
	 				add_genes.append(u)
	 				#print file, " : ", u
	 				if u in ggg:
	 					print file, " : ", u
final = genes + add_genes
H_1 = G.subgraph(final)
Gcc=sorted(nx.connected_component_subgraphs(H), key = len, reverse=True)
G0 = Gcc[0]
print list(G0.nodes())
pos=nx.spring_layout(G0)
nx.draw(G0, pos=pos)
nx.draw_networkx_labels(G0, pos=pos)

plt.show()



