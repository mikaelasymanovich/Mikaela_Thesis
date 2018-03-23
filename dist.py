import networkx as nx
import matplotlib.pyplot as plt
from decimal import Decimal
import numpy as np
import os
import random
import pylab as P
import math
import statistics

seconds = []
firsts = []

num_gt_cc = 0
num_gt_large_cc = 0
num_gt_cc_excluding_ones = 0
second_largest_num = 0


N = 10000
k = 190
disease_largest_cc = 35
second_largest = 14


fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
G = nx.read_edgelist(fh, delimiter='\t')
nodes_all = G.nodes()

for i in range(0,N):
	#os.system('gshuf -n 109 reactome_genes.txt  > random_nodes.txt')
	#g = open('random_nodes.txt', 'rb')
	genes = []

	#for line in g:
	#	x = line.split('\n')
	#	genes.append(x[0])

	#G = nx.read_edgelist(fh, delimiter='\t')
	genes = random.sample(nodes_all, k)

	H = G.subgraph(genes)
	nodes = nx.number_of_nodes(H)
	
	if (nodes != 0):
		Gcc =sorted(nx.connected_component_subgraphs(H), key = len, reverse=True)
		G0 =nx.number_of_nodes(Gcc[0])
		G1 = nx.number_of_nodes(Gcc[1])		

		if (G0 >= disease_largest_cc):
			firsts.append(G0)
		if (G1 >= second_largest):
			seconds.append(G1)

plt.hist(firsts)
plt.xlabel('Distribution of Largest Components')
plt.show()
plt.hist(seconds)
plt.xlabel('Distibution of Second Largest Components')
plt.show()

fh.close()
#g.close()




# for x in range(10):
#   print random.randint(1,101)