import networkx as nx
import matplotlib.pyplot as plt
from decimal import Decimal
import numpy as np
import os
import random
import pylab as P
import math


ratios = []
num_gt = 0

N = 1000
disease_ratio = 0.08

for i in range(0,N):
	fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
	os.system('gshuf -n 228 reactome_genes.txt  > random_nodes.txt')
	g = open('random_nodes.txt', 'rb')
	genes = []
	for line in g:
		x = line.split('\n')
		genes.append(x[0])

	G = nx.read_edgelist(fh, delimiter='\t')
	H = G.subgraph(genes)
	edges = nx.number_of_edges(H)
	nodes = nx.number_of_nodes(H)
	if (nodes != 0):
		ratio = float(edges)/float(math.factorial(nodes)/(2*math.factorial((nodes-2))))
		ratios.append(ratio)
		#Bonferroni correction on disease_ration?
		if (ratio > disease_ratio):
			num_gt+=1;
print("P(X>x) = ", float(num_gt)/float(N))

#plt.hist(ratios)
#plt.xlabel('Ratio of Edges / Nodes (for 150 Nodes)')

#plt.show()

#
# first create a single histogram
#
mu, sigma = 200, 25
x = mu + sigma*P.randn(10000)

# the histogram of the data with histtype='step'
n, bins, patches = P.hist(ratios, 10, normed=1, histtype='stepfilled')
P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

# add a line showing the expected distribution
y = P.normpdf(bins, mu, sigma)
l = P.plot(bins, y, 'k--', linewidth=1.5)

# the histogram of the data with histtype='step'
P.xlabel('Edge Density')
P.show()



fh.close()
g.close()




# for x in range(10):
#   print random.randint(1,101)