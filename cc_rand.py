import networkx as nx
import matplotlib.pyplot as plt
from decimal import Decimal
import numpy as np
import os
import random
import pylab as P
import math
import statistics

num_ccs = []
largest_ccs = []
num_ccs_excluding_ones = []
seconds = []

num_gt_cc = 0
num_gt_large_cc = 0
num_gt_cc_excluding_ones = 0
second_largest_num = 0


N = 1000
k = 127
disease_cc_number = 16
disease_largest_cc = 111
second_largest = 2
disease_cc_number_excluding_ones = 2


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
		cc_number = nx.number_connected_components(H)
		num_ccs.append(cc_number) #add number of CC's to list

		Gcc =sorted(nx.connected_component_subgraphs(H), key = len, reverse=True)
		G0 =nx.number_of_nodes(Gcc[0])
		G1 = nx.number_of_nodes(Gcc[1])
		largest_ccs.append(G0) #add largest CC to list

		cc_number_excuding_ones = 0
		for component in Gcc:
			if (nx.number_of_nodes(component) != 1):
				cc_number_excuding_ones += 1
		num_ccs_excluding_ones.append(cc_number_excuding_ones)

		if (cc_number > disease_cc_number):
			num_gt_cc += 1
		if (G0 > disease_largest_cc):
			num_gt_large_cc += 1
		if (cc_number_excuding_ones > disease_cc_number_excluding_ones):
			num_gt_cc_excluding_ones += 1
		if (G1 >= second_largest):
			second_largest_num += 1
			seconds.append(G1)


print("Number of Connected Components \n")
print("P(X>x) = ", float(num_gt_cc)/float(N))
#print ("The Average number of CCs is:", sum(num_ccs)/len(num_ccs))
#print ("The Median number of CCs is:", statistics.median(num_ccs))

print("Largest Connected Component \n")
print("P(X>x) = ", float(num_gt_large_cc)/float(N))
#print ("The Average biggest CC is:", sum(largest_ccs)/len(largest_ccs))
#print ("The Median largest CC is:", statistics.median(largest_ccs))

print("Number of CC's Excluding Single Nodes \n")
print("P(X>x) = ", float(num_gt_cc_excluding_ones)/float(N))

print("Probability that the second largest component \n")
print("P(X>x) = ", float(second_largest_num)/float(N))


plt.hist(num_ccs)
plt.xlabel('Number of Connected Components')
plt.show()
plt.hist(largest_ccs)
plt.xlabel('Largest Connected Component')
plt.show()
plt.hist(num_ccs_excluding_ones)
plt.xlabel('Number of Connected Components Excluding Single Nodes')
plt.show()

#first create a single histogram

# mu, sigma = 200, 25
# x = mu + sigma*P.randn(10000)

# # the histogram of the data with histtype='step'
# n, bins, patches = P.hist(ratios, 10, normed=1, histtype='stepfilled')
# P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

# # add a line showing the expected distribution
# y = P.normpdf(bins, mu, sigma)
# l = P.plot(bins, y, 'k--', linewidth=1.5)

# # the histogram of the data with histtype='step'
# P.xlabel('Edge Density')
# P.show()



fh.close()
g.close()




# for x in range(10):
#   print random.randint(1,101)