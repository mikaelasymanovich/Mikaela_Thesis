import networkx as nx
import matplotlib.pyplot as plt
from decimal import Decimal
import os
import random

fh = open(os.path.join('./PPI_Graphs', 'irefindex14.tsv'), 'rb')

count = 17103

g = open('random_nodes.txt', 'rb')
genes = []
for line in g:
	x = line.split('\n')
	print(x)
	genes.append(x[0])

# G = nx.read_edgelist(fh, delimiter='\t')
# H = G.subgraph(genes)
# edges = nx.number_of_edges(H)
# nodes = nx.number_of_nodes(H)

# print("Nodes", nodes, "Edges: ", edges)


fh.close()




# for x in range(10):
#   print random.randint(1,101)