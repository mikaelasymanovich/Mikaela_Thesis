import networkx as nx
import matplotlib.pyplot as plt
from decimal import Decimal
import os
import random

fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
count = 17103
xy = 162380

rands = []
index = 0
genes = []
#Generate 150 Random numbers
for j in range(150):
	num = random.randint(1,17103)
	rands.append(num)
#generate 150 random genes
for line in fh:
	x = line.split('\t')
	if index in rands:
		p = random.randint(0,1)
		genes.append(x[p])
	index+=1

G = nx.read_edgelist(fh, delimiter='\t')
H = G.subgraph(genes)
edges = nx.number_of_edges(H)
print("Edges: ", edges)

fh.close()




# for x in range(10):
#   print random.randint(1,101)