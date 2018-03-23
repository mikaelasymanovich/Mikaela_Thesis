import networkx as nx
import matplotlib.pyplot as plt
from decimal import Decimal
import os
import random


file = open('reactome_genes.txt','w')
iref = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
genes = []

for line in iref:
	x = line.split('\t')
	if x[0] not in genes:
		genes.append(x[0])
		file.write(x[0])
		file.write('\n')
	if x[1] not in genes:
		genes.append(x[1])
		file.write(x[1])

file.close()
iref.close()
