import networkx as nx
import matplotlib.pyplot as plt
from decimal import Decimal
import os
import random

fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
list1 = []
count = 0

for line in fh:
	gene = line.split('\t')
	if gene[0] not in list1:
		count+=1
		list1.append(gene[0])
	if gene[1] not in list1:
		count+=1
		list1.append(gene[1])

print(count)


# for x in range(10):
#   print random.randint(1,101)