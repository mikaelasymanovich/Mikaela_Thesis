## Check number of low, medium, high degree nodes
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal
import math
import os
import pandas
from collections import Counter
import random


genes = []
fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
lows = []; meds = []; highs = []

G = nx.read_edgelist(fh, delimiter='\t')

fh.close()

nodes = G.nodes()

for n in nodes:
	degree = G.degree(n)
	if degree > 150:
		
		elif degree > 50 and m < meds:
			genes.append(g)
			meds += 1
		elif degree > 0 and l < lows:
			genes.append(g)
			lows += 1
	f = open('RESULTS_most_recent!!.txt', 'w')
	f.write("Disease\tNodes\tEdges\tConnected_Components\tLargest_CC\tSecond_Largest_CC\tDensity\tClustering_coeff\tCC's no single nodes\tNumber of CCs as Ratio\tLow Degree Nodes\tMedium Degree Nodes\tHigh Degree Nodes" + "\n")

	#print all Values
	print "BEGIN"
	for file in os.listdir('./Diseases'):
	 	print(file)
	 	if (file == 'vogelstein.txt'):
	 		fx = open(os.path.join('./Diseases', file))
	 		genes = fx.read().splitlines()
	 		H = G.subgraph(genes)
	 		#H_ref = G.subgraph(genes)
	 		#print("REACTOME")
	 		printResults(G, H, f, file)
	 	elif (file != '.DS_Store'):
	 		genes = readGwasAssociationFile(file)
	 		f_new = open(file, 'w')
	 		H = G.subgraph(genes)
	 		nodes = list(H.nodes())
	 		for x in nodes:
	 			f_new.write(x + "\n")
	 		printResults(G, H, f, file)
	

	# nx.draw(H)
	# plt.show()


if  __name__ =='__main__':main()