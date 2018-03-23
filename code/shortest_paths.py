import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal
import math
import os
import sys
import graph

#compute the shortest paths from every node in component A to every node in component B
#Return the shortest path
def shortestPath(G, A, B):
	shortest_length = sys.maxint

	for source in A.nodes():
		for target in B.nodes():
			path = nx.shortest_path(G, source=source, target=target)
			if len(path) < shortest_length:
				shortest_length = len(path)
				shortest_path = path
	print length, path

def get_components(subgraph):
	
	Gcc=sorted(nx.connected_component_subgraphs(subgraph), key = len, reverse=True)
	G0 = Gcc[0]
	G1 = Gcc[1]

	return G0, G1


def main():

	reactome = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')

	G = nx.read_edgelist(reactome, delimiter='\t')

	reactome.close()

	print "BEGIN"
	genes = graph.readGwasAssociationFile('gwas_depression.tsv')
	H = G.subgraph(genes)
	first, second = get_components(H)
	shortestPath(G, first, second)
	print "DONE"


if  __name__ =='__main__':main()