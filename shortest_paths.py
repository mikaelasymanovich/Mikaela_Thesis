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
	paths = []

	for source in A.nodes():
		for target in B.nodes():
			path = nx.shortest_path(G, source=source, target=target)
			l = len(path)
			paths.append(l)

			if l < shortest_length:
				shortest_length = l
				shortest_path = path

	#print nx.average_clustering(A)
	#print nx.average_clustering(B)
	return shortest_length, shortest_path

def get_components(G, subgraph):
	
	Gcc=sorted(nx.connected_component_subgraphs(subgraph), key = len, reverse=True)
	G0 = Gcc[0]
	G1 = Gcc[1]

	show_components(G, G0, G1)
	return G0, G1

def check(node, disease):
	f = open(os.path.join('./Diseases', disease))

	first_line = 0

	for line in f:
		if (first_line == 1):
			x = line.split('\t')
			associated_gene = x[14]

			formatted_gene = associated_gene.replace('"', '').replace(',', '').replace(';', '').replace('-', '')
			mapped_genes = formatted_gene.split(' ')
			for gene in mapped_genes:
				if gene == node:
					print "true"
					p_value = Decimal(x[27])
					print node, " : ", p_value
					return

		first_line = 1
	print "node not found: ", node

def show_components(G, first, second):
	#nx.draw(G,  pos=nx.spring_layout(G), node_color='g')
	#nx.draw(first, pos=nx.spring_layout(first), with_labels=True)
	#plt.show()
	#nx.draw(second, pos=nx.spring_layout(second), with_labels=True)
	first_nodes = list(first.nodes)
	second_nodes = list(second.nodes)
	
	#nx.draw_networkx_nodes(first, pos=nx.spring_layout(first), nodelist=first_nodes, node_color='g')
	#nx.draw_networkx_nodes(second, pos=nx.spring_layout(second), nodelist=second_nodes, node_color='y')
	# also need to draw the egdes
	#nx.draw_networkx_edges(first, pos=nx.spring_layout(first))
	#nx.draw_networkx_edges(second, pos=nx.spring_layout(second))

	plt.show()


def main():

	reactome = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')

	G = nx.read_edgelist(reactome, delimiter='\t')

	reactome.close()
	file = 'gwas_bipolar_disorder.tsv'

	#For ONE file
	genes = graph.readGwasAssociationFile(file)
	H = G.subgraph(genes)
	first, second = get_components(G, H)
	print "Diameter 1: ", nx.diameter(first)
	print "Diameter 2: ", nx.diameter(second)
	length, path = shortestPath(G, first, second)
	# print "DONE"

	
	#print length, path
	#for n in path:
	#	check(n, file)

	#edges = G.edges(nbunch='AGRN')
	#print len(edges)
	#x = nx.algorithms.distance_measures.diameter(G)
	#print x
	#FOR ALL FILES
	# for file in os.listdir('./Diseases'):
	#  	print(file)
	#  	if (file == 'vogelstein.txt'):
	#  		fx = open(os.path.join('./Diseases', file))
	#  		genes = fx.read().splitlines()
	#  		H = G.subgraph(genes)
	#  		first, second = get_components(H)
	# 		shortestPath(G, first, second)

	#  	elif (file != '.DS_Store'):
	#  		genes = graph.readGwasAssociationFile(file)
	#  		H = G.subgraph(genes)
	#  		first, second = get_components(H)
	# 		shortestPath(G, first, second)

if  __name__ =='__main__':main()