import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal
import math
import os

def readDiseaseFile():
	genes = []

	f = open(os.path.join('./Diseases', 'ulcerative_colitis.txt'))

	for line in f:
		x = line.split('\r')
		for entry in x:
			y = entry.split('\t')
			if (Decimal(y[1]) < 0.05):
				z = y[0].split()
				for m in z:
					if (m != '-'):
						l = m.replace('"', '').replace(',', '').replace(';', '')
						#l2 = l.replace(',', '')
						genes.append(l)
	return genes

def readVegas():
	f = open('vegas.csv')
	genes = []
	count = 0
	for line in f:
		entry = line.split(',')
		if (count > 0):
			if (Decimal(entry[7]) < 1):
				genes.append(entry[1])
		count+=1
	return genes
				

def readGwasAssociationFile(disease):
	genes = []

	#f = open(os.path.join('./Diseases', 'ulcerative_colitis.tsv'))
	f = open(os.path.join('./Diseases', disease))

	first_line = 0
	num_genes = 0;
	num_assoc = 0;
	for line in f:
		if (first_line == 1):
			num_assoc += 1;
			x = line.split('\t')
			p_value = Decimal(x[27])
			associated_gene = x[14]
			if (p_value < 0.000005):
				formatted_gene = associated_gene.replace('"', '').replace(',', '').replace(';', '').replace('-', '')
				mapped_genes = formatted_gene.split(' ')
				for gene in mapped_genes:
					if gene: 
						genes.append(gene)
						num_genes += 1
		first_line = 1
	
	#print("Number of mapped genes: ", num_genes, "\nNumber of associations: ", num_assoc)
	return genes

def printResults(subgraph):
	edges = nx.number_of_edges(subgraph)
	nodes = nx.number_of_nodes(subgraph)
	components = sorted(nx.connected_components(subgraph), key = len, reverse=True)
	print("Edges: ", edges)
	print("Nodes: ", nodes)
	if (nodes != 0):
		print("Ratio: ", Decimal(edges)/ (math.factorial(nodes)/(2*math.factorial((nodes-2)))))
		cluster_coeff = nx.average_clustering(subgraph)
		#print("Clustering Coefficient: ", cluster_coeff)


	#print("components: ", components)


def main():

	fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
	gh = open(os.path.join('./PPI_Graphs', 'irefindex14.tsv'), 'rb')

	G = nx.read_edgelist(fh, delimiter='\t')
	G_ref = nx.read_edgelist(gh, delimiter='\t')

	fh.close()
	gh.close()

	genes = readGwasAssociationFile('gwas_depression.tsv')
	H = G.subgraph(genes)
	#H_ref = G_ref.subgraph(genes)
	nx.draw(H, node_color='#FF4500')
	#nx.draw(H, node_color='g')
	plt.show()
	printResults(H)

	#genes = readVegas()
	#print all Values
	# print "BEGIN"
	# for file in os.listdir('./Diseases'):
	# 	print(file)
	# 	if (file != '.DS_Store'):
	# 		genes = readGwasAssociationFile(file)
	# 		H = G.subgraph(genes)
	# 		H_ref = G_ref.subgraph(genes)
	# 		print("REACTOME")
	# 		printResults(H)
	# 		print("IREFINDEX14")
	# 		printResults(H_ref)

	# nx.draw(H)
	# plt.show()


if  __name__ =='__main__':main()