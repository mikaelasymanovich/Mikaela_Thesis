### Program that produces subnetworks for each disorder and prints statistics about that network
### Results from this program seen in RESULTS_.txt
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal
import math
import os
import pandas
from collections import Counter
import random

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

# Read directly from Vegas File and create list of genes below p-value threshold
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
		

# Read directly from GWAS downloaded file and create a list of genes
# disease = a file name
def readGwasAssociationFile(disease):
	genes = []
	insignificants = []
	total_genes = []
	#f = open(os.path.join('./Diseases', 'ulcerative_colitis.tsv'))
	f = open(os.path.join('./Diseases', disease))
	max = 0;
	first_line = 0
	num_genes = 0;
	num_assoc = 0;
	for line in f:
		if (first_line == 1):
			num_assoc += 1;
			x = line.split('\t')
			p_value = Decimal(x[27])
			if p_value > max:
				max = p_value
			associated_gene = x[14]

			if (disease == 'gwas_schizophrenia.tsv' or disease == 'gwas_Alzheimer.tsv'):
				if (p_value < 0.0000001):
					formatted_gene = associated_gene.replace('"', '').replace(',', '').replace(';', '').replace('-', '')
					mapped_genes = formatted_gene.split(' ')
					for gene in mapped_genes:
						if gene:
							if gene not in genes: 
								genes.append(gene)
								total_genes.append(gene)
								num_genes += 1
			elif (disease == 'gwas_Dementia.tsv'):
				if (p_value < 0.00000005):
					formatted_gene = associated_gene.replace('"', '').replace(',', '').replace(';', '').replace('-', '')
					mapped_genes = formatted_gene.split(' ')
					for gene in mapped_genes:
						if gene: 
							if gene not in genes: 
								genes.append(gene)
								total_genes.append(gene)
								num_genes += 1
			else:
				if (p_value < 0.000005):
					formatted_gene = associated_gene.replace('"', '').replace(',', '').replace(';', '').replace('-', '')
					mapped_genes = formatted_gene.split(' ')
					for gene in mapped_genes:
						if gene: 
							if gene not in genes:
								genes.append(gene)
								total_genes.append(gene)
								num_genes += 1

		first_line = 1
	
	print max
	#print("Number of mapped genes: ", num_genes, "\nNumber of associations: ", num_assoc)
	return genes

def sampleDistribution(disease):
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
			if (disease == 'gwas_Dementia.tsv' or disease == 'gwas_schizophrenia.tsv'):
				if (p_value < 0.00000005):
					formatted_gene = associated_gene.replace('"', '').replace(',', '').replace(';', '').replace('-', '')
					mapped_genes = formatted_gene.split(' ')
					for gene in mapped_genes:
						if gene: 
							genes.append(gene)
							num_genes += 1
			else:
				if (p_value < 0.0005):
					formatted_gene = associated_gene.replace('"', '').replace(',', '').replace(';', '').replace('-', '')
					mapped_genes = formatted_gene.split(' ')
					for gene in mapped_genes:
						if gene: 
							genes.append(gene)
							num_genes += 1
		first_line = 1
	
	#print("Number of mapped genes: ", num_genes, "\nNumber of associations: ", num_assoc)
	
	a = genes
	letter_counts = Counter(a)
	#descending = letter_counts.most_common()
	df = pandas.DataFrame.from_dict(letter_counts, orient='index')
	#df.sort()
	df.plot(kind='bar')
	plt.show()

def permutation_test(k, lows, meds, highs, file, G, cc_number, G0, G1, density, cluster_coeff):
	N = 1000
	print "permutation test for: ", file
	nodes_all = G.nodes()

	num_gt_cc = 0
	num_gt_large_cc = 0
	second_largest_num = 0
	density_pr = 0
	cluster_coeff_pr = 0


	for i in range(0,N):
		l = 0; m = 0; h = 0;
		genes = []
		while l < lows and m < meds and h < highs:
			g = random.sample(nodes_all, 1)
			g = g[0] #pick a random node
			degree = G.degree(g)
			if degree > 150 and h < highs:
				genes.append(g)
				highs += 1
			elif degree > 50 and m < meds:
				genes.append(g)
				meds += 1
			elif degree > 0 and l < lows:
				genes.append(g)
				lows += 1
		
		if len(genes) != k:
			print("error!")

		H = G.subgraph(genes)
		nodes = nx.number_of_nodes(H)
		edges = nx.number_of_edges(H)

		if (nodes != 0):
			cc = nx.number_connected_components(H)
			density_ = (2.0*(edges))/(nodes*(nodes-1.0))
			cluster_coeff_ = nx.average_clustering(H)
			Gcc = sorted(nx.connected_component_subgraphs(H), key = len, reverse=True)
			G0_ = nx.number_of_nodes(Gcc[0])
			G1_ = nx.number_of_nodes(Gcc[1])


			if (cc > cc_number):
				num_gt_cc += 1
			if (G0_ > G0):
				num_gt_large_cc += 1
			if (G1_ > G1):
				second_largest_num += 1
			if (density_ > density):
				density_pr += 1
			if (cluster_coeff_ > cluster_coeff):
				cluster_coeff_pr += 1

	print("Number of Connected Components")
	print("P(X>x) = ", float(num_gt_cc)/float(N), "\n")

	print("Largest Connected Component")
	print("P(X>x) = ", float(num_gt_large_cc)/float(N), "\n")

	print("Probability that the second largest component")
	print("P(X>x) = ", float(second_largest_num)/float(N))

# Writes results for one disease to a file
# subgraph is the subgraph of the PPI graph of the genes for a specific disease
# f is the file to write to
# file is the filename aka the disease
def printResults(G, subgraph, f, file):
	edges = nx.number_of_edges(subgraph)
	nodes = nx.number_of_nodes(subgraph)
	genes = subgraph.nodes()
	components = sorted(nx.connected_components(subgraph), key = len, reverse=True)
	
	if (nodes != 0):
		density = (2.0*(edges))/(nodes*(nodes-1.0))
		cluster_coeff = nx.average_clustering(subgraph)
		cc_number = nx.number_connected_components(subgraph)
		Gcc=sorted(nx.connected_component_subgraphs(subgraph), key = len, reverse=True)
		G0=nx.number_of_nodes(Gcc[0])
		G1=nx.number_of_nodes(Gcc[1])

		graph_results = []
		for comp in Gcc:
			graph_results.append(nx.number_of_nodes(comp))

		#plt.hist(graph_results, G0)
		#name = "Number of Connected Components Distribution in " + file
		#plt.xlabel(name)
		#plt.show()
		#print graph_results
		low_ = 50
		medium_ = 150 

		lows = 0
		meds = 0
		highs = 0

		for g in genes:
			degree = G.degree(g)
			if degree > medium_:
				highs += 1
			elif degree > low_:
				meds += 1
			else:
				lows += 1

		#permutation_test(nodes, lows, meds, highs, file, G, cc_number, G0, G1, density, cluster_coeff)

    	cc_number_excuding_ones = 0
    	for component in Gcc:
    		if (nx.number_of_nodes(component) != 1):
    			cc_number_excuding_ones += 1;

    	#if (file == 'gwas_major_depressive_disorder.tsv'):
    	#	print("largest CC for MDD is: ", Gcc[0].nodes())
    	#if (file == 'gwas_depression.tsv'):
    	#	print("largest CC for Depression is: ", Gcc[0].nodes())

    	ratio = Decimal(cc_number)/nodes

	f.write(file + "\t" + str(nodes) + "\t" + str(edges) + "\t" + str(cc_number) + "\t" + str(G0) + "\t" + str(G1) + "\t" + str(density) + "\t" + str(cluster_coeff) + "\t" + str(cc_number_excuding_ones) + "\t" + str(ratio) + "\t" + str(lows) + "\t" + str(meds)+ "\t" + str(highs) + "\n")


def main():

	fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
	#gh = open(os.path.join('./PPI_Graphs', 'irefindex14.tsv'), 'rb')

	G = nx.read_edgelist(fh, delimiter='\t')
	#G_ref = nx.read_edgelist(gh, delimiter='\t')

	fh.close()
	#gh.close()

	f = open('confused.txt', 'w')
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
	 		#printResults(G, H, f, file)
	 	elif (file != '.DS_Store'):
	 		genes = readGwasAssociationFile(file)
	 		f_new = open(file, 'w')
	 		H = G.subgraph(genes)
	 		nodes = list(H.nodes())
	 		for x in nodes:
	 			f_new.write(x + "\n")
	 		#printResults(G, H, f, file)
	

	# nx.draw(H)
	# plt.show()


if  __name__ =='__main__':main()