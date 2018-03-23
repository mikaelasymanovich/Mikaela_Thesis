import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal
import math
import os
import pandas
from collections import Counter

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

	first_line = 0
	num_genes = 0;
	num_assoc = 0;
	for line in f:
		if (first_line == 1):
			num_assoc += 1;
			x = line.split('\t')
			p_value = Decimal(x[27])
			associated_gene = x[14]

			if (disease == 'gwas_schizophrenia.tsv' or disease == 'gwas_Alzheimer.tsv'):
				if (p_value < 0.0000001):
					formatted_gene = associated_gene.replace('"', '').replace(',', '').replace(';', '').replace('-', '')
					mapped_genes = formatted_gene.split(' ')
					for gene in mapped_genes:
						if gene: 
							genes.append(gene)
							total_genes.append(gene)
							num_genes += 1
			elif (disease == 'gwas_Dementia.tsv'):
				if (p_value < 0.00000005):
					formatted_gene = associated_gene.replace('"', '').replace(',', '').replace(';', '').replace('-', '')
					mapped_genes = formatted_gene.split(' ')
					for gene in mapped_genes:
						if gene: 
							genes.append(gene)
							total_genes.append(gene)
							num_genes += 1
			else:
				if (p_value < 0.000005):
					formatted_gene = associated_gene.replace('"', '').replace(',', '').replace(';', '').replace('-', '')
					mapped_genes = formatted_gene.split(' ')
					for gene in mapped_genes:
						if gene: 
							genes.append(gene)
							total_genes.append(gene)
							num_genes += 1

		first_line = 1
	
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

# Writes results for one disease to a file
# subgraph is the subgraph of the PPI graph of the genes for a specific disease
# f is the file to write to
# file is the filename aka the disease
def printResults(subgraph, f, file):
	edges = nx.number_of_edges(subgraph)
	nodes = nx.number_of_nodes(subgraph)
	components = sorted(nx.connected_components(subgraph), key = len, reverse=True)
	
	if (nodes != 0):
		density = Decimal(edges)/ (math.factorial(nodes)/(2*math.factorial((nodes-2))))
		cluster_coeff = nx.average_clustering(subgraph)
		cc_number = nx.number_connected_components(subgraph)
		Gcc=sorted(nx.connected_component_subgraphs(subgraph), key = len, reverse=True)
		G0=nx.number_of_nodes(Gcc[0])
		G1=nx.number_of_nodes(Gcc[1])

		graph_results = []
		for comp in Gcc:
			graph_results.append(nx.number_of_nodes(comp))

		plt.hist(graph_results, G0)
		name = "Number of Connected Components Distribution in " + file
		plt.xlabel(name)
		plt.show()
		print graph_results


    	cc_number_excuding_ones = 0
    	for component in Gcc:
    		if (nx.number_of_nodes(component) != 1):
    			cc_number_excuding_ones += 1;

    	#if (file == 'gwas_major_depressive_disorder.tsv'):
    	#	print("largest CC for MDD is: ", Gcc[0].nodes())
    	#if (file == 'gwas_depression.tsv'):
    	#	print("largest CC for Depression is: ", Gcc[0].nodes())

    	ratio = Decimal(cc_number)/nodes

	f.write(file + "\t" + str(nodes) + "\t" + str(edges) + "\t" + str(cc_number) + "\t" + str(G0) + "\t" + str(G1) + "\t" + str(cluster_coeff) + "\t" + str(cc_number_excuding_ones) + "\t" + str(ratio) + "\n")


def main():

	fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
	gh = open(os.path.join('./PPI_Graphs', 'irefindex14.tsv'), 'rb')

	G = nx.read_edgelist(fh, delimiter='\t')
	G_ref = nx.read_edgelist(gh, delimiter='\t')

	fh.close()
	gh.close()

	#sampleDistribution('gwas_depression.tsv')
	#genes, insignificants, total_genes = readGwasAssociationFile('gwas_multiple_sclerosis.tsv')
	
	#H = G.subgraph(genes)
	#H_total = G.subgraph(total_genes)
	#H_insig = G.subgraph(insignificants)
	#H_ref = G_ref.subgraph(genes)
	#nx.draw_networkx_nodes(H_total, genes, node_color = 'r')
	#nx.draw_networkx_nodes(H_total, insignificants, node_color = 'b')
	#nx.draw(H, node_color='#FF4500')
	#plt.show()
	#printResults(H)

	f = open('results_irefindex.txt', 'w')
	f.write("Disease\tNodes\tEdges\tConnected_Components\tLargest_CC\tSecond_Largest_CC\tClustering_coeff\tCC's no single nodes\tNumber of CCs as Ratio" + "\n")
	#genes = readVegas()
	#print all Values
	print "BEGIN"
	for file in os.listdir('./Diseases'):
	 	print(file)
	 	if (file == 'vogelstein.txt'):
	 		fx = open(os.path.join('./Diseases', file))
	 		genes = fx.read().splitlines()
	 		#H = G.subgraph(genes)
	 		H_ref = G.subgraph(genes)
	 		print("REACTOME")
	 		printResults(H_ref, f, file)
	 	elif (file != '.DS_Store'):
	 		genes = readGwasAssociationFile(file)
	 		#H = G.subgraph(genes)
	 		H_ref = G.subgraph(genes)
	 		#print("REACTOME")
	 		#printResults(H, f, file)
	 		print("IREFINDEX14")
	 		printResults(H_ref, f, file)

	# nx.draw(H)
	# plt.show()


if  __name__ =='__main__':main()