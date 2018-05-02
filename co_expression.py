## Generate list of genes from gene IDs, print results
import networkx as nx
import numpy as np
from decimal import Decimal
import os
		

# Read directly from GWAS downloaded file and create a list of genes
# disease = a file name
def readGwasFile(disease):
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
			if x[15] and x[16]:
				upstream_downstream = x[15] + ' ' + x[16]
			if x[15] and not x[16]:
				upstream_downstream = x[15]
			if x[16] and not x[15]:
				upstream_downstream = x[16]
			else:
				upstream_downstream = ''
			#downstream_gene = x[16]
			if x[17]:
				middle_gene = x[17].replace(',', '')
			else:
				middle_gene = ''
			
			if upstream_downstream and middle_gene:
				_genes = upstream_downstream + ' ' + middle_gene
			if upstream_downstream and not middle_gene:
				_genes = upstream_downstream
			if middle_gene and not upstream_downstream:
				_genes = middle_gene
			else:
				_genes = ''

			mapped_genes = _genes.split(' ')

			if (disease == 'gwas_schizophrenia.tsv' or disease == 'gwas_Alzheimer.tsv'):
				if (p_value < 0.0000001):
					for gene in mapped_genes:
						if gene:
							if gene not in genes: 
								genes.append(gene)
								num_genes += 1
			elif (disease == 'gwas_Dementia.tsv'):
				if (p_value < 0.00000005):
					for gene in mapped_genes:
						if gene: 
							if gene not in genes: 
								genes.append(gene)
								num_genes += 1
			else:
				if (p_value < 0.000005):
					for gene in mapped_genes:
						if gene: 
							if gene not in genes:
								genes.append(gene)
								num_genes += 1

		first_line = 1
	
	#print("Number of mapped genes: ", num_genes, "\nNumber of associations: ", num_assoc)
	return genes

# Writes results for one disease to a file
# subgraph is the subgraph of the PPI graph of the genes for a specific disease
# f is the file to write to
# file is the filename aka the disease
def printResults(subgraph, f, file):
	edges = nx.number_of_edges(subgraph)
	nodes = nx.number_of_nodes(subgraph)
	components = sorted(nx.connected_components(subgraph), key = len, reverse=True)
	
	if (nodes != 0):
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

    	cc_number_excuding_ones = 0
    	for component in Gcc:
    		if (nx.number_of_nodes(component) != 1):
    			cc_number_excuding_ones += 1;

    	ratio = Decimal(cc_number)/nodes

	f.write(file + "\t" + str(nodes) + "\t" + str(edges) + "\t" + str(cc_number) + "\t" + str(G0) + "\t" + str(G1) + "\t" + str(cluster_coeff) + "\t" + str(cc_number_excuding_ones) + "\n")

def get_components(G, subgraph):
    Gcc=sorted(nx.connected_component_subgraphs(subgraph), key = len, reverse=True)
    G0 = Gcc[0]
    G1 = Gcc[1]

    #show_components(G, G0, G1)
    return G0, G1

def main():

	fh = open(os.path.join('./Tissue_Networks','front_03'), 'rb')

	G = nx.read_edgelist(fh, delimiter='\t')

	fh.close()

	f = open('_front_03.txt', 'w')
	f.write("Disease\tNodes\tEdges\tConnected_Components\tLargest_CC\tSecond_Largest_CC\tClustering_coeff\tCC's no single nodes" + "\n")
	print "BEGIN"
	for file in os.listdir('./Diseases'):
		print(file)
		if (file != '.DS_Store') and (file != 'vogelstein.txt'):
			genes = readGwasFile(file)
	 		H = G.subgraph(genes)
	 		printResults(H, f, file)



if  __name__ =='__main__':main()