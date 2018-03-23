import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal
import math
import os
import pandas
from collections import Counter
import matplotlib.patches as mpatches

def readGwasAssociationFile(disease):
	genes = []
	g_p = {}
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
							g_p[gene] = p_value
			elif (disease == 'gwas_Dementia.tsv'):
				if (p_value < 0.00000005):
					formatted_gene = associated_gene.replace('"', '').replace(',', '').replace(';', '').replace('-', '')
					mapped_genes = formatted_gene.split(' ')
					for gene in mapped_genes:
						if gene: 
							genes.append(gene)
							total_genes.append(gene)
							num_genes += 1
							g_p[gene] = p_value
			else:
				if (p_value < 0.000005):
					formatted_gene = associated_gene.replace('"', '').replace(',', '').replace(';', '').replace('-', '')
					mapped_genes = formatted_gene.split(' ')
					for gene in mapped_genes:
						if gene: 
							genes.append(gene)
							total_genes.append(gene)
							num_genes += 1
							g_p[gene] = p_value

		first_line = 1
	
	#print("Number of mapped genes: ", num_genes, "\nNumber of associations: ", num_assoc)
	return genes, g_p

def all_shortest_paths(G, A, B, genes):
    new_nodes = []
    for source in A.nodes():
        for target in B.nodes():
            # List of all lists of shortest paths from node A to node B
            gen = nx.all_shortest_paths(G, source=source, target=target)
            for path in gen:
                #this is a list
                for new_gene in path:
                    # if it is not already in the list of genes
                    if new_gene not in genes:
                        new_nodes.append(new_gene)
    return new_nodes

def get_components(G, subgraph):
    Gcc=sorted(nx.connected_component_subgraphs(subgraph), key = len, reverse=True)
    G0 = Gcc[0]
    G1 = Gcc[1]

    #show_components(G, G0, G1)
    return G0, G1

def main():

	fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
	G = nx.read_edgelist(fh, delimiter='\t')
	fh.close()


	file = 'gwas_depression.tsv'
	genes, g_p = readGwasAssociationFile(file)

	H = G.subgraph(genes)
	first, second = get_components(G, H)
	comp1 = first.nodes()
	comp2 = second.nodes()
	new_genes = all_shortest_paths(G, first, second, g_p)

	all_genes = list(comp1) + list(comp2) + new_genes
	Hx = G.subgraph(all_genes)

	pos=nx.spring_layout(Hx)

	nx.draw_networkx_nodes(Hx,pos,nodelist=comp1,node_color='#FF1A18',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(Hx,pos,nodelist=comp2,node_color='#1964F5',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(Hx,pos,nodelist=new_genes,node_color='#820EDB',node_size=300,alpha=0.9)

	nx.draw_networkx_edges(Hx,pos=pos)
	#nx.draw_networkx_labels(F, pos=pos)

	plt.show()


if  __name__ =='__main__':main()
