### Produces Graphs colored for degree for each disorder

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

def visualize(G, genes):
	sm = []
	sm_1 = []
	mid = []
	lg_1 = []
	lg_2 = []
	lg_3 = []

	H = G.subgraph(genes)
	Gcc=sorted(nx.connected_component_subgraphs(H), key = len, reverse=True)
	G0 = Gcc[0]
	G1 = Gcc[1]

	
	F = nx.compose(G0, G1)
	nodes = F.nodes()
	pos=nx.spring_layout(F)

	for n in nodes:
		# 10e-9
	    if G.degree(n) < 50:
	    	sm.append(n)
	    # 5 x 10e-9
	    elif G.degree(n) < 100:
	    	sm_1.append(n)
	    elif G.degree(n) < 150:
	        mid.append(n)
	    elif G.degree(n) < 200:
	    	lg_1.append(n)
	    elif G.degree(n) < 250:
			lg_2.append(n)
	    else:
	        lg_3.append(n)

	b0 = mpatches.Patch(color='#F8CDC1', label='d < 50')
	b1 = mpatches.Patch(color='#F9A58A', label='d < 100')
	b2 = mpatches.Patch(color='#F75927', label='d < 150')
	b3 = mpatches.Patch(color='#F73E04', label='d < 200')
	b4 = mpatches.Patch(color='#B02F07', label='d < 250')
	b5 = mpatches.Patch(color='#750504', label='d > 250')
	l = [b0, b1, b2, b3, b4, b5]
	fig, ax = plt.subplots()
	#ax.legend(l, ['d < 50', 'd < 100', 'd < 150', 'd < 200', "d < 250", "d > 250"])
	#plt.legend(handles=l)
	nx.draw_networkx_nodes(F,pos,nodelist=sm,node_color='#F8CDC1',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos,nodelist=sm_1,node_color='#F9A58A',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos,nodelist=mid,node_color='#F75927',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos,nodelist=lg_1,node_color='#F73E04',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos,nodelist=lg_2,node_color='#B02F07',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos,nodelist=lg_3,node_color='#750504',node_size=300,alpha=0.9)

	nx.draw_networkx_edges(F,pos=pos)
	#nx.draw_networkx_labels(F, pos=pos)
	#plt.legend(handles=[b0, b1, b2, b3, b4, b5], loc='upper left')
	#plt.colorbar(l)
	#plt.axis('off')
	ax.set_yticklabels([])
	ax.set_xticklabels([])
	ax.axis('off')
	plt.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off',
    top='off',
    bottom='off')
	plt.show()

def main():

	fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
	G = nx.read_edgelist(fh, delimiter='\t')
	fh.close()


	#file = 'gwas_depression.tsv'
	#genes, g_p = readGwasAssociationFile(file)
	#visualize(G, genes)
	for file in os.listdir('./Diseases'):
	 	print(file)
	 	if (file != '.DS_Store' and file != 'vogelstein.txt'):
	 		genes, g_p = readGwasAssociationFile(file)
	 		H = G.subgraph(genes)
	 		visualize(G, genes)

	# for file in os.listdir('./Diseases'):
	#  	print(file)
	#  	if (file == 'vogelstein.txt'):
	#  		fx = open(os.path.join('./Diseases', file))
	#  		genes = fx.read().splitlines()
	#  		visualize(G, genes)
	#  	elif (file != '.DS_Store'):
	#  		genes, g_p = readGwasAssociationFile(file)
	#  		visualize(G, genes)


if  __name__ =='__main__':main()
