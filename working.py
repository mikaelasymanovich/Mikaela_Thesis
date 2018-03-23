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


def main():

	fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
	G = nx.read_edgelist(fh, delimiter='\t')
	fh.close()

	# file = 'gwas_depression.tsv'
	# genes, g_p = readGwasAssociationFile(file)
	file = 'gwas_depression.tsv'
	mdd, g_p_mdd = readGwasAssociationFile(file)
	file = 'gwas_bipolar_disorder.tsv'
	bp, g_p_bp = readGwasAssociationFile(file)
	file = 'gwas_schizophrenia.tsv'
	schizo, g_p_schizo = readGwasAssociationFile(file)

	g_p = g_p_mdd.copy()   # start with x's keys and values
	g_p.update(g_p_bp)
	g_p.update(g_p_schizo)
	genes = mdd + bp + schizo

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
	    if g_p[n] < 0.000000001:
	    	sm.append(n)
	    # 5 x 10e-9
	    elif g_p[n] < 0.000000005:
	    	sm_1.append(n)
	    elif g_p[n] < 0.00000001:
	        mid.append(n)
	    elif g_p[n] < 0.00000005:
	    	lg_1.append(n)
	    elif g_p[n] < 0.0000001:
			lg_2.append(n)
	    else:
	        lg_3.append(n)

	b0 = mpatches.Patch(color='#f6fefd', label='< 10e-9')
	b1 = mpatches.Patch(color='#a6faf1', label='< 5 x 10e-9')
	b2 = mpatches.Patch(color='#6beaf9', label='< 10e-8')
	b3 = mpatches.Patch(color='#3ca3e8', label='< 5 x 10e-8')
	b4 = mpatches.Patch(color='#0E5FF6', label='< 10e-7')
	b5 = mpatches.Patch(color='#033579', label='> 10e-7')
	l = [b0, b1, b2, b3, b4, b5]
	fig, ax = plt.subplots()
	ax.legend(l, ['< 10e-9', '< 5 x 10e-9', '< 10e-8', '< 5 x 10e-8', "< 10e-7", "> 10e-7"])
	plt.legend(handles=l)
	nx.draw_networkx_nodes(F,pos,nodelist=sm,node_color='#f6fefd',node_size=300,alpha=1)
	nx.draw_networkx_nodes(F,pos,nodelist=sm_1,node_color='#a6faf1',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos,nodelist=mid,node_color='#6beaf9',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos,nodelist=lg_1,node_color='#3ca3e8',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos,nodelist=lg_2,node_color='#0E5FF6',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos,nodelist=lg_3,node_color='#033579',node_size=300,alpha=0.9)

	nx.draw_networkx_edges(F,pos=pos)
	#nx.draw_networkx_labels(F, pos=pos)
	#plt.legend(handles=[b0, b1, b2, b3, b4, b5], loc='upper left')
	#plt.colorbar(l)
	#plt.axis('off')
	plt.show()



if  __name__ =='__main__':main()
