#######################################################
# This file generates the vendiagrams for all 3 and 2 #
#######################################################


import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal
import math
import os
import pandas
from collections import Counter
import matplotlib.patches as mpatches
from matplotlib_venn import venn3
from matplotlib_venn import venn2



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
	

	H = G.subgraph(genes)
	Gcc=sorted(nx.connected_component_subgraphs(H), key = len, reverse=True)
	G0 = Gcc[0]
	G1 = Gcc[1]

	
	F = nx.compose(G0, G1)
	nodes = F.nodes()
	pos=nx.spring_layout(F)

	print "NUMBER OF NODES", len(nodes)
	mdd_bp = []
	bp_schizo = []
	schizo_1 = []
	all_disorders = []
	mdd_schizo = []
	bp_1 = []
	mdd_1 = []

	for n in nodes:
		if ((n in mdd) and (n in bp) and (n in schizo)):
			all_disorders.append(n)
		elif (n in mdd and n in bp):
			mdd_bp.append(n)
		elif (n in mdd and n in schizo):
			mdd_schizo.append(n)
		elif (n in bp and n in schizo):
			bp_schizo.append(n)
		elif n in schizo:
			schizo_1.append(n)
		elif n in mdd:
			mdd_1.append(n)
		elif n in bp:
			bp_1.append(n)
		else:
			print "fail"


	b0 = mpatches.Patch(color='#663300', label='GD & BP & SCZ')
	b1 = mpatches.Patch(color='#800080', label='GD & BP')
	b2 = mpatches.Patch(color='#33cc33', label='GD & SCZ')
	b3 = mpatches.Patch(color='#ff8c00', label='BP & SCZ')
	b4 = mpatches.Patch(color='#ffff00', label='SCZ')
	b5 = mpatches.Patch(color='#0000ff', label='GP')
	b6 = mpatches.Patch(color='#ff0000', label='BP')

	
	l = [b0, b1, b2, b3, b4, b5, b6]
	fig, ax = plt.subplots()
	ax.legend(l, ['GD, BP, SCZ', 'Depression & Bipolar Disorder', 'Depression & Schizophrenia', 'Schizophrenia & Bipolar Disorder', 'Schizophrenia Alone',
		'Depression Alone', 'Bipolar Disorder Alone'])
	plt.legend(handles=l)

	nx.draw_networkx_nodes(F,pos=pos,nodelist=all_disorders,node_color='#663300',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos=pos,nodelist=mdd_bp,node_color='#800080',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos=pos,nodelist=mdd_schizo,node_color='#33cc33',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos=pos,nodelist=bp_schizo,node_color='#ff8c00',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos=pos,nodelist=schizo_1,node_color='#ffff00',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos=pos,nodelist=mdd_1,node_color='#0000ff',node_size=300,alpha=0.9)
	nx.draw_networkx_nodes(F,pos=pos,nodelist=bp_1,node_color='#ff0000',node_size=300,alpha=0.9)

	nx.draw_networkx_edges(F,pos=pos)
	#nx.draw_networkx_labels(F, pos=pos)
	plt.legend()

	plt.legend(handles=[b0, b1, b2, b3, b4, b5, b6], loc='upper left')
	ax.set_yticklabels([])
	ax.set_xticklabels([])
	plt.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off',
    top='off',
    bottom='off')
	#plt.axis('off')
	# print "all", len(mdd_bp) + len(bp_schizo) + len(all_disorders) + len(mdd_schizo)
	# print "mdd_bp", len(mdd_bp)
	# print "bp_schizo", len(bp_schizo)
	# print "all", len(all_disorders)
	# print "mdd", len(mdd_1)
	# print "bp", len(bp_1)
	# print "schizo", len(schizo_1)

	# venn3(subsets = (len(bp_1), len(mdd_1), len(mdd_bp), len(schizo_1), len(bp_schizo), len(mdd_schizo), len(all_disorders)),
	# 	set_labels = ('Bipolar Disorder', 'Depression', 'Schizophrenia'))
	plt.show()

####################################################################################################
	
	genes = mdd + bp
	H = G.subgraph(genes)
	nodes = H.nodes() #224

	mdd_bp = []
	bp_1 = []
	mdd_1 = []

	for n in nodes:
		if (n in mdd and n in bp):
			mdd_bp.append(n)
		elif n in mdd:
			mdd_1.append(n)
		elif n in bp:
			bp_1.append(n)

	print "mdd_bp", len(mdd_bp)
	print "mdd", len(mdd_1)
	print "bp", len(bp_1)

	venn2(subsets = (len(bp_1), len(mdd_1), len(mdd_bp)),
		set_labels = ('Bipolar Disorder', 'Depression'))
	plt.show()



if  __name__ =='__main__':main()
