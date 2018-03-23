import networkx as nx
import matplotlib.pyplot as plt
from decimal import Decimal
import numpy as np
import os
import random
import pylab as P
import math
import sys

seconds = [] 
firsts = []

paths = []
c1_diam = []

c2_diam= []
c1d = 0
c2d = 0
c1dd = 0
c1ddd = 0
cc1 = []
cc2 = []
paths = []

N = 1000
k = 190
disease_largest_cc = 35
second_largest = 14

fh = open(os.path.join('./PPI_Graphs', 'reactomefi2015.tsv'), 'rb')
G = nx.read_edgelist(fh, delimiter='\t')
nodes_all = G.nodes()

for i in range(0,N):
    genes = []
    genes = random.sample(nodes_all, k)
    shortest_length = sys.maxint
    
    H = G.subgraph(genes)
    nodes = nx.number_of_nodes(H)

    if (nodes != 0):
        #Get All Components
        Gcc = sorted(nx.connected_component_subgraphs(H), key = len, reverse=True)
        first = Gcc[0]
        second = Gcc[1]
        G0 = nx.number_of_nodes(first)
        G1 = nx.number_of_nodes(second)
        
        #Diameter
        c1_diam.append(nx.diameter(first))
        if nx.diameter(first) > 12:
            c1d += 1
        if nx.diameter(second) > 5:
            c2d += 1
        if nx.diameter(first) > 7:
            c1dd += 1
        if nx.diameter(first) > 5:
            c1ddd += 1
   
        c2_diam.append(nx.diameter(second))
        #Clustering Coefficients
        cc1.append(nx.average_clustering(first))
        cc2.append(nx.average_clustering(second))

        #Shortest Path between Components
        for source in first.nodes():
            for target in second.nodes():
                path = nx.shortest_path(G, source=source, target=target)
                l = len(path)
                
                if l < shortest_length:
                    shortest_length = l
        paths.append(l) #add the shortest path to a running list of 1000 shortest paths

print("Average Shortest Path over 1000 trials: ", sum(paths)/1000.0)
print("Average Diameter of Component 1: ", sum(c1_diam)/1000.0)
print("P(X>x) Component 1: 12 ", float(c1d)/1000.0)
print("P(X>x) Component 1: 7 ", float(c1dd)/1000.0)
print("P(X>x) Component 1: 5 ", float(c1ddd)/1000.0)

print("Average Diameter of Component 2: ", sum(c2_diam)/1000.0)
print("P(X>x) Component 2: ", float(c2d)/1000.0)

print("Average CC of Component 1: ", sum(cc1)/1000.0)
print("Average CC of Component 2: ", sum(cc2)/1000.0)
    