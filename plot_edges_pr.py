### Plots Number of Edges as a Funtion of Probability
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal
import os
import matplotlib.ticker as ticker



#total_edges = 94797384
total_edges = 1048576

#fig, ax = plt.subplots()
x_values = []
y_values = []

with open(os.path.join('./Tissue_Networks','brain_sorted'), 'rb') as f:
	for line in f:
		line = line.rstrip()
		l = line.split('\t')
		y_values.append(total_edges)
		x_values.append(l[2])
		total_edges = total_edges-1
f.close()

# # Data for plotting
#t = np.arange(0.0, 1.0, 0.005)
x = np.array(x_values)
y = np.array(y_values)
# s = 1 + np.sin(2 * np.pi * t)

# # Note that using plt.subplots below is equivalent to using
# # fig = plt.figure() and then ax = fig.add_subplot(111)
fig, ax = plt.subplots(1)
#plt.ylim(0, total_edges)
ax.plot(x, y)
#plt.xticks(np.arange(0, 1.0, 0.1))
ax.xaxis.set_ticks(np.arange(0.0, 1.0, 0.1))
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

#plt.scatter(t,y_values)

ax.set(xlabel='Probability (p-value)', ylabel='Edges',
       title='Number of Edges as a Funtion of Probability')
ax.grid()

#fig.savefig("test.png")
plt.show()