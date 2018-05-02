#### Program to convert excel tables into Latex format
f = open("program.txt", "r")
print "Disorder & Nodes & Edges & Largest CC $(X)$ & Clustering Coefficient & $P(X \geq x)$ \\ [0.5ex] \hline\hline"
for line in f:
	x = line.split('\t')
	print x[0] + " & " + x[1] + " & " + x[2] + " & " + x[3] + " & " + x[4] + " & " + x[5] + "\\ \hline"
f.close()
