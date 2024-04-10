import networkx as nx
import matplotlib.pyplot as plt

g=nx.Graph()
g.add_edges_from([(1,2),(1,3),(1,4),(1,5),(2,3),(4,5)])
c1=nx.clustering(g,1)
d1=g.degree(1)
print ('D(1)= %d' %d1)
print ('C(1)= %.3f' %c1)
colors=['red']+4*['blue']
nx.draw(g,with_labels=True,node_color=colors)
plt.show()
