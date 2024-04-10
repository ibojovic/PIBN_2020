
#N=[1,2,3,4,5,6,7,8,9]
#E=[(1,2),(1,3),(1,4),(1,5),(1,6),(2,7),(2,8),(2,9)]

%matplotlib inline
import networkx as nx
import matplotlib.pyplot as plt

g=nx.Graph()
g.add_edges_from([(1,2),(1,3),(1,4),(1,5),(1,6),(2,7),(2,8),(2,9)])
s1=nx.betweenness_centrality_source(g,normalized=False)[1]
d1=g.degree(1)
print ('D(1)= %d' %d1)
print ('S(1)= %.0f' %s1)
colors=['red']+8*['blue']
nx.draw(g,with_labels=True,node_color=colors)
plt.show()
