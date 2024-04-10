import networkx as nx
import matplotlib.pyplot as plt 
import numpy as np

nodes=50
ba_net=nx.barabasi_albert_graph(nodes,1)

mean_deg=np.mean(list(dict(ba_net.degree()).values()))
list_bets=list(nx.betweenness_centrality_source(ba_net,normalized=False).values())
mean_bet=np.mean(list_bets)
mean_tra=np.mean(list(nx.clustering(ba_net).values()))
print ('<DEGREE>:',mean_deg)
print ('<BETWEENNESS>:',mean_bet)
print ('<TRANSITIVITY>: %.3f' %mean_tra)
nx.draw(ba_net,with_labels = True)
plt.show()
