import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

nodes=50
er_net=nx.erdos_renyi_graph(nodes,0.04)
mean_deg=np.mean(list(dict(er_net.degree()).values()))
list_bets=list(nx.betweenness_centrality_source(er_net,normalized=False).values())
mean_bet=np.mean(list_bets)
mean_tra=np.mean(list(nx.clustering(er_net).values()))
print ('<DEGREE>:',mean_deg)
print ('<BETWEENNESS>:',mean_bet)
print ('<TRANSITIVITY>: %.3f' %mean_tra)
nx.draw(er_net,with_labels = True)
plt.show()
