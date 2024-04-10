import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

nodes=50
ws_net=nx.watts_strogatz_graph(nodes,2,0.5)
mean_deg=np.mean(list(dict(ws_net.degree()).values()))
list_bets=list(nx.betweenness_centrality_source(ws_net,normalized=False).values())
mean_bet=np.mean(list_bets)
mean_tra=np.mean(list(nx.clustering(ws_net).values()))
print ('<DEGREE>:',mean_deg)
print ('<BETWEENNESS>:',mean_bet)
print ('<TRANSITIVITY>: %.3f' %mean_tra)
nx.draw(ws_net,with_labels = True)
plt.show()
