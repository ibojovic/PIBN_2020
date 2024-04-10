#!/usr/bin/env python3
import networkx as nx
import matplotlib.pyplot as plt
import itertools as it
import numpy as np
from scipy.stats import norm, shapiro


def count_isomorphisms(net,mot):
	n=list(net.nodes())
	k=mot.number_of_nodes()
	combs=list(it.combinations(n,k))
	i=0
	for comb in combs:
		sg=net.subgraph(comb)
		if nx.is_isomorphic(sg,mot): i=i+1
	return i


def match_isomorphisms(net,mot):
	diso={}
	ms=nx.isomorphism.DiGraphMatcher(net,mot)
	for m in ms.subgraph_isomorphisms_iter():
		k=list(m.keys())
		k.sort()
		diso[tuple(k)]=True
	i=len(list(diso.keys()))
	return i


def check_random(nn,ne,mot,k):
	occ=[]
	for i in range(k):
		g=nx.gnm_random_graph(nn,ne,directed=True)
		c=count_isomorphisms(g,mot)
		occ.append(c)
	return occ


def get_norm_statistics(net,mot,k=100):
	nn=net.number_of_nodes()
	ne=net.number_of_edges()
	# Here we used the maching function 
	# based on DiGraphMatcher replacing
	# n_match=count_isomorphisms(net,mot)
	# with
	n_match=match_isomorphisms(net,mot)
	# Generate k random networks and 
	# get their number of matches
	r_match=check_random(nn,ne,mot,k)
	# Calculate mean, stdev and z-score
	mean=np.mean(r_match)
	std=np.std(r_match)
	z=(n_match-mean)/std
	# Calculate the p-value assuming that  
	# r_match values have a normal distribution
	p=norm.sf(n_match,mean,std)
	# Shapiro test for testing tha r_match values
	# have a normal distribution
	pshapiro=shapiro(r_match)[1]
	# Calculate Empirical probability
	ep=len([i for i in r_match if i>n_match])/float(len(r_match))
	print ('Matches:',n_match)
	print ('Z-score: %.3f' %z)
	print ('P-value %.3e' %p)
	print ('Shapiro: %.3e' %pshapiro)
	print ('P(t>%d): %.3e\n' %(n_match,ep))
 

def main():
	k=1000
	# Define Feed Forward Loop
	ffl=nx.DiGraph()
	ffl.add_edges_from([('A','B'),('A','C'),('B','C')])
	# Define 3-node Cycle
	c3=nx.DiGraph()
	c3.add_edges_from([('A','B'),('B','C'),('C','A')])
	# Define the network
	g=nx.DiGraph()
	g.add_edges_from([(1,1),(2,1),(1,3),(2,3),(3,4),(4,2),(4,5),(4,6),(5,6)])
	print ('# Match Feed-Forvard Loop')
	get_norm_statistics(g,ffl,k)
	print ('# Match 3-Node Cycle')
	get_norm_statistics(g,c3,k)
	# Modify the network removing edge (1,1) and changing (3,4) with (4,3).
	g.remove_edge(1,1)
	g.remove_edge(3,4)
	g.add_edge(4,3)
 	# Test the significance of the occurrence of the Feed-Forvard Loop
	print ('# Match Feed-Forvard Loop')
	get_norm_statistics(g,ffl,k)


if __name__ == '__main__':
	main()
