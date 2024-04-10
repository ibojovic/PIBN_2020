#!/usr/bin/env python3
import sys
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

def parse_mitab(filename,id1=2,id2=3):
	gi={}
	f=open(filename,'r')
	for line in f:
		if line[0]=='#': continue
		line=line.rstrip()
		v=line.split('\t')
		gid1=v[id1].split('|')[1].split(':')[1]
		gid2=v[id2].split('|')[1].split(':')[1]
		gids=[gid1,gid2]
		# sort gene identifiers to generate undirected edges
		gids.sort()
		gi[tuple(gids)]=True
	return list(gi.keys())
	
	
def get_net(edges):
	g=nx.Graph()
	g.add_edges_from(edges)
	return g


def print_net_features(net):
	# Search for components
	comps=nx.connected_components(net)
	i=1
	for comp in comps:
		print ('Component:',i)
		# Select subgraph
		sg=net.subgraph(comp)
		print ('\t#Nodes:',sg.number_of_nodes())
		print ('\t#Edges:',sg.number_of_edges())
		i=i+1
	# Calculate network features
	dgs=net.degree()
	print ('MEAN DEGREE:',np.mean(list(dgs)))
	cs=nx.clustering(net)
	print ('MEAN CLUSTERING:',np.mean(list(cs.values())))
	# Betweenness for all nodes is too long to be clalculated. Consider only 100 nodes
	bs=nx.betweenness_centrality_source(net,100,normalized=False)
	print ('MEAN BETWEENNES: %.1f' %np.mean(list(bs.values())))
	# Find gene with highest degree
	list_dgs=[(v,k) for k,v in dgs.items()]
	list_dgs.sort()
	list_dgs.reverse()
	print ('MAX DEGREE:',list_dgs[0][1],list_dgs[0][0])



def fit_deg_distribution(degs_vals,nbin,outfile='figure'):
	hist,bins=np.histogram(degs_vals,nbin)
	# Transform the power law function in linear
	# calculating the log10 of x and y
	logx=[]
	logy=[]
	for i in range(len(hist)):
		if hist[i]==0: continue
		logy.append(np.log10(hist[i]))
		logx.append(np.log10((bins[i+1]+bins[i])/2))	
	reg=linregress(logx,logy)
	print ('CORRELATION: %.3f' %reg[2])
	print ('p-value: %.3e' %reg[3])
	x=[]
	yf=[]
	y=[]
	# Calculate the inverse fuction to transform
	# the linear function to power law.
	for i in range(len(logx)):
		x.append(10**logx[i])
		y.append(10**logy[i])
		yf.append(10**(logx[i]*reg[0]+reg[1]))
	plt.plot(x,y,'o')
	plt.plot(x,yf)
	plt.xscale('log')
	plt.yscale('log')
	plt.savefig(outfile+'.png')

	
if __name__ == '__main__':
	filename=sys.argv[1]
	nbin=int(sys.argv[2])
	gi=parse_mitab(filename)
	net=get_net(gi)
	print_net_features(net)
	degs_vals=list(net.degree().values())
	fit_deg_distribution(degs_vals,nbin)
