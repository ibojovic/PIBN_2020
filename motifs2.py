#!/usr/bin/env python3
import sys
import networkx as nx
from networkx.algorithms import isomorphism
import matplotlib.pyplot as plt

# Function for parsing the tsv file

def read_net(filename):
	list_edges=[]
	f=open(filename,'r')
	for line in f:
		# Do not read header lines
		if line[0]=='#': continue
		line=line.rstrip('\n')
		v=line.split('\t')
		# Select data with Strong supporting evidences
		# if v[-1]!='Strong': continue
		# Define the sign of the edges
		if v[4]=='-':
			sign='-'
		elif v[4]=='+':
			sign='+'
		else:
			continue
		list_edges.append((v[1],v[3],sign))
	return list_edges


# Define the regulatory network assigning
# attributes to node and edges

def get_net(list_edges):
	g=nx.MultiDiGraph()
	for e in list_edges:
		g.add_node(e[0],color='blue')
		g.add_node(e[1],color='red')
		g.add_edge(e[0],e[1],sign=e[2])
	return g


def motif_dpf():
	dpf=nx.MultiDiGraph()
	dpf.add_nodes_from([0,1],color='blue')
	dpf.add_nodes_from([2,3],color='red')
	dpf.add_edges_from([(0,2),(0,3),(1,2),(1,3)],sign='+')
	return dpf

	
def motif_mim():
	mim=nx.MultiDiGraph()
	mim.add_nodes_from([0,1],color='blue')
	mim.add_nodes_from([2,3,4],color='red')
	mim.add_edges_from([(0,2),(0,3),(0,4),(1,2),(1,3),(1,4)],sign='+')
	return mim


def show_graph(g):
	colors=[]
	for n in g.nodes:
		colors.append(g.nodes[n]['color'])
	nx.draw(g,node_color=colors)
	plt.show()


def get_max_degree(g,list_nodes):
	degs=[(g.degree(node),node) for node in list_nodes]
	degs.sort()
	degs.reverse()
	return degs[0][1],degs[0][0]
	
	
def get_elements(list_edges,pos):
	d={}
	for e in list_edges:
		d[e[pos]]=True
	return list(d.keys())
	

def get_effect(g,node,effect,out=True):
	i=0
	if out:
		edges=g.out_edges(node)
	else:
		edges=g.in_edges(node)
	edges=list(set(edges))
	for e in edges:
		# The next line can change depending 
		# on the version of networkx.
		se=g[e[0]][e[1]]
		for k in list(se.keys()):
			if se[k]['sign']==effect: i=i+1
	return i


def match_isomorphisms(net,mot):
	diso={}
	em=isomorphism.categorical_multiedge_match('sign','+')
	nm=isomorphism.categorical_node_match('color','blue')
	ms=nx.isomorphism.DiGraphMatcher(net,mot,edge_match=em, node_match=nm)
	i=0
	# Save sorted matching key in dictionary to exclude duplicates
	for m in ms.subgraph_isomorphisms_iter():
		k=list(m.keys())
		k.sort()
		diso[tuple(k)]=True
		i=i+1
	i=len(list(diso.keys()))
	return i


if __name__ == '__main__':
	filename=sys.argv[1]
	list_edges=read_net(filename)
	list_tf=get_elements(list_edges,0)
	list_gene=get_elements(list_edges,1)
	g=get_net(list_edges)
	print ('Regulation Network')
	print (' Number of TFs: %d' %len(list_tf))
	print (' Number of Genes: %d' %len(list_gene))
	print (' Number of Edges: %d' %len(list_edges))  
	max_tf=get_max_degree(g,list_tf)
	max_gene=get_max_degree(g,list_gene)
	act_tf=get_effect(g,max_tf[0],'+')
	rep_tf=get_effect(g,max_tf[0],'-')
	act_gene=get_effect(g,max_gene[0],'+',False)
	rep_gene=get_effect(g,max_gene[0],'-',False)
	print ('Max Degree TF: %s' %max_tf)
	print ('   Activation: %s' %act_tf)
	print ('   Repression: %s' %rep_tf)
	print ('Max Degree Gene: %s' %max_gene)
	print ('   Activation: %s' %act_gene)
	print ('   Repression: %s' %rep_gene)
	# Define the Double Positive Feedback Loop
	dpf=motif_dpf()
	n_dpf=match_isomorphisms(g,dpf)
	print ('Double Positive Feedback: %d' %n_dpf)
	# Define the Multi-Input Module
	mim=motif_mim()
	n_mim=match_isomorphisms(g,mim)
	print ('Multi-Input Module: %d' %n_mim)
	show_graph(g)



	
		
