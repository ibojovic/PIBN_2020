#!/usr/bin/env python3
import sys
import numpy as np


# Parsing secondary structure file containing 4 columns
# PDBID CHAIN SEQUENCE SECONDARY_STRUCTURE
def parse_ssdata(filename):
	# Parsing ss_data and generete 3 dicitionaries 
	# The dicitionaries are for counting the AA, SS and (AA,SS)	
	d_aa={}
	d_ss={}
	d_aass={}
	f=open(filename,'r')
	for line in f:
		line=line.rstrip()
		cols=line.split('\t')
		if (len(cols[2])!=len(cols[3])): print ('WARNING: Incorrect data.',file=sys.stderr)
		n=len(cols[2])
		for i in range(n):
			aa=cols[2][i]
			ss=cols[3][i]
			d_aa[aa]=d_aa.get(aa,0)+1
			d_ss[ss]=d_ss.get(ss,0)+1
			d_aass[(aa,ss)]=d_aass.get((aa,ss),0)+1
	return d_aa,d_ss,d_aass	


# Calculate the propensity scale for a given SS
def calculate_ssprop(d_aa,d_ss,d_aass,ss_type):
	d_prop={}
	aas=list(d_aa.keys())
	n=float(sum(list(d_aa.values())))
	for aa in aas:
		# Check for zero countings
		if d_aass.get((aa,ss_type),0)==0 or \
			d_aa.get(aa,0)==0 or \
			d_ss.get(ss_type,0)==0: continue
		d_prop[aa]=(d_aass[(aa,ss_type)]/n)/(d_aa[aa]/n*d_ss[ss_type]/n)
	return d_prop


if __name__ == '__main__':
	filename=sys.argv[1]
	ss_type=sys.argv[2]
	d_aa,d_ss,d_aass=parse_ssdata(filename)
	d_prop=calculate_ssprop(d_aa,d_ss,d_aass,ss_type)
  ks=list(d_prop.key())
  for k in ks:
		print ('%s\t%s:\t%.2f' %(ss_type,k,d_prop[k]))
 
