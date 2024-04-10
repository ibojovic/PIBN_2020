#!/usr/bin/env python3
import sys
import numpy as np

Norm_Acc= {"A" :106.0,  "B" :160.0,         # D or N
   "C" :135.0,  "D" :163.0,  "E" :194.0,
   "F" :197.0,  "G" : 84.0,  "H" :184.0,
   "I" :169.0,  "K" :205.0,  "L" :164.0,
   "M" :188.0,  "N" :157.0,  "P" :136.0,
   "Q" :198.0,  "R" :248.0,  "S" :130.0,
   "T" :142.0,  "V" :142.0,  "W" :227.0,
   "X" :180.0,         # undetermined (deliberate)
   "Y" :222.0,  "Z" :196.0}         # E or Q


# Function parsing the DSSP file
def parse_dssp(filename,chain):
  dssp_list=[]
	# Initialize state variable c=0
  c=0
  f=open(filename)
  for line in f:
    line=line.rstrip()
		# Check for the beginning of data
		# Change the state variable c=1
    if line.find('  #  RESIDUE')>-1: 
      c=1
      continue
		# Check if the state variable is 1
		# Select the correct chain
    if c==0 or line[11]!=chain: continue
		# Read all the variables resn, aa, acc, phi, psi
    aa=line[13]
    if aa.islower(): aa='C'
		if aa=='!': continue
    resn=int(line[5:10])
    ss=line[16]
    if ss==' ': ss='C'
    acc=float(line[34:38])
    phi=float(line[103:109])
    psi=float(line[109:115])
    aa_dssp=[resn,aa,ss,acc,phi,psi]
    dssp_list.append(aa_dssp)
  return dssp_list
		

# Get sequence and secondary structure
def get_ss(dssp_list):		
	seq=''
	ss=''
	for aa in dssp_list:
		ss=ss+aa[2]
		seq=seq+aa[1]
	return seq,ss
	

# Count residues with a given secondary structure
def count_ss(ss,ss_type):
	c=0
	for i in ss:
		if i==ss_type: c=c+1
	return c 
	

# Return the list of all phi and psi angles
def get_ss_angle(dssp_list,ss_type):
	angles=[]
	for aa in dssp_list:
		if aa[2]==ss_type: angles.append((aa[-2],aa[-1]))
	return angles


# Return the list of residues accessibility
def get_aa_acc(dssp_list,aa_type):
	accs=[]
	for aa in dssp_list:
		if aa[1]==aa_type: accs.append(aa[3])
	return accs


def get_chain_acc(dssp_list):
	accs=[]
	for aa in dssp_list:
		accs.append(aa[3])
	return accs



if __name__ == '__main__':
	filename=sys.argv[1]
	chain=sys.argv[2]
	dssp_list=parse_dssp(filename,chain)
	seq,ss=get_ss(dssp_list)
	print ('>SEQ\n%s' %seq)
	print ('>SS\n%s' %ss)
	print ('H count:',count_ss(ss,'H'))
	print ('E count:',count_ss(ss,'E'))
	h_angles=get_ss_angle(dssp_list,'H')
	e_angles=get_ss_angle(dssp_list,'E')
	print ('H angles:',np.mean([i[0] for i in h_angles]),np.mean([i[1] for i in h_angles]))
	print ('E angles:',np.mean([i[0] for i in e_angles]),np.mean([i[1] for i in e_angles]))
	print ('K acc:',np.mean(get_aa_acc(dssp_list,'K'))/Norm_Acc['K'])
	print ('V acc:',np.mean(get_aa_acc(dssp_list,'V'))/Norm_Acc['V'])
