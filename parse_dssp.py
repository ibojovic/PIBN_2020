#!/usr/bin/env python3
import sys
import numpy as np

# Function parsing the PDB file
def parse_pdb(filename, chain):
  pdb_coords={}
  f=open(filename,'r')
  for line in f:
      line=line.rstrip()
      # Check the field ATOM at the beginning of the line
      # Check the correct chain in column 22
      if line[:4]!='ATOM' or line[21]!=chain: continue
      resn=int(line[22:26].strip())
      atom=line[12:16].strip()
      x=float(line[30:38])
      y=float(line[38:46])
      z=float(line[46:54])
      coord=[x,y,z]
      # Initialize the dictionary with atoms coordinates of residue resn
      pdb_coords[resn]=pdb_coords.get(resn,{})
      # Add the atom's coordinates
      pdb_coords[resn][atom]=coord
  return pdb_coords


# Function calculating the distance between to points
def get_distance(coord1,coord2):
  return np.sqrt((coord1[0]-coord2[0])**2+\
                 (coord1[1]-coord2[1])**2+\
                 (coord1[2]-coord2[2])**2)


# Function returning the diatnces between all consecutives CA
def get_ca_dist(pdb_coords):
  ca_dists=[]
  keys=list(pdb_coords.keys())
  keys.sort()
  n=len(keys)
  for i in range(n-1):
    dist=get_distance(pdb_coords[keys[i]]['CA'],pdb_coords[keys[i+1]]['CA'])
    ca_dists.append(dist)
  return ca_dists


if __name__ == '__main__':
  filename=sys.argv[1]
  chain=sys.argv[2]
  pdb_coords=parse_pdb(filename,chain)
  ca_dists=get_ca_dist(pdb_coords)
  print ('Dist CA:',np.mean(ca_dists),np.std(ca_dists))))
  
