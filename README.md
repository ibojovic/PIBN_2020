# Protein Sequence and Structure Data

## Protein Sequence Data
From the UniProt FTP web site (ftp://ftp.expasy.org/databases/uniprot/) download the Human protein UP000005640_9606 in fasta format.

```bash
wget ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz

# Total numer of SwissProt protein
echo -e "Annotated Human Protein: \c"
zcat data/UP000005640_9606.fasta.gz |grep ">sp" | wc -l

# Total number of TrEMBL protein
echo -e "Human Protein TrEMBL: \c"
zcat data/UP000005640_9606.fasta.gz |grep ">tr" | wc -l 

# Download P53 file sequence
wget https://www.uniprot.org/uniprot/P04637.fasta

# Count residues
echo -e "P53 Sequence Length: \c"
grep -v ">" data/P04637.fasta | tr -d "\n" |wc -c

```

## Protein Structure Data
PDB file of the Ribonuclease A (PDB: 7RSA) from the web (http://ftp.rcsb.org/pub/pdb/data/structures/all/pdb/pdb7rsa.ent.gz)

```bash
# Download the PDB file
wget http://ftp.rcsb.org/pub/pdb/data/structures/all/pdb/pdb7rsa.ent.gz

# Count the CA atmes
echo -e "Residues in PDB Structure: \c"
zcat data/pdb7rsa.ent.gz |grep ^ATOM |grep -w CA |wc -l
```

```bash
# Download the DSSP file
wget ftp://ftp.cmbi.umcn.nl//pub/molbio/data/dssp/7rsa.dssp
```
```bash
# Run the script calculating the distance between the firt two CA atoms
python3 script/parse_dssp.py data/7rsa.dssp A
```
```bash
# Run the script calculating the distance between the firt two CA atoms
!python3 script/parse_pdb.py data/pdb7rsa.ent A
```

## Amino Acid Propensities
### Secondary structure propensity scale
Develop alpha helix propensity scale based on the data from **ss_data.tsv.gz**

```bash
# Run the script calculating the popensity scale
python3 script/ss_scale.py data/ss_data.tsv H
```
## Protein-Protein Interaction

Download the DSSP file of the Bacterial luciferase (Vibrio harveyi) from the PDB (code: 1BRL)
```bash
# Download the PDB file
wget http://ftp.rcsb.org/pub/pdb/data/structures/all/pdb/pdb1brl.ent.gz
mv pdb1brl.ent.gz data

# Unzip the pdb file
gunzip data/pdb1brl.ent.gz

# Extract the chain A
grep -w A data/pdb1brl.ent >data/1brlA.pdb

# Extract the chain B
grep -w B data/pdb1brl.ent >data/1brlB.pdb

# Merge chains A and B
cat data/1brlA.pdb data/1brlB.pdb >data/1brlAB.pdb

# Generate the DSSP file for chain A
script/dssp data/1brlA.pdb >data/1brlA.dssp

# Generate the DSSP file for chain B
script/dssp data/1brlB.pdb >data/1brlB.dssp

# Generate the DSSP file for the complex AB
script/dssp data/1brlAB.pdb >data/1brlAB.dssp
```

```bash
# Run the script on the 1brlA.dssp file and redirect the output
python3 script/parse_dssp_ppi.py data/1brlA.dssp A >data/m1brlA.acc

# Run the script on the 1brlB.dssp file and redirect the output
python3 script/parse_dssp_ppi.py data/1brlB.dssp B >data/m1brlB.acc

# Run the script on the 1brlAB.dssp file and get the accessibility of chain A
python3 script/parse_dssp_ppi.py data/1brlAB.dssp A >data/c1brlA.acc

# Run the script on the 1brlAB.dssp file and get the accessibility of chain B
python3 script/parse_dssp_ppi.py data/1brlAB.dssp B >data/c1brlB.acc
```

```bash
# Show the Total surface of the monomer A
tail -n 1 data/m1brlA.acc

# Show the Total surface of the monomer B
tail -n 1 data/m1brlB.acc

# Show the Total surface of the chain A in complex
tail -n 1 data/c1brlA.acc

# Show the Total surface of the chain B in complex
tail -n 1 data/c1brlB.acc
```

The surface of interaction is:

> $Interaction Surface = [S(A) + S(B) + S(AB)]/2$

where S(AB) is the sum of S(A) + S(B) obtaind for the structure of the protein complex (1brlAB.pdb)

```bash
# Compare chains A in monomeric and complex forms showing only residues with differecne in accessibility greater than 50
paste data/m1brlA.acc data/c1brlA.acc |awk -v OFS='\t' '{if ($4-$10>50) print $1,$2,$4,$10,$4-$10}' |sort -nrk 5  |grep -v TOTAL
```

```bash
# Compare chains B in monomeric and complex forms showing only residues with differecne in accessibility greater than 50
paste data/m1brlB.acc data/c1brlB.acc |awk -v OFS='\t' '{if ($4-$10>50) print $1,$2,$4,$10,$4-$10}' |sort -nrk 5 |grep -v TOTAL
```


## Analysis of Protein-Protein Interaction Network

```bash
# Download the IntAct PSI-MITAB file
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip
mv intact.zip data

# Search for the inteactions of MEKK1 protein using grep and count them
echo -e "Number of interactions of MEKK1: \c"
unzip -c data/intact.zip |grep 'uniprotkb:MEKK1' > data/mekk1.txt
cat data/mekk1.txt | wc -l

# From the previous search, count the number of unique interactions
echo -e "Unique interactions of MEKK1: \c"
cat data/mekk1.txt | awk -F '\t' '{print $1"\t"$2}' |sort -u |wc -l

# Remove interactions with drugs (chebi:) and other non proteins (intact:)
echo -e "Unique interactions of MEKK1 with proteins: \c"
cat data/mekk1.txt | awk -F '\t' '{print $1"\t"$2}' | grep -v 'chebi:' | grep -v 'intact:' | sort -u | wc -l

# Select only interactions correponiding to the human MEKK1 searching for the swissprot identifier Q13233
echo -e "Unique interactions of the human MEKK1 with proteins: \c"
grep  'uniprotkb:Q13233' data/mekk1.txt | awk -F '\t' '{print $1"\t"$2}' | grep -v 'chebi:' | grep -v 'intact:' | sort -u | wc -l

# Calculate the number of inteaction of BRAF (P15056)
echo -e "Unique interactions of BRAF with proteins: \c"
unzip -c data/intact.zip |grep 'uniprotkb:P15056' | awk -F '\t' '{print $1"\t"$2}' >data/braf.txt
sort -u data/braf.txt | wc -l

# Is BRAF inteacting with MEKK1? How many experiments are supporting this interactions?
echo -e "Number of interactions between MEKK1 and BRAF: \c"
cat data/braf.txt |grep 'uniprotkb:P15056' | grep  'uniprotkb:Q13233' | wc -l 

```

```bash
# Using the TaxID of the Homo sapiens (9606) calculate the numeber of total interactions between human proteins
unzip -c data/intact.zip | awk -F '\t' '{split($10,org1,"|"); split($11,org2,"|"); if (org1[1]=="taxid:9606(human)" && org2[1]=="taxid:9606(human)") print $0}' >data/human.txt
echo -e "Number of interactions between human proteins: \c"
cat data/human.txt | wc -l

# What is the number of unique interactions considering the problem of the directionality A->B is equal to B->A
echo -e "Unique interactions between human proteins: \c"
awk -F '\t' '{ if ($1>$2) {print $2"\t"$1} else {print $1"\t"$2} }'  data/human.txt | sort -u | wc -l

```

## Graph theory and Networks

networkx python library was used to generate the following networks:

* Random network (Erdos Renyi model)
* Small-world network (Watts-Strogatz model)
* Scale-free network (Barabasi-Albert model)

### Analysis of the Yeast interactome

```bash
# Download the zip file of all interactomes
%env web=https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.203/BIOGRID-ORGANISM-4.4.203.mitab.zip 
wget $web -P data

# unzip only the Saccharomyces file
%env yeast=BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.4.203.mitab.txt
unzip -p  data/BIOGRID-ORGANISM-4.4.203.mitab.zip $yeast >data/$yeast
```
```bash
# Run the above script giving in input 
# the mitab file and the number of bins.
python3 script/biogrid-net.py data/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.4.203.mitab.txt 37
```
```bash
from IPython.display import Image
Image('data/yeast.png')
```
























