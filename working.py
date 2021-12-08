#!/usr/bin/python3
from matplotlib.pyplot import table
import numpy as np
from pandas.io.pytables import Table
import requests, sys, os, re
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import json
#from Bio.SeqFeature import SeqFeature, FeatureLocation


#Step 1 - Get ensemble ID 
gene = "MC1R"  #redhair gene
server = "https://rest.ensembl.org"
ext = "/lookup/symbol/homo_sapiens/MC1R?format=condensed;db_type=core"

r = requests.get(server+ext, headers={ "Content-Type" : "text/xml"})
if not r.ok:
  r.raise_for_status()
  sys.exit()
print(r.text)

# Mygene id: ENSG00000258839


# Step 2. 
        #2A. Use the ensembl database endpoints to grab the nucleotide sequence data for MC1R. Write this sequence to a fasta file. 

        #Define endpoint
def fetch_endpoint (server, request, content_type): #from https://www.ebi.ac.uk/training/online/sites/ebi.ac.uk.training.online/files/u1218/8_REST_API_Emily.pdf
    r = requests.get(server+request, headers={"Content-Type" : content_type})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    if content_type == 'application/json':
        return r.json()
    else:
        return r.text

ext2 = "/xrefs/symbol/homo_sapiens/" + gene + "?"
        #https://rest.ensembl.org/documentation/info/xref_external
con = "application/json"


        #use ensembl ID to look up fasta file for redhair gene, MCIR
get_lookup = fetch_endpoint(server, ext2, con)
e_ID = get_lookup[0]['id']
print(e_ID)

        #Create fasta file name, pull sequence and save file as redhair.fasta
seq = "/sequence/id/" + e_ID + "?"
get_seq = fetch_endpoint(server, seq, "text/x-fasta")
with open('redhair.fasta', 'w') as text_file:
    text_file.write(get_seq) 
#print(get_seq)

    # 2B.Get the long ORF. 
    # https://stackoverflow.com/questions/31757876/python-find-longest-orf-in-dna-sequence

records = SeqIO.parse('redhair.fasta', 'fasta')

max_pro = ""
for record in records:
    for strand, seq in (1, record.seq), (-1, record.seq.reverse_complement()):
        for frame in range(3):
            index = frame
            while index < len(record) - 6:
                match = re.match('(ATG(?:\S{3})*?T(?:AG|AA|GA))', str(seq[index:]))
                if match:
                    orf = match.group()
                    index += len(orf)
                    if len(orf) > 1300:
                        pos = str(record.seq).find(orf) + 1 
                        long_orf = ("{}...{} - length {}, strand {}, frame {}, pos {}, name {}".format\
                           (orf[:6], orf[-3:], len(orf), strand, frame, pos, record.id))
                else: index += 3
print(long_orf) # prints the longest orf 
seq_orf = str(orf) # prints the entire seq of orf
print(seq_orf)    

        #2C. Translate DNA to amino acid

amino_acid = Seq.translate(orf, to_stop=True)
print(amino_acid)
# amino_acid = orf_seq.translate(to_stop=True)   #http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec25  

        #2D. Add longest ORF to fasta file

redhair = open("redhair.fasta", 'a')
redhair.write("Amino Acid Sequence:\n{}". format(amino_acid))
redhair.close()


#Step 3. Show homologous species for MC1R gene. 
     # https://rest.ensembl.org/documentation/info/homology_symbol
     # https://github.com/Ensembl/rest-api-jupyter-course
con = "application/json"
ext = "/homology/symbol/human/" + "MC1R" + "?"  
get_home = fetch_endpoint(server,ext,con)

hom = set()

for i in get_home['data'][0]['homologies']:
    if i['target']['species'] != 'homo_sapiens':
        hom.add(i['target']['species'])

#sort
homologous = "\n".join(sorted(hom))

#Step 3 : List homologous species of MC1R gene
hoy = open("homology.txt", 'w')
hoy.write("List of Homologous Species to MC1R:\n\n")
hoy.write(homologous)
hoy.close()
