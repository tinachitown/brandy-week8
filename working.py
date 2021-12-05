import requests, sys, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import json


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


# Step 2. Use the ensembl database endpoints to grab the nucleotide sequence data for MC1R. Write this sequence to a fasta file. 
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
server = "http://rest.ensembl.org"
ext2 = "/xrefs/symbol/homo_sapiens/" + gene + "?"
#https://rest.ensembl.org/documentation/info/xref_external
con = "application/json"


#get ensembl ID
get_lookup = fetch_endpoint(server, ext2, con)
e_ID = get_lookup[0]['id']
print(e_ID)


#ex1 = "/sequec/id/ENSG00000258839?type=genomic"  #ensembl id

# record = SeqIO.read("redhair.fasta", "fasta")
# table = 11
# min_pro_len = 100

# max_len_pro = 0
# max_pro = ""
# for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
#     for frame in range(3):
#       length = 3 * ((len(record)-frame) // 3) #Multiple of three
#       for pro in nuc[frame:frame+length].translate(table).split("*"): # to translate into amino acids
#         if len(pro) >= min_pro_len:
#           print("%s length %i, strand %i, frame %i" 
#             % (pro, len(pro), strand, frame)) 
#           if len(pro) > max_len_pro:
#                 max_len_pro = len(pro)
#                 max_pro = pro
# print(max_len_pro, max_pro)
# mpro = max_pro






# with open('redhair.fasta', 'w') as text_file:
#   text_file.write(r.text)


# #look up information for gene mc1r 

# record = SeqIO.read("red.fasta", "fasta")
# table = 11
# min_pro_len = 100



# # Step 3. Show homologous species for MC1R gene. 
 
# con = "application/json"
# ext = "/homology/symbol/human/" + "MC1R" + "?"
# get_home = endpoint(server,ext,con)

# hom = set()

# for i in get_home['data'][0]['homologies']:
#     if i['target']['species'] != 'homo_sapiens':
#         hom.add(i['target']['species'])

# #sort
# sort_hom = sorted(hom)
# str_hom = "\n".join(sort_hom)



# #Step 3 : List homologous species of MC1R gene
# hoy = open("homology.txt", 'w')
# hoy.write("List of Homologous Species to MC1R:\n\n")
# hoy.write(str_hom)
# hoy.close()
