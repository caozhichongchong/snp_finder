################################################### END ########################################################
################################################### SET PATH ########################################################
# summarize truncated genes
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

all_fasta_truc = '/scratch/users/anniz44/genomes/donor_species/IBD/vcf_round1/merge/summary/all.trunc.gene.faa'
Donor_species = dict()
for record in SeqIO.parse(all_fasta_truc, 'fasta'):
    record_id = str(record.id)
    cluster = record_id.split('__')[0].replace('CL','clustercluster')
    Donor_species.setdefault(cluster,0)
    Donor_species[cluster] += 1

Output = []
for cluster in Donor_species:
    Output.append('%s\t%s\n'%(cluster,Donor_species[cluster]))

foutput = open(all_fasta_truc + '.sum', 'w')
foutput.write('donor_species\tNo.trunc\n')
foutput.write(''.join(Output))
foutput.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
