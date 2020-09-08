# start
# simulate adaptive genes
# calculate gene number for each clonal population
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics

# set up path -> selected
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round*/*')+\
glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/round*/*')
output = os.path.join(input_script,'clonal_pop_gene_num.txt')

for folder in genome_dir:
    allfasta = glob.glob(os.path.join(folder,'*.fasta.corrected.fasta.faa'))
    if allfasta!= []:
        allfasta = allfasta[0]
        os.system('echo %s:: >> %s'%(os.path.split(folder)[-1],output))
        os.system('grep \">\" %s | wc -l >> %s'%(allfasta,output))

################################################### END ########################################################
