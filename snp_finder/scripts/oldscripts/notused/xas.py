# start
# core genome or flexible genome
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/pangenome'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
mutation_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome/summary'
mutation_fasta = os.path.join(mutation_dir,'all.selected.gene.faa.High_select2.faa')
mutation_fasta2 = os.path.join(mutation_dir,'all.denovo.gene.faa')
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/round*/*')+\
             glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round*/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome/pangenome'

try:
    os.mkdir(output_dir)
except IOError:
    pass

os.system('rm -rf %s'%(input_script_sub))
try:
    os.mkdir(input_script_sub)
except IOError:
    pass

cutoff = 50
cutoff2 = 80
genome_num = 0
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor_species_set = donor_species.split('_')
    donor = donor_species_set[0]
    cmds = ''
    allfasta = glob.glob(os.path.join(folder,'*.corrected.fasta'))
    for fasta in allfasta:
        fastaname = os.path.split(fasta)[-1]
        try:
            f1 = open(fasta + '.faa', 'r')
        except FileNotFoundError:
            os.system('prodigal -q -i %s -a %s' % (fasta, fasta + '.faa'))
        fasta = fasta + '.faa'
        cmds += (
                    "diamond blastp --query %s --db %s.dmnd --out %s.denovo.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                    % (fasta, mutation_fasta2, os.path.join(output_dir, fastaname), cutoff, cutoff2))
    if cmds!= '':
        f1 = open(os.path.join(input_script_sub, '%s.denovo.sh' % donor_species), 'a')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
        f1.close()

f1 = open(os.path.join(input_script, 'pangenome.denovo.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.denovo.sh')):
    if 'jobmit' in args.job:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    else:
        f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
