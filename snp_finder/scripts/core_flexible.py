# core genome or flexible genome
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of all WGS genomes",
                      type=str, default='.',
                      metavar='input/')
# optional output setup
optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='snp_output/',
                      metavar='snp_output/')
# optional search parameters
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 1)",
                      metavar="1 or more", action='store', default=1, type=int)
optional.add_argument('-job',
                      help="Optional: command to submit jobs",
                      metavar="nohup or customized",
                      action='store', default='jobmit', type=str)
# requirement for software calling
optional.add_argument('-pro',
                      help="Optional: complete path to prodigal if not in PATH",
                      metavar="/usr/local/bin/prodigal",
                      action='store', default='prodigal', type=str)
optional.add_argument('--dm',
                        help="Optional: complete path to diamond if not in PATH",
                        metavar="/usr/local/bin/diamond",
                        action='store', default='diamond', type=str)


################################################## Definition ########################################################
args = parser.parse_args()

input_script_sub = args.s +'/pangenome'
input_script = args.s
mutation_dir = args.o + '/summary'
mutation_fasta = os.path.join(mutation_dir,'all.selected.gene.faa.High_select2.faa')
mutation_fasta2 = os.path.join(mutation_dir,'all.denovo.gene.faa')
genome_dir = glob.glob(args.i + '/*')
output_dir = args.o + '/pangenome'

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
            os.system('%s -q -i %s -a %s' % (args.pro,fasta, fasta + '.faa'))
        fasta = fasta + '.faa'
        cmds += ( "%s blastp --query %s --db %s.dmnd --out %s.denovo.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                    % (args.dm, fasta, mutation_fasta2, os.path.join(output_dir, fastaname), cutoff, cutoff2))
    if cmds!= '':
        f1 = open(os.path.join(input_script_sub, '%s.denovo.sh' % donor_species), 'a')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
        f1.close()

f1 = open(os.path.join(input_script, 'pangenome.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.denovo.sh')):
    if 'jobmit' in args.job:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    else:
        f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/%s'%(input_script,'pangenome.sh'))