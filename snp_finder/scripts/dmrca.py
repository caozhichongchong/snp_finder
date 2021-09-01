# start
# after round 4 calculate dMRCA
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
import argparse
#import difflib
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
required.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='snp_output/',
                      metavar='snp_output/')
################################################## Definition ########################################################
args = parser.parse_args()
input_script = args.s
output_dir_merge = args.o
vcf_name = '.all.parsi.fasta'

def SNP_seq_diff(seq1, seq2,allSNP_loci):
    # diff to ancestral allele
    SNP_total = 0
    for i in allSNP_loci:
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
    return SNP_total

def SNP_seq_loci(seq1, seq2):
    SNP_loci = set()
    total_length = len(seq1)
    for i in range(0,total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_loci.add(i)
    return SNP_loci

def SNP_loci(all_seq):
    # polymorphic allele loci only
    allSNP_loci = []
    num_seq = len(all_seq)
    for i in range(0,num_seq-1):
        for j in range(i+1,num_seq):
            allSNP_loci += list(SNP_seq_loci(all_seq[i], all_seq[j]))
    return list(set(allSNP_loci))

# calculate dMRCA
all_vcf_file=glob.glob(os.path.join(output_dir_merge,'*%s'%(vcf_name)))
alloutput = []
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split('.all.parsi.fasta')[0]
    SNP_to_ref_all = []
    record_seq = ''
    Ref_seq = ''
    All_seq = dict()
    genome_num = 0
    all_seq = set()
    for record in open(vcf_file,'r'):
        if not record.startswith('   '):
            record_name,record_seq = record.split('\n')[0].split('    ')
            All_seq.setdefault(record_name, record_seq)
            all_seq.add(record_seq)
            genome_num += 1
    if 'Srefer' in All_seq:
        Ref_seq = All_seq['Srefer']
        All_seq.pop('Srefer')
        all_seq.remove(Ref_seq)
        genome_num -= 1
    else:
        Ref_seq = All_seq[record_name] # take the last one
    allSNP_loci = SNP_loci(list(all_seq))# polymorphic allele loci only
    for record_name in All_seq:
        SNP_to_ref_all.append(SNP_seq_diff(Ref_seq, All_seq[record_name],allSNP_loci))# diff to ancestral allele
    alloutput.append('%s\t%s\t%.3f\t%.3f\n'%(donor_species,genome_num,statistics.mean(SNP_to_ref_all),
                                               statistics.stdev(SNP_to_ref_all)))

f1 = open(os.path.join(output_dir_merge + '/summary/', 'alldmrca.txt'),'w')
f1.write('donor_species\tgenome_num\tdmrca\tdmrca_stdev\n' + ''.join(alloutput))
f1.close()

