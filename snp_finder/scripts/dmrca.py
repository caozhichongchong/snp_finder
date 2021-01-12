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
vcf_name = '.raw.vcf.filtered.vcf.final.removerec.fasta'

def SNP_seq_diff(seq1, seq2):
    SNP_total = 0
    j = 0
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            j = i
    return SNP_total

# calculate dMRCA
all_vcf_file=glob.glob(os.path.join(output_dir_merge,'*%s'%(vcf_name)))
alloutput = []
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split('.raw.vcf')[0].replace('.all','').replace('1_PB','1_PaDi')
    genome_num = 0
    SNP_to_ref_all = []
    record_seq = ''
    Ref_seq = ''
    for record in SeqIO.parse(vcf_file, 'fasta'):
        record_name = str(record.id)
        record_seq = str(record.seq)
        if 'reference' in record_name:
            Ref_seq = record_seq
        else:
            SNP_to_ref_all.append(SNP_seq_diff(Ref_seq, record_seq))
            genome_num += 1
    if genome_num > 0:
        alloutput.append('%s\t%s\t%.3f\t%.3f\n'%(donor_species,genome_num,statistics.mean(SNP_to_ref_all),
                                               statistics.stdev(SNP_to_ref_all)))

f1 = open(os.path.join(output_dir_merge + '/summary/', 'alldmrca.txt'),'w')
f1.write('donor_species\tgenome_num\tdmrca\tdmrca_stdev\n' + ''.join(alloutput))
f1.close()

