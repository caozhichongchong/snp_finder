# start
# after round 4 calculate NS ratio
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
required.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='snp_output/',
                      metavar='snp_output/')
################################################## Definition ########################################################
args = parser.parse_args()

output_dir_merge = args.o
sum_name = '.all.parsi.fasta.sum.txt'
vcf_name = '.raw.vcf.filtered.vcf.final.snp.txt'
outsum_name = '.all.parsi.fasta.linktrunc.sum.txt'
try:
    os.mkdir(output_dir_merge + '/nolinktrunc')
except IOError:
    pass
################################################### Function ########################################################
def load_snp(snpfile):
    trunc_set = dict()
    for lines in open(snpfile,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            genename, POS,N_or_S,AAchange = lines_set[-4:]
            if AAchange!='' and '*' == AAchange[1]:
                # truncation SNP
                trunc_set.setdefault(genename,POS)
    return trunc_set

def remove_linked_SNP(trunc_genomeset,temp_output,newoutput):
    for lines in temp_output:
        lines_set = lines.split('\n')[0].split('\t')
        Gene, Gene_POS, Genome_set_noSNP, Genome_SNP_set = lines_set[5:9]
        if Genome_SNP_set!=trunc_genomeset:
            newoutput.append(lines)
        else:
            print('remove linked SNP %s %s  %s'%(donor_species,Gene,Gene_POS))
    return newoutput

def load_sum(sumfile,trunc_set):
    newoutput = []
    trunc_genomeset = ''
    temp_output = []
    oldgene = ''
    for lines in open(sumfile, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        Gene, Gene_POS, Genome_set_noSNP, Genome_SNP_set = lines_set[5:9]
        if oldgene != Gene:
            # a new gene
            newoutput = remove_linked_SNP(trunc_genomeset,temp_output,newoutput)
            temp_output = []
        oldgene = Gene
        if Gene not in trunc_set:
            # genes with no truncation
            newoutput.append(lines)
        else:
            # genes with no truncation
            if Gene_POS == trunc_set[Gene]:
                # truncation site
                trunc_genomeset = Genome_SNP_set
                newoutput.append(lines)
            else:
                # other SNPs
                temp_output.append(lines)
    f1 = open(sumfile.replace(sum_name,outsum_name),'w')
    f1.write(''.join(newoutput))
    f1.close()
    os.system('mv %s %s/nolinktrunc/'%(sumfile,output_dir_merge))

################################################### Main ########################################################
all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0].replace('.all', '')
    outsum_file = glob.glob('%s/%s%s'%(output_dir_merge,donor_species,outsum_name))
    if outsum_file == []:
        sum_file = glob.glob('%s/%s%s'%(output_dir_merge,donor_species,sum_name))[0]
        trunc_set = load_snp(vcf_file)
        load_sum(sum_file, trunc_set)

################################################### END ########################################################
