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
required.add_argument("-i",
                      help="input fasta",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/WGS/vcf_round2/merge/summary/all.denovo.gene.fna',
                      metavar='all.denovo.gene.fna')

################################################## Definition ########################################################
args = parser.parse_args()

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
trunc_set = dict()
trunc_set['*']=0
trunc_set['']=1
Analyze = True
################################################### new class #########################################################
__metaclass__ = type

class SNP_gene:
    # create a class to store SNP_gene
    'a class to store SNP_gene'
    def init(self, gene):
        self.gene = gene
        self.NSratio = [0,0]
    def addSNP_pair(self, temp_NorS):
        self.temp_NorS = temp_NorS
        for codonchange in temp_NorS:
            self.NSratio[trunc_set[codonchange]] += 1
    def sumtrunc(self):
        temp_line = '%s\t%s\t%s\t%s\n'%(self.gene,self.NSratio[0],self.NSratio[1],self.NSratio[0]/(self.NSratio[0]+self.NSratio[1]))
        return temp_line

################################################### Function ########################################################

def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def dnORds(amino1, amino2):
    if amino1 != amino2 and amino2 == '*':
        return '*'
    else:
        return ''

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def expecttrunc(record_name,record_seq,position=0):
    Total = int(len(record_seq)/3)
    temp_SNP_gene = SNP_gene()
    temp_SNP_gene.init(record_name)
    for i in range(0, (Total - position)):
        codon = record_seq[(i * 3 + position):((i + 1) * 3 + position)]
        try:
            codon_NSratio = codontable_NSratio[codon]
            temp_SNP_gene.addSNP_pair(codon_NSratio.temp_NorS)
            SNP_gene_all.addSNP_pair(codon_NSratio.temp_NorS)
        except KeyError:
            pass
    return temp_SNP_gene

################################################### Prepare ########################################################

# calculate codon freq
codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
codontable_NSratio = dict()
for codon in codontable:
    SNP_gene_temp = SNP_gene()
    SNP_gene_temp.init(codon)
    codontable_NSratio.setdefault(codon, SNP_gene_temp)
    temp_NorS = []
    for position in range(0, 3):
        REF = codon[position]
        temp_NorS1 = []
        for ALT in Allels_order:
            if ALT != REF:
                new_codon = causeSNP(codon, position+1, ALT,0)
                temp_NorS1.append(dnORds(translate(codon)[0], translate(new_codon)[0]))
        if '*' in temp_NorS1:
            temp_NorS.append('*')
        else:
            temp_NorS.append('')
    SNP_gene_temp.addSNP_pair(temp_NorS)

################################################### Main ########################################################
if Analyze:
    # set up sum all seqs
    SNP_gene_all = SNP_gene() # all denovo mutation
    SNP_gene_all.init('all_denovo')
    Output = []
    # process all fasta
    for record in SeqIO.parse(args.i, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        SNP_gene_seq = expecttrunc(record_id,record_seq)
        Output.append(SNP_gene_seq.sumtrunc())
    Output.append(SNP_gene_all.sumtrunc())
    # output
    foutput = open(args.i + 'truncation.simulation.txt', 'w')
    foutput.write('genename\ttruncation_loci\tnon_truncation_loci\ttruncation_ratio\n')
    foutput.write(''.join(Output))
    foutput.close()


################################################### END ########################################################
