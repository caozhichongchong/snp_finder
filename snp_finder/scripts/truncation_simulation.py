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
import random
import numpy as np
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="input fasta",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/WGS/vcf_round2/merge/summary//all.genome.gene.fna',
                      metavar='all.genome.gene.fna')
required.add_argument("-snp",
                      help="input snp freq",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/MG/summary/all.truncsum.withdetails.countalleles.txt',
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
purines=['A','G']
pyrimidines=['C','T']
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
trunc_set = dict()
trunc_set['*']=0
trunc_set['']=1
trunc_set[0]='*'
trunc_set[1]=''
Analyze = True
simulation_time = 100
################################################### new class #########################################################
__metaclass__ = type

class SNP_gene:
    # create a class to store SNP_gene
    'a class to store SNP_gene'
    def init(self, gene, position):
        self.gene = gene
        self.position = position
        self.SNP_pair = {'A-T': set(),
                         'A-C': set(),
                         'G-C': set(),
                         'G-T': set(),
                         'A-G': set(),
                         'G-A': set()}
    def addSNP_pair(self, trunc,pair):
        self.SNP_pair[pair].add(trunc)
    def addSNP_pairposition(self, position,pair):
        self.SNP_pair[pair].add(position)
    def sumtrunc(self,gene_freq):
        temp_line = '%s\t%s\t'%(self.gene,gene_freq)
        temp_line_set = ''
        for i in range(0, simulation_time):
            positionset = set()
            for gene_freq1 in gene_freq:
                pair,pair_freq = gene_freq1
                random_set = set(random.choices(list(self.SNP_pair[pair]), k=pair_freq))
                positionset = positionset.union(random_set) #Random sampling with replacement
            temp_line_set += '%s\t%s\t%s\t%s\n'%(temp_line,i,len(positionset),list(positionset))
            allsum[self.gene].append(len(positionset))
        return temp_line_set

################################################### Function ########################################################

def transitions(REF,ALT):
    if REF in pyrimidines:
        REF = complement[REF]
        ALT = complement[ALT]
    return '%s-%s'%(REF,ALT)

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
    temp_SNP_gene.init(record_name,len(record_seq))
    j = 0
    for i in range(0, (Total - position)):
        codon = record_seq[(i * 3 + position):((i + 1) * 3 + position)]
        try:
            codon_NSratioset = codontable_NSratio[codon]
            for codon_NSratio in codon_NSratioset:
                for pair in codon_NSratio.SNP_pair:
                    if codon_NSratio.SNP_pair[pair] == {'*'}:
                        temp_SNP_gene.addSNP_pairposition(j,pair)
                j += 1
        except KeyError:
            print('no codon found for %s'%(codon))
    return temp_SNP_gene

def expecttruncadd(record_seq,temp_SNP_gene,position=0):
    Total = int(len(record_seq) / 3)
    j = 0
    for i in range(0, (Total - position)):
        codon = record_seq[(i * 3 + position):((i + 1) * 3 + position)]
        try:
            codon_NSratioset = codontable_NSratio[codon]
            for codon_NSratio in codon_NSratioset:
                for pair in codon_NSratio.SNP_pair:
                    if codon_NSratio.SNP_pair[pair] == {'*'}:
                        temp_SNP_gene.addSNP_pairposition(j, pair)
                j += 1
        except KeyError:
            print('no codon found for %s' % (codon))
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
    codontable_NSratio.setdefault(codon, [])
    for position in range(0, 3):
        SNP_gene_temp = SNP_gene()
        SNP_gene_temp.init(codon,position)
        codontable_NSratio[codon].append(SNP_gene_temp)
        REF = codon[position]
        for ALT in Allels_order:
            if ALT != REF:
                new_codon = causeSNP(codon, position+1, ALT,0)
                SNP_pair = transitions(REF, ALT)
                SNP_gene_temp.addSNP_pair(dnORds(translate(codon)[0], translate(new_codon)[0]),
                                          SNP_pair)

################################################### Main ########################################################
if Analyze:
    allsum = dict()
    Output = []
    Gene_freq = dict()
    Cluster = dict()
    # process real data
    for lines in open(args.snp,'r'):
        if not lines.startswith('gene_cluster'):
            lines_set = lines.split('\n')[0].split('\t')
            cluster, Major_ALT,Minor_ALT,freq = lines_set[0:4]
            Gene_freq.setdefault(cluster,[])
            Gene_freq[cluster].append([transitions(Major_ALT, Minor_ALT),int(freq)])
    for lines in open(args.i.replace('.fna','.faa.uc.short'),'r'):
        lines_set = lines.split('\n')[0].split('\t')
        cluster,gene = lines_set[0:2]
        if cluster in Gene_freq:
            Cluster.setdefault(gene,cluster)
    # process all fasta
    cluster_done = dict()
    for record in SeqIO.parse(args.i, 'fasta'):
        record_id = str(record.id)
        #record_id = record_id.split('.donor')[0]+'__' + record_id.split('__')[1]
        record_seq = str(record.seq)
        if record_id in Cluster:#Cluster
            cluster = Cluster[record_id]
            if cluster not in cluster_done:
                allsum.setdefault(cluster, [])
                cluster_done.setdefault(cluster, expecttrunc(cluster, record_seq))
            else:
                SNP_gene_seq = cluster_done[cluster]
                SNP_gene_seq = expecttruncadd(record_seq, SNP_gene_seq)
    for cluster in cluster_done:
        SNP_gene_seq = cluster_done[cluster]
        Output.append(SNP_gene_seq.sumtrunc(Gene_freq[cluster]))
    # output
    foutput = open(args.i + 'truncation.simulation.withmutfreq.txt', 'w')
    foutput.write('gene_cluster\tcluster_freq\tsimulation_time\ttruncation_loci_num\ttruncation_loci\n')
    foutput.write(''.join(Output))
    foutput.close()
    Output = []
    for cluster in allsum:
        trun_loci_set = np.array(allsum[cluster])
        Output.append('%s\t%s\t%s\t%s\t%s\n'%(cluster,
                                              np.percentile(trun_loci_set, 1),
                                              np.percentile(trun_loci_set, 50),np.percentile(trun_loci_set, 99),
                                              Gene_freq[cluster]))
    foutput = open(args.i + 'truncation.simulation.withmutfreq.sum.txt', 'w')
    foutput.write('gene_cluster\ttruncation_loci_num_1\ttruncation_loci_num_50\ttruncation_loci_num_99\tfreqset\n')
    foutput.write(''.join(Output))
    foutput.close()


################################################### END ########################################################
