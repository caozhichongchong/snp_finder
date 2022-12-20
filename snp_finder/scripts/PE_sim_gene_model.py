# start
# simulate PE
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
from scipy.stats import poisson

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-snp",
                      help="input folder of snp summary results",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/removerec_SNP/summary/',
                      metavar='summary/')
required.add_argument("-co",
                      help="a list of co-assemblies",
                      type=str, default='all.reference.list',
                      metavar='all.reference.list')
################################################## Definition ########################################################
args = parser.parse_args()
workingdir=os.path.abspath(os.path.dirname(__file__))
# set up path
summary_file = os.path.join('%s/all.lineagednds.txt'%(args.snp))
clonal_file = os.path.join('%s/clonal_genelength_new.txt'%(args.snp))
# set up parameters
simulation_round = 1000

################################################### Function ########################################################
def load_sum(summary_file):
    lineage_SNP = dict()
    for lines in open(summary_file,'r'):
        if not lines.startswith('lineage'):
            lines_set = lines.split('\t')
            lineage = lines_set[0].replace('.all','')
            SNP = int(lines_set[1])+int(lines_set[2]) # SNPs on ORFs
            lineage_SNP.setdefault(lineage,SNP)
    return lineage_SNP

def load_clonal(clonal_file):
    clonal = dict()
    for lines in open(clonal_file,'r'):
        if not lines.startswith('cluster'):
            lines_set = lines.split('\n')[0].split('\t')
            lineage = lines_set[0]
            ORFlength = int(lines_set[2])
            clonal.setdefault(lineage,ORFlength)
    return clonal

def load_coassembly_list(coassembly_list):
    coassembly_list_set = dict()
    for lines in open(coassembly_list, 'r'):
        database = lines.split('##reference=file:')[1].split('\n')[0] + '.fna'
        lineage = os.path.split(os.path.split(database)[0])[-1]
        coassembly_list_set.setdefault(lineage,database)
    return coassembly_list_set

def load_genes(assembly):
    gene_num = []
    gene_length = []
    for record in SeqIO.parse(assembly, 'fasta'):
        gene = str(record.id)
        gene = 'C_%s_G_%s' % (gene.split('_')[1], gene.split('_')[-1])
        gene = '%s__%s' % (lineage_new, gene)
        record_seq_len = len(str(record.seq))
        gene_num.append(gene)
        gene_length.append(record_seq_len)
    return [gene_num,gene_length]

def find_clonal(lineage_short):
    if lineage_short in clonal:
        return clonal[lineage_short]
    else:
        return clonal.get(lineage_short.split('_')[0],0)

def find_assemlby(lineage_short):
    if lineage_short in coassembly_list_set:
        return coassembly_list_set[lineage_short]
    else:
        return coassembly_list_set.get(lineage_short.split('_')[0],'')

def pvalue_mutgene(mut_rate,SNP_gene,this_gene_length):
    # pvalue to be >= observed value
    pvalue = 1-poisson.cdf(SNP_gene, mut_rate * this_gene_length)
    return pvalue

def mutation_sim():
    No_SNP = lineage_SNP[lineage]
    allsim_gene = dict()
    allsim_gene_mut = dict()
    allsim_gene_empty = dict()
    mut_rate = No_SNP / ORFlength
    for i in range(0,simulation_round):
        # simulate No_SNP on each gene
        gene_mut = random.choices(gene_num, weights = gene_length,k=No_SNP)
        # count SNPs on each gene
        for geneID in gene_num:
            gene_length_geneID = gene_length[gene_num.index(geneID)]
            No_mut = 0
            if geneID in gene_mut:
                No_mut = gene_mut.count(geneID)
                pvalue = pvalue_mutgene(mut_rate, No_mut, gene_length_geneID)
            elif geneID not in allsim_gene_empty:
                pvalue = pvalue_mutgene(mut_rate, No_mut, gene_length_geneID)
                allsim_gene_empty.setdefault(geneID,pvalue)
            else:
                pvalue = allsim_gene_empty[geneID]
            allsim_gene.setdefault(geneID,[])
            allsim_gene[geneID].append(pvalue)
            allsim_gene_mut.setdefault(geneID, [])
            allsim_gene_mut[geneID].append(No_mut)
            # record pvalue for each gene for each simulation
    # sum up all simulations
    for gene in allsim_gene:
        gene_length_geneID = gene_length[gene_num.index(gene)]
        pvalueset = allsim_gene.get(gene,[])
        No_mut_set = allsim_gene_mut.get(gene,[])
        pvalueset += [pvalue_mutgene(mut_rate, 0, gene_length_geneID)]*(simulation_round-len(pvalueset))
        No_mut_set += [0] * (simulation_round - len(No_mut_set))
        allsum.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (lineage,gene,statistics.mean(No_mut_set),
                                                np.percentile(pvalueset, 5),
                                                np.percentile(pvalueset, 50), np.percentile(pvalueset, 95)))

################################################### Main ########################################################
# load SNPs
lineage_SNP = load_sum(summary_file)

# load non ORF length
clonal = load_clonal(clonal_file)

# load reference genomes
coassembly_list_set = load_coassembly_list(args.co)
# simulation
allsum = []
allsum.append('lineage\tGene_ID\tnum_mutation\tlow\tmedium\thigh\n')
foutput = open(summary_file + '.simulation.genemodel.geneswithmut.sum.txt', 'w')
for lineage in lineage_SNP:
    lineage_short = lineage.split('.donor')[0]
    lineage_new = lineage.replace('_clustercluster', '_CL').replace('_PB_', '_PaDi_')
    ORFlength = find_clonal(lineage)
    print('process %s  ORF %sbp' % (lineage, ORFlength))
    if ORFlength == 0:
        print('no clonal infor for %s in %s' % (lineage, clonal_file))
    else:
        # load gene length
        assembly = find_assemlby(lineage_short)
        if assembly == '':
            print('no assembly for %s in %s' % (lineage, args.co))
        else:
            gene_num, gene_length = load_genes(assembly)
            # simulation
            mutation_sim()
            print('finish simulation %s' % (lineage))
    foutput.write(''.join(allsum))
    allsum = []

foutput.close()
################################################### END ########################################################
