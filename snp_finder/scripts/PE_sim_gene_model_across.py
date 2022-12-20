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
required.add_argument("-core",
                      help="core flexible genes",
                      type=str, default='None',
                      metavar='all.genome.gene.faa.uc.species.sum')

################################################## Definition ########################################################
args = parser.parse_args()
workingdir=os.path.abspath(os.path.dirname(__file__))
# set up path
summary_file = os.path.join('%s/all.speciesdnds.noparients.txt'%(args.snp))
clonal_file = os.path.join('%s/clonal_genelength_new.txt'%(args.snp))
core_file = args.core
# set up parameters
simulation_round = 1000

################################################### Function ########################################################
def load_sum(summary_file):
    lineage_SNP = dict()
    for lines in open(summary_file,'r'):
        if not lines.startswith('lineage'):
            lines_set = lines.split('\t')
            lineage = lines_set[0]
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

def load_cluster(cluster_file,cluster_fasta):
    Cluster_length = dict()
    Species_gene_length = dict()
    for record in SeqIO.parse(cluster_fasta, 'fasta'):
        Cluster_length.setdefault(str(record.id),len(str(record.seq)))
    if cluster_file!= 'None':
        for lines in open(cluster_file,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            cluster = lines_set[9]
            gene_name = lines_set[8]
            species = gene_name.split('_')[0]
            if cluster == '*':
                cluster = gene_name
            Species_gene_length.setdefault(species, set())
            Species_gene_length[species].add(cluster)
    return [Cluster_length,Species_gene_length]

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
    No_SNP = lineage_SNP[species]
    allsim_gene = dict()
    allsim_gene_mut = dict()
    allsim_gene_empty = dict()
    mut_rate = No_SNP / ORFlength
    gene_num = Species_gene_length[species]
    gene_length = [Cluster_length[x] for x in gene_num]
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
        allsum.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (species,gene,statistics.mean(No_mut_set),
                                                np.percentile(pvalueset, 5),
                                                np.percentile(pvalueset, 50), np.percentile(pvalueset, 95)))

################################################### Main ########################################################
# load cluster of all genes
Cluster_length,Species_gene_length = load_cluster(core_file.split('.species.sum')[0],
                                                   core_file.split('.uc.species.sum')[0]+'.cluster.aa')

# load SNPs
lineage_SNP = load_sum(summary_file)

# load non ORF length
clonal = load_clonal(clonal_file)

# load reference genomes
coassembly_list_set = load_coassembly_list(args.co)
# simulation
allsum = []
allsum.append('species\tGene_ID\tnum_mutation\tlow\tmedium\thigh\n')
foutput = open(summary_file + '.simulation.genemodelacross.geneswithmut.sum.txt', 'w')
for species in lineage_SNP:
    SNP = lineage_SNP[species]
    if SNP > 2:
        ORFlength = statistics.mean(clonal[species])
        # simulation
        mutation_sim()
        print('finish simulation %s' % (species))
        foutput.write(''.join(allsum))
        allsum = []

foutput.close()
################################################### END ########################################################
