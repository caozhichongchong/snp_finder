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
import numpy as np
from scipy.stats import poisson
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-snp",
                      help="input folder of snp summary results",
                      type=str, default='summary/',
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
genemut_file = os.path.join('%s/all.species.all.gene.txt'%(args.snp))
clonal_file = os.path.join('%s/clonal_genelength.txt'%(args.snp))
core_file = args.core
genus_num_cutoff = 3
mutation_start = 0 # needs to have 0 mutations
################################################### Function ########################################################
def load_core(core_file):
    Core = set()
    if core_file!= 'None':
        for lines in open(core_file,'r'):
            if not lines.startswith('record'):
                lines_set = lines.split('\t')
                genus_num = int(lines_set[3])
                if genus_num >= genus_num_cutoff:
                    Core.add(lines_set[0])
    return Core

def change_genename(gene,lineage_new):
    gene = 'C_%s_G_%s' % (gene.split('_')[1], gene.split('_')[-1])
    newgene = '%s__%s' % (lineage_new.split('.donor')[0], gene)
    return newgene

def load_genemut(genemut_file):
    lineage_SNP = dict()
    lineage_SNP.setdefault('all',[0,dict()])
    for lines in open(genemut_file,'r'):
        if not lines.startswith('species'):
            lines_set = lines.split('\t')
            if lines_set[4] != 'NA':
                lineage = lines_set[1]
                genename = lines_set[2]
                lineage_new = lineage.replace('_clustercluster', '_CL').replace('_PB_', '_PaDi_')
                if 'other' not in genename:
                    lineage_SNP.setdefault(lineage,0)
                    SNP = int(lines_set[6])
                    if SNP >= mutation_start:
                        # genes with mutation_start mutations
                        genename = change_genename(genename,lineage_new)
                        if genename not in Core and genename in Cluster:
                            cluster = Cluster[genename]
                            lineage_SNP[lineage] += SNP - mutation_start
                            lineage_SNP['all'][0] += SNP - mutation_start
                            lineage_SNP['all'][1].setdefault(cluster, SNP)
    return lineage_SNP

def load_clonal(clonal_file):
    clonal = dict()
    for lines in open(clonal_file,'r'):
        if not lines.startswith('species'):
            lines_set = lines.split('\n')[0].split('\t')
            lineage = lines_set[1].split('.all')[0]
            nonORFlength = int(lines_set[2])
            ORFlength = int(float(lines_set[3])*float(lines_set[4]))
            clonal.setdefault(lineage,[nonORFlength,ORFlength])
    return clonal

def load_coassembly_list(coassembly_list):
    coassembly_list_set = dict()
    for lines in open(coassembly_list, 'r'):
        database = lines.split('##reference=file:')[1].split('\n')[0] + '.fna'
        lineage = os.path.split(os.path.split(database)[0])[-1]
        coassembly_list_set.setdefault(lineage,database)
    return coassembly_list_set

def load_genes(assembly,gene_length,total_gene_length,lineage_new):
    for record in SeqIO.parse(assembly, 'fasta'):
        record_id = str(record.id)
        genename = change_genename(record_id, lineage_new)
        if genename not in Core and genename in Cluster:
            cluster = Cluster[genename]
            record_seq_len = len(str(record.seq))
            gene_length.setdefault(cluster, record_seq_len)
            total_gene_length += record_seq_len
    return [gene_length,total_gene_length]

def find_clonal(lineage_short):
    if lineage_short in clonal:
        return clonal[lineage_short]
    else:
        return clonal.get(lineage_short.split('_')[0],[0,0])

def find_assemlby(lineage_short):
    if lineage_short in coassembly_list_set:
        return coassembly_list_set[lineage_short]
    else:
        return coassembly_list_set.get(lineage_short.split('_')[0],'')

def pvalue_mutgene(SNP, ORFlength,gene_length,gene_set):
    mut_rate = float(SNP)/ORFlength
    for cluster in gene_length:
        SNP_gene = gene_set.get(cluster,0) - mutation_start # how many mutations beyond mutation_start
        this_gene_length = gene_length[cluster]
        pvalue = poisson.pmf(SNP_gene,mut_rate*this_gene_length)
        allsum_details.append('%s\t%s\t%s\t%s\t%s\t%.2f\n'%(lineage,cluster,pvalue,SNP_gene,this_gene_length,mut_rate*1000))

def load_cluster(cluster_file,Core):
    Cluster = dict()
    if cluster_file!= 'None':
        for lines in open(cluster_file,'r'):
            lines_set = lines.split('\t')
            cluster = lines_set[1]
            gene_name = lines_set[8]
            if gene_name not in Core:
                Cluster.setdefault(gene_name,cluster)
    return Cluster

################################################### Main ########################################################
# load core flexible genes
Core = load_core(core_file)
# load cluster of all genes
Cluster = load_cluster(core_file.split('.species.sum')[0],Core)
# load SNPs and genes with mutations
lineage_SNP = load_genemut(genemut_file)
# load non ORF length
clonal = load_clonal(clonal_file)

# load reference genomes
coassembly_list_set = load_coassembly_list(args.co)
# simulation
allsum_details = []
allsum_details.append('lineage\tgene_name\tpvalue\tSNP_gene\tgene_length\tmut_rate_1kbp\n')
total_gene_length = 0
gene_length = dict()
for lineage in lineage_SNP:
    if lineage != 'all':
        lineage_short = lineage.split('.donor')[0]
        lineage_new = lineage.replace('_clustercluster', '_CL').replace('_PB_', '_PaDi_')
        nonORFlength,ORFlength = find_clonal(lineage_short)
        SNP = lineage_SNP[lineage]
        if nonORFlength == 0:
            print('no clonal infor for %s in %s' % (lineage, clonal_file))
        else:
            # load gene length
            assembly = find_assemlby(lineage_short)
            if assembly == '':
                print('no assembly for %s in %s' % (lineage, args.co))
            else:
                # load gene length for all genes
                gene_length, total_gene_length = load_genes(assembly, gene_length, total_gene_length,lineage_new)
                print('load %s %s SNPs %s length of genes' % (lineage, SNP, total_gene_length))

# compute pvalue for all genes
lineage_set = ['all']
for lineage in lineage_set:
    SNP, gene_set = lineage_SNP[lineage]
    pvalue_mutgene(SNP,total_gene_length,gene_length,gene_set)
    print('finish compute poisson pvalue ')

foutput = open(genemut_file + '.poisson.pvalue.across.txt', 'w')
foutput.write(''.join(allsum_details))
foutput.close()
################################################### END ########################################################
