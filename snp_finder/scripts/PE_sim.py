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
optional.add_argument("-cutoff",
                      help="a file of cutoff of how many SNPs on a gene for each clonal population to call parallel evolution",
                      type=str, default='None',
                      metavar='total_SNP_cutoff.txt')
################################################## Definition ########################################################
args = parser.parse_args()
workingdir=os.path.abspath(os.path.dirname(__file__))
# set up path
summary_file = os.path.join('%s/all.species.lineage.dmrca.txt'%(args.snp))
clonal_file = os.path.join('%s/clonal_genelength.txt'%(args.snp))
# set up parameters
simulation_round = 101
Min_SNP_highselect_cutoff = 1/3000

################################################### Function ########################################################
def load_sum(summary_file):
    lineage_SNP = dict()
    for lines in open(summary_file,'r'):
        if not lines.startswith('X.donor_species'):
            lines_set = lines.split('\t')
            if lines_set[4] != 'NA':
                lineage = lines_set[0]
                SNP = int(lines_set[4])
                PE_SNP = int(lines_set[14])
                lineage_SNP.setdefault(lineage,[SNP,PE_SNP])
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

def load_genes(assembly):
    gene_num = []
    gene_length = []
    i = 0
    for record in SeqIO.parse(assembly, 'fasta'):
        record_seq_len = len(str(record.seq))
        gene_num.append(i)
        gene_length.append(record_seq_len)
        i+=1
    return [gene_num,gene_length]

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

def mutation_sim_orf(No_SNP):
    genome_mut = random.choices(['nonORF','ORF'], weights=[nonORFlength, ORFlength], k=No_SNP)
    return genome_mut.count('ORF')

def mutation_sim():
    No_SNP,NO_PE_SNP = lineage_SNP[lineage]
    allsum_dict = []
    for i in range(0,simulation_round):
        PE_num = 0
        gene_mut_sum = dict()
        # simulate No. SNPs on ORF region
        NO_SNP_ORF = mutation_sim_orf(No_SNP)
        # simulate NO_SNP_ORF on each gene
        gene_mut = random.choices(gene_num, weights = gene_length,k=NO_SNP_ORF)
        # count SNPs on each gene
        gene_mut_set = set(gene_mut)
        for geneID in gene_mut_set:
            No_mut = gene_mut.count(geneID)
            gene_length_geneID = gene_length[gene_num.index(geneID)]
            if No_mut/gene_length_geneID < Min_SNP_highselect_cutoff:
                # normalize against gene length
                No_mut = PE_cutoff - 1 # should not be a PE gene
            gene_mut_sum.setdefault(No_mut,0)
            gene_mut_sum[No_mut]+=1
        for No_mut in gene_mut_sum:
            allsum_details.append('%s\t%s\t%s\t%s\n'%(lineage,i,No_mut,gene_mut_sum[No_mut]))
            if No_mut >= PE_cutoff:
                PE_num += gene_mut_sum[No_mut]*No_mut # total number of PE SNPs
        allsum_dict.append(PE_num)
    # sum up all simulations
    realPE = NO_PE_SNP
    simPE = allsum_dict
    simPE.append(realPE)
    simPE.sort()
    pvalue = 1 - (simPE.index(realPE) + 1) / (simulation_round + 1)
    allsum.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (lineage,PE_cutoff,
                                            int(np.percentile(allsum_dict, 5)),
                                            int(np.percentile(allsum_dict, 50)), int(np.percentile(allsum_dict, 95)),pvalue))

################################################### Main ########################################################
# load cutoff
Donor_species = dict()
if args.cutoff != 'None':
    for lines in open(args.cutoff):
        lines_set = lines.replace('\n', '').replace('\r', '').split('\t')
        donor_species = lines_set[0]
        cutoff = int(lines_set[1])
        Donor_species.setdefault(donor_species,cutoff)

# load SNPs
lineage_SNP = load_sum(summary_file)

# load non ORF length
clonal = load_clonal(clonal_file)

# load reference genomes
coassembly_list_set = load_coassembly_list(args.co)
# simulation
allsum = []
allsum.append('lineage\tPE_SNP_cutoff\tlow\tmedium\thigh\n')
allsum_details = []
allsum_details.append('lineage\tsim_round\tNo_SNPs_per_gene\tNo_genes\tpvalue\n')
for lineage in lineage_SNP:
    lineage_short = lineage.split('.donor')[0]
    print(lineage_short)
    nonORFlength,ORFlength = find_clonal(lineage_short)
    PE_cutoff = Donor_species.get(lineage,2)
    print('process %s nonORF %sbp ORF %sbp, PE cutoff %s'%(lineage,nonORFlength,ORFlength,PE_cutoff))
    if nonORFlength == 0:
        print('no clonal infor for %s in %s'%(lineage,clonal_file))
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

foutput = open(summary_file + '.simulation.details.txt', 'w')
foutput.write(''.join(allsum_details))
foutput.close()

foutput = open(summary_file + '.simulation.sum.txt', 'w')
foutput.write(''.join(allsum))
foutput.close()
################################################### END ########################################################
