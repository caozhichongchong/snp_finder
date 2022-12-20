# start
import os,glob
from Bio import SeqIO
import argparse
from scipy.stats import poisson
from statistics import mean
import pandas as pd
import numpy as np
import random

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-snp",
                      help="input folder of snp summary results",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/removerec_SNP/',
                      metavar='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/removerec_SNP/')
required.add_argument("-MW",
                      help="path of folders of all vcfs",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/MV_cov',
                      metavar='input/')

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
allgenemut_file = glob.glob('%s/*.norecom.gooddepth.txt'%(args.snp))
clonal_file = ('%s/clonal_genelength_new.txt'%(args.snp))
core_file = args.core
genus_num_cutoff = 3 # core as genes shared by at least 3 genera
depth_cutoff = 6
min_good_alignment_samples = .4 # % of samples with bad alignment
simulation_round = 1000
pvalue_max = 0.01

try:
    os.mkdir(args.snp + '/summary')
except IOError:
    pass

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

def load_genemut(genemut_file,lineage_SNP):
    lineage = os.path.basename(genemut_file).split('.norecom')[0].replace('.all','')
    lineage_SNP.setdefault(lineage, [0, dict()]) # total SNPs on genes, gene SNPs
    for lines in open(genemut_file,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\t')
            if lines_set[6] != 'None':
                # SNPs on a gene
                genename = lines_set[6]
                lineage_SNP[lineage][0] += 1
                lineage_SNP[lineage][1].setdefault(genename, [0,set(),[],0,0])# No. SNPs + genome set with SNPs, No. isolates, No. N SNPs, No. truncation SNPs
                lineage_SNP[lineage][1][genename][0] += 1
                lineage_SNP[lineage][1][genename][1].add(lines_set[5])
                alleles = lines_set[4]
                alleles_total = len(alleles)
                lineage_SNP[lineage][1][genename][2].append(
                    (len(alleles) - len([x for x in alleles if x == '-'])) / alleles_total)
                if lines_set[8] in ['N','NN']:
                    # N SNPs
                    lineage_SNP[lineage][1][genename][-1] += 1
                    if '*' in lines_set[9]:
                        # truncation
                        lineage_SNP[lineage][1][genename][-2] += 1
    return lineage_SNP

def load_clonal(clonal_file):
    clonal = dict()
    for lines in open(clonal_file,'r'):
        if not lines.startswith('cluster'):
            lines_set = lines.split('\n')[0].split('\t')
            lineage = lines_set[0]
            nonORFlength = int(lines_set[1])-int(lines_set[2])
            ORFlength = int(lines_set[2])
            clonal.setdefault(lineage,[nonORFlength,ORFlength])
    return clonal

def load_coassembly_list(coassembly_list):
    coassembly_list_set = dict()
    for lines in open(coassembly_list, 'r'):
        database = lines.split('##reference=file:')[1].split('\n')[0] + '.fna'
        lineage = os.path.split(os.path.split(database)[0])[-1]
        coassembly_list_set.setdefault(lineage,database)
    return coassembly_list_set

def load_genes(assembly,allSNPscov):
    total_gene_length = 0
    gene_length = dict()
    for record in SeqIO.parse(assembly, 'fasta'):
        record_id = str(record.id)
        CHR = '_'.join(record_id.split('_')[:-1])
        record_des = str(record.description)
        startPOS = int(record_des.split(' # ')[1])-1000
        endPOS = int(record_des.split(' # ')[2])+1000
        if allSNPscov.loc[(allSNPscov['CHR']==CHR) & (allSNPscov['POS']<=endPOS) & (allSNPscov['POS']>=startPOS),:].shape[0] > 0:
            record_seq_len = len(str(record.seq))
            gene_length.setdefault(record_id, record_seq_len)
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
    SNP_genome_set_all = set()
    pvalueset = set()
    mut_rate = float(SNP)/ORFlength
    for gene in gene_set:
        if gene in gene_length:
            SNP_gene,SNP_genome_set,num_strains_with_mut,truncation_gene,N_SNP_gene = gene_set.get(gene,[0,set(),[1],0,0])
            num_strains_with_mut = mean(num_strains_with_mut)
            this_gene_length = gene_length[gene]
            # considering the prevalence of this gene among strains
            pvalue = 1-poisson.cdf(SNP_gene,mut_rate*this_gene_length*num_strains_with_mut) + \
                     poisson.pmf(SNP_gene,mut_rate*this_gene_length*num_strains_with_mut)# greater than and equal to
            gene = 'C_%s_G_%s' % (gene.split('_')[1], gene.split('_')[-1])
            newgene ='%s__%s'%(lineage_new,gene)
            gene = '%s__%s'%(lineage_new,gene)
            flexible = 'Flexible'
            if newgene in Core:
                # consider flexible genes
                flexible = 'Core'
            allsum_details.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%.5f\t%s\t%s\n'%(lineage,gene,pvalue,SNP_gene,N_SNP_gene,truncation_gene,this_gene_length,num_strains_with_mut,mut_rate*1000,flexible,len(SNP_genome_set)))
            if SNP_gene > 1 and len(SNP_genome_set) > 1 and pvalue <= pvalue_max:
                # potential PE
                pvalueset.add(pvalue)
                SNP_genome_set_all.add(len(SNP_genome_set))
        else:
            print('gene %s not qualified in lineage %s'%(gene,lineage))
    pvalueset = list(pvalueset)
    pvalueset.sort(reverse=True)
    SNP_genome_set_all = list(SNP_genome_set_all)
    SNP_genome_set_all.sort()
    return [pvalueset,SNP_genome_set_all]

def mutation_sim(SNP, ORFlength,gene_length,pvalueset,SNP_genome_set_all):
    gene_num = []
    gene_length_set = []
    num_strains_with_mut = 1 - min_good_alignment_samples# lower bound
    mut_rate = float(SNP)/ORFlength
    for gene,this_gene_length in gene_length.items():
        gene_num.append(gene)
        gene_length_set.append(this_gene_length)
    for i in range(0,simulation_round):
        gene_mut = random.choices(gene_num, weights=gene_length_set, k=SNP)
        allgenes_mut = set(gene_mut)
        allsim_gene = dict()
        for geneID in allgenes_mut:
            SNP_gene = gene_mut.count(geneID)
            if SNP_gene > 1:
                # FP PE
                this_gene_length = gene_length[geneID]
                pvalue = 1 - poisson.cdf(SNP_gene, mut_rate * this_gene_length * num_strains_with_mut) + \
                         poisson.pmf(SNP_gene,
                                     mut_rate * this_gene_length * num_strains_with_mut)  # greater than and equal to
                for SNP_gene_sub in range(2,SNP_gene):
                    allsim_gene.setdefault(SNP_gene_sub,[])
                    allsim_gene[SNP_gene_sub].append(pvalue)
        for SNP_gene_sub in SNP_genome_set_all:
            for pvalue in pvalueset:
                num_genes_pass_pvalue = len([x for x in allsim_gene.get(SNP_gene_sub,[]) if x <= pvalue])
                allsum.append('%s\t%s\t%s\t%s\t%s\n'%(lineage,SNP_gene_sub,pvalue,i,num_genes_pass_pvalue))
    return allsum

def depth_screening(allSNPscov):
    allSNPscov = allSNPscov[allSNPscov['POS'] > 100]  # not considering the depth of the first POS of each contig
    allSNPscov.index = range(0, allSNPscov.shape[0])
    # avg depth of non zero
    total_genomes = allSNPscov.shape[1] - 3
    for i in allSNPscov.index:
        non_zero_depth = [x for x in allSNPscov.iloc[i, 3:] if x > 0]
        if len(non_zero_depth) >= min_good_alignment_samples * total_genomes:
            # good prevalence
            allSNPscov.loc[i, 'avg_depth'] = np.mean(non_zero_depth)
        else:
            allSNPscov.loc[i, 'avg_depth'] = 0
    allSNPscov = allSNPscov[allSNPscov['avg_depth'] >= depth_cutoff] # good depth
    return allSNPscov

################################################### Main ########################################################
# load core flexible genes
Core = load_core(core_file)

# load SNPs and genes with mutations
lineage_SNP = dict()
for genemut_file in allgenemut_file:
    print('processing %s'%(genemut_file))
    lineage_SNP = load_genemut(genemut_file,lineage_SNP)

# load non ORF length
clonal = load_clonal(clonal_file)

# load reference genomes
coassembly_list_set = load_coassembly_list(args.co)
# simulation
allsum_details = []
allsum_details.append('lineage\tgene_name\tpvalue\tSNP_gene\tN_gene\ttruncation_gene\tgene_length\tprevalence\tmut_rate_1kbp\tflexible\tgenome_set\n')
allsum = []
allsum.append('lineage\tSNP_cutoff\tpvalue_cutoff\tsimulation_round\tFP\n')
for lineage in lineage_SNP:
    allSNPscov = pd.read_csv(os.path.join(args.MW, '%s.raw.vcf.filtered.cov.MW.txt') % (
        lineage.replace('.donor', '.all.donor')).split('_lineage')[0]
                             , sep='\t')
    allSNPscov = depth_screening(allSNPscov)
    lineage_short = lineage.split('.donor')[0]
    nonORFlength,ORFlength = find_clonal(lineage.split('_lineage')[0]) # callable genome length
    SNP, gene_set = lineage_SNP[lineage]
    if SNP > 1 and gene_set!=dict(): #SNPs on gene > 1
        lineage_new = lineage.replace('_clustercluster', '_CL').replace('_PB_', '_PaDi_').split('.donor')[0]
        if ORFlength == 0:
            print('no clonal infor for %s %s in %s'%(lineage,lineage.split('_lineage')[0],clonal_file))
        else:
            # load gene length
            assembly = find_assemlby(lineage_short)
            if assembly == '':
                print('no assembly for %s in %s' % (lineage, args.co))
            else:
                # load gene length for all genes
                gene_length,total_gene_length = load_genes(assembly,allSNPscov) #total_gene_length not used
                print('process %s %s SNPs %s ORF length' % (lineage, SNP, ORFlength))
                # compute pvalue for all genes
                pvalueset, SNP_genome_set_all = pvalue_mutgene(SNP,ORFlength,gene_length,gene_set)
                print('finish compute poisson pvalue %s PE pvalue set %s SNP genome set %s' % (lineage,pvalueset,SNP_genome_set_all))
                if len(pvalueset) > 0:
                    # simulation
                    mutation_sim(SNP, ORFlength,gene_length,pvalueset,SNP_genome_set_all)
                print('finish simulation %s' % (lineage))

foutput = open('%s/summary/allgenes.poisson.pvalue.within.txt'%(args.snp), 'w')
foutput.write(''.join(allsum_details))
foutput.close()
foutput = open('%s/summary/allgenes.poisson.pvalue.within.simulation.txt'%(args.snp), 'w')
foutput.write(''.join(allsum))
foutput.close()
################################################### END ########################################################
