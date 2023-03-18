# start
# simulate PE
import glob
import os
from Bio import SeqIO
import statistics
import argparse
from scipy.stats import poisson
from statistics import mean
import random
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-snp",
                      help="input folder of snp summary results",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/removerec_SNP/',
                      metavar='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/removerec_SNP/')
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
clonal_file = os.path.join('%s/clonal_genelength_new.txt'%(args.snp))
core_file = args.core
genus_num_cutoff = 3
simulation_round = 1000
pvalue_max = 0.05
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

def load_genemut(genemut_file,lineage_SNP):
    lineage_original = os.path.basename(genemut_file).split('.norecom')[0].replace('.all', ''
                                                                          )
    lineage = lineage_original.replace('BaFr_clustercluster7.donor.amnew', 'BaFr_clustercluster3.donor.am')
    species = lineage.split('_')[0]
    species_lineage.setdefault(species, set())
    species_lineage[species].add(lineage)
    lineage_new = lineage.replace('_clustercluster', '_CL').replace('_PB_', '_PaDi_').split('_lineage')[0]
    lineage_SNP.setdefault(species, [0, dict()]) # No. SNPs, gene dict
    for lines in open(genemut_file,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\t')
            if lines_set[6] != 'None':
                genename = lines_set[6]
                genename = change_genename(genename,lineage_new)
                if genename in Cluster:
                    cluster = Cluster[genename]
                    # No. SNPs per gene cluster
                    lineage_SNP[species][0] += 1
                    lineage_SNP[species][1].setdefault(cluster, [0,set(),0,0])  # No. SNPs + genome set with SNPs, No. N SNPs, No. truncation SNPs
                    lineage_SNP[species][1][cluster][0] += 1
                    lineage_SNP[species][1][cluster][1].add(lineage_original)
                    if lines_set[8] in ['N', 'NN']:
                        # N SNPs
                        lineage_SNP[species][1][cluster][-1] += 1
                        if '*' in lines_set[9]:
                            # truncation
                            lineage_SNP[species][1][cluster][-2] += 1
    return lineage_SNP

def load_clonal(clonal_file):
    clonal = dict()
    for lines in open(clonal_file,'r'):
        if not lines.startswith('cluster'):
            lines_set = lines.split('\n')[0].split('\t')
            species = lines_set[0].split('_')[0]
            ORFlength = int(lines_set[2])
            clonal.setdefault(species,[])
            clonal[species].append(ORFlength)
    return clonal

def load_coassembly_list(coassembly_list):
    coassembly_list_set = dict()
    for lines in open(coassembly_list, 'r'):
        database = lines.split('##reference=file:')[1].split('\n')[0] + '.fna'
        lineage = os.path.split(os.path.split(database)[0])[-1]
        coassembly_list_set.setdefault(lineage,database)
    return coassembly_list_set

def find_assemlby(lineage_short):
    if lineage_short in coassembly_list_set:
        return coassembly_list_set[lineage_short]
    else:
        return coassembly_list_set.get(lineage_short.split('_')[0],'')

def pvalue_mutgene(SNP, ORFlength,gene_set):
    SNP_genome_set_all = set()
    pvalueset = set()
    mut_rate = float(SNP)/ORFlength
    num_lineage_in_species = len(species_lineage[species])
    allgenepvalue = ['gene\tpvalue\n']
    for cluster in Cluster_length:
        if cluster in gene_set:
            this_gene_length = Cluster_length[cluster]*3 # aa to dna
            SNP_gene,lineage_num,truncation_gene,N_SNP_gene = gene_set.get(cluster,[0,set(),[1], 0,0]) # how many mutations
            num_strains_with_mut = len(Cluster_count[cluster][species])/num_lineage_in_species # num unique genes of the same cluster in a species / number of lineages
            # considering the prevalence of this gene among strains
            pvalue = 1-poisson.cdf(SNP_gene - 1,mut_rate*this_gene_length*num_strains_with_mut) # greater than and equal to
            flexible = 'Flexible'
            if cluster in Core_cluster:
                flexible = 'Core'
            allsum_details.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.5f\t%s\t%s\t%s\n'%(species,cluster,pvalue,SNP_gene,N_SNP_gene,truncation_gene,this_gene_length,num_strains_with_mut,mut_rate*1000,flexible,len(lineage_num)))
            allgenepvalue.append('%s\t%s\n' % (cluster, pvalue))
            if SNP_gene > 1 and len(lineage_num) > 1 and pvalue <= pvalue_max:
                # potential PE
                pvalueset.add(pvalue)
                SNP_genome_set_all.add(len(lineage_num))
        else:
            allgenepvalue.append('%s\t%s\n' % (cluster, 1))
    pvalueset = list(pvalueset)
    pvalueset.sort(reverse=True)
    SNP_genome_set_all = list(SNP_genome_set_all)
    SNP_genome_set_all.sort()
    if len(pvalueset) > 0:
        foutput = open('%s/summary/allgenes.poisson.pvalue.across.%s.details.txt' % (args.snp, species), 'w')
        foutput.write(''.join(allgenepvalue))
        foutput.close()
    return [pvalueset, SNP_genome_set_all]

def mutation_sim(SNP, ORFlength,pvalueset,SNP_genome_set_all):
    cluster_num = []
    cluster_length_set = []
    mut_rate = float(SNP) / ORFlength
    num_lineage_in_species = len(species_lineage[species])
    allgenepvalue = ['simulation\tpvalue\n']
    print(len(Cluster_species[species]))
    for cluster in Cluster_species[species]:
        cluster_num.append(cluster)
        cluster_length_set.append(Cluster_length[cluster])
    for i in range(0, simulation_round):
        cluster_mut = random.choices(cluster_num, weights=cluster_length_set, k=SNP)
        allclusters_mut = set(cluster_mut)
        allsim_cluster = dict()
        for clusterID in cluster_num:
            if clusterID in allclusters_mut:
                SNP_cluster = cluster_mut.count(clusterID)
                this_cluster_length = Cluster_length[clusterID]
                num_strains_with_mut = len(Cluster_count[clusterID][species]) / num_lineage_in_species
                pvalue = 1 - poisson.cdf(SNP_cluster - 1, mut_rate * this_cluster_length * num_strains_with_mut) # greater than and equal to
                allgenepvalue.append('%s\t%s\n' % (i, pvalue))
                if SNP_cluster > 1:
                    # FP PE
                    for SNP_cluster_sub in range(2, SNP_cluster):
                        allsim_cluster.setdefault(SNP_cluster_sub, [])
                        allsim_cluster[SNP_cluster_sub].append(pvalue)
            else:
                allgenepvalue.append('%s\t%s\n' % (i, 1))
        for SNP_gene_sub in SNP_genome_set_all:
            for pvalue in pvalueset:
                num_genes_pass_pvalue = len([x for x in allsim_cluster.get(SNP_gene_sub,[]) if x <= pvalue])
                allsum.append('%s\t%s\t%s\t%s\t%s\n'%(species,SNP_gene_sub,pvalue,i,num_genes_pass_pvalue))
    if len(pvalueset) > 0:
        foutput = open('%s/summary/allgenes.poisson.pvalue.across.simulation.%s.details.txt' % (args.snp, species), 'w')
        foutput.write(''.join(allgenepvalue))
        foutput.close()
    return allsum

def load_blastn(allblastnfiles):
    Cluster = dict()
    Cluster_length = dict()
    Cluster_count = dict()
    Cluster_species = dict()
    for blastnfile in allblastnfiles:
        blastnfilename = os.path.basename(blastnfile)
        species = blastnfilename.split('_')[0]
        lineage_new1 = blastnfilename.split('_%s'%(species))[0]
        lineage_new2 = '_%s'%(species) + blastnfilename.split('_%s'%(species))[1].split('.blasn.txt')[0]
        for lines in open(blastnfile,'r'):
            lines_set = lines.split('\t')
            Gene1, Gene2, Identity, Length = lines_set[0:4]
            Gene1 = change_genename(Gene1, lineage_new1)
            Gene2 = change_genename(Gene2, lineage_new2)
            Length = int(Length)
            Cluster.setdefault(Gene1, Gene2)
            Cluster.setdefault(Gene2, Gene2)
            Cluster_length.setdefault(Gene2, Length)
            Cluster_length[Gene2] = max(Cluster_length[Gene2],Length)
            Cluster_count.setdefault(Gene2, dict())
            Cluster_count[Gene2].setdefault(species, set())
            Cluster_count[Gene2][species].add(Gene2)
            Cluster_count[Gene2][species].add(Gene1)
            Cluster_species.setdefault(species, set())
            Cluster_species[species].add(Gene2)
    return [Cluster, Cluster_length, Cluster_count, Cluster_species]

################################################### Main ########################################################
# load core flexible genes
Core = load_core(core_file)
# load blastn results
allblastnfiles = glob.glob('%s/*clustercluster*clustercluster*.txt'%(assembly_folder))
Cluster,Cluster_length,Cluster_count,Cluster_species = load_blastn(allblastnfiles)

# load SNPs and genes with mutations
lineage_SNP = dict()
species_lineage = dict()
for genemut_file in allgenemut_file:
    print('processing %s'%(genemut_file))
    lineage_SNP = load_genemut(genemut_file,lineage_SNP)

# load non ORF length
clonal = load_clonal(clonal_file)

# load reference genomes
coassembly_list_set = load_coassembly_list(args.co)
# simulation
allsum_details = []
allsum_details.append('species\tgene_name\tpvalue\tSNP_gene\tN_gene\ttruncation_gene\tgene_length\tprevalence\tmut_rate_1kbp\tflexible\tlineage_set\n')
total_gene_length = 0
for species in lineage_SNP:
        SNP, gene_set = lineage_SNP[species]
        if SNP > 2: #SNPs on gene > 1
            ORFlength = statistics.mean(clonal[species])
            pvalueset,SNP_genome_set_all = pvalue_mutgene(SNP, ORFlength, gene_set)
            print('finish compute poisson pvalue for species %s with ORFlength %s and SNPs %s PE pvalue set %s SNP genome set %s' % (species,ORFlength,SNP,pvalueset,SNP_genome_set_all))
            allsum = []
            allsum.append('species\tSNP_cutoff\tpvalue_cutoff\tsimulation_round\tFP\n')
            if len(pvalueset) > 0:
            # simulation
                mutation_sim(SNP, ORFlength, pvalueset,SNP_genome_set_all)
            print('finish simulation %s' % (species))
            foutput = open('%s/summary/allgenes.poisson.pvalue.%s.across.simulation.txt' % (args.snp,species), 'w')
            foutput.write(''.join(allsum))
            foutput.close()

foutput = open('%s/summary/allgenes.poisson.pvalue.across.txt'%(args.snp), 'w')
foutput.write(''.join(allsum_details))
foutput.close()

################################################### END ########################################################
