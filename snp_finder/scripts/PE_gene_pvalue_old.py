# start
import os,glob
from Bio import SeqIO
import argparse
from scipy.stats import poisson
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
allgenemut_file = glob.glob('%s/EsCo*.norecom.txt'%(args.snp))
clonal_file = ('%s/clonal_genelength_new.txt'%(args.snp))
core_file = args.core
genus_num_cutoff = 3 # core as genes shared by at least 3 genera

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
    lineage = os.path.basename(genemut_file).split('.norecom.txt')[0].replace('.all','')
    lineage_SNP.setdefault(lineage, [0, dict()]) # total SNPs on genes, gene SNPs
    for lines in open(genemut_file,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\t')
            if lines_set[7] != 'None':
                # SNPs on a gene
                genename = lines_set[7]
                lineage_SNP[lineage][0] += 1
                lineage_SNP[lineage][1].setdefault(genename, 0)
                lineage_SNP[lineage][1][genename] += 1
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

def load_genes(assembly):
    total_gene_length = 0
    gene_length = dict()
    for record in SeqIO.parse(assembly, 'fasta'):
        record_id = str(record.id)
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
    mut_rate = float(SNP)/ORFlength
    for gene in gene_length:
        SNP_gene = gene_set.get(gene,0)
        this_gene_length = gene_length[gene]
        pvalue = 1-poisson.cdf(SNP_gene,mut_rate*this_gene_length)
        gene = 'C_%s_G_%s' % (gene.split('_')[1], gene.split('_')[-1])
        newgene ='%s__%s'%(lineage_new.split('.donor')[0],gene)
        gene = '%s__%s'%(lineage_new,gene)
        flexible = 'Flexible'
        if newgene in Core:
            # consider flexible genes
            flexible = 'Core'
        allsum_details.append('%s\t%s\t%s\t%s\t%s\t%.5f\t%s\n'%(lineage,gene,pvalue,SNP_gene,this_gene_length,mut_rate*1000,flexible))

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
allsum_details.append('lineage\tgene_name\tpvalue\tSNP_gene\tgene_length\tmut_rate_1kbp\tflexible\n')
for lineage in lineage_SNP:
    lineage_short = lineage.split('.donor')[0]
    print(lineage_short)
    nonORFlength,ORFlength = find_clonal(lineage)
    SNP, gene_set = lineage_SNP[lineage]
    if gene_set!=dict():
        lineage_new = lineage.replace('_clustercluster', '_CL').replace('_PB_', '_PaDi_')
        if ORFlength == 0:
            print('no clonal infor for %s in %s'%(lineage,clonal_file))
        else:
            # load gene length
            assembly = find_assemlby(lineage_short)
            if assembly == '':
                print('no assembly for %s in %s' % (lineage, args.co))
            else:
                # load gene length for all genes
                gene_length,total_gene_length = load_genes(assembly)
                print('process %s %s SNPs %s length of genes' % (lineage, SNP, total_gene_length))
                # compute pvalue for all genes
                pvalue_mutgene(SNP,total_gene_length,gene_length,gene_set)
                print('finish compute poisson pvalue %s' % (lineage))

foutput = open('%s/summary/allgenes.poisson.pvalue.PEonly.txt'%(args.snp), 'w')
foutput.write(''.join(allsum_details))
foutput.close()
################################################### END ########################################################
