# start
# simulate PE
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import argparse
import random
import numpy as np
import math
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-faa",
                      help="input fasta file of genes with snps",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.denovo.gene.faa',
                      metavar='all.denovo.gene.faa')
required.add_argument("-all",
                      help="input annotation of genes with snps",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.denovo.gene.faa.cluster.aa.all.eggnog.sum.species.sum',
                      metavar='all.denovo.gene.faa.cluster.aa.all.eggnog.sum.species.sum')
required.add_argument("-pe",
                      help="input annotation of PE",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.selected.gene.faa.cluster.aa.all.eggnog.sum.species.sum',
                      metavar='all.selected.gene.faa.cluster.aa.all.eggnog.sum.species.sum')
required.add_argument("-species",
                      help="species name",
                      type=str, default='None',
                      metavar='BA, BL, PaDi')

################################################## Definition ########################################################
args = parser.parse_args()
workingdir=os.path.abspath(os.path.dirname(__file__))
# set up parameters
simulation_round = 1000
summary_file = os.path.join(os.path.split(args.pe)[0],'all.PE.pathway.sim')
if args.species != 'None':
    summary_file += '.%s' % (args.species)
checkflexible = False
################################################### Function ########################################################

def load_anno(annnofile,CL_set):
    Anno = dict()
    CL_set_setup = False
    if CL_set == set():
        CL_set_setup = True
    for lines in open(annnofile,'r'):
        if not lines.startswith('cluster'):
            lines_set = lines.split('\n')[0].split('\t')
            if not checkflexible or lines_set[-2] == 'False' or lines_set[-1] == 'False':
                # flexible genes
                genename = lines_set[4]
                if args.species == 'None' or args.species in genename:
                    CL = genename.split('__')[0]
                    COG2 = lines_set[9]
                    if CL_set_setup:
                        CL_set.add(CL)
                        Anno.setdefault(genename, COG2)
                    elif CL in CL_set:
                        Anno.setdefault(genename, COG2)
    return [CL_set,Anno]

def load_length(fasta,Anno):
    genename_list = []
    gene_length = []
    for record in SeqIO.parse(fasta, 'fasta'):
        record_id = str(record.id)
        if record_id in Anno:
            record_seq_len = len(str(record.seq))
            # square root of length
            gene_length.append(math.sqrt(int(record_seq_len)))
            genename_list.append(record_id)
    return [genename_list,gene_length]

def simulate_PE(Anno_PE,Anno_all,genename_list,gene_length):
    num_PE = len(Anno_PE)
    gene_length=gene_length/np.sum(gene_length)
    allsimdetails = []
    allsimdetails.append('sim\tgene_name\tCOG2\n')
    allsimsum = []
    allsimsum.append('sim\tCOG2\tnum_PE\n')
    allsim = dict()
    allPE = dict()
    allsum = []
    allsum.append('COG2\texpectedPE\trealPE\tpvalue\n')
    for i in range(0, simulation_round):
        # random pick num_PE genes with weights = square root of gene_length
        gene_PE = np.random.choice(genename_list,num_PE,p = gene_length,replace=False)
        simsum = dict()
        # find COG for each gene
        for gene in gene_PE:
            COG2 = Anno_all[gene]
            allsimdetails.append('%s\t%s\t%s\n'%(i,gene,COG2))
            simsum.setdefault(COG2,0)
            simsum[COG2] += 1
        # sum up the COG
        for COG2 in simsum:
            allsimsum.append('%s\t%s\t%s\n' % (i, COG2, simsum[COG2]))
            allsim.setdefault(COG2,[])
            allsim[COG2].append(simsum[COG2])
    # compute real number of PE in pathways
    for gene in Anno_PE:
        COG2 = Anno_PE[gene]
        allPE.setdefault(COG2,0)
        allPE[COG2] += 1
    # compute pvalue for real number
    for COG2 in allPE:
        realPE = allPE[COG2]
        simPE = allsim[COG2]
        meanPE = statistics.mean(simPE)
        simPE.append(realPE)
        simPE.sort()
        pvalue = 1 - (simPE.index(realPE)+1)/(simulation_round+1)
        allsum.append('%s\t%s\t%s\t%s\n'%(COG2,meanPE,realPE,pvalue))
    # output to files
    foutput = open(summary_file + '.details.txt', 'w')
    foutput.write(''.join(allsimdetails))
    foutput.close()
    foutput = open(summary_file + '.simsum.txt', 'w')
    foutput.write(''.join(allsimsum))
    foutput.close()
    foutput = open(summary_file + '.pvaluesum.txt', 'w')
    foutput.write(''.join(allsum))
    foutput.close()

################################################### Main ########################################################
# load annotation
CL_set_PE = set()
# load PE genes
CL_set_PE,Anno_PE = load_anno(args.pe,CL_set_PE)
# load all genes with SNPs in PE lineages
CL_set_PE,Anno_all = load_anno(args.all,CL_set_PE)
# load all gene length
genename_list,gene_length=load_length(args.faa,Anno_all)

# simulation
simulate_PE(Anno_PE,Anno_all,genename_list,gene_length)

################################################### END ########################################################
