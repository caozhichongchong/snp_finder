# start
# sum PE mapping length + all gene mapping results
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
from statistics import stdev
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-vcf",
                      help="path of vcf file to process",
                      type=str, default='None',
                      metavar='input/')
required.add_argument("-snp",
                      help="path of snp files",
                      type=str, default='None',
                      metavar='/scratch/users/anniz44/genomes/donor_species/WGS/vcf_round1/merge/')
required.add_argument("-mfq",
                      help="file extension of metagenomes fastq #1 files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')
# optional output setup
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='snp_output/',
                      metavar='snp_output/')

################################################## Definition ########################################################
args = parser.parse_args()
# setup path all clonal population
output_dir = args.o + '/MG/bwa/'
try:
    os.mkdir(output_dir + '/finished')
except IOError:
    pass

try:
    os.mkdir(args.o + '/MG/truncgenesum/')
except IOError:
    pass

trunc_output  = args.o + '/MG/summary/all.lineage.PESNP.sum'

trunc_output_gene  = args.o + '/MG/summary/all.lineage.genemapping.sum'
################################################### Set up ########################################################
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

mapping_qual = 42
################################################### Function ########################################################
# set up functions
def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def contig_to_gene(CHR, POS):
    all_genes = Mapping_loci.get(CHR,[])
    Reverse_chr = 0
    for a_gene in all_genes:
        POS1, POS2, GENE = a_gene
        if POS >= POS1 and POS <= POS2:
            Ref_seq_chr = Ref_seq.get(GENE, 'None')
            Gene_length = len(Ref_seq_chr)
            if GENE in Reverse:  # reversed
                POS_gene = Gene_length-(int(POS-POS1))
                Reverse_chr = 1
            else:
                POS_gene = int(POS-POS1)+1
            codon_start = POS_gene - 1 - int((POS_gene - 1) % 3)
            return [GENE,POS_gene,codon_start,Ref_seq_chr,Reverse_chr]
    return []

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def get_snp(lines_set,No_PESNP,allgenelist):
    CHR,POS = lines_set[0:2]
    POS = int(POS)
    # check truncation
    gene_info = contig_to_gene(CHR, POS)
    if gene_info != []:
        Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
        withinHS = 'False'
        if Chr_gene in HS_gene_within_lineage:
            withinHS = 'True'
            No_PESNP += 1
        allgenelist.add('%s\t%s' % (Chr_gene, withinHS))
    return No_PESNP

def loadHS(sum_file):
    HS_gene = dict()
    for lines in open(sum_file, 'r'):
        if not lines.startswith("#donor_species"):
            lines_set = lines.split('\n')[0].split('\t')
            if lines_set[-1] == 'True':
                lineage, gene = lines_set[0:2]
                lineage = lineage.split('.donor')[0]
                HS_gene.setdefault(lineage,set())
                HS_gene[lineage].add(gene)
    return HS_gene

def loadmultiple(sum_file):
    HS_gene = dict()
    for lines in open(sum_file, 'r'):
        if not lines.startswith("Refgenome"):
            lines_set = lines.split('\n')[0].split('\t')
            lineage, gene,POS = lines_set[0:3]
            HS_gene.setdefault(lineage, set())
            HS_gene[lineage].add('%s\t%s'%(gene,POS))
    return HS_gene

def loaddatabase(database):
    # load database seq
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping_loci,Reverse]

def loadcluster(clusterfile):
    cluster_gene = dict()
    for lines in open(clusterfile, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        genecluster, genename = lines_set[0:2]
        cluster_gene.setdefault(genename, genecluster)
    return cluster_gene


def loadchangename(clusterfile):
    changegene = dict()
    for lines in open(clusterfile, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        lineage, lineageshort,genename,newname = lines_set[0:4]
        changegene.setdefault('%s\t%s'%(lineage,genename), newname)
    return changegene

def sum_snp(No_PESNP, allgenelist,samplename,donorspecies):
    allsum = []
    cluster_count = set()
    cluster_count_PE = set()
    gene_name_PE = set()
    gene_name_all = set()
    for a_gene_name in allgenelist:
        gene_name,withinHS = a_gene_name.split('\t')
        newgenename = changegene['%s\t%s'%(donorspecies,gene_name)]
        genecluster = cluster_gene[newgenename]
        if withinHS == 'True':
            cluster_count_PE.add(genecluster)
            gene_name_PE.add(gene_name)
        else:
            cluster_count.add(genecluster)
            gene_name_all.add(gene_name)
    allsum.append('%s\t%s\t%s\t%s\t%s\t%s\n'%(donorspecies,samplename,len(cluster_count_PE),len(cluster_count),len(gene_name_PE),len(gene_name_all)))
    f1 = open(trunc_output_gene, 'a')
    f1.write(''.join(allsum))
    f1.close()
    f1 = open(trunc_output, 'a')
    f1.write('%s\t%s\t%s\n'%(donorspecies,samplename,No_PESNP))
    f1.close()

################################################### Main ########################################################
# load HS genes
HS_gene_within = loadHS('%s/summary/all.species.txt'%(args.snp))

# load change gene info
changegene = loadchangename('%s/summary/all.genome.gene.faa.changename.txt'%(args.snp))
# load cluster info
cluster_gene = loadcluster('%s/summary/all.genome.gene.faa.uc.short'%(args.snp))
# run SNP truncation prediction
if args.vcf == 'None':
    all_vcf_file=glob.glob(os.path.join(output_dir,'*.raw.vcf'))
else:
    all_vcf_file = glob.glob(os.path.join(output_dir, '%s'%(args.vcf)))

for vcf_file in all_vcf_file:
    samplename = os.path.split(vcf_file)[-1].split(args.mfq)[0]
    donorspecies = os.path.split(vcf_file)[-1].split(args.mfq + '.')[1].split('.raw.vcf')[0]
    print('processing', samplename, donorspecies)
    assembly_file = glob.glob(
        '%s/../co-assembly/withPE/%s/%s.all.spades1.fasta.noHM.fasta.fna' % (args.snp, donorspecies, donorspecies)) + \
                    glob.glob('%s/../co-assembly/withPE_new/%s/%s.all.spades1.fasta.noHM.fasta.fna' % (
                        args.snp, donorspecies, donorspecies)) + \
                    glob.glob('%s/../co-assembly/%s/%s.all.spades1.fasta.noHM.fasta.fna' % (
                        args.snp, donorspecies, donorspecies))
    Ref_seq, Mapping_loci, Reverse = loaddatabase(assembly_file[0])
    HS_gene_within_lineage = HS_gene_within.get(donorspecies, set())
    No_PESNP = 0
    allgenelist = set()
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#") and 'INDEL' not in lines:
            lines_set = lines.split('\n')[0].split('\t')
            quality = float(lines_set[5])
            if quality >= mapping_qual:
                # check mapping quality and no >2 genotypes
                # map snp to snps found in a lineage
                No_PESNP = get_snp(lines_set,No_PESNP,allgenelist)
    # output snps found in a lineage
    sum_snp(No_PESNP, allgenelist,samplename,donorspecies)
    outputfolderzip = output_dir + '/finished/%s'%(donorspecies)
    try:
        os.mkdir(outputfolderzip)
    except IOError:
        pass
    try:
        f1 = open('%s/%s.zip'%(outputfolderzip,os.path.split(vcf_file)[-1]),'r')
    except IOError:
        os.system('zip -r -4 %s.zip %s'%(vcf_file,vcf_file))
        os.system('mv %s.zip %s/'%(vcf_file,outputfolderzip))
        os.system('rm %s' % (vcf_file))


################################################### END ########################################################
