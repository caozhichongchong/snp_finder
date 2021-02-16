# start
# filter vcf of metagenomes for clonal populations
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
################################################### Set up ########################################################
# set up cutoff
Gene_cov_cutoff = 0.75 # coverage of a gene to count as mapped
Gene_num_cutoff = 0.8 # minimum percentage of genes in a donor species
SNP_max_cutoff = 0.1 # maximum ratio of SNPs on a gene to count as mapped
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

checkdepth = True # 'MV', False, or True
################################################### Function ########################################################
# set up functions
def vcf_to_depth(lines_set):
    CHR = lines_set[0]
    POS = int(lines_set[1])
    Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
    Donorspecies_sample.setdefault(CHR,[])
    Donorspecies_sample[CHR].append(Depth)
    if CHR not in Length:
        try:
            total_length = CHR.split('size')[1]
        except IndexError:
            try:
                total_length = CHR.split('length_')[1].split('_cov')[0]
            except IndexError:
                total_length = POS
        total_length = int(total_length)
        Length.setdefault(CHR, total_length)
    if 'size' not in CHR and 'length_' not in CHR:
        Length[CHR] = POS



def vcf_to_depth_sum_moving_window(lines_set):
    CHR = lines_set[0]
    POS = int(lines_set[1])
    if CHR not in Length:
        try:
            total_length = CHR.split('size')[1]
        except IndexError:
            try:
                total_length = CHR.split('length_')[1].split('_cov')[0]
            except IndexError:
                total_length = POS
        total_length = int(total_length)
        Length.setdefault(CHR, total_length)
    if 'size' not in CHR and 'length_' not in CHR:
        Length[CHR] = POS
    total_length = Length[CHR]
    Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
    Donorspecies_sample_MV.setdefault(CHR, [0]*int(total_length/1000 + 1))
    try:
        Donorspecies_sample_MV[CHR][int(POS/1000)] += Depth
    except IndexError:
        Donorspecies_sample_MV[CHR].append(Depth)

def depth_sum_chr(alldepth,CHR):
    CHR_output = set()
    alldepth_sum = dict()
    for depth in alldepth:
        alldepth_sum.setdefault(depth, 0)
        alldepth_sum[depth] += 1
    for depth in alldepth_sum:
        CHR_output.add('%s\t%s\t%s\t\n'%(CHR,depth,alldepth_sum[depth]))
    return CHR_output

def depth_sum(Donorspecies_sample, Length,vcf_file):
    f1 = open('%s.depth.sum' % (vcf_file), 'w')
    f1.write('CHR\tdepth\tnum_loci\t\n')
    f2 = open('%s.coverage.sum' % (vcf_file), 'w')
    f2.write('CHR\tdepth_covered\ttotal_length\t\n')
    CHR_coverage = set()
    for CHR in Donorspecies_sample:
        total_length = Length[CHR]
        alldepth = Donorspecies_sample[CHR]
        CHR_output = depth_sum_chr(alldepth,CHR)
        f1.write(''.join(list(CHR_output)))
        CHR_coverage.add('%s\t%s\t%s\t\n'%(CHR,len(alldepth),total_length))
    f1.close()
    f2.write(''.join(list(CHR_coverage)))
    f2.close()
    f1.close()

def depth_sum_MV(Donorspecies_sample_MV, Length,vcf_file):
    f1 = open('%s.depth.MV.sum' % (vcf_file), 'w')
    f1.write('CHR\tdepth\tnum_loci\t\n')
    CHR_coverage = set()
    for CHR in Donorspecies_sample_MV:
        total_length = Length[CHR]
        alldepth = Donorspecies_sample_MV[CHR]
        CHR_output = depth_sum_chr(alldepth,CHR)
        f1.write(''.join(list(CHR_output)))
        CHR_coverage.add('%s\t%s\t%s\t\n'%(CHR,len(alldepth),total_length))
    f1.close()
    f1.close()

__metaclass__ = type

class SNP_lineage:
    # create a class to store SNP_lineage
    'a class to store SNP_lineage'
    def init(self, lineage):
        self.lineage = lineage
        self.position = dict()
    def addSNP(self,CHR_POS,Major_ALT,Minor_ALT,withinHS,allHS,genename):
        self.position.setdefault(CHR_POS,[Major_ALT,Minor_ALT,0,0,0,withinHS,allHS,genename])
    def mapSNP(self,CHR_POS,Major_ALTmeta,Minor_ALTmeta,Major_ALTmeta_freq, Minor_ALTmeta_freq):
        try:
            self.position[CHR_POS][self.position[CHR_POS].index(Major_ALTmeta) + 2] += Major_ALTmeta_freq
        except ValueError:
            self.position[CHR_POS][4] += Major_ALTmeta_freq
        if Minor_ALTmeta_freq != '.' or Minor_ALTmeta_freq != 0:
            try:
                self.position[CHR_POS][self.position[CHR_POS].index(Minor_ALTmeta) + 2] += Minor_ALTmeta_freq
            except ValueError:
                self.position[CHR_POS][4] += Minor_ALTmeta_freq
    def sum_snpmeta(self,samplename):
        allsum = []
        for CHR_POS in self.position:
            Major_ALT, Minor_ALT, Major_ALT_freq, Minor_ALT_freq, other_freq,withinHS,allHS,genename = self.position[CHR_POS]
            allsum.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n'%(CHR_POS,Major_ALT, Minor_ALT, Major_ALT_freq, Minor_ALT_freq, other_freq,withinHS,allHS,genename))
        f1 = open(os.path.join(args.o + '/MG/bwa/','%s.IN.%s.snp.sum'%(samplename,self.lineage)),'w')
        f1.write('CHR\tPOS\tMajor_ALT\tMinor_ALT\tMajor_ALT_freq\tMinor_ALT_freq\tother_freq\tHS_within\tHS_all\tgenename\t\n')
        f1.write(''.join(allsum))
        f1.close()


def load_snp(snp_folder):
    allsnp = []
    for snp_file in snp_folder:
        lineage = os.path.split(snp_file)[-1].split('.all.parsi.fasta.sum.txt')[0]
        temp_SNP_lineage = SNP_lineage()
        temp_SNP_lineage.init(lineage)
        for lines in open(snp_file,'r'):
            if not lines.startswith("CHR\tPOS"):
                lines_set = lines.split('\t')
                CHR,POS,Major_ALT,Minor_ALT = lines_set[0:4]
                CHR_POS = '%s\t%s'%(CHR,POS)
                genename = lines_set[5]
                withinHS = 'False'
                allHS = 'False'
                if genename in HS_gene_within.get(lineage,set()):
                    withinHS = 'True'
                if genename in HS_gene_all.get(lineage,set()):
                    allHS = 'True'
                temp_SNP_lineage.addSNP(CHR_POS,Major_ALT,Minor_ALT,withinHS,allHS,genename)
        allsnp.append(temp_SNP_lineage)
    return allsnp

def get_snp(allsnp,lines_set):
    CHR,POS = lines_set[0:2]
    CHR_POS = '%s\t%s' % (CHR, POS)
    for temp_SNP_lineage in allsnp:
        if CHR_POS in temp_SNP_lineage.position:
            Major_ALTmeta, Minor_ALTmeta = lines_set[3:5]
            Depth4 = lines_set[7].split('DP4=')[1].split(';')[0].split(',')
            Major_ALTmeta_freq = int(Depth4[0]) + int(Depth4[1])
            Minor_ALTmeta_freq = int(Depth4[2]) + int(Depth4[3])
            temp_SNP_lineage.mapSNP(CHR_POS, Major_ALTmeta, Minor_ALTmeta, Major_ALTmeta_freq, Minor_ALTmeta_freq)

def sum_snp(allsnp,samplename):
    for temp_SNP_lineage in allsnp:
        temp_SNP_lineage.sum_snpmeta(samplename)

def loadHS(sum_file):
    HS_gene = dict()
    for lines in open(sum_file, 'r'):
        if not lines.startswith("#donor_species"):
            lines_set = lines.split('\n')[0].split('\t')
            if lines_set[-1] == 'True':
                lineage, gene = lines_set[0:2]
                HS_gene.setdefault(lineage,set())
                HS_gene[lineage].add(gene)
    return HS_gene

################################################### Main ########################################################
# load HS genes
HS_gene_within = loadHS('%s/summary/all.species.txt'%(args.snp))
HS_gene_all = loadHS('%s/summary/all.species.txt.High_select2.txt'%(args.snp))
# run vcf depth and SNP summarizing
if args.vcf == 'None':
    all_vcf_file=glob.glob(os.path.join(output_dir,'*.raw.vcf'))
else:
    all_vcf_file = glob.glob(os.path.join(output_dir, '%s'%(args.vcf)))

for vcf_file in all_vcf_file:
    samplename = os.path.split(vcf_file)[-1].split(args.mfq)[0]
    donorspecies = os.path.split(vcf_file)[-1].split(args.mfq + '.')[1].split('.raw.vcf')[0]
    print('processing',samplename,donorspecies)
    snp_folder = glob.glob('%s/%s.*.all.parsi.fasta.sum.txt'%(args.snp,donorspecies))
    print(snp_folder)
    allsnp = load_snp(snp_folder)
    Donorspecies_sample = dict()
    Donorspecies_sample_MV = dict()
    Length = dict()
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            lines_set = lines.split('\n')[0].split('\t')
            if checkdepth:
                # check depth of moving window
                vcf_to_depth_sum_moving_window(lines_set)
                # check depth of each locus
                vcf_to_depth(lines_set)
            # map snp to snps found in a lineage
            get_snp(allsnp,lines_set)
    if checkdepth:
        # summarize depth
        depth_sum(Donorspecies_sample, Length,vcf_file)
        depth_sum_MV(Donorspecies_sample_MV, Length, vcf_file)
    # output snps found in a lineage
    sum_snp(allsnp,samplename)
    os.system('zip -r -4 %s.zip %s'%(vcf_file,vcf_file))
    os.system('mv %s.zip %s/'%(vcf_file,output_dir + '/finished'))
    os.system('rm %s' % (vcf_file))


################################################### END ########################################################
