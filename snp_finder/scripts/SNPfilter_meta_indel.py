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
try:
    os.mkdir(args.o + '/MG/covsum/')
except IOError:
    pass
try:
    os.mkdir(args.o + '/MG/indelsum/')
except IOError:
    pass

################################################### Set up ########################################################
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

checkdepth = False # 'MV', False, or True
mapping_qual = 42
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
    vcf_filename = os.path.split(vcf_file)[-1]
    lineage = vcf_filename.split(args.mfq + '.')[1].split('.raw.vcf')[0]
    outputfolder = args.o + '/MG/covsum/' + lineage
    try:
        os.mkdir(outputfolder)
    except IOError:
        pass
    f1 = open(os.path.join(outputfolder, '%s.depth.sum' % (vcf_filename)), 'w')
    f1.write('CHR\tdepth\tnum_loci\t\n')
    f2 = open(os.path.join(outputfolder, '%s.coverage.sum' % (vcf_filename)), 'w')
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
    vcf_filename = os.path.split(vcf_file)[-1]
    lineage = vcf_filename.split(args.mfq + '.')[1].split('.raw.vcf')[0]
    outputfolder = args.o + '/MG/covsum/' + lineage
    try:
        os.mkdir(outputfolder)
    except IOError:
        pass
    f1 = open(os.path.join(outputfolder,'%s.depth.MV.sum' % (vcf_filename)), 'w')
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
    def addSNP(self,CHR_POS):
        self.position.setdefault(CHR_POS,['','',0,0])
    def mapSNP(self,CHR_POS,Major_ALTmeta,Minor_ALTmeta,Major_ALTmeta_freq, Minor_ALTmeta_freq,gene,genePOS):
        self.position.setdefault(CHR_POS, [Major_ALTmeta,Minor_ALTmeta,Major_ALTmeta_freq, Minor_ALTmeta_freq,gene,genePOS])
    def sum_snpmeta(self,samplename):
        allsum = []
        for CHR_POS in self.position:
            Major_ALTmeta, Minor_ALTmeta, Major_ALTmeta_freq, Minor_ALTmeta_freq,gene,genePOS = self.position[CHR_POS]
            allsum.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n'%(CHR_POS,gene,genePOS,Major_ALTmeta,Minor_ALTmeta,Major_ALTmeta_freq, Minor_ALTmeta_freq))
        outputfolder = args.o + '/MG/indelsum/' + self.lineage
        try:
            os.mkdir(outputfolder)
        except IOError:
            pass
        f1 = open(os.path.join(outputfolder,'%s.IN.%s.snp.sum'%(samplename,self.lineage)),'w')
        f1.write('CHR\tPOS\tGene\tGenePOS\tMajor_ALT\tMinor_ALT\tMajor_ALT_freq\tMinor_ALT_freq\t\n')
        f1.write(''.join(allsum))
        f1.close()


def get_snp(lines_set,temp_SNP_lineage):
    CHR,POS = lines_set[0:2]
    CHR_POS = '%s\t%s' % (CHR, POS)
    gene = 'None'
    genePOS = POS
    POS = int(POS)
    if CHR in gene_loci:
        for genes in gene_loci[CHR]:
            genename,pos1,pos2 = genes
            if (POS >= pos1 and POS <= pos2):
                gene = genename
                genePOS = POS - pos1 + 1
    Major_ALTmeta, Minor_ALTmeta = lines_set[3:5]
    Depth4 = lines_set[7].split('DP4=')[1].split(';')[0].split(',')
    Major_ALTmeta_freq = int(Depth4[0]) + int(Depth4[1])
    Minor_ALTmeta_freq = int(Depth4[2]) + int(Depth4[3])
    temp_SNP_lineage.mapSNP(CHR_POS, Major_ALTmeta, Minor_ALTmeta, Major_ALTmeta_freq, Minor_ALTmeta_freq,gene,genePOS)

def sum_snp(temp_SNP_lineage,samplename):
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

def loadmultiple(sum_file):
    HS_gene = dict()
    for lines in open(sum_file, 'r'):
        if not lines.startswith("Refgenome"):
            lines_set = lines.split('\n')[0].split('\t')
            lineage, gene,POS = lines_set[0:3]
            HS_gene.setdefault(lineage, set())
            HS_gene[lineage].add('%s\t%s'%(gene,POS))
    return HS_gene

def loadgene(fasta):
    gene_loci = dict()
    for record in SeqIO.parse(fasta, 'fasta'):
        gene = str(record.id)
        CHR = '_'.join(gene.split('_')[:-1])
        pos1 = int(str(record.description).split('#')[1].replace(' ',''))
        pos2 = int(str(record.description).split('#')[2].replace(' ',''))
        gene_loci.setdefault(CHR,[]
                             )
        gene_loci[CHR].append([gene,pos1,pos2])
    return gene_loci


################################################### Main ########################################################
# run vcf depth and SNP summarizing
if args.vcf == 'None':
    all_vcf_file=glob.glob(os.path.join(output_dir,'*.raw.vcf'))
else:
    all_vcf_file = glob.glob(os.path.join(output_dir, '%s'%(args.vcf)))

for vcf_file in all_vcf_file:
    samplename = os.path.split(vcf_file)[-1].split(args.mfq)[0]
    donorspecies = os.path.split(vcf_file)[-1].split(args.mfq + '.')[1].split('.raw.vcf')[0]
    print('processing',samplename,donorspecies)
    #snp_folder = glob.glob('%s/%s.*.all.parsi.fasta.sum.txt'%(args.snp,donorspecies))
    #print(snp_folder)
    #allsnp = load_snp(snp_folder)
    assembly_file = glob.glob('%s/../co-assembly/withPE/%s/%s.all.spades1.fasta.noHM.fasta.faa'%(args.snp,donorspecies,donorspecies))+ \
                    glob.glob('%s/../co-assembly/withPE_new/%s/%s.all.spades1.fasta.noHM.fasta.faa' % (
                    args.snp, donorspecies, donorspecies))+ \
                    glob.glob('%s/../co-assembly/%s/%s.all.spades1.fasta.noHM.fasta.faa' % (
                        args.snp, donorspecies, donorspecies))
    print(assembly_file)
    gene_loci = loadgene(assembly_file[0])
    temp_SNP_lineage = SNP_lineage()
    temp_SNP_lineage.init(donorspecies)
    Donorspecies_sample = dict()
    Donorspecies_sample_MV = dict()
    Length = dict()
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#") and 'INDEL' in lines:
            lines_set = lines.split('\n')[0].split('\t')
            quality = float(lines_set[5])
            if quality >= mapping_qual:
                if checkdepth:
                    # check depth of moving window
                    vcf_to_depth_sum_moving_window(lines_set)
                    # check depth of each locus
                    vcf_to_depth(lines_set)
                if lines_set[4]!= '.' and len(lines_set[4]) > 1:
                    if lines_set[3] != '.':
                        codondis = len(lines_set[4])-len(lines_set[3])
                    else:
                        codondis = len(lines_set[4])
                    if not ',' in lines_set[4] and codondis%3!=0:
                        # check mapping quality and no >2 genotypes
                        # map snp to snps found in a lineage
                        get_snp(lines_set,temp_SNP_lineage)
    if checkdepth:
        # summarize depth
        depth_sum(Donorspecies_sample, Length,vcf_file)
        depth_sum_MV(Donorspecies_sample_MV, Length, vcf_file)
    # output snps found in a lineage
    sum_snp(temp_SNP_lineage,samplename)
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
