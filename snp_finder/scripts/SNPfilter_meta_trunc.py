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
optional.add_argument("-snplist",
                      help="a list of snps to look for (output SNPs in all MG)",
                      type=str, default='None',
                      metavar='snp_list.txt')
optional.add_argument("-trunclist",
                      help="a list of trunc to look for (output SNPs before truncation)",
                      type=str, default='None',
                      metavar='trunc_list.txt')
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

if args.snplist != 'None':
    try:
        os.mkdir(args.o + '/MG/truncsumHSnew/')
    except IOError:
        pass
elif args.trunclist != 'None':
    try:
        os.mkdir(args.o + '/MG/truncallsnps/')
    except IOError:
        pass
else:
    try:
        os.mkdir(args.o + '/MG/truncsumHS/')
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

mapping_qual = 42
min_separate_allele_freq = 5
################################################### Function ########################################################
# set up functions
__metaclass__ = type

class SNP_lineage:
    # create a class to store SNP_lineage
    'a class to store SNP_lineage'
    def init(self, lineage):
        self.lineage = lineage
        self.position = dict()
    def addSNP(self,CHR_POS):
        self.position.setdefault(CHR_POS,['','',0,0])
    def mapSNP(self,CHR_POS,Major_ALTmeta,Minor_ALTmeta,Major_ALTmeta_freq, Minor_ALTmeta_freq,gene,genePOS,Ref_seq_aa,SNP_seq_aa,withinHS):
        self.position.setdefault(CHR_POS, [Major_ALTmeta,Minor_ALTmeta,Major_ALTmeta_freq, Minor_ALTmeta_freq,gene,genePOS,Ref_seq_aa,SNP_seq_aa,withinHS])
    def sum_snpmeta(self,samplename):
        allsum = []
        for CHR_POS in self.position:
            Major_ALTmeta, Minor_ALTmeta, Major_ALTmeta_freq, Minor_ALTmeta_freq,gene,genePOS,Ref_seq_aa,SNP_seq_aa,withinHS = self.position[CHR_POS]
            allsum.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n'%(CHR_POS,gene,genePOS,withinHS,Major_ALTmeta,Minor_ALTmeta,Major_ALTmeta_freq, Minor_ALTmeta_freq,Ref_seq_aa,SNP_seq_aa))
        if args.snplist != 'None':
            outputfolder = args.o + '/MG/truncsumHSnew/' + self.lineage
        elif args.trunclist != 'None':
            outputfolder = args.o + '/MG/truncallsnps/' + self.lineage
        else:
            outputfolder = args.o + '/MG/truncsumHS/' + self.lineage
        try:
            os.mkdir(outputfolder)
        except IOError:
            pass
        f1 = open(os.path.join(outputfolder,'%s.IN.%s.snp.sum'%(samplename,self.lineage)),'w')
        f1.write('CHR\tPOS\tGene\tGenePOS\twithinHS\tMajor_ALT\tMinor_ALT\tMajor_ALT_freq\tMinor_ALT_freq\tRef_aa\tSNP_aa\n')
        f1.write(''.join(allsum))
        f1.close()

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

def get_snp(lines_set,temp_SNP_lineage):
    CHR,POS = lines_set[0:2]
    CHR_POS = '%s\t%s' % (CHR, POS)
    POS = int(POS)
    Major_ALTmeta, Minor_ALTmeta = lines_set[3:5]
    # check truncation
    gene_info = contig_to_gene(CHR, POS)
    if gene_info != []:
        Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
        if Ref_seq_chr != 'None':
            if codon_start <= POS_gene - 1:
                Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, Major_ALTmeta, Reverse_chr)
                Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                SNP_seq_chr = Ref_seq_chr
                if len(Ref_seq_codon) == 3:
                    Ref_seq_aa = translate(Ref_seq_codon)[0]
                    SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, Minor_ALTmeta, Reverse_chr)
                    SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                    SNP_seq_aa = translate(SNP_seq_codon)[0]
                    if SNP_seq_aa == '*' and Ref_seq_aa!= '*' :
                        # early stop codon
                        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0].split(',')
                        Major_ALTmeta_freq = int(Depth4[0]) + int(Depth4[1])
                        Minor_ALTmeta_freq = int(Depth4[2]) + int(Depth4[3])
                        withinHS = 'False'
                        if Chr_gene in HS_gene_within_lineage:
                            withinHS = 'True'
                        temp_SNP_lineage.mapSNP(CHR_POS, Major_ALTmeta, Minor_ALTmeta, Major_ALTmeta_freq,
                                                Minor_ALTmeta_freq, Chr_gene, POS_gene,Ref_seq_aa,SNP_seq_aa,withinHS)

def get_snp_certainsnp(lines_set,temp_SNP_lineage):
    CHR,POS = lines_set[0:2]
    CHR_POS = '%s\t%s' % (CHR, POS)
    CHRPOS = '%s %s' % (CHR, POS)
    Major_ALTmeta, Minor_ALTmeta = lines_set[3:5]
    Depth4 = lines_set[7].split('DP4=')[1].split(';')[0].split(',')
    Major_ALTmeta_freq = int(Depth4[0]) + int(Depth4[1])
    Minor_ALTmeta_freq = int(Depth4[2]) + int(Depth4[3])
    Minor, Gene, Genepos = SNP_list2[CHRPOS]
    if Minor == Minor_ALTmeta:
        SNP_seq_aa = '*'
    else:
        SNP_seq_aa = ''
    temp_SNP_lineage.mapSNP(CHR_POS, Major_ALTmeta, Minor_ALTmeta, Major_ALTmeta_freq,
                            Minor_ALTmeta_freq, Gene, Genepos, '', SNP_seq_aa, '')

def get_snp_certainsnpnotrunc(lines_set,temp_SNP_lineage,truncationsnp):
    CHR,POS = lines_set[0:2]
    CHR_POS = '%s\t%s' % (CHR, POS)
    Major_ALTmeta, Minor_ALTmeta = lines_set[3:5]
    Depth4 = lines_set[7].split('DP4=')[1].split(';')[0].split(',')
    Major_ALTmeta_freq = int(Depth4[0]) + int(Depth4[1])
    Minor_ALTmeta_freq = int(Depth4[2]) + int(Depth4[3])
    if Major_ALTmeta_freq > min_separate_allele_freq or Minor_ALTmeta_freq > min_separate_allele_freq:
        if truncationsnp:
            SNP_seq_aa = '*'
        else:
            SNP_seq_aa = ''
        temp_SNP_lineage.mapSNP(CHR_POS, Major_ALTmeta, Minor_ALTmeta, Major_ALTmeta_freq,
                                Minor_ALTmeta_freq, Gene, Genepos, '', SNP_seq_aa, '')

def sum_snp(temp_SNP_lineage,samplename):
    temp_SNP_lineage.sum_snpmeta(samplename)

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

def load_snp(snpfile):
    SNP_list = dict()
    SNP_list2 = dict()
    if snpfile != 'None':
        for lines in open(snpfile,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            lineage,CHRPOS,Major,Minor,Gene,Genepos = lines_set[:6]
            SNP_list.setdefault(lineage,set())
            SNP_list[lineage].add(CHRPOS)
            SNP_list2.setdefault(CHRPOS, [Minor,Gene,Genepos])
    return[SNP_list,SNP_list2]

################################################### Main ########################################################
# load HS genes
HS_gene_within = loadHS('%s/summary/all.species.txt'%(args.snp))

# load SNP list
SNP_list,SNP_list2 = load_snp(args.snplist)

Trunc_list = dict()
if args.trunclist!='None':
    for alltrunc in args.trunclist.split(';'):
        chr,pos,gene,genepos = alltrunc.split(' ')
        Trunc_list.setdefault(chr,[])
        Trunc_list[chr].append([int(pos),gene])

# run SNP truncation prediction
if args.vcf == 'None':
    all_vcf_file=glob.glob(os.path.join(output_dir,'*.raw.vcf'))
else:
    all_vcf_file = glob.glob(os.path.join(output_dir, '%s'%(args.vcf)))

for vcf_file in all_vcf_file:
    samplename = os.path.split(vcf_file)[-1].split(args.mfq)[0]
    donorspecies = os.path.split(vcf_file)[-1].split(args.mfq + '.')[1].split('.raw.vcf')[0]
    print('processing', samplename, donorspecies)
    if SNP_list!= dict():
        temp_SNP_lineage = SNP_lineage()
        temp_SNP_lineage.init(donorspecies)
        for lines in open(vcf_file, 'r'):
            if not lines.startswith("#") and 'INDEL' not in lines:
                lines_set = lines.split('\n')[0].split('\t')
                quality = float(lines_set[5])
                CHRPOS = '%s %s'%(lines_set[0],lines_set[1])
                if CHRPOS in SNP_list.get(donorspecies) and quality >= mapping_qual:
                    # map snp to snps found in a lineage
                    get_snp_certainsnp(lines_set, temp_SNP_lineage)
        # output snps found in a lineage
        sum_snp(temp_SNP_lineage, samplename)
    else:
        assembly_file = glob.glob('%s/../co-assembly/withPE/%s/%s.all.spades1.fasta.noHM.fasta.fna'%(args.snp,donorspecies,donorspecies))+ \
                        glob.glob('%s/../co-assembly/withPE_new/%s/%s.all.spades1.fasta.noHM.fasta.fna' % (
                        args.snp, donorspecies, donorspecies))+ \
                        glob.glob('%s/../co-assembly/%s/%s.all.spades1.fasta.noHM.fasta.fna' % (
                            args.snp, donorspecies, donorspecies))
        Ref_seq, Mapping_loci, Reverse = loaddatabase(assembly_file[0])
        HS_gene_within_lineage = HS_gene_within.get(donorspecies,set())
        if Trunc_list != dict():
            temp_SNP_lineage = SNP_lineage()
            temp_SNP_lineage.init(donorspecies)
            for lines in open(vcf_file, 'r'):
                if not lines.startswith("#") and 'INDEL' not in lines:
                    lines_set = lines.split('\n')[0].split('\t')
                    quality = float(lines_set[5])
                    CHR, POS = lines_set[0:2]
                    POS = int(POS)
                    gene_info = contig_to_gene(CHR, POS)
                    if gene_info != []:
                        Gene, Genepos, codon_start, Ref_seq_chr, Reverse_chr = gene_info
                        if CHR in Trunc_list:
                            for alltrunc in Trunc_list[CHR]:
                                pos, gene = alltrunc
                                if pos >= POS and gene == Gene:
                                    # same gene, pos before trunc
                                    if quality >= mapping_qual and lines_set[4]!= '.' and not ',' in lines_set[4]:
                                        # check mapping quality and no >2 genotypes
                                        get_snp_certainsnpnotrunc(lines_set, temp_SNP_lineage,pos==POS)
                                    break
            # output snps found in a lineage
            sum_snp(temp_SNP_lineage, samplename)
        else:
            temp_SNP_lineage = SNP_lineage()
            temp_SNP_lineage.init(donorspecies)
            for lines in open(vcf_file, 'r'):
                if not lines.startswith("#") and 'INDEL' not in lines:
                    lines_set = lines.split('\n')[0].split('\t')
                    quality = float(lines_set[5])
                    if quality >= mapping_qual and lines_set[4]!= '.' and not ',' in lines_set[4]:
                        # check mapping quality and no >2 genotypes
                        # map snp to snps found in a lineage
                        get_snp(lines_set,temp_SNP_lineage)
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
