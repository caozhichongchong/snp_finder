# start
# after round 4 calculate NS ratio
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of folders of WGS of each species",
                      type=str, default='.',
                      metavar='input/')
required.add_argument("-fq",
                      help="file extension of WGS fastq #1 files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')
# optional output setup

optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='snp_output/',
                      metavar='snp_output/')
# optional search parameters
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 1)",
                      metavar="1 or more", action='store', default=1, type=int)
optional.add_argument('-rd',
                      help="Round of SNP calling and filtering",
                      metavar="1-4", action='store', default=1, type=int)
optional.add_argument('-job',
                      help="Optional: command to submit jobs",
                      metavar="nohup or customized",
                      action='store', default='jobmit', type=str)
# requirement for software calling
optional.add_argument('-bw', '--bowtie',
                          help="Optional: complete path to bowtie if not in PATH",
                          metavar="/usr/local/bin/bowtie",
                          action='store', default='bowtie', type=str)
optional.add_argument('-sp', '--spades',
                          help="Optional: complete path to spades if not in PATH",
                          metavar="/usr/local/bin/spades",
                          action='store', default='spades', type=str)
optional.add_argument('-pro', '--prodigal',
                      help="Optional: complete path to prodigal if not in PATH, None for no prodigal (default)",
                      metavar="/usr/local/bin/prodigal",
                      action='store', default='None', type=str)
optional.add_argument('-bcf', '--bcftools',
                      help="Optional: complete path to bcftools if not in PATH",
                      metavar="/usr/local/bin/bcftools",
                      action='store', default='bcftools', type=str)
optional.add_argument('-sam', '--samtools',
                      help="Optional: complete path to bwa if not in PATH",
                      metavar="/usr/local/bin/samtools",
                      action='store', default='samtools', type=str)
optional.add_argument('-mini', '--minimap2',
                      help="Optional: complete path to minimap2 if not in PATH",
                      metavar="/usr/local/bin/minimap2",
                      action='store', default='minimap2', type=str)

################################################## Definition ########################################################
args = parser.parse_args()

# set up path
Round = args.rd
Cluster = True
Tree = True
Paircompare = False
Cov_dis = 20

input_script_sub = args.s + '/vcf_round%s_tree'%(Round)
input_script_sub_merge = args.s + '/vcf_round%s'%(Round)
input_script = args.s
genome_root = args.i + '/round*'
output_dir_merge = args.o +'/vcf_round%s/merge_genome/'%(Round)
vcf_name = '.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf'
ref_filename = '.all.spades*.fasta'
fasta_name = '.fasta.corrected.fasta'
fastq_name = '.sorted.bam'

try:
    os.mkdir(output_dir_merge + '/summary')
except IOError:
    pass

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# Set up N or S
N_S_set = dict()
N_S_set['N']=0
N_S_set['S']=1
purines=['A','G']
pyrimidines=['C','T']
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
# Set up NS ratio cutoff
NSratioobserve_cutoff = 1.0
Min_SNP_highselect_cutoff = 1/2000
Max_SNP_highselect_cutoff = 0.02
countOther = False #do not count other SNPs (non-gene SNPs) when calculating expected NS ratio
corecutoff = 0.9
################################################### new class #########################################################
__metaclass__ = type

class SNP_gene:
    # create a class to store SNP_gene
    'a class to store SNP_gene'
    def init(self, gene):
        self.gene = gene
        self.position = dict()
        self.position.setdefault(gene,set())
        # [[N,S],freq]
        # not observed but predicted SNP pair by codon freq
        self.SNP_pair = {'A-T': [0,0],
                         'A-C': [0,0],
                         'G-C': [0,0],
                         'G-T': [0,0],
                         'A-G': [0,0],
                         'G-A': [0,0]}
        self.SNP_pair_freq = {'A-T': 0,
                         'A-C': 0,
                         'G-C': 0,
                         'G-T': 0,
                         'A-G': 0,
                         'G-A': 0}
        self.depth = []
        self.cov = 0
        self.mutposition = dict()
        self.mutposition.setdefault(gene, set())
        self.NSratio = [0,0,0]
        self.protein = ''
        self.minor_freq = []
    def addmutposition(self, gene, position):
        self.mutposition.setdefault(gene, set())
        self.mutposition[gene].add(position)
    def deletemutposition(self, gene,position):
        self.mutposition.setdefault(gene, set())
        if position == 0 :
            self.mutposition.pop(gene, None)
        else:
            self.mutposition[gene].discard(position)
    def addposition(self, gene, position,depth):
        self.position.setdefault(gene, set())
        self.position[gene].add(position)
        self.cov += 1
        self.depth.append(depth)
    def deleteposition(self, gene, position):
        if position == 0 :
            self.position.pop(gene, None)
        else:
            self.position[gene].discard(position)
    def addprotein(self, aa):
        self.protein += aa
    def mutpositioncal(self):
        self.mutpositionsum1 = {0: 0,
                         1: 0,
                         2: 0}
        self.mutpositionsum2 = {0: 0,
                                1: 0,
                                2: 0}
        self.mutpositionsum = 0
        for Chr in self.mutposition:
            self.mutposition2 = list(self.mutposition[Chr])
            self.mutposition2.sort()
            Total = len(self.mutposition2)
            self.mutpositionsum += Total
            if Total >= 1:
                self.mutpositionsum1[abs(self.mutposition2[0]) % 3] += 1
                self.mutpositionsum2[(self.mutposition2[-1]) % 3] += 1
                if Total != 1:
                    for i in range(0,len(self.mutposition2)-1):
                        self.mutpositionsum1[abs(self.mutposition2[i+1] - self.mutposition2[i]) % 3] += 1
                        self.mutpositionsum2[(self.mutposition2[i]) % 3] += 1
    def addSNP_pair(self, pair, position, count, unique_snp_count,depth = 0):
        self.SNP_pair_freq[pair] += unique_snp_count
        self.NSratio[position] += unique_snp_count
        if position < 2 and depth == 0:
            # add to NS for each SNP pair of reference genes
            self.SNP_pair[pair][position] += unique_snp_count
    def addpredictSNP_pair(self,refSNP_pair_sum):
        for pair in refSNP_pair_sum:
            self.SNP_pair[pair][0] += refSNP_pair_sum[pair][0]
            self.SNP_pair[pair][1] += refSNP_pair_sum[pair][1]
    def addalt(self,AllALT_frq):
        self.minor_freq.append(AllALT_frq)
    def deleteSNP_pair(self,SNP_gene):
        for position in [0,1]:
            self.NSratio[position] -= SNP_gene.NSratio[position]
            for pair in SNP_gene.SNP_pair:
                self.SNP_pair[pair][position] -= SNP_gene.SNP_pair[pair][position]
                if position == 0:
                    self.SNP_pair_freq[pair] -= SNP_gene.SNP_pair_freq[pair]
    def sum_SNP_pair(self):
        self.SNP_pair_sum = {'A-T': [0, 0],
                         'A-C': [0, 0],
                         'G-C': [0, 0],
                         'G-T': [0, 0],
                         'A-G': [0, 0],
                         'G-A': [0, 0]}
        for pair in self.SNP_pair:
            self.SNP_pair_sum[pair][0] += self.SNP_pair[pair][0]
            self.SNP_pair_sum[pair][1] += self.SNP_pair[pair][1]
    def dN_dS(self,SNP_gene_all,normalize=0):
        self.expectNSratio = 'No_expect'
        expectNSratio = [0, 0]
        if normalize == 1:
            for pair in self.SNP_pair_freq:
                # use selected codon NS ratio (SNP pair) * all genes freq
                expectNSratio[0] += self.SNP_pair[pair][0] * SNP_gene_all.SNP_pair_freq[pair]
                expectNSratio[1] += self.SNP_pair[pair][1] * SNP_gene_all.SNP_pair_freq[pair]
        else:
            for pair in self.SNP_pair_freq:
                expectNSratio[0] += self.SNP_pair[pair][0] * self.SNP_pair_freq[pair]
                expectNSratio[1] += self.SNP_pair[pair][1] * self.SNP_pair_freq[pair]
        if expectNSratio[1] > 0:
            # no expect S
            self.expectNSratio = expectNSratio[0] / expectNSratio[1]
        elif self.NSratio[0] == 0:
            # only S observed
            self.expectNSratio = 'expect_None'
        else:
            # only N observed
            self.expectNSratio = 'expect_N_only'
        if self.NSratio[1] > 0:
            # S observed
            self.NSratiosum = self.NSratio[0] / self.NSratio[1]
        elif self.NSratio[0] == 0:
            # only S observed
            self.NSratiosum = 'observe_None'
        else:
            # only N observed
            self.NSratiosum = 'observe_N_only'
        self.dNdS = self.NSratiosum
        if type(self.expectNSratio) == float \
                and type(self.NSratiosum) == float \
                and self.expectNSratio > 0:
            # calculate dNdS only if expected NSratio and observed NSratio are effective
            self.dNdS = self.NSratiosum / self.expectNSratio
            self.dNdS = '%.3f' % (self.dNdS)
            self.NSratiosum = '%.3f' % (self.NSratiosum)
            self.expectNSratio = '%.3f' % (self.expectNSratio)
        elif type(self.expectNSratio) == float:
            self.expectNSratio = '%.3f' % (self.expectNSratio)
        elif type(self.NSratiosum) == float:
            self.NSratiosum = '%.3f' % (self.NSratiosum)

################################################### Function ########################################################

def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def dnORds(amino1, amino2):
    if amino1 == amino2:
        return 'S'
    else:
        return 'N'

def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        if ALT_frq > 0:
            ALT_set.setdefault(ALT_frq, set())
            ALT_set[ALT_frq].add(alleles)
            ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def transitions(REF,ALT):
    if REF in pyrimidines:
        REF = complement[REF]
        ALT = complement[ALT]
    return '%s-%s'%(REF,ALT)

def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = int(Subdepth[0]) + int(Subdepth[1])
    Subdepth_ALT = int(Subdepth[2]) + int(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return [allels_set[1],1]
    else:
        return [allels_set[0],0]

def expectNSsub(record_name,record_seq,position=0):
    Total = int(len(record_seq)/3)
    temp_SNP_gene = SNP_gene()
    temp_SNP_gene.init(record_name)
    for i in range(0, (Total - position)):
        codon = record_seq[(i * 3 + position):((i + 1) * 3 + position)]
        try:
            codon_NSratio = codontable_NSratio[codon]
            temp_SNP_gene.addprotein(codontable[codon])
            for pair in codon_NSratio.SNP_pair:
                temp_SNP_gene.addSNP_pair(pair, 0, codon_NSratio.SNP_pair[pair][0],codon_NSratio.SNP_pair[pair][0],0)
                temp_SNP_gene.addSNP_pair(pair, 1, codon_NSratio.SNP_pair[pair][1],codon_NSratio.SNP_pair[pair][1],0)
        except KeyError:
            pass
    temp_SNP_gene.sum_SNP_pair()
    return [temp_SNP_gene.SNP_pair_sum,temp_SNP_gene.protein,position]

def expectNS(record_name,record_seq):
    Total = int(len(record_seq) / 3)
    temp_result = expectNSsub(record_name, record_seq)
    if len(temp_result[1]) < 0.8 * Total:
        temp_result = expectNSsub(record_name, record_seq,1)
        if len(temp_result[1]) < 0.8 * Total:
            temp_result = expectNSsub(record_name, record_seq,2)
            if len(temp_result[1]) < 0.8 * Total:
                return 'None'
    return temp_result

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Ref_NSratio=dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        Ref_seq.setdefault(record_id, record_seq)
        Ref_NSratio.setdefault(record_id,
                               expectNS(record_id, record_seq))
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    foutput = open(database + '.ref.NS.ratio', 'w')
    foutput_list = []
    #for Ref in Ref_NSratio:
    #    foutput_list.append('%s\t%s\t\n' % (Ref, Ref_NSratio[Ref]))
    #foutput.write(''.join(foutput_list))
    #foutput.close()
    return [Ref_seq,Ref_NSratio,Mapping,Mapping_loci,Reverse]

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

def sumgene(SNP_gene_temp,genome_set_snp,donor_species,SNP_gene_species,Total,normalize = 1):
    High_select = False
    N_temp = SNP_gene_temp.NSratio[0]
    S_temp = SNP_gene_temp.NSratio[1]
    Other_temp = SNP_gene_temp.NSratio[2]
    N_S_sum = N_temp + S_temp
    total_SNP = N_S_sum + Other_temp
    Gene_length = 1000
    Chr = SNP_gene_temp.gene
    total_SNP_position = len(SNP_gene_temp.mutposition[Chr])
    if Chr in Mapping:
        Gene_length = Mapping[Chr]
    if Chr == 'allspecies':
        Gene_length = 0
        for allchr in SNP_gene_all.position:
            Gene_length += Mapping.get(allchr, 0)
    new_line = '%s\t%s\t%s\t%s\t%s' % (donor_species, Chr, Total,
                                       Gene_length, genome_set_snp)
    if total_SNP > 0:
        new_line += '\t%s\t%s' % (total_SNP,total_SNP_position)
        SNP_gene_temp.dN_dS(SNP_gene_species, normalize)  # normalized
        new_line += ('\t%s\t%s\t%s\t%s' % (N_temp, S_temp, Other_temp,SNP_gene_temp.NSratiosum))
        new_line += '\t%s\t%s' % (SNP_gene_temp.expectNSratio, SNP_gene_temp.dNdS)
        for pair in SNP_gene_temp.SNP_pair:
            pair_freq = SNP_gene_temp.SNP_pair_freq[pair]
            pair_N = SNP_gene_temp.SNP_pair[pair][0]
            pair_S = SNP_gene_temp.SNP_pair[pair][1]
            new_line += ('\t%s\t%s:%s' % ('%d' % pair_freq, pair_N, pair_S))
    if N_S_sum > 0 and genome_set_snp > 1 \
            and total_SNP_position >= 2 and \
            total_SNP_position / Gene_length >= Min_SNP_highselect_cutoff \
            and total_SNP_position / Gene_length <= Max_SNP_highselect_cutoff and\
        SNP_gene_temp.NSratio[0] > SNP_gene_temp.NSratio[1] * NSratioobserve_cutoff:
        High_select = True
    new_line += '\t%s\n'%(High_select)
    return [new_line,High_select]

def freq_call(vcf_file,Ref_seq, Ref_NSratio,SNP_gene_species,SNP_gene_all,SNP_gene_all_highselect,SNP_gene_all_flexible,SNP_gene_all_core,SNP_gene_species_highselect,Output2,donor_species):
    Output = []
    all_SNP_gene_temp = dict()
    Total = 0
    SNP_type = dict()
    for lines in open(vcf_file, 'r'):
        lines_set = lines.replace('\n','').replace('\r','').split('\t')
        # set up the basic
        if Total == 0:
            Total = len(lines_set) - 9
        Chr = lines_set[0]
        position = int(lines_set[1])
        gene_info = contig_to_gene(Chr, position)
        if gene_info != []:
            Chr, position, codon_start, Ref_seq_chr, Reverse_chr = gene_info# a gene
        else:
            Chr = Chr + '_other' # not a gene
        Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
        REF = lines_set[3]
        ALT_set = lines_set[4].split(',')
        allels_set = [REF] + ALT_set
        Total_alleles = len(allels_set)
        SNP_count_genome_count = [[0] * Total_alleles, '', '']
        SNP_type.setdefault(Chr, [''] * Total)
        New_gene = 0
        if Chr not in all_SNP_gene_temp:
            SNP_gene_temp = SNP_gene()
            SNP_gene_temp.init(Chr)
            all_SNP_gene_temp.setdefault(Chr, SNP_gene_temp)
            New_gene = 1
        # set up SNP_gene
        SNP_gene_temp = all_SNP_gene_temp[Chr]
        SNP_gene_temp.addposition(Chr, position, Depth)
        SNP_gene_all.addposition(Chr, position, Depth)
        SNP_gene_species.addposition(Chr, position, Depth)
        # count genome set with SNP
        genome_ID = 0
        for Subdepth_all in lines_set[9:]:
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '')
            Subdepth_set = Subdepth.split(',')
            Subdepth_set_int = []
            for sub_depth in Subdepth_set:
                Subdepth_set_int.append(int(sub_depth))
            major_alt_frq = max(Subdepth_set_int)
            major_alt_frq_index = Subdepth_set_int.index(major_alt_frq)
            major_alt = allels_set[major_alt_frq_index]
            SNP_count_genome_count[0][major_alt_frq_index] += 1
            SNP_count_genome_count[1] += major_alt
            SNP_type[Chr][genome_ID] += major_alt
            genome_ID += 1
            SNP_count_genome_count[1] += '\t'
        # currate REF and ALT
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        REF, REF_where = curate_REF(allels_set, Depth4)
        ALT_set = allels_set
        ALT_set.remove(REF)
        # calculate N or S
        refSNP_pair_sum_all = Ref_NSratio.get(Chr, 'None')
        if refSNP_pair_sum_all != 'None':
            # predicted NS store in SNP_pair
            if New_gene == 1:
                refSNP_pair = refSNP_pair_sum_all[0]
                SNP_gene_temp.addpredictSNP_pair(refSNP_pair)
                SNP_gene_species.addpredictSNP_pair(refSNP_pair)
                SNP_gene_all.addpredictSNP_pair(refSNP_pair)
            #  observed NS ratio calculated
            refSNP_condon_start = refSNP_pair_sum_all[-1]
            codon_start = position - 1 - int((position - 1) % 3) + refSNP_condon_start
            if codon_start <= position - 1:
                Ref_seq_chr = Ref_seq[Chr]
                SNP_seq_chr = Ref_seq_chr
                Ref_seq_chr = causeSNP(Ref_seq_chr, position, REF,Reverse_chr)
                Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                if len(Ref_seq_codon) == 3:
                    Ref_seq_aa = translate(Ref_seq_codon)[0]
                    ALT_num = 1
                    for ALT in ALT_set:
                        ALT_frq = SNP_count_genome_count[0][ALT_num]
                        SNP_seq_chr = causeSNP(SNP_seq_chr, position, ALT,Reverse_chr)
                        SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                        SNP_seq_aa = translate(SNP_seq_codon)[0]
                        temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                        SNP_count_genome_count[2] += temp_NorS
                        SNP_pair = transitions(REF, ALT)
                        ALT_num += 1
                        SNP_gene_temp.addSNP_pair(SNP_pair, N_S_set[temp_NorS],
                                                  ALT_frq, 1, Depth)
                        SNP_gene_temp.addmutposition(Chr, position)
                        SNP_gene_all.addSNP_pair(SNP_pair, N_S_set[temp_NorS],
                                                 ALT_frq, 1, Depth)
                        SNP_gene_all.addmutposition('allspecies', position)
                        SNP_gene_species.addSNP_pair(SNP_pair, N_S_set[temp_NorS],
                                                     ALT_frq, 1, Depth)
                        SNP_gene_species.addmutposition(donor_species, position)
        else:
            # not a gene
            ALT_num = 1
            for ALT in ALT_set:
                ALT_frq = SNP_count_genome_count[0][ALT_num]
                SNP_pair = transitions(REF, ALT)
                # add to P
                SNP_gene_temp.addSNP_pair(SNP_pair, 2,
                                          ALT_frq, 1, Depth)
                SNP_gene_temp.addmutposition(Chr, position)
                if countOther:
                    SNP_gene_all.addSNP_pair(SNP_pair, 2,
                                             ALT_frq, 1, Depth)
                    SNP_gene_all.addmutposition('allspecies', position)
                SNP_gene_species.addSNP_pair(SNP_pair, 2,
                                             ALT_frq, 1, Depth)
                SNP_gene_species.addmutposition(donor_species, position)
                ALT_num += 1
        # output genome SNP
        Output.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n' % (Chr, position, REF, ALT,
                                                              SNP_count_genome_count[0][0],
                                                              SNP_count_genome_count[0][1],
                                                              SNP_count_genome_count[2], SNP_count_genome_count[1]))
    if all_SNP_gene_temp!= dict():
        for Chr in all_SNP_gene_temp:
            SNP_gene_temp = all_SNP_gene_temp[Chr]
            total_genome_set = len(set(SNP_type[Chr])) - 1
            sumgene_line,High_select = sumgene(SNP_gene_temp,total_genome_set,donor_species,SNP_gene_species,Total,1)
            Output2.append(sumgene_line)
            if High_select:
                SNP_gene_all_highselect.addpredictSNP_pair(SNP_gene_temp.SNP_pair)
                SNP_gene_species_highselect.addpredictSNP_pair(SNP_gene_temp.SNP_pair)
                for pair in SNP_gene_temp.SNP_pair_freq:
                    # use selected genes frequency * all genes NS ratio (codon NS sum of all genes)
                    SNP_gene_all_highselect.SNP_pair[pair][0] += SNP_gene_temp.SNP_pair[pair][0]
                    SNP_gene_all_highselect.SNP_pair[pair][1] += SNP_gene_temp.SNP_pair[pair][1]
                    SNP_gene_all_highselect.SNP_pair_freq[pair] += SNP_gene_temp.SNP_pair_freq[pair]
                    SNP_gene_species_highselect.SNP_pair[pair][0] += SNP_gene_temp.SNP_pair[pair][0]
                    SNP_gene_species_highselect.SNP_pair[pair][1] += SNP_gene_temp.SNP_pair[pair][1]
                    SNP_gene_species_highselect.SNP_pair_freq[pair] += SNP_gene_temp.SNP_pair_freq[pair]
                SNP_gene_all_highselect.NSratio[0] += SNP_gene_temp.NSratio[0]
                SNP_gene_all_highselect.NSratio[1] += SNP_gene_temp.NSratio[1]
                SNP_gene_species_highselect.NSratio[0] += SNP_gene_temp.NSratio[0]
                SNP_gene_species_highselect.NSratio[1] += SNP_gene_temp.NSratio[1]
            if '%s:%s'%(donor_species,Chr) in Core:
                core = Core['%s:%s'%(donor_species,Chr)]
                if core == 'all_core':
                    SNP_gene_all_core.addpredictSNP_pair(SNP_gene_temp.SNP_pair)
                    for pair in SNP_gene_temp.SNP_pair_freq:
                        # use selected genes frequency * all genes NS ratio (codon NS sum of all genes)
                        SNP_gene_all_core.SNP_pair[pair][0] += SNP_gene_temp.SNP_pair[pair][0]
                        SNP_gene_all_core.SNP_pair[pair][1] += SNP_gene_temp.SNP_pair[pair][1]
                        SNP_gene_all_core.SNP_pair_freq[pair] += SNP_gene_temp.SNP_pair_freq[pair]
                    SNP_gene_all_core.NSratio[0] += SNP_gene_temp.NSratio[0]
                    SNP_gene_all_core.NSratio[1] += SNP_gene_temp.NSratio[1]
                elif core == 'species_flexible':
                    SNP_gene_all_flexible.addpredictSNP_pair(SNP_gene_temp.SNP_pair)
                    for pair in SNP_gene_temp.SNP_pair_freq:
                        # use selected genes frequency * all genes NS ratio (codon NS sum of all genes)
                        SNP_gene_all_flexible.SNP_pair[pair][0] += SNP_gene_temp.SNP_pair[pair][0]
                        SNP_gene_all_flexible.SNP_pair[pair][1] += SNP_gene_temp.SNP_pair[pair][1]
                        SNP_gene_all_flexible.SNP_pair_freq[pair] += SNP_gene_temp.SNP_pair_freq[pair]
                    SNP_gene_all_flexible.NSratio[0] += SNP_gene_temp.NSratio[0]
                    SNP_gene_all_flexible.NSratio[1] += SNP_gene_temp.NSratio[1]
            else:
                print('missing genes %s in %s'%(Chr, donor_species))
        foutput = open(vcf_file + '.frq.snp', 'w')
        foutput.write('#CHR\tPOS\tMajor_ALT\tMinor_ALT\tMajor_ALT_frq\tMinor_Alt_frq\tN_or_S\tgenotype_allgenomes\n')
        foutput.write(''.join(Output))
        foutput.close()
    return Total

# calculate codon freq
codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
codontable_NSratio = dict()
for codon in codontable:
    SNP_gene_temp = SNP_gene()
    SNP_gene_temp.init(codon)
    codontable_NSratio.setdefault(codon, SNP_gene_temp)
    for position in range(0, 3):
        REF = codon[position]
        for ALT in Allels_order:
            if ALT != REF:
                new_codon = causeSNP(codon, position+1, ALT,0)
                temp_NorS = dnORds(translate(codon)[0], translate(new_codon)[0])
                SNP_pair = transitions(REF, ALT)
                SNP_gene_temp.addSNP_pair(SNP_pair, N_S_set[temp_NorS], 1, 1, 0)

# load core/flexible genes
Core = dict()
for lines in open(os.path.join(input_script,'all.denovo.gene.faa.allpangenome.sum.allsum.species.multispecies.txt')):
    lines_set = lines.replace('\n', '').replace('\r', '').split('\t')
    donor_species = lines_set[-3]
    geneID = lines_set[-1]
    core = lines_set[-4]
    Core.setdefault('%s:%s'%(donor_species,geneID),core)

# set up sum all species
SNP_gene_all = SNP_gene() # all denovo mutation
SNP_gene_all.init('allspecies')
SNP_gene_all_highselect = SNP_gene() # all mutations of highly selected genes
SNP_gene_all_highselect.init('allspecies_highselect')
SNP_gene_all_flexible = SNP_gene() # all mutations of flexible genes
SNP_gene_all_flexible.init('allspecies_flexible')
SNP_gene_all_core = SNP_gene() # all mutations of flexible genes
SNP_gene_all_core.init('allspecies_core')
# process each vcf file
Output2 = []
all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))

for vcf_file in all_vcf_file:
    print(vcf_file)
    vcf_ref_file_name = os.path.split(vcf_file)[-1].split('.all.')[0]
    donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0].split('.flt.snp.vcf')[0]
    database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
    if len(database) > 1:
        print(vcf_file, database)
    database = database[0]
    ref_dir, ref_name = os.path.split(database)
    database_file = database.replace('.fasta', '.fna')
    print('running %s' % donor_species)
    Ref_seq, Ref_NSratio, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
    SNP_gene_species = SNP_gene()  # all mutations of a species
    SNP_gene_species.init(donor_species)
    SNP_gene_species_highselect = SNP_gene()  # all mutations of highly selected genes in a species
    SNP_gene_species_highselect.init(donor_species + '_highselect')
    Total = freq_call(vcf_file, Ref_seq, Ref_NSratio, SNP_gene_species, SNP_gene_all, SNP_gene_all_highselect,
                      SNP_gene_all_flexible,SNP_gene_all_core,SNP_gene_species_highselect,
                          Output2, donor_species)
    if sum(SNP_gene_species.NSratio) > 0:
        # there's a SNP
        sumgene_line, High_select = sumgene(SNP_gene_species, 1, donor_species, SNP_gene_species, Total, 0)
        Output2.append(sumgene_line)
        print(
            SNP_gene_species.expectNSratio, SNP_gene_species.dNdS, SNP_gene_species.NSratio,
            SNP_gene_species.NSratiosum)
        if sum(SNP_gene_species_highselect.NSratio) > 0:
            sumgene_line, High_select = sumgene(SNP_gene_species_highselect, 1, donor_species, SNP_gene_species, Total, 1)
            Output2.append(sumgene_line)
            print(
                SNP_gene_species_highselect.expectNSratio, SNP_gene_species_highselect.dNdS, SNP_gene_species_highselect.NSratio,
                SNP_gene_species_highselect.NSratiosum)
    else:
        Output2.append('%s\t0\t\n'%(donor_species))

# sum all species dNdS
sumgene_line,High_select = sumgene(SNP_gene_all,1,'allspecies',SNP_gene_all,'None',0)
Output2.append(sumgene_line)
# HS genes
sumgene_line,High_select = sumgene(SNP_gene_all_highselect,1,'allspecies',SNP_gene_all,'None',1)
Output2.append(sumgene_line)
# flexible genes
sumgene_line,High_select = sumgene(SNP_gene_all_flexible,1,'allspecies',SNP_gene_all,'None',1)
Output2.append(sumgene_line)
# core genes
sumgene_line,High_select = sumgene(SNP_gene_all_core,1,'allspecies',SNP_gene_all,'None',1)
Output2.append(sumgene_line)

# output
foutput = open(output_dir_merge + '/summary/all.donor.species.dnds.txt', 'w')
foutput.write('#donor_species\tgene\tNo.genome\tgene_length\ttotal_SNP_genomeset\tNo.SNP\tNo.SNP_position\tN\tS\tOther\tobserved_ratio\texpected_ratio\tdNdS\t' +\
            'A-T_freq\tA-T_N:S\tA-C_freq\tA-C_N:S\tG-C_freq\tG-C_N:S\tG-T_freq\tG-T_N:S\tA-G_freq\tA-G_N:S\tG-A_freq\tG-A_N:S\tHigh_selected\n')
foutput.write(''.join(Output2))
foutput.close()

print(SNP_gene_all.expectNSratio,SNP_gene_all.dNdS,
      SNP_gene_all.NSratio,SNP_gene_all.NSratiosum)
print(SNP_gene_all_highselect.expectNSratio,SNP_gene_all_highselect.dNdS,
      SNP_gene_all_highselect.NSratio,SNP_gene_all_highselect.NSratiosum)
print(SNP_gene_all_flexible.expectNSratio,SNP_gene_all_flexible.dNdS,
      SNP_gene_all_flexible.NSratio,SNP_gene_all_flexible.NSratiosum)
print(SNP_gene_all_core.expectNSratio,SNP_gene_all_core.dNdS,
      SNP_gene_all_core.NSratio,SNP_gene_all_core.NSratiosum)

################################################### END ########################################################
