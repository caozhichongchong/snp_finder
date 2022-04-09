################################################### END ########################################################
################################################### SET PATH ########################################################
# After round 4 filter results of WGS
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import itertools
import random
from datetime import datetime
# set up path
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path to all vcf files",
                      type=str, default='.',
                      metavar='input/')
required.add_argument("-vcf",
                      help="file extension of vcfs with only SNPs(.flt.snp.vcf)",
                      type=str, default='.flt.snp.vcf',
                      metavar='.flt.snp.vcf')
# optional input genome
optional.add_argument("-cluster",
                      help="a cluster to run, default is all clusters",
                      type=str, default='',
                      metavar='cluster1')
# optional output setup
optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
optional.add_argument("-smp",
                      help="a folder to store all scripts for mapping",
                      type=str, default='scripts/mapping',
                      metavar='scripts/mapping')
# requirement for software calling
optional.add_argument('-pro',
                          help="Optional: complete path to prodigal if not in PATH",
                          metavar="/usr/local/bin/prodigal",
                          action='store', default='prodigal', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
# set up path
input_script = args.s
vcf_name = '.raw.vcf'
ref_filename = '.fasta'
bowtie_filter = True
################################################### Set up ########################################################
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels['-']=4
Allels_order = ['A','T','G','C','-']

min_qual_for_call = 20 #Remove sample*candidate that has lower than this quality
min_maf_for_call = .9 #Remove sample*candidate
min_cov = 3 # at least 3 reads mapped to POS
min_cov_for_call_per_strand_round = 3 #Remove sample*candidate
pair_strand = False # min_cov_for_call_per_strand for each straind

################################################### Function ########################################################
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 5):
        ALT_frq = int(Allels_count[alleles])
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


def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = float(Subdepth[0]) + float(Subdepth[1])
    Subdepth_ALT = float(Subdepth[2]) + float(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return [allels_set[1],1]
    else:
        return [allels_set[0],0]

def outputvcf(output_name,vcf_file_list_freq,vcf_file_list_vcf,Sample_name):
    vcf_file_filtered = open(vcf_file + '.%s.snpfreq.txt' % (output_name), 'w')
    vcf_file_filtered.write('CHR\tPOS\tMajor_ALT\tMinor_ALT\tGenome_set_noSNP\tGenome_set\tQuality\tGene\tGene_POS\tN_or_S\tAA_change\t%s\n'%('\t'.join(Sample_name))\
                            +''.join(vcf_file_list_freq))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.vcf' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list_vcf))
    vcf_file_filtered.close()


def SNP_seq(seq1, seq2, POS_info,POS_info_CHR,POS_info_CHR_LEN,POS_info_output,G1,G2):
    SNP_total = 0
    j = 0
    POS_DIS = []
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            CHR = POS_info_CHR[i]
            POS = POS_info[i]
            LEN = POS_info_CHR_LEN[CHR]
            if CHR == POS_info_CHR[j]:  # same CHR
                DIS = abs(POS - POS_info[j])
                POS_DIS.append(DIS)  # POS diff
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, DIS, LEN))
            else:  # new CHR
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, 0, LEN))
            j = i
    return SNP_total

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
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

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

def SNP_check_fq(lines_set,vcf_file_list,vcf_file_list_freq,vcf_file_list_vcf,vcf_file_POS_candidate, min_cov_for_call_per_strand):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_frq2 = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    REF = lines_set[3]
    allels_set = [REF]
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    SNP_quality = lines_set[5]
    allels_set += lines_set[4].split(',')
    Total_alleles = len(allels_set)
    Subdepth_all = lines_set_sub[0]
    Qual = float(SNP_quality)
    if Qual >= min_qual_for_call:
        # Keep SNP that has higher than this quality
        Allels_frq = [0, 0, 0, 0, 0]
        Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        Allels_frq_sub[0] = Qual
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
        total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
        Subdepth_forward = [int(i) for i in Subdepth_all.split(':')[-3].split(',')]
        Subdepth_reverse = [int(i) for i in Subdepth_all.split(':')[-2].split(',')]
        for num_allels in range(0, min(len(Subdepth), Total_alleles)):
            allels = allels_set[num_allels]
            #Subdepth_alleles = float(Subdepth[num_allels])
            if allels in Allels:
                forward = Subdepth_forward[num_allels]
                reverse = Subdepth_reverse[num_allels]
                Allels_frq[Allels[allels]] += reverse + forward
                Allels_frq_sub[Allels[allels] * 2 + 1] += forward
                Allels_frq_sub[Allels[allels] * 2 + 2] += reverse
            else:
                pass
        # find major alt and calculate frequency
        if sum(Allels_frq) > 0:
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            MAF = Major_ALT[1] / total_sub_depth
            if Major_ALT[0] not in [REF,'-']:
                if (MAF == 1 and total_sub_depth >= min_cov) or (MAF >= min_maf_for_call and
                                                                 total_sub_depth >= min_cov_for_call_per_strand * 2):
                    if not pair_strand or \
                            sum(Subdepth_forward) >= min_cov_for_call_per_strand and \
                            sum(Subdepth_reverse) >= min_cov_for_call_per_strand:
                            temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
                            temp_snp_line_frq2.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub[1:]))
                            # a potential SNP
                            # calculate NS
                            gene_info = contig_to_gene(CHR, POS)
                            if gene_info != []:
                                Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
                                if Ref_seq_chr != 'None':
                                    #  observed NS ratio calculated
                                    temp_snp_line_NS = [Chr_gene, str(POS_gene), '']
                                    if codon_start <= POS_gene - 1:
                                        Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                                        Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                                        SNP_seq_chr = Ref_seq_chr
                                        if len(Ref_seq_codon) == 3:
                                            Ref_seq_aa = translate(Ref_seq_codon)[0]
                                            temp_snp_line_AA += Ref_seq_aa
                                            SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, Major_ALT[0], Reverse_chr)
                                            SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                            SNP_seq_aa = translate(SNP_seq_codon)[0]
                                            temp_snp_line_AA += SNP_seq_aa
                                            temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                            temp_snp_line_NS[-1] += temp_NorS
                            # output lines and output major alt
                            temp_snp_line.append(CHR)
                            temp_snp_line.append(str(POS))
                            temp_snp_line.append(REF)
                            temp_snp_line.append(Major_ALT[0])
                            vcf_file_list.append(
                                '\t'.join(temp_snp_line) + '\t' + '\t'.join(temp_snp_line_frq2) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                                    '1', 'PASS', '\t'.join(temp_snp_line_NS),
                                    temp_snp_line_AA))
                            vcf_file_list_freq.append(
                                '\t'.join(temp_snp_line) + '\t\"%s\"\t\"%s\"\t%s\t%s\t%s\t%s\n' % (
                                    'None',
                                    '1',
                                    Qual, '\t'.join(temp_snp_line_NS),
                                    temp_snp_line_AA, '\t'.join(temp_snp_line_frq)))
                            vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
                            vcf_file_list_vcf.append('\t'.join(lines_set[0:9]) + '\t' + '\t'.join(lines_set_sub) + '\n')
    return [vcf_file_list,vcf_file_list_freq,vcf_file_list_vcf,vcf_file_POS_candidate]

def SNP_check_mapper(lines_set,vcf_file_list,vcf_file_list_freq,vcf_file_list_vcf,vcf_file_POS_candidate, min_cov_for_call_per_strand):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_frq2 = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    REF = lines_set[2]
    allels_set = [REF]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    allels_set += lines_set[3].split(',')
    Total_alleles = len(allels_set)
    Subdepth_all = lines_set[5]
    # Keep SNP that has higher than this quality
    Allels_frq = [0, 0, 0, 0, 0]
    Allels_frq_sub = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    Subdepth = Subdepth_all.split(';')
    for num_allels in range(0, min(len(Subdepth), Total_alleles)):
        allels = allels_set[num_allels]
        if allels in Allels:
            forward = float(Subdepth[num_allels].split(',')[0])
            reverse = float(Subdepth[num_allels].split(',')[1])
            Allels_frq[Allels[allels]] += forward + reverse
            Allels_frq_sub[Allels[allels] * 2] += forward
            Allels_frq_sub[Allels[allels] * 2 + 1] += reverse
        else:
            pass
    total_sub_depth = sum(Allels_frq)
    # find major alt and calculate frequency
    if sum(Allels_frq) > 0:
        Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
        MAF = Major_ALT[1] / total_sub_depth
        if Major_ALT[0] not in [REF, '-']:
            if (MAF == 1 and total_sub_depth >= min_cov) or (MAF >= min_maf_for_call and
                                                             total_sub_depth >= min_cov_for_call_per_strand * 2):
                if not pair_strand or \
                        sum(Subdepth_forward) >= min_cov_for_call_per_strand and \
                        sum(Subdepth_reverse) >= min_cov_for_call_per_strand:
                    temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
                    temp_snp_line_frq2.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub[1:]))
                    # a potential SNP
                    # calculate NS
                    gene_info = contig_to_gene(CHR, POS)
                    if False and gene_info != []: # disabled
                        Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
                        if Ref_seq_chr != 'None':
                            #  observed NS ratio calculated
                            temp_snp_line_NS = [Chr_gene, str(POS_gene), '']
                            if codon_start <= POS_gene - 1:
                                Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                                Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                                SNP_seq_chr = Ref_seq_chr
                                if len(Ref_seq_codon) == 3:
                                    Ref_seq_aa = translate(Ref_seq_codon)[0]
                                    temp_snp_line_AA += Ref_seq_aa
                                    SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, Major_ALT[0], Reverse_chr)
                                    SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                    SNP_seq_aa = translate(SNP_seq_codon)[0]
                                    temp_snp_line_AA += SNP_seq_aa
                                    temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                    temp_snp_line_NS[-1] += temp_NorS
                    # output lines and output major alt
                    temp_snp_line.append(CHR)
                    temp_snp_line.append(str(POS))
                    temp_snp_line.append(REF)
                    temp_snp_line.append(Major_ALT[0])
                    vcf_file_list.append(
                        '\t'.join(temp_snp_line) + '\t' + '\t'.join(temp_snp_line_frq2) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                            '1', 'PASS', '\t'.join(temp_snp_line_NS),
                            temp_snp_line_AA))
                    vcf_file_list_freq.append(
                        '\t'.join(temp_snp_line) + '\t\"%s\"\t\"%s\"\t%s\t%s\t%s\t%s\n' % (
                            'None',
                            '1',
                            'PASS', '\t'.join(temp_snp_line_NS),
                            temp_snp_line_AA, '\t'.join(temp_snp_line_frq)))
                    vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
                    vcf_file_list_vcf.append('\t'.join(lines_set) + '\n')
    return [vcf_file_list,vcf_file_list_freq,vcf_file_list_vcf,vcf_file_POS_candidate]

def load_sample(vcf_file):
    Sample_name = []
    for lines in open(os.path.join(args.smp,os.path.split(vcf_file)[-1].split(args.vcf)[0])+'.vcf.sh', 'r'):
        if ' -x ' in lines:
            # set database for bowtie
            database_file = lines.split(' -x ')[1].split(' -1 ')[0]
            Sample_name.append(os.path.split(lines.split(' -1 ')[1].split(' -2 ')[0])[-1])
            break
        if ' -ax sr -N 10 -p 0.9 -t 40 ' in lines:
            # set database for minimap2
            database_file = lines.split(' -ax sr -N 10 -p 0.9 -t 40 ')[1].split('.mmi')[0]
            Sample_name.append(os.path.split(lines.split('.mmi ')[1].split(' ')[0])[-1])
            break
        if 'bwa mem ' in lines:
            # set database for bwa
            database_file = lines.split('bwa mem ')[1].split(' ')[0]
            Sample_name.append(os.path.split(lines.split('bwa mem ')[1].split(' ')[1])[-1])
            break
    print(database_file, Sample_name)
    if database_file.split('.')[-1] != '.fna':
        # not gene file
        try:
            f1 = open(database_file + '.fna', 'r')
        except FileNotFoundError:
            print('%s -q -i %s -d %s.fna' % (args.pro, database_file, database_file))
            os.system('%s -q -i %s -d %s.fna' % (args.pro, database_file, database_file))
        database_file = database_file + '.fna'
    return [database_file,Sample_name]

def load_sample_mapper(vcf_file):
    Sample_name = []
    for lines in open(os.path.join(args.smp,os.path.split(vcf_file)[-1])+'.sh', 'r'):
        if ' --reference ' in lines:
            # set database
            database_file = lines.split(' --reference ')[1].split(' --queries')[0]
            Sample_name.append(os.path.split(lines.split(' --queries ')[1].split(' --out-vcf')[0])[-1])
            break
    print(database_file, Sample_name)
    if database_file.split('.')[-1] != '.fna':
        # not gene file
        try:
            f1 = open(database_file + '.fna', 'r')
        except FileNotFoundError:
            print('%s -q -i %s -d %s.fna' % (args.pro, database_file, database_file))
            os.system('%s -q -i %s -d %s.fna' % (args.pro, database_file, database_file))
        database_file = database_file + '.fna'
    return [database_file,Sample_name]

def contig_length(CHR):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        try:
            total_length = CHR.split('length_')[1].split('_cov')[0]
        except IndexError:
            total_length = 100000
    return int(total_length)

def SNP_filter(vcf_file,Sample_name,output_name,min_cov_for_call_per_strand, bowtie = True):
    vcf_file_list = []
    vcf_file_list_freq = []
    vcf_file_list_vcf = []
    vcf_file_POS_candidate = set()
    m = 0
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#") and not lines.startswith("CHR"):
            lines_set = lines.split('\n')[0].split('\t')
            if bowtie:
                if lines_set[3] != 'N' and '.' not in lines_set[4]:
                    # potential SNP
                    Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                    if (Depth >= min_cov):
                        # Remove candidate locations have lower than this coverage
                        m += 1
                        if m % 1000 == 0:
                            print('%s processed %s SNPs' % (datetime.now(), m))
                        vcf_file_list, vcf_file_list_freq, vcf_file_list_vcf, vcf_file_POS_candidate = \
                            SNP_check_fq(lines_set, vcf_file_list,
                                         vcf_file_list_freq, vcf_file_list_vcf, vcf_file_POS_candidate,
                                         min_cov_for_call_per_strand)
            else:  # mapper1
                if lines_set[2] not in ['N','-'] and lines_set[3] not in ['','-']:
                    # potential SNP
                    Depth = float(lines_set[4])
                    if (Depth >= min_cov):
                        # Remove candidate locations have lower than this coverage
                        m += 1
                        if m % 100000 == 0:
                            print('%s processed %s loci' % (datetime.now(), m))
                        vcf_file_list, vcf_file_list_freq, vcf_file_list_vcf, vcf_file_POS_candidate = \
                            SNP_check_mapper(lines_set, vcf_file_list,
                                             vcf_file_list_freq, vcf_file_list_vcf, vcf_file_POS_candidate,
                                             min_cov_for_call_per_strand)
    outputvcf(output_name,vcf_file_list_freq,vcf_file_list_vcf,Sample_name)
    print('%s finished output %s' % (datetime.now(), donor_species))

################################################### Main ########################################################
# run vcf filtering for bowtie
output_name = 'final'
if bowtie_filter:
    if args.cluster != 'None':
        allvcf_file = glob.glob(os.path.join(args.i, '%s*%s' % (args.cluster, args.vcf)))
    else:
        allvcf_file = glob.glob(os.path.join(args.i, '*%s' % (args.vcf)))
    print(allvcf_file)
    for vcf_file in allvcf_file:
        try:
            f1 = open(vcf_file + '.%s.snpfreq.txt' % (output_name), 'r')
        except IOError:
            donor_species = os.path.split(vcf_file)[-1]
            print('%s start processing %s' % (datetime.now(), donor_species))
            Ref_seq = dict()
            donor_species_sub = os.path.split(vcf_file)[-1].split(vcf_name)[0]
            database_file, Sample_name = load_sample(vcf_file)
            if Ref_seq == dict():
                Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
            print('%s load database %s' % (datetime.now(), database_file))
            Total = len(Sample_name)
            print('%s running %s genomes in %s' % (datetime.now(), Total, donor_species_sub))
            # SNP filtering
            print('%s start filtering SNPs %s' % (datetime.now(), donor_species_sub))
            SNP_filter(vcf_file, Sample_name, output_name,
                       min_cov_for_call_per_strand_round, True)

# run vcf filtering for mapper
if args.cluster != 'None':
    allvcf_file = glob.glob(os.path.join(args.i, '%s*mapper1.vcf' % (args.cluster)))
else:
    allvcf_file = glob.glob(os.path.join(args.i, '*mapper1.vcf'))

print(allvcf_file)
for vcf_file in allvcf_file:
    os.system('grep \';\' %s > %s.snp'%(vcf_file,vcf_file))
    try:
        f1 = open(vcf_file + '.%s.snpfreq.txt'%(output_name),'r')
    except IOError:
        donor_species = os.path.split(vcf_file)[-1]
        print('%s start processing %s' % (datetime.now(), donor_species))
        Ref_seq = dict()
        donor_species_sub = os.path.split(vcf_file)[-1].split(vcf_name)[0]
        database_file, Sample_name = load_sample_mapper(vcf_file)
        if Ref_seq == dict():
            Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
        print('%s load database %s' % (datetime.now(), database_file))
        Total = len(Sample_name)
        print('%s running %s genomes in %s' % (datetime.now(), Total, donor_species_sub))
        # SNP filtering
        print('%s start filtering SNPs %s' % (datetime.now(), donor_species_sub))
        SNP_filter(vcf_file + '.snp', Sample_name, output_name,
                   min_cov_for_call_per_strand_round, False)


