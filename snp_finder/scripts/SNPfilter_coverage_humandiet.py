################################################### END ########################################################
################################################### SET PATH ########################################################
# After round 4 filter results of WGS
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from statistics import median
from statistics import mean
from datetime import datetime
import numpy as np
# set up path
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path to all vcf files",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/xq_data/merge/',
                      metavar='input/')
# optional output setup
optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='/scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/xq_data/',
                      metavar='scripts/')

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
################################################### Set up ########################################################
# set up steps
SNP_cluster = dict()
cluster_set = set()
separate_donor_genome = []
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels['-']=4
Allels['+A']=5
Allels['+T']=6
Allels['+G']=7
Allels['+C']=8
Allels_order = ['A','T','G','C','-','+A','+T','+G','+C']
nomal_alleles = ['A','T','G','C']

min_qual_for_call = 40 #Remove sample*candidate that has lower than this quality
neighbour_cov_range = 100 # 100bp
################################################### Function ########################################################
# set up functions
def all_ALT(Allels_count,REF,total_depth):
    allALT = []
    ALT_freq = 0
    ALT_set = dict()
    ALT_frq_set = set()
    Depth_fm = sum([Allels_count[i * 4] for i in range(0,9)])
    Depth_fe = sum([Allels_count[i * 4 + 1] for i in range(0,9)])
    Depth_rm = sum([Allels_count[i * 4 + 2] for i in range(0,9)])
    Depth_re = sum([Allels_count[i * 4 + 3] for i in range(0,9)])
    ALT_freq_fm = 0
    ALT_freq_fe = 0
    ALT_freq_rm = 0
    ALT_freq_re = 0
    for alleles in range(0, 9):
        forwardm = Allels_count[alleles * 4]
        forwarde = Allels_count[alleles * 4 + 1]
        reversem = Allels_count[alleles * 4 + 2]
        reversee = Allels_count[alleles * 4 + 3]
        ALT_frq = forwardm + forwarde + reversem + reversee
        if ALT_frq >= depth_ALT_cutoff:
            ALT_set.setdefault(ALT_frq, [])
            ALT_set[ALT_frq].append(alleles)
            ALT_frq_set.add(ALT_frq)
            # SNP
            ALT = Allels_order[alleles]
            if ALT != REF:
                if ALT_frq / total_depth >= ALT_freq_cutoff:
                    qualify_allele = True
                    # ALT freq cutoff
                    if reversem + forwardm < depth_ALT_cutoff:
                        # low middle reads
                        expectm = 0
                        if Depth_fe + Depth_re > 0:
                            expectm = (Depth_fm + Depth_rm) / (Depth_fe + Depth_re) * (reversee + forwarde)
                        if reversem + forwardm < int(expectm * 0.75):
                            # not qualified
                            qualify_allele = False
                            # end freq cutoff
                    if qualify_allele:
                        ALT_freq += ALT_frq
                        ALT_freq_fm += forwardm
                        ALT_freq_fe += forwarde
                        ALT_freq_rm += reversem
                        ALT_freq_re += reversee
    # order by major to minor
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            ALT = Allels_order[alleles]
            allALT.append(ALT)
    return [allALT,ALT_freq,Depth_fm,Depth_fe,Depth_rm,Depth_re,ALT_freq_fm,ALT_freq_fe,ALT_freq_rm,ALT_freq_re]

def outputvcf(output_name,vcf_file_list_freq,vcf_file_list_freq_snp,Sample_name):
    vcf_file_filtered = open(vcf_file + '.%s.allfreq.txt' % (output_name), 'w')
    vcf_file_filtered.write(
        'CHR\tPOS\tREF\tALT\tDepth\tDepth_median\tDepth_ALT\tDepth_fm\tDepth_fe\tDepth_rm\tDepth_re\tDepth_ALT_fm\tDepth_ALT_fe\tDepth_ALT_rm\tDepth_ALT_re\tQuality\tGene\tGene_POS\tN_or_S\tAA_change\t'+
                            'reads_sub\n'
        + ''.join(vcf_file_list_freq))
    vcf_file_filtered.close()

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

def SNP_check_fq(lines_set,vcf_file_list_freq, vcf_file_list_freq_snp,depth_ALT_cutoff, depth_cutoff,median_depth,oldPOS,oldPOS_out,bowtie=True):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    if bowtie:
        REF = lines_set[3]
        allels_set = [REF]
        CHR = lines_set[0]
        POS = int(lines_set[1])
        SNP_quality = lines_set[5]
        allels_set += lines_set[4].split(',')
        allels_set = [x for x in allels_set if x in Allels]
        Total_alleles = len(allels_set)
        Subdepth_all = lines_set[9]
        Qual = float(SNP_quality)
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
        Subdepth_forward = [int(i) for i in Subdepth_all.split(':')[-3].split(',')]
        Subdepth_reverse = [int(i) for i in Subdepth_all.split(':')[-2].split(',')]
    else:
        REF = lines_set[2]
        allels_set = [REF]
        CHR = lines_set[0].split(' ')[0]
        POS = int(lines_set[1])
        allels_set += lines_set[3].split(',')
        Total_alleles = len(allels_set)
        Subdepth_all = lines_set[5]
        # Keep SNP that has higher than this quality
        Subdepth = Subdepth_all.split(';')
        Qual = min_qual_for_call  # not quality cutoff for mapper
    if POS == -1 and oldPOS_out!=[]:
        # insertion
        Allels_frq_sub,allels_setall,REF = oldPOS_out
        allels_set = [REF]
        allels_set += ['+'+i for i in lines_set[3].split(',') if i != '']
        Total_alleles = len(allels_set)
    elif oldPOS!= 0:
        # output last POS
        Allels_frq_sub,allels_setall,REFold = oldPOS_out
        total_depth = sum(Allels_frq_sub)
        # find major alt and calculate frequency
        if total_depth > 0:
            allALT, ALT_freq, Depth_fm,Depth_fe,Depth_rm,Depth_re,ALT_freq_fm,ALT_freq_fe,ALT_freq_rm,ALT_freq_re = all_ALT(Allels_frq_sub, REFold,total_depth)
            temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
            # calculate NS
            gene_info = contig_to_gene(CHR, oldPOS)
            if gene_info != []:
                Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
                if Ref_seq_chr != 'None':
                    #  observed NS ratio calculated
                    temp_snp_line_NS = [Chr_gene, str(POS_gene), '']
                    if codon_start <= POS_gene - 1:
                        Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REFold, Reverse_chr)
                        Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                        SNP_seq_chr = Ref_seq_chr
                        if len(Ref_seq_codon) == 3:
                            Ref_seq_aa = translate(Ref_seq_codon)[0]
                            temp_snp_line_AA += Ref_seq_aa
                            for ALT in allels_setall:
                                if ALT != REFold and ALT in nomal_alleles:
                                    SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                    SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                    SNP_seq_aa = translate(SNP_seq_codon)[0]
                                    temp_snp_line_AA += SNP_seq_aa
                                    temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                    temp_snp_line_NS[-1] += temp_NorS
            # output lines and output major alt
            temp_snp_line.append(CHR)
            temp_snp_line.append(str(oldPOS))
            temp_snp_line.append(REFold)
            temp_snp_line.append(','.join(allALT))
            vcf_file_list_freq.append(
                '\t'.join(
                    temp_snp_line) + '\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%s\t%s\t%s\t%s\n' % (
                    total_depth, median_depth, ALT_freq, Depth_fm, Depth_fe, Depth_rm, Depth_re, ALT_freq_fm,
                    ALT_freq_fe, ALT_freq_rm, ALT_freq_re,
                    Qual, '\t'.join(temp_snp_line_NS),
                    temp_snp_line_AA, '\t'.join(temp_snp_line_frq)))
            if False and total_depth >= depth_cutoff and ALT_freq >= depth_ALT_cutoff and allALT != [REFold]:
                # for total depth between depth_ratio_cutoff1 and depth_ratio_cutoff2 * median depth of neighbour
                # a SNP
                vcf_file_list_freq_snp.append(
                    '\t'.join(
                        temp_snp_line) + '\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%s\t%s\t%s\t%s\n' % (
                        total_depth, median_depth, ALT_freq, Depth_fm, Depth_fe, Depth_rm, Depth_re, ALT_freq_fm,
                        ALT_freq_fe, ALT_freq_rm, ALT_freq_re,
                        Qual, '\t'.join(temp_snp_line_NS),
                        temp_snp_line_AA, '\t'.join(temp_snp_line_frq)))
        # update OLD POS
        oldPOS = POS
        Allels_frq_sub = [0] * (4 * len(Allels_order))
        oldPOS_out = [Allels_frq_sub,allels_setall,REF]
    else:
        # update OLD POS
        oldPOS = POS
        Allels_frq_sub = [0] * (4 * len(Allels_order))
        oldPOS_out = [Allels_frq_sub,set(),'']
    if Qual >= min_qual_for_call:
        # Keep SNP that has higher than this quality
        for num_allels in range(0, min(len(Subdepth), Total_alleles)):
            allels = allels_set[num_allels]
            if bowtie:
                forwardm = Subdepth_forward[num_allels]
                forwarde = 0
                reversem = Subdepth_reverse[num_allels]
                reversee = 0
                #alldepth = forwardm + reversem
            else:
                forwardm, forwarde = Subdepth[num_allels].split(',')[0].split('/')
                reversem, reversee = Subdepth[num_allels].split(',')[1].split('/')
                forwardm = float(forwardm)
                forwarde =  float(forwarde)
                reversem = float(reversem)
                reversee = float(reversee)
                #alldepth = forwardm + forwarde + reversem + reversee
            Allels_frq_sub[Allels[allels] * 4] += forwardm
            Allels_frq_sub[Allels[allels] * 4 + 1] += forwarde
            Allels_frq_sub[Allels[allels] * 4 + 2] += reversem
            Allels_frq_sub[Allels[allels] * 4 + 3] += reversee
        allels_setall = oldPOS_out[-1]
        allels_setall = list(allels_setall)
        allels_setall += allels_set
        allels_setall = set(allels_setall)
        oldPOS_out = [Allels_frq_sub,allels_setall,REF]
    return [vcf_file_list_freq,vcf_file_list_freq_snp,oldPOS,oldPOS_out]

def load_sample(vcf_file):
    Sample_name = []
    for lines in open(vcf_file, 'r'):
        if lines.startswith('##reference=file:'):
            # set database
            database_file = lines.split('##reference=file:')[1].split('\n')[0]
        if lines.startswith('#CHROM'):
            Sample_name.append(os.path.split(lines.split('\t')[9].split('\n')[0])[-1].split('.sorted.bam')[0])
            break
    if database_file.split('.')[-1] != '.fna':
        # not gene file
        try:
            f1 = open(database_file + '.fna', 'r')
        except FileNotFoundError:
            os.system('%s -q -i %s -d %s.fna' % (args.pro, database_file, database_file))
        database_file = database_file + '.fna'
    return [database_file,Sample_name]

def SNP_filter(vcf_file,Sample_name,database):
    alldepth = [[],0]
    # read all depth
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#") and not lines.startswith("CHR"):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0].split(' ')[0]
            Depth = float(lines_set[4])
            Qual = min_qual_for_call
            if Qual >= min_qual_for_call:
                if CHR in ['NODE_7_length_83618_cov_574.879053', 'NODE_21_length_39475_cov_506.843088',
                                    'NODE_13_length_54051_cov_1504.877960']:
                    POS = int(lines_set[1])
                    gene_info = contig_to_gene(CHR, POS)
                    if gene_info != []:
                        Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
                        if Chr_gene in ['NODE_7_length_83618_cov_574.879053_45', 'NODE_21_length_39475_cov_506.843088_11',
                                        'NODE_13_length_54051_cov_1504.877960_15']:
                            print(Chr_gene, POS_gene)
                            # target gene
                            REF = lines_set[2]
                            allels_set = lines_set[3].split(',')
                            temp_snp_line_NS = ['', '']
                            if allels_set != []:
                                #  observed NS ratio calculated
                                if codon_start <= POS_gene - 1:
                                    Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                                    Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                                    SNP_seq_chr = Ref_seq_chr
                                    if len(Ref_seq_codon) == 3:
                                        Ref_seq_aa = translate(Ref_seq_codon)[0]
                                        temp_snp_line_NS[0] += Ref_seq_aa
                                        for ALT in allels_set:
                                            if ALT != REF and ALT in nomal_alleles:
                                                SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                                SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                                SNP_seq_aa = translate(SNP_seq_codon)[0]
                                                temp_snp_line_NS[0] += SNP_seq_aa
                                                temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                                temp_snp_line_NS[1] += temp_NorS
                            alloutput2.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(
                                Sample_name[0],database,Chr_gene,POS_gene,REF,lines_set[3],temp_snp_line_NS[0],temp_snp_line_NS[1]))
                alldepth[0].append(Depth)
                if lines_set[2] != '-':
                    alldepth[1]+=1
    genomedepth = alldepth[0]
    alloutput.append('%s\t%s\t%.1f\t%.1f\t%s\t%s\n'%(
        Sample_name[0],database,sum(genomedepth),mean(genomedepth),'\t'.join(
                [str(i) for i in np.quantile(genomedepth, [0.1, 0.5, 0.9])]
            ),
            alldepth[1]))


def load_sample_mapper(vcf_file):
    Sample_name = []
    sample = os.path.split(vcf_file)[-1].split('.mapper')[0].split('_')[-1]
    Sample_name.append(sample)
    database = '_'.join(os.path.split(vcf_file)[-1].split('.mapper')[0].split('_')[:-1])
    if database == 'BL_CL27':
        database_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/BL_clustercluster27/BL_clustercluster27.all.spades2.fasta.noHM.fasta'
    if database == 'BL_CL33':
        database_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/BL_clustercluster33/BL_clustercluster33.all.spades2.fasta.noHM.fasta'
    if database == 'BiPs_CL1':
        database_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/BiPs_clustercluster1/BiPs_clustercluster1.all.spades2.fasta.noHM.fasta'
    if database_file.split('.')[-1] != '.fna':
        # not gene file
        try:
            f1 = open(database_file + '.fna', 'r')
        except FileNotFoundError:
            os.system('%s -q -i %s -d %s.fna' % (args.pro, database_file, database_file))
        database_file = database_file + '.fna'
    return [database_file,Sample_name,database]

################################################### Main ########################################################

# run vcf filtering for mapper
depth_ratio_cutoff1 = 0 # for SNP at least depth_ratio_cutoff1 * median depth of neighbour_cov_range bp
depth_ALT_cutoff = 0  # for each REF + ALT freq
ALT_freq_cutoff =0 # for SNP depth of one ALT / total depth
output_name = 'final'
allvcf_file = glob.glob(os.path.join(args.i, '*mapper1.vcf'))
print(allvcf_file)
alloutput = []
alloutput.append('sample\tgenome\tsumdepth\tmean_depth\tdepth_0.1\tmedian_depth\tdepth_0.9\tgenome_coverage\n')
alloutput2 = []
alloutput2.append('sample\tgenome\tgene\tgene_POS\tREF\tALT\tAAchange\tN_S\n')

for vcf_file in allvcf_file:
    try:
        f1 = open(vcf_file + '.%s.allfreqsum.txt' % (output_name), 'r')
    except IOError:
        sample = os.path.split(vcf_file)[-1].split('mapper1.vcf')[0]
        print('%s load database for %s' % (datetime.now(), sample))
        Ref_seq = dict()
        database_file, Sample_name,database = load_sample_mapper(vcf_file)
        if Ref_seq == dict():
            Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
        # SNP filtering
        print('%s start filtering SNPs %s' % (datetime.now(), sample))
        SNP_filter(vcf_file, Sample_name, database)
        print('%s finished output %s' % (datetime.now(), sample))

f1 = open(args.i + 'all.allfreqsum.txt', 'w')
f1.write(
         ''.join(alloutput))
f1.close()
f1 = open(args.i + 'all.target.gene.txt', 'w')
f1.write(
         ''.join(alloutput2))
f1.close()