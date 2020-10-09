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
                      help="file extension of vcfs with only SNPs(.removerec.vcf)",
                      type=str, default='.removerec.vcf',
                      metavar='.filtered.vcf')
required.add_argument("-fq",
                      help="file extension of fastq #1 files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')
# optional output setup
optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
# requirement for software calling
optional.add_argument('-pro',
                          help="Optional: complete path to prodigal if not in PATH",
                          metavar="/usr/local/bin/prodigal",
                          action='store', default='prodigal', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
# set up path
Cluster = True
Tree = True
Cov_dis = 20
Cov_dis_overall = 1000 # calculate coverage per 1000 bp
input_script = args.s
vcf_name = '.raw.vcf'
ref_filename = '.fasta'
fastq_name = args.fq

################################################### Set up ########################################################
# set up steps
SNP_cluster = dict()
cluster_set = set()
reference_set = ['reference']

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

#end_cutoff = 50 # contig end no SNP calling
min_average_coverage_to_include_sample = 5 #Filter out samples have lower than that number of coverage
max_fraction_ambigious_samples = .4 #If more than % of the samples have ambiguous NT, discard the candidate location
min_median_coverage_position = 3 #Remove candidate locations have lower than this coverage
min_qual_for_call = 50 #Remove sample*candidate that has lower than this quality
min_maf_for_call = .9 #Remove sample*candidate
min_cov_for_call_per_strand = 3 #Remove sample*candidate

try:
    os.mkdir(args.i + '/tree')
except IOError:
    pass
input_script_sub = input_script + '/tree'
try:
    os.mkdir(input_script_sub)
except IOError:
    pass
################################################### Function ########################################################
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
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

def vcf_to_txt(lines,output_list):
    lines_set = lines.split('\n')[0].split('\t')
    if len(lines_set) >9:
        CHR = lines_set[0]
        POS = int(lines_set[1])
        temp_line = []
        temp_line.append(CHR)
        temp_line.append(str(POS))
        i = 9
        for Subdepth_all in lines_set[9:]:
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
            total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
            temp_line.append(str(total_sub_depth))
            i += 1
        output_list.append('\t'.join(temp_line)+'\n')
    else:
        print(lines)

def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = int(Subdepth[0]) + int(Subdepth[1])
    Subdepth_ALT = int(Subdepth[2]) + int(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return [allels_set[1],1]
    else:
        return [allels_set[0],0]

def outputvcf(output_name):
    vcf_file_filtered = open(vcf_file + '.%s.snpfreq.txt' % (output_name), 'w')
    vcf_file_filtered.write('CHR\tPOS\tMajor_ALT\tMinor_ALT\tGenome_set_noSNP\tGenome_set\tQuality\tGene\tGene_POS\tN_or_S\tAA_change\t%s\n'%('\t'.join(Sample_name))\
                            +''.join(vcf_file_list_freq))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.snp.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.samplename.txt' % (output_name), 'w')
    vcf_file_filtered.write('\t'.join(Sample_name) + '\n')
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.POS.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_POS))
    vcf_file_filtered.close()

def outputtree(output_name,runparsi = False):
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        newgenomename = genomename
        if len(genomename) > 8:
            newgenomename = genomename[0:4] + '_' + genomename[-4:]
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, SNP_alignment[genomename]))
            SNP_alignment_output_parsi.append('S%s    %s\n' % (newgenomename, SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.%s.vcf' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list_vcf))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.fasta' % (output_name), 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    if runparsi:
        vcf_file_filtered = open(vcf_file + '.%s.parsi.fasta' %(output_name), 'w')
        vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
        vcf_file_filtered.close()
        run_parsi(vcf_file + '.%s.parsi.fasta' % (output_name))

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

def outputcovwindow(MV_cov):
    cov_output = []
    for lines in open(vcf_file.split(vcf_name)[0] + vcf_name, 'r'):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = int(lines_set[1])
            if POS%Cov_dis_overall == 1:
                # subset this CHR
                temp_cov = []
                lines_set_sub = lines_set[9:]
                for Subdepth_all in lines_set_sub:
                    Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                    total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                    temp_cov.append(total_sub_depth)
                cov_output.append('%s\t%s\t%s\t%s\n'%(CHR,POS,statistics.mean(temp_cov),
                                                      '\t'.join(str(cov) for cov in temp_cov)))
    vcf_file_filtered = open(MV_cov, 'w')
    vcf_file_filtered.write('CHR\tPOS\tavg_depth\t%s\n'%('\t'.join(Sample_name))+''.join(cov_output))
    vcf_file_filtered.close()
    return 'finished sampling coverage for %s'%(vcf_file)

def coverage_sample(MV_cov):
    Depth = dict()
    Badsample = []
    try:
        f1 = open(MV_cov, 'r')
    except IOError:
        outputcovwindow(MV_cov)
    for lines in open(MV_cov,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        if not lines.startswith("CHR"):
            i = 3
            for samplename in Depth:
                Depth[samplename].append(int(lines_set[i]))
                i+=1
        else:
            for samplename in lines_set[3:]:
                Depth.setdefault(samplename,[])
    for samplename in Depth:
        depth_mean = statistics.mean(Depth[samplename])
        if depth_mean < min_average_coverage_to_include_sample:
            Badsample.append(samplename)
            print('excluding sample %s for low coverage %s'%(samplename,depth_mean))
    return Badsample

def SNP_check_all_fq(lines_set,CHR_old,POS_old,reference_name,Total,Badsamplenum):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_frq2 = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    SNP = set()
    NOSNP = set()
    BADSNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    lines_set_sub = lines_set[9:]
    if Badsamplenum!= []:
        lines_set_sub = [lines_set[i] for i in range(9,len(lines_set)) if i not in Badsamplenum]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    SNP_quality = lines_set[5]
    if '.' not in lines_set[4]:
        allels_set += lines_set[4].split(',')
    Total_alleles = len(allels_set)
    genome_order = 0
    Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
    REF, REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
    sample_num = 9
    for Subdepth_all in lines_set_sub:
        Qual = float(Subdepth_all.split(':')[0])
        genome_order += 1
        Allels_frq = [0, 0, 0, 0]
        Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        Allels_frq_sub[0] = Qual
        SNP_seq.append(REF)  # set as reference
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
        total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
        Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
        Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
        for num_allels in range(0, min(len(Subdepth), Total_alleles)):
            allels = allels_set[num_allels]
            Subdepth_alleles = int(Subdepth[num_allels])
            if allels in Allels:
                Allels_frq[Allels[allels]] += Subdepth_alleles
                Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_forward[num_allels])
                Allels_frq_sub[Allels[allels] * 2 + 2] += int(Subdepth_reverse[num_allels])
            else:
                pass
        if Qual >= min_qual_for_call:
            # Keep sample*candidate that has lower than this quality
            # find major alt and calculate frequency
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            if total_sub_depth >= min_cov_for_call_per_strand * 2 and Major_ALT[1]/total_sub_depth >= min_maf_for_call:
                # Keep sample*candidate
                if Major_ALT[0] != REF:
                    SNP_seq[-1] = Major_ALT[0]
                    SNP.add(genome_order)
                else:
                    NOSNP.add(genome_order)
            else:
                BADSNP.add(genome_order)
                NOSNP.add(genome_order)
        else:
            BADSNP.add(genome_order)
            NOSNP.add(genome_order)
        temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
        temp_snp_line_frq2.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub[1:]))
        sample_num += 1
    if SNP!= set() and NOSNP != set() and len(BADSNP) <= max_fraction_ambigious_samples*Total:
        # If less than X% bad samples have ambiguous NT, keep the candidate location
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
                        ALT_set = allels_set
                        for ALT in ALT_set:
                            if ALT != REF:
                                SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                SNP_seq_aa = translate(SNP_seq_codon)[0]
                                temp_snp_line_AA += SNP_seq_aa
                                temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                temp_snp_line_NS[-1] += temp_NorS
        # output lines and output major alt
        temp_snp_line_pass = 'PASS'
        if CHR == CHR_old:
            # same CHR
            POS_DIS = abs(POS - POS_old)
            vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
        else:
            # diff CHR first SNP
            vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
        POS_old = POS
        CHR_old = CHR
        if len(SNP) > len(NOSNP):
            NOSNP2 = NOSNP
            NOSNP = SNP
            SNP = NOSNP2
            REF = allels_set[1-REF_where]
        temp_snp_line.append(CHR)
        temp_snp_line.append(str(POS))
        temp_snp_line.append(REF)
        temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
        vcf_file_list.append(
            '\t'.join(temp_snp_line) + '\t' + '\t'.join(temp_snp_line_frq2) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in sorted(SNP)), temp_snp_line_pass, '\t'.join(temp_snp_line_NS),
                temp_snp_line_AA))
        vcf_file_list_freq.append(
            '\t'.join(temp_snp_line) + '\t\"%s\"\t\"%s\"\t%s\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in sorted(NOSNP)),';'.join(str(genome_order) for genome_order in sorted(SNP)),
                SNP_quality, '\t'.join(temp_snp_line_NS),
                temp_snp_line_AA, '\t'.join(temp_snp_line_frq)))
        vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
        vcf_file_list_vcf.append('\t'.join(lines_set[0:9]) + '\t' + '\t'.join(lines_set_sub) + '\n')
        i = 9
        j = 0
        SNP_alignment[reference_name] += REF
        for genomename in SNP_alignment:
            if genomename != reference_name:
                SNP_alignment[genomename] += SNP_seq[j]
                j += 1
                i += 1
    return [CHR_old,POS_old]

def contig_end(CHR,POS):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        try:
            total_length = CHR.split('length_')[1].split('_cov')[0]
        except IndexError:
            return False
    total_length = int(total_length)
    if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
        return True
    else:
        return False

def run_parsi(a_parsi_file):
    if 'RAxML_parsimonyTree' not in a_parsi_file:
        os.system('rm -rf %s %s' % (a_parsi_file + '.out.txt',
                                    a_parsi_file + '.out.tree'))
        SNP_tree_cmd3 = ('%s\n5\nV\n1\ny\n' % (a_parsi_file))
        f1 = open(os.path.join(input_script_sub, 'parsi.optionfile.txt'), 'w')
        f1.write(SNP_tree_cmd3)
        f1.close()
        os.system('rm -rf outfile outtree')
        os.system('dnapars < %s/parsi.optionfile.txt > %s/%s.parsi.output\n' % (
            input_script_sub, input_script_sub, os.path.split(a_parsi_file)[-1]))
        os.system('mv outfile %s' % (a_parsi_file + '.out.txt'))
        os.system('mv outtree %s' % (a_parsi_file + '.out.tree'))

def outputtree_parsi(output_file):
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        newgenomename = genomename
        if len(genomename) > 8:
            newgenomename = genomename[0:4] + '_' + genomename[-4:]
        if seq_len > 0:
            SNP_alignment_output_parsi.append('S%s    %s\n' % (newgenomename, SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(output_file, 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()
    vcf_file_filtered = open(output_file + '.sum.txt', 'w')
    vcf_file_filtered.write(''.join(alloutput))
    vcf_file_filtered.close()

def ALT_freq_sub(Allels_subcount):
    Allels_count = [Allels_subcount[0] + Allels_subcount[1],
                    Allels_subcount[2] + Allels_subcount[3],
                    Allels_subcount[4] + Allels_subcount[5],
                    Allels_subcount[6] + Allels_subcount[7]]
    return Allels_order[Allels_count.index(max(Allels_count))]

def findSNP_sample(SNP_list,samplelist):
    samplelist_genotype = []
    genotype_set = []
    i = 0
    for SNP in SNP_list:
        SNP_sub = [int(i) for i in SNP.split(';')[1:]]
        samplelist_genotype.append(samplelist[i])
        if sum(SNP_sub) > 0:
            ALT = ALT_freq_sub(SNP_sub)
            genotype_set.append(ALT)
        else:
            # EMPTY
            genotype_set.append('-')
        i += 1
    return [samplelist_genotype,genotype_set]

def find_major(genotype_set):
    allALT = list(set(genotype_set))
    if '-' in allALT:
        allALT.remove('-')
    countALT = []
    for ALT in allALT:
        countALT.append(genotype_set.count(ALT))
    Major = allALT[countALT.index(max(countALT))]
    return [Major,','.join([i for i in allALT if i != Major])]

################################################### Main ########################################################
# run vcf filtering
allvcf_file = glob.glob(os.path.join(args.i,'*%s'%(args.vcf)))
print(allvcf_file)
reference_name = reference_set[0]
output_name = 'final'
for vcf_file in allvcf_file:
    filesize = 0
    try:
        filesize = int(os.path.getsize(vcf_file + '.%s.snpfreq.txt' % (output_name)))
    except FileNotFoundError:
        pass
    if filesize == 0:
        donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0]
        vcf_file_list = []
        vcf_file_list_freq = []
        vcf_file_list_vcf = []
        Sample_name = []
        vcf_file_POS = []
        vcf_file_POS_candidate = set()
        Ref_seq = dict()
        Mapping = dict()
        Mapping_loci = dict()
        CHR_old = ''
        POS_old = 0
        # filter out samples with low coverage
        print('remove samples with low coverage %s'%(donor_species))
        Badsample = coverage_sample(vcf_file.split(vcf_name)[0] + vcf_name + '.filtered.cov.MW.txt')
        Badsamplenum = []
        for lines in open(vcf_file.split(vcf_name)[0] + vcf_name, 'r'):
            if lines.startswith('##bcftoolsCommand=mpileup '):
                # setup samples
                sample_set = lines.split(ref_filename + ' ')[1].split('\n')[0].split(' ')
                samplenum = 9
                for samples in sample_set:
                    genomename = os.path.split(samples)[-1].split(fastq_name)[0].split('all')[0].split('.sorted.bam')[0]
                    if genomename not in Badsample:
                        Sample_name.append(genomename.replace('.', ''))
                    else:
                        Badsamplenum.append(samplenum)
                    samplenum += 1
            if lines.startswith('##reference=file:'):
                database_file = lines.split('##reference=file:')[1].split('\n')[0]
                break
        print('running %s' % (donor_species))
        if database_file.split('.')[-1] != '.fna':
            # not gene file
            try:
                f1 = open(database_file + '.fna', 'r')
            except FileNotFoundError:
                os.system('%s -q -i %s -d %s.fna' % (args.pro, database_file, database_file))
            database_file = database_file + '.fna'
        Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
        vcf_file_list = []
        vcf_file_list_freq = []
        vcf_file_list_vcf = []
        vcf_file_POS = []
        vcf_file_POS_candidate = set()
        SNP_alignment = dict()
        SNP_alignment.setdefault(reference_name, '')
        CHR_old = ''
        POS_old = 0
        Total = len(Sample_name)
        for genomename in Sample_name:
            SNP_alignment.setdefault(genomename, '')
        for lines in open(vcf_file + '', 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\n')[0].split('\t')
                CHR = lines_set[0]
                POS = int(lines_set[1])
                Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                if Depth/Total >= min_median_coverage_position:
                    # Remove candidate locations have lower than this coverage
                    CHRPOS = '%s\t%s' % (CHR, POS)
                    CHR_old, POS_old = SNP_check_all_fq(lines_set,
                                                            CHR_old, POS_old, reference_name,Total,Badsamplenum)
        outputvcf(output_name)
        outputtree(output_name,True)

# merge parsi file
allvcf_file = glob.glob(os.path.join(args.i,'*%s.%s.snpfreq.txt'%(args.vcf,output_name)))
for vcf_file in allvcf_file:
    donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
    print('processing %s %s' % (donor_species, vcf_file))
    parsi_output = os.path.join(os.path.split(vcf_file)[0],'%s.all.parsi.fasta'%(donor_species))
    filesize = 0
    try:
        filesize = int(os.path.getsize(parsi_output))
    except FileNotFoundError:
        pass
    if filesize == 0:
        print('start processing %s' % (donor_species))
        allvcf_filesub = [vcf_file for vcf_file in allvcf_file if donor_species in vcf_file]
        vcf_sample = dict()
        SNP_alignment = dict()
        allsamplename = []
        # load all samplename
        for a_vcf_file in allvcf_filesub:
            samplefile = a_vcf_file.replace('.snpfreq.txt','.samplename.txt')
            vcf_sample.setdefault(a_vcf_file,[])
            for lines in open(samplefile,'r'):
                sample_set = lines.split('\n')[0].split('\t')
                vcf_sample[a_vcf_file] += sample_set
                allsamplename += sample_set
                for genomename in sample_set:
                    SNP_alignment.setdefault(genomename, '')
        # load all SNPs and POSs
        SNP_info = dict()
        SNP_sample = dict()
        for a_vcf_file in allvcf_filesub:
            for lines in open(a_vcf_file.replace('.snpfreq.txt','.vcf')):
                lines_set = lines.split('\n')[0].split('\t')
                CHR, POS, NOTUSE, Ref, ALT = lines_set[0:5]
                CHRPOS = '%s\t%s' % (CHR, POS)
                SNP_info.setdefault(CHRPOS, [Ref, '', set()])
            for lines in open(a_vcf_file, 'r'):
                if not lines.startswith('CHR'):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR, POS, Ref, ALT = lines_set[0:4]
                    ALT = ALT.split(',')
                    CHRPOS = '%s\t%s' % (CHR, POS)
                    Gene_info = '\t'.join(lines_set[7:9])
                    N_or_S = lines_set[9]
                    SNP_info[CHRPOS][1] = Gene_info
                    SNP_sample.setdefault(CHRPOS, [[], []])
                    samplelist_genotype, genotype_set = findSNP_sample(lines_set[11:], vcf_sample[a_vcf_file])
                    for minor_ALT in list(set(genotype_set)):
                        if minor_ALT != Ref and minor_ALT != '-':
                            SNP_info[CHRPOS][-1].add(N_or_S[ALT.index(minor_ALT)])
                    SNP_sample[CHRPOS][0] += samplelist_genotype
                    SNP_sample[CHRPOS][1] += genotype_set
        # summarize parsi seq
        alloutput = []
        alloutput.append('CHR\tPOS\tMajor_ALT\tMinor_ALT\tN_or_S\tGene\tGene_POS\tGenome_set_noSNP\tGenome_SNP_set\t%s\n'%('\t'.join(allsamplename)))
        for CHRPOS in SNP_info:
            Ref, Gene_info,N_or_S = SNP_info[CHRPOS]
            samplelist_genotype, genotype_set = SNP_sample[CHRPOS]
            Major, Minor = find_major(genotype_set + [Ref]* (len(SNP_alignment)-len(samplelist_genotype)))
            Genome_set_noSNP = []
            Genome_set_SNP = []
            temp_line = '%s\t%s\t%s\t%s\t%s'%(CHRPOS,Major,Minor,','.join(list(N_or_S)),Gene_info)
            temp_line2 = ''
            i = 1
            for genomename in SNP_alignment:
                if genomename in samplelist_genotype:
                    ALT = genotype_set[samplelist_genotype.index(genomename)]
                else:
                    # not significantly different than Ref
                    ALT = Ref
                if ALT in [Major, '-']:
                    ALT = Major
                    Genome_set_noSNP.append(str(i))
                else:
                    Genome_set_SNP.append(str(i))
                temp_line2 += '\t%s' % (ALT)
                SNP_alignment[genomename] += ALT
                i += 1
            alloutput.append(temp_line + '\t%s\t%s'%(';'.join(Genome_set_noSNP),';'.join(Genome_set_SNP)) + temp_line2 + '\n')
        outputtree_parsi(parsi_output)
        print('finished processing %s' % (donor_species))
        print('run parsi tree for %s' % (donor_species))
        run_parsi(parsi_output)

# run parsi tree
os.system('mv %s/*.parsi.fasta.out* %s/tree/' % (
    args.i, args.i))


