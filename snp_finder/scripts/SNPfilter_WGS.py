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
genome_root = args.i + '/round%s'%(Round)
genome_dir = glob.glob(args.i + '/round%s/*'%(Round))
output_dir = args.o + '/vcf_round%s/bwa/0/'%(Round)
output_dir_merge = args.o +'/vcf_round%s/merge'%(Round)
vcf_name = '.all.flt.snp.vcf'
ref_filename = '.all.spades%s.fasta'%(Round)
fasta_name = '.fasta.corrected.fasta'
fastq_name = '.sorted.bam'
deleting_file = []
Species_replace = dict()
if Round == 4:
    genome_root = args.i + '/round*'
    output_dir_merge = argso + '/vcf_round%s/merge' % (Round)

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

def vcf_to_txt(lines,output_list,cluster_sub=[]):
    lines_set = lines.split('\n')[0].split('\t')
    if len(lines_set) >9:
        CHR = lines_set[0]
        POS = int(lines_set[1])
        temp_line = []
        temp_line.append(CHR)
        temp_line.append(str(POS))
        i = 9
        for Subdepth_all in lines_set[9:]:
            if (cluster_sub==[] and i not in deleting_set) or i in cluster_sub:
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
    vcf_file_filtered = open(vcf_file + '.%s.snp.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.samplename.txt' % (output_name), 'w')
    vcf_file_filtered.write('\t'.join(Sample_name) + '\n')
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.POS.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_POS))
    vcf_file_filtered.close()

def outputcov(output_name,vcf_file_POS_candidate,cluster_sub=[]):
    if len(vcf_file_list) > 0:
        vcf_file_POS_candidate = '\n'.join(vcf_file_POS_candidate)
        vcf_file_POS_candidate_output = ('%s' % (vcf_file_POS_candidate))
        f1 = open(os.path.join(input_script,'grep.temp.txt'),'w')
        f1.write(vcf_file_POS_candidate_output)
        f1.close()
        os.system('grep -%s -f %s %s --no-group-separator > %s'% (
            Cov_dis, os.path.join(input_script,'grep.temp.txt'),
            vcf_file.split('.flt.snp.vcf')[0] + '.raw.vcf',
            vcf_file + '.%s.cov.temp' % (output_name)))
        os.system('cat %s | sort | uniq > %s' % (
            vcf_file + '.%s.cov.temp' % (output_name),
            vcf_file + '.%s.uniqcov.temp' % (output_name)))
        for lines in open(vcf_file + '.%s.uniqcov.temp' % (output_name), 'r'):
            if not lines.startswith("#"):
                vcf_to_txt(lines, cov_file_list,cluster_sub)
        os.system('rm -rf %s %s' % (vcf_file + '.%s.cov.temp' % (output_name),
                                    vcf_file + '.%s.uniqcov.temp' % (output_name)))
        vcf_file_filtered = open(vcf_file + '.%s.cov.txt' % (output_name), 'w')
        vcf_file_filtered.write(''.join(cov_file_list))
        vcf_file_filtered.close()

def outputtree(output_name):
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, SNP_alignment[genomename]))
            SNP_alignment_output_parsi.append('S%s    %s\n' % (genomename[-8:], SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.%s.vcf' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list_vcf))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.fasta' % (output_name), 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.parsi.fasta' %(output_name), 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
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

def SNP_distance_correct(distanceNEW, distanceREF, SNPREF):
    return int((distanceNEW/distanceREF)*SNPREF)

def find_neighbor(Cluster_SNP,neighbor,Cluster_SNP_set,cluster,Cluster_SNP_set_added):
    if neighbor != []:
        for record_name in neighbor:
            if record_name not in Cluster_SNP_set_added:
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    subneighbor = Cluster_SNP.get(record_name,[])
                    find_neighbor(Cluster_SNP, subneighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)

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

def SNP_check_all(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name,SNP_presence_cutoff,SNP_presence_sample_cutoff,cluster_sub=[]):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    Total_qualify = 0
    Total_qualify_SNP = 0
    Total_qualify_notSNP = 0
    Total_unqualify_alt_freq = 0
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    if cluster_sub!= []:
        lines_set_sub = [lines_set[i] for i in cluster_sub]
        Total_subsample = len(cluster_sub)
        if Total_subsample >= 15:
            SNP_presence_cutoff = 0.33  # for a large group of samples
        if Total_subsample <= 3:
            SNP_presence_cutoff = 1  # for a small group of samples
            SNP_presence_sample_cutoff = 2
    else:
        cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        REF, REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
                Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
                total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
                total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                        Allels_frq_sub[Allels[allels] * 2] += int(Subdepth_forward[num_allels])
                        Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_reverse[num_allels])
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    MLF = Major_ALT[1] / total_sub_depth
                    if total_sub_depth_forward >= Sample_depth_cutoff and \
                            total_sub_depth_reverse >= Sample_depth_cutoff:
                        # forward and reverse cutoff POS detected
                        if MLF >= Major_alt_freq_cutoff:
                        # major alt frequency cutoff
                            Total_qualify += 1
                            # check for qualified SNP
                            if Major_ALT[0] != REF:
                                SNP_seq[-1] = Major_ALT[0]  # unqualified SNP also include in alignment
                                # qualified SNP
                                if (total_sub_depth_forward - int(Subdepth_forward[REF_where]) >= SNP_depth_cutoff and \
                                        total_sub_depth_reverse - int(Subdepth_reverse[REF_where]) >= SNP_depth_cutoff) and \
                                        Major_ALT[1] >= SNP_depth_cutoff:
                                    Total_qualify_SNP += 1
                                    SNP.add(genome_order)  # only take qualified SNP as valid SNP
                                    SNP_seq[-1] = Major_ALT[0]
                            else:
                                Total_qualify_notSNP += 1
                        else:
                            # major alt frequency low
                            Total_unqualify_alt_freq += 1
            sample_num += 1
        if Total_qualify / Total_subsample >= SNP_presence_cutoff and \
                Total_unqualify_alt_freq / Total_subsample <= Poor_MLF_freq_cutoff and\
                Total_qualify >= SNP_presence_sample_cutoff and \
                Total_qualify_SNP >= 1 and Total_qualify_SNP < Total_qualify and\
                Total_qualify_notSNP > 0:
            # -> qualified SNP
            # qualified samples cutoff + unqualified samples cutoff
            # at least 1 qualified SNP
            # calculate NS
            gene_info = contig_to_gene(CHR, POS)
            if gene_info!= []:
                Chr_gene, POS_gene,codon_start,Ref_seq_chr,Reverse_chr  = gene_info
                if Ref_seq_chr != 'None':
                    #  observed NS ratio calculated
                    temp_snp_line_NS= [Chr_gene,str(POS_gene),'']
                    if codon_start <= POS_gene - 1:
                        Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                        Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                        SNP_seq_chr = Ref_seq_chr
                        if len(Ref_seq_codon) == 3:
                            Ref_seq_aa = translate(Ref_seq_codon)[0]
                            temp_snp_line_AA += Ref_seq_aa
                            ALT_set = allels_set
                            ALT_set.remove(REF)
                            for ALT in ALT_set:
                                SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                SNP_seq_aa = translate(SNP_seq_codon)[0]
                                temp_snp_line_AA += SNP_seq_aa
                                temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                temp_snp_line_NS[-1]+=temp_NorS
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
            temp_snp_line.append(CHR)
            temp_snp_line.append(str(POS))
            temp_snp_line.append(REF)
            temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
            vcf_file_list.append('\t'.join(temp_snp_line)+ '\t' +'\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,'\t'.join(temp_snp_line_NS),temp_snp_line_AA))
            vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
            vcf_file_list_vcf.append('\t'.join(lines_set[0:9])+'\t'+'\t'.join(lines_set_sub)+'\n')
            i = 9
            j = 0
            SNP_alignment[reference_name] += REF
            for genomename in SNP_alignment:
                if genomename != reference_name:
                    if i in cluster_sub:
                        SNP_alignment[genomename] += SNP_seq[j]
                        j += 1
                    i += 1
    return [CHR_old,POS_old]

def SNP_check_all_fq(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        REF, REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
                Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                        Allels_frq_sub[Allels[allels] * 2] += int(Subdepth_forward[num_allels])
                        Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_reverse[num_allels])
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    if Major_ALT[0] != REF:
                        SNP_seq[-1] = Major_ALT[0]  # unqualified SNP also include in alignment
                        SNP.add(genome_order)
            sample_num += 1
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
                        ALT_set.remove(REF)
                        for ALT in ALT_set:
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
        temp_snp_line.append(CHR)
        temp_snp_line.append(str(POS))
        temp_snp_line.append(REF)
        temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
        vcf_file_list.append(
            '\t'.join(temp_snp_line) + '\t' + '\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass, '\t'.join(temp_snp_line_NS),
                temp_snp_line_AA))
        vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
        vcf_file_list_vcf.append('\t'.join(lines_set[0:9]) + '\t' + '\t'.join(lines_set_sub) + '\n')
        i = 9
        j = 0
        SNP_alignment[reference_name] += REF
        for genomename in SNP_alignment:
            if genomename != reference_name:
                if i in cluster_sub:
                    SNP_alignment[genomename] += SNP_seq[j]
                    j += 1
                i += 1
    return [CHR_old,POS_old]

def load_ref_vcf(ref_vcf_file):
    ref_chr = dict()
    for files in ref_vcf_file:
        Set_length = False
        for lines in open(files,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS, Notused, REF, ALT = lines_set[0:5]
            CHR_POS = '%s__%s'%(CHR, POS)
            ref_chr.setdefault(CHR_POS,[])
            ref_chr[CHR_POS]=[REF,ALT]
    return ref_chr

def contig_end(CHR,POS):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        total_length = CHR.split('length_')[1].split('_cov')[0]
    total_length = int(total_length)
    if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
        return True
    else:
        return False
################################################### Set up ########################################################
# set up steps
SNP_cluster = dict()
cluster_set = set()
reference_set = ['reference']
outputname_set = ['filtered']

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

end_cutoff = 70 # contig end no SNP calling
################################################### Main ########################################################
# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir_merge,'*.flt.snp.vcf'))
reference_name = reference_set[0]
output_name = outputname_set[0]
for vcf_file in all_vcf_file:
    filesize = 0
    try:
        filesize = int(os.path.getsize(vcf_file + '.%s.vcf'%(output_name)))
    except FileNotFoundError:
        pass
    if filesize == 0:
        Total = 0
        donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
        SNP_tree_cmd = []
        SNP_tree_cmd2 = []
        vcf_file_list = []
        vcf_file_list_vcf = []
        Sample_name = []
        deleting_set = []
        vcf_file_POS = []
        vcf_file_POS_candidate = set()
        SNP_alignment = dict()
        SNP_alignment.setdefault(reference_name, '')
        cov_file_list = []
        Ref_seq = dict()
        Mapping = dict()
        Mapping_loci = dict()
        CHR_old = ''
        POS_old = 0
        SNP_cluster_donor_species = dict()
        try:
            # genome mapping file
            vcf_genome = glob.glob(
                output_dir_merge + '/../merge_genome/' + '_'.join(donor_species.split('_')[:-1]) + '*.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf')[0]
            ref_chr = load_ref_vcf([vcf_genome])
            for cluster_type in cluster_set:
                SNP_cluster_donor_species.setdefault(cluster_type, [])
            for lines in open(os.path.join(input_script_sub_merge, '%s.vcf.sh' % (donor_species)), 'r'):
                if lines.startswith('bcftools mpileup '):
                    # setup samples
                    sample_set = lines.split('.fasta ')[1].split('\n')[0].split(' | ')[0].split(' ')
                    samplenum = 9
                    for samples in sample_set:
                        genomename = os.path.split(samples)[-1].split(fastq_name)[0].split('all')[0].split('.sorted.bam')[0]
                        Sample_name.append(genomename.replace('.', ''))
                        if genomename in deleting_file:
                            deleting_set.append(samplenum)
                        else:
                            SNP_alignment.setdefault(genomename, '')
                            if SNP_cluster != dict() and genomename in SNP_cluster:
                                SNP_cluster_donor_species[SNP_cluster[genomename]].append(samplenum)
                        samplenum += 1
            print('running %s' % (donor_species))
            # load database
            database_file = glob.glob(os.path.join(genome_root,
                                                   '%s/*all.spades*.fna' % (donor_species)))
            if database_file == []:
                database_file = glob.glob(os.path.join(genome_root,
                                                       '%s/*all.spades*.fna' % ('_'.join(donor_species.split('_')[:-1]))))
            if database_file != []:
                Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file[0])
            for lines in open(vcf_file, 'r'):
                if not lines.startswith("#"):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR = lines_set[0]
                    POS = int(lines_set[1])
                    # a SNP confirmed in WGS mapping
                    CHR_POS = '%s__%s' % (CHR, POS)
                    if CHR_POS in ref_chr and not contig_end(CHR, POS):
                        Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                        if Total == 0:
                            Total = len(lines_set) - 9 - len(deleting_set)
                        if "INDEL" not in lines_set[7] \
                                and (lines_set[6] != 'LowQual'):
                            CHR_old, POS_old = SNP_check_all_fq(lines_set, '',
                                                                CHR_old, POS_old, reference_name)
            outputvcf(output_name)
            outputtree(output_name)
            try:
                f1 = open(vcf_file + '.%s.cov.txt' % (output_name), 'r')
            except IOError:
                outputcov(output_name, list(vcf_file_POS_candidate), [])
        except IndexError:
            print('missing genome mapping file %s' % (
                    output_dir_merge + '/../merge_genome/vcf/' + '_'.join(
                donor_species.split('_')[:-1]) + '*.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf'))
