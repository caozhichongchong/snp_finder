# start
# filtering genome SNP
# round 1-2 clonal population selection vcf filering, and round 3 SNP calling, round 4 depth filtering
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
    output_dir_merge = argso + '/vcf_round%s/merge_genome' % (Round)

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
        if '.fna.flt.snp.vcf' in vcf_file:
            cov_file = vcf_file.split('.fna.flt.snp.vcf')[0] + '.fna.raw.vcf'
        else:
            cov_file = vcf_file.split('.flt.snp.vcf')[0] + '.raw.vcf'
        os.system('grep -%s -f %s %s --no-group-separator > %s'% (
            Cov_dis, os.path.join(input_script,'grep.temp.txt'),
            cov_file,
            vcf_file + '.%s.cov.temp' % (output_name)))
        os.system('cat %s | sort | uniq > %s' % (
            vcf_file + '.%s.cov.temp' % (output_name),
            vcf_file + '.%s.uniqcov.temp' % (output_name)))
        for lines in open(vcf_file + '.%s.uniqcov.temp' % (output_name), 'r'):
            if not lines.startswith("#"):
                vcf_to_txt(lines, cov_file_list,cluster_sub)
        os.system('rm -rf %s %s %s' % (vcf_file + '.%s.cov.temp' % (output_name),
                                    vcf_file + '.%s.uniqcov.temp' % (output_name),
                                       os.path.join(input_script, 'grep.temp.txt')))
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
    if Tree:
        vcf_file_filtered = open(vcf_file + '.%s.parsi.fasta' % (output_name), 'w')
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
    if all_genes == []:
        # database is a gene database
        codon_start = POS - 1 - int((POS - 1) % 3)
        Ref_seq_chr = Ref_seq.get(CHR, 'None')
        return [CHR, POS, codon_start, Ref_seq_chr, Reverse_chr]
    else:
        # database is a contig database
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

def depthcheck(vcf_genome,vcf_fq):
    os.system('cat %s | cut -f 1,2 > %s.temp' %(vcf_genome,vcf_genome))
    os.system('grep -f %s.temp %s --no-group-separator > %s.temp.depth' % (
        vcf_genome,
        vcf_fq,
        vcf_genome))
    Length = dict()
    Depth_set = dict()
    Total = 0
    for lines in open(vcf_genome + '.temp.depth'):
        lines_set = lines.split('\n')[0].split('\t')
        if Total == 0:
            Total = len(lines_set) - 9
        CHR = lines_set[0]
        if CHR not in Length:
            try:
                total_length = CHR.split('size')[1]
            except IndexError:
                total_length = CHR.split('length_')[1].split('_cov')[0]
            total_length = int(total_length)
            Length.setdefault(CHR, total_length)
        total_length = Length[CHR]
        if total_length >= Length_cutoff:
            POS = lines_set[1]
            CHRPOS = '%s\t%s' % (CHR, POS)
            lines_set_sub = lines_set[9:]
            for i in range(0,Total):
                Subdepth_all = lines_set_sub[i]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                if total_sub_depth >= Depth_cutoff:
                    Depth_set.setdefault(CHRPOS, [])
                    Depth_set[CHRPOS].append(i)
    os.system('rm -rf %s.temp*'%(vcf_genome))
    return Depth_set

def SNP_check_all(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name,SNP_presence_cutoff,SNP_presence_sample_cutoff,no_SNP_cutoff,Depth_set,cluster_sub=[]):
    CHR = lines_set[0]
    POS = int(lines_set[1])
    CHRPOS = '%s\t%s'%(CHR,POS)
    if Depth_set == {} or CHRPOS in Depth_set:
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
        REF_where=0
        if cluster_sub!= []:
            lines_set_sub = [lines_set[i] for i in cluster_sub]
            Total_subsample = len(cluster_sub)
            if Total_subsample >= 15:
                SNP_presence_cutoff = 0.33  # for a large group of samples
            elif Total_subsample in [3,4]:
                SNP_presence_cutoff = 1  # for a small group of samples
                SNP_presence_sample_cutoff = 2
            elif Total_subsample in [1,2]:
                SNP_presence_cutoff = 1  # for a small group of samples
                SNP_presence_sample_cutoff = 1
                no_SNP_cutoff = 0
        else:
            cluster_sub = list(range(9, len(lines_set)))
        if Total_subsample > 0:
            if '.' not in lines_set[4]:
                allels_set += lines_set[4].split(',')
            Total_alleles = len(allels_set)
            genome_order = 0
            Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
            if Total_subsample > 2:
                REF,REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
            sample_num = 9
            for Subdepth_all in lines_set_sub:
                if sample_num not in deleting_set:
                    genome_order += 1
                    Allels_frq = [0, 0, 0, 0]
                    Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                    total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                    for num_allels in range(0, Total_alleles):
                        allels = allels_set[num_allels]
                        Subdepth_alleles = int(Subdepth[num_allels])
                        if allels in Allels:
                            Allels_frq[Allels[allels]] += Subdepth_alleles
                        else:
                            pass
                    # find major alt and calculate frequency
                    Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                    temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq))
                    SNP_seq.append(REF)  # set as reference
                    if total_sub_depth > 0:
                        MLF = Major_ALT[1] / total_sub_depth
                        if MLF >= Major_alt_freq_cutoff:
                            # major alt frequency cutoff
                            Total_qualify += 1
                            # check for qualified SNP
                            if Major_ALT[0] != REF:
                                if (Depth_set == {} or sample_num - 9 in Depth_set[CHRPOS]): # depth cutoff
                                    # a qualified SNP
                                    temp_snp_line_pass += 'PASS'
                                    Total_qualify_SNP += 1
                                    SNP.add(genome_order)  # only take qualified SNP as valid SNP
                                    SNP_seq[-1] = Major_ALT[0] # only qualified SNP include in alignment
                                else:
                                    print('deleting sample %s of CHRPOS %s' % (sample_num - 9, CHRPOS))
                            else:
                                Total_qualify_notSNP += 1
                        else:
                            # major alt frequency low
                            Total_unqualify_alt_freq += 1
                sample_num += 1
            if Total_qualify / Total_subsample >= SNP_presence_cutoff and \
                    Total_unqualify_alt_freq < Poor_MLF_freq_cutoff and\
                    Total_qualify >= SNP_presence_sample_cutoff and \
                    Total_qualify_SNP >= 1 and Total_qualify_SNP <= Total_qualify - no_SNP_cutoff and\
                    Total_qualify_notSNP >= no_SNP_cutoff:
                # -> qualified SNP, qualified samples cutoff + unqualified samples cutoff, at least 1 qualified SNP, calculate NS
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
                if 'PASS' in temp_snp_line_pass:
                    temp_snp_line_pass = 'PASS'
                else:
                    temp_snp_line_pass = 'NOPASS'
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
    else:
        print('deleting CHRPOS %s'%(CHRPOS))
    return [CHR_old,POS_old]

def SNP_check_output(lines_set,CHR_old,POS_old,reference_name):
    CHR = lines_set[0]
    POS = int(lines_set[1])
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
    REF_where = 0
    cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        if Total_subsample > 2:
            REF, REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    qualify_loci = 0
                    MLF = Major_ALT[1] / total_sub_depth
                    if Major_ALT[0] != REF and MLF >= Major_alt_freq_cutoff or qualify_loci == 1:
                        # check for qualified SNP
                        SNP.add(genome_order)  # only take qualified SNP as valid SNP
                        SNP_seq[-1] = Major_ALT[0]  # only qualified SNP include in alignment
            sample_num += 1
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
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,
                '\t'.join(temp_snp_line_NS), temp_snp_line_AA))
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

def dis_pos(current_CHRPOS,pre_CHRPOS):
    POS1 = int(pre_CHRPOS.split('\t')[1])
    POS2 = int(current_CHRPOS.split('\t')[1])
    return POS2 - POS1

def compare_set(Ge_pre,Ge_cu):
    Ge_preset = Ge_pre.replace('\"','').split(';')
    Ge_cuset = Ge_cu.replace('\"', '').split(';')
    if (len(Ge_preset) >= 3 and len(Ge_cuset) >= 2) or (len(Ge_preset) >= 2 and len(Ge_cuset) >= 3):
        return all(elem in Ge_cuset for elem in Ge_preset) or all(elem in Ge_preset for elem in Ge_cuset)
    else:
        return False

def N_ratio(allrec,SNP_N_set):
    N = 0
    for CHRPOS in allrec:
        if SNP_N_set[CHRPOS] == 'N':
            N += 0
    return N / len(allrec)

def cluster_rec(CHR_set,SNP_gen_set,SNP_N_set,CHRPOS_set):
    for CHR in CHR_set:
        potential_rec = dict()
        allCHRPOS = CHR_set[CHR]
        for i in range(1,len(allCHRPOS)):
            current_CHRPOS = allCHRPOS[i]
            for j in reversed(range(0, i)):
                pre_CHRPOS = allCHRPOS[j]
                Ge_pre = SNP_gen_set[pre_CHRPOS]
                Ge_cu = SNP_gen_set[current_CHRPOS]
                Dis = dis_pos(current_CHRPOS,pre_CHRPOS)
                if Dis <= Rec_length_cutoff:
                    if Ge_pre == Ge_cu or compare_set(Ge_pre,Ge_cu):
                        # cluster rec sites
                        potential_rec.setdefault(Ge_pre, set())
                        potential_rec[Ge_pre].add(current_CHRPOS)
                        potential_rec[Ge_pre].add(pre_CHRPOS)
                        potential_rec.setdefault(Ge_cu, set())
                        potential_rec[Ge_cu].add(current_CHRPOS)
                        potential_rec[Ge_cu].add(pre_CHRPOS)
                        potential_rec[Ge_cu].update(list(potential_rec[Ge_pre]))
                        potential_rec[Ge_pre].update(list(potential_rec[Ge_cu]))
                        break
                else:
                    if j == i-1:
                        # output recombination
                        delete_set = []
                        for Ge_pre in potential_rec:
                            allrec = potential_rec[Ge_pre]
                            if len(allrec) >= Rec_SNP_cutoff:
                                CHRPOS_set += allrec
                            delete_set.append(Ge_pre)
                        for Ge_pre in delete_set:
                            potential_rec.pop(Ge_pre, 'None')
                    break
        # output the last recombination
        for Ge_pre in potential_rec:
            allrec = potential_rec[Ge_pre]
            if len(allrec) >= Rec_SNP_cutoff:
                CHRPOS_set += allrec
    return CHRPOS_set

def remove_rec(SNP_file):
    CHRPOS_set = []
    CHR_set = dict()
    SNP_gen_set = dict()
    SNP_N_set = dict()
    # import SNP info
    for lines in open(SNP_file,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        CHR = lines_set[0]
        POS = lines_set[1]
        CHRPOS = '%s\t%s' % (CHR, POS)
        Ge = lines_set[-6]
        NS = lines_set[-2]
        if 'NN' in NS:
            NS = 'N'
        SNP_gen_set.setdefault(CHRPOS,Ge)
        SNP_N_set.setdefault(CHRPOS, NS)
        CHR_set.setdefault(CHR, [])
        CHR_set[CHR].append(CHRPOS)
    CHRPOS_set = cluster_rec(CHR_set,SNP_gen_set,SNP_N_set,CHRPOS_set)
    return CHRPOS_set

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

# set up output
if Tree:
    try:
        os.mkdir(output_dir_merge + '/tree')
    except IOError:
        pass
    try:
        os.mkdir(input_script_sub)
    except IOError:
        pass

# set up cutoff
reference_set = ['reference']
outputname_set = ['filtered']
SNP_presence_cutoff = 0.66  # avg presence in all samples
SNP_presence_sample_cutoff = 3  # num of samples passing the above criteria
Major_alt_freq_cutoff = 0.9 # major alt freq in a genome, do not allow multiple homolougs genes
no_SNP_cutoff = 1
Poor_MLF_freq_cutoff = 1 # no sample should have homologous genes (low major alt freq)
# set up strict cutoff
SNP_presence_cutoff2 = 0.66 # avg coverage in all samples
SNP_presence_sample_cutoff2 = 3  # num of samples passing the above criteria
Poor_MLF_freq_cutoff2 = 1 # no sample should have homologous genes (low major alt freq)
# cluster cutoff Step 3
SNP_total_cutoff_2 = 100
cluster_cutoff = 2
# Depth and recombination cutoff Round4
Depth_cutoff = 10 # covered by 10 reads
Length_cutoff = 2000 # minimum ref contig length
Rec_length_cutoff = 1000 # maximum distance between recombination sites
Rec_SNP_cutoff = 4 # minumum no. of SNPs grouped/clustered as a recombination
end_cutoff = 50 # contig end no SNP calling

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

# run vcf filtering
if Round == 4:
    vcf_name = '.all.flt.snp.vcf.filtered.vcf'
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    outputname_set = ['final']
    ref_filename = '.all.spades*.fasta'
else:
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))

reference_name = reference_set[0]
output_name = outputname_set[0]
for vcf_file in all_vcf_file:
    filesize = 0
    try:
        filesize = int(os.path.getsize(vcf_file + '.%s.vcf'%(output_name)))
    except FileNotFoundError:
        pass
    if filesize == 0:
        SNP_presence_cutoff = SNP_presence_cutoff2  # for group of samples
        SNP_presence_sample_cutoff = SNP_presence_sample_cutoff2
        no_SNP_cutoff = 1
        print(vcf_file)
        Total = 0
        # filter depth
        Depth_set = dict()
        ref_chr = dict()
        if Round == 4:
            vcf_ref_file_name = os.path.split(vcf_file)[-1].split('.all.')[0]
            vcf_file_raw = \
            glob.glob(output_dir_merge + '/../../vcf_round*/merge/' + vcf_ref_file_name + '.all.raw.vcf')[0]
            try:
                # WGS
                vcf_fq = glob.glob(output_dir_merge + '/../merge/' + vcf_ref_file_name + '*.all.fq.flt.snp.vcf')[0]
                #Depth_set = depthcheck(vcf_file, vcf_fq)
                ref_chr = load_ref_vcf([vcf_fq])
            except IndexError:
                pass
        else:
            vcf_file_raw = vcf_file.replace('.flt.snp.vcf', '.raw.vcf')
        donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0].split('.flt.snp.vcf')[0]
        database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
        if len(database) > 1:
            print(vcf_file,database)
        database = database[0]
        ref_dir, ref_name = os.path.split(database)
        ref_fna = database.replace('.fasta', '.fna')
        try:
            f1 = open(ref_fna, 'r')
        except FileNotFoundError:
            os.system('prodigal -q -i %s -d %s' % (database, ref_fna))
        Sample_name = []
        deleting_set = []
        Ref_seq = dict()
        Mapping = dict()
        Mapping_loci = dict()
        for lines in open(vcf_file_raw, 'r'):
            if lines.startswith('##bcftoolsCommand=mpileup '):
                # setup samples
                sample_set = lines.split(ref_name + ' ')[1].split('\n')[0].split(' |')[0].split(' ')
                samplenum = 9
                for samples in sample_set:
                    genomename = os.path.split(samples)[-1].split(fastq_name)[0]
                    Sample_name.append(genomename.replace('.', ''))
                    if genomename in deleting_file:
                        deleting_set.append(samplenum)
                    samplenum += 1
                break
        print('running %s' % (donor_species))
        # load database
        database_file = ref_fna
        Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
        SNP_tree_cmd = []
        SNP_tree_cmd2 = []
        vcf_file_list = []
        vcf_file_list_vcf = []
        vcf_file_POS = []
        vcf_file_POS_candidate = set()
        SNP_alignment = dict()
        SNP_alignment.setdefault(reference_name, '')
        cov_file_list = []
        CHR_old = ''
        POS_old = 0
        for genomename in Sample_name:
            SNP_alignment.setdefault(genomename, '')
        for lines in open(vcf_file, 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\n')[0].split('\t')
                CHR = lines_set[0]
                POS = int(lines_set[1])
                # a SNP confirmed in WGS mapping
                CHR_POS = '%s__%s' % (CHR, POS)
                if (Round < 4 or CHR_POS in ref_chr) and not contig_end(CHR, POS):
                    Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                    if Total == 0:
                        Total = len(lines_set) - 9 - len(deleting_set)
                        if Total >= 15:
                            SNP_presence_cutoff = 0.33  # for a large group of genomes
                        elif Total in [3, 4]:
                            SNP_presence_cutoff = 1  # for a small group of genomes
                            SNP_presence_sample_cutoff = 2
                        elif Total in [1, 2]:
                            SNP_presence_cutoff = 1  # for only 1 or 2 samples, compare to ref
                            SNP_presence_sample_cutoff = 1
                            no_SNP_cutoff = 0
                    if Depth / Total >= SNP_presence_cutoff:
                        # average depth in all samples cutoff
                        if "INDEL" not in lines_set[7] \
                                and (lines_set[6] != 'LowQual'):
                            CHR_old, POS_old = SNP_check_all(lines_set, '',
                                                             CHR_old, POS_old, reference_name,
                                                             SNP_presence_cutoff,
                                                             SNP_presence_sample_cutoff, no_SNP_cutoff, Depth_set)
        outputvcf(output_name)
        outputtree(output_name)

# remove recombination
if Round == 4:
    vcf_name = '.all.flt.snp.vcf.filtered.vcf.final.vcf'
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    print(all_vcf_file)
    outputname_set = ['removerec']
    ref_filename = '.all.spades*.fasta'
    reference_name = reference_set[0]
    output_name = outputname_set[0]
    for vcf_file in all_vcf_file:
        filesize = 0
        try:
            filesize = int(os.path.getsize(vcf_file + '.%s.vcf' % (output_name)))
        except FileNotFoundError:
            pass
        if filesize == 0:
            print(vcf_file)
            Total = 0
            # filter recombination
            vcf_ref_file_name = os.path.split(vcf_file)[-1].split('.all.')[0]
            vcf_file_raw = \
                glob.glob(output_dir_merge + '/../../vcf_round*/merge/' + vcf_ref_file_name + '.all.raw.vcf')[0]
            SNP_file = vcf_file.replace('.final.vcf', '.final.snp.txt')
            CHRPOS_set = remove_rec(SNP_file)
            donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0].split('.flt.snp.vcf')[0]
            database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
            if len(database) > 1:
                print(vcf_file, database)
            database = database[0]
            ref_dir, ref_name = os.path.split(database)
            ref_fna = database.replace('.fasta', '.fna')
            try:
                f1 = open(ref_fna, 'r')
            except FileNotFoundError:
                os.system('prodigal -q -i %s -d %s' % (database, ref_fna))
            Sample_name = []
            deleting_set = []
            Ref_seq = dict()
            Mapping = dict()
            Mapping_loci = dict()
            for lines in open(vcf_file_raw, 'r'):
                if lines.startswith('##bcftoolsCommand=mpileup '):
                    # setup samples
                    sample_set = lines.split(ref_name + ' ')[1].split('\n')[0].split(' |')[0].split(' ')
                    samplenum = 9
                    for samples in sample_set:
                        genomename = os.path.split(samples)[-1].split(fastq_name)[0]
                        Sample_name.append(genomename.replace('.', ''))
                        if genomename in deleting_file:
                            deleting_set.append(samplenum)
                        samplenum += 1
                    break
            print('running %s' % (donor_species))
            # load database
            database_file = ref_fna
            Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
            SNP_tree_cmd = []
            SNP_tree_cmd2 = []
            vcf_file_list = []
            vcf_file_list_vcf = []
            vcf_file_POS = []
            vcf_file_POS_candidate = set()
            SNP_alignment = dict()
            SNP_alignment.setdefault(reference_name, '')
            cov_file_list = []
            CHR_old = ''
            POS_old = 0
            for genomename in Sample_name:
                SNP_alignment.setdefault(genomename, '')
            for lines in open(vcf_file, 'r'):
                if not lines.startswith("#"):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR = lines_set[0]
                    POS = int(lines_set[1])
                    CHRPOS = '%s\t%s' % (CHR, POS)
                    if Total == 0:
                        Total = len(lines_set) - 9 - len(deleting_set)
                    if CHRPOS not in CHRPOS_set and not contig_end(CHR, POS):
                        CHR_old, POS_old = SNP_check_output(lines_set,
                                                            CHR_old, POS_old, reference_name)
            outputvcf(output_name)
            outputtree(output_name)

# run parsi tree
if Tree:
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    for vcf_file in all_vcf_file:
        a_parsi_file = vcf_file + '.%s.parsi.fasta'%(output_name)
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
    os.system('mv %s/*.parsi* %s/tree' % (
        output_dir_merge, output_dir_merge))

# cluster cutoff
SNP_total_cutoff_2 = 1000
cluster_cutoff = 3
# run clustering
second_strain = dict()
output_name = outputname_set[0]
if Cluster:
    # sum up SNP
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    f1 = open(os.path.join(input_script, 'SNP_round%s.sum' % (Round)), 'w')
    f1.write('donor_species\tGenome\tSNP_total\tcluster\tsubcluster\t\n')
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round%s.allpair.sum' % (Round)), 'w')
    f1.write('Genome1\tGenome2\tSNP_total\t\n')
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round%s.sum.cutoff' % (Round)), 'w')
    f1.write('donor_species\tmax_cluster_diff\tSNP_cutoff\tSNP_total_len\t\n')
    f1.close()
    for vcf_file in all_vcf_file:
        fasta = vcf_file + '.%s.fasta'%(output_name)
        SNP_pair = []
        POS_info_output = []
        tree_distance_output = []
        SNP_cutoff = []
        fasta_name = os.path.split(fasta)[-1]
        donor_species = fasta_name.split('.all.flt.snp.vcf')[0]
        POS_file = glob.glob(os.path.join(output_dir_merge, donor_species + '.*.%s.POS.txt'%(output_name)))[0]
        print(donor_species)
        tree_distance = dict()
        Seq_list = dict()
        Ref = ''
        max_cluster_diff = 0
        POS_info = []
        POS_info_CHR = []
        POS_info_CHR_LEN = dict()
        Cluster_SNP = dict()
        Cluster_SNP_set = dict()
        Cluster_SNP_set_added = set()
        # load SNP POS info
        for lines in open(POS_file, 'r'):
            CHR = lines.split('\t')[0]
            POS_info.append(int(lines.split('\t')[1]))
            POS_info_CHR.append(CHR)
            POS_info_CHR_LEN.setdefault(CHR, 0)
            POS_info_CHR_LEN[CHR] = max(POS_info_CHR_LEN[CHR], int(lines.split('\t')[1]))
        # load genome SNP fasta and calculate pair-wise SNPs
        for record in SeqIO.parse(fasta, 'fasta'):
            record_name = str(record.id)
            if 'reference' not in record_name:
                new_center = 1
                record_seq = str(record.seq)
                if Ref == '':
                    # set upt the first seq as ref
                    Ref = record_seq
                    REF_name = record_name
                    new_center = 0
                    SNP_total = 0
                    SNP_total_length = len(Ref)
                    SNP_total_cutoff = SNP_total_cutoff_2
                else:
                    for record_before in Seq_list:
                        SNP_total = SNP_seq(Seq_list[record_before], record_seq, POS_info, POS_info_CHR,
                                            POS_info_CHR_LEN,
                                            POS_info_output, record_before, record_name)
                        SNP_pair.append('%s\t%s\t%s\t\n' % (record_before, record_name, SNP_total))
                        if SNP_total <= SNP_total_cutoff:
                            Cluster_SNP.setdefault(record_before, [])
                            Cluster_SNP[record_before].append(record_name)
                            max_cluster_diff = max(max_cluster_diff, SNP_total)
                    SNP_total = SNP_seq(Ref, record_seq, POS_info, POS_info_CHR, POS_info_CHR_LEN,
                                        POS_info_output, REF_name, record_name)
                    SNP_pair.append('%s\t%s\t%s\t\n' % (REF_name, record_name, SNP_total))
                    if SNP_total <= SNP_total_cutoff:
                        Cluster_SNP.setdefault(REF_name,[])
                        Cluster_SNP[REF_name].append(record_name)
                        max_cluster_diff = max(max_cluster_diff,SNP_total)
                tree_distance.setdefault(record_name, SNP_total)
                Seq_list.setdefault(record_name, record_seq)
        cluster = 0
        # cluster genomes by SNP distance
        for record_name in Cluster_SNP:
            neighbor = Cluster_SNP.get(record_name,[])
            if neighbor != [] and record_name not in Cluster_SNP_set_added:
                cluster += 1
                Cluster_SNP_set.setdefault(cluster,set())
                Cluster_SNP_set[cluster].add(record_name)
                Cluster_SNP_set_added.add(record_name)
                find_neighbor(Cluster_SNP, neighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)
        # output single genome
        for record_name in Seq_list:
            if record_name not in Cluster_SNP_set_added:
                cluster += 1
                Cluster_SNP_set.setdefault(cluster, set())
                Cluster_SNP_set[cluster].add(record_name)
                Cluster_SNP_set_added.add(record_name)
        Sub_cluster = len(Cluster_SNP_set)
        for cluster in Cluster_SNP_set:
            tree_name_list = Cluster_SNP_set[cluster]
            tree_SNP_count = 'cluster%s' % (cluster)
            second_strain.setdefault('%s_%s'%(donor_species,tree_SNP_count),[])
            for tree_name in tree_name_list:
                second_strain['%s_%s'%(donor_species,tree_SNP_count)].append(tree_name)
                tree_distance_output.append(
                    '%s\t%s\t%s\t%s\t%s\t\n' % (
                        donor_species, tree_name, tree_distance[tree_name], tree_SNP_count, Sub_cluster))
        SNP_cutoff.append('%s\t%s\t%s\t%s\t\n' % (donor_species, max_cluster_diff, SNP_total_cutoff, SNP_total_length))
        if max_cluster_diff == 0:
            # no SNPs in a cluster
            for lines in open(fasta.replace('.all.flt.snp.vcf.filtered.fasta','.all.flt.snp.vcf.filtered.samplename.txt')):
                genomenames=lines.split('\n')[0].split('\t')
                for tree_name in genomenames:
                    tree_distance_output.append(
                        '%s\t%s\t%s\t%s\t%s\t\n' % (
                            donor_species, tree_name, 0, 'cluster1', 1))
                    for tree_name2 in genomenames:
                        if tree_name != tree_name2:
                            SNP_pair.append('%s\t%s\t%s\t\n' % (tree_name, tree_name2, 0))
        f1 = open(os.path.join(input_script, 'SNP_round%s.sum'%(Round)), 'a')
        f1.write('%s' % (''.join(tree_distance_output)))
        f1.close()
        f1 = open(os.path.join(input_script, 'SNP_round%s.allpair.sum' % (Round)), 'a')
        f1.write('%s' % (''.join(SNP_pair)))
        f1.close()
        f1 = open(os.path.join(input_script, 'SNP_round%s.sum.cutoff'%(Round)), 'a')
        f1.write('%s' % (''.join(SNP_cutoff)))
        f1.close()
        # f1 = open(os.path.join(input_script, 'SNP_round%s.POS.sum'%(Round)), 'w')
        # f1.write('%s' % (''.join(POS_info_output)))
        # f1.close()

# move
if Round < 3:
    cmd_move = ''
    fastq_name = '_1.fastq'
    fastq_name_2 = fastq_name.replace('1','2')
    for donor_species in second_strain:
        try:
            os.mkdir(genome_root + '/../round%s'%(Round + 1))
        except IOError:
            pass
        folders_dir = os.path.join(genome_root + '/../round%s'%(Round + 1), donor_species)
        try:
            os.mkdir(folders_dir)
        except IOError:
            pass
        try:
            os.mkdir(folders_dir + '/fastq')
        except IOError:
            pass
        donor_species_original = donor_species.split('_cluster')[0]
        for genome_files in second_strain[donor_species]:
            genome_files = genome_files.split('sortedbam')[0].split('fasta')[0].replace('af_Pseudoflavonifractor_sp_',
                                                                                        'af_Pseudoflavonifractor_sp._')
            cmd_move += ('mv %s/../round*/%s*/%s.* %s/\n'%(genome_root,donor_species_original,genome_files,folders_dir))
            cmd_move += ('mv %s/../round*/%s*/fastq/%s%s %s/fastq/\n' % (genome_root, donor_species_original,
                                                              genome_files, fastq_name, folders_dir))
            cmd_move += ('mv %s/../round*/%s*/fastq/%s%s %s/fastq/\n' % (genome_root, donor_species_original,
                                                              genome_files, fastq_name.replace('1','2'), folders_dir))
    f1 = open(os.path.join(input_script, 'SNP_round%s.move.sh'%(Round)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmd_move)))
    f1.close()

# compare genomes of different clusters
if Round == 1:
    genome_subsample = 2 # each cluster randomly pickup 2 genomes
    def run_vcf(genome_file,database,tempbamoutput,sumfile):
        # generate code
        # for curated genome
        cmds = 'minimap2 -d %s.mmi %s\n' % (database, database)
        cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
            min(40, args.t), database, genome_file, 'samtools', min(40, args.t),
            tempbamoutput, 'samtools', min(40, args.t), tempbamoutput, tempbamoutput, 'samtools', min(40, args.t),
            tempbamoutput)
        cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s.sorted.bam | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
            'bcftools', min(40, args.t), database,
            tempbamoutput, 'bcftools', min(40, args.t), tempbamoutput)
        cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
            'bcftools', min(40, args.t), tempbamoutput, tempbamoutput)
        cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
            'bcftools', tempbamoutput, tempbamoutput)
        cmds += 'wc -l %s.flt.snp.vcf >> %s\n'%(tempbamoutput,sumfile)
        cmds += 'rm -rf %s* %s.mmi\n' % (tempbamoutput,database)
        return cmds
    def pairgenome(allgenome, note,tag):
        allcmds = ''
        i = 0
        cluster_num = len(allgenome)
        for allpairs in itertools.combinations(range(0,cluster_num), 2):
            genomeset1 = allgenome[allpairs[0]]
            genomeset2 = allgenome[allpairs[1]]
            allgenomepairs = [[a, b] for a in genomeset1 for b in genomeset2 if a != b]
            for genomepair in allgenomepairs:
                genome1, genome2 = genomepair
                if '.all.' not in genome1 and '.all.' not in genome2:
                    tempbamoutput = os.path.join(input_script_temp,
                                                 '%s__%s' % (
                                                     os.path.split(genome1)[-1].split(fasta_name)[0],
                                                     os.path.split(genome2)[-1].split(fasta_name)[0]
                                                 ))
                    i += 1
                    allcmds += run_vcf(genome2, genome1, tempbamoutput,
                                       os.path.join(input_script_temp,
                                                    '%s.%s.sum.txt'%(tag,int(i / 100))))
                    if i % 100 == 0 and allcmds != '':
                        f1 = open(os.path.join(input_script_temp, '%s.%s.sh' % (tag,int(i / 100))), 'a')
                        f1.write('#!/bin/bash\nsource ~/.bashrc\n')
                        f1.write(''.join(allcmds))
                        f1.close()
                        allcmds = ''
        if allcmds!= '':
            print(note)
            f1 = open(os.path.join(input_script_temp, '%s.%s.sh' % (tag, int(i / 100))), 'a')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n')
            f1.write(''.join(allcmds))
            f1.close()
    def group_genome(genomelist):
        Preset = dict()
        Preset2 = dict()
        Cluster = dict()
        # input cluster
        for lines in open(os.path.join(input_script, 'SNP_round%s.all.sum'%(Round)), 'r'):
            lines_set = lines.split('\t')
            # genomename, subcluster
            Cluster.setdefault(lines_set[1],lines_set[3])
        # group genomes into subcluster
        for genome in genomelist:
            genomename = os.path.split(genome)[-1]
            donor = genomename.split('_')[0]
            species = '_'.join(genomename.split('_')[1:-1])
            cluster = Cluster.get(genomename, 'cluster1')
            donor_species = '%s_%s_%s' % (donor, species, cluster)
            if species in Species_replace:
                species = Species_replace[species]
            Preset.setdefault(species, set())
            Preset[species].add(donor_species)
            Preset2.setdefault(donor_species, set())
            Preset2[donor_species].add(genome)
        for setname in Preset:
            # onsidering genomes of diff donors
            if Preset2!= dict():
                allgenome2 = []
                num_genome = 0
                for donor_species in Preset[setname]:
                    allgenome = list(Preset2[donor_species])
                    # no need to pairing genomes in a donor species
                    # use SNP_round1.allpair.sum
                    # sub set genome_subsample genomes for a donor
                    allgenome2.append(random.choices(allgenome, k=genome_subsample))
                    num_genome += genome_subsample
                if len(Preset[setname]) > 1:
                    # at least 2 donors
                    note = ('pairing genome for %s %s clusters: %s genomes' % (
                        setname,len(Preset[setname]),num_genome))
                    print(Preset[setname])
                    print(note)
                    pairgenome(allgenome2, note, setname)
    os.system('#cat %s %s > %s'
              %(os.path.join(input_script, 'SNP_round%s.sum'%(Round)),
    os.path.join(input_script2, 'SNP_round%s.sum'%(Round)),
    os.path.join(input_script, 'SNP_round%s.all.sum'%(Round))))
    # run clustering for pairwise comparison below
    if Paircompare:
        # sum up SNP
        SNP_cutoff = []
        Fasta_SNP = glob.glob(os.path.join(output_dir_merge, '*.filtered.fasta'))
        f1 = open(os.path.join(input_script, 'SNP_round%s.%s.sum' % (Round, SNP_total_cutoff_2)), 'w')
        f1.write('donor_species\tGenome\tSNP_total\tcluster1\tsubcluster\t\n')
        f1.close()
        f1 = open(os.path.join(input_script, 'SNP_round%s.%s.allpair.sum' % (Round, SNP_total_cutoff_2)), 'w')
        f1.write('Genome1\tGenome2\tSNP_total\t\n')
        f1.close()
        for fasta in Fasta_SNP:
            SNP_pair = []
            POS_info_output = []
            tree_distance_output = []
            fasta_name = os.path.split(fasta)[-1]
            donor_species = fasta_name.split('.all.flt.snp.vcf.filtered.fasta')[0]
            POS_file = glob.glob(os.path.join(output_dir_merge, donor_species + '.all.flt.snp.vcf.filtered.POS.txt'))[0]
            print(donor_species)
            tree_distance = dict()
            Seq_list = dict()
            Ref = ''
            max_cluster_diff = 0
            POS_info = []
            POS_info_CHR = []
            POS_info_CHR_LEN = dict()
            Cluster_SNP = dict()
            Cluster_SNP_set = dict()
            Cluster_SNP_set_added = set()
            # load SNP POS info
            for lines in open(POS_file, 'r'):
                CHR = lines.split('\t')[0]
                POS_info.append(int(lines.split('\t')[1]))
                POS_info_CHR.append(CHR)
                POS_info_CHR_LEN.setdefault(CHR, 0)
                POS_info_CHR_LEN[CHR] = max(POS_info_CHR_LEN[CHR], int(lines.split('\t')[1]))
            # load genome SNP fasta and calculate pair-wise SNPs
            for record in SeqIO.parse(fasta, 'fasta'):
                record_name = str(record.id)
                if 'reference' not in record_name:
                    new_center = 1
                    record_seq = str(record.seq)
                    if Ref == '':
                        # set upt the first seq as ref
                        Ref = record_seq
                        REF_name = record_name
                        new_center = 0
                        SNP_total = 0
                        SNP_total_length = len(Ref)
                        SNP_total_cutoff = SNP_total_cutoff_2
                    else:
                        for record_before in Seq_list:
                            SNP_total = SNP_seq(Seq_list[record_before], record_seq, POS_info, POS_info_CHR,
                                                POS_info_CHR_LEN,
                                                POS_info_output, record_before, record_name)
                            SNP_pair.append('%s\t%s\t%s\t\n' % (record_before, record_name, SNP_total))
                            if SNP_total <= SNP_total_cutoff:
                                Cluster_SNP.setdefault(record_before, [])
                                Cluster_SNP[record_before].append(record_name)
                                max_cluster_diff = max(max_cluster_diff, SNP_total)
                        SNP_total = SNP_seq(Ref, record_seq, POS_info, POS_info_CHR, POS_info_CHR_LEN,
                                            POS_info_output, REF_name, record_name)
                        SNP_pair.append('%s\t%s\t%s\t\n' % (REF_name, record_name, SNP_total))
                        if SNP_total <= SNP_total_cutoff:
                            Cluster_SNP.setdefault(REF_name, [])
                            Cluster_SNP[REF_name].append(record_name)
                            max_cluster_diff = max(max_cluster_diff, SNP_total)
                    tree_distance.setdefault(record_name, SNP_total)
                    Seq_list.setdefault(record_name, record_seq)
            cluster = 0
            # cluster genomes by SNP distance
            for record_name in Cluster_SNP:
                neighbor = Cluster_SNP.get(record_name, [])
                if neighbor != [] and record_name not in Cluster_SNP_set_added:
                    cluster += 1
                    Cluster_SNP_set.setdefault(cluster, set())
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    find_neighbor(Cluster_SNP, neighbor, Cluster_SNP_set, cluster, Cluster_SNP_set_added)
            # output single genome
            for record_name in Seq_list:
                if record_name not in Cluster_SNP_set_added:
                    cluster += 1
                    Cluster_SNP_set.setdefault(cluster, set())
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
            Sub_cluster = len(Cluster_SNP_set)
            for cluster in Cluster_SNP_set:
                tree_name_list = Cluster_SNP_set[cluster]
                tree_SNP_count = 'cluster%s' % (cluster)
                for tree_name in tree_name_list:
                    tree_distance_output.append(
                        '%s\t%s\t%s\t%s\t%s\t\n' % (
                            donor_species, tree_name, tree_distance[tree_name], tree_SNP_count, Sub_cluster))
            SNP_cutoff.append(
                '%s\t%s\t%s\t%s\t\n' % (donor_species, max_cluster_diff, SNP_total_cutoff, SNP_total_length))
            if max_cluster_diff == 0:  # no SNPs in a cluster
                for lines in open(
                        fasta.replace('.all.flt.snp.vcf.filtered.fasta', '.all.flt.snp.vcf.filtered.samplename.txt')):
                    genomenames = lines.split('\n')[0].split('\t')
                    for tree_name in genomenames:
                        tree_distance_output.append(
                            '%s\t%s\t%s\t%s\t%s\t\n' % (
                                donor_species, tree_name, 0, 'cluster1', 1))
                        for tree_name2 in genomenames:
                            if tree_name != tree_name2:
                                SNP_pair.append('%s\t%s\t%s\t\n' % (tree_name, tree_name2, 0))
            f1 = open(os.path.join(input_script, 'SNP_round%s.%s.sum' % (Round, SNP_total_cutoff_2)), 'a')
            f1.write('%s' % (''.join(tree_distance_output)))
            f1.close()
            f1 = open(os.path.join(input_script, 'SNP_round%s.%s.allpair.sum' % (Round, SNP_total_cutoff_2)), 'a')
            f1.write('%s' % (''.join(SNP_pair)))
            f1.close()
            # f1 = open(os.path.join(input_script, 'SNP_round%s.%s.POS.sum'%(Round,SNP_total_cutoff_2)), 'w')
            # f1.write('%s' % (''.join(POS_info_output)))
            # f1.close()
        f1 = open(os.path.join(input_script, 'SNP_round%s.%s.sum.cutoff' % (Round, SNP_total_cutoff_2)), 'w')
        f1.write('donor_species\tmax_cluster_diff\tSNP_cutoff\tSNP_total_len\t\n%s' % (''.join(SNP_cutoff)))
        f1.close()
    if Paircompare:
        os.system('rm -rf %s'%(input_script_temp))
        try:
            os.mkdir(input_script_temp)
        except IOError:
            pass
        allgenome = glob.glob('%s/*/*%s'%(genome_root, fasta_name))+\
        glob.glob('%s/*/*%s' % (genome_root2, fasta_name))
        group_genome(allgenome)
        f1 = open(os.path.join(input_script, 'allpaircmds.sh'), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n')
        for sub_scripts in glob.glob(os.path.join(input_script_temp, '*.sh')):
            if 'jobmit' in args.job:
                f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
            else:
                f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
        f1.close()
    # after finished summarize inter donor compare
    os.system('cat %s > %s'
              %(os.path.join(input_script_temp, '*.sum.txt'),
                  os.path.join(input_script, 'SNP_round%s.allpair.all.diffcluster.sum'%(Round)))
              )
    os.system('#cat %s %s > %s'
              %(os.path.join(input_script2, 'SNP_round%s.allpair.sum'%(Round)),
                  os.path.join(input_script, 'SNP_round%s.allpair.sum'%(Round)),
                os.path.join(input_script, 'SNP_round%s.allpair.all.sum' % (Round)))
              )
    def annogenome(G1):
        G1_set = G1.split('_')
        donor = G1_set[0]
        species = '_'.join(G1_set[1:-1])
        if species in Species_replace:
            species = Species_replace[species]
        return [donor,species]
    def diff(item1, item2):
        if item1 == item2:
            return 'same'
        else:
            return 'diff'
    def splitgepair(G_pair):
        genomename = os.path.split(G_pair)[-1].split('.flt.snp.vcf')[0]
        G1,G2 = genomename.split('__')
        donor1, species1 = annogenome(G1)
        donor2, species2 = annogenome(G2)
        donor_anno = '%s_donor'%(diff(donor1, donor2))
        if diff(species1, species2)=='diff':
            # check whether 2 genomes are of the same species
            print('wrong species pair: %s %s'%(G1,G2))
        if diff(donor1, donor2)=='sane':
            # check whether 2 genomes are of diff donors
            print('wrong donor pair: %s %s'%(G1,G2))
        return '\t'.join([species1,
            G1,donor1,
            G2,donor2])
    Sumoutput = []
    Sumoutput.append('SNP_count\tspecies\tgenome1\tdonor1\tgenome2\tdonor2\n')
    for lines in open(os.path.join(input_script, 'SNP_round%s.allpair.all.diffcluster.sum'%(Round)),'r'):
        lines_set = lines.split(' ')
        SNP_count, G_pair = lines_set[0:2]
        Sumoutput.append('%s\t%s\n'%(SNP_count, splitgepair(G_pair)))
    f1 = open(os.path.join(input_script, 'SNP_round%s.allpair.all.diffcluster.sumnew.txt'%(Round)),'w')
    f1.write(''.join(Sumoutput))
    f1.close()

################################################### END ########################################################
