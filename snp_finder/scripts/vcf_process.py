################################################### END ########################################################
################################################### SET PATH ########################################################
# Filter results of WGS
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
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='WGS/',
                      metavar='WGS/')
# optional search parameters
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 1)",
                      metavar="1 or more", action='store', default=40, type=int)
optional.add_argument('-rd',
                      help="Round of SNP calling and filtering",
                      metavar="1-4", action='store', default=2, type=int)
optional.add_argument('-addqual',
                      help="add quality of SNP mapping",
                       action='store', default='True', type=str)
optional.add_argument('-job',
                      help="Optional: command to submit jobs",
                      metavar="nohup or customized",
                      action='store', default='jobmit', type=str)
# requirement for software calling
optional.add_argument('-pro',
                          help="Optional: complete path to prodigal if not in PATH",
                          metavar="/usr/local/bin/prodigal",
                          action='store', default='prodigal', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
# set up path
Round = args.rd
Cluster = True
Tree = True
Cov_dis = 20
Cov_dis_overall = 1000 # calculate coverage per 1000 bp
input_script = args.s
output_dir = args.o + '/vcf_round%s/bwa'%(Round)
output_dir_merge = args.o +'/vcf_round%s/merge'%(Round)
vcf_name = '.all*.raw.vcf'
ref_filename = '.noHM.fasta'
fastq_name = args.fq

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

CHR_length_cutoff = 2000
end_cutoff = 50 # contig end no SNP calling -> not used
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

def outputcovwindow(output_name):
    cov_output = []
    for lines in open(vcf_file, 'r'):
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
    vcf_file_filtered = open(vcf_file + '.%s.cov.MW.txt' % (output_name), 'w')
    vcf_file_filtered.write('CHR\tPOS\tavg_depth\t%s\n'%('\t'.join(Sample_name))+''.join(cov_output))
    vcf_file_filtered.close()
    return 'finished sampling coverage for %s'%(vcf_file)

def outputcov(output_name,vcf_file_POS_candidate):
    if len(vcf_file_list) > 0:
        vcf_file_POS_candidate = '\n'.join(vcf_file_POS_candidate)
        vcf_file_POS_candidate_output = ('%s' % (vcf_file_POS_candidate))
        f1 = open(os.path.join(input_script,'grep.temp.txt'),'w')
        f1.write(vcf_file_POS_candidate_output)
        f1.close()
        os.system('grep -%s -f %s %s --no-group-separator > %s'% (
            Cov_dis, os.path.join(input_script,'grep.temp.txt'),
            vcf_file,
            vcf_file + '.%s.cov.temp' % (output_name)))
        os.system('cat %s | sort | uniq > %s' % (
            vcf_file + '.%s.cov.temp' % (output_name),
            vcf_file + '.%s.uniqcov.temp' % (output_name)))
        for lines in open(vcf_file + '.%s.uniqcov.temp' % (output_name), 'r'):
            if not lines.startswith("#"):
                vcf_to_txt(lines, cov_file_list)
        os.system('rm -rf %s %s' % (vcf_file + '.%s.cov.temp' % (output_name),
                                    vcf_file + '.%s.uniqcov.temp' % (output_name)))
        vcf_file_filtered = open(vcf_file + '.%s.cov.txt' % (output_name), 'w')
        vcf_file_filtered.write(''.join(cov_file_list))
        vcf_file_filtered.close()
        return 'finished sampling coverage around SNPs for %s' % (vcf_file)

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

def SNP_check_all_fq(lines_set,CHR_old,POS_old,reference_name):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_frq2 = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    SNP = set()
    NOSNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    SNP_quality = lines_set[5]
    cluster_sub = list(range(9, len(lines_set)))
    if '.' not in lines_set[4]:
        allels_set += lines_set[4].split(',')
    Total_alleles = len(allels_set)
    genome_order = 0
    Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
    REF, REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
    sample_num = 9
    for Subdepth_all in lines_set_sub:
        genome_order += 1
        Allels_frq = [0, 0, 0, 0]
        Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        Allels_frq_sub[0]= Subdepth_all.split(':')[0]
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
        total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
        Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
        Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
        for num_allels in range(0, min(len(Subdepth),Total_alleles)):
            allels = allels_set[num_allels]
            Subdepth_alleles = int(Subdepth[num_allels])
            if allels in Allels:
                Allels_frq[Allels[allels]] += Subdepth_alleles
                Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_forward[num_allels])
                Allels_frq_sub[Allels[allels] * 2 + 2] += int(Subdepth_reverse[num_allels])
            else:
                pass
        # find major alt and calculate frequency
        Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
        temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
        temp_snp_line_frq2.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub[1:]))
        SNP_seq.append(REF)  # set as reference
        if total_sub_depth > 0:
            if Major_ALT[0] != REF:
                SNP_seq[-1] = Major_ALT[0]
                SNP.add(genome_order)
            else:
                NOSNP.add(genome_order)
        sample_num += 1
    if SNP!= set() and NOSNP != set():
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
                if i in cluster_sub:
                    SNP_alignment[genomename] += SNP_seq[j]
                    j += 1
                i += 1
    return [CHR_old,POS_old]

def short_contig_end(CHR,POS):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        try:
            total_length = CHR.split('length_')[1].split('_cov')[0]
        except IndexError:
            return False
    total_length = int(total_length)
    if total_length >= CHR_length_cutoff:
        if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
            # not short but end
            return True
        else:
            # not short not end
            return False
    else:
        # short
        return True

def mplieup_fastq(all_subvcffiles):
    try:
        f1=open(snp_file,'r')
    except FileNotFoundError:
        snp_file_old = snp_file.replace('.flt.snp.addqual.vcf', '.flt.snp.vcf')
        os.system('cat %s | cut -f 1,2 > %s.CHRPOS'%(snp_file_old,snp_file_old))
        allchrpos = dict()
        newoutput = []
        for lines in open(snp_file_old, 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\n')[0].split('\t')
                CHR = lines_set[0]
                POS = int(lines_set[1])
                if not short_contig_end(CHR, POS):
                    CHR_POS = '%s__%s' % (CHR, POS)
                    allchrpos.setdefault(CHR_POS,['0:0/0:0,0,0:0,0:0,0:0,0']*(len(all_subvcffiles)+1))
                    allchrpos[CHR_POS][0] = '\t'.join(lines_set[0:8]) + '\tQual:GT:PL:ADF:ADR:AD'
        i = 1
        for onesnp_file in all_subvcffiles:
            cmds = '#!/bin/bash\nsource ~/.bashrc\n'
            try:
                f1 = open(onesnp_file + '.filter', 'r')
            except FileNotFoundError:
                cmds += 'bcftools filter -T %s.CHRPOS %s > %s.filter\n'%(snp_file_old,onesnp_file,onesnp_file)
            f1 = open(snp_file_old + '.sh','w')
            f1.write('%s'%(cmds))
            f1.close()
            os.system('sh %s.sh'%(snp_file_old))
            print('sh %s.sh'%(snp_file_old))
            for lines in open(onesnp_file + '.filter', 'r'):
                if not lines.startswith("#"):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR = lines_set[0]
                    POS = lines_set[1]
                    CHR_POS = '%s__%s' % (CHR, POS)
                    if CHR_POS in allchrpos:
                        SNP_quality = lines_set[5]
                        allchrpos[CHR_POS][i]=('%s:%s'%(SNP_quality,lines_set[9]))
            i += 1
        for CHR_POS in allchrpos:
            newoutput.append('\t'.join(allchrpos[CHR_POS]) + '\n')
        f1 = open(snp_file,'w')
        f1.write(''.join(newoutput))
        f1.close()

################################################### Main ########################################################
# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir_merge,'%s*%s'%(args.cluster,vcf_name)))
reference_name = reference_set[0]
output_name = outputname_set[0]
for vcf_file in all_vcf_file:
    filesize = 0
    try:
        filesize = int(os.path.getsize(vcf_file + '.%s.snpfreq.txt'%(output_name)))
    except FileNotFoundError:
        pass
    if filesize == 0:
        Total = 0
        donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
        donor_species_folder = vcf_file.split('.all')[0]
        SNP_tree_cmd = []
        SNP_tree_cmd2 = []
        vcf_file_list = []
        vcf_file_list_freq = []
        vcf_file_list_vcf = []
        Sample_name = []
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
        all_subvcffiles = []
        for lines in open(vcf_file, 'r'):
            if lines.startswith('##bcftoolsCommand=mpileup '):
                # setup samples
                sample_set = lines.split(ref_filename + ' ')[1].split('\n')[0].split(' ')
                samplenum = 9
                for samples in sample_set:
                    genomename = os.path.split(samples)[-1].split(fastq_name)[0].split('.all')[0].split('.sorted.bam')[0]
                    Sample_name.append(genomename.replace('.', ''))
                    SNP_alignment.setdefault(genomename, '')
                    samplenum += 1
                    if args.addqual =='True':
                        all_subvcffiles.append((glob.glob(donor_species_folder + '.%s%s.raw.vcf'%(genomename,fastq_name))+ \
                                               glob.glob(donor_species_folder + '.%s.raw.vcf' % (genomename))
                                                )[0])
            if lines.startswith('##reference=file:'):
                database_file = lines.split('##reference=file:')[1].split('\n')[0]
                break
        print('running %s' % (donor_species))
        snp_file = vcf_file.replace('.raw.vcf', '.flt.snp.vcf')
        if args.addqual =='True':
            mplieup_fastq(all_subvcffiles)
            snp_file = vcf_file.replace('.raw.vcf', '.flt.snp.addqual.vcf')
        if database_file.split('.')[-1] != '.fna':
            # not gene file
            try:
                f1 = open(database_file + '.fna','r')
            except FileNotFoundError:
                os.system('%s -q -i %s -d %s.fna' % (args.pro, database_file, database_file))
            database_file = database_file + '.fna'
        Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
        for lines in open(snp_file, 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\n')[0].split('\t')
                CHR = lines_set[0]
                POS = int(lines_set[1])
                CHR_POS = '%s__%s' % (CHR, POS)
                Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                if Total == 0:
                    Total = len(lines_set) - 9
                if "INDEL" not in lines_set[7] \
                        and (lines_set[6] != 'LowQual'):
                    CHR_old, POS_old = SNP_check_all_fq(lines_set,
                                                        CHR_old, POS_old, reference_name)
        outputvcf(output_name)
        outputtree(output_name)
        if args.addqual =='True':
            try:
                f1 = open(vcf_file + '.%s.cov.txt' % (output_name), 'r')
            except IOError:
                pass
                #outputcov(output_name, list(vcf_file_POS_candidate))
            try:
                f1 = open(vcf_file + '.%s.cov.MW.txt' % (output_name), 'r')
            except IOError:
                outputcovwindow(output_name)
