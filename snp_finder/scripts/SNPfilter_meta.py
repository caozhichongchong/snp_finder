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

optional.add_argument("-m",
                      help="path of metagenomes for tracking strains",
                      type=str, default='meta/',
                      metavar='meta/')
optional.add_argument("-mfq",
                      help="file extension of metagenomes fastq #1 files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')
# optional search parameters
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 1)",
                      metavar="1 or more", action='store', default=1, type=int)
optional.add_argument('-rd',
                      help="Round of SNP calling and filtering",
                      metavar="2-4", action='store', default=2, type=int)
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
                      action='store', default=args.bcf, type=str)
optional.add_argument('-sam', '--samtools',
                      help="Optional: complete path to bwa if not in PATH",
                      metavar="/usr/local/bin/samtools",
                      action='store', default=args.sam, type=str)
optional.add_argument('-mini', '--minimap2',
                      help="Optional: complete path to minimap2 if not in PATH",
                      metavar="/usr/local/bin/minimap2",
                      action='store', default='minimap2', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
# setup path all clonal population

input_script = args.s + '/MG'
output_dir = args.o + '/vcf_round4/MG'
fastq_name = args.mfq
ref_merge_vcf_dir = args.o + '/vcf_round4/merge_genome'
cluster_set = {}

################################################### Set up ########################################################

# set up cutoff
Sample_depth_cutoff = 2  # both forward and reverse reads cutoff in a sample
# SNPs merge into strain cutoff
MLF_cutoff = 0.1 # minimum MLF cutoff
MLF_variant_cutoff = 0.1 # maximum variance of MLF to be grouped as one genome seq
# strains merge into representative strains cutoff
SNP_cluster_cutoff = 10 # maximum SNP cutoff to be a cluster
top_relative = 50 # maximum nest to cluster
uniq_seq_count_cutoff = 4 # minimum presence of a seq and its downstream relatives to output
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

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

def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = int(Subdepth[0]) + int(Subdepth[1])
    Subdepth_ALT = int(Subdepth[2]) + int(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return [allels_set[1],1]
    else:
        return [allels_set[0],0]

def create_SNP_seq(total_length):
    return 'N'*total_length

def change_SNP_seq(Sample_alignment_seq,subname,ALT,SNP_pos):
    Seq = Sample_alignment_seq[subname]
    Sample_alignment_seq[subname] = Seq[:SNP_pos] + ALT + Seq[SNP_pos + 1:]

def add_SNP_seq(genomename,ALT,MLF,SNP_pos,Main = True):
    if MLF > MLF_cutoff:
        Set = False
        if Main == True:
            # major alt stores as genomename.1
            change_SNP_seq(Sample_alignment_seq,genomename + '.1', ALT, SNP_pos)
            Sample_alignment_freq[genomename + '.1'].append(MLF)
        else:
            # minor alt stores as genomename.i
            for genome_sub in Sample_alignment[genomename]:
                MLF_genome = statistics.mean(Sample_alignment_freq[genome_sub])
                if MLF <= MLF_variant_cutoff +  MLF_genome and MLF >= MLF_genome - MLF_variant_cutoff:
                    # find the right sub genome sub
                    Sample_alignment_freq[genome_sub].append(MLF)
                    change_SNP_seq(Sample_alignment_seq, genome_sub, ALT, SNP_pos)
                    Set = True
                    break
            if Set == False:
                # the right sub genome sub not found
                # create another genome sub
                i = int(Sample_alignment[genomename][-1].split('.')[-1])
                newsubname = genomename + '.%s'%(i + 1)
                Sample_alignment[genomename].append(newsubname)
                Sample_alignment_seq.setdefault(newsubname, Sample_alignment_seq[genomename + '.1'])
                change_SNP_seq(Sample_alignment_seq,newsubname, ALT, SNP_pos)
                Sample_alignment_freq.setdefault(newsubname, [MLF])

def SNP_check_all(lines_set,SNP_pos):
    REF = lines_set[3]
    allels_set = [REF]
    if '.' not in lines_set[4]:
        allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        REF, REF_where = curate_REF(allels_set, Depth4)
        for Subdepth_all in lines_set[9:]:
            genomename = Sample_name[genome_order]
            genome_order += 1
            Allels_frq = [0, 0, 0, 0]
            Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0]
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
            total_sub_depth = 0
            Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
            Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
            for num_allels in range(0, Total_alleles):
                allels = allels_set[num_allels]
                Subdepth_alleles = int(Subdepth[num_allels])
                if allels in Allels:
                    forward_alt = int(Subdepth_forward[num_allels])
                    reverse_alt = int(Subdepth_reverse[num_allels])
                    if forward_alt >= Sample_depth_cutoff and \
                            reverse_alt >= Sample_depth_cutoff:
                        # forward and reverse depth cutoff to call
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                        total_sub_depth += Subdepth_alleles
                        Allels_frq_sub[Allels[allels] * 2] += int(Subdepth_forward[num_allels])
                        Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_reverse[num_allels])
                else:
                    pass
            # find major alt and calculate frequency
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            # change to major ALT
            for genome_sub in Sample_alignment[genomename]:
                change_SNP_seq(Sample_alignment_seq, genome_sub, REF, SNP_pos)
            Sample_alignment_cov[genomename].append(total_sub_depth/MGsize[genomename])
            if total_sub_depth > 0:
                add_SNP_seq(genomename, Major_ALT[0], Major_ALT[1] / total_sub_depth, SNP_pos)
                if Minor_ALT!= []:
                    for sub_alt in Minor_ALT:
                        ALT, ALT_MLF = sub_alt
                        add_SNP_seq(genomename, ALT, ALT_MLF / total_sub_depth, SNP_pos, False)

def load_ref_vcf(ref_vcf_file):
    ref_chr = dict()
    ref_chr_order = dict()
    i = 0
    total_length = 0
    for files in ref_vcf_file:
        Set_length = False
        for lines in open(files,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS, Notused, REF, ALT = lines_set[0:5]
            CHR_POS = '%s__%s'%(CHR, POS)
            ref_chr.setdefault(CHR_POS,[])
            ref_chr[CHR_POS]=[REF,ALT]
            ref_chr_order.setdefault(CHR_POS,i)
            i += 1
        fasta_file = '.'.join(files.split('.')[:-1])+'.fasta'
        for record in SeqIO.parse(fasta_file, 'fasta'):
            record_name = str(record.id)
            if not Set_length:
                total_length += len(str(record.seq))
                Set_length = True
    return [ref_chr,ref_chr_order,total_length]

def outputtree(uniq_seq):
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in Sample_alignment_seq:
        seq_len = len(Sample_alignment_seq[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, Sample_alignment_seq[genomename]))
            SNP_alignment_output_parsi.append('%s    %s\n' % (genomename[-9:], Sample_alignment_seq[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.fasta', 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.parsi.fasta', 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in uniq_seq:
        seq_len = len(uniq_seq[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, uniq_seq[genomename]))
            SNP_alignment_output_parsi.append('%s    %s\n' % (genomename[-9:], uniq_seq[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max, seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.uniq.fasta', 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.uniq.parsi.fasta', 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()

def SNP_seq(seq1, seq2, total_length):
    SNP_total = 0
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
    return SNP_total

def find_close_relative(seq,uniq_seq,uniq_ID):
    SNP_set = dict()
    # look for closest relative from the closest timepoint
    allsample_order = sorted(uniq_seq,reverse=True)
    for samplenamesub in allsample_order:
        if samplenamesub != uniq_ID:
            meta_seq = uniq_seq[samplenamesub]
            SNP_total = SNP_seq(seq, meta_seq, total_length)
            SNP_set.setdefault(SNP_total,samplenamesub)
    SNP_set2 = sorted(SNP_set)
    print(SNP_set2,SNP_set)
    return [SNP_set2[0],SNP_set[SNP_set2[0]]]

def find_close_3_relative(close_relative_set,uniq_seq_count,pre_result,top_relative = 3):
    temp_result = close_relative_set.get(pre_result[1], [])
    closest_relative_ID = pre_result[1]
    if 'First' in closest_relative_ID:
        closest_relative_ID = find_first_uniqID(close_relative_set, closest_relative_ID)
    uniq_seq_count[closest_relative_ID] += 1
    print(temp_result,top_relative,pre_result)
    if pre_result[1] in cluster_set.get(donor_species,[]):
        # a pre set cluster
        return pre_result
    if temp_result != [] and top_relative > 0 and temp_result[0] <= SNP_cluster_cutoff:
        top_relative -= 1
        pre_result = find_close_3_relative(close_relative_set, uniq_seq_count, temp_result, top_relative)
    return pre_result

def find_first_uniqID(close_relative_set,closest_relative):
    for uniqID in close_relative_set:
        if close_relative_set[uniqID][1] == closest_relative:
            return uniqID

def outputfreq():
    timepoint = sorted(Sample_alignment_seq)
    cluster = True
    meta_seq_uniq = dict()
    uniq_seq = dict()
    # calculate uniq seq number
    for samplenamesub in timepoint:
        if Sample_alignment_freq[samplenamesub] != []:
            # not empty mapping
            meta_seq = Sample_alignment_seq[samplenamesub]
            SNP_short, closest_relative = [0, 'First']
            # uniq seq
            if meta_seq not in meta_seq_uniq:
                meta_seq_uniq.setdefault(meta_seq, [])
                uniq_seq.setdefault(samplenamesub, meta_seq)
    if len(uniq_seq) <= 10:
        cluster = False
    meta_seq_uniq = dict()
    uniq_seq = dict()
    output_freq = []
    close_relative_set = dict()
    close_relative_set_output = dict()
    uniq_seq_count = dict()
    # unique meta seq
    for samplenamesub in timepoint:
        if Sample_alignment_freq[samplenamesub] != []:
            # not empty mapping
            meta_seq = Sample_alignment_seq[samplenamesub]
            SNP_short, closest_relative = [0,'First']
            # uniq seq
            if meta_seq not in meta_seq_uniq:
                meta_seq_uniq.setdefault(meta_seq, [])
                uniq_seq.setdefault(samplenamesub,meta_seq)
                uniq_seq_count.setdefault(samplenamesub, 0)
            meta_seq_uniq[meta_seq].append(samplenamesub)
            uniq_ID = meta_seq_uniq[meta_seq][0]
            uniq_seq_count[uniq_ID] += 1
            if len(uniq_seq) == 1:
                # set up the first timepoint
                close_relative_set.setdefault(uniq_ID, [SNP_short,closest_relative])
                close_relative_set_output.setdefault(uniq_ID, [SNP_short, closest_relative])
                uniq_ID_secondstrain = uniq_ID.replace('.1', '.2')
                if uniq_ID_secondstrain in timepoint:
                    close_relative_set.setdefault(uniq_ID_secondstrain, [0, 'First2'])
                    close_relative_set_output.setdefault(uniq_ID_secondstrain, [0, 'First2'])
                uniq_ID_secondstrain = uniq_ID.replace('.1', '.3')
                if uniq_ID_secondstrain in timepoint:
                    close_relative_set.setdefault(uniq_ID_secondstrain, [0, 'First3'])
                    close_relative_set_output.setdefault(uniq_ID_secondstrain, [0, 'First3'])
            # find closest relative
            if uniq_ID not in close_relative_set:
                SNP_short, closest_relative = find_close_relative(meta_seq, uniq_seq,uniq_ID)
                close_relative_set.setdefault(uniq_ID, [SNP_short, closest_relative])
                if cluster:
                    SNP_short, closest_relative = find_close_3_relative(close_relative_set, uniq_seq_count,[SNP_short, closest_relative], top_relative)
                    closest_relative_ID = closest_relative
                    if 'First' in closest_relative:
                        closest_relative_ID = find_first_uniqID(close_relative_set,closest_relative)
                    SNP_short = SNP_seq(uniq_seq[closest_relative_ID], uniq_seq[uniq_ID], total_length)
                close_relative_set_output.setdefault(uniq_ID,[SNP_short,closest_relative])
            else:
                SNP_short, closest_relative = close_relative_set_output[uniq_ID]
            # output freq
            output_freq.append('%s\t%d\t%.3f\t%s\t%s\t\n'%(uniq_ID,
                                                   int(samplenamesub.replace('am','').split('.')[0]),
                                                   statistics.mean(Sample_alignment_freq[samplenamesub]),
                                                           closest_relative,SNP_short))
    outputtree(uniq_seq)
    vcf_file_filtered = open(vcf_file + '.freq', 'w')
    vcf_file_filtered.write('#seq_type\ttimepoint\tfreq\tclosest_relative\tSNP\n' + ''.join(output_freq))
    vcf_file_filtered.close()

################################################### Main ########################################################
# read MG size
MGsize = dict()
for lines in open(os.path.join(input_script + '/..','MGsize.txt')):
    lines_set = lines.split('\n')[0].split('\t')
    MGsize.setdefault(lines_set[0].split(fastq_name)[0],2*int(lines_set[-1]))

# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir + '/merge','*.all.flt.snp.vcf'))
Coverage = []
Coverage.append('donor_species\ttimepoint\tavg_depth\tstdev_depth\tSNP_retrieved\ttotal_refSNP\t\n')
for vcf_file in all_vcf_file:
    try:
        vcf_file_filtered = open(vcf_file + '.freq', 'r')
    except FileNotFoundError:
        donor_species = os.path.split(vcf_file)[-1].split('.all.flt.snp.vcf')[0]
        Sample_depth_cutoff = 2
        # always round 2 result
        ref_vcf_file = glob.glob(os.path.join(ref_merge_vcf_dir, '%s.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf' % (donor_species)))
        print(vcf_file,ref_vcf_file)
        Sample_name = []
        Sample_alignment = dict()
        Sample_alignment_seq = dict()
        Sample_alignment_freq = dict()
        Sample_alignment_cov = dict()
        Mapped_snp = 0
        ref_chr, ref_chr_order, total_length = load_ref_vcf(ref_vcf_file)
        for lines in open(os.path.join(input_script_merge_sub, '%s.sh' % (donor_species)), 'r'):
            if lines.startswith('bcftools mpileup '):
                # setup samples
                sample_set = lines.split('.fasta ')[1].split('\n')[0].split('  |')[0].split(' ')
                for samples in sample_set:
                    genomename = os.path.split(samples)[-1].split(fastq_name)[0]
                    Sample_name.append(genomename)
                    Sample_alignment_cov.setdefault(genomename,[])
                    Sample_alignment.setdefault(genomename,[genomename + '.1'])
                    Sample_alignment_seq.setdefault(genomename + '.1', create_SNP_seq(total_length))
                    Sample_alignment_freq.setdefault(genomename + '.1', [])
        print('running %s' % (donor_species))
        for lines in open(vcf_file, 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\n')[0].split('\t')
                CHR = lines_set[0]
                POS = int(lines_set[1])
                # a SNP confirmed in genome analysis
                CHR_POS = '%s__%s' % (CHR, POS)
                if CHR_POS in ref_chr:
                    REF_ref, ALT_ref = ref_chr[CHR_POS]
                    REF_meta, ALT_meta = lines_set[3:5]
                    # same ALTs
                    if [REF_ref, ALT_ref].sort() == [REF_meta, ALT_meta].sort():
                        # add to coverage
                        Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                        Mapped_snp += 1
                        # which ALT in SNP seq
                        SNP_pos = ref_chr_order[CHR_POS]
                        SNP_check_all(lines_set, SNP_pos)
        outputfreq()
        for genomename in Sample_alignment_cov:
            cov_genome = Sample_alignment_cov[genomename]
            cov_genome_notzero = [i for i in cov_genome if i > 0]
            if len(cov_genome) > 1:
                Coverage.append('%s\t%s\t%.1e\t%.1e\t%s\t%s\t\n'%(donor_species,
                                                              int(genomename.replace('am','').split('.')[0]),
                                                              statistics.mean(cov_genome),
                                                        stdev(cov_genome),len(cov_genome_notzero),
                                                                  total_length))
            elif len(cov_genome) == 1:
                Coverage.append('%s\t%s\t%.1e\t%.1e\t%s\t%s\t\n' % (donor_species,
                                                                    int(genomename.replace('am', '').split('.')[0]),
                                                                    cov_genome[0],
                                                                    0, len(cov_genome_notzero),
                                                                    total_length))

f1 = open(os.path.join(output_dir + '/merge', 'allcov.sum.txt'),'w')
f1.write(''.join(Coverage))
f1.close()

################################################### END ########################################################
