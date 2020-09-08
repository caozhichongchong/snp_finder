# start
# run allcurate
# step 2 check SNPs
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
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

input_script = args.s
genome_root = args.i
genome_dir = glob.glob(args.i + '/*')
output_dir = args.o + '/vcf_round1/SNP_curate'
fastq_name = '.curate.flt.snp.vcf'

# set up cutoff
# good mapping
Major_alt_freq_cutoff = 0.7 # major alt freq in a sample
Sample_depth_cutoff = 5 # both forward and reverse reads cutoff in a sample

# good coverage
total_coverage_cutoff = 0.8 # at least X reads map to its original genome
genome_avg_coverage_cutoff = 10 # genome average coverage cutoff
# reasonable curation
Major_alt_freq_cutoff2 = 0.8 # major alt freq in a sample
Sample_depth_cutoff2 = Sample_depth_cutoff # both forward and reverse reads cutoff in a sample

Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
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

def coverage_check(sh_file,Coverage):
    Coverage1 = []
    Coverage2 = []
    os.system('rm -rf %s/temp.sh %s/temp.err'%(input_script_split_sub,input_script_split_sub))
    os.system('grep \"bowtie2\" %s > %s/temp.sh' % (sh_file, input_script_split_sub))
    os.system('grep \"overall alignment rate\" %s.err > %s/temp.err' % (sh_file, input_script_split_sub))
    # load fastq name
    for lines in open('%s/temp.sh' % (input_script_split_sub), 'r'):
        try:
            filename = lines.split('>')[1].split('\n')[0].split(fastq_name + '.bam')[0]
            filename = os.path.split(filename)[-1]
            Coverage1.append(filename)
        except IndexError:
            pass
    # load reads alignment rate
    for lines in open('%s/temp.err' % (input_script_split_sub), 'r'):
        Coverage2.append(float(lines.split('%')[0]))
    i = 0
    for filename in Coverage1:
        # load average coverage of ref genome
        tempbamoutput = os.path.join(output_dir, filename + fastq_name)
        coverage_file = glob.glob(tempbamoutput + '.sorted.bam.avgcov')
        if coverage_file == []:
            cmds = '%s depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
                args.sam, tempbamoutput, tempbamoutput)
            os.system(cmds)
            coverage_file = tempbamoutput + '.sorted.bam.avgcov'
        else:
            coverage_file = coverage_file[0]
        length_total = 0
        coverage_total = 0
        donor_species = filename.split(genome_split)[0]
        Donor_species.setdefault(donor_species, [[], 0])
        try:
            for lines in open(coverage_file, 'r'):
                length_contig = float(lines.split('\t')[1])
                length_total += length_contig
                coverage_contig = float(lines.split('\t')[2])
                coverage_total += coverage_contig * length_contig
            coverage_num = Coverage2[i] / 100
            avg_coverage = coverage_total / length_total
            i += 1
        except IOError:
            coverage_num = total_coverage_cutoff - 0.1
            avg_coverage = genome_avg_coverage_cutoff - 1
        temp_line = ('%s\t%.2f\t%.1f\t%.1f' % (filename, coverage_num, length_total, avg_coverage))
        if coverage_num < total_coverage_cutoff or avg_coverage < genome_avg_coverage_cutoff:
            # not qualified
            print(filename, coverage_num, avg_coverage)
            Donor_species[donor_species][0].append(filename)
            temp_line += ('\tnot_qualified\t\n')
        else:
            # qualified
            Donor_species[donor_species][1] += 1
            temp_line += ('\tqualified\t\n')
        Coverage.append(temp_line)

def SNP_check(lines,donor_species,vcf_file_list):
    # CHR, POS, REF, ALT, good assembly, qualified mapping
    lines_set = lines.split('\n')[0].split('\t')
    report_line = ''
    temp_report_line = ['T','T']
    need_curation = 'F'
    REF = lines_set[3]
    allels_set = [REF]
    if '.' not in lines_set[4]:
        allels_set += lines_set[4].split(',')
    Total_alleles = len(allels_set)
    for Subdepth_all in lines_set[9:]:
            Allels_frq = [0, 0, 0, 0]
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
                else:
                    pass
            # find major alt and calculate frequency
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            MLF = Major_ALT[1] / total_sub_depth
            if REF != Major_ALT[0]:
                # wrong assembly
                temp_report_line[0] = 'F' # F: not good assembly
                if total_sub_depth_forward >= Sample_depth_cutoff2 and \
                        total_sub_depth_reverse >= Sample_depth_cutoff2 and \
                        MLF >= Major_alt_freq_cutoff2:
                    # can be curated
                    need_curation = 'T'  # T: need curation
            if total_sub_depth_forward < Sample_depth_cutoff or \
                    total_sub_depth_reverse < Sample_depth_cutoff or \
                    MLF < Major_alt_freq_cutoff:
                # unqualified mapping
                temp_report_line[1] = 'F' # F: bad mapping
    if temp_report_line != ['T','T']:
        report_line += '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t\n' %(donor_species,lines_set[0],lines_set[1],
                                                    REF,','.join([ALT for ALT in allels_set if ALT != REF]),
                                                    temp_report_line[0],temp_report_line[1],need_curation,
                                                     Major_ALT[0], MLF,
                                                     total_sub_depth_forward, total_sub_depth_reverse)
        vcf_file_list.append(donor_species + '\t' + lines)
    return report_line

# check major alt
all_vcf_file=glob.glob(os.path.join(output_dir,'*%s'%(fastq_name)))
vcf_file_report = []
vcf_file_list = []
vcf_file_report.append('donor_species\tCHR\tPOS\tREF\tALT\tAssembly\tMapping_quality\tNeed_curation\tMajor_ALT\tMajor_ALT_frq\tDepth_F\tDepth_R\n')
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split(fastq_name)[0]
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            vcf_file_report.append(SNP_check(lines,donor_species,vcf_file_list))

f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.sum'), 'w')
f1.write(''.join(vcf_file_report))
f1.close()

f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.vcf'), 'w')
f1.write(''.join(vcf_file_list))
f1.close()

################################################### END ########################################################
