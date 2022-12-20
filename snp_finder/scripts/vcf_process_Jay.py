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
                      help="a folder to store all output",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/moredetails/BaFr_jay/',
                      metavar='WGS/')
required.add_argument("-fq",
                      help="file extension of WGS fastq #1 files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')
required.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='/scratch/users/anniz44/scripts/1MG/donor_species/assembly/',
                      metavar='scripts/')


################################################## Definition ########################################################
args = parser.parse_args()
# set up path
input_script = args.s  + '/vcf_filter.sh'
output_dir = args.i
vcf_name = '*.raw.vcf'
ref_filename = '.noHM.fasta'
fastq_name = args.fq
Cov_dis_overall = 1000
MW_cov_dir = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/MV_cov'

def filter_vcfs(allvcfs,lineage):
    CHRPOS_withSNPs = os.path.join(output_dir,'%s.CHRPOS'%(lineage))
    for vcf in allvcfs:
        snp_file = vcf.replace('.raw.vcf', '.flt.snp.vcf')
        os.system('cat %s | cut -f 1,2 >> %s' % (snp_file, CHRPOS_withSNPs))
    os.system('cat %s | sort | uniq > %s.unique'%(CHRPOS_withSNPs,CHRPOS_withSNPs))
    cmds = '#!/bin/bash\nsource ~/.bashrc\npy39\n'
    for vcf in allvcfs:
        cmds += 'bcftools filter -T %s.unique %s > %s.filter\n' % (CHRPOS_withSNPs, vcf, vcf)
    f1 = open(os.path.join(input_script) , 'w')
    f1.write(cmds)
    f1.close()

def outputcovwindow(allvcfs,lineage):
    cov_genome = dict()
    i = 0
    sample_len = len(allvcfs)
    for vcf in allvcfs:
        for lines in open(vcf, 'r'):
            if not lines.startswith('#'):
                lines_set = lines.split('\n')[0].split('\t')
                CHR = lines_set[0]
                POS = int(lines_set[1])
                if POS % Cov_dis_overall == 500:
                    CHRPOS = '%s\t%s' % (CHR, POS)
                    cov_genome.setdefault(CHRPOS, [0] * sample_len)
                    # subset this CHR
                    Subdepth_all = lines_set[9]
                    Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                    total_sub_depth = sum([int(Subdepth_sub) for Subdepth_sub in Subdepth])
                    cov_genome[CHRPOS][i] = total_sub_depth
        i += 1
    cov_output = []
    allCHRPOS = list(cov_genome.keys())
    allCHRPOS.sort()
    for CHRPOS in allCHRPOS:
        temp_cov = cov_genome[CHRPOS]
        cov_output.append('%s\t%s\t%s\n' % (CHRPOS, statistics.mean(temp_cov),
                                                            '\t'.join(str(cov) for cov in temp_cov)))
    vcf_file_filtered = open(MW_cov_dir + '/%s.all.donor.amnew.raw.vcf.filtered.cov.MW.txt' % (lineage), 'w')
    vcf_file_filtered.write('CHR\tPOS\tavg_depth\t%s\n' % ('\t'.join([os.path.basename(vcf).split(fastq_name)[0].split('.')[1] for vcf in allvcfs])) + ''.join(cov_output))
    vcf_file_filtered.close()

################################################### Main ########################################################
# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir,'*%s'%(vcf_name)))

# cluster donor species
all_vcf_file_set = dict()
for vcf_file in all_vcf_file:
    lineage = os.path.basename(vcf_file).split('.')[0]
    all_vcf_file_set.setdefault(lineage, [])
    all_vcf_file_set[lineage].append(vcf_file)

for lineage in all_vcf_file_set:
    allvcfs = all_vcf_file_set[lineage]
    allvcfs.sort()
    #filter_vcfs(allvcfs, lineage)
    outputcovwindow(allvcfs,lineage)