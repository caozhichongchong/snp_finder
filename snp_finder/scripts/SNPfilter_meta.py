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
                      help="path of co-assembly of all populations",
                      type=str, default='.',
                      metavar='input/')
required.add_argument("-m",
                      help="path of metagenomes for tracking strains",
                      type=str, default='meta/',
                      metavar='meta/')
required.add_argument("-mfq",
                      help="file extension of metagenomes fastq #1 files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')
optional.add_argument("-d",
                      help="database for AMPHORA (protein)",
                      type=str, default='/scratch/users/anniz44/scripts/database/31_marker.fas',
                      metavar='./31_marker.fas')
# optional output setup
optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='snp_output/',
                      metavar='snp_output/')

################################################## Definition ########################################################
args = parser.parse_args()
# setup path all clonal population
output_dir = args.o + '/MG/bwa'
try:
    os.mkdir(output_dir + '/finished')
except IOError:
    pass
################################################### Set up ########################################################
# set up cutoff
Gene_cov_cutoff = 0.75 # coverage of a gene to count as mapped
Gene_num_cutoff = 0.8 # minimum percentage of genes in a donor species
SNP_max_cutoff = 0.1 # maximum ratio of SNPs on a gene to count as mapped
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

################################################### Function ########################################################
# set up functions
def vcf_to_depth(lines):
    lines_set = lines.split('\n')[0].split('\t')
    CHR = lines_set[0]
    Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
    Donorspecies_sample.setdefault(CHR,[[],[]])
    Donorspecies_sample[CHR][0].append(Depth)
    if lines_set[4] != '.':
        # a SNP
        Donorspecies_sample[CHR][1].append(Depth)

def process_donor_species(Donorspecies_sample,donorspecies):
    Donorspecies_sample2 = dict()
    for CHR in Donorspecies_sample:
        Depth_all, Depth_SNP = Donorspecies_sample[CHR]
        newCHR = '%s__%s' % (donorspecies, CHR)
        if len(Depth_all) >= Gene_cov_cutoff*Length[newCHR]:
            # enough coverage of a gene to count as mapped
            if len(Depth_SNP) < SNP_max_cutoff*Length[newCHR]:
                # less than max ratio of SNPs on a gene to count as mapped
                Donorspecies_sample2.setdefault(donorspecies, set())
                Donorspecies_sample2[donorspecies].add(CHR)
            else:
                print('not qualified for max SNPs %s %s %s'% (len(Depth_SNP),Length[newCHR],CHR))
        else:
            print('not qualified for enough coverage %s %s' % (len(Depth_all)/Length[newCHR], CHR))
    for donorspecies in Donorspecies:
        if donorspecies in Donorspecies_sample2:
            allCHR = Donorspecies_sample2[donorspecies]
            if len(allCHR) >= Gene_num_cutoff * Donorspecies[donorspecies]:
                # enough percentage of all genes in a donor species
                Depth_all = []
                Depth_all_SNP = []
                for CHR in allCHR:
                    Depth_all += Donorspecies_sample[CHR][0]
                    Depth_all_SNP += Donorspecies_sample[CHR][1]
                max_depth = statistics.mean(Depth_all) + 1.5*statistics.stdev(Depth_all)
                Depth_all_filter = [i for i in Depth_all if i < max_depth]
                Coverage.append('%s\t%s\t%.1f\t%.3f\t%.2f\t\n'%(donorspecies,samplename,statistics.mean(Depth_all_filter),
                                                            statistics.stdev(Depth_all_filter),len(Depth_all_SNP)/len(Depth_all)))
            else:
                print('not qualified for enough genes %s %s' % (len(allCHR)/Donorspecies[donorspecies],donorspecies))
                Coverage.append(
                    '%s\t%s\t%.1f\t%.3f\t%.2f\t\n' % (donorspecies, samplename, 0,0,0))

################################################### Main ########################################################
# load gene length
Length = dict()
Donorspecies = dict()
for record in SeqIO.parse(open(os.path.join(args.i,'all.31marker.fna'), 'r'), 'fasta'):
    record_id = str(record.id)
    donorspecies = record_id.split('__')[0]
    Donorspecies.setdefault(donorspecies,0)
    Donorspecies[donorspecies] += 1
    Length.setdefault(record_id, len(str(record.seq)))

# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir,'*.raw.vcf'))
Coverage = []
Coverage.append('donor_species\tsample\tavg_depth\tstdev_depth\tSNP_num\t\n')

for vcf_file in all_vcf_file:
    samplename = os.path.split(vcf_file)[-1].split(args.mfq)[0]
    donorspecies = os.path.split(vcf_file)[-1].split(args.mfq + '.')[1].split('.raw.vcf')[0]
    print('processing',samplename,donorspecies)
    Donorspecies_sample = dict()
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            vcf_to_depth(lines)
    process_donor_species(Donorspecies_sample,donorspecies)
    os.system('mv %s %s/finished/'%(vcf_file,output_dir))

f1 = open(os.path.join(output_dir + '/../', 'MG.abu.sum.txt'),'a')
f1.write(''.join(Coverage))
f1.close()

################################################### END ########################################################
