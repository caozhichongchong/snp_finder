# start
# after round 4 set cutoff for recombination windows to calculate NS ratio
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics

Round = 4
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')+\
glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/*')
genome_root = '/scratch/users/anniz44/genomes/donor_species/*/round*'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s'%(Round)
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/merge_genome/'%(Round)
output_dir_merge2 = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/merge_genome/'%(Round)
vcf_name = '.all.flt.snp.vcf.filtered.vcf.final.snp.txt'
ref_filename = '.all.spades*.fasta'
fasta_name = '.fasta.corrected.fasta'
fastq_name = '.sorted.bam'
windowsize_set = [810000, 729000, 656100, 590490, 531441, 478296, 430466, 387419, 348677, 313809,
                  282428, 254185, 228766, 205889, 185300, 166770, 150093, 135083, 121574, 109416,
                  98474, 88626, 79763, 71786, 64607, 58146, 52331, 47097, 42387, 38148, 34333, 30899,
                  27809, 25028, 22525, 20272, 18244, 16419, 14777, 13299, 11969, 10772, 9694, 8724,
                  7851, 7065, 6358, 5722, 5149, 4634, 4170, 3753, 3377, 3039, 2735, 2461, 2214, 1992,
                  1792, 1612, 1450, 1305, 1174, 1056, 950, 855, 769, 692, 622, 559, 503, 452, 406, 365,
                  328, 295, 265, 238, 214, 192, 172, 154, 138, 124, 111, 99, 89, 80, 72, 64, 57, 51, 45,
                  40, 36, 32, 28, 25, 22, 19, 17, 15, 13, 11, 9, 8, 7, 6, 5, 4, 3, 2, 1]

try:
    os.mkdir(output_dir_merge + '/summary')
except IOError:
    pass

def windowing(Seq_N):
    for windowsize in windowsize_set:
        N_S_set = [0, 0, 0]
        total_NS = ['']
        POS_old = 0
        for CHR in Seq_N:
            # calculate interval
            try:
                total_length = CHR.split('size')[1]
            except IndexError:
                total_length = CHR.split('length_')[1].split('_cov')[0]
            total_length = int(total_length)
            if POS_old + total_length >= windowsize:
                total_interval = int(total_length / windowsize) + 1
                total_NS += [''] * (total_interval)
            # windowing SNPs
            POS_set = Seq_N[CHR][0]
            NS_set = Seq_N[CHR][1]
            for i in range(0, len(POS_set)):
                POS = POS_set[i] + POS_old
                loci_POS = int(POS / windowsize)
                if total_NS[loci_POS] == '':
                    total_NS[loci_POS] = NS_set[i]
            POS_old += total_length
        N_S_set[0] += total_NS.count('N')
        N_S_set[1] += total_NS.count('S')
        try:
            N_S_set[2] = N_S_set[0] / N_S_set[1]
        except ZeroDivisionError:
            N_S_set[2] = 'N_only'
        Output.append('%s\t%s\t%s\t%s\t%s\t\n'%(donor_species,windowsize,
                                                N_S_set[0],N_S_set[1],N_S_set[2]))

def readSNPfile(vcf_file):
    Seq_N = dict()
    for lines in open(vcf_file,'r'):
        lines_set = lines.replace('\n', '').replace('\r', '').split('\t')
        CHR = lines_set[0]
        POS = int(lines_set[1])
        N_S = lines_set[-2]
        if 'None' not in N_S and 'S' not in N_S:
            N_S = 'N'
        Seq_N.setdefault(CHR,[[],[]])
        Seq_N[CHR][0].append(POS)
        Seq_N[CHR][1].append(N_S)
    return Seq_N

# before remove rec
Output = []
all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))+\
glob.glob(os.path.join(output_dir_merge2, '*%s' % (vcf_name)))
for vcf_file in all_vcf_file:
    print(vcf_file)
    vcf_ref_file_name = os.path.split(vcf_file)[-1].split('.all.')[0]
    donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0].split('.flt.snp.vcf')[0]
    Seq_N = readSNPfile(vcf_file)
    windowing(Seq_N)

foutput = open(output_dir_merge + '/summary/all.donor.species.NSratio.new.txt', 'w')
foutput.write('donor_species\twindowsize\tN\tS\tNS_ratio\t\n')
foutput.write(''.join(Output))
foutput.close()

# after remove rec
vcf_name = '.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.snp.txt'
Output = []
all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))+\
glob.glob(os.path.join(output_dir_merge2, '*%s' % (vcf_name)))
for vcf_file in all_vcf_file:
    print(vcf_file)
    vcf_ref_file_name = os.path.split(vcf_file)[-1].split('.all.')[0]
    donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0].split('.flt.snp.vcf')[0]
    Seq_N = readSNPfile(vcf_file)
    windowing(Seq_N)

foutput = open(output_dir_merge + '/summary/all.donor.species.NSratio.removerec.new.txt', 'w')
foutput.write('donor_species\twindowsize\tN\tS\tNS_ratio\t\n')
foutput.write(''.join(Output))
foutput.close()

################################################### END ########################################################
