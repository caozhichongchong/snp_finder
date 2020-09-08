# start
# after round 4 calculate dMRCA
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics

Round = 4
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')
genome_root = '/scratch/users/anniz44/genomes/donor_species/*/round*'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s'%(Round)
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/merge_genome/vcf/'%(Round)
vcf_name = '.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf.frq.snp'

def groupSNP(CHRPOS,genotypes,alltimepoint,timegroup,timegroupSNP):
    for timepoint in alltimepoint:
        samples = timegroup[timepoint]
        allSNPs = set()
        for sample in samples:
            allSNPs.add(genotypes[samplenames.index(sample)])
        if len(allSNPs) > 1:
            # a polymorphic SNP at this timepoint
            timegroupSNP[timepoint].add(CHRPOS)
    return timegroupSNP

def sumgroupSNP(timegroupSNP,timegroup):
    Output1 = []
    Output2 = []
    for timepoint in alltimepoint:
        if timepoint == alltimepoint[0]:
            oldset = list(timegroupSNP[timepoint])
        newset = list(timegroupSNP[timepoint])
        newSNP = [i for i in newset if i not in oldset]
        Output1.append('%s\t%s\t%s\t%.3f\t\n'%(timepoint - alltimepoint[0],
                                               len(newSNP),len(timegroup[timepoint]),
                                             len(newSNP)/len(timegroup[timepoint])))
        for CHRPOS in newSNP:
            Output2.append('%s\t%s\t\n' % (timepoint, CHRPOS))
    f1 = open(vcf_file + '.dmrca', 'w')
    f1.write('timepoint_diff\tnewSNP\ttotalgenome_thistime\tavgnewSNP\t\n'+ ''.join(Output1))
    f1.close()
    f1 = open(vcf_file + '.dmrca.chr.pos', 'w')
    f1.write('timepoint\tCHR\tPOS\t\n' + ''.join(Output2))
    f1.close()

# read time tag
Time = dict()
for lines in open(os.path.join(input_script,'BN10_WGS_newname_meta_multitime.txt'),'r'):
    if not lines.startswith("oldname"):
        lines_set = lines.split('\n')[0].split('\t')
        Time.setdefault(lines_set[2],int(lines_set[-1]))

# calculate dMRCA
all_vcf_file=glob.glob(os.path.join(output_dir_merge,'am*%s'%(vcf_name)))
for vcf_file in all_vcf_file:
    try:
        vcf_file_filtered = open(vcf_file + '.dmrca2', 'r')
    except FileNotFoundError:
        donor_species = os.path.split(vcf_file)[-1].split('.all.flt.snp.vcf')[0]
        timegroup = dict()
        timegroupSNP = dict()
        for lines in open(vcf_file.replace('.vcf.frq.snp','.samplename.txt'), 'r'):
            samplenames = lines.split('\n')[0].split('\t')
            samplenum = len(samplenames)
            for samplename in samplenames:
                if samplename in Time:
                    timepoint = Time[samplename]
                    timegroup.setdefault(timepoint,[])
                    timegroupSNP.setdefault(timepoint, set())
                    timegroup[timepoint].append(samplename)
        if len(timegroup) > 1:
            # at least 2 time points
            alltimepoint = sorted(timegroup)
            for lines in open(vcf_file,'r'):
                if not lines.startswith("#"):
                    # find polymorphic SNPs at each timepoint
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR = lines_set[0]
                    POS = int(lines_set[1])
                    CHRPOS = '%s\t%s'%(CHR, POS)
                    genotypes = lines_set[7:]
                    timegroupSNP = groupSNP(CHRPOS, genotypes, alltimepoint, timegroup, timegroupSNP)
            sumgroupSNP(timegroupSNP,timegroup)

all_dmrca=glob.glob(os.path.join(output_dir_merge,'am*.dmrca'))
alloutput = []
for dmrca in all_dmrca:
    donor_species = os.path.split(dmrca)[-1].split('.all.flt.snp.vcf')[0]
    species = '_'.join(donor_species.split('_')[1:-1])
    for lines in open(dmrca,'r'):
        alloutput.append('%s\t%s\t%s'%(donor_species,species,lines))

f1 = open(os.path.join(output_dir_merge, 'alldmrca.txt'),'w')
f1.write(''.join(alloutput))
f1.close()

