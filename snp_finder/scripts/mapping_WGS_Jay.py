import glob
import os
import statistics
from Bio import SeqIO
from Bio.Seq import Seq
# set up path
import argparse
from scipy.stats import poisson

fastqs = glob.glob('/scratch/users/mit_alm/gutevo/2016_09_20_Bfragilis_TS1/Bfrag_*/*_1.fastq') +\
glob.glob('/scratch/users/mit_alm/gutevo/2016_09_20_Bfragilis_TS1/S1_*/*_1.fastq') +\
glob.glob('/scratch/users/mit_alm/gutevo/2016_09_20_Bfragilis_TS1/D2_*/*_1.fastq')

reference = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/BaFr_clustercluster7/BaFr_clustercluster7.all.spades2.fasta.noHM.fasta'
outputfolder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/moredetails/BaFr_jay/'
input_script_vcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/BaFr_jay'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
genomesize_range = [0.9,1.1]  #gemove genomes out of genome size * X-Y
assembly_genome_path = outputfolder

def run_vcf_WGS(files,files2,database,tempbamoutput):
    # generate code
    cmds = ''
    try:
        f1 = open('%s.sorted.bam' % (tempbamoutput),'r')
    except IOError:
        cmds = 'bowtie2' + ' --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
            40, database, files, files2, 'samtools', 40,
            tempbamoutput, 'samtools', 40, tempbamoutput, tempbamoutput, 'samtools', 40,
            tempbamoutput)
        cmds += 'rm -rf %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def runspades_single(file1,file2,output_name):
    temp_output = output_name.split('.fasta')[0]
    cmds = '%s --careful -1 %s -2 %s -o %s --threads %s --memory 100 --cov-cutoff 7\n' % \
            ('spades.py',file1, file2, temp_output,40)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_output, output_name)
    cmds += 'rm -rf %s\n' % (temp_output)
    return cmds

def merge_sample(database,vcfoutput,allsam):
    cmds = ''
    try:
        f1 = open('%s.raw.vcf' % (vcfoutput), 'r')
    except FileNotFoundError:
        cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
            'bcftools', 40, database,
            ' '.join(allsam), 'bcftools', 40, vcfoutput)
    try:
        f1 = open('%s.flt.snp.vcf' % (vcfoutput))
    except FileNotFoundError:
        cmds += '%s view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
            'bcftools', vcfoutput, vcfoutput)
    return cmds

def check_genome_size(fastqs):
    allgenomesize = []
    for fastq_file in fastqs:
        fastq_file_name = os.path.basename(fastq_file)
        genomename = os.path.join(assembly_genome_path, fastq_file_name.replace('_1.fastq', '.fasta'))
        try:
            filesize = int(os.path.getsize(genomename))
            allgenomesize.append(filesize)
        except FileNotFoundError:
            allgenomesize.append(0)

    # filter out genomesize that's above the size cutoff
    avg_genome_size = statistics.mean(allgenomesize)
    min_genome_size = genomesize_range[0]*avg_genome_size#poisson.ppf(genomesize_range[0],avg_genome_size)
    max_genome_size = genomesize_range[1]*avg_genome_size#poisson.ppf(genomesize_range[1],avg_genome_size)
    print(avg_genome_size*1.2,min_genome_size,max_genome_size)
    qualify_fastqs = [fastqs[i] for i in range(0, len(allgenomesize)) if
                     allgenomesize[i] <= max_genome_size and allgenomesize[i] >= min_genome_size]
    return qualify_fastqs

# check genome size
qualify_fastqs = check_genome_size(fastqs)
print('not qualified samples by genome size',[x for x in fastqs if x not in qualify_fastqs])

# run mapping
for fastq_file in qualify_fastqs:
    cmds = '#!/bin/bash\nsource ~/.bashrc\npy39\n'
    fastq_file2 = fastq_file.replace('_1.fastq','_2.fastq')
    fastq_file_name = os.path.basename(fastq_file)
    results = run_vcf_WGS(fastq_file, fastq_file2,
                          reference,
                          os.path.join(outputfolder,
                                       fastq_file_name))
    cmds += results[0]
    cmds += merge_sample(reference, os.path.join(outputfolder,
                                       'BaFr_clustercluster7.' + fastq_file_name), [results[1]])
    f1 = open(os.path.join(input_script_vcf, '%s.vcf.sh' % (fastq_file_name)), 'w')
    f1.write(cmds)
    f1.close()
    cmds = runspades_single(fastq_file, fastq_file2, '%s/%s'%(outputfolder,fastq_file_name.replace('_1.fastq','.fasta')))
    f1 = open(os.path.join(input_script_vcf, '%s.assembly.sh' % (fastq_file_name)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s' % (''.join(cmds)))
    f1.close()

# sum all codes
f1 = open(os.path.join(input_script, 'allWGS.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.vcf.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run: sh %s/allWGS.sh'%(input_script))

f1 = open(os.path.join(input_script, 'allassembly.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.assembly.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run: sh %s/allassembly.sh'%(input_script))