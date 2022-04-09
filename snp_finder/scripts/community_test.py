import glob
import os
input_fastq = '/scratch/users/xy43/mapper_testing/'
target_folder = ['diff_bc_fq','same_cc_fq','test_reads']
reference_dir = '/scratch/users/anniz44/genomes/donor_species/SNP_curate/xq_data/reference'
script_dir = '/scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/xq_data'
output_dir = '/scratch/users/anniz44/genomes/donor_species/SNP_curate/xq_data/'
mapper_dir = '/scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/'
try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(script_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/merge')
except IOError:
    pass

def run_mapper(files,files2,database,tempbamoutput):
    # NEED CHANGE IF TIME EVALUATION
    cmds = '#time java -Xms900g -Xmx900g -jar %s/mapper1.19.jar --max-num-matches 1 --max-penalty 0.05 --num-threads 40 --reference %s --queries %s  --queries %s --out-vcf %s.vcf\n' % (mapper_dir,database, files, files2, tempbamoutput)
    return cmds

def run_vcf_WGS(files,files2,database,tempbamoutput):
    # generate code
    cmds = '#bowtie2-build %s %s\n'%(database,database)
    try:
        f1 = open('%s.sorted.bam' % (tempbamoutput),'r')
    except IOError:# -0.5 = 8 mutations, -0.3125 = 5 mutations = 95% similarity for 100 bp
        cmds += 'bowtie2 --score-min L,0,-0.3125 -k 1 --threads %s -x %s -1 %s -2 %s |samtools view -@ %s -S -b >%s.bam\nsamtools sort -@ %s %s.bam -o %s.sorted.bam\nsamtools index -@ %s %s.sorted.bam\n' % (
            40, database, files, files2,  40,
            tempbamoutput,  40, tempbamoutput, tempbamoutput, 40,
            tempbamoutput)
        cmds += 'rm -rf %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def merge_sample(database,vcfoutput,allsam):
    cmds = ''
    try:
        f1 = open('%s.raw.vcf' % (vcfoutput), 'r')
    except FileNotFoundError:
        cmds += 'bcftools mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -Ou -B -d300000 -f %s %s | bcftools call -c -Ov -A --threads %s > %s.raw.vcf\n' % (
             40, database,
            ' '.join(allsam), 40, vcfoutput)
    try:
        f1 = open('%s.flt.snp.vcf' % (vcfoutput))
    except FileNotFoundError:
        cmds += 'bcftools view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
            vcfoutput, vcfoutput)
    return cmds


for folder in target_folder:
    reference_genome = '%s/%s.fa'%(reference_dir,folder)
    print('%s/%s/test*_1.fq'%(input_fastq,folder))
    allfastq = glob.glob('%s/%s/test*_1.fq'%(input_fastq,folder))
    for fastq_file in allfastq:
        fastq_file_name = os.path.split(fastq_file)[-1].split('_1.fq')[0]
        fastq_file2 = fastq_file.replace('_1.fq','_2.fq')
        # call SNPs by time bowtie2
        results = run_vcf_WGS(fastq_file, fastq_file2,
                              reference_genome,
                              os.path.join(output_dir + '/bwa',
                                           fastq_file_name + '.bowtie'))
        outputvcf = os.path.join(output_dir + '/merge',
                                 fastq_file_name + '.bowtie')
        cmds = results[0]
        cmds += merge_sample(reference_genome, outputvcf, [results[1]])
        f1 = open(os.path.join(script_dir, '%s.bowtie.vcf.sh' % (fastq_file_name)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\npy39\n%s' % (''.join(cmds)))
        f1.close()
        # call SNPs by our mapper
        # run bowtie and mapper on the same node
        # NEED CHANGE IF TIME EVALUATION
        cmds = 'time sh %s\n' % (os.path.join(script_dir, '%s.bowtie.vcf.sh' % (fastq_file_name)))
        cmds += run_mapper(fastq_file, fastq_file2, reference_genome, os.path.join(output_dir + '/merge',
                                                                                 fastq_file_name + '.mapper1'))
        f1 = open(os.path.join(script_dir, '%s.mapper1.vcf.sh' % (fastq_file_name)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
        f1.close()


f1 = open(os.path.join(script_dir, '../allxqdata.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for m in range(0,1):
    try:
        os.mkdir(script_dir + '/%s'%(m))
    except IOError:
        pass
    for sub_scripts in glob.glob(os.path.join(script_dir, '*.mapper1.vcf.sh')):
        sub_scripts_name = os.path.split(sub_scripts)[-1]
        os.system('cp %s %s/%s/%s'%(sub_scripts,script_dir,m,sub_scripts_name))
        f1.write('jobmit %s/%s/%s %s%s c7\n' % (script_dir,m,sub_scripts_name,sub_scripts_name,m))

f1.close()
