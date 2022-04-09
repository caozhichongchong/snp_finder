import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import random
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="input raw reads to check",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/human/SNP_model_noindel/merge/selected.reads2.fa',
                      metavar='snp.unique.fa')
required.add_argument("-o",
                      help="output folder",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/human/SNP_model_noindel/SNP_compare/',
                      metavar='SNP_compare_output')
required.add_argument("-r",
                      help="reference fasta",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/human//SNP_model_noindel/data/SRR2842672.fasta.0.SNP.fasta',
                      metavar='reference.fasta')

################################################## Definition ########################################################
args = parser.parse_args()

try:
    os.mkdir(args.o)
except IOError:
    pass

output_dir = args.o + '/SNP_compare_output'

try:
    os.mkdir(output_dir)
except IOError:
    pass

script_dir = args.o + '/SNP_compare_script'
try:
    os.mkdir(script_dir)
except IOError:
    pass

def run_vcf_WGS(files,database,tempbamoutput):
    # generate code
    cmds = '#bowtie2-build %s %s\n'%(database,database)
    try:
        f1 = open('%s.sorted.bam' % (tempbamoutput),'r')
    except IOError:
        cmds += 'bowtie2 --threads %s -x %s -U %s |samtools view -@ %s -S -b >%s.bam\nsamtools sort -@ %s %s.bam -o %s.sorted.bam\nsamtools index -@ %s %s.sorted.bam\n' % (
            40, database, files,  40,
            tempbamoutput,  40, tempbamoutput, tempbamoutput, 40,
            tempbamoutput)
        cmds += 'rm -r %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def run_mapper(files,database,tempbamoutput):
    cmds = ''
    try:
        f1 = open('%s.vcf' % (tempbamoutput), 'r')
    except FileNotFoundError:
        cmds = 'time java -Xms220g -Xmx220g -jar /scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/mapper1.31.jar --num-threads 40 --reference %s --queries %s --out-vcf %s.vcf\n' % (database, files, tempbamoutput)
    return cmds

def merge_sample(database,vcfoutput,allsam):
    cmds = ''
    try:
        f1 = open('%s.raw.vcf' % (vcfoutput), 'r')
    except FileNotFoundError:
        cmds += 'bcftools mpileup --ff UNMAP,QCFAIL,SECONDARY --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -Ou -d3000 -A -f %s %s | bcftools call -Ov -A -M -c --threads %s > %s.raw.vcf\n' % (
             40, database,
            ' '.join(allsam), 40, vcfoutput)
        cmds += 'grep -v \'#\' %s.raw.vcf > %s.flt.snp.vcf\n'%(vcfoutput,vcfoutput)
        cmds += 'rm -r %s.sorted.bam %s.sorted.bam.bai\n' % (vcfoutput, vcfoutput)
    return cmds


def generating_rawdata(seq):
    lengthseq = len(seq)
    return '@SRR2842672.1.1 1 length=%s\n'%(lengthseq)+\
    seq +'\n+SRR2842672.1.1 1 length=%s\n%s\n'%(lengthseq,'F'*(lengthseq))

# call SNPs by time bowtie2
i = 0
if True:
    for lines in open(args.i,'r'):
        if i <= 200:
            outputfilename = 'seq%s'%(i)
            seq = lines.split('\n')[0]
            # generate raw data
            fastq = os.path.join(output_dir, '%s.fq' % (outputfilename))
            f1 = open(fastq, 'w')
            f1.write(generating_rawdata(seq))
            f1.close()
            # call SNPs by bowtie
            results = run_vcf_WGS(fastq, args.r, os.path.join(output_dir,
                                               outputfilename + '.bowtie'))
            outputvcf = os.path.join(output_dir,
                                     outputfilename + '.bowtie')
            cmds = results[0]
            cmds += merge_sample(args.r, outputvcf, [results[1]])
            f1 = open(os.path.join(script_dir, '%s.bowtie.vcf.sh' % (outputfilename)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\npy39\n%s' % (''.join(cmds)))
            f1.close()
            # call SNPs by our mapper
            # run bowtie and mapper on the same node
            cmds = run_mapper(fastq, args.r, os.path.join(output_dir, outputfilename + '.mapper1'))
            cmds += 'time sh %s\n' % (os.path.join(script_dir, '%s.bowtie.vcf.sh' % (outputfilename)))
            f1 = open(os.path.join(script_dir, '%s.mapper1.vcf.sh' % (outputfilename)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            f1.close()
            i+=1

allscripts = glob.glob('%s/*mapper1.vcf.sh'%(script_dir))
f1 = open('%s/../allscripts.sh'%(script_dir),'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n#time PYTHONUNBUFFERED=1 python /scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/SNP_compare_mappertobowtiesum.py -i /scratch/users/anniz44/genomes/donor_species/SNP_curate/human/SNP_model_noindel/SNP_compare/SNP_compare_output/\n')
for script in allscripts:
    try:
        f1 = open(script + '.err','r')
    except IOError:
        f1.write('jobmit %s %s big\n'%(script,script))
f1.close()
