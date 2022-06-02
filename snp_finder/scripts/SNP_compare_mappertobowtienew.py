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
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/SNP_model/SNP_model/merge/',
                      metavar='snp.unique.fa')

################################################## Definition ########################################################
args = parser.parse_args()
output_dir_all = args.i + '/../SNP_compare/'
latestmapper = glob.glob('/scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/mapper1.*.jar')[0]
maxseq = 200
try:
    os.mkdir(output_dir_all)
except IOError:
    pass

output_dir = output_dir_all + '/SNP_compare_output'

try:
    os.mkdir(output_dir)
except IOError:
    pass

script_dir = output_dir_all + '/SNP_compare_script'
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

def run_minimap(files,database,tempbamoutput):
    cmds = '#minimap2 -d %s.mmi %s \n' % (database, database)
    cmds += 'minimap2 -ax sr -N 10 -p 0.9 -t %s %s.mmi %s >%s.sam\npy39\nsamtools view -@ %s -S -b %s.sam >%s.bam\nsamtools sort -@ %s %s.bam -o %s.sorted.bam\nsamtools index -@ %s %s.sorted.bam\n' % (
        40, database, files, tempbamoutput, 40,tempbamoutput,
            tempbamoutput,  40, tempbamoutput, tempbamoutput, 40,
            tempbamoutput)
    cmds += 'rm -r %s.sam %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def run_bwa(files,database,tempbamoutput):
    cmds = '#bwa index %s\n'%(database)
    cmds += 'bwa mem %s %s > %s.sam\npy39\nsamtools view -@ %s -S -b %s.sam >%s.bam\nsamtools sort -@ %s %s.bam -o %s.sorted.bam\nsamtools index -@ %s %s.sorted.bam\n' % (
         database, files, tempbamoutput, 40,tempbamoutput,
            tempbamoutput,  40, tempbamoutput, tempbamoutput, 40,
            tempbamoutput)
    cmds += 'rm -r %s.sam %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def run_mapper(files,database,tempbamoutput):
    cmds = ''
    try:
        f1 = open('%s.vcf' % (tempbamoutput), 'r')
    except FileNotFoundError:
        max_penalty = 0.1  # default
        if '200000' in database:
            max_penalty = 0.2
        if 'human' in args.i:
            cmds = '/usr/bin/time -v java -Xms850g -Xmx850g -jar %s  --max-penalty %s --num-threads 40 --reference %s --queries %s --out-vcf %s.vcf\n' % (
                latestmapper, max_penalty, database, files, tempbamoutput)
        else:
            cmds = '/usr/bin/time -v java -Xms10g -Xmx10g -jar %s --max-penalty %s --num-threads 40 --reference %s --queries %s --out-vcf %s.vcf\n' % (
                latestmapper, max_penalty, database, files, tempbamoutput)
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
        cmds += 'rm -r %s.raw.vcf\n'%(vcfoutput)
    return cmds


def generating_rawdata(seq):
    lengthseq = len(seq)
    return '@SRR2842672.1.1 1 length=%s\n'%(lengthseq)+\
    seq +'\n+SRR2842672.1.1 1 length=%s\n%s\n'%(lengthseq,'F'*(lengthseq))

def loadsnpseq(mapperbaselinefile):
    snpfile = mapperbaselinefile.split('.final.vcf')[0] + '.final.txt'
    ALTset = []
    seqall = set()
    for lines in open(snpfile, 'r'):
        lines_set = lines.split('\t')
        ALTset.append(lines_set[3])
    ALTset.remove('ALT')
    i = 0
    for lines in open(mapperbaselinefile, 'r'):
        lines_set = lines.split('\t')
        seqset = lines_set[-1].split('\n')[0].split(',')
        ALT = '[%s]'%(ALTset[i])
        for seq in seqset:
            if ALT in seq:
                seqall.add(seq.replace('[','').replace(']',''))
                break
        i += 1
    return list(seqall)

def compare_mapper_bowtie(seq,outputfilename,reference):
    # generate raw data
    fastq = os.path.join(output_dir, '%s.fq' % (outputfilename))
    f1 = open(fastq, 'w')
    f1.write(generating_rawdata(seq))
    f1.close()
    # call SNPs by bowtie
    results = run_vcf_WGS(fastq, reference, os.path.join(output_dir,
                                                         outputfilename + '.bowtie'))
    outputvcf = os.path.join(output_dir,
                             outputfilename + '.bowtie')
    cmds = results[0]
    cmds += merge_sample(reference, outputvcf, [results[1]])
    f1 = open(os.path.join(script_dir, '%s.bowtie.vcf.sh' % (outputfilename)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\npy39\n%s' % (''.join(cmds)))
    f1.close()
    # call SNPs by minimap
    results = run_minimap(fastq, reference, os.path.join(output_dir,
                                                         outputfilename + '.minimap'))
    outputvcf = os.path.join(output_dir,
                             outputfilename + '.minimap')
    cmds = results[0]
    cmds += merge_sample(reference, outputvcf, [results[1]])
    f1 = open(os.path.join(script_dir, '%s.minimap.vcf.sh' % (outputfilename)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\npy39\n%s' % (''.join(cmds)))
    f1.close()
    # call SNPs by bwa
    results = run_bwa(fastq, reference, os.path.join(output_dir,
                                                         outputfilename + '.bwa'))
    outputvcf = os.path.join(output_dir,
                             outputfilename + '.bwa')
    cmds = results[0]
    cmds += merge_sample(reference, outputvcf, [results[1]])
    f1 = open(os.path.join(script_dir, '%s.bwa.vcf.sh' % (outputfilename)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s' % (''.join(cmds)))
    f1.close()
    # call SNPs by our mapper
    # run mapper and others on the same node
    cmds = run_mapper(fastq, reference, os.path.join(output_dir, outputfilename + '.mapper1'))
    cmds += '/usr/bin/time -v  sh %s\n' % (os.path.join(script_dir, '%s.bowtie.vcf.sh' % (outputfilename)))
    cmds += '/usr/bin/time -v  sh %s\n' % (os.path.join(script_dir, '%s.minimap.vcf.sh' % (outputfilename)))
    if 'human' in args.i:
        cmds += '/usr/bin/time -v  #sh %s\n' % (os.path.join(script_dir, '%s.bwa.vcf.sh' % (outputfilename)))
    else:
        cmds += '/usr/bin/time -v  sh %s\n' % (os.path.join(script_dir, '%s.bwa.vcf.sh' % (outputfilename)))
    f1 = open(os.path.join(script_dir, '%s.mapper1.vcf.sh' % (outputfilename)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
    f1.close()

# find diff in baseline genome
allmapperbaseline = glob.glob('%s/*.0.SNP.fasta.mapper1.vcf.final.vcf'%(args.i))
#allmapperbaseline = glob.glob('%s/*.mapper1.vcf.final.vcf.FP'%(args.i))
for mapperbaselinefile in allmapperbaseline:
    try:
        seqall = loadsnpseq(mapperbaselinefile)
        i = 0
        genomename = os.path.split(mapperbaselinefile)[-1].split('.mapper1.vcf.final.vcf')[0]
        reference = glob.glob('%s/../data/%s'%(args.i,genomename))[0]
        for seq in seqall[:maxseq]:
            compare_mapper_bowtie(seq, genomename + '.seq%s'%(i), reference)
            i += 1
    except FileNotFoundError:
        pass

allscripts = glob.glob('%s/*mapper1.vcf.sh'%(script_dir))
f1 = open('%s/../allscripts.sh'%(script_dir),'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n#PYTHONUNBUFFERED=1 python /scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/SNP_compare_mappertobowtiesum.py -i %s\n'%(output_dir))
for script in allscripts:
    if 'human' in args.i:
        f1.write('jobmit %s %s big\n' % (script,script))
    else:
        f1.write('jobmit %s\n'%(script))
f1.close()
