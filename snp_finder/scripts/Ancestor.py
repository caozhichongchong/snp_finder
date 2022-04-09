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
                      help="path of folders of all ref",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/ancester/test_data/',
                      metavar='input/')
required.add_argument("-fa",
                      help="file extension of fasta files",
                      type=str, default='_final.scaffolds.fasta',
                      metavar='_final.scaffolds.fasta')
required.add_argument("-fq",
                      help="file extension of fastq files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')
required.add_argument("-vcf",
                      help="vcf file of ancestor",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/ancester/merge/ancester.flt.snp.vcf',
                      metavar='ancester.flt.snp.vcf')
required.add_argument("-ref",
                      help="ref file for ancestor alignment",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/ancester/test_data/am_BaOv_g0001_final.scaffolds.fasta',
                      metavar='ref.fasta')
# optional output setup
optional.add_argument('-sum',
                      help="summarize the results (default: False)",
                      metavar="True or False", action='store', default='False', type=str)
optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='/scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/',
                      metavar='ancestor/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/ancester/mapper_output/',
                      metavar='mapper_output/')
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 40)",
                      metavar="1 or more", action='store', default=40, type=int)


################################################## Definition ########################################################
args = parser.parse_args()
input_script = args.s
genome_root = args.i
output_dir = args.o
genome_name = args.fa
fastq_name=args.fq
fastq_name2=args.fq.replace('1','2')
input_script_sub = '%s/ancestor/'%(input_script)
regenerate_ancestor = False
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# Set up N or S
N_S_set = dict()
N_S_set['N']=0
N_S_set['S']=1
purines=['A','G']
pyrimidines=['C','T']
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

try:
    os.mkdir(output_dir)
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

# function
def run_vcf_WGS(files,files2,database,tempbamoutput):
    # generate code
    cmds = 'bowtie2-build %s %s\n'%(database,database)
    try:
        f1 = open('%s.sorted.bam' % (tempbamoutput),'r')
    except IOError:
        cmds += 'bowtie2 --threads %s -x %s -1 %s -2 %s |samtools view -@ %s -S -b >%s.bam\nsamtools sort -@ %s %s.bam -o %s.sorted.bam\nsamtools index -@ %s %s.sorted.bam\n' % (
            min(40, args.t), database, files, files2,  min(40, args.t),
            tempbamoutput,  min(40, args.t), tempbamoutput, tempbamoutput, min(40, args.t),
            tempbamoutput)
        cmds += 'rm -r %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def merge_sample(database,vcfoutput,allsam):
    cmds = ''
    try:
        f1 = open('%s.raw.vcf' % (vcfoutput), 'r')
    except FileNotFoundError:
        cmds += 'bcftools mpileup --ff UNMAP,QCFAIL,SECONDARY --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -Ou -d3000 -A -f  %s %s | bcftools call -Ov -A -M -c --threads %s > %s.raw.vcf\n' % (
             min(40, args.t), database,
            ' '.join(allsam), min(40, args.t), vcfoutput)
    try:
        f1 = open('%s.flt.snp.vcf' % (vcfoutput))
    except FileNotFoundError:
        cmds += 'bcftools view -H -v snps -i \'QUAL>=20 && MIN(DP)>=3\' %s.raw.vcf > %s.flt.snp.vcf \n' % (
            vcfoutput, vcfoutput)
    cmds += 'rm %s\n'%(' '.join(allsam))
    if '.0.SNP.fasta.bowtie' not in vcfoutput:
        cmds += '#rm %s.raw.vcf\n' % (vcfoutput)
    return cmds

def run_minimap(files,files2,database,tempbamoutput):
    cmds = 'minimap2 -d %s.mmi %s \n' % (database, database)
    cmds += 'minimap2 -ax sr -N 10 -p 0.9 -t %s %s.mmi %s %s >%s.sam\npy39\nsamtools view -@ %s -S -b %s.sam >%s.bam\nsamtools sort -@ %s %s.bam -o %s.sorted.bam\nsamtools index -@ %s %s.sorted.bam\n' % (
        min(40, args.t), database, files, files2, tempbamoutput, min(40, args.t),tempbamoutput,
            tempbamoutput,  min(40, args.t), tempbamoutput, tempbamoutput, min(40, args.t),
            tempbamoutput)
    cmds += 'rm -r %s.sam %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def run_bwa(files,files2,database,tempbamoutput):
    cmds = 'bwa index %s\n'%(database)
    cmds += 'bwa mem %s %s %s > %s.sam\npy39\nsamtools view -@ %s -S -b %s.sam >%s.bam\nsamtools sort -@ %s %s.bam -o %s.sorted.bam\nsamtools index -@ %s %s.sorted.bam\n' % (
         database, files, files2, tempbamoutput, min(40, args.t),tempbamoutput,
            tempbamoutput,  min(40, args.t), tempbamoutput, tempbamoutput, min(40, args.t),
            tempbamoutput)
    cmds += 'rm -r %s.sam %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def run_mapper(files,files2,databaseset,tempbamoutput):
    cmds = 'time java -Xmx200g -jar %s/mapper_mat1.1.jar --max-num-matches 254 --num-threads 40 --reference %s --queries %s  --queries %s --out-refs-map-count %s.txt  --out-vcf %s.vcf\n' % (args.s,' --reference '.join(databaseset), files, files2, tempbamoutput,tempbamoutput)
    return cmds

def generate_ancestor(ref,vcffile):
    if regenerate_ancestor:
        SNPset = dict()
        vcffilename = os.path.split(vcffile)[-1].split('.flt.snp.vcf')[0]
        print(vcffilename)
        for lines in open(vcffile, 'r'):
            lines_set = lines.split('\t')
            CON, POS, nothing, REF = lines_set[:4]
            SNPset.setdefault(CON, [])
            SNPset[CON].append([int(POS) - 1, REF])
        alloutput = []
        for record in SeqIO.parse(ref, 'fasta'):
            CON = str(record.id)
            seq = str(record.seq)
            if CON in SNPset:
                allSNP = SNPset[CON]
                seq = list(seq)
                for POS, REF in allSNP:
                    if seq[POS]!=REF:
                        print('wrong POS in %s %s %s!=%s'%(CON,POS,REF,seq[POS]))
                    seq[POS] = 'N'
                seq = ''.join(seq)
            alloutput.append('>%s\n%s\n' % (CON, seq))
        f1 = open('%s/%s%s'%(genome_root,vcffilename,genome_name),'w')
        f1.write(''.join(alloutput))
        f1.close()

def load_allfasta(allgenome):
    if regenerate_ancestor:
        alloutput = []
        for genome in allgenome:
            genomename = os.path.split(genome)[-1].split(genome_name)[0]
            i = 0
            for record in SeqIO.parse(genome, 'fasta'):
                seq = str(record.seq)
                alloutput.append('>%s_CON%s\n%s\n' % (genomename, i, seq))
                i+=1
        f1 = open('%s/../allref.fasta' % (genome_root), 'w')
        f1.write(''.join(alloutput))
        f1.close()
    return '%s/../allref.fasta' % (genome_root)

def process_mapperresult(vcf,alloutput):
    vcffile = os.path.split(vcf)[-1]
    genomesum = dict()
    for lines in open(vcf,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\t')
            CHR = lines_set[0]
            genome = CHR.split('_CON')[0].replace('ancester2','ancester')
            if genome not in ['ancester','ancestor']:
                genome = 'daughter'
            depth = float(lines_set[4])
            genomesum.setdefault(genome,[0,0,0])
            genomesum[genome][0] += depth
            genomesum[genome][1] += 1
            if lines_set[3] != '' and lines_set[2] != 'N':
                genomesum[genome][2] += 1
    for genome in genomesum:
        vcfref = '%s\t%s\t%s'%(vcffile,'mapper',genome)
        alloutput.setdefault(vcfref,
                             [0, 0, 0])
        alloutput[vcfref][0] += genomesum[genome][0]
        alloutput[vcfref][1] += genomesum[genome][1]
        alloutput[vcfref][2] += genomesum[genome][2]
    return alloutput

if args.sum == 'False':
    os.system('rm -r %s' % (input_script_sub))
    try:
        os.mkdir(input_script_sub)
    except IOError:
        pass
    # generate ancestor genome
    ref = args.ref
    vcffile = args.vcf
    generate_ancestor(ref,vcffile)
    # load allfasta
    allgenome = glob.glob('%s/*%s'%(genome_root,genome_name))
    allref = load_allfasta(allgenome)
    # map fastq files
    allfastq = glob.glob('%s/*%s'%(genome_root,fastq_name))
    for fastq_file in allfastq:
        # find fastq
        fastq_file2 = fastq_file.replace(fastq_name,fastq_name2)
        fastq_filename = os.path.split(fastq_file)[-1].split(fastq_name)[0]
        # call SNPs by time bowtie2
        results = run_vcf_WGS(fastq_file, fastq_file2,
                              allref,
                              os.path.join(output_dir + '/bwa',
                                           fastq_filename + '.bowtie'))
        outputvcf = os.path.join(output_dir + '/merge',
                                 fastq_filename + '.bowtie')
        cmds = results[0]
        cmds += merge_sample(allref, outputvcf, [results[1]])
        f1 = open(os.path.join(input_script_sub, '%s.bowtie.vcf.sh' % (fastq_filename)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\npy39\n%s' % (''.join(cmds)))
        f1.close()
        # call SNPs by time bwa
        results = run_bwa(fastq_file, fastq_file2,
                          allref,
                          os.path.join(output_dir + '/bwa',
                                       fastq_filename + '.bwa'))
        outputvcf = os.path.join(output_dir + '/merge',
                                 fastq_filename + '.bwa')
        cmds = results[0]
        cmds += merge_sample(allref, outputvcf, [results[1]])
        f1 = open(os.path.join(input_script_sub, '%s.bwa.vcf.sh' % (fastq_filename)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s' % (''.join(cmds)))
        f1.close()
        # call SNPs by time minimap
        results = run_minimap(fastq_file, fastq_file2,
                              allref,
                              os.path.join(output_dir + '/bwa',
                                           fastq_filename + '.minimap'))
        outputvcf = os.path.join(output_dir + '/merge',
                                 fastq_filename + '.minimap')
        cmds = results[0]
        cmds += merge_sample(allref, outputvcf, [results[1]])
        f1 = open(os.path.join(input_script_sub, '%s.minimap.vcf.sh' % (fastq_filename)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
        f1.close()
        # call SNPs by our mapper
        # run bowtie and mapper on the same node
        cmds = run_mapper(fastq_file, fastq_file2, allgenome, os.path.join(output_dir + '/merge',
                                                                        fastq_filename + '.mapper1'))
        cmds += '#time sh %s\n' % (os.path.join(input_script_sub, '%s.bowtie.vcf.sh' % (fastq_filename)))
        cmds += '#time sh %s\n' % (os.path.join(input_script_sub, '%s.minimap.vcf.sh' % (fastq_filename)))
        cmds += '#time sh %s\n' % (os.path.join(input_script_sub, '%s.bwa.vcf.sh' % (fastq_filename)))
        f1 = open(os.path.join(input_script_sub, '%s.mapper1.vcf.sh' % (fastq_filename)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
        f1.close()
    # generate bash file for all commands
    f1 = open(os.path.join(input_script, 'allancestor.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.mapper1.vcf.sh')):
        f1.write('jobmit %s %s small\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()
else:
    alloutput = dict()
    alloutputsum = []
    allmappervcf = glob.glob('%s/merge/*.vcf'%(output_dir))
    for mapperfile in allmappervcf:
        print('start processing %s'%(mapperfile))
        alloutput = process_mapperresult(mapperfile,alloutput)
    print('start summarizing')
    for vcfref in alloutput:
        alloutputsum.append('%s\t%s\t%s\t%s\n'%(vcfref,alloutput[vcfref][0],alloutput[vcfref][1],alloutput[vcfref][2]))
    f1 = open('%s/allvcfsum.txt'%(output_dir),'w')
    f1.write('vcf\ttool\tgenome\tdepth\tcov\tSNPs\n')
    f1.write(''.join(alloutputsum))
    f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
