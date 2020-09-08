# start
# round 1 run genome assembly and map genomes to a reference genome
import glob
import os
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
                      action='store', default=args.sam, type=str)
optional.add_argument('-mini', '--minimap2',
                      help="Optional: complete path to minimap2 if not in PATH",
                      metavar="/usr/local/bin/minimap2",
                      action='store', default='minimap2', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
input_script_sub = args.s + '/SNP_curate'
input_script_vcf = args.s + '/vcf_round1'
input_script = args.s
genome_dir = glob.glob(os.path.join(args.i,'/*'))
output_dir = args.o + '/vcf_round1'
fastq_name = args.fq
genome_name = '.fasta'

try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa/0')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/merge')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/SNP_curate')
except IOError:
    pass

try:
    os.mkdir(input_script)
except IOError:
    pass

os.system('rm -rf %s'%(input_script_sub))
os.system('rm -rf %s'%(input_script_vcf))

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(input_script_vcf)
except IOError:
    pass


def subset(file1,file2,output_name1,output_name2):
    cmds = 'head -500000 %s >> %s\n' % (
        file1, output_name1)
    cmds += 'tail -500000 %s >> %s\n' % (
        file1, output_name1)
    cmds += 'head -500000 %s >> %s\n' % (
        file2, output_name2)
    cmds += 'tail -500000 %s >> %s\n' % (
        file2, output_name2)
    return cmds

def runspades(file1,file2,temp_output,output_name):
    cmds = '%s --careful -1 %s -2 %s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' % \
            (args.sp,file1, file2, temp_output)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_output, output_name)
    cmds += 'rm -rf %s %s %s\n' % (temp_output, file1, file2)
    return cmds

def run_vcf(genome_file,database,tempbamoutput):
    # generate code
    # for curated genome
    genome_file = genome_file + '.corrected.fasta'
    cmds = args.mini + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, args.t), database, genome_file, args.sam, min(40, args.t),
        tempbamoutput, args.sam, min(40, args.t), tempbamoutput, tempbamoutput, args.sam, min(40, args.t),
        tempbamoutput)
    cmds += '#%s depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
        args.sam,tempbamoutput, tempbamoutput)
    sample_set3 = ['%s.sorted.bam' % (tempbamoutput)]
    cmds += 'rm -rf %s.bam %s.bam.bai\n' %(tempbamoutput,tempbamoutput)
    return [cmds,sample_set3]

def run_vcf_curate(files,files2,database,tempbamoutput):
    # generate code
    cmds = ''
    Runcode = False
    try:
        f1 = open(tempbamoutput + '.flt.snp.vcf')
        filesize =int(os.path.getsize(tempbamoutput + '.flt.snp.vcf'))
        if filesize == 0:
            Runcode = True
    except IOError:
        Runcode = True
    if Runcode:
        cmds = 'rm -rf %s.fai\n' % (database)
        cmds += os.path.join(os.path.split('args.bw')[0],'bowtie2-build') + ' %s %s\n' % (database, database)
        cmds += args.bw + ' -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
            min(40, args.t), database, files, files2, args.sam, min(40, args.t),
            tempbamoutput, args.sam, min(40, args.t), tempbamoutput, tempbamoutput, args.sam, min(40, args.t),
            tempbamoutput)
        cmds += '#%s depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
            args.sam,tempbamoutput, tempbamoutput)
        cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s.sorted.bam  | %s call --threads %s -m > %s.raw.vcf\n' % (
            args.bcf, min(40, args.t), database,
            tempbamoutput, 'bcftools', min(40, args.t), tempbamoutput)
        cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
            args.bcf, min(40, args.t), tempbamoutput, tempbamoutput)
        cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
            args.bcf, tempbamoutput, tempbamoutput)
        cmds += 'rm -rf %s.bam %s.sorted.bam %s.sorted.bam.bai\n' % (
            tempbamoutput,tempbamoutput,tempbamoutput)
        cmds += 'rm -rf %s.flt.vcf %s.raw.vcf\n' % (tempbamoutput,tempbamoutput)
    else:
        os.system('rm -rf %s.fai %s.*.bt2' % (database,database))
    return cmds

def merge_sample(database,tempbamoutputGenome,samplesetGenomecorrect):
    # corrected genomes
    cmds = '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
        args.bcf, min(40, args.t), database,
        ' '.join(samplesetGenomecorrect), 'bcftools', min(40, args.t), tempbamoutputGenome)
    cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
        args.bcf, min(40, args.t), tempbamoutputGenome, tempbamoutputGenome)
    cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
        args.bcf, tempbamoutputGenome, tempbamoutputGenome)
    cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutputGenome)
    return cmds

# generate code
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    donor_species_genome = glob.glob(os.path.join(folder, '*.fasta'))
    sub_samples = []
    donor_species_fastqall = glob.glob(os.path.join(folder, 'fastq/*' + fastq_name)) + \
                             glob.glob(os.path.join(folder, '*' + fastq_name))
    if donor_species_fastqall != []:
        cmds_sub = ''
        cmds = ''
        cmds2 = ''
        i = 0
        split_bash = 1
        donor_species_folder_all = os.path.join(folder, donor_species + '_allspades1')
        donor_species_genomename = os.path.join(folder, donor_species + '.all.spades1.fasta')
        donor_species_fastq = os.path.join(folder, donor_species + '.all.spades1' + fastq_name)
        donor_species_fastq2 = os.path.join(folder, donor_species + '.all.spades1' + fastq_name.replace('1', '2'))
        tempbamoutputGenome = os.path.join(output_dir + '/merge', donor_species + '.all')
        samplesetGenomecorrect = []
        for fastq_file in donor_species_fastqall:
            if '.all' + fastq_name not in fastq_file and '.all.spades1' + fastq_name not in fastq_file:
                fastq_file_name = os.path.split(fastq_file)[-1]
                filename = fastq_file.split(fastq_name)[0]
                fastq_file2 = filename + fastq_name.replace('1', '2')
                genome_file = os.path.join(folder, fastq_file_name.split(fastq_name)[0] + genome_name)
                if glob.glob(genome_file) == []:
                    # no assembly
                    cmds_sub += runspades(fastq_file, fastq_file2, filename + '_spades', genome_file)
                    print('run spades for ', genome_file)
                # subset fastq
                cmds += subset(fastq_file, fastq_file2, donor_species_fastq, donor_species_fastq2)
                # mapping assembly to WGS to curate genome
                cmds_sub += run_vcf_curate(fastq_file, fastq_file2, genome_file,
                                           os.path.join(output_dir + '/SNP_curate',
                                                        fastq_file_name + '.curate'))
                # run each mapping corrected genome to pangenome
                results = run_vcf(genome_file, donor_species_genomename,
                                  os.path.join(output_dir + '/bwa/0',
                                               fastq_file_name)
                                  )
                cmds2 += results[0]
                samplesetGenomecorrect += results[1]
                i += 1
                if i % 10 == 0 and cmds_sub != '':
                    f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species, int(i / 10))), 'a')
                    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds_sub)))
                    f1.close()
                    cmds_sub = ''
        # pan genome
        cmds += runspades(donor_species_fastq, donor_species_fastq2, donor_species_folder_all, donor_species_genomename)
        # then run mapping corrected genome to pangenome
        cmds += 'minimap2 -d %s.mmi %s\n' % (donor_species_genomename, donor_species_genomename)
        cmds += 'rm -rf %s.fai\n' % (donor_species_genomename)
        cmds += os.path.join(os.path.split('args.bw')[0],'bowtie2-build') +' %s %s\n' % (donor_species_genomename, donor_species_genomename)
        cmds += cmds2
        cmds += merge_sample(donor_species_genomename, tempbamoutputGenome, samplesetGenomecorrect)
        if cmds_sub != '':
            f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species, int(i / 10))), 'a')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds_sub)))
            f1.close()
        f1 = open(os.path.join(input_script_vcf, '%s.vcf.sh' % (donor_species)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
        f1.close()

f1 = open(os.path.join(input_script, 'allcurate.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    if 'jobmit' in args.job:
        f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))
    else:
        f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'allround1.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.vcf.sh')):
    if 'jobmit' in args.job:
        f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))
    else:
        f1.write('nohup sh %s %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
