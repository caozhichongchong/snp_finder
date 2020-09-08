# start
# metagenomes mapping to a reference genome
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

optional.add_argument("-m",
                      help="path of metagenomes for tracking strains",
                      type=str, default='meta/',
                      metavar='meta/')
optional.add_argument("-mfq",
                      help="file extension of metagenomes fastq #1 files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')
# optional search parameters
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 1)",
                      metavar="1 or more", action='store', default=1, type=int)
optional.add_argument('-rd',
                      help="Round of SNP calling and filtering",
                      metavar="2-4", action='store', default=2, type=int)
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
                      action='store', default=args.bcf, type=str)
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

input_script_sub = args.s + '/MG/vcf'
input_script_merge = args.s + '/MG/vcf_merge'
input_script = args.s + '/MG'
genome_dir = glob.glob(args.i + '/round3/*')
metagenome_dir = args.m
output_dir = args.o + '/vcf_round4/MG'
fastq_name = args.mfq

try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir + '/merge')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/bwa/0')
except IOError:
    pass

try:
    os.mkdir(input_script)
except IOError:
    pass

os.system('rm -rf %s' %(input_script_sub))

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(input_script_merge)
except IOError:
    pass

allfastq = glob.glob(os.path.join(metagenome_dir,'*%s'%(fastq_name)))
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    donor_species_genome = glob.glob(os.path.join(folder, '*.all.spades*.fasta'))
    sub_samples = []
    if donor_species_genome != []:
        cmds = ''
        i = 0
        database = donor_species_genome[0]
        try:
            f1 = open(database + '.mmi')
        except IOError:
            os.system('minimap2 -d %s.mmi %s \n' % (database, database))
        for files in allfastq:
            donor_species_dir_file = os.path.split(files)[-1]
            if donor_species_dir_file.startswith('am'):
                print(donor_species_dir_file)
                tempbamoutput = os.path.join(output_dir + '/bwa/0', donor_species_dir_file + '.' + donor_species)
                try:
                    f1 = open(tempbamoutput + '.sorted.bam')
                except IOError:
                    files2 = files.replace(fastq_name, fastq_name.replace('1', '2'))
                    if 'fasta' in fastq_name:
                        cmds += args.mini + ' -ax sr -N 1000 -p 0.8 -t %s %s.mmi %s %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                            min(40, args.t), database, files, files2, args.sam, min(40, args.t),
                            tempbamoutput, args.sam, min(40, args.t), tempbamoutput, tempbamoutput, args.sam,
                            min(40, args.t),
                            tempbamoutput)
                    else:
                        cmds += args.mini + ' -ax sr -N 1000 -p 0.8 -t %s %s.mmi %s %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                            min(40, args.t), database, files, files2, args.sam, min(40, args.t),
                            tempbamoutput, args.sam, min(40, args.t), tempbamoutput, tempbamoutput, args.sam,
                            min(40, args.t),
                            tempbamoutput)
                    cmds += 'rm -rf %s.bam\n' % (tempbamoutput)
                    sub_samples.append(tempbamoutput + '.sorted.bam')
                    cmds += 'rm -rf %s.*.unaligned\n' % (tempbamoutput)
                    i += 1
                    if i % 25 == 0:
                        f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species, int(i / 25))), 'a')
                        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                        f1.close()
                        cmds = ''
        f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species, int(i / 25))), 'a')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
        f1.close()
        # mileup
        cmds = ''
        tempbamoutput = os.path.join(output_dir + '/merge', donor_species + '.all')
        try:
            f1 = open(tempbamoutput + '.sorted.bam')
        except IOError:
            if sub_samples != []:
                cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s  | %s call  --threads %s -m > %s.raw.vcf\n' % (
                    args.bcf, min(40, args.t), database,
                    ' '.join(sub_samples), args.bcf, min(40, args.t), tempbamoutput)
                cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
                    args.bcf, min(40, args.t), tempbamoutput, tempbamoutput)
                cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
                    args.bcf, tempbamoutput, tempbamoutput)
                cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
        f1 = open(os.path.join(input_script_merge, '%s.sh' % donor_species), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
        f1.close()

f1 = open(os.path.join(input_script, 'allMGvcf.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    if 'jobmit' in args.job:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    else:
        f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'allMGmergevcf.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_merge, '*.sh')):
    if 'jobmit' in args.job:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    else:
        f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
