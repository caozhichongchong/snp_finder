# start
# metagenomes mapping to a reference genome
from Bio import SeqIO
from Bio.Seq import Seq
import glob
import os
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
required.add_argument("-snp",
                      help="path of snp files",
                      type=str, default='None',
                      metavar='/scratch/users/anniz44/genomes/donor_species/WGS/vcf_round1/merge/')
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
                      metavar="1 or more", action='store', default=40, type=int)
optional.add_argument('-job',
                      help="Optional: command to submit jobs",
                      metavar="nohup or customized",
                      action='store', default='jobmit', type=str)
# requirement for software calling
optional.add_argument('-dm',
                          help="Optional: complete path to diamond if not in PATH",
                          metavar="/usr/local/bin/diamond",
                          action='store', default='diamond', type=str)
optional.add_argument('-bcf',
                      help="Optional: complete path to bcftools if not in PATH",
                      metavar="/usr/local/bin/bcftools",
                      action='store', default='bcftools', type=str)
optional.add_argument('-sam',
                      help="Optional: complete path to samtools if not in PATH",
                      metavar="/usr/local/bin/samtools",
                      action='store', default='samtools', type=str)
optional.add_argument('-mini',
                      help="Optional: complete path to minimap2 if not in PATH",
                      metavar="/usr/local/bin/minimap2",
                      action='store', default='minimap2', type=str)

################################################## Definition ########################################################
args = parser.parse_args()

input_script = args.s
input_script_sub = input_script + '/vcf'
genome_dir = glob.glob(args.i + '/*')
metagenome_dir = args.m
output_dir = args.o + '/MG'
fastq_name = args.mfq
database_AMPHORA = args.d
try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir + '/bwa')
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

allfastq = glob.glob(os.path.join(metagenome_dir,'*%s'%(fastq_name)))
print(allfastq)
alldatabase = []
# find all database
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    donor_species_genome = glob.glob(os.path.join(folder, '*.all.spades*.fasta.noHM.fasta'))
    sub_samples = []
    if donor_species_genome != []:
        database = donor_species_genome[0]
        alldatabase.append(database)
Total_job = max(int(len(allfastq)*len(alldatabase)/250),1)

# map metagenome to each genome
i = 0
cmds = ''
cmdssum = ''
for database in alldatabase:
    donor_species = os.path.split(database)[-1].split('.all')[0]
    try:
        f1 = open(database + '.mmi')
    except IOError:
        os.system('minimap2 -d %s.mmi %s \n' % (database, database))
    for files in allfastq:
        donor_species_dir_file = os.path.split(files)[-1]
        tempbamoutput = os.path.join(output_dir + '/bwa', donor_species_dir_file + '.' + donor_species)
        try:
            f1 = open(tempbamoutput + '.raw.vcf')
            try:
                f1 = open(tempbamoutput + '.raw.vcf.depth.MV.sum')
            except IOError:
                cmdssum += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta.py -vcf %s -snp %s -mfq %s -o %s \n' % (
                tempbamoutput + '.raw.vcf',
                args.snp, args.mfq, args.o)
        except IOError:
            try:
                tempbamoutputfinish = os.path.join(output_dir + '/bwa/finished', donor_species_dir_file + '.' + donor_species)
                f1 = open(tempbamoutputfinish + '.raw.vcf.zip')
            except IOError:
                if 'fasta' in fastq_name:
                    cmds += args.mini + ' -ax sr -N 1000 -p 0.99 -t %s %s.mmi %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                        min(40, args.t), database, files, args.sam, min(40, args.t),
                        tempbamoutput, args.sam, min(40, args.t), tempbamoutput, tempbamoutput, args.sam,
                        min(40, args.t),
                        tempbamoutput)
                else:
                    cmds += args.mini + ' -ax sr -N 1000 -p 0.99 -t %s %s.mmi %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                        min(40, args.t), database, files, args.sam, min(40, args.t),
                        tempbamoutput, args.sam, min(40, args.t), tempbamoutput, tempbamoutput, args.sam,
                        min(40, args.t),
                        tempbamoutput)
                cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q10 -Ou -B -d3000 -f %s %s.sorted.bam | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
                    args.bcf, min(40, args.t), database,
                    tempbamoutput, args.bcf, min(40, args.t), tempbamoutput)
                cmds += 'rm -rf %s.bam %s.bam.bai %s.sorted.bam %s.sorted.bam.bai\n' % (
                    tempbamoutput, tempbamoutput, tempbamoutput, tempbamoutput)
                try:
                    f1 = open(tempbamoutput + '.raw.vcf.depth.MV.sum')
                except IOError:
                    cmds += 'py37\npython /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta.py -vcf %s -snp %s -mfq %s -o %s \n'%(tempbamoutput + '.raw.vcf',
                                                                                                                                               args.snp, args.mfq, args.o)
                i += 1
                if i % Total_job == 0:
                    f1 = open(os.path.join(input_script_sub, '%s.vcf.sh' % int(i / Total_job)), 'a')
                    f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s' % (''.join(cmds)))
                    f1.close()
                    cmds = ''
    f1 = open(os.path.join(input_script_sub, '%s.vcf.sh' % int(i / Total_job)), 'a')
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

f1 = open(os.path.join(input_script, 'allMGsum.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n')
f1.write(''.join(cmdssum))
f1.close()
################################################### END ########################################################
