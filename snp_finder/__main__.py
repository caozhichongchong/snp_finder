import os
import argparse
import glob
import snp_finder

################################################### Decalration #######################################################
print ("\
------------------------------------------------------------------------\n\
snp_finder identifies and summarizes SNPs in clonal populations\n\
input: folders of WGS (whole genome sequences) of each species\n\
requirement: samtools, bcftools, spades, minimap2, usearch \n\n\
Copyright: An-Ni Zhang, Prof. Eric Alm, MIT\n\n\
Citation:\n\
Contact anniz44@mit.edu\n\
------------------------------------------------------------------------\n\
")

def main():
    usage = ("usage: snp_finder -t your.otu.table -s your.otu.seqs")
    version_string = 'snp_finder {v}, on Python {pyv[0]}.{pyv[1]}.{pyv[2]}'.format(
        v=snp_finder.__version__,
        pyv=sys.version_info,
    )
    ############################################ Arguments and declarations ##############################################
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('command',
                        help="snp_finder prepare for preparing and correcting assembled genomes for each WGS, \
                        snp_finder snp for identifying and sumarizing SNPs, \
                        snp_finder select for calculating dN/dS ratio and selecting/extracting genes under parallel evolution, \
                        snp_finder anno for annotating genes under parallel evolution and genes with de novo mutaitons,\
                        snp_finder meta for tracking strains in metagenomes",
                        type=str,
                        default='prepare',
                        choices=['prepare','snp', 'select',
                        'anno','meta'],
                        action='store',
                        metavar='snp_finder command')
    required.add_argument("-i",
                        help="path of folders of WGS of each species",
                        type=str, default='.',
                        metavar='input/')
    required.add_argument("-fq",
                          help="file extension of WGS fastq #1 files",
                          type=str, default='_1.fastq',
                          metavar='_1.fastq')
    # optional input/output setup
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
    # requirement for software calling
    optional.add_argument('-sp', '--spades',
                          help="Optional: complete path to spades if not in PATH",
                          metavar="/usr/local/bin/spades",
                          action='store', default='spades', type=str)
    optional.add_argument('-bw', '--bowtie',
                          help="Optional: complete path to bowtie if not in PATH",
                          metavar="/usr/local/bin/bowtie",
                          action='store', default='bowtie', type=str)
    optional.add_argument('-pro','--prodigal',
                        help="Optional: complete path to prodigal if not in PATH",
                        metavar="/usr/local/bin/prodigal",
                        action='store', default='prodigal', type=str)
    optional.add_argument('-bcf','--bcftools',
                        help="Optional: complete path to bcftools if not in PATH",
                        metavar="/usr/local/bin/bcftools",
                        action='store', default='bcftools', type=str)
    optional.add_argument('-sam','--samtools',
                        help="Optional: complete path to bwa if not in PATH",
                        metavar="/usr/local/bin/samtools",
                        action='store', default='samtools', type=str)
    optional.add_argument('-mini', '--minimap2',
                          help="Optional: complete path to minimap2 if not in PATH",
                          metavar="/usr/local/bin/minimap2",
                          action='store', default='minimap2', type=str)
    optional.add_argument('--u','--usearch',
                        help="Optional: cluster genes with SNPs",
                        metavar="usearch",
                        action='store', default='usearch', type=str)

    ################################################## Definition ########################################################
    args = parser.parse_args()
    workingdir=os.path.abspath(os.path.dirname(__file__))
    jobmit = 'jobmit'
    ################################################### Programme #######################################################
    f1 = open (args.s + '/snp_finder.sh','w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    thread = int(args.t)
    # run traits finding
    if args.command == 'prepare':
        cmd = ('python '+workingdir+'/scripts/prepare_genome.py -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s\n'
               % (args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf, args.sam, args.mini))
        f1.write(cmd)
        f1.write('sh %s\n'%(os.path.join(input_script, 'allcurate.sh')))
        cmd = ('python ' + workingdir + '/scripts/genome_curate.py -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s\n'
                    % (
                    args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf, args.sam,
                    args.mini))
        f1.write(cmd)
        cmd = ('python ' + workingdir + '/scripts/genome_curate2.py -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini))
        f1.write(cmd)
        f1.write('sh %s\n' % (os.path.join(input_script, 'cleanoldgenome.sh')))
        f1.write('sh %s\n'%(os.path.join(input_script, 'allcurate.sh')))
    if args.command == 'snp':
        cmd = ('python ' + workingdir + '/scripts/SNPfilter_genome.py -rd 1 -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s\n'
                    % (
                    args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf, args.sam,
                    args.mini))
        f1.write(cmd)
        f1.write('sh %s\n' % (os.path.join(input_script, 'SNP_round1.move.sh')))
        cmd = ('python ' + workingdir + '/scripts/mapping_genome.py -rd 2 -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini))
        f1.write(cmd)
        f1.write('sh %s\n' % (os.path.join(input_script, 'allround2.sh')))
        cmd = ('python ' + workingdir + '/scripts/SNPfilter_genome.py -rd 2 -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini))
        f1.write(cmd)
        f1.write('sh %s\n' % (os.path.join(input_script, 'SNP_round2.move.sh')))
        cmd = ('python ' + workingdir + '/scripts/mapping_genome.py -rd 3 -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini))
        f1.write(cmd)
        f1.write('sh %s\n' % (os.path.join(input_script, 'allround3.sh')))
        cmd = ('python ' + workingdir + '/scripts/SNPfilter_genome.py -rd 3 -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini))
        f1.write(cmd)
        f1.write('sh %s\n' % (os.path.join(input_script, 'SNP_round3.move.sh')))
        cmd = ('python ' + workingdir + '/scripts/mapping_genome.py -rd 4 -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s\n'
                % (
                    args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                    args.sam,
                    args.mini))
        f1.write(cmd)
        f1.write('sh %s\n' % (os.path.join(input_script, 'allround4.sh')))
        cmd = ('python ' + workingdir + '/scripts/SNPfilter_genome.py -rd 4 -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini))
        f1.write(cmd)
        cmd = ('python ' + workingdir + '/scripts/SNPfilter_WGS.py -rd 4 -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini))
        f1.write(cmd)
    if args.command == 'select':
        cmd = ('python ' + workingdir + '/scripts/dnds.py -rd 4 -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini))
        f1.write(cmd)
        cmd = ('python ' + workingdir + '/scripts/extract.py -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s --u %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini, args.u))
        f1.write(cmd)
        cmd = ('python ' + workingdir + '/scripts/parallel_evolution.py -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s --u %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini, args.u))
        f1.write(cmd)
    if args.command == 'meta':
        cmd = ('python ' + workingdir + '/scripts/mapping_meta.py -i %s -fq %s -s %s -o %s -m %s -mfq %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s --u %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.m, args.mfq, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini, args.u))
        f1.write(cmd)
        f1.write('sh %s\n' % (os.path.join(input_script, 'allMGvcf.sh')))
        f1.write('sh %s\n' % (os.path.join(input_script, 'allMGmergevcf.sh')))
        cmd = ('python ' + workingdir + '/scripts/SNPfilter_meta.py -i %s -fq %s -s %s -o %s -m %s -mfq %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s --u %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.m, args.mfq, args.t, args.job, args.bw, args.sp, args.pro,
                        args.bcf,
                        args.sam,
                        args.mini, args.u))
        f1.write(cmd)
    if args.command == 'anno':
        cmd = ('python ' + workingdir + '/scripts/annotate.py -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s --u %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini, args.u))
        f1.write(cmd)
        f1.write('sh %s\n' % (os.path.join(input_script, 'allannotate.sh')))
        f1.write('sh %s\n' % (os.path.join(input_script, 'allannotate_all.sh')))
        cmd = ('python ' + workingdir + '/scripts/annotate_sum.py -i %s -fq %s -s %s -o %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s --u %s\n'
                    % (
                        args.i, args.fq, args.s, args.o, args.t, args.job, args.bw, args.sp, args.pro, args.bcf,
                        args.sam,
                        args.mini, args.u))
        f1.write(cmd)

################################################## Function ########################################################

if __name__ == '__main__':
    main()
