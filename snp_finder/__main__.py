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
    usage = ("usage: snp_finder auto -i input_folder")
    version_string = 'snp_finder {v}, on Python {pyv[0]}.{pyv[1]}.{pyv[2]}'.format(
        v=snp_finder.__version__,
        pyv=sys.version_info,
    )
    ############################################ Arguments and declarations ##############################################
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('command',
                        help="snp_finder auto for automatically running all scripts, \
                        snp_finder manual for manually running the scripts (fast and flexible)",
                        type=str,
                        default='manual',
                        choices=['manual','auto'],
                        action='store',
                        metavar='snp_finder command')
    required.add_argument("-i",
                          help="path of folders of WGS of each species",
                          type=str, default='.',
                          metavar='input/')
    required.add_argument("-fq",
                          help="file extension of WGS fastq #1 files",
                          type=str, default='_1.fq',
                          metavar='_1.fastq')
    # optional output setup
    optional.add_argument("-s",
                          help="a folder to store all scripts",
                          type=str, default='snpfinder_scripts/',
                          metavar='snpfinder_scripts/')
    optional.add_argument("-o",
                          help="a folder to store all output",
                          type=str, default='snpfinder_output/',
                          metavar='snpfinder_output/')
    # optional input support
    optional.add_argument("-fa",
                          help="file extension of assembled genomes if assembled",
                          type=str, default='.fasta',
                          metavar='.fasta')
    optional.add_argument("-ref",
                          help="a list of reference genome if not using co-assembly\n in the format of \'folder_name\'\t\'reference_genome_path\'",
                          type=str, default='None',
                          metavar='reference.genome.txt')
    optional.add_argument("-trunc",
                          help="whether to extract truncated genes caused by SNPs",
                          type=str, default='False',
                          metavar='False or True')
    # optional PE setup
    optional.add_argument("-cutoff",
                          help="a file to set cutoff for how many SNPs on a gene to call PE in lineages (default: 2 SNPs)",
                          type=str, default='None',
                          metavar='SNP_cutoff_lineage.txt')
    optional.add_argument("-cutoffsp",
                          help="a file to set cutoff for how many SNPs on a gene to call PE in species (default: 2 SNPs)",
                          type=str, default='None',
                          metavar='SNP_cutoff_species.txt')
    # optional search parameters
    optional.add_argument('-t',
                          help="Optional: set the thread number assigned for running XXX (default 1)",
                          metavar="1 or more", action='store', default=1, type=int)
    optional.add_argument('-job',
                          help="Optional: command to submit jobs if not using nohup",
                          metavar="nohup or customized",
                          action='store', default='nohup', type=str)
    # requirement for software calling
    optional.add_argument('-bw',
                          help="Optional: complete path to bowtie2 if not in PATH",
                          metavar="/usr/local/bin/bowtie2",
                          action='store', default='bowtie2', type=str)
    optional.add_argument('-sp',
                          help="Optional: complete path to spades.py if not in PATH",
                          metavar="/usr/local/bin/spades.py",
                          action='store', default='spades.py', type=str)
    optional.add_argument('-pro',
                          help="Optional: complete path to prodigal if not in PATH, None for no prodigal (default)",
                          metavar="/usr/local/bin/prodigal",
                          action='store', default='None', type=str)
    optional.add_argument('-bcf',
                          help="Optional: complete path to bcftools if not in PATH",
                          metavar="/usr/local/bin/bcftools",
                          action='store', default='bcftools', type=str)
    optional.add_argument('-sam',
                          help="Optional: complete path to bwa if not in PATH",
                          metavar="/usr/local/bin/samtools",
                          action='store', default='samtools', type=str)
    optional.add_argument('-mini',
                          help="Optional: complete path to minimap2 if not in PATH",
                          metavar="/usr/local/bin/minimap2",
                          action='store', default='minimap2', type=str)
    optional.add_argument('-blastn',
                          help="Optional: complete path to blastn if not in PATH",
                          metavar="/usr/local/bin/blastn",
                          action='store', default='blastn', type=str)
    optional.add_argument('-u',
                          help="Optional: complete path to usearch if not in PATH",
                          metavar="/usr/local/bin/usearch",
                          action='store', default='usearch', type=str)
    optional.add_argument('-roary',
                          help="Optional: complete path to roary if not in PATH",
                          metavar="/usr/local/bin/roary",
                          action='store', default='roary', type=str)
    optional.add_argument('-snp',
                          help="Optional: complete path to snp-sites if not in PATH",
                          metavar="/usr/local/bin/snp-sites",
                          action='store', default='snp-sites', type=str)
    optional.add_argument('-prokka',
                          help="Optional: complete path to prokka if not in PATH",
                          metavar="/usr/local/bin/prokka",
                          action='store', default='prokka', type=str)

    ################################################## Definition ########################################################
    args = parser.parse_args()
    workingdir=os.path.abspath(os.path.dirname(__file__))
    ################################################### Programme #######################################################
    f1 = open (args.s + '/snp_finder.sh','w')
    f1.write('#!/bin/bash\n')
    job = args.job
    # run snp finder
    if args.command == 'auto':
        job = 'None'
        cmd = ''
    else:
        cmd = 'Please run the following scripts step by step\nPlease make sure all tasks of one step are finished before running the next step\n'
    # step 1 clonal population clustering
    cmd += '#step 1 clonal population clustering\n'
    input_script = args.s
    genome_root = args.i
    output_dir = args.o
    cmd += ('python '+workingdir+'/bin/pangenome.py -i %s -s %s -o %s -t %s -job %s -roary %s -snp %s -prokka %s -fq %s\n' % (
        genome_root, input_script, output_dir, args.t, job, args.roary, args.snp, args.prokka, args.fq))
    cmd += ('sh %s/allprepangenome.sh\n' % (input_script))
    cmd += ('python ' + workingdir + '/bin/pangenome.py -i %s -s %s -o %s -t %s -job %s -roary %s -snp %s -prokka %s\n' % (
            genome_root, input_script, output_dir, args.t, job, args.roary, args.snp, args.prokka))
    cmd += ('sh %s/allpangenome.sh\n' % (input_script))
    cmd += ('python '+workingdir+'/bin/cluster_CP.py -i %s -s %s -o %s -t %s -roary %s -snp %s\n' % (
        genome_root, input_script, output_dir, args.t, args.roary, args.snp))
    cmd += ('python '+workingdir+'/bin/cluster_CP.py -i %s -s %s -o %s -clustering 2 -t %s -roary %s -snp %s\n' % (
        genome_root, input_script, output_dir, args.t, args.roary, args.snp))
    # step 2 co-assembly and mapping
    cmd += '\n#step 2 co-assembly and mapping\n'
    pangenome_dir = output_dir + '/pangenome/'
    cmd += ('python '+workingdir+'/bin/assembly_genome.py -i %s -s %s -o %s -fa %s -ref %s -cl %s -fq %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s -blastn %s\n' % (
        genome_root, input_script,output_dir,args.fa, args.ref,pangenome_dir, args.fq, args.t, job, args.bw, args.sp, args.pro, args.bcf, args.sam, args.mini, args.blastn))
    cmd += ('sh %s/allassembly.sh\n' % (input_script))
    cmd +=('python '+workingdir+'/bin/mapping_WGS.py -i %s -s %s -o %s -fa %s -ref %s -cl %s -fq %s -t %s -job %s -bw %s -sp %s -pro %s -bcf %s -sam %s -mini %s -blastn %s\n' % (
        genome_root, input_script,output_dir, args.fa, args.ref,pangenome_dir,  args.fq, args.t, job, args.bw, args.sp, args.pro, args.bcf, args.sam, args.mini, args.blastn))
    cmd +=('sh %s/allWGS.sh\n' % (input_script))
    cmd +=('python '+workingdir+'/bin/vcf_process.py -i %s -s %s -o %s -fq %s -t %s -pro %s\n' % (
        genome_root, input_script, output_dir, args.fq, args.t, args.pro))
    # step 3 remove recombination and low-quality SNPs
    cmd += '\n#step 3 remove recombination and low-quality SNPs\n'
    vcf_folder = output_dir + '/vcf_round1/merge'
    vcf_format = '.filtered.vcf'
    cmd += ('python SNPfilter_WGS.py -i %s -vcf %s -s %s -fq %s -pro %s\n' % (
        vcf_folder, vcf_format, input_script, args.fq, args.pro))
    # step 4 calculate dnds
    cmd += '\n#step 4 calculate dnds\n'
    output_dir_merge = output_dir + '/vcf_round1/merge'
    cmd += ('python dnds.py -s %s -o %s -cutoff %s\n' % (
        input_script, output_dir_merge, args.cutoff))
    # step 5 screen for parallel evolution
    cmd += '\n#step 5 screen for parallel evolution\n'
    co_assembly_dir = output_dir + '/vcf_round1/co-assembly/'
    output_dir_merge = output_dir + '/vcf_round1/merge/'
    cmd += ('python parallel_evolution.py -i %s -s %s -o %s -cutoff %s -trunc %s -t %s -pro %s -u %s\n' % (
    co_assembly_dir, input_script, output_dir_merge, args.cutoffsp, args.trunc, args.t, args.pro, args.u))
    cmd += ('sh %s/%s\n' % (input_script, 'allannotate.sh'))
    cmd += ('sh %s/%s\n' % (input_script, 'allannotate_all.sh'))
    f1.write(cmds)
    f1.close()
    print('please run %s'%(args.s + '/snp_finder.sh'))
    if args.command == 'manual':
        print('Please run the following scripts step by step\nPlease make sure all tasks of one step are finished before running the next step\n')

################################################## Function ########################################################

if __name__ == '__main__':
    main()
