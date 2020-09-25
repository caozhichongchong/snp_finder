# start
import glob
import os

# set up path
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of folders of all genomes",
                      type=str, default='.',
                      metavar='input/')
# optional output setup

optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='pan-genome/',
                      metavar='pan-genome/')
# optional search parameters
optional.add_argument('-id',
                      help="Optional: set the identity cutoff for gene clustering",
                      metavar="95", action='store', default=95, type=int)
optional.add_argument('-cd',
                      help="Optional: set the cutoff for core gene",
                      metavar="99", action='store', default=99, type=int)
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 1)",
                      metavar="1 or more", action='store', default=40, type=int)
optional.add_argument('-job',
                      help="Optional: command to submit jobs",
                      metavar="nohup or customized",
                      action='store', default='jobmit', type=str)
# requirement for software calling
optional.add_argument('-roary',
                          help="Optional: complete path to roary if not in PATH",
                          metavar="/usr/local/bin/roary",
                          action='store', default='roary', type=str)
optional.add_argument('-snp',
                          help="Optional: complete path to snp-sites if not in PATH",
                          metavar="/usr/local/bin/snp-sites",
                          action='store', default='snp-sites', type=str)


################################################## Definition ########################################################
args = parser.parse_args()

input_script_pan = args.s + '/pangenome'
input_script = args.s
genome_root = args.i
genome_dir = glob.glob('%s/*'%(args.i))
output_dir = args.o + '/pangenome'

try:
    os.mkdir(output_dir)
except IOError:
    pass

os.system('rm -rf %s'%(input_script_pan))

try:
    os.mkdir(input_script_pan)
except IOError:
    pass


################################################## Function ########################################################
def pangenome(input_dir,species,output_dir):
    cmds = ''
    try:
        f1 = open('%s/roary_%s/core_gene_alignment.snp_sites.aln'%(output_dir,species),'r')
    except FileNotFoundError:
        cmds = '%s -p %s -o roary_%s -e -n --mafft -f %s/roary_%s -i %s -cd %s %s/*.gff\n' % (
            args.roary,args.t, species, output_dir, species, args.id, args.cd, input_dir)
        cmds += 'py37\n%s %s/roary_%s/core_gene_alignment.aln > %s/roary_%s/core_gene_alignment.snp_sites.aln\n' %(
            args.snp, output_dir, species, output_dir, species)
    return cmds

################################################## Main ########################################################
for genome_folder in genome_dir:
    #species = os.path.split(genome_folder)[-1]
    species = '%s_%s' % (os.path.split(genome_folder)[0].split('/')[-1], os.path.split(genome_folder)[-1])
    cmds = pangenome(genome_folder,species,output_dir)
    if cmds != '':
        f1 = open(os.path.join(input_script_pan, '%s.sh' % (species)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
        f1.close()

f1 = open(os.path.join(input_script, 'allpangenome.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_pan, '*.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    if 'jobmit' in args.job:
            f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    else:
            f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run: sh %s/allpangenome.sh'%(input_script))