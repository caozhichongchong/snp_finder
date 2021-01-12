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
optional.add_argument('-prokka',
                      help="Optional: complete path to prokka if not in PATH",
                      metavar="/usr/local/bin/prokka",
                      action='store', default='prokka', type=str)

################################################## Definition ########################################################
args = parser.parse_args()

input_script_pan = args.s + '/pangenome'
input_script = args.s
genome_root = args.i
species_folder = glob.glob('%s/*'%(args.i))
genome_dir = glob.glob('%s/*'%(args.i))
output_dir = args.o + '/pangenome'
output_dir_species_gff = output_dir + '/species'
fastq_name = args.fq
genome_name = args.fa
input_script_vcf = args.s + '/preassembly'

try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir_species_gff)
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
        cmds += '%s -c %s/roary_%s/core_gene_alignment.aln > %s/roary_%s/core_gene_alignment.snp_sites.aln\n' %(
            args.snp, output_dir, species, output_dir, species)
    return cmds

def prokka(output_dir, species, genome):
    genome_filename = os.path.split(genome)[-1].replace(fasta_name,'.gff')
    cmd = '%s --kingdom Bacteria --outdir %s/prokka_%s  --locustag Bacter %s\n' % \
           (args.prokka, output_dir, species, genome)
    cmd += 'mv %s/prokka_%s/*.gff %s/%s\n' % (
        output_dir, species, output_dir, genome_filename)
    cmd += 'rm -rf %s/prokka_%s\n' % (output_dir, species)
    return cmd

def runspades_single(file1,file2,output_name):
    temp_output = output_name.split(genome_name)[0]
    cmds = '%s --careful -1 %s -2 %s -o %s --threads %s --memory 100 --cov-cutoff 7\n' % \
            (args.sp,file1, file2, temp_output,args.t)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_output, output_name)
    cmds += 'rm -rf %s\n' % (temp_output)
    return cmds

def assembly_genomes(donor_species,donor_species_fastqall):
    species = donor_species.split('_')[1]
    species_folder = output_dir_species_gff + '/' + species
    try:
        os.mkdir(species_folder)
    except IOError:
        pass
    for fastq_file in donor_species_fastqall:
        if '.all' + fastq_name not in fastq_file and '.all.spades' not in fastq_file:
            original_folder, fastq_file_name = os.path.split(fastq_file)
            filename = fastq_file.split(fastq_name)[0]
            fastq_file2 = filename + fastq_name.replace('1', '2')
            newgnome_filename = fastq_file_name.split(fastq_name)[0] + genome_name
            genome_file = glob.glob(
                os.path.join(original_folder, newgnome_filename)) + \
                          glob.glob(os.path.join(original_folder + '/../',
                                                 newgnome_filename))
            prokka_output = genome_file = glob.glob(
                os.path.join(species_folder, newgnome_filename))
            cmds = ''
            if genome_file == []:
                genome_file = os.path.join(original_folder, newgnome_filename)
                # need assemble the genome first
                print('assemblying genome %s' % (genome_file))
                cmds += runspades_single(fastq_file, fastq_file2, genome_file)
            genome_file = os.path.join(original_folder, newgnome_filename)
            if prokka_output == []:
                cmds += prokka(species_folder, newgnome_filename,
                               genome_file)
            if cmds != '':
                f1 = open(os.path.join(input_script_vcf, '%s.assembly.sh' % (os.path.split(genome_file)[-1])), 'w')
                f1.write('#!/bin/bash\n%s' % (''.join(cmds)))
                f1.close()

################################################## Main ########################################################
# assembly and gff annotation
for folder in species_folder:
    donor_species = os.path.split(folder)[-1]
    donor_species_fastqall = glob.glob(os.path.join(folder, 'fastq/*' + fastq_name)) + \
                             glob.glob(os.path.join(folder, '*' + fastq_name))
    assembly_genomes(donor_species, donor_species_fastqall)

f1 = open(os.path.join(input_script, 'allpreassembly.sh'), 'w')
f1.write('#!/bin/bash\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.assembly.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    if 'jobmit' in args.job:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    elif 'nohup' in args.job:
        f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    else:
        f1.write('sh %s > %s.out\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run: sh %s/allpreassembly.sh'%(input_script))

# run roary
for genome_folder in output_dir_species_gff:
    species = os.path.split(genome_folder)[-1]
    cmds = pangenome(genome_folder,species,output_dir)
    if cmds != '':
        f1 = open(os.path.join(input_script_pan, '%s.sh' % (species)), 'w')
        f1.write('#!/bin/bash\n%s' % (''.join(cmds)))
        f1.close()

f1 = open(os.path.join(input_script, 'allpangenome.sh'), 'w')
f1.write('#!/bin/bash\n')
for sub_scripts in glob.glob(os.path.join(input_script_pan, '*.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    if 'jobmit' in args.job:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    elif 'nohup' in args.job:
        f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    else:
        f1.write('sh %s > %s.out\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run: sh %s/allpangenome.sh'%(input_script))