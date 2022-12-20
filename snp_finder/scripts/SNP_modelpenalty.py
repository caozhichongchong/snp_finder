# step 1 modelling SNPs and test 2 methods
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
                      help="path of folders of WGS of each species",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/test_data/penalty_test/',
                      metavar='input/')
# optional output setup
optional.add_argument("-fa",
                      help="file extension of fasta files",
                      type=str, default='.fasta',
                      metavar='.corrected.fasta')
optional.add_argument("-fq",
                      help="file extension of fastq files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')

optional.add_argument("-s",
                      help="a folder for your mapper",
                      type=str, default='/scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/',
                      metavar='mapper_folder/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/',
                      metavar='.')
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 40)",
                      metavar="1 or more", action='store', default=40, type=int)

################################################## Definition ########################################################
# setting match score as 0 for all
args = parser.parse_args()
penaltyset = {
    'bowtie':[6,5,3,0,1],
    'minimap':[6,4,2,0,1],
    'bwa':[5,6,1,0,1],
    'mapper':[1*10,2*10,int(0.5*10),0*10,int(0.1*10)]
}
input_script = args.s
genome_root = args.i
output_dir = args.o + '/SNP_model_penalty'
genome_name = args.fa
fastq_name=args.fq
fastq_name2=args.fq.replace('1','2')
input_script_sub = '%s/SNP_model_penalty'%(input_script)
latest_mapper = glob.glob('%s/mapper-1*.jar'%(args.s))[0]

try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir + '/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/minimap')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/bowtie')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/mapper')
except IOError:
    pass

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

def remove_unmappedreads(sam):
    cmds = 'py37\nsamtools view -F 4,256 %s > %s.mapped\n'%(sam,sam)# remove unmapped and secondary alignments
    cmds += 'mv %s.mapped %s\n'%(sam,sam)
    return cmds

def run_bowtie(files,files2,database,tempbamoutput):
    # generate code
    cmds = '#bowtie2-build %s %s\n'%(database,database)
    # not using --ma for end to end model, --ignore-quals for mismatch quality, always use the highest penalty
    cmds += 'bowtie2 --threads %s --mp %s,%s --rfg %s,%s --rdg %s,%s --np %s --ignore-quals -x %s -U %s -S %s.sam\n' % (
            min(40, args.t), penalty[0],penalty[0],penalty[1],penalty[2],penalty[1],penalty[2],penalty[4], database, files, tempbamoutput)
    cmds += remove_unmappedreads(tempbamoutput + '.sam')
    return cmds


def run_minimap(files,files2,database,tempbamoutput):
    cmds = '#minimap2 -d %s.mmi %s \n' % (database, database)
    # -A match score set as 2, -B -= 2 and --score-N += 2
    if penalty[0] > 2:
        cmds += 'minimap2 -ax sr -t %s -B %s -O %s -E %s -A %s --score-N %s %s.mmi %s >%s.sam\n' % (
                min(40, args.t), penalty[0] - 2,penalty[1],penalty[2],penalty[3] + 2 ,penalty[3]-penalty[4] + 2, database, files, tempbamoutput)
    else:
        cmds += 'minimap2 -ax sr -t %s -B %s -O %s -E %s -A %s --score-N %s %s.mmi %s >%s.sam\n' % (
                min(40, args.t), penalty[0] - 0.5,penalty[1],penalty[2],penalty[3] + 0.5 ,penalty[3]-penalty[4] + 0.5, database, files, tempbamoutput)
    cmds += remove_unmappedreads(tempbamoutput + '.sam')
    return cmds

def run_bwa(files,files2,database,tempbamoutput):
    cmds = '#bwa index %s\n'%(database)
    # -A match score set as 2, -B -= 2
    if penalty[0] > 2:
        cmds += 'bwa mem -t %s -B %s -O %s -E %s -A %s %s %s > %s.sam\n' % (
                min(40, args.t), penalty[0] - 2,penalty[1],penalty[2],penalty[3]+ 2, database, files, tempbamoutput)
    else:
        cmds += 'bwa mem -t %s -B %s -O %s -E %s -A %s %s %s > %s.sam\n' % (
            min(40, args.t), penalty[0] - 0.5, penalty[1], penalty[2], penalty[3] + 0.5, database, files,
            tempbamoutput)
    cmds += remove_unmappedreads(tempbamoutput + '.sam')
    return cmds

def run_mapper(files,files2,database,tempbamoutput):
    penalty2 = [penalty[3]+penalty[0],penalty[1],penalty[2]]
    max_penalty = penalty2[0]*0.1
    cmds = '/usr/bin/time -v java -Xms10g -Xmx10g -jar %s --max-penalty %.2f --snp-penalty %s --new-indel-penalty %s --extend-indel-penalty %s --num-threads %s --reference %s --queries %s   --out-sam %s.sam\n' % (
                latest_mapper, max_penalty, penalty2[0],penalty2[1],penalty2[2],min(40, args.t),database, files, tempbamoutput)
    return cmds

################################################## Main ########################################################
# find all files
allgenome = glob.glob('%s/*%s'%(genome_root,fastq_name))
print(allgenome)
for fastq_file in allgenome:
    databaseset = ['%s/%s'%(genome_root,os.path.split(fastq_file)[-1].replace(fastq_name,genome_name))]
    for database in databaseset:
        sample_name = os.path.split(fastq_file)[-1]
        database_name = os.path.split(database)[-1]
        # find fastq
        fastq_file2 = '%s/%s'%(genome_root,sample_name.replace(genome_name,fastq_name2))
        for tool in penaltyset:
            penalty = penaltyset[tool]
            # mapper
            tempbamoutput = '%s/%s/%s_%s.mapper1' % (output_dir, tool, sample_name,database_name)
            cmds = '#!/bin/bash\nsource ~/.bashrc\n'
            cmds += run_mapper(fastq_file, fastq_file2, database, tempbamoutput)
            # minimap
            tempbamoutput = '%s/%s/%s_%s.minimap' % (output_dir, tool, sample_name,database_name)
            cmds += run_minimap(fastq_file, fastq_file2, database, tempbamoutput)
            # bwa
            tempbamoutput = '%s/%s/%s_%s.bwa' % (output_dir, tool, sample_name,database_name)
            cmds += 'py37\n'
            cmds += run_bwa(fastq_file, fastq_file2, database, tempbamoutput)
            # bowtie
            tempbamoutput = '%s/%s/%s_%s.bowtie'%(output_dir,tool,sample_name,database_name)
            cmds += 'py39\n'
            cmds += run_bowtie(fastq_file, fastq_file2, database, tempbamoutput)
            f1 = open(os.path.join(input_script_sub, '%s.%s.%s.sh' % (sample_name,database_name,tool)), 'w')
            f1.write(cmds)
            f1.close()

f1 = open(os.path.join(input_script, 'allsnpmodel.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.sh')):
    f1.write('jobmit %s %s small\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
f1.close()

################################################### END ########################################################
