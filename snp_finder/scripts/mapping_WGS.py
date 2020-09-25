# start
# round 1-3 run genome assembly and map genomes to a reference genome
import glob
import os
import statistics
from Bio import SeqIO
from Bio.Seq import Seq
# set up path
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
required.add_argument("-cl",
                      help="path to output of step 1 clustering clonal population (None if there's no clustering)",
                      type=str, default='.',
                      metavar='pan-genome/')
# optional output setup

optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='WGS/',
                      metavar='WGS/')
optional.add_argument("-fa",
                      help="file extension of assembled genomes",
                      type=str, default='.fasta.corrected.fasta',
                      metavar='.fasta.corrected.fasta')
optional.add_argument("-ref",
                      help="reference genome list if not co-assembly \n in the format of \'folder of species\'\t\'reference_genome_path\'",
                      type=str, default='None',
                      metavar='reference.genome.txt')
# optional search parameters
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 1)",
                      metavar="1 or more", action='store', default=40, type=int)
optional.add_argument('-rd',
                      help="Round of SNP calling and filtering",
                      metavar="1-4", action='store', default=1, type=int)
optional.add_argument('-job',
                      help="Optional: command to submit jobs",
                      metavar="nohup or customized",
                      action='store', default='jobmit', type=str)
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

################################################## Definition ########################################################
args = parser.parse_args()
Round = args.rd
input_script_vcf = args.s + '/vcf_round%s'%(Round)
input_script = args.s
genome_root = args.i
species_folder = glob.glob('%s/*'%(args.i))
output_dir = args.o + '/vcf_round%s'%(Round)
fastq_name = args.fq
genome_name = args.fa
# co-assembly cutoff
cluster_cutoff = 3 # at least 3 genomes for SNP calling
genomesize_cutoff = 1.2 # remove genomes with >=120% the mean size
# homologous cutoff
length_cutoff = 500 # 500 bp for homologous region
identity_cutoff = 95 # 95% identity for homologous region
CHR_length_cutoff = 2000 # minimum contig lengths for reference genome
workingdir=os.path.abspath(os.path.dirname(__file__))

try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/co-assembly')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/merge')
except IOError:
    pass

try:
    os.mkdir(input_script)
except IOError:
    pass

os.system('rm -rf %s'%(input_script_vcf))

try:
    os.mkdir(input_script_vcf)
except IOError:
    pass
################################################## Function ########################################################

def subset(file1,file2,output_name1,output_name2):
    cmds = 'head -400000 %s >> %s\n' % (
        file1, output_name1)
    cmds += 'tail -400000 %s >> %s\n' % (
        file1, output_name1)
    cmds += 'head -400000 %s >> %s\n' % (
        file2, output_name2)
    cmds += 'tail -400000 %s >> %s\n' % (
        file2, output_name2)
    return cmds

def runspades(file1,file2,temp_output,output_name):
    cmds = '%s --careful -1 %s -2 %s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' % \
            (args.sp,file1, file2, temp_output)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_output, output_name)
    cmds += 'rm -rf %s %s %s\n' % (temp_output, file1, file2)
    return cmds

def runspades_single(file1,file2,output_name):
    temp_output = output_name.split(genome_name)[0]
    cmds = '%s --careful -1 %s -2 %s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' % \
            (args.sp,file1, file2, temp_output)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_output, output_name)
    cmds += 'rm -rf %s\n' % (temp_output)
    return cmds

def run_vcf_WGS(files,files2,database,tempbamoutput):
    # generate code
    cmds = args.bw + ' --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, args.t), database, files, files2, args.sam, min(40, args.t),
        tempbamoutput, args.sam, min(40, args.t), tempbamoutput, tempbamoutput, args.sam, min(40, args.t),
        tempbamoutput)
    cmds += 'rm -rf %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    sample_set3 = ['%s.sorted.bam' % (tempbamoutput)]
    return [cmds, sample_set3]

def run_vcf(genome_file,database,tempbamoutput):
    # generate code
    # for curated genome
    cmds = args.mini + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, args.t), database, genome_file, args.sam, min(40, args.t),
        tempbamoutput, args.sam, min(40, args.t), tempbamoutput, tempbamoutput, args.sam, min(40, args.t),
        tempbamoutput)
    cmds += '#samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
        tempbamoutput, tempbamoutput)
    sample_set3 = ['%s.sorted.bam' % (tempbamoutput)]
    cmds += 'rm -rf %s.bam %s.bam.bai\n' %(tempbamoutput,tempbamoutput)
    return [cmds,sample_set3]

def merge_sample(database,tempbamoutputGenome,samplesetGenomecorrect):
    # corrected genomes
    cmds = '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
        args.bcf, min(40, args.t), database,
        ' '.join(samplesetGenomecorrect), args.bcf, min(40, args.t), tempbamoutputGenome)
    cmds += '%s view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
        args.bcf, tempbamoutputGenome, tempbamoutputGenome)
    return cmds

def remove_homologous(genome):
    cmds = 'py37\npython %s/remove_homologous.py -i %s\n'%(workingdir,genome)
    return ['%s.noHM.fasta' %(genome),cmds]

def mapping_WGS(donor_species,donor_species_fastqall):
    try:
        f1 = open(os.path.join(output_dir + '/merge', donor_species + '.all.flt.snp.vcf'), 'r')
    except FileNotFoundError:
        print('generate mapping code for %s' % (donor_species))
        folder = output_dir + '/co-assembly/%s' % (donor_species)
        try:
            os.mkdir(folder)
        except IOError:
            pass
        cmds = ''
        if len(donor_species_fastqall) >= cluster_cutoff:
            # check assmebly size
            allgenomesize = []
            for fastq_file in donor_species_fastqall:
                if '.all' + fastq_name not in fastq_file and '.all.spades' not in fastq_file:
                    original_folder, fastq_file_name = os.path.split(fastq_file)
                    filename = fastq_file.split(fastq_name)[0]
                    fastq_file2 = filename + fastq_name.replace('1', '2')
                    genome_file = glob.glob(
                        os.path.join(original_folder, fastq_file_name.split(fastq_name)[0] + genome_name)) + \
                                  glob.glob(
                                      os.path.join(original_folder + '/../',
                                                   fastq_file_name.split(fastq_name)[0] + genome_name))
                    if genome_file == []:
                        # need assemble the genome first
                        print('assemblying genome %s' % (genome_file))
                        os.system(runspades_single(fastq_file, fastq_file2, genome_file))
                    genome_file = genome_file[0]
                    filesize = int(os.path.getsize(genome_file))
                    allgenomesize.append(filesize)
                else:
                    allgenomesize.append(0)
            # filter out genomesize that's above the size cutoff
            maxgenome = statistics.mean(allgenomesize) * genomesize_cutoff
            qualifygenome = [donor_species_fastqall[i] for i in range(0, len(allgenomesize)) if
                             allgenomesize[i] <= maxgenome and allgenomesize[i] > 0]
            if donor_species not in Ref_set:
                # co-assembly
                donor_species_folder_all = os.path.join(folder, donor_species + '_allspades%s' % (Round))
                donor_species_genomename = os.path.join(folder, donor_species + '.all.spades%s.fasta' % (Round))
                donor_species_fastq = os.path.join(folder, donor_species + '.all.spades%s' % (Round) + fastq_name)
                donor_species_fastq2 = os.path.join(folder,
                                                    donor_species + '.all.spades%s' % (Round) + fastq_name.replace('1',
                                                                                                                   '2'))
                try:
                    f1 = open(donor_species_genomename, 'r')
                except FileNotFoundError:
                    print('run subset fastq for %s' % (donor_species_genomename))
                    for fastq_file in qualifygenome:
                        filename = fastq_file.split(fastq_name)[0]
                        fastq_file2 = filename + fastq_name.replace('1', '2')
                        # subset fastq
                        cmds += subset(fastq_file, fastq_file2, donor_species_fastq, donor_species_fastq2)
                    # pan genome
                    cmds += runspades(donor_species_fastq, donor_species_fastq2, donor_species_folder_all,
                                     donor_species_genomename)
                    print('run co-assembly for %s' % (donor_species_genomename))
            else:
                # use reference genome
                donor_species_genomename = Ref_set[donor_species]
            if 'noHM.fasta' not in donor_species_genomename:
                donor_species_genomename,cmds2 = remove_homologous(donor_species_genomename)
                cmds += cmds2
            tempbamoutputGenome = os.path.join(output_dir + '/merge', donor_species + '.all')
            vcf_output = []
            # then run library
            cmds += args.mini + ' -d %s.mmi %s\n' % (donor_species_genomename, donor_species_genomename)
            cmds += 'rm -rf %s.fai\n' % (donor_species_genomename)
            cmds += os.path.join(os.path.split('args.bw')[0], 'bowtie2-build') + ' %s %s\n' % (
                donor_species_genomename, donor_species_genomename)
            for fastq_file in qualifygenome:
                fastq_file_name = os.path.split(fastq_file)[-1]
                filename = fastq_file.split(fastq_name)[0]
                fastq_file2 = filename + fastq_name.replace('1', '2')
                # map each WGS to pangenome
                results = run_vcf_WGS(fastq_file, fastq_file2,
                                      donor_species_genomename,
                                      os.path.join(output_dir + '/bwa',
                                                   fastq_file_name))
                vcf_output += results[1]
                cmds += results[0]
            cmds += merge_sample(donor_species_genomename, tempbamoutputGenome, vcf_output)
            f1 = open(os.path.join(input_script_vcf, '%s.vcf.sh' % (donor_species)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            f1.close()

################################################## Main ########################################################
# load reference genomes if provided
Ref_set = dict()
if args.ref != 'None':
    for lines in open(args.ref,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        folder, ref = lines_set[0:2]
        Ref_set.setdefault(folder, ref)

# load clustering results
if args.cl != 'None':
    Cluster = dict()
    allclustering = glob.glob(args.cl + '/clonal_population/*.genome.cluster.txt') + \
                    glob.glob(args.cl + '/*.genome.cluster.txt')
    for clusterfile in allclustering:
        species = os.path.split(clusterfile)[-1].split('.genome.cluster.txt')[0]
        print('find all fastq for each cluster of %s' % (species))
        for lines in open(clusterfile,'r'):
            if not lines.startswith('species'):
                lines_set = lines.split('\n')[0].split('\t')
                genome, cluster, total_cluster = lines_set[1:4]
                Cluster.setdefault(species, dict())
                Cluster[species].setdefault(cluster, [])
                fastq_file = glob.glob(os.path.join(args.i + '/*/fastq/', '%s%s' % (genome, fastq_name)))# + \
                #glob.glob(os.path.join(args.i + '/*/', '%s%s' % (genome, fastq_name)))
                if fastq_file!=[]:
                    Cluster[species][cluster].append(fastq_file[0])
    # generate code
    for species in Cluster:
        for cluster in Cluster[species]:
            donor_species = '%s_cluster%s'%(species,cluster)
            donor_species_fastqall = Cluster[species][cluster]
            print(donor_species_fastqall)
            mapping_WGS(donor_species, donor_species_fastqall)
else:
    # generate code based on donor-species
    for folder in species_folder:
        donor_species = os.path.split(folder)[-1]
        donor_species_fastqall = glob.glob(os.path.join(folder, 'fastq/*' + fastq_name)) + \
                                 glob.glob(os.path.join(folder, '*' + fastq_name))
        mapping_WGS(donor_species, donor_species_fastqall)

# sum all codes
f1 = open(os.path.join(input_script, 'allWGS.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.vcf.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    if 'jobmit' in args.job:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    else:
        f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run: sh %s/allWGS.sh'%(input_script))
################################################### END ########################################################
