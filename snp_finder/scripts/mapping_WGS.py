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
                      metavar="1-4", action='store', default=2, type=int)
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
optional.add_argument('--prokka',
                        help="Optional: complete path to prokka if not in PATH",
                        metavar="/usr/local/bin/prokka",
                        action='store', default='prokka', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
Folder_num = 7 # 7 for PB, BA, BL; 6 for BaFr
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
    os.mkdir(args.o)
except IOError:
    pass

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
    cmds = ''
    try:
        f1 = open('%s.sorted.bam' % (tempbamoutput),'r')
    except IOError:
        cmds = args.bw + ' --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
            min(40, args.t), database, files, files2, args.sam, min(40, args.t),
            tempbamoutput, args.sam, min(40, args.t), tempbamoutput, tempbamoutput, args.sam, min(40, args.t),
            tempbamoutput)
        cmds += 'rm -rf %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def merge_sample(database,vcfoutput,allsam,coverage_only = False):
    cmds = ''
    try:
        f1 = open('%s.raw.vcf' % (vcfoutput), 'r')
    except FileNotFoundError:
        cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
            args.bcf, min(40, args.t), database,
            ' '.join(allsam), args.bcf, min(40, args.t), vcfoutput)
    if not coverage_only:
        try:
            f1 = open('%s.flt.snp.vcf' % (vcfoutput))
        except FileNotFoundError:
            cmds += '%s view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
                args.bcf, vcfoutput, vcfoutput)
    return cmds

def merge_sample_subcluster(database,vcfoutput,allsam,coverage_only = False):
    cmds = ''
    allsam = list(set(allsam))
    if len(allsam) <= 50:
        try:
            f1 = open('%s.raw.vcf'%(vcfoutput),'r')
        except FileNotFoundError:
            cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
                args.bcf, min(40, args.t), database,
                ' '.join(allsam), args.bcf, min(40, args.t), vcfoutput)
        if not coverage_only:
            try:
                f1 = open('%s.flt.snp.vcf' % (vcfoutput))
            except FileNotFoundError:
                cmds += '%s view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
                    args.bcf, vcfoutput, vcfoutput)
    else:
        total_folder = int(len(allsam) / 50)
        print('creating %s subclusters for %s samples'%(total_folder + 1,len(allsam)))
        for i in range(0, total_folder + 1):
            vcfoutput2 = vcfoutput + '_subcluster%s'%(i)
            if i == total_folder:
                try:
                    f1 = open('%s.raw.vcf' % (vcfoutput2),'r')
                except FileNotFoundError:
                    cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
                        args.bcf, min(40, args.t), database,
                        ' '.join(allsam[i * 50:]), args.bcf, min(40, args.t), vcfoutput2)
                if not coverage_only:
                    try:
                        f1 = open('%s.flt.snp.vcf' % (vcfoutput2))
                    except FileNotFoundError:
                        cmds += '%s view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
                            args.bcf, vcfoutput2, vcfoutput2)
            else:
                try:
                    f1 = open('%s.raw.vcf' % (vcfoutput2))
                except FileNotFoundError:
                    cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
                        args.bcf, min(40, args.t), database,
                        ' '.join(allsam[i * 50:(i + 1) * 50]), args.bcf, min(40, args.t), vcfoutput2)
                if not coverage_only:
                    try:
                        f1 = open('%s.flt.snp.vcf' % (vcfoutput2))
                    except FileNotFoundError:
                        cmds += '%s view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
                            args.bcf, vcfoutput2, vcfoutput2)
    return cmds

def remove_homologous(genome):
    cmds = ''
    genome_noHM = '%s.noHM.fasta' %(genome)
    try:
        f1 = open(genome_noHM,'r')
    except FileNotFoundError:
        cmds = 'py37\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:$LD_LIBRARY_PATH\npython %s/remove_homologous.py -i %s\n'%(workingdir,genome)
        # then run library
        #cmds += args.mini + ' -d %s.mmi %s\n' % (genome_noHM, genome_noHM)
        cmds += 'rm -rf %s.fai\n' % (genome_noHM)
        cmds += 'source deactivate py37\n'
        cmds += os.path.join(os.path.split('args.bw')[0], 'bowtie2-build') + ' %s %s\n' % (
            genome_noHM, genome_noHM)
    return [genome_noHM,cmds]

def shortname(genome,dna = False):
    newoutput = []
    donor_species = os.path.split(genome)[-1].split('.all.spades')[0]
    donor_species_new = donor_species.replace('_clustercluster', '_CL').replace('_PB_', '_PaDi_')
    change_name = set()
    for record in SeqIO.parse(genome, 'fasta'):
        record_id = str(record.id)
        if dna:
            record_id = 'C_%s' % (record_id.split('_')[1])
        else:
            record_id = 'C_%s_G_%s' % (record_id.split('_')[1], record_id.split('_')[-1])
        if len('%s__%s'%(donor_species_new, record_id)) > 32:
            donor_species_new = donor_species_new.replace('_IBD','')
        newoutput.append('>%s__%s\n%s\n' % (donor_species_new, record_id, str(record.seq)))
        temp_line = '%s\t%s\t%s__%s\t\n' % (donor_species, record_id,donor_species_new, record_id)
        change_name.add(temp_line)
    f1 = open(genome + '.prokka.fasta', 'w')
    f1.write(''.join(newoutput))
    f1.close()
    f1 = open(genome + '.changename.txt', 'w')
    f1.write(''.join(list(change_name)))
    f1.close()


def runprokka(genome):
    cmds = ''
    if 'prokka_' not in genome:
        genomefolder,filename = os.path.split(genome)
        output_dir = '%s/prokka_%s'%(genomefolder,filename)
        if len(glob.glob('%s/PROKKA_*.tsv'%(output_dir))) == 0:
            try:
                f1 = open(genome + '.fna', 'r')
            except FileNotFoundError:
                os.system('%s -q -i %s -d %s.fna -a %s.faa' % (args.pro, genome, genome,genome))
            shortname(genome,True)
            cmds = '%s --kingdom Bacteria --force --outdir %s --locustag Bacter %s\n' % \
                         (args.prokka, output_dir,
                          genome + '.prokka.fasta')
    return cmds



def mapping_WGS(donor_species,donor_species_fastqall):
    filesize = 0
    try:
        filesize = int(os.path.getsize(os.path.join(output_dir + '/merge/', donor_species + '.all.flt.snp.vcf')))
    except FileNotFoundError:
        pass
    if filesize == 0:
        print('generate mapping code for %s' % (donor_species))
        folder = output_dir + '/co-assembly/%s' % (donor_species)
        try:
            os.mkdir(folder)
        except IOError:
            pass
        cmds = ''
        # check assmebly size
        allgenomesize = []
        species = donor_species.split('_')[0]
        for fastq_file in donor_species_fastqall:
            if '.all' + fastq_name not in fastq_file and '.all.spades' not in fastq_file:
                original_folder, fastq_file_name = os.path.split(fastq_file)
                fastq_file_name = os.path.split(fastq_file)[-1]
                filename = fastq_file.split(fastq_name)[0]
                fastq_file2 = filename + fastq_name.replace('1', '2')
                if fastq_file_name == 'filter_reads_1.fastq':
                    # IBD
                    genome_file = glob.glob('/scratch/users/mit_alm/IBD_Evo/%s/Assembly_for_gene_flow/%s/contigs.fasta'%(species,
                                     original_folder.split('/')[Folder_num])) #+\
                                  #glob.glob(os.path.join(args.i + '/*/', '%s%s' % (original_folder.split('/')[Folder_num], genome_name)))
                    #if 'PB' in original_folder:
                    #    genome_file = glob.glob(os.path.join(args.i + '/*/', '%s%s' % (
                    #        '%s_PB_%s' % (original_folder.split('/')[Folder_num].split('_')[0],
                    #                      original_folder.split('/')[Folder_num].split('_')[1]), genome_name)))
                else:
                    genome_file = glob.glob(
                        os.path.join(original_folder, fastq_file_name.split(fastq_name)[0] + genome_name)) + \
                                  glob.glob(os.path.join(original_folder + '/../',
                                                         fastq_file_name.split(fastq_name)[0] + genome_name))+ \
                                  glob.glob(os.path.join(original_folder + '/../fasta/',
                                                         fastq_file_name.split(fastq_name)[0] + '_final.scaffolds.fasta'))
                if genome_file == []:
                    # need assemble the genome first
                    if fastq_file_name == 'filter_reads_1.fastq':
                        # IBD
                        genome_file = os.path.join(args.i + '/%s' % (donor_species),
                                                   '%s%s' % (original_folder.split('/')[Folder_num], genome_name))
                        if 'PB' in original_folder:
                            genome_file = os.path.join(args.i + '/%s' % (donor_species), '%s%s' % (
                                '%s_PB_%s' % (original_folder.split('/')[Folder_num].split('_')[0],
                                              original_folder.split('/')[Folder_num].split('_')[1]), genome_name))
                    else:
                        genome_file = os.path.join(original_folder, fastq_file_name.split(fastq_name)[0] + genome_name)
                    print('assemblying genome %s' % (genome_file))
                    #cmds = runspades_single(fastq_file, fastq_file2, genome_file)
                    #os.system(cmds)
                    #cmds = '' # NEED CHANGE SKIP SPADES
                    #filesize = int(os.path.getsize(genome_file))
                    filesize = 1
                else:
                    genome_file = genome_file[0]
                    filesize = int(os.path.getsize(genome_file))
                allgenomesize.append(filesize)
            else:
                allgenomesize.append(0)
        # filter out genomesize that's above the size cutoff
        maxgenome = statistics.mean(allgenomesize) * genomesize_cutoff
        qualifygenome = [donor_species_fastqall[i] for i in range(0, len(allgenomesize)) if
                         allgenomesize[i] <= maxgenome and allgenomesize[i] > 0]
        if species not in Ref_set:
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
                try:
                    f1 = open(donor_species_fastq, 'r')
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
            donor_species_genomename = Ref_set[species]
        if 'noHM.fasta' not in donor_species_genomename:
            donor_species_genomename, cmds2 = remove_homologous(donor_species_genomename)
            cmds += cmds2
        #cmds += runprokka(donor_species_genomename)
        alloutputvcf = os.path.join(output_dir + '/merge', donor_species + '.all')
        allsam = []
        cmds += 'source deactivate py37\n'
        for fastq_file in donor_species_fastqall:
            if '.all' + fastq_name not in fastq_file and '.all.spades' not in fastq_file:
                fastq_file_name = os.path.split(fastq_file)[-1]
                if fastq_file_name == 'filter_reads_1.fastq':
                    # 'IBD'
                    original_folder, fastq_file_name = os.path.split(fastq_file)
                    fastq_file_name = original_folder.split('/')[Folder_num]
                    if 'PB' in original_folder:
                        fastq_file_name = '%s_PB_%s' % (fastq_file_name.split('_')[0], fastq_file_name.split('_')[1])
                filename = fastq_file.split(fastq_name)[0]
                fastq_file2 = filename + fastq_name.replace('1', '2')
                # map each WGS to pangenome
                if fastq_file == ref_fastq:
                    # set as reference genome
                    results = run_vcf_WGS(fastq_file, fastq_file2,
                                      donor_species_genomename,
                                      os.path.join(output_dir + '/bwa',
                                                   fastq_file_name + '.ref'))
                    outputvcf = os.path.join(output_dir + '/merge', '%s.%s' % (donor_species, fastq_file_name + '.ref'))
                else:
                    results = run_vcf_WGS(fastq_file, fastq_file2,
                                          donor_species_genomename,
                                          os.path.join(output_dir + '/bwa',
                                                       fastq_file_name))
                    outputvcf = os.path.join(output_dir + '/merge', '%s.%s' % (donor_species, fastq_file_name))
                    allsam += [results[1]]
                cmds += results[0]
                cmds += merge_sample(donor_species_genomename, outputvcf, [results[1]], True)
        print(allsam)
        cmds += merge_sample(donor_species_genomename, alloutputvcf, allsam, False)
        f1 = open(os.path.join(input_script_vcf, '%s.vcf.sh' % (donor_species)), 'w')
        # for blast
        f1.write(
            '#!/bin/bash\nsource ~/.bashrc\npy37\n%s' % (
                ''.join(cmds)))
        f1.close()

################################################## Main ########################################################
# load reference genomes if provided
Ref_set = dict()
if args.ref != 'None':
    for lines in open(args.ref,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        donor_species, ref = lines_set[0:2]
        Ref_set.setdefault(donor_species, ref)

# load clustering results
if args.cl != 'None':
    Cluster = dict()
    allclustering = glob.glob(args.cl + '/clonal_population/*.genome.cluster.txt')
    for clusterfile in allclustering:
        species = os.path.split(clusterfile)[-1].split('.genome.cluster.txt')[0]
        print('find all fastq for each cluster of %s' % (species))
        for lines in open(clusterfile,'r'):
            if not lines.startswith('species'):
                lines_set = lines.split('\n')[0].split('\t')
                genome, cluster, total_cluster = lines_set[1:4]
                Cluster.setdefault(species, dict())
                Cluster[species].setdefault(cluster, [])
                if genome.startswith('P') or genome.startswith('H') or genome.startswith('D'):
                    fastq_file = glob.glob('/scratch/users/mit_alm/IBD_Evo/BA/*/%s/sickle2050/filter_reads_1.fastq'%(genome))+ \
                                 glob.glob('/scratch/users/mit_alm/IBD_Evo/BL/*/%s/sickle2050/filter_reads_1.fastq' % (genome))+ \
                                 glob.glob('/scratch/users/mit_alm/IBD_Evo/PB/*/%s_%s/sickle2050/filter_reads_1.fastq' % (
                                     genome.split('_')[0],genome.split('_')[-1]))
                    #fastq_file = glob.glob('/scratch/users/mit_alm/gutevo/2016_09_20_Bfragilis_TS1/%s/sickle2050/filter_reads_1.fastq' % (genome))
                else:
                    fastq_file = glob.glob(os.path.join(args.i + '/*/fastq/', '%s%s' % (genome, fastq_name))) + \
                                 glob.glob(os.path.join(args.i + '/*/', '%s%s' % (genome, fastq_name)))
                if fastq_file!=[]:
                    Cluster[species][cluster].append(fastq_file[0])
    # generate code
    for species in Cluster:
        for cluster in Cluster[species]:
            donor_species = '%s_cluster%s'%(species,cluster)
            donor_species_fastqall = Cluster[species][cluster]
            if len(donor_species_fastqall) >= cluster_cutoff:
                try:
                    ref_fastq = Cluster[species]['cluster%s'%(int(cluster.split('cluster')[1]) + 1)][0]
                except KeyError:
                    ref_fastq = Cluster[species]['cluster1'][0]
                if ref_fastq not in donor_species_fastqall:
                    donor_species_fastqall.append(ref_fastq)
                else:
                    ref_fastq = ''
                print('map genomes for %s' % (donor_species))
                mapping_WGS(donor_species, donor_species_fastqall)
else:
    # generate code based on donor-species
    for folder in species_folder:
        donor_species = os.path.split(folder)[-1]
        donor_species_fastqall = glob.glob(os.path.join(folder, 'fastq/*' + fastq_name)) + \
                                 glob.glob(os.path.join(folder, '*' + fastq_name))
        ref_fastq = ''
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
