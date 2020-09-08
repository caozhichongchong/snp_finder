# start
# Jay's folder
# round 1-3 run genome assembly and map genomes to a reference genome
import glob
import os

Round = 4
input_script_vcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round%s'%(Round)
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay'
old_genome_dir = '/scratch/users/mit_alm/IBD_Evo/'
genome_root = '/scratch/users/anniz44/genomes/donor_species/jay'
genome_dir = glob.glob('%s/round%s/*'%(genome_root,Round))
if Round == 4:
    genome_dir = glob.glob('%s/round*/*' % (genome_root))

output_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s'%(Round)
output_dir_old = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s'%(Round-1)
fastq_name = '_1.fastq'
genome_name = '.fasta.corrected.fasta'

try:
    os.mkdir(input_script)
except IOError:
    pass

os.system('rm -rf %s' % (input_script_vcf))

try:
    os.mkdir(input_script_vcf)
except IOError:
    pass

try:
    os.mkdir(output_dir)
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
    os.mkdir(output_dir + '/merge')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/merge_genome')
except IOError:
    pass

def subset(file1, file2, output_name1, output_name2):
    cmds = 'head -500000 %s >> %s\n' % (
        file1, output_name1)
    cmds += 'tail -500000 %s >> %s\n' % (
        file1, output_name1)
    cmds += 'head -500000 %s >> %s\n' % (
        file2, output_name2)
    cmds += 'tail -500000 %s >> %s\n' % (
        file2, output_name2)
    return cmds

def runspades(file1, file2, temp_output, output_name):
    cmds = 'spades.py --careful -1 %s -2 %s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' % \
           (file1, file2, temp_output)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_output, output_name)
    cmds += 'rm -rf %s %s %s\n' % (temp_output,file1, file2)
    return cmds

def run_vcf(genome_file, database, tempbamoutput):
    # generate code
    # for curated genome
    cmds = 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, genome_file, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    sample_set3 = ['%s.sorted.bam' % (tempbamoutput)]
    cmds += 'rm -rf %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    return [cmds, sample_set3]

def run_vcf_WGS(files,files2,database,tempbamoutput):
    # generate code
    cmds = 'rm -rf %s.fai\n' % (database)
    cmds += 'bowtie2' + ' -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, files, files2, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    cmds += 'rm -rf %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    sample_set3 = ['%s.sorted.bam' % (tempbamoutput)]
    return [cmds, sample_set3]

def merge_sample(database, tempbamoutputGenome, samplesetGenomecorrect):
    # corrected genomes
    cmds = '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
        'bcftools', min(40, 40), database,
        ' '.join(samplesetGenomecorrect), 'bcftools', min(40, 40), tempbamoutputGenome)
    cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
        'bcftools', min(40, 40), tempbamoutputGenome, tempbamoutputGenome)
    cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
        'bcftools', tempbamoutputGenome, tempbamoutputGenome)
    cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutputGenome)
    return cmds

# generate code

for new_folder in genome_dir:
    donor, species = os.path.split(new_folder)[-1].split('_cluster')[0].split('_')
    donor_species = os.path.split(new_folder)[-1]
    allgenome = glob.glob('%s/*%s'%(new_folder,genome_name))
    cmds = ''
    cmds2 = ''
    # set up spades
    donor_species_folder_all = os.path.join(new_folder, donor_species + '_allspades%s' % (Round))
    donor_species_genomename = os.path.join(new_folder, donor_species + '.all.spades%s.fasta' % (Round))
    donor_species_fastq = os.path.join(new_folder, donor_species + '.all.spades%s' % (Round) + fastq_name)
    donor_species_fastq2 = os.path.join(new_folder,
                                        donor_species + '.all.spades%s' % (Round) + fastq_name.replace(
                                            '1', '2'))
    tempbamoutputGenome = os.path.join(output_dir + '/merge', donor_species + '.all')
    samplesetGenomecorrect = []
    if len(allgenome) >= 3:
        if Round == 4:
            donor_species_genomename = glob.glob(os.path.join(new_folder, donor_species + '.all.spades*.fasta'))
            if donor_species_genomename == []:
                donor_species_old = '_cluster'.join(donor_species.split('_cluster')[0:-1])
                donor_species_genomename = glob.glob(os.path.join(genome_root + '/round*/%s' % (donor_species_old),
                                                                  donor_species_old + '.all.spades*.fasta'))[0]
                cmds += 'cp %s/vcf_round*/merge/%s.*filtered* %s/merge_genome/\n' % (genome_root, donor_species_old,
                                                                              output_dir)
            else:
                donor_species_genomename = donor_species_genomename[0]
                cmds += 'cp %s/vcf_round*/merge/%s.*filtered* %s/merge_genome/\n' % (genome_root, donor_species,
                                                                              output_dir)
            print(donor_species_genomename)
            tempbamoutputGenome = tempbamoutputGenome + '.fq'
            print(new_folder,tempbamoutputGenome)
        for genome_file in allgenome:
            fastq_file_name = os.path.split(genome_file)[-1].split(genome_name)[0]
            genomeID = fastq_file_name.split('_')[-1]
            fastq_file = glob.glob('%s/%s/*%s/%s*_%s/sickle2050/*%s'%(old_genome_dir,species,donor,
                                                                      donor,genomeID,fastq_name))
            if fastq_file!= []:
                fastq_file = fastq_file[0]
                fastq_file2 = fastq_file.replace(fastq_name,
                                                 fastq_name.replace('1', '2'))
                if Round <= 3:
                    # subset fastq
                    cmds += subset(fastq_file, fastq_file2, donor_species_fastq, donor_species_fastq2)
                    # run each mapping corrected genome to pangenome
                    results = run_vcf(genome_file, donor_species_genomename,
                                      os.path.join(output_dir + '/bwa/0',
                                                   fastq_file_name)
                                      )
                else:
                    # WGS for confirmation
                    results = run_vcf_WGS(fastq_file, fastq_file2,
                                donor_species_genomename,
                                          os.path.join(output_dir + '/bwa/0',
                                                   fastq_file_name))
                cmds2 += results[0]
                samplesetGenomecorrect += results[1]
            else:
                print('no fastq in %s' % (
                        '%s/%s/*%s/%s/sickle2050/*%s' % (old_genome_dir, species, donor,
                                                         fastq_file_name, fastq_name)))
        if Round <= 3:
            # pan-genome
            cmds += runspades(donor_species_fastq, donor_species_fastq2, donor_species_folder_all,
                              donor_species_genomename)  # pan genome
            # then run mapping corrected genome to pangenome
            cmds += 'minimap2 -d %s.mmi %s\n' % (donor_species_genomename, donor_species_genomename)
            cmds += 'rm -rf %s.fai\n' % (donor_species_genomename)
            cmds += 'bowtie2-build %s %s\n' % (donor_species_genomename, donor_species_genomename)
        cmds += cmds2
        cmds += merge_sample(donor_species_genomename, tempbamoutputGenome, samplesetGenomecorrect)
    f1 = open(os.path.join(input_script_vcf, '%s.vcf.sh' % (donor_species)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
    f1.close()

f1 = open(os.path.join(input_script, 'allround%s.sh'%(Round)), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.vcf.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    donor_species_old = '_cluster'.join(sub_scripts_name.split('_cluster')[0:-1])
    if Round <= 3:
        if len(glob.glob(os.path.join(input_script_vcf,donor_species_old+'*'))) > 1:
            # more than 1 subcluster
            f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
        else:
            # only 1 cluster, copy previous round results
            f1.write('cp %s/merge/%s.*filtered* %s/merge/\n' % (output_dir_old, donor_species_old,
                                                       output_dir))
    else:
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
