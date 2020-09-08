# start
# Jay's folder
# round 1 run genome assembly and map genomes to a reference genome
import glob
import os
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/SNP_curate'
input_script_vcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round1'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay'
old_genome_dir = glob.glob('/scratch/users/mit_alm/IBD_Evo/*')
genome_dir = '/scratch/users/anniz44/genomes/donor_species/jay/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round1'
fastq_name = '_1.fastq'
genome_name = '.fasta'

os.system('rm -rf %s'%(input_script_sub))
os.system('rm -rf %s'%(input_script_vcf))

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(input_script_vcf)
except IOError:
    pass

try:
    os.mkdir(genome_dir)
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
    os.mkdir(output_dir+'/bwa/0')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/merge')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/SNP_curate')
except IOError:
    pass

try:
    os.mkdir(input_script)
except IOError:
    pass

def subset(file1,file2,output_name1,output_name2):
    cmds = 'head -500000 %s >> %s\n' % (
        file1, output_name1)
    cmds += 'tail -500000 %s >> %s\n' % (
        file1, output_name1)
    cmds += 'head -500000 %s >> %s\n' % (
        file2, output_name2)
    cmds += 'tail -500000 %s >> %s\n' % (
        file2, output_name2)
    return cmds

def runspades(file1,file2,temp_output,output_name):
    cmds = 'spades.py --careful -1 %s -2 %s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' % \
            (file1, file2, temp_output)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_output, output_name)
    cmds += 'rm -rf %s %s %s\n' % (temp_output, file1, file2)
    return cmds

def run_vcf(genome_file,database,tempbamoutput):
    # generate code
    # for curated genome
    genome_file = genome_file + '.corrected.fasta'
    cmds = 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, genome_file, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    sample_set3 = ['%s.sorted.bam' % (tempbamoutput)]
    cmds += 'rm -rf %s.bam %s.bam.bai\n' %(tempbamoutput,tempbamoutput)
    return [cmds,sample_set3]

def run_vcf_curate(files,files2,database,tempbamoutput):
    # generate code
    cmds = ''
    Runcode = False
    try:
        f1 = open(tempbamoutput + '.flt.snp.vcf')
        filesize =int(os.path.getsize(tempbamoutput + '.flt.snp.vcf'))
        if filesize == 0:
            Runcode = True
    except IOError:
        Runcode = True
    if Runcode:
        cmds = 'rm -rf %s.fai\n' % (database)
        cmds += 'bowtie2-build %s %s\n' % (database, database)
        cmds += 'bowtie2' + ' -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
            min(40, 40), database, files, files2, 'samtools', min(40, 40),
            tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
            tempbamoutput)
        cmds += '#samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
            tempbamoutput, tempbamoutput)
        cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s.sorted.bam  | %s call --threads %s -m > %s.raw.vcf\n' % (
            'bcftools', min(40, 40), database,
            tempbamoutput, 'bcftools', min(40, 40), tempbamoutput)
        cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
            'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
        cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
            'bcftools', tempbamoutput, tempbamoutput)
        cmds += 'rm -rf %s.bam %s.sorted.bam %s.sorted.bam.bai\n' % (
            tempbamoutput,tempbamoutput,tempbamoutput)
        cmds += 'rm -rf %s.flt.vcf %s.raw.vcf\n' % (tempbamoutput,tempbamoutput)
    else:
        os.system('rm -rf %s.fai %s.*.bt2' % (database,database))
    return cmds

def merge_sample(database,tempbamoutputGenome,samplesetGenomecorrect):
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
cmds_cp = ''
for folder in old_genome_dir:
    species = os.path.split(folder)[-1]
    if len(species) == 2 and species in ['BA','BL','PB']:
        alldonors = glob.glob('%s/*' % (folder))
        for donorfolder in alldonors:
            donor = os.path.split(donorfolder)[-1].split('_')[-1]
            if donor!= 'TOT':
                new_folder = '%s/%s_%s' %(genome_dir,donor,species)
                # generate codes
                donor_species = '%s_%s' %(donor,species)
                allfastq = glob.glob('%s/*/sickle2050/*%s' % (donorfolder, fastq_name))
                if allfastq != []:
                    cmds_sub = ''
                    cmds = ''
                    cmds2 = ''
                    i = 0
                    split_bash = 1
                    # set up spades
                    donor_species_folder_all = os.path.join(new_folder, donor_species + '_allspades1')
                    donor_species_genomename = os.path.join(new_folder, donor_species + '.all.spades1.fasta')
                    donor_species_fastq = os.path.join(new_folder, donor_species + '.all.spades1' + fastq_name)
                    donor_species_fastq2 = os.path.join(new_folder,
                                                        donor_species + '.all.spades1' + fastq_name.replace('1', '2'))
                    tempbamoutputGenome = os.path.join(output_dir + '/merge', donor_species + '.all')
                    samplesetGenomecorrect = []
                    try:
                        os.mkdir(new_folder)
                    except IOError:
                        pass
                    for fastq_file in allfastq:
                        if '.all' + fastq_name not in fastq_file and '.all.spades1' + fastq_name not in fastq_file:
                            genomename = '%s_%s_%s%s' % (donor, species,
                                                         fastq_file.split('/sickle2050')[0].split('_')[-1], genome_name)
                            genome = glob.glob('%s/Assembly_for_gene_flow/%s/scaffolds.fasta' % (folder,
                                                                                                 fastq_file.split(
                                                                                                     '/sickle2050')[0].split(
                                                                                                     '/')[-1]))
                            genome_file = '%s/%s' %(new_folder,genomename)
                            fastq_file2 = fastq_file.replace(fastq_name,
                                                             fastq_name.replace('1', '2'))
                            fastq_file_name = os.path.split(genome_file)[-1]
                            if genome!= []:
                                # copy genome to genome_folder
                                genome = genome[0]
                                cmds_cp += ('cp %s %s/%s\n' % (genome, new_folder, genomename))
                                # subset fastq
                                cmds += subset(fastq_file, fastq_file2, donor_species_fastq, donor_species_fastq2)
                                # mapping assembly to simulated WGS to curate genome
                                cmds_sub += run_vcf_curate(fastq_file, fastq_file2, genome_file,
                                                           os.path.join(output_dir + '/SNP_curate',
                                                                        fastq_file_name + '.curate'))
                                # run each mapping corrected genome to pangenome
                                results = run_vcf(genome_file, donor_species_genomename,
                                                  os.path.join(output_dir + '/bwa/0',
                                                               fastq_file_name)
                                                  )
                                cmds2 += results[0]
                                samplesetGenomecorrect += results[1]
                                i += 1
                                if i % 10 == 0 and cmds_sub != '':
                                    f1 = open(
                                        os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species, int(i / 10))),
                                        'a')
                                    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds_sub)))
                                    f1.close()
                                    cmds_sub = ''
                            else:
                                # run spades
                                #cmds_sub += runspades(fastq_file, fastq_file2, genome_file.split(genome_name)[0],
                                #                      genome_file)
                                print('missing genome file for %s'%(fastq_file))
                    cmds += runspades(donor_species_fastq, donor_species_fastq2, donor_species_folder_all,
                                      donor_species_genomename) # pan genome
                    # then run mapping corrected genome to pangenome
                    cmds += 'minimap2 -d %s.mmi %s\n' % (donor_species_genomename, donor_species_genomename)
                    cmds += 'rm -rf %s.fai\n' % (donor_species_genomename)
                    cmds += 'bowtie2-build %s %s\n' % (donor_species_genomename, donor_species_genomename)
                    cmds += cmds2
                    cmds += merge_sample(donor_species_genomename, tempbamoutputGenome, samplesetGenomecorrect)
                    if cmds_sub != '':
                        f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species, int(i / 10))), 'a')
                        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds_sub)))
                        f1.close()
                    f1 = open(os.path.join(input_script_vcf, '%s.vcf.sh' % (donor_species)), 'w')
                    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                    f1.close()
                else:
                    print('no fastq in %s'%('%s/*/sickle2050/*%s' % (donorfolder, fastq_name)))

#f1 = open(os.path.join(input_script, 'cpgenome.sh'), 'w')
#f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds_cp)))
#f1.close()

f1 = open(os.path.join(input_script, 'allcurate.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'allround1.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
