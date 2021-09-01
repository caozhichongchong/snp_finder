################################################### SET PATH ########################################################
import glob
import os
# select genomes of same person of same species
metadata='/scratch/users/anniz44/scripts/1MG/metadata/GTDB_taxon_CG_GMC.brief.habitat.species.selected.all'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/roary2'
input_script_subgubbbins = '/scratch/users/anniz44/scripts/1MG/donor_species/gubbins'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
genome_dir = '/scratch/users/anniz44/genomes/GHM/newgenomes/*_Genome_Library_*'
output_dir = '/scratch/users/anniz44/genomes/donor_species'
#os.system('rm -rf %s' %(output_dir))
#os.system('rm -rf %s' %('.'))
try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(input_script_subgubbbins)
except IOError:
    pass

try:
    os.mkdir(output_dir + '/gubbins')
except IOError:
    pass

Set1 = dict()
Set2 = dict()
Species = dict()
Species_correct = dict()
for lines in open(metadata,'r'):
    if not lines.startswith('GCA_'):
        genome = lines.split('\t')[0]
        species = lines.split('\t')[1].replace(' ', '_')
        donor = genome.split('_')[0]
        #try:
        #    if species.split('_')[1].startswith('sp'):
        #        print(species, species.split('_')[0])
        #        species = species.split('_')[0]
        #except IndexError:
        #    pass
        donor_species = '%s_%s'%(donor,species)
        Set1.setdefault(donor_species, [])
        Set1[donor_species].append(genome)
        Set2.setdefault(donor_species, 0)
        Set2[donor_species]+=1

Set2_list=[]
cmdscopy = ''
cmdsgubbins = ''
i=0
donor_list =["cx","bk","av","ao","an","am","af","aa",
             "8306RM","7652PC","6937AY","6587MW","6383DG",
             "4117FE","3790QQ","3554DX","3321JW","3219KR"]
for donor_species in Set2:
    donor = donor_species.split('_')[0]
    print(donor,donor in donor_list)
    print(Set2[donor_species] )
    #if Set2[donor_species] >=50 and len(donor_species.split('_'))>2:
    if donor in donor_list and Set2[donor_species] >=10 and Set2[donor_species] < 50:
        new_line = '%s\t%s\t'%(donor_species,Set2[donor_species])
        try:
            output_subdir = os.path.join(output_dir, donor_species)
            cmdsroary = ''
            cmdsroary += 'roary -p 40 -o roary_%s -e -n --mafft -f %s/roary_%s  -i 80 -cd 95 %s/*.gff\n' % (
            donor_species, output_dir,donor_species, output_subdir)
            cmdsgubbins += 'run_gubbins.py --threads 40  --tree_builder hybrid --use_time_stamp --prefix %s/gubbins/gubbins_%s_core_gene_alignment --verbose %s/roary_%s/core_gene_alignment.aln\n' % (
                output_dir, donor_species, output_dir, donor_species)
            f1 = open(os.path.join(input_script_sub, '%s.roary.sh'%(i)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmdsroary)))
            f1.close()
            i+=1
            os.mkdir(output_subdir)
        except OSError:
            pass
        for genomes in Set1[donor_species]:
            new_line += '%s\t'%(genomes)
            cmdscopy += 'cp %s/%s %s\n' % (
            genome_dir, genomes + '_*.fasta', os.path.join(output_dir, donor_species + '/' + genomes + '.fasta'))
            cmdscopy += 'cp %s/%s %s\n' % (
                genome_dir, genomes + '_*.gff', os.path.join(output_dir, donor_species + '/' + genomes + '.gff'))
            cmdscopy += 'cp %s/%s %s\n' % (
                genome_dir, genomes + '_*.faa', os.path.join(output_dir, donor_species + '/' + genomes + '.faa'))
        new_line += '\n'
        Set2_list.append(new_line)

f1 = open(os.path.join(input_script,'species.donor.selected.txt'),'a')
f1.write(''.join(Set2_list))
f1.close()

f1 = open(os.path.join(input_script_subgubbbins, 'gubbins.sh'), 'a')
f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s' % (''.join(cmdsgubbins)))
f1.close()

f1 = open(os.path.join(input_script,'copy2.sh'),'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(''.join(cmdscopy)))
f1.close()

f1 = open(os.path.join(input_script, 'allroary2.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.roary.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

# put whole genome fastq into the right folder
import glob
import os
input_fastq = '/scratch/users/anniz44/genomes/GHM/newgenomes/BN10_WG'
input_genome = '/scratch/users/anniz44/genomes/donor_species/species'
species = []
for lines in open(os.path.join(input_fastq,'table_clusternames_taxo.sort.txt'),'r'):
    lines_set = lines.replace('\n','').replace('\r','').split('\t')
    cluster = lines_set[0]
    a_species = lines_set[1]
    os.system('mv %s %s'%(os.path.join(input_fastq,cluster),
                          os.path.join(input_fastq, a_species)))
    species.append(a_species)

for folders in glob.glob(os.path.join(input_genome,'*_cluster*')):
    donor_species = os.path.split(folders)[-1].split('_cluster')[0]
    print(donor_species)
    if glob.glob(os.path.join(input_genome, donor_species)) != []:
        os.system('mv %s/* %s' % (folders, os.path.join(input_genome, donor_species)))
        os.system('rm -rf %s' % (folders))
    else:
        print(os.path.join(input_genome, donor_species))

for a_species in species:
    input_fastq_subfolder = os.path.join(input_fastq,a_species)
    allfastq = glob.glob(input_fastq_subfolder + '/*.fastq')
    for a_fastq in allfastq:
        donor = os.path.split(a_fastq)[-1].split('_')[0]
        donor_species = '%s_%s' %(donor,a_species)
        if glob.glob(os.path.join(input_genome,donor_species + '/*.*')) !=[]:
            os.system('mv %s %s/' % (a_fastq,os.path.join(input_genome,donor_species)))
        else:
            print(os.path.join(input_genome,donor_species),input_fastq_subfolder)

################################################### END ########################################################
################################################### SET PATH ########################################################
# rename genome and whole genome fastq
# put whole genome fastq into the right folder
import glob
import os
input_genome = '/scratch/users/anniz44/genomes/donor_species/species'
all_folder = glob.glob(os.path.join(input_genome,'*'))
file_format_list = ['.faa','.gff','R2.concat.trim.fastq','.fasta']
Change_file_list = []
for donor_species_folder in all_folder:
    donor_species = os.path.split(donor_species_folder)[-1]
    donor=donor_species.split('_')[0]
    if len(donor) == 2: #BN10 data
        allfasta =  glob.glob(os.path.join(donor_species_folder,'*R1.concat.trim.fastq'))
        i = 0
        for a_fasta in allfasta:
            a_fasta_name = os.path.split(a_fasta)[-1].replace('R1.concat.trim.fastq','')
            i+=1
            Change_file_list.append('%s\t%s_g%s\t\n'%(a_fasta_name,donor_species,i))
            cmd = ('mv %s %s//%s_g%s%s' % (a_fasta,donor_species_folder,donor_species,i,'R1.concat.trim.fastq'))
            for file_format in file_format_list:
                cmd += (' && mv %s %s//%s_g%s%s' % (a_fasta.replace('R1.concat.trim.fastq',file_format), donor_species_folder,donor_species, i, file_format))
            os.system(cmd)

Change_file_list_file = open(os.path.join(input_genome,'file.change.name.txt'),'a')
Change_file_list_file.write(''.join(Change_file_list))
Change_file_list_file.close()

# remove wrong genomes
import glob
import os
input_genome = '/scratch/users/anniz44/genomes/donor_species/species'
not_used = '/scratch/users/anniz44/genomes/donor_species/species/no_fastq'
all_folder = glob.glob(os.path.join(input_genome,'*'))

try:
    os.mkdir(not_used)
except IOError:
    pass

Change_file_list_file=dict()
for lines in open(os.path.join(input_genome,'file.change.name.txt'),'r'):
    Change_file_list_file.setdefault(lines.split('\t')[1],lines.split('\t')[0])

file_format_list = ['.faa','.gff','R2.concat.trim.fastq','.fasta']
for donor_species_folder in all_folder:
    donor_species = os.path.split(donor_species_folder)[-1]
    donor = donor_species.split('_')[0]
    if len(donor) == 2:  # BN10 data
        allR1 = glob.glob(os.path.join(donor_species_folder, '*R1.concat.trim.fastq'))
        for a_fasta in allR1:
            try:
                f1 = open(a_fasta.replace('R1.concat.trim.fastq','R2.concat.trim.fastq'),'r')
            except IOError:
                print(a_fasta.replace('R1.concat.trim.fastq','R2.concat.trim.fastq'),'not found')
                a_fasta_folder, a_fasta_name = os.path.split(a_fasta)
                a_fasta_name = a_fasta_name.split('R1.concat.trim.fastq')[0]
                if a_fasta_name in Change_file_list_file:
                    for file_format in file_format_list:
                        os.system('mv %s/%s%s %s/%s%s' % (a_fasta_folder,Change_file_list_file[a_fasta_name],file_format,
                                                      a_fasta_folder,a_fasta_name,file_format))
                    try:
                        f1 = open(a_fasta.replace('R1.concat.trim.fastq', 'R2.concat.trim.fastq'), 'r')
                    except IOError:
                        print(a_fasta.replace('R1.concat.trim.fastq', 'R2.concat.trim.fastq'), 'not found')
                        os.system('mv %s %s/' % (a_fasta, not_used))
        allfasta = glob.glob(os.path.join(donor_species_folder, '*.fasta'))
        for a_fasta in allfasta:
            try:
                f1 = open(a_fasta.replace('.fasta','R1.concat.trim.fastq'),'r')
            except IOError:
                print(a_fasta)
                os.system('mv %s %s/'%(a_fasta.replace('.fasta','.*'),not_used))
                os.system('mv %s %s/' % (a_fasta.replace('.fasta', 'R2.concat.trim.fastq'), not_used))

# add names to fastq without fasta
import glob
import os
Change_file_list = []
input_genome = '/scratch/users/anniz44/genomes/donor_species/species'
all_folder = glob.glob(os.path.join(input_genome,'*'))
file_format_list = ['.faa','.gff','R2.concat.trim.fastq','.fasta']
for donor_species_folder in all_folder:
    donor_species = os.path.split(donor_species_folder)[-1]
    fastq = glob.glob(os.path.join(donor_species_folder,'*R1.concat.trim.fastq'))
    named = []
    not_named = []
    if fastq != []:
        for fastq_file in fastq:
            fastq_name = os.path.split(fastq_file)[-1]
            if donor_species in fastq_name:
                named.append(int(fastq_name.split('R1.concat.trim.fastq')[0].split(donor_species+'_g')[1]))
            else:
                not_named.append(fastq_file)
        try:
            i = max(named)
        except ValueError:
            i = 0
        print(i,named)
        print(not_named)
        if not_named !=[]:
            for fastq_file in not_named:
                i += 1
                fastq_name = os.path.split(fastq_file)[-1].split('R1.concat.trim.fastq')[0]
                Change_file_list.append('%s\t%s_g%04d\t\n' % (fastq_name, donor_species, i))
                cmd = ('mv %s %s//%s_g%04d%s\n' % (fastq_file, donor_species_folder, donor_species, i, 'R1.concat.trim.fastq'))
                for file_format in file_format_list:
                    cmd += ('mv %s %s//%s_g%04d%s\n' % (fastq_file.replace('R1.concat.trim.fastq',
                                                                       file_format),
                                                    donor_species_folder, donor_species, i, file_format))
                os.system(cmd)

Change_file_list_file = open(os.path.join(input_genome,'file.change.name.txt'),'a')
Change_file_list_file.write(''.join(Change_file_list))
Change_file_list_file.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# Jay's folder
# round 1 run genome assembly and map genomes to a reference genome
import glob
import os
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/SNP_curate'
input_script_vcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay'
old_genome_dir = glob.glob('/scratch/users/mit_alm/IBD_Evo/*')
genome_dir = '/scratch/users/anniz44/genomes/donor_species/jay/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round2'
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
        cmds += 'bowtie2' + ' --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
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
################################################### SET PATH ########################################################
# round 1 run genome assembly and map genomes to a reference genome
import glob
import os
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNP_curate'
input_script_vcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_round2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/af_Bifidobacterium_longum')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round2'
fastq_name = '_1.fastq'
genome_name = '.fasta'

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

os.system('#rm -rf %s'%(input_script_sub))
os.system('$rm -rf %s'%(input_script_vcf))

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(input_script_vcf)
except IOError:
    pass

# put clonal populations back -> done
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    if len(donor) == 2 and 'cluster' in donor_species: # BN10 donors
        origina_folder = folder.split('_newcluster')[0].split('_cluster')[0]
        if glob.glob(origina_folder)!=[]:
            os.system('mv %s/* %s'%(folder,origina_folder))
            os.system('mv %s/fastq/* %s' % (folder, origina_folder))
        else:
            print(origina_folder, ' not found')

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
    cmds += '#samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
        tempbamoutput, tempbamoutput)
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
        cmds += 'bowtie2' + ' --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
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
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    if len(donor) == 2:
        # BN10 donors
        donor_species_genome = glob.glob(os.path.join(folder,'*.fasta'))
        sub_samples = []
        donor_species_fastqall = glob.glob(os.path.join(folder, 'fastq/*' + fastq_name)) +\
                              glob.glob(os.path.join(folder, '*' + fastq_name))
        if donor_species_fastqall != []:
            cmds_sub = ''
            cmds = ''
            cmds2 = ''
            i = 0
            split_bash = 1
            donor_species_folder_all = os.path.join(folder, donor_species + '_allspades1')
            donor_species_genomename = os.path.join(folder, donor_species + '.all.spades1.fasta')
            donor_species_fastq = os.path.join(folder, donor_species + '.all.spades1' + fastq_name)
            donor_species_fastq2 = os.path.join(folder, donor_species + '.all.spades1' + fastq_name.replace('1', '2'))
            tempbamoutputGenome = os.path.join(output_dir + '/merge', donor_species + '.all')
            samplesetGenomecorrect = []
            for fastq_file in donor_species_fastqall:
                if '.all' + fastq_name not in fastq_file and '.all.spades1' + fastq_name not in fastq_file:
                    fastq_file_name = os.path.split(fastq_file)[-1]
                    filename = fastq_file.split(fastq_name)[0]
                    fastq_file2 = filename + fastq_name.replace('1', '2')
                    genome_file = os.path.join(folder, fastq_file_name.split(fastq_name)[0] + genome_name)
                    if glob.glob(genome_file) == []:
                        # no assembly
                        cmds_sub += runspades(fastq_file, fastq_file2, filename + '_spades', genome_file)
                        print('run spades for ', genome_file)
                    # subset fastq
                    cmds +=subset(fastq_file,fastq_file2,donor_species_fastq,donor_species_fastq2)
                    # mapping assembly to simulated WGS to curate genome
                    cmds_sub += run_vcf_curate(fastq_file, fastq_file2, genome_file,
                                               os.path.join(output_dir+'/SNP_curate',
                                                            fastq_file_name+'.curate'))
                    # run each mapping corrected genome to pangenome
                    results = run_vcf(genome_file, donor_species_genomename,
                                      os.path.join(output_dir + '/bwa/0',
                                                   fastq_file_name)
                                      )
                    cmds2 += results[0]
                    samplesetGenomecorrect += results[1]
                    i += 1
                    if i % 10 == 0 and cmds_sub != '':
                        f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species, int(i / 10))), 'a')
                        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds_sub)))
                        f1.close()
                        cmds_sub = ''
            # pan genome
            cmds += runspades(donor_species_fastq, donor_species_fastq2, donor_species_folder_all, donor_species_genomename)
            # then run mapping corrected genome to pangenome
            cmds += 'minimap2 -d %s.mmi %s\n' % (donor_species_genomename, donor_species_genomename)
            cmds += 'rm -rf %s.fai\n' % (donor_species_genomename)
            cmds += 'bowtie2-build %s %s\n' % (donor_species_genomename, donor_species_genomename)
            cmds += cmds2
            cmds += merge_sample(donor_species_genomename, tempbamoutputGenome, samplesetGenomecorrect)
            if cmds_sub != '':
                f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species,int(i/10))), 'a')
                f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds_sub)))
                f1.close()
            f1 = open(os.path.join(input_script_vcf, '%s.vcf.sh' % (donor_species)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            f1.close()

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
################################################### SET PATH ########################################################
# run allcurate
# step 2 check SNPs
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
# set up path
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round2/SNP_curate'
fastq_name = '.curate.flt.snp.vcf'

input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay'
genome_root = '/scratch/users/anniz44/genomes/donor_species/jay/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/')
output_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round2/SNP_curate'
fastq_name = '.curate.flt.snp.vcf'
# set up cutoff
# good mapping
Major_alt_freq_cutoff = 0.7 # major alt freq in a sample
Sample_depth_cutoff = 5 # both forward and reverse reads cutoff in a sample

# good coverage
total_coverage_cutoff = 0.8 # at least X reads map to its original genome
genome_avg_coverage_cutoff = 10 # genome average coverage cutoff
# reasonable curation
Major_alt_freq_cutoff2 = 0.8 # major alt freq in a sample
Sample_depth_cutoff2 = Sample_depth_cutoff # both forward and reverse reads cutoff in a sample

Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def coverage_check(sh_file,Coverage):
    Coverage1 = []
    Coverage2 = []
    os.system('rm -rf %s/temp.sh %s/temp.err'%(input_script_split_sub,input_script_split_sub))
    os.system('grep \"bowtie2\" %s > %s/temp.sh' % (sh_file, input_script_split_sub))
    os.system('grep \"overall alignment rate\" %s.err > %s/temp.err' % (sh_file, input_script_split_sub))
    # load fastq name
    for lines in open('%s/temp.sh' % (input_script_split_sub), 'r'):
        try:
            filename = lines.split('>')[1].split('\n')[0].split(fastq_name + '.bam')[0]
            filename = os.path.split(filename)[-1]
            Coverage1.append(filename)
        except IndexError:
            pass
    # load reads alignment rate
    for lines in open('%s/temp.err' % (input_script_split_sub), 'r'):
        Coverage2.append(float(lines.split('%')[0]))
    i = 0
    for filename in Coverage1:
        # load average coverage of ref genome
        tempbamoutput = os.path.join(output_dir, filename + fastq_name)
        coverage_file = glob.glob(tempbamoutput + '.sorted.bam.avgcov')
        if coverage_file == []:
            cmds = 'samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
                tempbamoutput, tempbamoutput)
            os.system(cmds)
            coverage_file = tempbamoutput + '.sorted.bam.avgcov'
        else:
            coverage_file = coverage_file[0]
        length_total = 0
        coverage_total = 0
        donor_species = filename.split(genome_split)[0]
        Donor_species.setdefault(donor_species, [[], 0])
        try:
            for lines in open(coverage_file, 'r'):
                length_contig = float(lines.split('\t')[1])
                length_total += length_contig
                coverage_contig = float(lines.split('\t')[2])
                coverage_total += coverage_contig * length_contig
            coverage_num = Coverage2[i] / 100
            avg_coverage = coverage_total / length_total
            i += 1
        except IOError:
            coverage_num = total_coverage_cutoff - 0.1
            avg_coverage = genome_avg_coverage_cutoff - 1
        temp_line = ('%s\t%.2f\t%.1f\t%.1f' % (filename, coverage_num, length_total, avg_coverage))
        if coverage_num < total_coverage_cutoff or avg_coverage < genome_avg_coverage_cutoff:
            # not qualified
            print(filename, coverage_num, avg_coverage)
            Donor_species[donor_species][0].append(filename)
            temp_line += ('\tnot_qualified\t\n')
        else:
            # qualified
            Donor_species[donor_species][1] += 1
            temp_line += ('\tqualified\t\n')
        Coverage.append(temp_line)

def SNP_check(lines,donor_species,vcf_file_list):
    # CHR, POS, REF, ALT, good assembly, qualified mapping
    lines_set = lines.split('\n')[0].split('\t')
    report_line = ''
    temp_report_line = ['T','T']
    need_curation = 'F'
    REF = lines_set[3]
    allels_set = [REF]
    if '.' not in lines_set[4]:
        allels_set += lines_set[4].split(',')
    Total_alleles = len(allels_set)
    for Subdepth_all in lines_set[9:]:
            Allels_frq = [0, 0, 0, 0]
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
            total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
            Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
            Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
            total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
            total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
            for num_allels in range(0, Total_alleles):
                allels = allels_set[num_allels]
                Subdepth_alleles = int(Subdepth[num_allels])
                if allels in Allels:
                    Allels_frq[Allels[allels]] += Subdepth_alleles
                else:
                    pass
            # find major alt and calculate frequency
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            MLF = Major_ALT[1] / total_sub_depth
            if REF != Major_ALT[0]:
                # wrong assembly
                temp_report_line[0] = 'F' # F: not good assembly
                if total_sub_depth_forward >= Sample_depth_cutoff2 and \
                        total_sub_depth_reverse >= Sample_depth_cutoff2 and \
                        MLF >= Major_alt_freq_cutoff2:
                    # can be curated
                    need_curation = 'T'  # T: need curation
            if total_sub_depth_forward < Sample_depth_cutoff or \
                    total_sub_depth_reverse < Sample_depth_cutoff or \
                    MLF < Major_alt_freq_cutoff:
                # unqualified mapping
                temp_report_line[1] = 'F' # F: bad mapping
    if temp_report_line != ['T','T']:
        report_line += '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t\n' %(donor_species,lines_set[0],lines_set[1],
                                                    REF,','.join([ALT for ALT in allels_set if ALT != REF]),
                                                    temp_report_line[0],temp_report_line[1],need_curation,
                                                     Major_ALT[0], MLF,
                                                     total_sub_depth_forward, total_sub_depth_reverse)
        vcf_file_list.append(donor_species + '\t' + lines)
    return report_line

# check major alt
all_vcf_file=glob.glob(os.path.join(output_dir,'*%s'%(fastq_name)))
vcf_file_report = []
vcf_file_list = []
vcf_file_report.append('donor_species\tCHR\tPOS\tREF\tALT\tAssembly\tMapping_quality\tNeed_curation\tMajor_ALT\tMajor_ALT_frq\tDepth_F\tDepth_R\n')
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split(fastq_name)[0]
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            vcf_file_report.append(SNP_check(lines,donor_species,vcf_file_list))

f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.sum'), 'w')
f1.write(''.join(vcf_file_report))
f1.close()

f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.vcf'), 'w')
f1.write(''.join(vcf_file_list))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# step 3 curate genome
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
# set up path
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/*/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round2/SNP_curate'
fastq_name = '.curate.flt.snp.vcf'
samplename = '_1.fastq'
genome_name = '.fasta'
#species_select = 'am_Escherichia_coli'

input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay'
genome_root = '/scratch/users/anniz44/genomes/donor_species/jay/*/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/')
output_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round2/SNP_curate'
fastq_name = '.curate.flt.snp.vcf'
samplename = '_1.fastq'
genome_name = '.fasta'
species_select = ''

# set up cutoff
end_cutoff = 10 # don't curate the first and last end_cutoff bp of a contig

# function

def loadreport(vcf_file_report,coverage_report):
    assembly_quality = dict()
    assembly_quality_chr = dict()
    report_curate = dict()
    coverage = dict()
    remove_list = []
    try:
        for lines in open(coverage_report,'r'):
            lines_set = lines.split('\t')
            genome = lines_set[0]
            quality = lines_set[4]
            coverage.setdefault(genome,quality)
    except FileNotFoundError:
        pass
    for lines in open(vcf_file_report,'r'):
        if species_select == '' or species_select in lines:
            lines_set = lines.split('\t')
            genome = lines_set[0]
            if coverage.get(genome,'') == 'qualified' or coverage == {}:
                assembly_quality.setdefault(genome,[])
                CHR = lines_set[1]
                POS = lines_set[2]
                REF = lines_set[3]
                ALT_major = lines_set[8]
                Assembly = lines_set[5]
                Mapping_quality = lines_set[6]
                Need_curation = lines_set[7]
                genome_chr = '%s_%s'%(genome,CHR)
                report_curate['%s_%s_%s' % (genome, CHR, POS)] = lines
                if Need_curation == 'T':
                    # can be curated
                    assembly_quality_chr.setdefault(genome_chr,set())
                    assembly_quality_chr[genome_chr].add('__'.join([POS, REF, ALT_major]))
                    assembly_quality[genome].append(CHR)
                elif Assembly == 'F' or Mapping_quality == 'F':
                    # wrong assembly or bad mapping quality
                    assembly_quality_chr.setdefault(genome_chr, set())
                    assembly_quality_chr[genome_chr].add('__'.join([POS, REF, 'N'])) # use ambiguous characters instead
                    assembly_quality[genome].append(CHR)
                else:
                    print(lines)
            else:
                remove_list.append(genome)
    return [assembly_quality,assembly_quality_chr,report_curate]

def checkREF(seq,position,REF):
    return seq[position - 1].upper() == REF.upper()

def causeSNP(seq,position,ALT):
    seq = list(seq)
    seq[position - 1] = ALT
    return ''.join(seq)

def correct_genome(database,assembly_quality_genome,assembly_quality_chr,report_curate):
    Output = []
    Output_report = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        total_length = len(record_seq)
        if record_id in assembly_quality_genome:
            genome_chr = '%s_%s' % (genome, record_id)
            allposition = assembly_quality_chr.get(genome_chr,set())
            for a_position in allposition:
                POS, REF, ALT = a_position.split('__')
                POS = int(POS)
                if POS > end_cutoff and POS < total_length-end_cutoff + 1:
                    # don't curate the first and last end_cutoff bp of a contig
                    if not checkREF(record_seq,int(POS),REF):
                        # wrong POS that needs to be checked
                        print(genome_chr, allposition)
                        print('wrong SNP POS %s %s for %s %s in %s'%(POS,REF,record_seq[POS - 1],
                                                                  record_id,database))
                    record_seq = causeSNP(record_seq, int(POS), ALT)
                    Output_report.append(report_curate['%s_%s_%s' % (genome, record_id, POS)])
        Output.append('>%s\n%s\n'%(record_id,record_seq))
    f1 = open(database + '.corrected.fasta','w')
    f1.write(''.join(Output))
    f1.close()
    return Output_report

# load report
assembly_quality,assembly_quality_chr,report_curate = loadreport(os.path.join(input_script, 'SNP_currate1.assembly.sum'),
           os.path.join(input_script, 'SNP_currate1.coverage.sum'))

# output curated loci
#f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.sum.curated'), 'w')
#f1.write('donor_species\tCHR\tPOS\tREF\tALT\tAssembly\tMapping_quality\tNeed_curation\tMajor_ALT\tMajor_ALT_frq\tDepth_F\tDepth_R\n')
#f1.close()

# correct genomes
cmds = ''
Output_report = []
for genome in assembly_quality:
    assembly_quality_genome = assembly_quality.get(genome, [])
    if assembly_quality_genome != []:
        database = glob.glob(os.path.join(genome_root, genome.split(fastq_name)[0].replace(samplename, genome_name)))
        if database != []:
            database = database[0]
            Output_report += correct_genome(database, assembly_quality_genome, assembly_quality_chr, report_curate)
            # clean up old genomes
            cmds += 'rm -rf %s\n' % (database)
        else:
            print('missing input fasta for %s' % (genome))
    else:
        os.system('cp %s %s.corrected.fasta'%(genome,genome))

f1 = open(os.path.join(input_script, 'cleanoldgenome.sh'), 'a')
f1.write('#!/bin/bash\nsource ~/.bashrc\n' + 'rm -rf %s\n'%(output_dir))
f1.write(''.join(cmds))
f1.close()

f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.sum.curated'), 'a')
f1.write(''.join(Output_report))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# CHANGED filtering genome SNP
# CHANGED ALL curate_REF, outputvcf
# round 1-2 clonal population selection vcf filering, and round 3 SNP calling, round 4 depth filtering
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import itertools
import random
# set up path
Round = 4
Cluster = True
Tree = True
Paircompare = False
Cov_dis = 20

#input_script_temp = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round%s_allpair'%(Round)
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round%s_tree'%(Round)
input_script_sub_merge = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round%s'%(Round)
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay'
genome_root = '/scratch/users/anniz44/genomes/donor_species/jay/round%s'%(Round)
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/round%s/*'%(Round))
output_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/bwa/0/'%(Round)
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/merge'%(Round)
vcf_name = '.all.flt.snp.vcf'
ref_filename = '.all.spades%s.fasta'%(Round)
fasta_name = '.fasta.corrected.fasta'
fastq_name = '.sorted.bam'
deleting_file = []
if Round == 4:
    genome_root = '/scratch/users/anniz44/genomes/donor_species/jay/round*'
    output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/merge_genome/vcf' % (Round)

#input_script_temp = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_round%s_allpair'%(Round)
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_round%s_tree'%(Round)
input_script_sub_merge = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_round%s'%(Round)
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
input_script2 = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round%s'%(Round)
genome_root2 = '/scratch/users/anniz44/genomes/donor_species/jay'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round%s/*'%(Round))
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/bwa/0/'%(Round)
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/merge'%(Round)
vcf_name = '.all.flt.snp.vcf'
ref_filename = '.all.spades%s.fasta'%(Round)
fasta_name = '.fasta.corrected.fasta'
deleting_file = []
fastq_name = '_1.fastq'
Species_replace = dict()
Species_replace.setdefault('BA','Bifidobacterium_adolescentis')
Species_replace.setdefault('BL','Bifidobacterium_longum')
Species_replace.setdefault('PB','Parabacteroides_butyrate')
if Round == 4:
    genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round*'
    output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/merge_genome/vcf' % (
        Round)

# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def vcf_to_txt(lines,output_list,cluster_sub=[]):
    lines_set = lines.split('\n')[0].split('\t')
    if len(lines_set) >9:
        CHR = lines_set[0]
        POS = int(lines_set[1])
        temp_line = []
        temp_line.append(CHR)
        temp_line.append(str(POS))
        i = 9
        for Subdepth_all in lines_set[9:]:
            if (cluster_sub==[] and i not in deleting_set) or i in cluster_sub:
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                temp_line.append(str(total_sub_depth))
            i += 1
        output_list.append('\t'.join(temp_line)+'\n')
    else:
        print(lines)

def curate_REF(allels_set,lines_set_sub):
    SNP_count_genome_count = [0,0,0,0]
    for Subdepth_all in lines_set_sub:
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '')
        Subdepth_set = Subdepth.split(',')
        Subdepth_set_int = []
        for sub_depth in Subdepth_set:
            Subdepth_set_int.append(int(sub_depth))
        if sum(Subdepth_set_int) > 0:
            major_alt_frq = max(Subdepth_set_int)
            major_alt_frq_index = Subdepth_set_int.index(major_alt_frq)
            major_alt = allels_set[major_alt_frq_index]
            SNP_count_genome_count[major_alt_frq_index] += 1
    major_alt_frq = max(SNP_count_genome_count)
    major_alt_frq_index = SNP_count_genome_count.index(major_alt_frq)
    major_alt = allels_set[major_alt_frq_index]
    return [major_alt, major_alt_frq_index]

def outputvcf(output_name):
    vcf_file_filtered = open(vcf_file + '.%s.vcf' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list_vcf))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.snp.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.samplename.txt' % (output_name), 'w')
    vcf_file_filtered.write('\t'.join(Sample_name) + '\n')
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.POS.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_POS))
    vcf_file_filtered.close()

def outputcov(output_name,vcf_file_POS_candidate,cluster_sub=[]):
    if len(vcf_file_list) > 0:
        vcf_file_POS_candidate = '\n'.join(vcf_file_POS_candidate)
        vcf_file_POS_candidate_output = ('%s' % (vcf_file_POS_candidate))
        f1 = open(os.path.join(input_script,'grep.temp.txt'),'w')
        f1.write(vcf_file_POS_candidate_output)
        f1.close()
        if '.fna.flt.snp.vcf' in vcf_file:
            cov_file = vcf_file.split('.fna.flt.snp.vcf')[0] + '.fna.raw.vcf'
        else:
            cov_file = vcf_file.split('.flt.snp.vcf')[0] + '.raw.vcf'
        os.system('grep -%s -f %s %s --no-group-separator > %s'% (
            Cov_dis, os.path.join(input_script,'grep.temp.txt'),
            cov_file,
            vcf_file + '.%s.cov.temp' % (output_name)))
        os.system('cat %s | sort | uniq > %s' % (
            vcf_file + '.%s.cov.temp' % (output_name),
            vcf_file + '.%s.uniqcov.temp' % (output_name)))
        for lines in open(vcf_file + '.%s.uniqcov.temp' % (output_name), 'r'):
            if not lines.startswith("#"):
                vcf_to_txt(lines, cov_file_list,cluster_sub)
        os.system('rm -rf %s %s %s' % (vcf_file + '.%s.cov.temp' % (output_name),
                                    vcf_file + '.%s.uniqcov.temp' % (output_name),
                                       os.path.join(input_script, 'grep.temp.txt')))
        vcf_file_filtered = open(vcf_file + '.%s.cov.txt' % (output_name), 'w')
        vcf_file_filtered.write(''.join(cov_file_list))
        vcf_file_filtered.close()

def outputtree(output_name):
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, SNP_alignment[genomename]))
            SNP_alignment_output_parsi.append('S%s    %s\n' % (genomename[-8:], SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.%s.fasta' % (output_name), 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    if Tree:
        vcf_file_filtered = open(vcf_file + '.%s.parsi.fasta' % (output_name), 'w')
        vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
        vcf_file_filtered.close()

def SNP_seq(seq1, seq2, POS_info,POS_info_CHR,POS_info_CHR_LEN,POS_info_output,G1,G2):
    SNP_total = 0
    j = 0
    POS_DIS = []
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            CHR = POS_info_CHR[i]
            POS = POS_info[i]
            LEN = POS_info_CHR_LEN[CHR]
            if CHR == POS_info_CHR[j]:  # same CHR
                DIS = abs(POS - POS_info[j])
                POS_DIS.append(DIS)  # POS diff
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, DIS, LEN))
            else:  # new CHR
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, 0, LEN))
            j = i
    return SNP_total

def SNP_distance_correct(distanceNEW, distanceREF, SNPREF):
    return int((distanceNEW/distanceREF)*SNPREF)

def find_neighbor(Cluster_SNP,neighbor,Cluster_SNP_set,cluster,Cluster_SNP_set_added):
    if neighbor != []:
        for record_name in neighbor:
            if record_name not in Cluster_SNP_set_added:
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    subneighbor = Cluster_SNP.get(record_name,[])
                    find_neighbor(Cluster_SNP, subneighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)

def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def dnORds(amino1, amino2):
    if amino1 == amino2:
        return 'S'
    else:
        return 'N'

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

def contig_to_gene(CHR, POS):
    all_genes = Mapping_loci.get(CHR,[])
    Reverse_chr = 0
    if all_genes == []:
        # database is a gene database
        codon_start = POS - 1 - int((POS - 1) % 3)
        Ref_seq_chr = Ref_seq.get(CHR, 'None')
        return [CHR, POS, codon_start, Ref_seq_chr, Reverse_chr]
    else:
        # database is a contig database
        for a_gene in all_genes:
            POS1, POS2, GENE = a_gene
            if POS >= POS1 and POS <= POS2:
                Ref_seq_chr = Ref_seq.get(GENE, 'None')
                Gene_length = len(Ref_seq_chr)
                if GENE in Reverse:  # reversed
                    POS_gene = Gene_length-(int(POS-POS1))
                    Reverse_chr = 1
                else:
                    POS_gene = int(POS-POS1)+1
                codon_start = POS_gene - 1 - int((POS_gene - 1) % 3)
                return [GENE,POS_gene,codon_start,Ref_seq_chr,Reverse_chr]
    return []

def contig_end(CHR,POS):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        total_length = CHR.split('length_')[1].split('_cov')[0]
    total_length = int(total_length)
    if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
        return True
    else:
        return False

def SNP_check_all(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name,SNP_presence_cutoff,SNP_presence_sample_cutoff,no_SNP_cutoff,Depth_set,cluster_sub=[]):
    CHR = lines_set[0]
    POS = int(lines_set[1])
    CHRPOS = '%s\t%s'%(CHR,POS)
    if Depth_set == {} or CHRPOS in Depth_set:
        temp_snp_line = []
        temp_snp_line_frq = []
        temp_snp_line_NS = ['None', 'None', 'None']
        temp_snp_line_AA = ''
        Total_qualify = 0
        Total_qualify_SNP = 0
        Total_qualify_notSNP = 0
        Total_unqualify_alt_freq = 0
        SNP = set()
        SNP_seq = []
        REF = lines_set[3]
        allels_set = [REF]
        Total_subsample = Total
        lines_set_sub = lines_set[9:]
        REF_where=0
        if cluster_sub!= []:
            lines_set_sub = [lines_set[i] for i in cluster_sub]
            Total_subsample = len(cluster_sub)
            if Total_subsample >= 15:
                SNP_presence_cutoff = 0.33  # for a large group of samples
            elif Total_subsample in [3,4]:
                SNP_presence_cutoff = 1  # for a small group of samples
                SNP_presence_sample_cutoff = 2
            elif Total_subsample in [1,2]:
                SNP_presence_cutoff = 1  # for a small group of samples
                SNP_presence_sample_cutoff = 1
                no_SNP_cutoff = 0
        else:
            cluster_sub = list(range(9, len(lines_set)))
        if Total_subsample > 0:
            if '.' not in lines_set[4]:
                allels_set += lines_set[4].split(',')
            Total_alleles = len(allels_set)
            genome_order = 0
            Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
            if Total_subsample > 2:
                REF,REF_where = curate_REF(allels_set, lines_set_sub)  # as the major alt in the population
            sample_num = 9
            for Subdepth_all in lines_set_sub:
                if sample_num not in deleting_set:
                    genome_order += 1
                    Allels_frq = [0, 0, 0, 0]
                    Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                    total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                    for num_allels in range(0, Total_alleles):
                        allels = allels_set[num_allels]
                        Subdepth_alleles = int(Subdepth[num_allels])
                        if allels in Allels:
                            Allels_frq[Allels[allels]] += Subdepth_alleles
                        else:
                            pass
                    # find major alt and calculate frequency
                    Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                    temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq))
                    SNP_seq.append(REF)  # set as reference
                    if total_sub_depth > 0:
                        MLF = Major_ALT[1] / total_sub_depth
                        if MLF >= Major_alt_freq_cutoff:
                            # major alt frequency cutoff
                            Total_qualify += 1
                            # check for qualified SNP
                            if Major_ALT[0] != REF:
                                if (Depth_set == {} or sample_num - 9 in Depth_set[CHRPOS]): # depth cutoff
                                    # a qualified SNP
                                    temp_snp_line_pass += 'PASS'
                                    Total_qualify_SNP += 1
                                    SNP.add(genome_order)  # only take qualified SNP as valid SNP
                                    SNP_seq[-1] = Major_ALT[0] # only qualified SNP include in alignment
                                else:
                                    print('deleting sample %s of CHRPOS %s' % (sample_num - 9, CHRPOS))
                            else:
                                Total_qualify_notSNP += 1
                        else:
                            # major alt frequency low
                            Total_unqualify_alt_freq += 1
                sample_num += 1
            if Total_qualify / Total_subsample >= SNP_presence_cutoff and \
                    Total_unqualify_alt_freq < Poor_MLF_freq_cutoff and\
                    Total_qualify >= SNP_presence_sample_cutoff and \
                    Total_qualify_SNP >= 1 and Total_qualify_SNP <= Total_qualify - no_SNP_cutoff and\
                    Total_qualify_notSNP >= no_SNP_cutoff:
                # -> qualified SNP, qualified samples cutoff + unqualified samples cutoff, at least 1 qualified SNP, calculate NS
                gene_info = contig_to_gene(CHR, POS)
                if gene_info!= []:
                    Chr_gene, POS_gene,codon_start,Ref_seq_chr,Reverse_chr  = gene_info
                    if Ref_seq_chr != 'None':
                        #  observed NS ratio calculated
                        temp_snp_line_NS= [Chr_gene,str(POS_gene),'']
                        if codon_start <= POS_gene - 1:
                            Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                            Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                            SNP_seq_chr = Ref_seq_chr
                            if len(Ref_seq_codon) == 3:
                                Ref_seq_aa = translate(Ref_seq_codon)[0]
                                temp_snp_line_AA += Ref_seq_aa
                                ALT_set = allels_set
                                ALT_set.remove(REF)
                                for ALT in ALT_set:
                                    SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                    SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                    SNP_seq_aa = translate(SNP_seq_codon)[0]
                                    temp_snp_line_AA += SNP_seq_aa
                                    temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                    temp_snp_line_NS[-1]+=temp_NorS
                # output lines and output major alt
                if 'PASS' in temp_snp_line_pass:
                    temp_snp_line_pass = 'PASS'
                else:
                    temp_snp_line_pass = 'NOPASS'
                if CHR == CHR_old:
                    # same CHR
                    POS_DIS = abs(POS - POS_old)
                    vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
                else:
                    # diff CHR first SNP
                    vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
                POS_old = POS
                CHR_old = CHR
                temp_snp_line.append(CHR)
                temp_snp_line.append(str(POS))
                temp_snp_line.append(REF)
                temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
                vcf_file_list.append('\t'.join(temp_snp_line)+ '\t' +'\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                    ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,'\t'.join(temp_snp_line_NS),temp_snp_line_AA))
                vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
                vcf_file_list_vcf.append('\t'.join(lines_set[0:9])+'\t'+'\t'.join(lines_set_sub)+'\n')
                i = 9
                j = 0
                SNP_alignment[reference_name] += REF
                for genomename in SNP_alignment:
                    if genomename != reference_name:
                        if i in cluster_sub:
                            SNP_alignment[genomename] += SNP_seq[j]
                            j += 1
                        i += 1
    else:
        print('deleting CHRPOS %s'%(CHRPOS))
    return [CHR_old,POS_old]

def SNP_check_output(lines_set,CHR_old,POS_old,reference_name):
    CHR = lines_set[0]
    POS = int(lines_set[1])
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    REF_where = 0
    cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        if Total_subsample > 2:
            REF, REF_where = curate_REF(allels_set, lines_set_sub)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    qualify_loci = 0
                    MLF = Major_ALT[1] / total_sub_depth
                    if Major_ALT[0] != REF and MLF >= Major_alt_freq_cutoff or qualify_loci == 1:
                        # check for qualified SNP
                        SNP.add(genome_order)  # only take qualified SNP as valid SNP
                        SNP_seq[-1] = Major_ALT[0]  # only qualified SNP include in alignment
            sample_num += 1
        gene_info = contig_to_gene(CHR, POS)
        if gene_info != []:
            Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
            if Ref_seq_chr != 'None':
                #  observed NS ratio calculated
                temp_snp_line_NS = [Chr_gene, str(POS_gene), '']
                if codon_start <= POS_gene - 1:
                    Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                    Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                    SNP_seq_chr = Ref_seq_chr
                    if len(Ref_seq_codon) == 3:
                        Ref_seq_aa = translate(Ref_seq_codon)[0]
                        temp_snp_line_AA += Ref_seq_aa
                        ALT_set = allels_set
                        ALT_set.remove(REF)
                        for ALT in ALT_set:
                            SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                            SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                            SNP_seq_aa = translate(SNP_seq_codon)[0]
                            temp_snp_line_AA += SNP_seq_aa
                            temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                            temp_snp_line_NS[-1] += temp_NorS
        temp_snp_line_pass = 'PASS'
        if CHR == CHR_old:
            # same CHR
            POS_DIS = abs(POS - POS_old)
            vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
        else:
            # diff CHR first SNP
            vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
        POS_old = POS
        CHR_old = CHR
        temp_snp_line.append(CHR)
        temp_snp_line.append(str(POS))
        temp_snp_line.append(REF)
        temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
        vcf_file_list.append(
            '\t'.join(temp_snp_line) + '\t' + '\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,
                '\t'.join(temp_snp_line_NS), temp_snp_line_AA))
        vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
        vcf_file_list_vcf.append('\t'.join(lines_set[0:9]) + '\t' + '\t'.join(lines_set_sub) + '\n')
        i = 9
        j = 0
        SNP_alignment[reference_name] += REF
        for genomename in SNP_alignment:
            if genomename != reference_name:
                if i in cluster_sub:
                    SNP_alignment[genomename] += SNP_seq[j]
                    j += 1
                i += 1
    return [CHR_old,POS_old]

def dis_pos(current_CHRPOS,pre_CHRPOS):
    POS1 = int(pre_CHRPOS.split('\t')[1])
    POS2 = int(current_CHRPOS.split('\t')[1])
    return POS2 - POS1

def compare_set(Ge_pre,Ge_cu):
    Ge_preset = Ge_pre.replace('\"','').split(';')
    Ge_cuset = Ge_cu.replace('\"', '').split(';')
    if (len(Ge_preset) >= 3 and len(Ge_cuset) >= 2) or (len(Ge_preset) >= 2 and len(Ge_cuset) >= 3):
        #return all(elem in Ge_cuset for elem in Ge_preset) or all(elem in Ge_preset for elem in Ge_cuset)
        # allowing one mismatch
        return len([elem for elem in Ge_cuset if elem in Ge_preset]) >= len(Ge_cuset) - 1 or \
               len([elem for elem in Ge_preset if elem in Ge_cuset]) >= len(Ge_preset) - 1
    else:
        return False

def cluster_rec(CHR_set,SNP_gen_set,SNP_N_set,CHRPOS_set):
    for CHR in CHR_set:
        potential_rec = dict()
        allCHRPOS = CHR_set[CHR]
        for i in range(1,len(allCHRPOS)):
            current_CHRPOS = allCHRPOS[i]
            for j in reversed(range(0, i)):
                pre_CHRPOS = allCHRPOS[j]
                Ge_pre = SNP_gen_set[pre_CHRPOS]
                Ge_cu = SNP_gen_set[current_CHRPOS]
                Dis = dis_pos(current_CHRPOS,pre_CHRPOS)
                if Dis <= Rec_length_cutoff:
                    if Ge_pre == Ge_cu or compare_set(Ge_pre,Ge_cu):
                        # cluster rec sites
                        potential_rec.setdefault(Ge_pre, set())
                        potential_rec[Ge_pre].add(current_CHRPOS)
                        potential_rec[Ge_pre].add(pre_CHRPOS)
                        potential_rec.setdefault(Ge_cu, set())
                        potential_rec[Ge_cu].add(current_CHRPOS)
                        potential_rec[Ge_cu].add(pre_CHRPOS)
                        potential_rec[Ge_cu].update(list(potential_rec[Ge_pre]))
                        potential_rec[Ge_pre].update(list(potential_rec[Ge_cu]))
                        break
                else:
                    if j == i-1:
                        # output recombination
                        delete_set = []
                        for Ge_pre in potential_rec:
                            allrec = potential_rec[Ge_pre]
                            if len(allrec) >= Rec_SNP_cutoff:
                                CHRPOS_set += allrec
                            delete_set.append(Ge_pre)
                        for Ge_pre in delete_set:
                            potential_rec.pop(Ge_pre, 'None')
                    break
        # output the last recombination
        for Ge_pre in potential_rec:
            allrec = potential_rec[Ge_pre]
            if len(allrec) >= Rec_SNP_cutoff:
                CHRPOS_set += allrec
    return CHRPOS_set

def remove_rec(SNP_file):
    CHRPOS_set = []
    CHR_set = dict()
    SNP_gen_set = dict()
    SNP_N_set = dict()
    # import SNP info
    for lines in open(SNP_file,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        CHR = lines_set[0]
        POS = lines_set[1]
        CHRPOS = '%s\t%s' % (CHR, POS)
        Ge = lines_set[-6]
        NS = lines_set[-2]
        if 'NN' in NS:
            NS = 'N'
        SNP_gen_set.setdefault(CHRPOS,Ge)
        SNP_N_set.setdefault(CHRPOS, NS)
        CHR_set.setdefault(CHR, [])
        CHR_set[CHR].append(CHRPOS)
    CHRPOS_set = cluster_rec(CHR_set,SNP_gen_set,SNP_N_set,CHRPOS_set)
    return CHRPOS_set

def cluster_event(CHR_set, SNP_gen_set, SNP_N_set, CHRPOS_set):
    for CHRgene in CHR_set:
        allCHRPOS = CHR_set[CHRgene]
        allge = SNP_gen_set[CHRgene]
        allN = SNP_N_set[CHRgene]
        uniqge = set(allge)
        for ge in uniqge:
            if allge.count(ge) > 1:
                # at least 2 SNPs in the same genome set
                samegenomeset = [i for i in range(0, len(allge)) if allge[i] == ge]
                removeCHRPOS = [allCHRPOS[i] for i in samegenomeset]
                N_set = [allN[i] for i in samegenomeset]
                try:
                    keepCHRPOS = removeCHRPOS[N_set.index('N')]
                except ValueError:
                    keepCHRPOS = removeCHRPOS[0]
                removeCHRPOS.remove(keepCHRPOS)
                CHRPOS_set += removeCHRPOS
    return CHRPOS_set

def remove_event(SNP_file):
    CHRPOS_set = []
    CHR_set = dict()
    SNP_gen_set = dict()
    SNP_N_set = dict()
    # import SNP info
    for lines in open(SNP_file,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        CHR = lines_set[0]
        POS = lines_set[1]
        CHRPOS = '%s\t%s' % (CHR, POS)
        Ge = lines_set[-6]
        NS = lines_set[-2]
        if 'NN' in NS:
            NS = 'N'
        CHRgene = lines_set[-4]  # gene based
        if CHRgene != 'None':
            CHR_set.setdefault(CHRgene, [])
            CHR_set[CHRgene].append(CHRPOS)
            SNP_gen_set.setdefault(CHRgene,[])
            SNP_gen_set[CHRgene].append(Ge)
            SNP_N_set.setdefault(CHRgene, [])
            SNP_N_set[CHRgene].append(NS)
    CHRPOS_set = cluster_event(CHR_set, SNP_gen_set, SNP_N_set, CHRPOS_set)
    return CHRPOS_set

def curate_REF_fq(allels_set,lines_set_sub):
    SNP_count_genome_count = [0,0,0,0]
    Good_REF = False
    Good_ALT = False
    REF = allels_set[0]
    for Subdepth_all in lines_set_sub:
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '')
        Subdepth_set = Subdepth.split(',')
        Subdepth_set_int = []
        for sub_depth in Subdepth_set:
            Subdepth_set_int.append(int(sub_depth))
        total_sub_depth = sum(Subdepth_set_int)
        if total_sub_depth > 0:
            major_alt_frq = max(Subdepth_set_int)
            major_alt_frq_index = Subdepth_set_int.index(major_alt_frq)
            major_alt = allels_set[major_alt_frq_index]
            SNP_count_genome_count[major_alt_frq_index] += 1
            MLF = major_alt_frq / total_sub_depth
            if MLF >= MLF_cutoff and total_sub_depth >= Depth_cutoff:
                # MLF cutoff and depth cutoff
                if major_alt!=REF:
                    Good_ALT = True
                else:
                    Good_REF = True
    major_alt_frq = max(SNP_count_genome_count)
    major_alt_frq_index = SNP_count_genome_count.index(major_alt_frq)
    major_alt = allels_set[major_alt_frq_index]
    return [major_alt, major_alt_frq_index,Good_REF,Good_ALT]

def load_ref_vcf(ref_vcf_file):
    ref_chr = dict()
    Length = dict()
    for files in ref_vcf_file:
        for lines in open(files,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS, Notused, REF, ALT = lines_set[0:5]
            CHR_POS = '%s__%s'%(CHR, POS)
            if CHR not in Length:
                try:
                    total_length = CHR.split('size')[1]
                except IndexError:
                    total_length = CHR.split('length_')[1].split('_cov')[0]
                total_length = int(total_length)
                Length.setdefault(CHR, total_length)
            total_length = Length[CHR]
            if total_length >= Length_cutoff:
                # CONTIG LENGTH cutoff
                lines_set_sub = lines_set[9:]
                allels_set = [REF]
                if '.' not in lines_set[4]:
                    allels_set += lines_set[4].split(',')
                REF, REF_where,Good_REF,Good_ALT = curate_REF_fq(allels_set, lines_set_sub)  # as the major alt in the population
                if Good_REF and Good_ALT:
                    ref_chr.setdefault(CHR_POS, [])
                    ref_chr[CHR_POS] = [REF, ALT]
    return ref_chr

# set up output
if Tree:
    try:
        os.mkdir(output_dir_merge + '/tree')
    except IOError:
        pass
    try:
        os.mkdir(input_script_sub)
    except IOError:
        pass

# set up cutoff
reference_set = ['reference']
outputname_set = ['filtered']
SNP_presence_cutoff = 0.66  # avg presence in all samples
SNP_presence_sample_cutoff = 3  # num of samples passing the above criteria
Major_alt_freq_cutoff = 0.9 # major alt freq in a genome, do not allow multiple homolougs genes
no_SNP_cutoff = 1
Poor_MLF_freq_cutoff = 1 # no sample should have homologous genes (low major alt freq)
# set up strict cutoff
SNP_presence_cutoff2 = 0.66 # avg coverage in all samples
SNP_presence_sample_cutoff2 = 3  # num of samples passing the above criteria
Poor_MLF_freq_cutoff2 = 1 # no sample should have homologous genes (low major alt freq)
# cluster cutoff Step 3
SNP_total_cutoff_2 = 100
cluster_cutoff = 2
# Depth and recombination cutoff Round4
Depth_cutoff = 6 # covered by 6 reads
MLF_cutoff = 0.7 # MLF cutoff for WGS
Length_cutoff = 2000 # minimum ref contig length
Rec_length_cutoff = 5000 # maximum distance between recombination sites
Rec_SNP_cutoff = 3 # minumum no. of SNPs grouped/clustered as a recombination
end_cutoff = 70 # contig end no SNP calling

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

# run vcf filtering
if Round == 4:
    vcf_name = '.all.flt.snp.vcf.filtered.vcf'
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    outputname_set = ['final']
    ref_filename = '.all.spades*.fasta'
else:
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))

reference_name = reference_set[0]
output_name = outputname_set[0]
for vcf_file in all_vcf_file:
    filesize = 0
    try:
        filesize = int(os.path.getsize(vcf_file + '.%s.vcf2'%(output_name)))
    except FileNotFoundError:
        pass
    if filesize == 0:
        SNP_presence_cutoff = SNP_presence_cutoff2  # for group of samples
        SNP_presence_sample_cutoff = SNP_presence_sample_cutoff2
        no_SNP_cutoff = 1
        print(vcf_file)
        Total = 0
        # filter depth
        Depth_set = dict()
        ref_chr = dict()
        if Round == 4:
            vcf_ref_file_name = os.path.split(vcf_file)[-1].split('.all.')[0]
            vcf_file_raw = \
            glob.glob(output_dir_merge + '/../../../vcf_round*/merge/' + vcf_ref_file_name + '.all.raw.vcf')[0]
            try:
                # WGS
                vcf_fq = glob.glob(output_dir_merge + '/../../merge/' + vcf_ref_file_name + '*.all.fq.flt.snp.vcf')
                ref_chr = load_ref_vcf(vcf_fq)
            except IndexError:
                pass
        else:
            vcf_file_raw = vcf_file.replace('.flt.snp.vcf', '.raw.vcf')
        donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0].split('.flt.snp.vcf')[0]
        database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
        if len(database) > 1:
            print(vcf_file,database)
        database = database[0]
        ref_dir, ref_name = os.path.split(database)
        ref_fna = database.replace('.fasta', '.fna')
        try:
            f1 = open(ref_fna, 'r')
        except FileNotFoundError:
            os.system('prodigal -q -i %s -d %s' % (database, ref_fna))
        Sample_name = []
        deleting_set = []
        Ref_seq = dict()
        Mapping = dict()
        Mapping_loci = dict()
        for lines in open(vcf_file_raw, 'r'):
            if lines.startswith('##bcftoolsCommand=mpileup '):
                # setup samples
                sample_set = lines.split(ref_name + ' ')[1].split('\n')[0].split(' |')[0].split(' ')
                samplenum = 9
                for samples in sample_set:
                    genomename = os.path.split(samples)[-1].split(fastq_name)[0]
                    Sample_name.append(genomename.replace('.', ''))
                    if genomename in deleting_file:
                        deleting_set.append(samplenum)
                    samplenum += 1
                break
        print('running %s' % (donor_species))
        # load database
        database_file = ref_fna
        Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
        SNP_tree_cmd = []
        SNP_tree_cmd2 = []
        vcf_file_list = []
        vcf_file_list_vcf = []
        vcf_file_POS = []
        vcf_file_POS_candidate = set()
        SNP_alignment = dict()
        SNP_alignment.setdefault(reference_name, '')
        cov_file_list = []
        CHR_old = ''
        POS_old = 0
        for genomename in Sample_name:
            SNP_alignment.setdefault(genomename, '')
        for lines in open(vcf_file, 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\n')[0].split('\t')
                CHR = lines_set[0]
                POS = int(lines_set[1])
                # a SNP confirmed in WGS mapping
                CHR_POS = '%s__%s' % (CHR, POS)
                if (Round < 4 or CHR_POS in ref_chr) and not contig_end(CHR, POS):
                    Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                    if Total == 0:
                        Total = len(lines_set) - 9 - len(deleting_set)
                        if Total >= 15:
                            SNP_presence_cutoff = 0.33  # for a large group of genomes
                        elif Total in [3, 4]:
                            SNP_presence_cutoff = 1  # for a small group of genomes
                            SNP_presence_sample_cutoff = 2
                        elif Total in [1, 2]:
                            SNP_presence_cutoff = 1  # for only 1 or 2 samples, compare to ref
                            SNP_presence_sample_cutoff = 1
                            no_SNP_cutoff = 0
                    if Depth / Total >= SNP_presence_cutoff:
                        # average depth in all samples cutoff
                        if "INDEL" not in lines_set[7] \
                                and (lines_set[6] != 'LowQual'):
                            CHR_old, POS_old = SNP_check_all(lines_set, '',
                                                             CHR_old, POS_old, reference_name,
                                                             SNP_presence_cutoff,
                                                             SNP_presence_sample_cutoff, no_SNP_cutoff, Depth_set)
        outputvcf(output_name)
        if Round < 4:
            outputtree(output_name)

# remove recombination
if Round == 4:
    vcf_name = '.all.flt.snp.vcf.filtered.vcf.final.vcf'
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    outputname_set = ['removerec']
    ref_filename = '.all.spades*.fasta'
    reference_name = reference_set[0]
    output_name = outputname_set[0]
    for vcf_file in all_vcf_file:
        filesize = 0
        try:
            filesize = int(os.path.getsize(vcf_file + '.%s.vcf2' % (output_name)))
        except FileNotFoundError:
            pass
        if filesize == 0:
            print(vcf_file)
            Total = 0
            # filter recombination
            vcf_ref_file_name = os.path.split(vcf_file)[-1].split('.all.')[0]
            vcf_file_raw = \
                glob.glob(output_dir_merge + '/../../../vcf_round*/merge/' + vcf_ref_file_name + '.all.raw.vcf')[0]
            SNP_file = vcf_file.replace('.final.vcf', '.final.snp.txt')
            CHRPOS_set = remove_rec(SNP_file)
            donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0].split('.flt.snp.vcf')[0]
            database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
            if len(database) > 1:
                print(vcf_file, database)
            database = database[0]
            ref_dir, ref_name = os.path.split(database)
            ref_fna = database.replace('.fasta', '.fna')
            try:
                f1 = open(ref_fna, 'r')
            except FileNotFoundError:
                os.system('prodigal -q -i %s -d %s' % (database, ref_fna))
            Sample_name = []
            deleting_set = []
            Ref_seq = dict()
            Mapping = dict()
            Mapping_loci = dict()
            for lines in open(vcf_file_raw, 'r'):
                if lines.startswith('##bcftoolsCommand=mpileup '):
                    # setup samples
                    sample_set = lines.split(ref_name + ' ')[1].split('\n')[0].split(' |')[0].split(' ')
                    samplenum = 9
                    for samples in sample_set:
                        genomename = os.path.split(samples)[-1].split(fastq_name)[0]
                        Sample_name.append(genomename.replace('.', ''))
                        if genomename in deleting_file:
                            deleting_set.append(samplenum)
                        samplenum += 1
                    break
            print('running %s' % (donor_species))
            # load database
            database_file = ref_fna
            Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
            SNP_tree_cmd = []
            SNP_tree_cmd2 = []
            vcf_file_list = []
            vcf_file_list_vcf = []
            vcf_file_POS = []
            vcf_file_POS_candidate = set()
            SNP_alignment = dict()
            SNP_alignment.setdefault(reference_name, '')
            cov_file_list = []
            CHR_old = ''
            POS_old = 0
            for genomename in Sample_name:
                SNP_alignment.setdefault(genomename, '')
            for lines in open(vcf_file, 'r'):
                if not lines.startswith("#"):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR = lines_set[0]
                    POS = int(lines_set[1])
                    CHRPOS = '%s\t%s' % (CHR, POS)
                    if Total == 0:
                        Total = len(lines_set) - 9 - len(deleting_set)
                    if CHRPOS not in CHRPOS_set and not contig_end(CHR, POS):
                        CHR_old, POS_old = SNP_check_output(lines_set,
                                                            CHR_old, POS_old, reference_name)
            outputvcf(output_name)
            outputtree(output_name)

# run parsi tree
if Tree:
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    for vcf_file in all_vcf_file:
        a_parsi_file = vcf_file + '.%s.parsi.fasta'%(output_name)
        if 'RAxML_parsimonyTree' not in a_parsi_file:
            os.system('rm -rf %s %s' % (a_parsi_file + '.out.txt',
                                        a_parsi_file + '.out.tree'))
            SNP_tree_cmd3 = ('%s\n5\nV\n1\ny\n' % (a_parsi_file))
            f1 = open(os.path.join(input_script_sub, 'parsi.optionfile.txt'), 'w')
            f1.write(SNP_tree_cmd3)
            f1.close()
            os.system('rm -rf outfile outtree')
            os.system('dnapars < %s/parsi.optionfile.txt > %s/%s.parsi.output\n' % (
                input_script_sub, input_script_sub, os.path.split(a_parsi_file)[-1]))
            os.system('mv outfile %s' % (a_parsi_file + '.out.txt'))
            os.system('mv outtree %s' % (a_parsi_file + '.out.tree'))
    os.system('mv %s/*.parsi* %s/tree' % (
        output_dir_merge, output_dir_merge))

# cluster historical events based on genes
if Round == 4:
    vcf_name = '.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf'
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    outputname_set = ['singleevent']
    ref_filename = '.all.spades*.fasta'
    reference_name = reference_set[0]
    output_name = outputname_set[0]
    for vcf_file in all_vcf_file:
        filesize = 0
        try:
            filesize = int(os.path.getsize(vcf_file + '.%s.vcf2' % (output_name)))
        except FileNotFoundError:
            pass
        if filesize == 0:
            print(vcf_file)
            Total = 0
            # filter recombination
            vcf_ref_file_name = os.path.split(vcf_file)[-1].split('.all.')[0]
            vcf_file_raw = \
                glob.glob(output_dir_merge + '/../../../vcf_round*/merge/' + vcf_ref_file_name + '.all.raw.vcf')[0]
            SNP_file = vcf_file.replace('.removerec.vcf', '.removerec.snp.txt')
            CHRPOS_set = remove_event(SNP_file)
            donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0].split('.flt.snp.vcf')[0]
            database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
            if len(database) > 1:
                print(vcf_file, database)
            database = database[0]
            ref_dir, ref_name = os.path.split(database)
            ref_fna = database.replace('.fasta', '.fna')
            try:
                f1 = open(ref_fna, 'r')
            except FileNotFoundError:
                os.system('prodigal -q -i %s -d %s' % (database, ref_fna))
            Sample_name = []
            deleting_set = []
            Ref_seq = dict()
            Mapping = dict()
            Mapping_loci = dict()
            for lines in open(vcf_file_raw, 'r'):
                if lines.startswith('##bcftoolsCommand=mpileup '):
                    # setup samples
                    sample_set = lines.split(ref_name + ' ')[1].split('\n')[0].split(' |')[0].split(' ')
                    samplenum = 9
                    for samples in sample_set:
                        genomename = os.path.split(samples)[-1].split(fastq_name)[0]
                        Sample_name.append(genomename.replace('.', ''))
                        if genomename in deleting_file:
                            deleting_set.append(samplenum)
                        samplenum += 1
                    break
            print('running %s' % (donor_species))
            # load database
            database_file = ref_fna
            Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
            SNP_tree_cmd = []
            SNP_tree_cmd2 = []
            vcf_file_list = []
            vcf_file_list_vcf = []
            vcf_file_POS = []
            vcf_file_POS_candidate = set()
            SNP_alignment = dict()
            SNP_alignment.setdefault(reference_name, '')
            cov_file_list = []
            CHR_old = ''
            POS_old = 0
            for genomename in Sample_name:
                SNP_alignment.setdefault(genomename, '')
            for lines in open(vcf_file, 'r'):
                if not lines.startswith("#"):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR = lines_set[0]
                    POS = int(lines_set[1])
                    CHRPOS = '%s\t%s' % (CHR, POS)
                    if Total == 0:
                        Total = len(lines_set) - 9 - len(deleting_set)
                    if CHRPOS not in CHRPOS_set and not contig_end(CHR, POS):
                        CHR_old, POS_old = SNP_check_output(lines_set,
                                                            CHR_old, POS_old, reference_name)
            outputvcf(output_name)
            outputtree(output_name)

# run parsi tree
if Tree:
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    for vcf_file in all_vcf_file:
        a_parsi_file = vcf_file + '.%s.parsi.fasta'%(output_name)
        if 'RAxML_parsimonyTree' not in a_parsi_file:
            os.system('rm -rf %s %s' % (a_parsi_file + '.out.txt',
                                        a_parsi_file + '.out.tree'))
            SNP_tree_cmd3 = ('%s\n5\nV\n1\ny\n' % (a_parsi_file))
            f1 = open(os.path.join(input_script_sub, 'parsi.optionfile.txt'), 'w')
            f1.write(SNP_tree_cmd3)
            f1.close()
            os.system('rm -rf outfile outtree')
            os.system('dnapars < %s/parsi.optionfile.txt > %s/%s.parsi.output\n' % (
                input_script_sub, input_script_sub, os.path.split(a_parsi_file)[-1]))
            os.system('mv outfile %s' % (a_parsi_file + '.out.txt'))
            os.system('mv outtree %s' % (a_parsi_file + '.out.tree'))
    os.system('mv %s/*.parsi* %s/tree' % (
        output_dir_merge, output_dir_merge))

# cluster cutoff
SNP_total_cutoff_2 = 1000
cluster_cutoff = 3
# run clustering
second_strain = dict()
output_name = outputname_set[0]
if Cluster:
    # sum up SNP
    all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
    f1 = open(os.path.join(input_script, 'SNP_round%s.sum' % (Round)), 'w')
    f1.write('donor_species\tGenome\tSNP_total\tcluster\tsubcluster\t\n')
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round%s.allpair.sum' % (Round)), 'w')
    f1.write('Genome1\tGenome2\tSNP_total\t\n')
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round%s.sum.cutoff' % (Round)), 'w')
    f1.write('donor_species\tmax_cluster_diff\tSNP_cutoff\tSNP_total_len\t\n')
    f1.close()
    for vcf_file in all_vcf_file:
        fasta = vcf_file + '.%s.fasta'%(output_name)
        SNP_pair = []
        POS_info_output = []
        tree_distance_output = []
        SNP_cutoff = []
        fasta_name = os.path.split(fasta)[-1]
        donor_species = fasta_name.split('.all.flt.snp.vcf')[0]
        POS_file = glob.glob(os.path.join(output_dir_merge, donor_species + '.*.%s.POS.txt'%(output_name)))[0]
        print(donor_species)
        tree_distance = dict()
        Seq_list = dict()
        Ref = ''
        max_cluster_diff = 0
        POS_info = []
        POS_info_CHR = []
        POS_info_CHR_LEN = dict()
        Cluster_SNP = dict()
        Cluster_SNP_set = dict()
        Cluster_SNP_set_added = set()
        # load SNP POS info
        for lines in open(POS_file, 'r'):
            CHR = lines.split('\t')[0]
            POS_info.append(int(lines.split('\t')[1]))
            POS_info_CHR.append(CHR)
            POS_info_CHR_LEN.setdefault(CHR, 0)
            POS_info_CHR_LEN[CHR] = max(POS_info_CHR_LEN[CHR], int(lines.split('\t')[1]))
        # load genome SNP fasta and calculate pair-wise SNPs
        for record in SeqIO.parse(fasta, 'fasta'):
            record_name = str(record.id)
            if 'reference' not in record_name:
                new_center = 1
                record_seq = str(record.seq)
                if Ref == '':
                    # set upt the first seq as ref
                    Ref = record_seq
                    REF_name = record_name
                    new_center = 0
                    SNP_total = 0
                    SNP_total_length = len(Ref)
                    SNP_total_cutoff = SNP_total_cutoff_2
                else:
                    for record_before in Seq_list:
                        SNP_total = SNP_seq(Seq_list[record_before], record_seq, POS_info, POS_info_CHR,
                                            POS_info_CHR_LEN,
                                            POS_info_output, record_before, record_name)
                        SNP_pair.append('%s\t%s\t%s\t\n' % (record_before, record_name, SNP_total))
                        if SNP_total <= SNP_total_cutoff:
                            Cluster_SNP.setdefault(record_before, [])
                            Cluster_SNP[record_before].append(record_name)
                            max_cluster_diff = max(max_cluster_diff, SNP_total)
                    SNP_total = SNP_seq(Ref, record_seq, POS_info, POS_info_CHR, POS_info_CHR_LEN,
                                        POS_info_output, REF_name, record_name)
                    SNP_pair.append('%s\t%s\t%s\t\n' % (REF_name, record_name, SNP_total))
                    if SNP_total <= SNP_total_cutoff:
                        Cluster_SNP.setdefault(REF_name,[])
                        Cluster_SNP[REF_name].append(record_name)
                        max_cluster_diff = max(max_cluster_diff,SNP_total)
                tree_distance.setdefault(record_name, SNP_total)
                Seq_list.setdefault(record_name, record_seq)
        cluster = 0
        # cluster genomes by SNP distance
        for record_name in Cluster_SNP:
            neighbor = Cluster_SNP.get(record_name,[])
            if neighbor != [] and record_name not in Cluster_SNP_set_added:
                cluster += 1
                Cluster_SNP_set.setdefault(cluster,set())
                Cluster_SNP_set[cluster].add(record_name)
                Cluster_SNP_set_added.add(record_name)
                find_neighbor(Cluster_SNP, neighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)
        # output single genome
        for record_name in Seq_list:
            if record_name not in Cluster_SNP_set_added:
                cluster += 1
                Cluster_SNP_set.setdefault(cluster, set())
                Cluster_SNP_set[cluster].add(record_name)
                Cluster_SNP_set_added.add(record_name)
        Sub_cluster = len(Cluster_SNP_set)
        for cluster in Cluster_SNP_set:
            tree_name_list = Cluster_SNP_set[cluster]
            tree_SNP_count = 'cluster%s' % (cluster)
            second_strain.setdefault('%s_%s'%(donor_species,tree_SNP_count),[])
            for tree_name in tree_name_list:
                second_strain['%s_%s'%(donor_species,tree_SNP_count)].append(tree_name)
                tree_distance_output.append(
                    '%s\t%s\t%s\t%s\t%s\t\n' % (
                        donor_species, tree_name, tree_distance[tree_name], tree_SNP_count, Sub_cluster))
        SNP_cutoff.append('%s\t%s\t%s\t%s\t\n' % (donor_species, max_cluster_diff, SNP_total_cutoff, SNP_total_length))
        if max_cluster_diff == 0:
            # no SNPs in a cluster
            for lines in open(fasta.replace('.all.flt.snp.vcf.filtered.fasta','.all.flt.snp.vcf.filtered.samplename.txt')):
                genomenames=lines.split('\n')[0].split('\t')
                for tree_name in genomenames:
                    tree_distance_output.append(
                        '%s\t%s\t%s\t%s\t%s\t\n' % (
                            donor_species, tree_name, 0, 'cluster1', 1))
                    for tree_name2 in genomenames:
                        if tree_name != tree_name2:
                            SNP_pair.append('%s\t%s\t%s\t\n' % (tree_name, tree_name2, 0))
        f1 = open(os.path.join(input_script, 'SNP_round%s.sum'%(Round)), 'a')
        f1.write('%s' % (''.join(tree_distance_output)))
        f1.close()
        f1 = open(os.path.join(input_script, 'SNP_round%s.allpair.sum' % (Round)), 'a')
        f1.write('%s' % (''.join(SNP_pair)))
        f1.close()
        f1 = open(os.path.join(input_script, 'SNP_round%s.sum.cutoff'%(Round)), 'a')
        f1.write('%s' % (''.join(SNP_cutoff)))
        f1.close()
        # f1 = open(os.path.join(input_script, 'SNP_round%s.POS.sum'%(Round)), 'w')
        # f1.write('%s' % (''.join(POS_info_output)))
        # f1.close()

# move
if Round < 3:
    cmd_move = ''
    fastq_name = '_1.fastq'
    fastq_name_2 = fastq_name.replace('1','2')
    for donor_species in second_strain:
        try:
            os.mkdir(genome_root + '/../round%s'%(Round + 1))
        except IOError:
            pass
        folders_dir = os.path.join(genome_root + '/../round%s'%(Round + 1), donor_species)
        try:
            os.mkdir(folders_dir)
        except IOError:
            pass
        try:
            os.mkdir(folders_dir + '/fastq')
        except IOError:
            pass
        donor_species_original = donor_species.split('_cluster')[0]
        for genome_files in second_strain[donor_species]:
            genome_files = genome_files.split('sortedbam')[0].split('fasta')[0].replace('af_Pseudoflavonifractor_sp_',
                                                                                        'af_Pseudoflavonifractor_sp._')
            cmd_move += ('mv %s/../round*/%s*/%s.* %s/\n'%(genome_root,donor_species_original,genome_files,folders_dir))
            cmd_move += ('mv %s/../round*/%s*/fastq/%s%s %s/fastq/\n' % (genome_root, donor_species_original,
                                                              genome_files, fastq_name, folders_dir))
            cmd_move += ('mv %s/../round*/%s*/fastq/%s%s %s/fastq/\n' % (genome_root, donor_species_original,
                                                              genome_files, fastq_name.replace('1','2'), folders_dir))
    f1 = open(os.path.join(input_script, 'SNP_round%s.move.sh'%(Round)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmd_move)))
    f1.close()

# compare genomes of different clusters
if Round == 1:
    genome_subsample = 2 # each cluster randomly pickup 2 genomes
    def run_vcf(genome_file,database,tempbamoutput,sumfile):
        # generate code
        # for curated genome
        cmds = 'minimap2 -d %s.mmi %s\n' % (database, database)
        cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
            min(40, 40), database, genome_file, 'samtools', min(40, 40),
            tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
            tempbamoutput)
        cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s.sorted.bam | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
            'bcftools', min(40, 40), database,
            tempbamoutput, 'bcftools', min(40, 40), tempbamoutput)
        cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
            'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
        cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
            'bcftools', tempbamoutput, tempbamoutput)
        cmds += 'wc -l %s.flt.snp.vcf >> %s\n'%(tempbamoutput,sumfile)
        cmds += 'rm -rf %s* %s.mmi\n' % (tempbamoutput,database)
        return cmds
    def pairgenome(allgenome, note,tag):
        allcmds = ''
        i = 0
        cluster_num = len(allgenome)
        for allpairs in itertools.combinations(range(0,cluster_num), 2):
            genomeset1 = allgenome[allpairs[0]]
            genomeset2 = allgenome[allpairs[1]]
            allgenomepairs = [[a, b] for a in genomeset1 for b in genomeset2 if a != b]
            for genomepair in allgenomepairs:
                genome1, genome2 = genomepair
                if '.all.' not in genome1 and '.all.' not in genome2:
                    tempbamoutput = os.path.join(input_script_temp,
                                                 '%s__%s' % (
                                                     os.path.split(genome1)[-1].split(fasta_name)[0],
                                                     os.path.split(genome2)[-1].split(fasta_name)[0]
                                                 ))
                    i += 1
                    allcmds += run_vcf(genome2, genome1, tempbamoutput,
                                       os.path.join(input_script_temp,
                                                    '%s.%s.sum.txt'%(tag,int(i / 100))))
                    if i % 100 == 0 and allcmds != '':
                        f1 = open(os.path.join(input_script_temp, '%s.%s.sh' % (tag,int(i / 100))), 'a')
                        f1.write('#!/bin/bash\nsource ~/.bashrc\n')
                        f1.write(''.join(allcmds))
                        f1.close()
                        allcmds = ''
        if allcmds!= '':
            print(note)
            f1 = open(os.path.join(input_script_temp, '%s.%s.sh' % (tag, int(i / 100))), 'a')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n')
            f1.write(''.join(allcmds))
            f1.close()
    def group_genome(genomelist):
        Preset = dict()
        Preset2 = dict()
        Cluster = dict()
        # input cluster
        for lines in open(os.path.join(input_script, 'SNP_round%s.all.sum'%(Round)), 'r'):
            lines_set = lines.split('\t')
            # genomename, subcluster
            Cluster.setdefault(lines_set[1],lines_set[3])
        # group genomes into subcluster
        for genome in genomelist:
            genomename = os.path.split(genome)[-1]
            donor = genomename.split('_')[0]
            species = '_'.join(genomename.split('_')[1:-1])
            cluster = Cluster.get(genomename, 'cluster1')
            donor_species = '%s_%s_%s' % (donor, species, cluster)
            if species in Species_replace:
                species = Species_replace[species]
            Preset.setdefault(species, set())
            Preset[species].add(donor_species)
            Preset2.setdefault(donor_species, set())
            Preset2[donor_species].add(genome)
        for setname in Preset:
            # onsidering genomes of diff donors
            if Preset2!= dict():
                allgenome2 = []
                num_genome = 0
                for donor_species in Preset[setname]:
                    allgenome = list(Preset2[donor_species])
                    # no need to pairing genomes in a donor species
                    # use SNP_round1.allpair.sum
                    # sub set genome_subsample genomes for a donor
                    allgenome2.append(random.choices(allgenome, k=genome_subsample))
                    num_genome += genome_subsample
                if len(Preset[setname]) > 1:
                    # at least 2 donors
                    note = ('pairing genome for %s %s clusters: %s genomes' % (
                        setname,len(Preset[setname]),num_genome))
                    print(Preset[setname])
                    print(note)
                    pairgenome(allgenome2, note, setname)
    os.system('#cat %s %s > %s'
              %(os.path.join(input_script, 'SNP_round%s.sum'%(Round)),
    os.path.join(input_script2, 'SNP_round%s.sum'%(Round)),
    os.path.join(input_script, 'SNP_round%s.all.sum'%(Round))))
    # run clustering for pairwise comparison below
    if Paircompare:
        # sum up SNP
        SNP_cutoff = []
        Fasta_SNP = glob.glob(os.path.join(output_dir_merge, '*.filtered.fasta'))
        f1 = open(os.path.join(input_script, 'SNP_round%s.%s.sum' % (Round, SNP_total_cutoff_2)), 'w')
        f1.write('donor_species\tGenome\tSNP_total\tcluster1\tsubcluster\t\n')
        f1.close()
        f1 = open(os.path.join(input_script, 'SNP_round%s.%s.allpair.sum' % (Round, SNP_total_cutoff_2)), 'w')
        f1.write('Genome1\tGenome2\tSNP_total\t\n')
        f1.close()
        for fasta in Fasta_SNP:
            SNP_pair = []
            POS_info_output = []
            tree_distance_output = []
            fasta_name = os.path.split(fasta)[-1]
            donor_species = fasta_name.split('.all.flt.snp.vcf.filtered.fasta')[0]
            POS_file = glob.glob(os.path.join(output_dir_merge, donor_species + '.all.flt.snp.vcf.filtered.POS.txt'))[0]
            print(donor_species)
            tree_distance = dict()
            Seq_list = dict()
            Ref = ''
            max_cluster_diff = 0
            POS_info = []
            POS_info_CHR = []
            POS_info_CHR_LEN = dict()
            Cluster_SNP = dict()
            Cluster_SNP_set = dict()
            Cluster_SNP_set_added = set()
            # load SNP POS info
            for lines in open(POS_file, 'r'):
                CHR = lines.split('\t')[0]
                POS_info.append(int(lines.split('\t')[1]))
                POS_info_CHR.append(CHR)
                POS_info_CHR_LEN.setdefault(CHR, 0)
                POS_info_CHR_LEN[CHR] = max(POS_info_CHR_LEN[CHR], int(lines.split('\t')[1]))
            # load genome SNP fasta and calculate pair-wise SNPs
            for record in SeqIO.parse(fasta, 'fasta'):
                record_name = str(record.id)
                if 'reference' not in record_name:
                    new_center = 1
                    record_seq = str(record.seq)
                    if Ref == '':
                        # set upt the first seq as ref
                        Ref = record_seq
                        REF_name = record_name
                        new_center = 0
                        SNP_total = 0
                        SNP_total_length = len(Ref)
                        SNP_total_cutoff = SNP_total_cutoff_2
                    else:
                        for record_before in Seq_list:
                            SNP_total = SNP_seq(Seq_list[record_before], record_seq, POS_info, POS_info_CHR,
                                                POS_info_CHR_LEN,
                                                POS_info_output, record_before, record_name)
                            SNP_pair.append('%s\t%s\t%s\t\n' % (record_before, record_name, SNP_total))
                            if SNP_total <= SNP_total_cutoff:
                                Cluster_SNP.setdefault(record_before, [])
                                Cluster_SNP[record_before].append(record_name)
                                max_cluster_diff = max(max_cluster_diff, SNP_total)
                        SNP_total = SNP_seq(Ref, record_seq, POS_info, POS_info_CHR, POS_info_CHR_LEN,
                                            POS_info_output, REF_name, record_name)
                        SNP_pair.append('%s\t%s\t%s\t\n' % (REF_name, record_name, SNP_total))
                        if SNP_total <= SNP_total_cutoff:
                            Cluster_SNP.setdefault(REF_name, [])
                            Cluster_SNP[REF_name].append(record_name)
                            max_cluster_diff = max(max_cluster_diff, SNP_total)
                    tree_distance.setdefault(record_name, SNP_total)
                    Seq_list.setdefault(record_name, record_seq)
            cluster = 0
            # cluster genomes by SNP distance
            for record_name in Cluster_SNP:
                neighbor = Cluster_SNP.get(record_name, [])
                if neighbor != [] and record_name not in Cluster_SNP_set_added:
                    cluster += 1
                    Cluster_SNP_set.setdefault(cluster, set())
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    find_neighbor(Cluster_SNP, neighbor, Cluster_SNP_set, cluster, Cluster_SNP_set_added)
            # output single genome
            for record_name in Seq_list:
                if record_name not in Cluster_SNP_set_added:
                    cluster += 1
                    Cluster_SNP_set.setdefault(cluster, set())
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
            Sub_cluster = len(Cluster_SNP_set)
            for cluster in Cluster_SNP_set:
                tree_name_list = Cluster_SNP_set[cluster]
                tree_SNP_count = 'cluster%s' % (cluster)
                for tree_name in tree_name_list:
                    tree_distance_output.append(
                        '%s\t%s\t%s\t%s\t%s\t\n' % (
                            donor_species, tree_name, tree_distance[tree_name], tree_SNP_count, Sub_cluster))
            SNP_cutoff.append(
                '%s\t%s\t%s\t%s\t\n' % (donor_species, max_cluster_diff, SNP_total_cutoff, SNP_total_length))
            if max_cluster_diff == 0:  # no SNPs in a cluster
                for lines in open(
                        fasta.replace('.all.flt.snp.vcf.filtered.fasta', '.all.flt.snp.vcf.filtered.samplename.txt')):
                    genomenames = lines.split('\n')[0].split('\t')
                    for tree_name in genomenames:
                        tree_distance_output.append(
                            '%s\t%s\t%s\t%s\t%s\t\n' % (
                                donor_species, tree_name, 0, 'cluster1', 1))
                        for tree_name2 in genomenames:
                            if tree_name != tree_name2:
                                SNP_pair.append('%s\t%s\t%s\t\n' % (tree_name, tree_name2, 0))
            f1 = open(os.path.join(input_script, 'SNP_round%s.%s.sum' % (Round, SNP_total_cutoff_2)), 'a')
            f1.write('%s' % (''.join(tree_distance_output)))
            f1.close()
            f1 = open(os.path.join(input_script, 'SNP_round%s.%s.allpair.sum' % (Round, SNP_total_cutoff_2)), 'a')
            f1.write('%s' % (''.join(SNP_pair)))
            f1.close()
            # f1 = open(os.path.join(input_script, 'SNP_round%s.%s.POS.sum'%(Round,SNP_total_cutoff_2)), 'w')
            # f1.write('%s' % (''.join(POS_info_output)))
            # f1.close()
        f1 = open(os.path.join(input_script, 'SNP_round%s.%s.sum.cutoff' % (Round, SNP_total_cutoff_2)), 'w')
        f1.write('donor_species\tmax_cluster_diff\tSNP_cutoff\tSNP_total_len\t\n%s' % (''.join(SNP_cutoff)))
        f1.close()
    if Paircompare:
        os.system('rm -rf %s'%(input_script_temp))
        try:
            os.mkdir(input_script_temp)
        except IOError:
            pass
        allgenome = glob.glob('%s/*/*%s'%(genome_root, fasta_name))+\
        glob.glob('%s/*/*%s' % (genome_root2, fasta_name))
        group_genome(allgenome)
        f1 = open(os.path.join(input_script, 'allpaircmds.sh'), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n')
        for sub_scripts in glob.glob(os.path.join(input_script_temp, '*.sh')):
            f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
        f1.close()
    # after finished summarize inter donor compare
    os.system('cat %s > %s'
              %(os.path.join(input_script_temp, '*.sum.txt'),
                  os.path.join(input_script, 'SNP_round%s.allpair.all.diffcluster.sum'%(Round)))
              )
    os.system('#cat %s %s > %s'
              %(os.path.join(input_script2, 'SNP_round%s.allpair.sum'%(Round)),
                  os.path.join(input_script, 'SNP_round%s.allpair.sum'%(Round)),
                os.path.join(input_script, 'SNP_round%s.allpair.all.sum' % (Round)))
              )
    def annogenome(G1):
        G1_set = G1.split('_')
        donor = G1_set[0]
        species = '_'.join(G1_set[1:-1])
        if species in Species_replace:
            species = Species_replace[species]
        return [donor,species]
    def diff(item1, item2):
        if item1 == item2:
            return 'same'
        else:
            return 'diff'
    def splitgepair(G_pair):
        genomename = os.path.split(G_pair)[-1].split('.flt.snp.vcf')[0]
        G1,G2 = genomename.split('__')
        donor1, species1 = annogenome(G1)
        donor2, species2 = annogenome(G2)
        donor_anno = '%s_donor'%(diff(donor1, donor2))
        if diff(species1, species2)=='diff':
            # check whether 2 genomes are of the same species
            print('wrong species pair: %s %s'%(G1,G2))
        if diff(donor1, donor2)=='sane':
            # check whether 2 genomes are of diff donors
            print('wrong donor pair: %s %s'%(G1,G2))
        return '\t'.join([species1,
            G1,donor1,
            G2,donor2])
    Sumoutput = []
    Sumoutput.append('SNP_count\tspecies\tgenome1\tdonor1\tgenome2\tdonor2\n')
    for lines in open(os.path.join(input_script, 'SNP_round%s.allpair.all.diffcluster.sum'%(Round)),'r'):
        lines_set = lines.split(' ')
        SNP_count, G_pair = lines_set[0:2]
        Sumoutput.append('%s\t%s\n'%(SNP_count, splitgepair(G_pair)))
    f1 = open(os.path.join(input_script, 'SNP_round%s.allpair.all.diffcluster.sumnew.txt'%(Round)),'w')
    f1.write(''.join(Sumoutput))
    f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
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
    cmds += 'bowtie2' + ' --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
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
################################################### SET PATH ########################################################
# round 1-3 run genome assembly and map genomes to a reference genome
import glob
import os
Round = 4
input_script_vcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_round%s'%(Round)
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
genome_dir = glob.glob('%s/round%s/*'%(genome_root,Round))
if Round == 4:
    genome_dir = glob.glob('%s/round*/*' % (genome_root))

output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s'%(Round)
output_dir_old = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s'%(Round-1)
fastq_name = '_1.fastq'
genome_name = '.fasta.corrected.fasta'

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
    os.mkdir(input_script)
except IOError:
    pass

os.system('rm -rf %s'%(input_script_vcf))

try:
    os.mkdir(input_script_vcf)
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

def run_vcf_WGS(files,files2,database,tempbamoutput):
    # generate code
    cmds = 'rm -rf %s.fai\n' % (database)
    cmds += 'bowtie2' + ' --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, files, files2, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    cmds += 'rm -rf %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    sample_set3 = ['%s.sorted.bam' % (tempbamoutput)]
    return [cmds, sample_set3]

def run_vcf(genome_file,database,tempbamoutput):
    # generate code
    # for curated genome
    cmds = 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, genome_file, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    cmds += '#samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
        tempbamoutput, tempbamoutput)
    sample_set3 = ['%s.sorted.bam' % (tempbamoutput)]
    cmds += 'rm -rf %s.bam %s.bam.bai\n' %(tempbamoutput,tempbamoutput)
    return [cmds,sample_set3]

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
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    if len(donor) == 2:
        # BN10 donors
        donor_species_genome = glob.glob(os.path.join(folder,genome_name))
        sub_samples = []
        donor_species_fastqall = glob.glob(os.path.join(folder, 'fastq/*' + fastq_name)) +\
                              glob.glob(os.path.join(folder, '*' + fastq_name))
        cmds = ''
        cmds2 = ''
        if len(donor_species_fastqall) >= 3:
            donor_species_folder_all = os.path.join(folder, donor_species + '_allspades%s'%(Round))
            donor_species_genomename = os.path.join(folder, donor_species + '.all.spades%s.fasta'%(Round))
            donor_species_fastq = os.path.join(folder, donor_species + '.all.spades%s'%(Round) + fastq_name)
            donor_species_fastq2 = os.path.join(folder, donor_species + '.all.spades%s'%(Round) + fastq_name.replace('1', '2'))
            tempbamoutputGenome = os.path.join(output_dir + '/merge', donor_species + '.all')
            samplesetGenomecorrect = []
            if Round == 4:
                donor_species_genomename = glob.glob(os.path.join(folder, donor_species + '.all.spades*.fasta'))
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
                print(folder, tempbamoutputGenome)
            for fastq_file in donor_species_fastqall:
                if '.all' + fastq_name not in fastq_file and '.all.spades' not in fastq_file:
                    fastq_file_name = os.path.split(fastq_file)[-1]
                    filename = fastq_file.split(fastq_name)[0]
                    fastq_file2 = filename + fastq_name.replace('1', '2')
                    genome_file = os.path.join(folder, fastq_file_name.split(fastq_name)[0] + genome_name)
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
            if Round <= 3:
                # pan genome
                cmds += runspades(donor_species_fastq, donor_species_fastq2, donor_species_folder_all, donor_species_genomename)
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
################################################### SET PATH ########################################################
# CHANGED After round 4 filter results of WGS
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import itertools
import random
# set up path
Cov_dis = 20
Round = 4
end_cutoff = 70 # 10 bp at the ends of a contig, separate

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round%s_tree'%(Round)
input_script_sub_merge = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/vcf_round%s'%(Round)
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay'
genome_root = '/scratch/users/anniz44/genomes/donor_species/jay/round*/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/round%s/*'%(Round))
output_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/bwa/0/'%(Round)
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/merge'%(Round)
vcf_name = '.all.flt.snp.vcf'
ref_filename = '.all.spades%s.fasta'%(Round)
fasta_name = '.fasta.corrected.fasta'
fastq_name = '.sorted.bam'
deleting_file = []


input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_round%s_tree'%(Round)
input_script_sub_merge = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_round%s'%(Round)
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
input_script2 = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round*/'
genome_root2 = '/scratch/users/anniz44/genomes/donor_species/jay'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round%s/*'%(Round))
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/bwa/0/'%(Round)
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/merge'%(Round)
vcf_name = '.all.flt.snp.vcf'
ref_filename = '.all.spades%s.fasta'%(Round)
fasta_name = '.fasta.corrected.fasta'
deleting_file = []
fastq_name = '_1.fastq'

################################################### Function ########################################################
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def contig_end(CHR,POS):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        total_length = CHR.split('length_')[1].split('_cov')[0]
    total_length = int(total_length)
    if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
        return True
    else:
        return False

def vcf_to_txt(lines,output_list,cluster_sub=[]):
    lines_set = lines.split('\n')[0].split('\t')
    if len(lines_set) >9:
        CHR = lines_set[0]
        POS = int(lines_set[1])
        temp_line = []
        temp_line.append(CHR)
        temp_line.append(str(POS))
        i = 9
        for Subdepth_all in lines_set[9:]:
            if (cluster_sub==[] and i not in deleting_set) or i in cluster_sub:
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                temp_line.append(str(total_sub_depth))
            i += 1
        output_list.append('\t'.join(temp_line)+'\n')
    else:
        print(lines)

def curate_REF(allels_set,lines_set_sub):
    SNP_count_genome_count = [0,0,0,0]
    for Subdepth_all in lines_set_sub:
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '')
        Subdepth_set = Subdepth.split(',')
        Subdepth_set_int = []
        for sub_depth in Subdepth_set:
            Subdepth_set_int.append(int(sub_depth))
        if sum(Subdepth_set_int) > 0:
            major_alt_frq = max(Subdepth_set_int)
            major_alt_frq_index = Subdepth_set_int.index(major_alt_frq)
            major_alt = allels_set[major_alt_frq_index]
            SNP_count_genome_count[major_alt_frq_index] += 1
    major_alt_frq = max(SNP_count_genome_count)
    major_alt_frq_index = SNP_count_genome_count.index(major_alt_frq)
    major_alt = allels_set[major_alt_frq_index]
    return [major_alt, major_alt_frq_index]

def outputvcf(output_name):
    vcf_file_filtered = open(vcf_file + '.%s.vcf' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list_vcf))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.snp.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.samplename.txt' % (output_name), 'w')
    vcf_file_filtered.write('\t'.join(Sample_name) + '\n')
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.POS.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_POS))
    vcf_file_filtered.close()

def outputcov(output_name,vcf_file_POS_candidate,cluster_sub=[]):
    if len(vcf_file_list) > 0:
        vcf_file_POS_candidate = '\n'.join(vcf_file_POS_candidate)
        vcf_file_POS_candidate_output = ('%s' % (vcf_file_POS_candidate))
        f1 = open(os.path.join(input_script,'grep.temp.txt'),'w')
        f1.write(vcf_file_POS_candidate_output)
        f1.close()
        os.system('grep -%s -f %s %s --no-group-separator > %s'% (
            Cov_dis, os.path.join(input_script,'grep.temp.txt'),
            vcf_file.split('.all.fq.flt.snp.vcf')[0] + '.all.fq.raw.vcf',
            vcf_file + '.%s.cov.temp' % (output_name)))
        os.system('cat %s | sort | uniq > %s' % (
            vcf_file + '.%s.cov.temp' % (output_name),
            vcf_file + '.%s.uniqcov.temp' % (output_name)))
        for lines in open(vcf_file + '.%s.uniqcov.temp' % (output_name), 'r'):
            if not lines.startswith("#"):
                vcf_to_txt(lines, cov_file_list,cluster_sub)
        os.system('rm -rf %s %s' % (vcf_file + '.%s.cov.temp' % (output_name),
                                    vcf_file + '.%s.uniqcov.temp' % (output_name)))
        vcf_file_filtered = open(vcf_file + '.%s.cov.txt' % (output_name), 'w')
        vcf_file_filtered.write(''.join(cov_file_list))
        vcf_file_filtered.close()

def outputtree(output_name):
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, SNP_alignment[genomename]))
            SNP_alignment_output_parsi.append('S%s    %s\n' % (genomename[-8:], SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.%s.fasta' % (output_name), 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.parsi.fasta' %(output_name), 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()

def SNP_seq(seq1, seq2, POS_info,POS_info_CHR,POS_info_CHR_LEN,POS_info_output,G1,G2):
    SNP_total = 0
    j = 0
    POS_DIS = []
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            CHR = POS_info_CHR[i]
            POS = POS_info[i]
            LEN = POS_info_CHR_LEN[CHR]
            if CHR == POS_info_CHR[j]:  # same CHR
                DIS = abs(POS - POS_info[j])
                POS_DIS.append(DIS)  # POS diff
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, DIS, LEN))
            else:  # new CHR
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, 0, LEN))
            j = i
    return SNP_total

def SNP_distance_correct(distanceNEW, distanceREF, SNPREF):
    return int((distanceNEW/distanceREF)*SNPREF)

def find_neighbor(Cluster_SNP,neighbor,Cluster_SNP_set,cluster,Cluster_SNP_set_added):
    if neighbor != []:
        for record_name in neighbor:
            if record_name not in Cluster_SNP_set_added:
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    subneighbor = Cluster_SNP.get(record_name,[])
                    find_neighbor(Cluster_SNP, subneighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)

def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def dnORds(amino1, amino2):
    if amino1 == amino2:
        return 'S'
    else:
        return 'N'

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

def contig_to_gene(CHR, POS):
    all_genes = Mapping_loci.get(CHR,[])
    Reverse_chr = 0
    for a_gene in all_genes:
        POS1, POS2, GENE = a_gene
        if POS >= POS1 and POS <= POS2:
            Ref_seq_chr = Ref_seq.get(GENE, 'None')
            Gene_length = len(Ref_seq_chr)
            if GENE in Reverse:  # reversed
                POS_gene = Gene_length-(int(POS-POS1))
                Reverse_chr = 1
            else:
                POS_gene = int(POS-POS1)+1
            codon_start = POS_gene - 1 - int((POS_gene - 1) % 3)
            return [GENE,POS_gene,codon_start,Ref_seq_chr,Reverse_chr]
    return []

def SNP_check_all(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name,SNP_presence_cutoff,SNP_presence_sample_cutoff,cluster_sub=[]):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    Total_qualify = 0
    Total_qualify_SNP = 0
    Total_qualify_notSNP = 0
    Total_unqualify_alt_freq = 0
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    if cluster_sub!= []:
        lines_set_sub = [lines_set[i] for i in cluster_sub]
        Total_subsample = len(cluster_sub)
        if Total_subsample >= 15:
            SNP_presence_cutoff = 0.33  # for a large group of samples
        if Total_subsample <= 3:
            SNP_presence_cutoff = 1  # for a small group of samples
            SNP_presence_sample_cutoff = 2
    else:
        cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        REF, REF_where = curate_REF(allels_set, lines_set_sub)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
                Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
                total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
                total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                        Allels_frq_sub[Allels[allels] * 2] += int(Subdepth_forward[num_allels])
                        Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_reverse[num_allels])
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    MLF = Major_ALT[1] / total_sub_depth
                    if total_sub_depth_forward >= Sample_depth_cutoff and \
                            total_sub_depth_reverse >= Sample_depth_cutoff:
                        # forward and reverse cutoff POS detected
                        if MLF >= Major_alt_freq_cutoff:
                        # major alt frequency cutoff
                            Total_qualify += 1
                            # check for qualified SNP
                            if Major_ALT[0] != REF:
                                SNP_seq[-1] = Major_ALT[0]  # unqualified SNP also include in alignment
                                # qualified SNP
                                if (total_sub_depth_forward - int(Subdepth_forward[REF_where]) >= SNP_depth_cutoff and \
                                        total_sub_depth_reverse - int(Subdepth_reverse[REF_where]) >= SNP_depth_cutoff) and \
                                        Major_ALT[1] >= SNP_depth_cutoff:
                                    Total_qualify_SNP += 1
                                    SNP.add(genome_order)  # only take qualified SNP as valid SNP
                                    SNP_seq[-1] = Major_ALT[0]
                            else:
                                Total_qualify_notSNP += 1
                        else:
                            # major alt frequency low
                            Total_unqualify_alt_freq += 1
            sample_num += 1
        if Total_qualify / Total_subsample >= SNP_presence_cutoff and \
                Total_unqualify_alt_freq / Total_subsample <= Poor_MLF_freq_cutoff and\
                Total_qualify >= SNP_presence_sample_cutoff and \
                Total_qualify_SNP >= 1 and Total_qualify_SNP < Total_qualify and\
                Total_qualify_notSNP > 0:
            # -> qualified SNP
            # qualified samples cutoff + unqualified samples cutoff
            # at least 1 qualified SNP
            # calculate NS
            gene_info = contig_to_gene(CHR, POS)
            if gene_info!= []:
                Chr_gene, POS_gene,codon_start,Ref_seq_chr,Reverse_chr  = gene_info
                if Ref_seq_chr != 'None':
                    #  observed NS ratio calculated
                    temp_snp_line_NS= [Chr_gene,str(POS_gene),'']
                    if codon_start <= POS_gene - 1:
                        Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                        Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                        SNP_seq_chr = Ref_seq_chr
                        if len(Ref_seq_codon) == 3:
                            Ref_seq_aa = translate(Ref_seq_codon)[0]
                            temp_snp_line_AA += Ref_seq_aa
                            ALT_set = allels_set
                            ALT_set.remove(REF)
                            for ALT in ALT_set:
                                SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                SNP_seq_aa = translate(SNP_seq_codon)[0]
                                temp_snp_line_AA += SNP_seq_aa
                                temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                temp_snp_line_NS[-1]+=temp_NorS
            # output lines and output major alt
            temp_snp_line_pass = 'PASS'
            if CHR == CHR_old:
                # same CHR
                POS_DIS = abs(POS - POS_old)
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
            else:
                # diff CHR first SNP
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
            POS_old = POS
            CHR_old = CHR
            temp_snp_line.append(CHR)
            temp_snp_line.append(str(POS))
            temp_snp_line.append(REF)
            temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
            vcf_file_list.append('\t'.join(temp_snp_line)+ '\t' +'\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,'\t'.join(temp_snp_line_NS),temp_snp_line_AA))
            vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
            vcf_file_list_vcf.append('\t'.join(lines_set[0:9])+'\t'+'\t'.join(lines_set_sub)+'\n')
            i = 9
            j = 0
            SNP_alignment[reference_name] += REF
            for genomename in SNP_alignment:
                if genomename != reference_name:
                    if i in cluster_sub:
                        SNP_alignment[genomename] += SNP_seq[j]
                        j += 1
                    i += 1
    return [CHR_old,POS_old]

def curate_REF_fq(allels_set,lines_set_sub):
    SNP_count_genome_count = [0,0,0,0]
    for Subdepth_all in lines_set_sub:
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '')
        Subdepth_set = Subdepth.split(',')
        Subdepth_set_int = []
        for sub_depth in Subdepth_set:
            Subdepth_set_int.append(int(sub_depth))
        total_sub_depth = sum(Subdepth_set_int)
        if total_sub_depth > 0:
            major_alt_frq = max(Subdepth_set_int)
            major_alt_frq_index = Subdepth_set_int.index(major_alt_frq)
            SNP_count_genome_count[major_alt_frq_index] += 1
    major_alt_frq = max(SNP_count_genome_count)
    major_alt_frq_index = SNP_count_genome_count.index(major_alt_frq)
    major_alt = allels_set[major_alt_frq_index]
    return [major_alt, major_alt_frq_index]

def SNP_check_all_fq(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    Total_qualify = 0
    Total_qualify_SNP = 0
    Total_qualify_notSNP = 0
    Total_unqualify_alt_freq = 0
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        genome_order = 0
        REF, REF_where = curate_REF(allels_set, lines_set_sub)  # as the major alt in the population
        sample_num = 9
        Total_alleles = len(allels_set)
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
                Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
                total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
                total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                        Allels_frq_sub[Allels[allels] * 2] += int(Subdepth_forward[num_allels])
                        Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_reverse[num_allels])
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    MLF = Major_ALT[1] / total_sub_depth
                    if Major_ALT[0] != REF:
                        SNP.add(genome_order)  # only take qualified SNP as valid SNP
                        SNP_seq[-1] = Major_ALT[0]
        # calculate NS
        gene_info = contig_to_gene(CHR, POS)
        if gene_info != []:
            Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
            if Ref_seq_chr != 'None':
                #  observed NS ratio calculated
                temp_snp_line_NS = [Chr_gene, str(POS_gene), '']
                if codon_start <= POS_gene - 1:
                    Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                    Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                    SNP_seq_chr = Ref_seq_chr
                    if len(Ref_seq_codon) == 3:
                        Ref_seq_aa = translate(Ref_seq_codon)[0]
                        temp_snp_line_AA += Ref_seq_aa
                        ALT_set = allels_set
                        ALT_set.remove(REF)
                        for ALT in ALT_set:
                            SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                            SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                            SNP_seq_aa = translate(SNP_seq_codon)[0]
                            temp_snp_line_AA += SNP_seq_aa
                            temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                            temp_snp_line_NS[-1] += temp_NorS
        # output lines and output major alt
        temp_snp_line_pass = 'PASS'
        if CHR == CHR_old:
            # same CHR
            POS_DIS = abs(POS - POS_old)
            vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
        else:
            # diff CHR first SNP
            vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
        POS_old = POS
        CHR_old = CHR
        temp_snp_line.append(CHR)
        temp_snp_line.append(str(POS))
        temp_snp_line.append(REF)
        temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
        vcf_file_list.append(
            '\t'.join(temp_snp_line) + '\t' + '\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass, '\t'.join(temp_snp_line_NS),
                temp_snp_line_AA))
        vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
        vcf_file_list_vcf.append('\t'.join(lines_set[0:9]) + '\t' + '\t'.join(lines_set_sub) + '\n')
        i = 9
        j = 0
        SNP_alignment[reference_name] += REF
        for genomename in SNP_alignment:
            if genomename != reference_name:
                if i in cluster_sub:
                    SNP_alignment[genomename] += SNP_seq[j]
                    j += 1
                i += 1
    return [CHR_old,POS_old]

def load_ref_vcf(ref_vcf_file):
    ref_chr = dict()
    for files in ref_vcf_file:
        for lines in open(files,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS, Notused, REF, ALT = lines_set[0:5]
            CHR_POS = '%s__%s'%(CHR, POS)
            ref_chr.setdefault(CHR_POS,[])
            ref_chr[CHR_POS]=[REF,ALT]
    return ref_chr

################################################### Set up ########################################################
# set up steps
SNP_cluster = dict()
cluster_set = set()
reference_set = ['reference']
outputname_set = ['filtered']

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

################################################### Main ########################################################
# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir_merge,'*.flt.snp.vcf'))
for set_num in range(0,len(reference_set)):
    reference_name = reference_set[set_num]
    output_name = outputname_set[set_num]
    for vcf_file in all_vcf_file:
        Total = 0
        donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
        SNP_tree_cmd = []
        SNP_tree_cmd2 = []
        vcf_file_list = []
        vcf_file_list_vcf = []
        Sample_name = []
        deleting_set = []
        vcf_file_POS = []
        vcf_file_POS_candidate = set()
        SNP_alignment = dict()
        SNP_alignment.setdefault(reference_name, '')
        cov_file_list = []
        Ref_seq = dict()
        Mapping = dict()
        Mapping_loci = dict()
        CHR_old = ''
        POS_old = 0
        SNP_cluster_donor_species = dict()
        try:
            # genome mapping file
            vcf_ref_file_name = os.path.split(vcf_file)[-1].split('.all.')[0]
            # need change later
            vcf_genome = glob.glob(output_dir_merge + '/../merge_genome/vcf/' + vcf_ref_file_name + '*.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf.singleevent.vcf')
            if vcf_genome == []:
                vcf_ref_file_name = '_cluster'.join(vcf_ref_file_name.split('_cluster')[:-1])
                vcf_genome = glob.glob(
                    output_dir_merge + '/../merge_genome/vcf/' + vcf_ref_file_name + '*.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf.singleevent.vcf')
            ref_chr = load_ref_vcf(vcf_genome)
            for cluster_type in cluster_set:
                SNP_cluster_donor_species.setdefault(cluster_type,[])
            for lines in open(os.path.join(input_script_sub_merge, '%s.vcf.sh' % (donor_species)), 'r'):
                if lines.startswith('bcftools mpileup '):
                    # setup samples
                    sample_set = lines.split('.fasta ')[1].split('\n')[0].split(' | ')[0].split(' ')
                    samplenum = 9
                    for samples in sample_set:
                        genomename = os.path.split(samples)[-1].split(fastq_name)[0].split('all')[0].split('.sorted.bam')[0]
                        Sample_name.append(genomename.replace('.', ''))
                        if genomename in deleting_file:
                            deleting_set.append(samplenum)
                        else:
                            SNP_alignment.setdefault(genomename, '')
                            if SNP_cluster!= dict() and genomename in SNP_cluster:
                                SNP_cluster_donor_species[SNP_cluster[genomename]].append(samplenum)
                        samplenum += 1
            print('running %s' % (donor_species))
            print(Sample_name)
            # load database
            database_file = glob.glob(os.path.join(genome_root + '/%s/'%(vcf_ref_file_name),'*all.spades*.fna' ))
            if database_file != []:
                Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file[0])
            for lines in open(vcf_file, 'r'):
                if not lines.startswith("#"):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR = lines_set[0]
                    POS = int(lines_set[1])
                    # a SNP confirmed in WGS mapping
                    CHR_POS = '%s__%s' % (CHR, POS)
                    if CHR_POS in ref_chr and not contig_end(CHR, POS):
                        Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                        if Total == 0:
                            Total = len(lines_set) - 9 - len(deleting_set)
                        if "INDEL" not in lines_set[7] \
                                and (lines_set[6] != 'LowQual'):
                            CHR_old, POS_old = SNP_check_all_fq(lines_set, '',
                                                                CHR_old, POS_old, reference_name)
            outputvcf(output_name)
            try:
                f1 = open(vcf_file + '.%s.cov.txt' % (output_name), 'r')
            except IOError:
                outputcov(output_name, list(vcf_file_POS_candidate), [])
        except IndexError:
            print('missing genome mapping file %s'%(output_dir_merge + '/../merge_genome/vcf/' + vcf_ref_file_name + '*.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf'))

################################################## END ########################################################
################################################### SET PATH ########################################################
# compare genome call SNPs VS WGS call SNPs
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import random
import argparse

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="subsample * 20", type=int,
                    default=0,metavar='0 to 11')
################################################## Definition ########################################################
args = parser.parse_args()
i = int(args.i)
# set up path
output_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round4/merge/'
ref_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round4/merge_genome/'
# genome
vcf_ref = glob.glob(ref_dir + '/*.all.flt.snp.vcf.filtered.snp.txt')
# set up cutoff
end_cutoff = 70 # 10 bp at the ends of a contig, separate

# function
def load_vcf_ref(vcf_file,end_count = False):
    end_set = ['end', 'notend']
    vcf_count = dict()
    vcf_count.setdefault('end', 0)
    vcf_count.setdefault('notend', 0)
    vcf_input = []
    vcf_ref = dict()
    for lines in open(vcf_file):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = lines_set[1]
            CHRPOS = '%s\t%s\t\n'%(CHR,POS)
            Gene = lines_set[-2]
            if Gene == 'None':
                Gene = Gene.replace('None','Other')
            else:
                Gene = 'Gene'
            vcf_input.append(CHRPOS)
            vcf_ref.setdefault(CHRPOS,Gene)
            if end_count:
                vcf_count[contig_end(CHR,POS)] += 1
    if end_count:
        summary_file_output.append('%s\t%s\t%s\t0\t0\n'%(os.path.split(vcf_file)[-1],vcf_count['end'],
                                                               vcf_count['notend']))
        print(vcf_count)
    return [vcf_input,vcf_ref]

def load_vcf(vcf_file,end_count = False):
    end_set = ['end', 'notend']
    vcf_count = dict()
    vcf_count.setdefault('end', 0)
    vcf_count.setdefault('notend', 0)
    vcf_input = []
    for lines in open(vcf_file):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = lines_set[1]
            vcf_input.append('%s\t%s\t\n'%(CHR,POS))
            if end_count:
                vcf_count[contig_end(CHR,POS)] += 1
    if end_count:
        summary_file_output.append('%s\t%s\t%s\t0\t0\n'%(os.path.split(vcf_file)[-1],vcf_count['end'],
                                                               vcf_count['notend']))
        print(vcf_count)
    return vcf_input

def contig_end(CHR,POS):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        total_length = CHR.split('length_')[1].split('_cov')[0]
    total_length = int(total_length)
    if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
        return 'end'
    else:
        return 'notend'

def compare_vcf(vcf_input1, vcf_input2, vcf_file1, vcf_file2, output_file):
    #end_set = ['end','notend']
    end_set = ['notend']
    vcf_1_diff = dict()
    vcf_1_diff.setdefault('end',[])
    vcf_1_diff.setdefault('notend', [])
    vcf_1_diff.setdefault('Gene', 0)
    vcf_1_diff.setdefault('Other', 0)
    # CHRPOS in 1 not in 2
    for CHRPOS in vcf_input1:
        if CHRPOS not in vcf_input2:
            CHR,POS = CHRPOS.split('\t')[0:2]
            contig_end_tag = contig_end(CHR,POS)
            vcf_1_diff[contig_end_tag].append(CHRPOS)
            if contig_end_tag == 'notend' and CHRPOS in ref_gene:
                vcf_1_diff[ref_gene[CHRPOS]] += 1
    for contig_end_tag in end_set:
        if len(vcf_1_diff[contig_end_tag]) > 0:
            print(len(vcf_1_diff[contig_end_tag]))
            temp_output = os.path.join(output_dir, 'grep.temp.%s.txt' % (contig_end_tag))
            print(temp_output)
            f1 = open(temp_output, 'w')
            f1.write(''.join(vcf_1_diff[contig_end_tag]))
            f1.close()
            #os.system('grep -T -f %s %s %s --no-group-separator > %s' % (
            #    temp_output,
            #    vcf_file1,vcf_file2,
            #    output_file + '.temp'))
            #os.system('sort -k3 -n %s | sort -k2 > %s' %
            #          (output_file + '.temp', output_file + '.' + contig_end_tag)
            #          )
            #os.system('rm -rf %s' % (output_file+ '.temp'))
        else:
            os.system('rm -rf %s'%(output_file + '.' + contig_end_tag))
    return vcf_1_diff

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def compare_vcf_all(vcf_1,vcf_2 = ''):
    vcf_input1 = load_vcf(vcf_1)
    vcf_1_all = glob.glob(vcf_1.split('.flt.snp.vcf')[0] + '.raw.vcf')
    if vcf_1_all != []:
        vcf_1_all = vcf_1_all[0]
    else:
        vcf_1_all = vcf_1
    if vcf_2 != '':
        vcf_input1 = intersection(vcf_input1, load_vcf(vcf_2))
        vcf_input1 = list(set(vcf_input1))
        vcf_1 = vcf_1 + '.all'
    FN_diff = compare_vcf(vcf_inputref, vcf_input1, vcf_ref_file, vcf_1_all, vcf_1 + '.FN') # vcf diff in ref not in 1, FN
    FP_diff = compare_vcf(vcf_input1, vcf_inputref, vcf_1, vcf_ref_file, vcf_1 + '.FP') # vcf diff in 1 not in ref, FP
    summary_file_output.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(
        os.path.split(vcf_1)[-1],
        len(FN_diff['end']),len(FN_diff['notend']),
        len(FP_diff['end']), len(FP_diff['notend']),
                                 FN_diff['Gene'],FN_diff['Other']
    ))

def depthcheck(vcf_1,vcf_2,vcf_3):
    vcf_input2 = load_vcf(vcf_2)
    vcf_input3 = load_vcf(vcf_3)
    Length = dict()
    Total = 0
    for lines in open(vcf_1):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
            if Total == 0:
                Total = len(lines_set) - 9
            if Depth / Total <= 100:
                CHR = lines_set[0]
                if CHR not in Length:
                    try:
                        total_length = CHR.split('size')[1]
                    except IndexError:
                        total_length = CHR.split('length_')[1].split('_cov')[0]
                    total_length = int(total_length)
                    Length.setdefault(CHR,total_length)
                total_length = Length[CHR]
                if total_length >= 5000:
                    POS = lines_set[1]
                    CHRPOS = '%s\t%s\t\n' % (CHR, POS)
                    lines_set_sub = lines_set[9:]
                    for Subdepth_all in lines_set_sub:
                        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                        total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                        # Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
                        # Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
                        # forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
                        # reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
                        if total_sub_depth > 0:
                            Depth_set.setdefault(total_sub_depth, [0, 0, 0, 0])
                            if CHRPOS in vcf_inputref:
                                Depth_set[total_sub_depth][1] += 1
                                if CHRPOS in vcf_input3:
                                    Depth_set[total_sub_depth][3] += 1
                            else:
                                Depth_set[total_sub_depth][0] += 1
                            if CHRPOS in vcf_input2:
                                Depth_set[total_sub_depth][2] += 1
    return 'done'

# WGS VS Genome
summary_file = output_dir + '/diff.sum.nofilter.txt'
summary_file_output = []
summary_file_output.append('sample\tFN_end\tFN_notend\tFP_end\tFP_notend\tFN_notend_gene\tFN_notend_other\n')
for vcf_ref_file in vcf_ref:
    vcf_ref_file_name = os.path.split(vcf_ref_file)[-1].split('.all.')[0]
    try:
        vcf_inputref, ref_gene = load_vcf_ref(vcf_ref_file,True)
        # WGS filtered
        vcf_1 = glob.glob(output_dir + vcf_ref_file_name + '*.all.fq.flt.snp.vcf.filtered.snp.txt')[0]
        compare_vcf_all(vcf_1)
        # WGS not filtered
        vcf_1 = glob.glob(output_dir + vcf_ref_file_name + '*.all.fq.flt.snp.vcf')[0]
        compare_vcf_all(vcf_1)
    except IndexError:
        pass

f1=open(summary_file,'w')
f1.write(''.join(summary_file_output))
f1.close()
os.system('rm -rf %s'%(os.path.join(output_dir, 'grep.temp.*.txt')))

# Depth check of Genome mutations
summary_file = output_dir + '/SNP.depth.sum.nofilter%s.txt' %(i)
f1=open(summary_file,'w')
f1.write('Depth\tGe_noSNP\tGe_SNP\tWGS_filterSNP\tWGS_Ge_SNP\n')
f1.close()
for vcf_ref_file in vcf_ref[i*20:(i+1)*20]:
    vcf_ref_file_name = os.path.split(vcf_ref_file)[-1].split('.all.')[0]
    Depth_set = dict()
    summary_file_output = []
    try:
        # WGS
        vcf_1 = glob.glob(output_dir + vcf_ref_file_name + '*.all.fq.raw.vcf')[0]
        vcf_2 = glob.glob(output_dir + vcf_ref_file_name + '*.all.fq.flt.snp.vcf.filtered.snp.txt')[0]
        vcf_3 = glob.glob(output_dir + vcf_ref_file_name + '*.all.fq.flt.snp.vcf')[0]
        vcf_inputref = load_vcf(vcf_ref_file,True)
        depthcheck(vcf_1, vcf_2, vcf_3)
        summary_file_output = []
        for Depth in Depth_set:
            summary_file_output.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (vcf_ref_file_name,Depth,
                                                                 Depth_set[Depth][0], Depth_set[Depth][1],
                                                                 Depth_set[Depth][2], Depth_set[Depth][3]))

        f1 = open(summary_file, 'a')
        f1.write(''.join(summary_file_output))
        f1.close()
        print(Depth_set.get(3, 'None'), Depth_set.get(10, 'None'))
    except IndexError:
        pass

################################################### END ########################################################
################################################### SET PATH ########################################################
# after round 4 set cutoff for recombination windows to calculate NS ratio
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics

input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge'
vcf_name = '.raw.vcf.filtered.snp.txt'
windowsize_set = [810000, 729000, 656100, 590490, 531441, 478296, 430466, 387419, 348677, 313809,
                  282428, 254185, 228766, 205889, 185300, 166770, 150093, 135083, 121574, 109416,
                  98474, 88626, 79763, 71786, 64607, 58146, 52331, 47097, 42387, 38148, 34333, 30899,
                  27809, 25028, 22525, 20272, 18244, 16419, 14777, 13299, 11969, 10772, 9694, 8724,
                  7851, 7065, 6358, 5722, 5149, 4634, 4170, 3753, 3377, 3039, 2735, 2461, 2214, 1992,
                  1792, 1612, 1450, 1305, 1174, 1056, 950, 855, 769, 692, 622, 559, 503, 452, 406, 365,
                  328, 295, 265, 238, 214, 192, 172, 154, 138, 124, 111, 99, 89, 80, 72, 64, 57, 51, 45,
                  40, 36, 32, 28, 25, 22, 19, 17, 15, 13, 11, 9, 8, 7, 6, 5, 4, 3, 2, 1]

try:
    os.mkdir(output_dir_merge + '/summary')
except IOError:
    pass

def windowing_concat(Seq_N):
    for windowsize in windowsize_set:
        N_S_set = [0, 0, 0]
        total_NS = ['']
        POS_old = 0
        for CHR in Seq_N:
            # calculate interval
            try:
                total_length = CHR.split('size')[1]
            except IndexError:
                total_length = CHR.split('length_')[1].split('_cov')[0]
            total_length = int(total_length)
            if POS_old + total_length >= windowsize:
                total_interval = int(total_length / windowsize) + 1
                total_NS += [''] * (total_interval)
            # windowing SNPs
            POS_set = Seq_N[CHR][0]
            NS_set = Seq_N[CHR][1]
            for i in range(0, len(POS_set)):
                POS = POS_set[i] + POS_old
                loci_POS = int(POS / windowsize)
                if total_NS[loci_POS] == '':
                    total_NS[loci_POS] = NS_set[i]
            POS_old += total_length
        N_S_set[0] += total_NS.count('N')
        N_S_set[1] += total_NS.count('S')
        try:
            N_S_set[2] = N_S_set[0] / N_S_set[1]
        except ZeroDivisionError:
            N_S_set[2] = 'N_only'
        Output.append('%s\t%s\t%s\t%s\t%s\t\n'%(donor_species,windowsize,
                                                N_S_set[0],N_S_set[1],N_S_set[2]))

def windowing(Seq_N):
    for windowsize in windowsize_set:
        N_S_set = [0, 0, 0]
        total_NS = ['']
        for CHR in Seq_N:
            # calculate interval
            try:
                total_length = CHR.split('size')[1]
            except IndexError:
                total_length = CHR.split('length_')[1].split('_cov')[0]
            total_length = int(total_length)
            if total_length >= windowsize:
                total_interval = int(total_length / windowsize) + 1
                total_NS += [''] * (total_interval)
            # windowing SNPs
            POS_set = Seq_N[CHR][0]
            NS_set = Seq_N[CHR][1]
            for i in range(0, len(POS_set)):
                POS = POS_set[i]
                loci_POS = int(POS / windowsize)
                if total_NS[loci_POS] == '':
                    total_NS[loci_POS] = NS_set[i]
        N_S_set[0] += total_NS.count('N')
        N_S_set[1] += total_NS.count('S')
        try:
            N_S_set[2] = N_S_set[0] / N_S_set[1]
        except ZeroDivisionError:
            N_S_set[2] = 'N_only'
        Output.append('%s\t%s\t%s\t%s\t%s\t\n'%(donor_species,windowsize,
                                                N_S_set[0],N_S_set[1],N_S_set[2]))

def readSNPfile(vcf_file):
    Seq_N = dict()
    for lines in open(vcf_file,'r'):
        lines_set = lines.replace('\n', '').replace('\r', '').split('\t')
        CHR = lines_set[0]
        POS = int(lines_set[1])
        N_S = lines_set[-2]
        if 'None' not in N_S and 'S' not in N_S:
            N_S = 'N'
        Seq_N.setdefault(CHR,[[],[]])
        Seq_N[CHR][0].append(POS)
        Seq_N[CHR][1].append(N_S)
    return Seq_N

# before remove rec
Output = []
all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
for vcf_file in all_vcf_file:
    if 'BaFr_clustercluster2' not in vcf_file:
        print(vcf_file)
        donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0]
        Seq_N = readSNPfile(vcf_file)
        windowing(Seq_N)

foutput = open(output_dir_merge + '/summary/all.donor.species.NSratio.txt', 'w')
foutput.write('donor_species\twindowsize\tN\tS\tNS_ratio\t\n')
foutput.write(''.join(Output))
foutput.close()

# after remove rec
vcf_name = '.raw.vcf.filtered.vcf.final.snp.txt'
Output = []
all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
for vcf_file in all_vcf_file:
    if 'BaFr_clustercluster2' not in vcf_file:
        print(vcf_file)
        donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0]
        Seq_N = readSNPfile(vcf_file)
        windowing(Seq_N)

foutput = open(output_dir_merge + '/summary/all.donor.species.NSratio.removerec.txt', 'w')
foutput.write('donor_species\twindowsize\tN\tS\tNS_ratio\t\n')
foutput.write(''.join(Output))
foutput.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# CHANGED after round 4 calculate NS ratio
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics

Round = 4
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')+\
glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/*')
genome_root = '/scratch/users/anniz44/genomes/donor_species/*/round*'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s'%(Round)
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/merge_genome/'%(Round)
output_dir_merge2 = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/merge_genome/'%(Round)
vcf_name = '.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf'
ref_filename = '.all.spades*.fasta'
fasta_name = '.fasta.corrected.fasta'
fastq_name = '.sorted.bam'

try:
    os.mkdir(output_dir_merge + '/summary')
except IOError:
    pass

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# Set up N or S
N_S_set = dict()
N_S_set['N']=0
N_S_set['S']=1
purines=['A','G']
pyrimidines=['C','T']
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
# Set up NS ratio cutoff
NSratioobserve_cutoff = 0
Min_SNP_highselect_cutoff = 1/2000
Max_SNP_highselect_cutoff = 0.02
total_SNP_position_cutoff = 2 # at least Xbp
singleevent = True # curate historical events
curate_empty = False # remove SNP type that is a subset of another SNP type
countOther = True #do not count other SNPs (non-gene SNPs) when calculating expected NS ratio
corecutoff = 0.9
if singleevent:
    vcf_name = '.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf.singleevent.vcf'
################################################### new class #########################################################
__metaclass__ = type

class SNP_gene:
    # create a class to store SNP_gene
    'a class to store SNP_gene'
    def init(self, gene):
        self.gene = gene
        self.position = dict()
        self.position.setdefault(gene,set())
        # [[N,S],freq]
        # not observed but predicted SNP pair by codon freq
        self.SNP_pair = {'A-T': [0,0],
                         'A-C': [0,0],
                         'G-C': [0,0],
                         'G-T': [0,0],
                         'A-G': [0,0],
                         'G-A': [0,0]}
        self.SNP_pair_freq = {'A-T': 0,
                         'A-C': 0,
                         'G-C': 0,
                         'G-T': 0,
                         'A-G': 0,
                         'G-A': 0}
        self.depth = []
        self.cov = 0
        self.mutposition = dict()
        self.mutposition.setdefault(gene, set())
        self.NSratio = [0,0,0]
        self.protein = ''
        self.minor_freq = []
    def addmutposition(self, gene, position):
        self.mutposition.setdefault(gene, set())
        self.mutposition[gene].add(position)
    def deletemutposition(self, gene,position):
        self.mutposition.setdefault(gene, set())
        if position == 0 :
            self.mutposition.pop(gene, None)
        else:
            self.mutposition[gene].discard(position)
    def addposition(self, gene, position,depth):
        self.position.setdefault(gene, set())
        self.position[gene].add(position)
        self.cov += 1
        self.depth.append(depth)
    def deleteposition(self, gene, position):
        if position == 0 :
            self.position.pop(gene, None)
        else:
            self.position[gene].discard(position)
    def addprotein(self, aa):
        self.protein += aa
    def mutpositioncal(self):
        self.mutpositionsum1 = {0: 0,
                         1: 0,
                         2: 0}
        self.mutpositionsum2 = {0: 0,
                                1: 0,
                                2: 0}
        self.mutpositionsum = 0
        for Chr in self.mutposition:
            self.mutposition2 = list(self.mutposition[Chr])
            self.mutposition2.sort()
            Total = len(self.mutposition2)
            self.mutpositionsum += Total
            if Total >= 1:
                self.mutpositionsum1[abs(self.mutposition2[0]) % 3] += 1
                self.mutpositionsum2[(self.mutposition2[-1]) % 3] += 1
                if Total != 1:
                    for i in range(0,len(self.mutposition2)-1):
                        self.mutpositionsum1[abs(self.mutposition2[i+1] - self.mutposition2[i]) % 3] += 1
                        self.mutpositionsum2[(self.mutposition2[i]) % 3] += 1
    def addSNP_pair(self, pair, position, count, unique_snp_count,depth = 0):
        self.SNP_pair_freq[pair] += unique_snp_count
        self.NSratio[position] += unique_snp_count
        if position < 2 and depth == 0:
            # add to NS for each SNP pair of reference genes
            self.SNP_pair[pair][position] += unique_snp_count
    def addpredictSNP_pair(self,refSNP_pair_sum):
        for pair in refSNP_pair_sum:
            self.SNP_pair[pair][0] += refSNP_pair_sum[pair][0]
            self.SNP_pair[pair][1] += refSNP_pair_sum[pair][1]
    def addalt(self,AllALT_frq):
        self.minor_freq.append(AllALT_frq)
    def deleteSNP_pair(self,SNP_gene):
        for position in [0,1]:
            self.NSratio[position] -= SNP_gene.NSratio[position]
            for pair in SNP_gene.SNP_pair:
                self.SNP_pair[pair][position] -= SNP_gene.SNP_pair[pair][position]
                if position == 0:
                    self.SNP_pair_freq[pair] -= SNP_gene.SNP_pair_freq[pair]
    def sum_SNP_pair(self):
        self.SNP_pair_sum = {'A-T': [0, 0],
                         'A-C': [0, 0],
                         'G-C': [0, 0],
                         'G-T': [0, 0],
                         'A-G': [0, 0],
                         'G-A': [0, 0]}
        for pair in self.SNP_pair:
            self.SNP_pair_sum[pair][0] += self.SNP_pair[pair][0]
            self.SNP_pair_sum[pair][1] += self.SNP_pair[pair][1]
    def dN_dS(self,SNP_gene_all,normalize=0):
        self.expectNSratio = 'No_expect'
        expectNSratio = [0, 0]
        if normalize == 1:
            for pair in self.SNP_pair_freq:
                # use selected codon NS ratio (SNP pair) * all genes freq
                expectNSratio[0] += self.SNP_pair[pair][0] * SNP_gene_all.SNP_pair_freq[pair]
                expectNSratio[1] += self.SNP_pair[pair][1] * SNP_gene_all.SNP_pair_freq[pair]
        else:
            for pair in self.SNP_pair_freq:
                expectNSratio[0] += self.SNP_pair[pair][0] * self.SNP_pair_freq[pair]
                expectNSratio[1] += self.SNP_pair[pair][1] * self.SNP_pair_freq[pair]
        if expectNSratio[1] > 0:
            # no expect S
            self.expectNSratio = expectNSratio[0] / expectNSratio[1]
        elif self.NSratio[0] == 0:
            # only S observed
            self.expectNSratio = 'expect_None'
        else:
            # only N observed
            self.expectNSratio = 'expect_N_only'
        if self.NSratio[1] > 0:
            # S observed
            self.NSratiosum = self.NSratio[0] / self.NSratio[1]
        elif self.NSratio[0] == 0:
            # only S observed
            self.NSratiosum = 'observe_None'
        else:
            # only N observed
            self.NSratiosum = 'observe_N_only'
        self.dNdS = self.NSratiosum
        if type(self.expectNSratio) == float \
                and type(self.NSratiosum) == float \
                and self.expectNSratio > 0:
            # calculate dNdS only if expected NSratio and observed NSratio are effective
            self.dNdS = self.NSratiosum / self.expectNSratio
            self.dNdS = '%.3f' % (self.dNdS)
            self.NSratiosum = '%.3f' % (self.NSratiosum)
            self.expectNSratio = '%.3f' % (self.expectNSratio)
        elif type(self.expectNSratio) == float:
            self.expectNSratio = '%.3f' % (self.expectNSratio)
        elif type(self.NSratiosum) == float:
            self.NSratiosum = '%.3f' % (self.NSratiosum)

################################################### Function ########################################################

def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def dnORds(amino1, amino2):
    if amino1 == amino2:
        return 'S'
    else:
        return 'N'

def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        if ALT_frq > 0:
            ALT_set.setdefault(ALT_frq, set())
            ALT_set[ALT_frq].add(alleles)
            ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def transitions(REF,ALT):
    if REF in pyrimidines:
        REF = complement[REF]
        ALT = complement[ALT]
    return '%s-%s'%(REF,ALT)

def curate_REF(allels_set,lines_set_sub):
    SNP_count_genome_count = [0,0,0,0]
    for Subdepth_all in lines_set_sub:
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '')
        Subdepth_set = Subdepth.split(',')
        Subdepth_set_int = []
        for sub_depth in Subdepth_set:
            Subdepth_set_int.append(int(sub_depth))
        if sum(Subdepth_set_int) > 0:
            major_alt_frq = max(Subdepth_set_int)
            major_alt_frq_index = Subdepth_set_int.index(major_alt_frq)
            major_alt = allels_set[major_alt_frq_index]
            SNP_count_genome_count[major_alt_frq_index] += 1
    major_alt_frq = max(SNP_count_genome_count)
    major_alt_frq_index = SNP_count_genome_count.index(major_alt_frq)
    major_alt = allels_set[major_alt_frq_index]
    return [major_alt, major_alt_frq_index]

def expectNSsub(record_name,record_seq,position=0):
    Total = int(len(record_seq)/3)
    temp_SNP_gene = SNP_gene()
    temp_SNP_gene.init(record_name)
    for i in range(0, (Total - position)):
        codon = record_seq[(i * 3 + position):((i + 1) * 3 + position)]
        try:
            codon_NSratio = codontable_NSratio[codon]
            temp_SNP_gene.addprotein(codontable[codon])
            for pair in codon_NSratio.SNP_pair:
                temp_SNP_gene.addSNP_pair(pair, 0, codon_NSratio.SNP_pair[pair][0],codon_NSratio.SNP_pair[pair][0],0)
                temp_SNP_gene.addSNP_pair(pair, 1, codon_NSratio.SNP_pair[pair][1],codon_NSratio.SNP_pair[pair][1],0)
        except KeyError:
            pass
    temp_SNP_gene.sum_SNP_pair()
    return [temp_SNP_gene.SNP_pair_sum,temp_SNP_gene.protein,position]

def expectNS(record_name,record_seq):
    Total = int(len(record_seq) / 3)
    temp_result = expectNSsub(record_name, record_seq)
    if len(temp_result[1]) < 0.8 * Total:
        temp_result = expectNSsub(record_name, record_seq,1)
        if len(temp_result[1]) < 0.8 * Total:
            temp_result = expectNSsub(record_name, record_seq,2)
            if len(temp_result[1]) < 0.8 * Total:
                return 'None'
    return temp_result

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Ref_NSratio=dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        Ref_seq.setdefault(record_id, record_seq)
        Ref_NSratio.setdefault(record_id,
                               expectNS(record_id, record_seq))
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    foutput = open(database + '.ref.NS.ratio', 'w')
    foutput_list = []
    #for Ref in Ref_NSratio:
    #    foutput_list.append('%s\t%s\t\n' % (Ref, Ref_NSratio[Ref]))
    #foutput.write(''.join(foutput_list))
    #foutput.close()
    return [Ref_seq,Ref_NSratio,Mapping,Mapping_loci,Reverse]

def contig_to_gene(CHR, POS):
    all_genes = Mapping_loci.get(CHR,[])
    Reverse_chr = 0
    for a_gene in all_genes:
        POS1, POS2, GENE = a_gene
        if POS >= POS1 and POS <= POS2:
            Ref_seq_chr = Ref_seq.get(GENE, 'None')
            Gene_length = len(Ref_seq_chr)
            if GENE in Reverse:  # reversed
                POS_gene = Gene_length-(int(POS-POS1))
                Reverse_chr = 1
            else:
                POS_gene = int(POS-POS1)+1
            codon_start = POS_gene - 1 - int((POS_gene - 1) % 3)
            return [GENE,POS_gene,codon_start,Ref_seq_chr,Reverse_chr]
    return []

def sumgene(SNP_gene_temp,genome_set_snp,donor_species,SNP_gene_species,Total,normalize = 1):
    High_select = False
    High_select2 = False
    N_temp = SNP_gene_temp.NSratio[0]
    S_temp = SNP_gene_temp.NSratio[1]
    Other_temp = SNP_gene_temp.NSratio[2]
    N_S_sum = N_temp + S_temp
    total_SNP = N_S_sum + Other_temp
    Gene_length = 1000
    Chr = SNP_gene_temp.gene
    total_SNP_position = len(SNP_gene_temp.mutposition[Chr])
    if Chr in Mapping:
        Gene_length = Mapping[Chr]
    if Chr == 'allspecies':
        Gene_length = 0
        for allchr in SNP_gene_all.position:
            Gene_length += Mapping.get(allchr, 0)
    new_line = '%s\t%s\t%s\t%s\t%s' % (donor_species, Chr, Total,
                                       Gene_length, genome_set_snp)
    if total_SNP > 0:
        new_line += '\t%s\t%s' % (total_SNP,total_SNP_position)
        SNP_gene_temp.dN_dS(SNP_gene_species, normalize)  # normalized
        new_line += ('\t%s\t%s\t%s\t%s' % (N_temp, S_temp, Other_temp,SNP_gene_temp.NSratiosum))
        new_line += '\t%s\t%s' % (SNP_gene_temp.expectNSratio, SNP_gene_temp.dNdS)
        for pair in SNP_gene_temp.SNP_pair:
            pair_freq = SNP_gene_temp.SNP_pair_freq[pair]
            pair_N = SNP_gene_temp.SNP_pair[pair][0]
            pair_S = SNP_gene_temp.SNP_pair[pair][1]
            new_line += ('\t%s\t%s:%s' % ('%d' % pair_freq, pair_N, pair_S))
    if N_S_sum > 0 and genome_set_snp > 1 \
            and total_SNP_position >= total_SNP_position_cutoff and \
            total_SNP_position / Gene_length >= Min_SNP_highselect_cutoff \
            and total_SNP_position / Gene_length <= Max_SNP_highselect_cutoff:
        High_select2 = True
        if SNP_gene_temp.NSratio[0] > SNP_gene_temp.NSratio[1] * NSratioobserve_cutoff:
            High_select = True
    new_line += '\t%s\n'%(High_select)
    return [new_line,High_select,High_select2]

def curate_SNP_type(SNP_genome_set):
    # remove SNP type that is a subset of another SNP type
    Major_alt_set = SNP_genome_set[-1]
    SNP_genome_set = list(set(SNP_genome_set[:-1]))
    SNP_total_length = len(Major_alt_set)
    SNP_genome_set_diff = set()
    for j in range(0, len(SNP_genome_set)):
        newSNP = SNP_genome_set[j]
        diff = ['%s %s'%(i, newSNP[i]) for i in range(0, SNP_total_length) if newSNP[i] != Major_alt_set[i]]
        SNP_genome_set_diff.add('\t'.join(diff))
    total_count = 0
    total_diff_length = len(SNP_genome_set_diff)
    SNP_genome_set_diff = list(SNP_genome_set_diff)
    for diff in SNP_genome_set_diff:
        diff_SNP_position = [i for i in range(0, total_diff_length) if all(elem in SNP_genome_set_diff[i] for elem in diff)]
        if len(diff_SNP_position) < 2:
            # no other SNP type contains this SNP type
            total_count += 1
    return total_count

def addSNP(SNP_gene_all_temp,SNP_gene_temp):
    SNP_gene_all_temp.addpredictSNP_pair(SNP_gene_temp.SNP_pair)
    SNP_gene_all_temp.NSratio[0] += SNP_gene_temp.NSratio[0]
    SNP_gene_all_temp.NSratio[1] += SNP_gene_temp.NSratio[1]
    for pair in SNP_gene_temp.SNP_pair_freq:
        # use selected genes frequency * all genes NS ratio (codon NS sum of all genes)
        SNP_gene_all_temp.SNP_pair[pair][0] += SNP_gene_temp.SNP_pair[pair][0]
        SNP_gene_all_temp.SNP_pair[pair][1] += SNP_gene_temp.SNP_pair[pair][1]
        SNP_gene_all_temp.SNP_pair_freq[pair] += SNP_gene_temp.SNP_pair_freq[pair]

def freq_call(vcf_file,Ref_seq, Ref_NSratio,SNP_gene_species,SNP_gene_all,SNP_gene_all_highselect,SNP_gene_all_highselect2,SNP_gene_all_flexible,SNP_gene_all_core,SNP_gene_species_highselect,Output2,donor_species):
    Output = []
    all_SNP_gene_temp = dict()
    Total = 0
    SNP_type = dict()
    for lines in open(vcf_file, 'r'):
        lines_set = lines.replace('\n','').replace('\r','').split('\t')
        # set up the basic
        if Total == 0:
            Total = len(lines_set) - 9
        Chr = lines_set[0]
        position = int(lines_set[1])
        gene_info = contig_to_gene(Chr, position)
        if gene_info != []:
            Chr, position, codon_start, Ref_seq_chr, Reverse_chr = gene_info# a gene
        else:
            Chr = Chr + '_other' # not a gene
        Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
        REF = lines_set[3]
        ALT_set = lines_set[4].split(',')
        allels_set = [REF] + ALT_set
        Total_alleles = len(allels_set)
        SNP_count_genome_count = [[0] * Total_alleles, '', '']
        SNP_type.setdefault(Chr, [''] * (Total+1)) # one more for major alt set
        New_gene = 0
        if Chr not in all_SNP_gene_temp:
            SNP_gene_temp = SNP_gene()
            SNP_gene_temp.init(Chr)
            all_SNP_gene_temp.setdefault(Chr, SNP_gene_temp)
            New_gene = 1
        # set up SNP_gene
        SNP_gene_temp = all_SNP_gene_temp[Chr]
        SNP_gene_temp.addposition(Chr, position, Depth)
        SNP_gene_all.addposition(Chr, position, Depth)
        SNP_gene_species.addposition(Chr, position, Depth)
        # count genome set with SNP
        genome_ID = 0
        for Subdepth_all in lines_set[9:]:
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '')
            Subdepth_set = Subdepth.split(',')
            Subdepth_set_int = []
            for sub_depth in Subdepth_set:
                Subdepth_set_int.append(int(sub_depth))
            major_alt_frq = max(Subdepth_set_int)
            major_alt_frq_index = Subdepth_set_int.index(major_alt_frq)
            major_alt = allels_set[major_alt_frq_index]
            SNP_count_genome_count[0][major_alt_frq_index] += 1
            SNP_count_genome_count[1] += major_alt
            SNP_type[Chr][genome_ID] += major_alt
            genome_ID += 1
            SNP_count_genome_count[1] += '\t'
        # currate REF and ALT
        REF, REF_where = curate_REF(allels_set, lines_set[9:])
        SNP_type[Chr][genome_ID] += REF # major alt set
        ALT_set = allels_set
        ALT_set.remove(REF)
        # calculate N or S
        refSNP_pair_sum_all = Ref_NSratio.get(Chr, 'None')
        if refSNP_pair_sum_all != 'None':
            # predicted NS store in SNP_pair
            if New_gene == 1:
                refSNP_pair = refSNP_pair_sum_all[0]
                SNP_gene_temp.addpredictSNP_pair(refSNP_pair)
                SNP_gene_species.addpredictSNP_pair(refSNP_pair)
                SNP_gene_all.addpredictSNP_pair(refSNP_pair)
            #  observed NS ratio calculated
            refSNP_condon_start = refSNP_pair_sum_all[-1]
            codon_start = position - 1 - int((position - 1) % 3) + refSNP_condon_start
            if codon_start <= position - 1:
                Ref_seq_chr = Ref_seq[Chr]
                SNP_seq_chr = Ref_seq_chr
                Ref_seq_chr = causeSNP(Ref_seq_chr, position, REF,Reverse_chr)
                Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                if len(Ref_seq_codon) == 3:
                    Ref_seq_aa = translate(Ref_seq_codon)[0]
                    ALT_num = 1
                    for ALT in ALT_set:
                        ALT_frq = SNP_count_genome_count[0][ALT_num]
                        SNP_seq_chr = causeSNP(SNP_seq_chr, position, ALT,Reverse_chr)
                        SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                        SNP_seq_aa = translate(SNP_seq_codon)[0]
                        temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                        SNP_count_genome_count[2] += temp_NorS
                        SNP_pair = transitions(REF, ALT)
                        ALT_num += 1
                        SNP_gene_temp.addSNP_pair(SNP_pair, N_S_set[temp_NorS],
                                                  ALT_frq, 1, Depth)
                        SNP_gene_temp.addmutposition(Chr, position)
                        SNP_gene_all.addSNP_pair(SNP_pair, N_S_set[temp_NorS],
                                                 ALT_frq, 1, Depth)
                        SNP_gene_all.addmutposition('allspecies', position)
                        SNP_gene_species.addSNP_pair(SNP_pair, N_S_set[temp_NorS],
                                                     ALT_frq, 1, Depth)
                        SNP_gene_species.addmutposition(donor_species, position)
        else:
            # not a gene
            ALT_num = 1
            for ALT in ALT_set:
                ALT_frq = SNP_count_genome_count[0][ALT_num]
                SNP_pair = transitions(REF, ALT)
                # add to P
                SNP_gene_temp.addSNP_pair(SNP_pair, 2,
                                          ALT_frq, 1, Depth)
                SNP_gene_temp.addmutposition(Chr, position)
                if countOther:
                    SNP_gene_all.addSNP_pair(SNP_pair, 2,
                                             ALT_frq, 1, Depth)
                    SNP_gene_all.addmutposition('allspecies', position)
                SNP_gene_species.addSNP_pair(SNP_pair, 2,
                                             ALT_frq, 1, Depth)
                SNP_gene_species.addmutposition(donor_species, position)
                ALT_num += 1
        # output genome SNP
        Output.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n' % (Chr, position, REF, ALT,
                                                              SNP_count_genome_count[0][0],
                                                              SNP_count_genome_count[0][1],
                                                              SNP_count_genome_count[2], SNP_count_genome_count[1]))
    if all_SNP_gene_temp!= dict():
        for Chr in all_SNP_gene_temp:
            SNP_gene_temp = all_SNP_gene_temp[Chr]
            estimate_total_genome_set = len(set(SNP_type[Chr][:-1]))-1
            if curate_empty and estimate_total_genome_set > 1:
                # remove SNP type that is a subset of another SNP type
                total_genome_set = curate_SNP_type(SNP_type[Chr]) # not including the major alt alt
                print(total_genome_set,estimate_total_genome_set)
                print(set(SNP_type[Chr]))
                estimate_total_genome_set = total_genome_set
            sumgene_line,High_select,High_select2 = sumgene(SNP_gene_temp,estimate_total_genome_set,donor_species,SNP_gene_species,Total,1)
            Output2.append(sumgene_line)
            if High_select:
                # use selected genes frequency * all genes NS ratio (codon NS sum of all genes)
                addSNP(SNP_gene_all_highselect, SNP_gene_temp)
            if High_select2:
                addSNP(SNP_gene_all_highselect2, SNP_gene_temp)
                addSNP(SNP_gene_species_highselect, SNP_gene_temp)
            if '%s:%s'%(donor_species,Chr) in Core:
                core = Core['%s:%s'%(donor_species,Chr)]
                if core == 'all_core':
                    addSNP(SNP_gene_all_core, SNP_gene_temp)
                elif core == 'species_flexible':
                    addSNP(SNP_gene_all_flexible, SNP_gene_temp)
            else:
                print('missing genes %s in %s'%(Chr, donor_species))
        foutput = open(vcf_file + '.frq.snp', 'w')
        foutput.write('#CHR\tPOS\tMajor_ALT\tMinor_ALT\tMajor_ALT_frq\tMinor_Alt_frq\tN_or_S\tgenotype_allgenomes\n')
        foutput.write(''.join(Output))
        foutput.close()
    return Total

# calculate codon freq
codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
codontable_NSratio = dict()
for codon in codontable:
    SNP_gene_temp = SNP_gene()
    SNP_gene_temp.init(codon)
    codontable_NSratio.setdefault(codon, SNP_gene_temp)
    for position in range(0, 3):
        REF = codon[position]
        for ALT in Allels_order:
            if ALT != REF:
                new_codon = causeSNP(codon, position+1, ALT,0)
                temp_NorS = dnORds(translate(codon)[0], translate(new_codon)[0])
                SNP_pair = transitions(REF, ALT)
                SNP_gene_temp.addSNP_pair(SNP_pair, N_S_set[temp_NorS], 1, 1, 0)

# load core/flexible genes
Core = dict()
for lines in open(os.path.join(input_script,'all.denovo.gene.faa.allpangenome.sum.allsum.species.multispecies.txt')):
    lines_set = lines.replace('\n', '').replace('\r', '').split('\t')
    donor_species = lines_set[-3]
    geneID = lines_set[-1]
    core = lines_set[-4]
    Core.setdefault('%s:%s'%(donor_species,geneID),core)

# set up sum all species
SNP_gene_all = SNP_gene() # all denovo mutation
SNP_gene_all.init('allspecies')
SNP_gene_all_highselect = SNP_gene() # all mutations of highly selected genes
SNP_gene_all_highselect.init('allspecies_highselect')
SNP_gene_all_highselect2 = SNP_gene() # all mutations of highly selected genes no NS ratio cutoff
SNP_gene_all_highselect2.init('allspecies_highselect_noNS')
SNP_gene_all_flexible = SNP_gene() # all mutations of flexible genes
SNP_gene_all_flexible.init('allspecies_flexible')
SNP_gene_all_core = SNP_gene() # all mutations of flexible genes
SNP_gene_all_core.init('allspecies_core')
# process each vcf file
Output2 = []
all_vcf_file = glob.glob(os.path.join(output_dir_merge + '/vcf', '*%s' % (vcf_name)))+\
glob.glob(os.path.join(output_dir_merge2 + '/vcf', 'D*%s' % (vcf_name)))
#all_vcf_file = glob.glob(os.path.join(output_dir_merge2, 'H*%s' % (vcf_name)))+\
#glob.glob(os.path.join(output_dir_merge2, 'P*%s' % (vcf_name)))
for vcf_file in all_vcf_file:
    print(vcf_file)
    vcf_ref_file_name = os.path.split(vcf_file)[-1].split('.all.')[0]
    donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0].split('.flt.snp.vcf')[0]
    database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
    if len(database) > 1:
        print(vcf_file, database)
    database = database[0]
    ref_dir, ref_name = os.path.split(database)
    database_file = database.replace('.fasta', '.fna')
    print('running %s cutoff %s' % (donor_species,total_SNP_position_cutoff))
    Ref_seq, Ref_NSratio, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
    SNP_gene_species = SNP_gene()  # all mutations of a species
    SNP_gene_species.init(donor_species)
    SNP_gene_species_highselect = SNP_gene()  # all mutations of highly selected genes in a species
    SNP_gene_species_highselect.init(donor_species + '_highselect')
    Total = freq_call(vcf_file, Ref_seq, Ref_NSratio, SNP_gene_species, SNP_gene_all, SNP_gene_all_highselect,
                      SNP_gene_all_highselect2,
                      SNP_gene_all_flexible,SNP_gene_all_core,SNP_gene_species_highselect,
                          Output2, donor_species)
    if sum(SNP_gene_species.NSratio) > 0:
        # there's a SNP
        sumgene_line, High_select, High_select2 = sumgene(SNP_gene_species, 1, donor_species, SNP_gene_species, Total, 0)
        Output2.append(sumgene_line)
        print('allspecies',SNP_gene_species.dNdS, SNP_gene_species.NSratio,
            SNP_gene_species.NSratiosum)
        if sum(SNP_gene_species_highselect.NSratio) > 0:
            sumgene_line, High_select, High_select2 = sumgene(SNP_gene_species_highselect, 1, donor_species, SNP_gene_species, Total, 1)
            Output2.append(sumgene_line)
            print('HSspecies',SNP_gene_species_highselect.dNdS, SNP_gene_species_highselect.NSratio,
                SNP_gene_species_highselect.NSratiosum)
    else:
        Output2.append('%s\t0\t\n'%(donor_species))

# sum all species dNdS
sumgene_line,High_select,High_select2 = sumgene(SNP_gene_all,1,'allspecies',SNP_gene_all,'None',0)
Output2.append(sumgene_line)
# HS genes
sumgene_line,High_select,High_select2 = sumgene(SNP_gene_all_highselect,1,'allspecies',SNP_gene_all,'None',1)
Output2.append(sumgene_line)
# HS genes no NS ratio cutoff
sumgene_line,High_select,High_select2 = sumgene(SNP_gene_all_highselect2,1,'allspecies',SNP_gene_all,'None',1)
Output2.append(sumgene_line)
# flexible genes
sumgene_line,High_select,High_select2 = sumgene(SNP_gene_all_flexible,1,'allspecies',SNP_gene_all,'None',1)
Output2.append(sumgene_line)
# core genes
sumgene_line,High_select,High_select2 = sumgene(SNP_gene_all_core,1,'allspecies',SNP_gene_all,'None',1)
Output2.append(sumgene_line)

# output
foutput = open(output_dir_merge + '/summary/all.WGSgenome.donor.species.removerec.sum.noother.species.txt', 'w')
foutput.write('#donor_species\tgene\tNo.genome\tgene_length\ttotal_SNP_genomeset\tNo.SNP\tNo.SNP_position\tN\tS\tOther\tobserved_ratio\texpected_ratio\tdNdS\t' +\
            'A-T_freq\tA-T_N:S\tA-C_freq\tA-C_N:S\tG-C_freq\tG-C_N:S\tG-T_freq\tG-T_N:S\tA-G_freq\tA-G_N:S\tG-A_freq\tG-A_N:S\tHigh_selected\n')
foutput.write(''.join(Output2))
foutput.close()

print('all',SNP_gene_all.dNdS)
print('HS',SNP_gene_all_highselect.dNdS)
print('HS no NS cutoff',SNP_gene_all_highselect2.dNdS)
print('flexible',SNP_gene_all_flexible.dNdS)
print('core',SNP_gene_all_core.dNdS)

################################################### END ########################################################
################################################### SET PATH ########################################################
# after calculate NS ratio, extract genes and run annotation
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics

# set up path -> selected
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/annotate_all'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round*')+\
glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/round*')
genome_root = '/scratch/users/anniz44/genomes/donor_species/*/round*'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome'
ref_filename = '.all.spades*.fasta'
input_summary = output_dir_merge + '/summary/all.WGSgenome.donor.species.removerec.sum.noother.species.txt'
output_gene = output_dir_merge + '/summary/all.selected.gene.faa'
output_gene_dna = output_dir_merge + '/summary/all.selected.gene.fna'

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# function
def sum_donor_species(input_summary,selected = 1):
    Donor_species = dict()
    for lines in open(input_summary,'r'):
        lines_set = lines.replace('\n','').split('\t')
        if not lines.startswith('#'):
            if selected == 0 or lines_set[-1] == 'True':
                # highly selected genes or all genes
                gene_name = lines_set[1]
                donor_species = lines_set[0]
                Donor_species.setdefault(donor_species,[])
                Donor_species[donor_species].append(gene_name)
    return Donor_species

def sum_vcf(vcf,Donor_species):
    donor_species = os.path.split(vcf)[-1].split('.all.flt.snp.vcf')[0]
    Donor_species.setdefault(donor_species, [])
    for lines in open(vcf,'r'):
        lines_set = lines.replace('\n','').split('\t')
        if not lines.startswith('#'):
            if '*' in lines_set[-1]:
                # truncated genes
                gene_name = lines_set[-4]
                Donor_species[donor_species].append(gene_name)
    return Donor_species

def extract_donor_species(donor_species,input_fasta,gene_name_list,output_fasta):
    gene_name_extract = []
    donor_species_set = donor_species.split('_')
    record_id_set = set()
    try:
        donor_species = '%s_%s_%s' % (donor_species_set[0],
                                      donor_species_set[1][0:min(6, len(donor_species_set[1]))],
                                      donor_species_set[2][0:min(6, len(donor_species_set[2]))])
    except IndexError:
        donor_species = '%s_%s_%s' % (donor_species_set[0],
                                      donor_species_set[1],
                                      donor_species_set[1])
    if 'cluster' in donor_species_set[-1]:
        try:
            donor_species += '_CL' + donor_species_set[-2].split('cluster')[1]
        except IndexError:
            donor_species += '_CL' + donor_species_set[-1].split('cluster')[1]
    gene_name_list = set(gene_name_list)
    for record in SeqIO.parse(input_fasta, 'fasta'):
        record_id = str(record.id)
        temp_line = '%s\t%s\t%s\t'%('_'.join(donor_species_set),donor_species,record_id)
        if record_id in gene_name_list and record_id not in record_id_set:
            gene_name_extract.append(record_id)
            record_id = 'C_%s_G_%s'%(record_id.split('_')[1], record_id.split('_')[-1])
            output_fasta.append('>%s__%s\n%s\n'%(donor_species, record_id, str(record.seq)))
            record_id_set.add(str(record.id))
            temp_line += '%s__%s\t\n'%(donor_species,record_id)
            change_name.add(temp_line)
    if len(gene_name_extract) < len(gene_name_list):
        print('missing genes in files %s extracted %s genes of %s genes'%(input_fasta,len(gene_name_extract),len(gene_name_list)))
        print('missing genes ' + ' '.join([gene_name for gene_name in gene_name_list if gene_name not in gene_name_extract]))
    return output_fasta

def annotation(all_filter_gene_fasta_file,script_name):
    # run prokka
    cmdsprokka = 'py37\nprokka --kingdom Bacteria --outdir %s/prokka_%s  --protein %s --locustag Bacter %s/%s\n' % \
                 (output_dir_merge, os.path.split(all_filter_gene_fasta_file)[-1],
                  all_filter_gene_fasta_file,
                  output_dir_merge,
                  os.path.split(all_filter_gene_fasta_file)[-1].replace('.faa', '.fna'))
    f1 = open(os.path.join(input_script_sub, 'prokka.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmdsprokka))
    f1.close()
    # run cluster
    cutoff = 0.7
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                   % ('usearch', all_filter_gene_fasta_file, cutoff, all_filter_gene_fasta_file,
                      all_filter_gene_fasta_file, 40))
    all_filter_gene_fasta_file = all_filter_gene_fasta_file + '.cluster.aa'
    # run metacyc
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/mit_alm/database/metacyc/protseq.fsa'
    cmds = ("diamond blastp --query %s --db %s.dmnd --out %s.metacyc.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    f1 = open(os.path.join(input_script_sub, 'metacyc.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    # run eggnog
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/eggnog/xaa.hmm'
    cmds = ('hmmsearch --tblout %s.eggnog.1.txt --cpu 40 -E %s %s %s\n') %(all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.1.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xab.hmm'
    cmds = ('hmmsearch --tblout %s.eggnog.2.txt --cpu 40 -E %s %s %s\n') % (
    all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.2.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xac.hmm'
    cmds = ('hmmsearch --tblout %s.eggnog.3.txt --cpu 40 -E %s %s %s\n') % (
    all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.3.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmds))
    f1.close()
    # run kegg
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/kegg/kofam/profiles/prokaryote/prokaryote.hmm'
    cmds = ('hmmsearch --tblout %s.kegg.txt --cpu 40 -E %s %s %s\n') %(all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'kegg.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    # run customed database
    cutoff = 80
    cutoff2 = 80
    cmds = ''
    database = '/scratch/users/anniz44/scripts/database/SARG.db.fasta'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.SARG.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 50
    cutoff2 = 50
    database = '/scratch/users/anniz44/scripts/database/AHR.aa.db'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.AHR.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 60
    cutoff2 = 80
    database = '/scratch/users/anniz44/scripts/database/Butyrate.pro.aa'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.buty.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/anniz44/scripts/database/IntI1_database.fasta'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.int.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/anniz44/scripts/database/SRB.AA'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.SRB.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 0.01
    database = '/scratch/users/anniz44/scripts/database/NR.hmm'
    cmds += ('hmmsearch --tblout %s.NR.txt --cpu 40 -E %s %s %s\n') %(all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'customed.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    # all scripts
    f1 = open(os.path.join(input_script, script_name), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    f1.write(cmd_cluster)
    for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.sh')):
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()

def annotation_dna(all_filter_gene_fasta_file):
    # run cluster
    cutoff = 0.7
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                   % ('usearch', all_filter_gene_fasta_file, cutoff, all_filter_gene_fasta_file,
                      all_filter_gene_fasta_file, 40))
    all_filter_gene_fasta_file = all_filter_gene_fasta_file + '.cluster.aa'
    # run metacyc
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/mit_alm/database/metacyc/protseq.fsa'
    cmds = ("diamond blastx --query %s --db %s.dmnd --out %s.metacyc.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    f1 = open(os.path.join(input_script_sub, 'metacyc.sh'), 'a')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    # run eggnog
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/eggnog/xaa.hmm'
    cmds = ('#hmmsearch --tblout %s.eggnog.1.txt --cpu 40 -E %s %s %s\n') %(all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.1.sh'), 'a')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xab.hmm'
    cmds = ('#hmmsearch --tblout %s.eggnog.2.txt --cpu 40 -E %s %s %s\n') % (
    all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.2.sh'), 'a')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xac.hmm'
    cmds = ('#hmmsearch --tblout %s.eggnog.3.txt --cpu 40 -E %s %s %s\n') % (
    all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.3.sh'), 'a')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmds))
    f1.close()
    # run kegg
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/kegg/kofam/profiles/prokaryote/prokaryote.hmm'
    cmds = ('#hmmsearch --tblout %s.kegg.txt --cpu 40 -E %s %s %s\n') %(all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'kegg.sh'), 'a')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    # run customed database
    cutoff = 80
    cutoff2 = 80
    cmds = ''
    database = '/scratch/users/anniz44/scripts/database/SARG.db.fasta'
    cmds += ("diamond blastx --query %s --db %s.dmnd --out %s.SARG.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 50
    cutoff2 = 50
    database = '/scratch/users/anniz44/scripts/database/AHR.aa.db'
    cmds += ("diamond blastx --query %s --db %s.dmnd --out %s.AHR.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 60
    cutoff2 = 80
    database = '/scratch/users/anniz44/scripts/database/Butyrate.pro.aa'
    cmds += ("diamond blastx --query %s --db %s.dmnd --out %s.buty.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/anniz44/scripts/database/IntI1_database.fasta'
    cmds += ("diamond blastx --query %s --db %s.dmnd --out %s.int.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/anniz44/scripts/database/SRB.AA'
    cmds += ("diamond blastx --query %s --db %s.dmnd --out %s.SRB.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 0.01
    database = '/scratch/users/anniz44/scripts/database/NR.dna.new.hmm'
    cmds += ('hmmsearch --tblout %s.NR.txt --cpu 40 -E %s %s %s\n') %(all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'customed.sh'), 'a')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    # all scripts
    f1 = open(os.path.join(input_script, 'allannotate.dna.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    f1.write(cmd_cluster)
    for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.sh')):
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()

# extract highly selected sequences
output_fasta = []
output_fasta_dna = []
change_name = set()
Donor_species = sum_donor_species(input_summary)
for donor_species in Donor_species:
    print(donor_species)
    database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
    database = database[0]
    ref_dir, ref_name = os.path.split(database)
    input_fasta = database.replace('.fasta', '.faa')
    input_fasta_dna = database.replace('.fasta', '.faa')
    gene_name_list = Donor_species[donor_species]
    try:
        f1 = open(input_fasta, 'r')
    except FileNotFoundError:
        os.system('prodigal -q -i %s -a %s' % (database, input_fasta))
    try:
        f1 = open(input_fasta_dna, 'r')
    except FileNotFoundError:
        os.system('prodigal -q -i %s -d %s' % (database, input_fasta_dna))
    if gene_name_list != []:
        if input_fasta != [] and input_fasta_dna != []:
            output_fasta = extract_donor_species(donor_species,input_fasta, gene_name_list, output_fasta)
            output_fasta_dna = extract_donor_species(donor_species,input_fasta_dna, gene_name_list, output_fasta_dna)
        else:
            print('missing files for %s' % (donor_species))
    else:
        print('missing genes for %s' % (donor_species))

f1 = open(output_gene, 'w')
f1.write(''.join(output_fasta))
f1.close()
f1 = open(output_gene + '.changename.txt', 'w')
f1.write(''.join(list(change_name)))
f1.close()
f1 = open(output_gene_dna, 'w')
f1.write(''.join(output_fasta_dna))
f1.close()

# run clustering
#annotation(output_gene)
#annotation_dna(output_gene_dna)

# extract all sequences
output_gene = output_dir_merge + '/summary/all.denovo.gene.faa'
output_gene_dna = output_dir_merge + '/summary/all.denovo.gene.fna'
output_fasta = []
output_fasta_dna = []
Donor_species = sum_donor_species(input_summary,0)
change_name = set()
for donor_species in Donor_species:
    print(donor_species)
    database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
    database = database[0]
    ref_dir, ref_name = os.path.split(database)
    input_fasta = database.replace('.fasta', '.faa')
    input_fasta_dna = database.replace('.fasta', '.faa')
    gene_name_list = Donor_species[donor_species]
    try:
        f1 = open(input_fasta, 'r')
    except FileNotFoundError:
        os.system('prodigal -q -i %s -a %s' % (database, input_fasta))
    try:
        f1 = open(input_fasta_dna, 'r')
    except FileNotFoundError:
        os.system('prodigal -q -i %s -d %s' % (database, input_fasta_dna))
    if gene_name_list != []:
        if input_fasta != [] and input_fasta_dna != []:
            output_fasta = extract_donor_species(donor_species,input_fasta, gene_name_list, output_fasta)
            output_fasta_dna = extract_donor_species(donor_species,input_fasta_dna, gene_name_list, output_fasta_dna)
        else:
            print('missing files for %s' % (donor_species))
    else:
        print('missing genes for %s' % (donor_species))

f1 = open(output_gene, 'w')
f1.write(''.join(output_fasta))
f1.close()
f1 = open(output_gene + '.changename.txt', 'w')
f1.write(''.join(list(change_name)))
f1.close()
f1 = open(output_gene_dna, 'w')
f1.write(''.join(output_fasta_dna))
f1.close()

annotation(output_gene,'allannotate.sh')

# run cluster
cutoff = 0.7
cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                   % ('usearch', output_gene, cutoff, output_gene,
                      output_gene, 40))
os.system(cmd_cluster)

# extract truncated sequences
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome/'
output_dir_merge2 = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round4/merge_genome/'
vcf_name = '.final.vcf.removerec.snp.txt'
output_gene = output_dir_merge + '/summary/all.trunc.gene.faa'
output_gene_dna = output_dir_merge + '/summary/all.trunc.gene.fna'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/annotate_trunc'
try:
    os.mkdir(input_script_sub)
except IOError:
    pass

output_fasta = []
output_fasta_dna = []
all_vcf_file = glob.glob(os.path.join(output_dir_merge + '/vcf/', '*%s' % (vcf_name)))+\
glob.glob(os.path.join(output_dir_merge2 + '/vcf/', 'D*%s' % (vcf_name)))
Donor_species = dict()
for vcf in all_vcf_file:
    sum_vcf(vcf,Donor_species)

change_name = set()
for donor_species in Donor_species:
    print(donor_species)
    gene_name_list = Donor_species[donor_species]
    if gene_name_list!= []:
        database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
        database = database[0]
        ref_dir, ref_name = os.path.split(database)
        input_fasta = database.replace('.fasta', '.faa')
        input_fasta_dna = database.replace('.fasta', '.faa')
        try:
            f1 = open(input_fasta, 'r')
        except FileNotFoundError:
            os.system('prodigal -q -i %s -a %s' % (database, input_fasta))
        try:
            f1 = open(input_fasta_dna, 'r')
        except FileNotFoundError:
            os.system('prodigal -q -i %s -d %s' % (database, input_fasta_dna))
        if gene_name_list != []:
            if input_fasta != [] and input_fasta_dna != []:
                output_fasta = extract_donor_species(donor_species,input_fasta, gene_name_list, output_fasta)
                output_fasta_dna = extract_donor_species(donor_species,input_fasta_dna, gene_name_list, output_fasta_dna)
            else:
                print('missing files for %s' % (donor_species))
        else:
            print('missing genes for %s' % (donor_species))

f1 = open(output_gene, 'w')
f1.write(''.join(output_fasta))
f1.close()
f1 = open(output_gene + '.changename.txt', 'w')
f1.write(''.join(list(change_name)))
f1.close()
f1 = open(output_gene_dna, 'w')
f1.write(''.join(output_fasta_dna))
f1.close()

annotation(output_gene,'allannotate.trunc.sh')

################################################### END ########################################################
################################################### SET PATH ########################################################
# CHANGED find highly selected genes in one species across populations and donors
import os
import glob
import copy
from Bio import SeqIO
from Bio.Seq import Seq

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/annotate'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome/'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
all_fasta = os.path.join(output_dir + '/summary', 'all.denovo.gene.faa')
all_fasta_HS = os.path.join(output_dir + '/summary', 'all.selected.gene.faa')
input_summary = output_dir_merge + '/summary/all.WGSgenome.donor.species.removerec.sum.noother.species.txt'
input_summary_nosingleevent = output_dir_merge + '/summary/all.WGSgenome.donor.species.removerec.sum.noother.species.nosingleevent.txt'
pre_cluster = ''
NSratioobserve_cutoff = 0
SNP_cutoff = ['Bifidobacterium_longum','Bifidobacterium_adolescentis','Bifidobacterium_pseudocatenulatum'] # SNP cutoff set as 3 for HS genes

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# functions
# clustering
def cluster_uc(cluster_input):
    Clusters = dict()
    High_select2 = set()
    High_select2_3SNP = dict()
    High_select2_output = []
    for lines in open(cluster_input, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        cluster = line_set[1]
        record_name = line_set[8]
        record_name2 = line_set[9]
        Clusters.setdefault(record_name, cluster)
        if record_name2!= '*':
            donor = record_name.split('_')[0]
            donor_species = record_name.split('__')[0]
            species = donor_species.replace(donor + '_', '')
            donor2 = record_name2.split('_')[0]
            donor_species2 = record_name2.split('__')[0]
            species2 = donor_species2.replace(donor2 + '_', '')
            if species == species2:
                if species not in SNP_cutoff:
                    High_select2.add(record_name2)
                    High_select2.add(record_name)
                else:
                    High_select2_3SNP.setdefault(record_name,set())
                    High_select2_3SNP[record_name].add(record_name2)
                    High_select2_3SNP.setdefault(record_name2, set())
                    High_select2_3SNP[record_name2].add(record_name)
    if High_select2_3SNP!= dict():
        for record in High_select2_3SNP:
            if len(High_select2_3SNP[record]) >= 3:
                High_select2.append(record)
    return [Clusters,High_select2_output,High_select2]

# correct highly selected genes dNdS
def format_freq(freq):
    freq = freq.split(':')
    return [int(freq[0]),int(freq[1])]

def init_highselect_notHS(line_set):
    allspecies_highselect = line_set
    for i in [3,5,6,7,8,-13,-11,-9,-7,-5,-3]:
        allspecies_highselect[i] = int(allspecies_highselect[i])
    for i in [-12,-10,-8,-6,-4,-2]:
        allspecies_highselect[i] = format_freq(line_set[i])
    return allspecies_highselect

def init_highselect_empty(line_set):
    allspecies_highselect = line_set
    for i in [3]:
        allspecies_highselect[i] = int(allspecies_highselect[i])
    for i in [5,6,7,8,9,10,11,12,-13,-11,-9,-7,-5,-3]:
        allspecies_highselect[i] = 0
    for i in [-12,-10,-8,-6,-4,-2]:
        allspecies_highselect[i] = format_freq('0:0')
    return allspecies_highselect

def add_freq(freq_new,freq_old):
    temp_result = format_freq(freq_new)
    freq_old[0] += temp_result[0]
    freq_old[1] += temp_result[1]

def add_highselect(line_set,allspecies_highselect):
    for i in [3,5,6,7,8,-13,-11,-9,-7,-5,-3]:
        allspecies_highselect[i] += int(line_set[i])
    for i in [-12,-10,-8,-6,-4,-2]:
        add_freq(line_set[i], allspecies_highselect[i])
    return allspecies_highselect

def calculate_NS(allspecies,allspecies_highselect):
    # expect NS
    # expected NS ratio
    tempNS = [0, 0]
    i = -12
    for freq in allspecies:
        tempNS[0] += freq * allspecies_highselect[i][0]
        tempNS[1] += freq * allspecies_highselect[i][1]
        i += 2
    if tempNS[1] > 0:
        allspecies_highselect[11] = tempNS[0] / tempNS[1] # expect
    else:
        allspecies_highselect[11] = 'expect_N_only'
    # NS ratio
    if allspecies_highselect[8] > 0:
        allspecies_highselect[10] = allspecies_highselect[7]/allspecies_highselect[8] # observe
        if allspecies_highselect[11] != 'expect_N_only':
            # dNdS
            allspecies_highselect[12] = allspecies_highselect[10] / allspecies_highselect[11]# dnds
    elif allspecies_highselect[7] > 0:
        allspecies_highselect[10] = 'observe_N_only'
        allspecies_highselect[12] = 'observe_N_only'
    return allspecies_highselect

def output_highlight(allspecies_highselect,allspecies,newoutput):
    allspecies_highselect = calculate_NS(allspecies, allspecies_highselect)
    for i in range(0, len(allspecies_highselect)):
        if i in [14, 16, 18, 20, 22, 24]:
            allspecies_highselect[i] = '%s:%s' % (allspecies_highselect[i][0],
                                                  allspecies_highselect[i][1])
        else:
            allspecies_highselect[i] = str(allspecies_highselect[i])
    newoutput.append('\t'.join(allspecies_highselect) + '\n')
    return newoutput

def genus_species(donor_species):
    genus = donor_species.split('_')[1].replace('BA', 'Bifidobacterium').replace('BL',
                                                                                 'Bifidobacterium').replace(
        'PB', 'Parabacteroides')
    if any(item in donor_species for item in ['BA','BL','PB']):
        species = donor_species.split('_')[1].replace('BA', 'Bifidobacterium_adolescentis').replace('BL',
                                                                                                      'Bifidobacterium_longum').replace(
            'PB', 'Parabacteroides_butyrate')
    else:
        species = (donor_species.split('_')[1] + '_' +
                   donor_species.split('_')[2])
    return [genus,species]

def add_new_selection(input_summary,High_select2):
    Donor_species = dict()
    Donor_species_notHS = dict()
    Genus = dict()
    Genus_notHS = dict()
    Species = dict()
    Species_notHS = dict()
    newoutput = []
    # use selected codon NS ratio (SNP pair) * all genes freq of no single event
    for lines in open(input_summary_nosingleevent, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        Selected = line_set[-1]
        if Selected == 'False':
            donor_species = line_set[0]
            record_id = line_set[1]
            # all species
            if record_id == 'allspecies':
                # all genes freq
                allspecies = [int(line_set[-13]), int(line_set[-11]), int(line_set[-9]),
                              int(line_set[-7]), int(line_set[-5]), int(line_set[-3])]
                allspecies_highselect = init_highselect_empty(line_set)
                allspecies_highselect[1] = 'allspecies_highselect'
                newoutput.append(lines)
            elif record_id in ['allspecies_flexible', 'allspecies_core']:
                newoutput.append(lines)
            elif donor_species == record_id:
                # donor_species not high select set up
                line_set = lines.split('\n')[0].split('\t')
                donor_species_highselect = init_highselect_notHS(line_set)
                Donor_species_notHS.setdefault(donor_species,
                                               donor_species_highselect)
                line_set = lines.split('\n')[0].split('\t')
                donor_species_highselect = init_highselect_empty(line_set)
                donor_species_highselect[1] = donor_species + '_highselect'
                Donor_species.setdefault(donor_species,
                                         donor_species_highselect)
                genus, species = genus_species(donor_species)
                line_set = lines.split('\n')[0].split('\t')
                if genus not in Genus_notHS:
                    genus_highselect = init_highselect_notHS(line_set)
                    genus_highselect[0] = genus
                    genus_highselect[1] = genus
                    Genus_notHS.setdefault(genus, genus_highselect)
                    line_set = lines.split('\n')[0].split('\t')
                    genus_highselect = init_highselect_empty(line_set)
                    genus_highselect[0] = genus
                    genus_highselect[1] = genus + '_highselect'
                    Genus.setdefault(genus, genus_highselect)
                else:
                    genus_highselect = Genus_notHS[genus]
                    Genus_notHS[genus] = add_highselect(line_set, genus_highselect)
                line_set = lines.split('\n')[0].split('\t')
                if species not in Species_notHS:
                    species_highselect = init_highselect_notHS(line_set)
                    species_highselect[0] = species
                    species_highselect[1] = species
                    Species_notHS.setdefault(species, species_highselect)
                    line_set = lines.split('\n')[0].split('\t')
                    species_highselect = init_highselect_empty(line_set)
                    species_highselect[0] = species
                    species_highselect[1] = species + '_highselect'
                    Species.setdefault(species, species_highselect)
                else:
                    species_highselect = Species_notHS[species]
                    Species_notHS[species] = add_highselect(line_set, species_highselect)
    # calculate NS for new HS genes
    for lines in open(input_summary, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        Selected = line_set[-1]
        donor_species = line_set[0]
        record_id = line_set[1]
        # species or genes
        if Selected == 'False':
            if record_id not in ['0', 'gene', 'allspecies', 'allspecies_highselect',
                                 'allspecies_core', 'allspecies_flexible',
                                 'allspecies_highselect_noNS']:
                genus, species = genus_species(donor_species)
                donor_species_set = donor_species.split('_')
                try:
                    donor_species = '%s_%s_%s' % (donor_species_set[0],
                                                  donor_species_set[1][0:min(6, len(donor_species_set[1]))],
                                                  donor_species_set[2][0:min(6, len(donor_species_set[2]))])
                except IndexError:
                    donor_species = '%s_%s_%s' % (donor_species_set[0],
                                                  donor_species_set[1],
                                                  donor_species_set[1])
                if 'cluster' in donor_species_set[-1]:
                    try:
                        donor_species += '_CL' + donor_species_set[-2].split('cluster')[1]
                    except IndexError:
                        donor_species += '_CL' + donor_species_set[-1].split('cluster')[1]
                record_id = '%s__C_%s_G_%s' % (donor_species, record_id.split('_')[1], record_id.split('_')[-1])
                # highly selected genes
                if record_id in High_select2:
                    if (line_set[10] == 'observe_N_only' or float(line_set[10]) > NSratioobserve_cutoff):
                        newoutput.append('\t'.join(line_set[:-1]) + '\tTrue\n')
                        allspecies_highselect = add_highselect(line_set, allspecies_highselect)
                        Donor_species[line_set[0]] = add_highselect(line_set, Donor_species[line_set[0]])
                        Genus[genus] = add_highselect(line_set, Genus[genus])
                        Species[species] = add_highselect(line_set, Species[species])
                    else:
                        newoutput.append('\t'.join(line_set[:-1]) + '\tFalse\n')
                        High_select2.remove(record_id)
                # other genes and species
                elif line_set[0] not in line_set[1]:
                    newoutput.append(lines)
                else:
                    print('not output %s' % (line_set[1]))
        # original highly selected
        elif Selected == 'True':
            newoutput.append(lines)
            allspecies_highselect = add_highselect(line_set, allspecies_highselect)
            Donor_species[line_set[0]] = add_highselect(line_set, Donor_species[line_set[0]])
            Genus[genus] = add_highselect(line_set, Genus[genus])
            Species[species] = add_highselect(line_set, Species[species])
        else:
            newoutput.append(lines)
    # output new highselect
    newoutput = output_highlight(allspecies_highselect, allspecies, newoutput)
    for donor_species in Donor_species:
        line_set = Donor_species_notHS[donor_species]
        line_set_sub = [int(line_set[-13]), int(line_set[-11]), int(line_set[-9]),
                        int(line_set[-7]), int(line_set[-5]), int(line_set[-3])]
        newoutput = output_highlight(Donor_species[donor_species],
                                     line_set_sub, newoutput)
        newoutput = output_highlight(Donor_species_notHS[donor_species],
                                     line_set_sub, newoutput)
    for genus in Genus:
        line_set = Genus_notHS[genus]
        line_set_sub = [int(line_set[-13]), int(line_set[-11]), int(line_set[-9]),
                        int(line_set[-7]), int(line_set[-5]), int(line_set[-3])]
        newoutput = output_highlight(Genus[genus],
                                     line_set_sub, newoutput)
        newoutput = output_highlight(Genus_notHS[genus],
                                     line_set_sub, newoutput)
    for species in Species:
        line_set = Species_notHS[species]
        line_set_sub = [int(line_set[-13]), int(line_set[-11]), int(line_set[-9]),
                        int(line_set[-7]), int(line_set[-5]), int(line_set[-3])]
        newoutput = output_highlight(Species[species],
                                     line_set_sub, newoutput)
        newoutput = output_highlight(Species_notHS[species],
                                     line_set_sub, newoutput)
    # output new summary
    f1 = open(input_summary + '.High_select2.txt', 'w')
    f1.write(''.join(newoutput))
    f1.close()
    return High_select2

def sum_gene(input_summary,High_select2, selected = 1):
    genelist = []
    genelist += High_select2
    for lines in open(input_summary,'r'):
        lines_set = lines.replace('\n','').split('\t')
        if not lines.startswith('#'):
            if selected == 0 or lines_set[-1] == 'True':
                # highly selected genes or all genes
                gene_name = lines_set[1]
                genelist.append(gene_name)
    print(len(genelist))
    genelist = list(set(genelist))
    print(len(genelist))
    Output = []
    output_set = []
    for record in SeqIO.parse(all_fasta_HS, 'fasta'):
        Output.append('>%s\n%s\n' % (str(record.id), str(record.seq)))
        output_set.append(str(record.id))
    for record in SeqIO.parse(all_fasta, 'fasta'):
        if str(record.id) in genelist and str(record.id) not in output_set:
            Output.append('>%s\n%s\n' % (str(record.id), str(record.seq)))
            output_set.append(str(record.id))
    f1 = open(all_fasta_HS + '.High_select2.faa' , 'w')
    f1.write(''.join(Output))
    f1.close()

def annotation(all_filter_gene_fasta_file,pre_cluster = ''):
    # run cluster
    cutoff = 0.7
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                   % ('usearch', all_filter_gene_fasta_file, cutoff, all_filter_gene_fasta_file,
                      all_filter_gene_fasta_file, 40))
    os.system(cmd_cluster)
    all_filter_gene_fasta_file = all_filter_gene_fasta_file + '.cluster.aa'
    if pre_cluster!= '':
        os.system('#%s -makeudb_usearch %s -output %s.udb' %
                  ('usearch', pre_cluster, pre_cluster))
        os.system('%s -ublast %s -db %s.udb  -evalue 1e-2 -accel 0.5 -blast6out %s -threads 2'%
                  ('usearch', all_filter_gene_fasta_file,pre_cluster, all_filter_gene_fasta_file + '.ref.out.txt'))
    # run prokka
    cmdsprokka = 'py37\nprokka --kingdom Bacteria --outdir %s/prokka_%s  --protein %s --locustag Bacter %s/%s\n' % \
                 (output_dir_merge + '/summary', os.path.split(all_filter_gene_fasta_file)[-1],
                  all_filter_gene_fasta_file,
                  output_dir_merge + '/summary',
                  os.path.split(all_filter_gene_fasta_file)[-1].replace('.faa', '.fna'))
    f1 = open(os.path.join(input_script_sub, 'prokka.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmdsprokka))
    f1.close()
    # run metacyc
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/mit_alm/database/metacyc/protseq.fsa'
    cmds = ("diamond blastp --query %s --db %s.dmnd --out %s.metacyc.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    f1 = open(os.path.join(input_script_sub, 'metacyc.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    # run eggnog
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/eggnog/xaa.hmm'
    cmds = ('hmmsearch --tblout %s.eggnog.1.txt --cpu 40 -E %s %s %s\n') %(all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.1.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xab.hmm'
    cmds = ('hmmsearch --tblout %s.eggnog.2.txt --cpu 40 -E %s %s %s\n') % (
    all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.2.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xac.hmm'
    cmds = ('hmmsearch --tblout %s.eggnog.3.txt --cpu 40 -E %s %s %s\n') % (
    all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.3.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmds))
    f1.close()
    # run kegg
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/kegg/kofam/profiles/prokaryote/prokaryote.hmm'
    cmds = ('hmmsearch --tblout %s.kegg.txt --cpu 40 -E %s %s %s\n') %(all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'kegg.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    # run customed database
    cutoff = 80
    cutoff2 = 80
    cmds = ''
    database = '/scratch/users/anniz44/scripts/database/SARG.db.fasta'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.SARG.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 50
    cutoff2 = 50
    database = '/scratch/users/anniz44/scripts/database/AHR.aa.db'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.AHR.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 60
    cutoff2 = 80
    database = '/scratch/users/anniz44/scripts/database/Butyrate.pro.aa'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.buty.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/anniz44/scripts/database/IntI1_database.fasta'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.int.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/anniz44/scripts/database/SRB.AA'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.SRB.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 0.01
    database = '/scratch/users/anniz44/scripts/database/NR.hmm'
    cmds += ('hmmsearch --tblout %s.NR.txt --cpu 40 -E %s %s %s\n') %(all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'customed.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    # all scripts
    f1 = open(os.path.join(input_script, 'allannotate.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.sh')):
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()

# clustering
Clusters_gene, High_select2_output,High_select2 = cluster_uc(all_fasta + '.uc')
f1 = open(all_fasta + '.High_select2.txt', 'w')
f1.write(''.join(High_select2_output))
f1.close()
# correcting
High_select2 = add_new_selection(input_summary,High_select2)

# run clustering
sum_gene(input_summary,High_select2,1)
annotation(all_fasta_HS + '.High_select2.faa',pre_cluster)

################################################### END ########################################################
################################################### SET PATH ########################################################
# CHANGED gene annotation summary
import os
import glob
import copy
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome/'
#input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly_jay/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
database_metacyc = '/scratch/users/mit_alm/database/metacyc/genes.col'
database_eggnog = '/scratch/users/mit_alm/database/eggnog/2_annotations.tsv'
database_eggnog2 = '/scratch/users/mit_alm/database/eggnog/COG_catogary.txt'
database_kegg = '/scratch/users/anniz44/scripts/database/Kegg/ko_formated'
#all_fasta_set = glob.glob(os.path.join(output_dir + '/summary_jay', '*.faa'))
all_fasta_set = glob.glob(os.path.join(output_dir + '/summary/', '*.faa'))

# function annotation
def best_hit(Blast_output,small = 0):
    for gene_name in Blast_output:
        db_name,identity = Blast_output[gene_name]
        if len(identity) > 2:
            identity2 = copy.deepcopy(identity)
            if small == 1:
                identity2.sort()
            else:
                identity2.sort(reverse=True)
            top1 = identity2[0]
            top2 = identity2[1]
            Blast_output[gene_name]=[db_name[identity.index(top1)],
                                     db_name[identity.index(top2)]]
        else:
            Blast_output[gene_name]=db_name
    return Blast_output

def annotate_metacyc(blast_search):
    Blast_output = dict()
    DB_name = dict()
    for lines in open(blast_search, 'r'):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            gene_name = lines_set[0]
            db_name = lines_set[1].split('|')[-1].split('-MONOMER')[0]
            identity = float(lines_set[2])
            Blast_output.setdefault(gene_name, [[],[]])
            Blast_output[gene_name][0].append(db_name)
            Blast_output[gene_name][1].append(identity)
            DB_name.setdefault(db_name, ['',''])
    for database in [database_metacyc]:
        for lines in open(database, 'r', encoding="utf8", errors='ignore'):
            if not lines.startswith('#'):
                lines_set = lines.split('\n')[0].split('\t')
                db_name = lines_set[0].split('-MONOMER')[0]
                function_name = lines_set[1]
                annotation_fun = lines_set[2]
                if db_name in DB_name:
                    DB_name[db_name][0] = function_name
                    DB_name[db_name][1] = annotation_fun
    best_hit(Blast_output, 0)
    return [Blast_output,DB_name]

def annotate_eggnog(blast_search):
    Anno_eggnog = dict()
    for lines in open(database_eggnog2,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        Anno_eggnog.setdefault(lines_set[0],lines.split('\n')[0])
    Blast_output = dict()
    DB_name = dict()
    DB_name_2 = dict()
    for blast_search_file in blast_search:
        for lines in open(blast_search_file, 'r'):
            if not lines.startswith('#'):
                db_name = ''
                identity = 0
                lines_set = lines.split('\n')[0].split(' ')
                gene_name = lines_set[0]
                for sub_line in lines_set[1:]:
                    if sub_line !='' and sub_line !='-':
                        if db_name == '':
                            db_name = sub_line.split('.')[0]
                        elif identity == 0:
                            identity = float(sub_line)
                            break
                Blast_output.setdefault(gene_name, [[], []])
                Blast_output[gene_name][0].append(db_name)
                Blast_output[gene_name][1].append(identity)
                DB_name.setdefault(db_name, ['', ''])
                DB_name_2.setdefault(db_name, '')
    for database in [database_eggnog]:
        for lines in open(database, 'r'):
            if not lines.startswith('#'):
                lines_set = lines.split('\n')[0].split('\t')
                db_name = lines_set[1]
                function_name = ''
                annotation_catogary = ''
                for COG in lines_set[2]:
                    annotation_catogary += '\t' + Anno_eggnog[COG]
                annotation_fun = lines_set[3]
                if db_name in DB_name:
                    DB_name[db_name][0] = function_name
                    DB_name[db_name][1] = annotation_fun
                    DB_name_2[db_name] = '%s\t%s%s'%(function_name,annotation_fun,annotation_catogary)
    best_hit(Blast_output, 1)
    return [Blast_output,DB_name,DB_name_2]

def annotate_kegg(blast_search):
    Blast_output = dict()
    DB_name = dict()
    DB_name2 = dict()
    for lines in open(blast_search, 'r'):
        if not lines.startswith('#'):
            db_name = ''
            identity = 0
            lines_set = lines.split('\n')[0].split(' ')
            gene_name = lines_set[0]
            for sub_line in lines_set[1:]:
                if sub_line != '' and sub_line != '-':
                    if db_name == '':
                        db_name = sub_line
                    elif identity == 0:
                        identity = float(sub_line)
                        break
            Blast_output.setdefault(gene_name, [[], []])
            Blast_output[gene_name][0].append(db_name)
            Blast_output[gene_name][1].append(identity)
            DB_name.setdefault(db_name, ['', ''])
            DB_name2.setdefault(db_name, '')
    for database in [database_kegg]:
        for lines in open(database, 'r'):
            if not lines.startswith('#'):
                lines_set = lines.split('\n')[0].split('\t')
                db_name = lines_set[0]
                function_name = lines_set[-6]
                annotation_fun = lines_set[-5]
                if db_name in DB_name:
                    DB_name[db_name][0] = function_name
                    DB_name[db_name][1] = annotation_fun
                    if 'Drug resistance' in lines_set[5] or ('Cellular community - eukaryotes' not in lines_set[5] and 'Overview' not in lines_set[5] \
                            and 'Human Diseases' not in lines_set[4] and 'Organismal Systems' not in lines_set[4]):
                        DB_name2[db_name] = '\t'.join(lines_set[4:])
    best_hit(Blast_output, 1)
    return [Blast_output,DB_name,DB_name2]

def annotate_prokka(fasta_input):
    cutoff = 80
    cutoff2 = 80
    cmds = ("diamond makedb --in %s -d %s.dmnd\ndiamond blastp --query %s --db %s.dmnd --out %s.prokka.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                % (database_prokka_fasta,database_prokka_fasta,
                   fasta_input, database_prokka_fasta, fasta_input, cutoff, cutoff2))
    os.system(cmds)
    blast_search = '%s.prokka.txt' %(fasta_input)
    Blast_output = dict()
    DB_name = dict()
    for lines in open(blast_search, 'r'):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            gene_name = lines_set[0]
            db_name = lines_set[1]
            identity = float(lines_set[2])
            Blast_output.setdefault(gene_name, [[], []])
            Blast_output[gene_name][0].append(db_name)
            Blast_output[gene_name][1].append(identity)
            DB_name.setdefault(db_name, ['', ''])
    for database in [database_prokka]:
        for lines in open(database, 'r'):
            if not lines.startswith('#'):
                lines_set = lines.split('\n')[0].split('\t')
                db_name = lines_set[0]
                function_name = lines_set[3].split('_')[0]
                annotation_fun = lines_set[6]
                if db_name in DB_name:
                    DB_name[db_name][0] = function_name
                    DB_name[db_name][1] = annotation_fun
    best_hit(Blast_output, 1)
    return [Blast_output,DB_name]

def cluster_uc(cluster_input):
    Clusters = dict()
    for lines in open(cluster_input, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        cluster = line_set[1]
        record_name = line_set[8]
        Clusters.setdefault(record_name, cluster)
    return Clusters

def annotate_gene(Blast_output1,DB_name1,All_annotation):
    for gene_name in Blast_output1:
        All_annotation.setdefault(gene_name, [[], set(), set(),[], [], []])
        for db_name in Blast_output1[gene_name]:
            if db_name in DB_name1:
                fun, anno = DB_name1[db_name]
                All_annotation[gene_name][1].add(fun)
                All_annotation[gene_name][2].add(anno)
    return All_annotation

def sum_kegg_eggnog(Blast_outputkegg,DB_namekegg,Blast_outputegg,DB_nameegg,Clusters):
    Output = []
    for gene_name in Blast_outputkegg:
        cluster = Clusters.get(gene_name, '')
        donor = gene_name.split('_')[0]
        donor_species = gene_name.split('__')[0].split('_CL')[0].replace('BA_BA', 'Bifido_adoles').replace('BL_BL',
                                                                                                           'Bifido_longum').replace(
            'PB_PB', 'Paraba_butyra').replace('BA_cluste', 'Bifido_adoles').replace('BL_cluste',
                                                                                    'Bifido_longum').replace(
            'PB_cluste', 'Paraba_butyra')
        species = donor_species.replace(donor + '_', '')
        for db_name in Blast_outputkegg[gene_name]:
            if db_name in DB_namekegg:
                Output.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(cluster,donor,species,donor_species,gene_name,db_name,DB_namekegg[db_name]))
    f1 = open(all_fasta + '.cluster.aa.kegg.sum','w')
    f1.write('cluster\tdonor\tspecies\tdonor_species\tgene_name\tKO\tBRITE_KO1\tBRITE_KO2\tBRITE_KO3\n' + ''.join(Output))
    f1.close()
    Output = []
    for gene_name in Blast_outputegg:
        cluster = Clusters.get(gene_name, '')
        donor = gene_name.split('_')[0]
        donor_species = gene_name.split('__')[0].split('_CL')[0].replace('BA_BA', 'Bifido_adoles').replace('BL_BL',
                                                                                                           'Bifido_longum').replace(
            'PB_PB', 'Paraba_butyra').replace('BA_cluste', 'Bifido_adoles').replace('BL_cluste',
                                                                                    'Bifido_longum').replace(
            'PB_cluste', 'Paraba_butyra')
        species = donor_species.replace(donor + '_', '')
        for db_name in Blast_outputegg[gene_name]:
            if db_name in DB_nameegg:
                Output.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                cluster, donor, species, donor_species, gene_name, db_name, DB_nameegg[db_name]))
    f1 = open(all_fasta + '.cluster.aa.eggnog.sum', 'w')
    f1.write(
        'cluster\tdonor\tspecies\tdonor_species\tgene_name\tEGGNOG\tannotation\tCOG\tCOG1\tCOG2\n' + ''.join(Output))
    f1.close()

def cluster_genes(All_annotation,Clusters):
    for gene_name in Clusters:
        cluster = Clusters.get(gene_name,'')
        donor = gene_name.split('_')[0]
        donor_species = gene_name.split('__')[0].split('_CL')[0].replace('BA_BA','Bifido_adoles').replace('BL_BL','Bifido_longum').replace('PB_PB','Paraba_butyra').replace('BA_cluste','Bifido_adoles').replace('BL_cluste','Bifido_longum').replace('PB_cluste','Paraba_butyra')
        species = donor_species.replace(donor + '_', '')
        All_annotation.setdefault(gene_name,[[], set(), set(),[], [], []])
        All_annotation[gene_name][3].append(donor)
        All_annotation[gene_name][4].append(species)
        All_annotation[gene_name][5].append(donor_species)
        if cluster!= []:
            All_annotation[gene_name][0].append(cluster)
        else:
            print('missing cluster for gene %s'%(gene_name))
    return All_annotation

def cluster_species(All_annotation):
    Clusters = dict()
    Species = dict()
    Donorspecies = dict()
    Donor = dict()
    output_cluster = []
    Species_species_count = dict()
    output_sum = []
    output_cluster.append('cluster\tTag_species\tTag_donor\tdonor\tspecies\tdonor_species\tfun\tanno\t\n')
    output_sum.append('gene_name\tcluster\tdonor\tspecies\tdonor_species\tfun\tanno\t\n')
    for gene_name in All_annotation:
        clusters, fun, anno, donor, species, donor_species = All_annotation[gene_name]
        All_annotation[gene_name][0] = ';'.join(list(set(clusters)))
        for cluster in clusters:
            Clusters.setdefault(cluster, [[], [], [], set(), set()])
            Clusters[cluster][0].append(donor[0])
            Clusters[cluster][1].append(species[0])
            Clusters[cluster][2].append(donor_species[0])
            All_annotation[gene_name][1] = set()
            All_annotation[gene_name][2] = set()
            for a_fun in fun:
                if a_fun != '':
                    All_annotation[gene_name][1].add(a_fun)
                    Clusters[cluster][3].add(a_fun)
            for a_anno in anno:
                if a_anno != '':
                    All_annotation[gene_name][2].add(a_anno)
                    Clusters[cluster][4].add(a_anno)
        cluster, fun, anno, donor, species, donor_species = All_annotation[gene_name]
        output_sum.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n' % (gene_name,
                                                              cluster, ';'.join(donor), ';'.join(species),
                                                              ';'.join(donor_species), ';'.join(fun), ';'.join(anno)
                                                              ))
    Clusters_sum = dict()
    for cluster in Clusters:
        Tag_species = 'unique_species'
        Tag_donor = 'unique_donor'
        species_set = set()
        genus_set = set()
        alldonor, allspecies, alldonorspecies, fun, anno = Clusters[cluster]
        for donor in alldonor:
            Donor.setdefault(donor, [])
            Donor[donor].append(cluster)
        for donorspecies in alldonorspecies:
            Donorspecies.setdefault(donorspecies, [])
            Donorspecies[donorspecies].append(cluster)
        for species in allspecies:
            new_species = species.split('_cluste')[0].split('_newcluster')[0].split('_CL')[0]
            species_set.add(new_species)
            Species.setdefault(new_species, [])
            Species[new_species].append(cluster)
            new_genus = new_species.split('_')[0]
            genus_set.add(new_genus)
            Species_species_count.setdefault(new_species, [])
        if len(allspecies) > 1:
            if len(alldonor) > 1:
                Tag_donor = 'multi_donor'
            if len(species_set) > 1:
                Tag_species = 'multi_species'
            if len(genus_set) > 1:
                Tag_species = 'multi_genera'
        Clusters_sum.setdefault(cluster, [Tag_species, Tag_donor, ';'.join(list(fun)), ';'.join(list(anno))])
        Tag_species, Tag_donor, fun, anno = Clusters_sum[cluster]
        output_cluster.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n' % (cluster, Tag_species, Tag_donor,
                                                         ';'.join(alldonor), ';'.join(list(species_set)),';'.join(alldonorspecies),
                                                          fun, anno))
        for new_species in species_set:
            for new_species2 in species_set:
                if new_species < new_species2:
                    Species_species_count[new_species].append(new_species2)
                    Species_species_count[new_species2].append(new_species)
    Species_sum = []
    Species_sum.append('#species\tgenus\tcluster\tTag_species\tTag_donor\t\n')
    for new_species in Species:
        new_genus = new_species.split('_')[0]
        clusters = list(set(Species[new_species]))
        for cluster in clusters:
            Species_sum.append('%s\t%s\t%s\t%s\t%s\t\n' % (
                new_species, new_genus, cluster, Clusters_sum[cluster][0], Clusters_sum[cluster][1]))
    Donor_sum = []
    Donor_sum.append('#donor\tcluster\tTag_species\tTag_donor\t\n')
    for donor in Donor:
        clusters = list(set(Donor[donor]))
        for cluster in clusters:
            Donor_sum.append(
                '%s\t%s\t%s\t%s\t\n' % (donor, cluster, Clusters_sum[cluster][0], Clusters_sum[cluster][1]))
    Donorspecies_sum = []
    Donorspecies_sum.append('#donorspecies\tcluster\tTag_species\tTag_donor\t\n')
    for donorspecies in Donorspecies:
        clusters = list(set(Donorspecies[donorspecies]))
        for cluster in clusters:
            Donorspecies_sum.append(
                '%s\t%s\t%s\t%s\t\n' % (donorspecies, cluster, Clusters_sum[cluster][0], Clusters_sum[cluster][1]))
    f1 = open(all_fasta + '.gene.sum', 'w')
    f1.write(''.join(output_sum))
    f1.close()
    f1 = open(all_fasta + '.cluster.sum', 'w')
    f1.write(''.join(output_cluster))
    f1.close()
    f1 = open(all_fasta + '.species.sum', 'w')
    f1.write(''.join(Species_sum))
    f1.close()
    f1 = open(all_fasta + '.donor.sum', 'w')
    f1.write(''.join(Donor_sum))
    f1.close()
    f1 = open(all_fasta + '.donor.species.sum', 'w')
    f1.write(''.join(Donorspecies_sum))
    f1.close()
    Output = []
    Output.append('#species1\tspecies2\tgenus1\tgenus2\tnum_cluster_shared\t\n')
    for new_species in Species_species_count:
        all_species = Species_species_count[new_species]
        all_species_unique = set(all_species)
        new_genus1 = new_species.split('_')[0]
        for new_species2 in all_species_unique:
            new_genus2 = new_species2.split('_')[0]
            Output.append(
                '%s\t%s\t%s\t%s\t%s\t\n' % (new_species, new_species2, new_genus1,new_genus2, all_species.count(new_species2)))
    f1 = open(all_fasta + '.speciesTospecies.sum', 'w')
    f1.write(''.join(Output))
    f1.close()

def sum_annotation(Blast_output1,DB_name1,Blast_output2,DB_name2,Blast_output3,DB_name3,Blast_output4,DB_name4,Clusters_gene):
    All_annotation = dict()
    All_annotation = annotate_gene(Blast_output1, DB_name1, All_annotation)
    All_annotation = annotate_gene(Blast_output2, DB_name2, All_annotation)
    All_annotation = annotate_gene(Blast_output3, DB_name3, All_annotation)
    All_annotation = annotate_gene(Blast_output4, DB_name4, All_annotation)
    All_annotation = cluster_genes(All_annotation, Clusters_gene)
    cluster_species(All_annotation)

for all_fasta in all_fasta_set:
    all_fasta_folder, all_fasta_file = os.path.split(all_fasta)
    try:
        database_prokka = glob.glob(os.path.join(all_fasta_folder, 'prokka_%s/*.tsv' % (all_fasta_file)))[0]
        database_prokka_fasta = glob.glob(os.path.join(all_fasta_folder, 'prokka_%s/*.faa' % (all_fasta_file)))[0]
        Blast_output4, DB_name4 = annotate_prokka(all_fasta)
    except IndexError:
        Blast_output4 = dict()
        DB_name4 = dict()
    # set up gene annotation and clustering
    Blast_output1, DB_name1, DB_name1_2 = annotate_kegg(all_fasta + '.cluster.aa.kegg.txt')
    Blast_output2, DB_name2 = annotate_metacyc(all_fasta + '.cluster.aa.metacyc.txt')
    Blast_output3, DB_name3, DB_name3_2 = annotate_eggnog(glob.glob(all_fasta + '.cluster.aa.eggnog.*.txt'))
    # sum up
    Clusters_gene = cluster_uc(all_fasta + '.uc')
    sum_kegg_eggnog(Blast_output1,DB_name1_2,Blast_output3,DB_name3_2,Clusters_gene)
    sum_annotation(Blast_output1, DB_name1, Blast_output2, DB_name2, Blast_output3, DB_name3, Blast_output4, DB_name4,
                   Clusters_gene)

################################################### END ########################################################
################################################### SET PATH ########################################################
# simulate adaptive genes
# calculate gene number for each clonal population
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics

# set up path -> selected
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round*/*')+\
glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/round*/*')
output = os.path.join(input_script,'clonal_pop_gene_num.txt')

for folder in genome_dir:
    allfasta = glob.glob(os.path.join(folder,'*.fasta.corrected.fasta.faa'))
    if allfasta!= []:
        allfasta = allfasta[0]
        os.system('echo %s:: >> %s'%(os.path.split(folder)[-1],output))
        os.system('grep \">\" %s | wc -l >> %s'%(allfasta,output))

################################################### END ########################################################
################################################### SET PATH ########################################################
# use annotated genes
import glob
allanno = glob.glob('*.sum')

# load annotated genes
annotated = dict()
for lines in open('annotated.genes.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    try:
        geneset = lines_set[-2].split(';')
        for gene in geneset:
            annotated.setdefault(gene.replace(' ',''),lines.split('\n')[0])
    except IndexError:
        pass

# output annotation
confirmed = ['lolD']
for anno in allanno:
    Output = []
    for lines in open(anno,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        geneset = lines_set[6].replace(',',';').split(';')
        funset = lines_set[7].replace(',',';').split(';')
        temp_line = ''
        for gene in geneset:
            gene = gene.replace(' ','')
            if gene in annotated and gene != '':
                newline = annotated[gene]
                funnew = newline.split('\t')[-1]
                if any(fun in funnew for fun in funset) or gene in confirmed:
                    if temp_line == '':
                        temp_line = newline
                    elif temp_line!= newline:
                        if newline.split('\t')[0] != temp_line.split('\t')[0]:
                            Output.append('%s\t%s\n' % (lines.split('\n')[0], temp_line))
                            temp_line = newline
                        elif len(newline) > len(temp_line):
                            temp_line = newline
                else:
                    print(gene,funset)
                    print(funnew + '\n')
        Output.append('%s\t%s\n'%(lines.split('\n')[0],temp_line))
    f1=open(anno + '.new.txt','w')
    f1.write(''.join(Output))
    f1.close()

# add new kegg ko -> rarely used
Output = set()
temp = ['','','']
Not_output = 0
for lines in open('kegg_new','r'):
    lines_set = lines.split('\n')[0].lstrip().split(' ')
    KO = lines_set[0]
    Brite = ' '.join(lines_set)
    if not KO.startswith('K'):
        if Not_output == 0:
            temp[-1] = Brite
        elif Not_output == 1:
            temp[-2] = temp[-1]
            temp[-1] = Brite
        elif Not_output == 2:
            temp[-3] = temp[-2]
            temp[-2] = temp[-1]
            temp[-1] = Brite
        Not_output += 1
    else:
        print(temp)
        Not_output = 0
        Output.add('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(
        KO,' '.join(lines_set[1:]).split(';')[0],
        ' '.join(lines_set[1:]),'',' '.join(temp[-3].split(' ')[1:]),
        ' '.join(temp[-2].split(' ')[1:]),temp[-1].split(' [')[0]))

f1 = open('ko_format_new','w')
f1.write(''.join(list(Output)))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# core genome or flexible genome
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/pangenome'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
mutation_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome/summary'
mutation_fasta = os.path.join(mutation_dir,'all.selected.gene.faa.High_select2.faa')
mutation_fasta2 = os.path.join(mutation_dir,'all.denovo.gene.faa')
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/round*/*')+\
             glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round*/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome/pangenome'

try:
    os.mkdir(output_dir)
except IOError:
    pass

os.system('rm -rf %s'%(input_script_sub))
try:
    os.mkdir(input_script_sub)
except IOError:
    pass

cutoff = 50
cutoff2 = 80
genome_num = 0
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor_species_set = donor_species.split('_')
    donor = donor_species_set[0]
    cmds = ''
    allfasta = glob.glob(os.path.join(folder,'*.corrected.fasta'))
    for fasta in allfasta:
        fastaname = os.path.split(fasta)[-1]
        try:
            f1 = open(fasta + '.faa', 'r')
        except FileNotFoundError:
            os.system('prodigal -q -i %s -a %s' % (fasta, fasta + '.faa'))
        fasta = fasta + '.faa'
        cmds += (
                    "diamond blastp --query %s --db %s.dmnd --out %s.denovo.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                    % (fasta, mutation_fasta2, os.path.join(output_dir, fastaname), cutoff, cutoff2))
    if cmds!= '':
        f1 = open(os.path.join(input_script_sub, '%s.denovo.sh' % donor_species), 'a')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
        f1.close()

f1 = open(os.path.join(input_script, 'pangenome.denovo.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.denovo.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# summary core genome or flexible genome
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/pangenome'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
mutation_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome/summary'
mutation_select_fasta = os.path.join(mutation_dir,'all.selected.gene.faa.High_select2.faa')
mutation_fasta = os.path.join(mutation_dir,'all.denovo.gene.faa')
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/round*/*')+\
             glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round*/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome/pangenome'
allresult = glob.glob(os.path.join(output_dir,'*.denovo.txt'))
core_cutoff = 0.8
# load diamond results
def load_result(file,speciesname,genomename,Gene,Gene_copy,alloutput):
    resultset = []
    tempspeciesset = dict()
    for lines in open(file,'r'):
        lines_set = lines.split('\t')
        geneID = lines_set[1]
        resultset.append(geneID)
        tempspeciesset.setdefault(geneID,set())
        tempspeciesset[geneID].add(speciesname)
        alloutput.append('%s\t%s\t%s\t\n'%(geneID,genomename,speciesname))
        Gene_copy[geneID].append(speciesname)
    for geneID in tempspeciesset:
        Gene[geneID]+=list(tempspeciesset[geneID])
    return [Gene,Gene_copy,alloutput]

# selected genes
# load selected genes
Gene_select = []
for record in SeqIO.parse(mutation_select_fasta, 'fasta'):
    geneID = str(record.id)
    Gene_select.append(geneID)

# load all genes
Gene = dict()
Gene_copy = dict()
Gene_summary = dict()
Gene_length = dict()
for record in SeqIO.parse(mutation_fasta, 'fasta'):
    geneID = str(record.id)
    gene_length = len(str(record.seq))
    Gene.setdefault(geneID,[])
    Gene_copy.setdefault(geneID, [])
    Gene_summary.setdefault(geneID, set())
    Gene_length.setdefault(geneID,gene_length)

# check all genes
for geneID in Gene_select:
    if geneID not in Gene:
        print('missing geneID %s in %s'%(geneID,mutation_fasta))

# summarize all genomes for a species
Species = dict()
try:
    f1 = open(mutation_fasta + '.allpangenome.txt','r')
    Gene_temp = dict()
    for lines in f1:
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            geneID,genomename,speciesname = lines_set[0:3]
            if '_g' in speciesname and '.denovo.txt' in speciesname: # fix Parabacteroides
                speciesname = speciesname.split('_')[1]
            Species.setdefault(speciesname, set())
            Species[speciesname].add(genomename)
            Gene_temp.setdefault(geneID, set())
            Gene_temp[geneID].add(genomename)
            Gene_copy[geneID].append(speciesname)
    for geneID in Gene_temp:
        for genomename in Gene_temp[geneID]:
            genomename_set = genomename.split('_')
            speciesname = genomename_set[1] + '_' + genomename_set[2]
            if genomename_set[1] in ['BA', 'BL', 'PB']:
                if genomename_set[1] == 'BA':
                    speciesname = 'Bifidobacterium_adolescentis'
                elif genomename_set[1] == 'BL':
                    speciesname = 'Bifidobacterium_longum'
                elif genomename_set[1] == 'PB':
                    speciesname = 'Parabacteroides_butyrate'
            if genomename_set[2].startswith('_g'):
                speciesname = genomename_set[1]
            Gene[geneID].append(speciesname)
except IOError:
    alloutput = []
    alloutput.append('#geneID\tgenome\tspecies\t\n')
    for resultfile in allresult:
        genomename = os.path.split(resultfile)[-1].split('.faa.selected.txt')[0]
        genomename_set = genomename.split('_')
        speciesname = genomename_set[1]+'_'+genomename_set[2]
        if genomename_set[1] in ['BA','BL','PB']:
            if genomename_set[1] == 'BA':
                speciesname = 'Bifidobacterium_adolescentis'
            elif genomename_set[1] == 'BL':
                speciesname = 'Bifidobacterium_longum'
            elif genomename_set[1] == 'PB':
                speciesname = 'Parabacteroides_butyrate'
        if genomename_set[2].startswith('_g'):
            speciesname = genomename_set[1]
        Species.setdefault(speciesname,set())
        Species[speciesname].add(genomename)
        Gene,Gene_copy,alloutput = load_result(resultfile,speciesname,genomename,Gene,Gene_copy,alloutput)
    f1 = open(mutation_fasta + '.allpangenome.txt','w')
    f1.write(''.join(alloutput))
    f1.close()

# summarize gene distribution in a species
#allsummary = []
allsummary2 = set()
allsummary2.add('#geneID\tspecies\tgene_length\tgenome_percentage\tavg_copy_num\tselected\t\n')
geneID_list = 'Species\tTotalgenome'
for geneID in Gene:
    geneID_list += '\t%s'%(geneID)

geneID_list += '\n'
#allsummary.append(geneID_list)
for speciesname in Species:
    Totalgenome = len(Species[speciesname])
    geneID_list = '%s\t%s'%(speciesname,Totalgenome)
    for geneID in Gene:
        Totalgenome_geneID = Gene[geneID].count(speciesname)
        Totalcopy_geneID = Gene_copy[geneID].count(speciesname)
        Percentage = Totalgenome_geneID/Totalgenome
        Percentage_copy = Totalcopy_geneID/Totalgenome
        geneID_list += '\t%.3f' % (Percentage)
        if Percentage >= core_cutoff:
            Gene_summary[geneID].add('core')
        elif Percentage > 0:
                Gene_summary[geneID].add('flexible')
        if geneID in Gene_select:
            tag = 'selected'
        else:
            tag = 'nonselected'
        if Percentage > 0:
            allsummary2.add('%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n'%(geneID,speciesname,Gene_length[geneID],Totalgenome_geneID,
                                                          Percentage,Percentage_copy,tag))
    geneID_list += '\n'
    #allsummary.append(geneID_list)

Gene_summary_output=[]
for geneID in Gene:
    allspecies = Gene[geneID]
    allspecies_unique = list(set(allspecies))
    Totalspecies = len(allspecies_unique)
    tempoutput = ('%s\t%s\t'%(geneID,Totalspecies))
    if geneID in Gene_select:
        tempoutput += 'selected\t'
    else:
        tempoutput += 'nonselected\t'
    if Totalspecies > 1:
        tempoutput += ('multispecies\t')
    else:
        tempoutput += ('singlespecies\t')
    tempoutput += ';'.join(allspecies_unique) + '\t'
    for tag in Gene_summary[geneID]:
        tempoutput += '%s\t'%(tag)
    Gene_summary_output.append(tempoutput)

# output results
#f1 = open(mutation_fasta + '.allpangenome.sum.2.txt','w')
#f1.write(''.join(allsummary))
#f1.close()
f1 = open(mutation_fasta + '.allpangenome.sum.txt','w')
f1.write('\n'.join(Gene_summary_output))
f1.close()
f1 = open(mutation_fasta + '.allpangenome.sum.3.txt','w')
f1.write(''.join(list(allsummary2)))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# annitator results
annotation = dict()
for lines in open('newannotation','r'):
    lines_set = lines.split('\n')[0].split('\t')
    annotation.setdefault(lines_set[0],lines)

# merge multiple gene names
Output = []
annotation_output = []
for lines in open('all.selected.gene.faa.High_select2.faa.cluster.sum.new.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    fun = lines_set[1].replace(' ','').replace(',',';')
    fun_set = fun.split(';')
    temp_line = ''
    for fun in fun_set:
        if fun in annotation:
            annotation_output.append(fun)
            temp_line += '%s\t%s'%(lines_set[1],annotation[fun])
    if temp_line != '':
        Output.append(temp_line + '\n')

for lines in open('all.trunc.gene.faa.cluster.sum.new.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    fun = lines_set[1].replace(' ','').replace(',',';')
    fun_set = fun.split(';')
    temp_line = ''
    for fun in fun_set:
        if fun in annotation:
            annotation_output.append(fun)
            temp_line += '%s\t%s'%(lines_set[1],annotation[fun])
    if temp_line != '':
        Output.append(temp_line + '\n')

notoutput = [fun for fun in annotation if fun not in annotation_output]
print(notoutput)

f1 = open('newannotation.order.txt','w')
f1.write(''.join(list(Output)))
f1.close()

# merge multiple gene names by uniprot
annotation = dict()
Output = []
for lines in open('newannotation.order.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    geneset = lines_set[0]
    gene = lines_set[2]
    protein = lines_set[3]
    geneprotein = '%s\t%s'%(gene,protein)
    if gene!='None':
        if geneprotein not in annotation:
            annotation.setdefault(geneprotein,lines.split('\n')[0])
        else:
            if annotation[geneprotein].split('\t')[0] != geneset:
                if lines_set[1] != annotation[geneprotein].split('\t')[1]:
                    annotation[geneprotein] = geneset + ';' + annotation[geneprotein] + '\t' + lines.split('\n')[0]
                else:
                    annotation[geneprotein] = geneset + ';' + annotation[geneprotein]
    else:
        Output.append('\t%s\n' % (lines.split('\n')[0]))

Output = []
for geneprotein in annotation:
    Output.append('%s\t%s\n'%(geneprotein,annotation[geneprotein]))

f1 = open('newannotation.order.txt','w')
f1.write(''.join(list(Output)))
f1.close()

# delete annotated genes
Output = []
annotated = dict()
notright = ['lacI','abfA','asa1','recG','tetM','cbpA']
for lines in open('annotated.genes.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    try:
        geneset = lines_set[-2].replace(' ','').replace(',',';').split(';')
        for gene in geneset:
            annotated.setdefault(gene.replace(' ',''),lines.split('\n')[0])
    except IndexError:
        pass

for lines in open('newannotation.order.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    gene = lines_set[2]
    geneset = gene.replace(' ','').replace(',',';').split(';')
    alreadyin = False
    for gene in geneset:
        if gene in annotated and gene not in notright:
            print(gene, annotated[gene])
            print(geneset, lines)
            alreadyin = True
            break
    if not alreadyin:
        Output.append(lines)

f1 = open('newannotation.order.2.txt','w')
f1.write(''.join(list(Output)))
f1.close()

#>> annotate and update annotated.genes.txt
# use annotated genes
import glob
allanno = glob.glob('*.sum')

# load annotated genes
annotated = dict()
for lines in open('annotated.genes.txt','r'):
    lines_set = lines.split('\n')[0].split('\t')
    try:
        geneset = lines_set[-2].split(';')
        for gene in geneset:
            annotated.setdefault(gene.replace(' ',''),lines.split('\n')[0])
    except IndexError:
        pass

# output annotation
confirmed = ['lolD','FAEPRAA2165_RS06535','FAEPRAA2165_RS06280','crtI','cbpA','valS','pepO']
Output = []
for anno in allanno[3:4]:
    tag = 'annotation'
    if 'High_select' in anno:
        tag = 'High_select'
    elif 'trunc' in anno:
        tag = 'Trunc'
    for lines in open(anno,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        geneset = lines_set[6].replace(',',';').split(';')
        funset = lines_set[7].replace(',',';').split(';')
        temp_line = ''
        for gene in geneset:
            gene = gene.replace(' ','')
            if gene in annotated and gene != '':
                newline = annotated[gene]
                funnew = newline.split('\t')[-1]
                if any(fun in funnew for fun in funset) or gene in confirmed:
                    if temp_line == '':
                        temp_line = newline
                    elif temp_line!= newline:
                        if newline.split('\t')[0] != temp_line.split('\t')[0]:
                            Output.append('%s\t%s\t%s\n' % ('\t'.join(lines_set[0:8]), tag, temp_line))
                            temp_line = newline
                        elif len(newline) > len(temp_line):
                            temp_line = newline
                else:
                    print(gene,funset)
                    print(funnew + '\n')
        Output.append('%s\t%s\t%s\n'%('\t'.join(lines_set[0:8]),tag,temp_line))

f1=open(anno + '.all.new.txt','w')
f1.write(''.join(Output))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# metagenomes mapping to a reference genome
import glob
import os

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/MGvcf'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round2/am_*')
metagenome_dir = '/scratch/users/anniz44/Metagenomes/public_metagenomes/fecal_human'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/MG'
fastq_name = '_1.fasta'
genus_select = ['am_Bifidobacterium','am_Bacteroides',
                'am_Collinsella', 'am_Escherichia_coli',
                'am_Turicibacter']

try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir + '/merge')
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
    os.mkdir(input_script)
except IOError:
    pass

os.system('rm -rf %s' %(input_script_sub))

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(input_script_merge)
except IOError:
    pass

allfastq = glob.glob(os.path.join(metagenome_dir,'am*%s'%(fastq_name)))
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    if donor == 'am' and not any(genus in donor_species for genus in genus_select):
        # BN10 donors and selected genus
        donor_species_genome = glob.glob(os.path.join(folder, '*.all.spades*.fasta'))
        sub_samples = []
        if donor_species_genome!= []:
            cmds = ''
            i = 0
            database = donor_species_genome[0]
            try:
                f1 = open(database + '.mmi')
            except IOError:
                os.system('minimap2 -d %s.mmi %s \n' % (database, database))
            for files in allfastq:
                donor_species_dir_file = os.path.split(files)[-1]
                if donor_species_dir_file.startswith('am'):
                    print(donor_species_dir_file)
                    tempbamoutput = os.path.join(output_dir + '/bwa/0', donor_species_dir_file + '.' + donor_species)
                    try:
                        f1 = open(tempbamoutput + '.sorted.bam')
                    except IOError:
                        files2 = files.replace(fastq_name, fastq_name.replace('1', '2'))
                        if 'fasta' in fastq_name:
                            cmds += 'minimap2' + ' -ax sr -N 1000 -p 0.8 -t %s %s.mmi %s %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                                min(40, 40), database, files, files2, 'samtools', min(40, 40),
                                tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools',
                                min(40, 40),
                                tempbamoutput)
                        else:
                            cmds += 'minimap2' + ' -ax sr -N 1000 -p 0.8 -t %s %s.mmi %s %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                                min(40, 40), database, files, files2, 'samtools', min(40, 40),
                                tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools',
                                min(40, 40),
                                tempbamoutput)
                        cmds += 'rm -rf %s.bam\n' % (tempbamoutput)
                        sub_samples.append(tempbamoutput + '.sorted.bam')
                        cmds += 'rm -rf %s.*.unaligned\n' % (tempbamoutput)
                        i += 1
                        if i % 25 == 0:
                            f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species, int(i / 25))), 'a')
                            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                            f1.close()
                            cmds = ''
            f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species, int(i / 25))), 'a')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            f1.close()
            # mileup
            cmds = ''
            tempbamoutput = os.path.join(output_dir + '/merge', donor_species + '.all')
            try:
                f1 = open(tempbamoutput + '.sorted.bam')
            except IOError:
                if sub_samples != []:
                    cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s  | %s call  --threads %s -m > %s.raw.vcf\n' % (
                        'bcftools', min(40, 40), database,
                        ' '.join(sub_samples), 'bcftools', min(40, 40), tempbamoutput)
                    cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
                        'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
                    cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
                        'bcftools', tempbamoutput, tempbamoutput)
                    cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
            f1 = open(os.path.join(input_script_merge, '%s.sh' % donor_species), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            f1.close()

f1 = open(os.path.join(input_script, 'allnewvcf.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'allnewmergevcf.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_merge, '*.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# clean up
import os
os.system('rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/MG/bwa')
os.system('rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/MG/merge/*.raw.vcf')
os.system('mv /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/MG/merge/*.all.flt.snp.vcf /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/MG/merge/rawdata/')

################################################### END ########################################################
################################################### SET PATH ########################################################
# filter vcf of metagenomes for clonal populations
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
from statistics import stdev

# setup path all clonal population
input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/MG/vcf'
input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/MG/vcf_merge'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/MG'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/MG/assembly/MG/tree'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round2/am_*')
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
ref_merge_vcf_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome/vcf'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/MG'
fastq_name = '_1.fasta'
lose_set = []
cluster_set = {}

################################################### Set up ########################################################

# set up cutoff
Sample_depth_cutoff = 2  # both forward and reverse reads cutoff in a sample
# SNPs merge into strain cutoff
MLF_cutoff = 0.1 # minimum MLF cutoff
MLF_variant_cutoff = 0.1 # maximum variance of MLF to be grouped as one genome seq
# strains merge into representative strains cutoff
SNP_cluster_cutoff = 10 # maximum SNP cutoff to be a cluster
top_relative = 50 # maximum nest to cluster
uniq_seq_count_cutoff = 4 # minimum presence of a seq and its downstream relatives to output
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

################################################### Function ########################################################
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def curate_REF(allels_set,lines_set_sub):
    SNP_count_genome_count = [0,0,0,0]
    for Subdepth_all in lines_set_sub:
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '')
        Subdepth_set = Subdepth.split(',')
        Subdepth_set_int = []
        for sub_depth in Subdepth_set:
            Subdepth_set_int.append(int(sub_depth))
        if sum(Subdepth_set_int) > 0:
            major_alt_frq = max(Subdepth_set_int)
            major_alt_frq_index = Subdepth_set_int.index(major_alt_frq)
            major_alt = allels_set[major_alt_frq_index]
            SNP_count_genome_count[major_alt_frq_index] += 1
    major_alt_frq = max(SNP_count_genome_count)
    major_alt_frq_index = SNP_count_genome_count.index(major_alt_frq)
    major_alt = allels_set[major_alt_frq_index]
    return [major_alt, major_alt_frq_index]


def create_SNP_seq(total_length):
    return 'N'*total_length

def change_SNP_seq(Sample_alignment_seq,subname,ALT,SNP_pos):
    Seq = Sample_alignment_seq[subname]
    Sample_alignment_seq[subname] = Seq[:SNP_pos] + ALT + Seq[SNP_pos + 1:]

def add_SNP_seq(genomename,ALT,MLF,SNP_pos,Main = True):
    if MLF > MLF_cutoff:
        Set = False
        if Main == True:
            # major alt stores as genomename.1
            change_SNP_seq(Sample_alignment_seq,genomename + '.1', ALT, SNP_pos)
            Sample_alignment_freq[genomename + '.1'].append(MLF)
        else:
            # minor alt stores as genomename.i
            for genome_sub in Sample_alignment[genomename]:
                MLF_genome = statistics.mean(Sample_alignment_freq[genome_sub])
                if MLF <= MLF_variant_cutoff +  MLF_genome and MLF >= MLF_genome - MLF_variant_cutoff:
                    # find the right sub genome sub
                    Sample_alignment_freq[genome_sub].append(MLF)
                    change_SNP_seq(Sample_alignment_seq, genome_sub, ALT, SNP_pos)
                    Set = True
                    break
            if Set == False:
                # the right sub genome sub not found
                # create another genome sub
                i = int(Sample_alignment[genomename][-1].split('.')[-1])
                newsubname = genomename + '.%s'%(i + 1)
                Sample_alignment[genomename].append(newsubname)
                Sample_alignment_seq.setdefault(newsubname, Sample_alignment_seq[genomename + '.1'])
                change_SNP_seq(Sample_alignment_seq,newsubname, ALT, SNP_pos)
                Sample_alignment_freq.setdefault(newsubname, [MLF])

def SNP_check_all(lines_set,SNP_pos):
    REF = lines_set[3]
    allels_set = [REF]
    if '.' not in lines_set[4]:
        allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        REF, REF_where = curate_REF(allels_set, lines_set[9:])
        for Subdepth_all in lines_set[9:]:
            genomename = Sample_name[genome_order]
            genome_order += 1
            Allels_frq = [0, 0, 0, 0]
            Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0]
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
            total_sub_depth = 0
            Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
            Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
            for num_allels in range(0, Total_alleles):
                allels = allels_set[num_allels]
                Subdepth_alleles = int(Subdepth[num_allels])
                if allels in Allels:
                    forward_alt = int(Subdepth_forward[num_allels])
                    reverse_alt = int(Subdepth_reverse[num_allels])
                    if forward_alt >= Sample_depth_cutoff and \
                            reverse_alt >= Sample_depth_cutoff:
                        # forward and reverse depth cutoff to call
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                        total_sub_depth += Subdepth_alleles
                        Allels_frq_sub[Allels[allels] * 2] += int(Subdepth_forward[num_allels])
                        Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_reverse[num_allels])
                else:
                    pass
            # find major alt and calculate frequency
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            # change to major ALT
            for genome_sub in Sample_alignment[genomename]:
                change_SNP_seq(Sample_alignment_seq, genome_sub, REF, SNP_pos)
            Sample_alignment_cov[genomename].append(total_sub_depth/MGsize[genomename])
            if total_sub_depth > 0:
                add_SNP_seq(genomename, Major_ALT[0], Major_ALT[1] / total_sub_depth, SNP_pos)
                if Minor_ALT!= []:
                    for sub_alt in Minor_ALT:
                        ALT, ALT_MLF = sub_alt
                        add_SNP_seq(genomename, ALT, ALT_MLF / total_sub_depth, SNP_pos, False)

def load_ref_vcf(ref_vcf_file):
    ref_chr = dict()
    ref_chr_order = dict()
    i = 0
    total_length = 0
    for files in ref_vcf_file:
        Set_length = False
        for lines in open(files,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS, Notused, REF, ALT = lines_set[0:5]
            CHR_POS = '%s__%s'%(CHR, POS)
            ref_chr.setdefault(CHR_POS,[])
            ref_chr[CHR_POS]=[REF,ALT]
            ref_chr_order.setdefault(CHR_POS,i)
            i += 1
        fasta_file = '.'.join(files.split('.')[:-1])+'.fasta'
        for record in SeqIO.parse(fasta_file, 'fasta'):
            record_name = str(record.id)
            if not Set_length:
                total_length += len(str(record.seq))
                Set_length = True
    return [ref_chr,ref_chr_order,total_length]

def outputtree(uniq_seq):
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in Sample_alignment_seq:
        seq_len = len(Sample_alignment_seq[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, Sample_alignment_seq[genomename]))
            SNP_alignment_output_parsi.append('%s    %s\n' % (genomename[-9:], Sample_alignment_seq[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.fasta', 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.parsi.fasta', 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in uniq_seq:
        seq_len = len(uniq_seq[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, uniq_seq[genomename]))
            SNP_alignment_output_parsi.append('%s    %s\n' % (genomename[-9:], uniq_seq[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max, seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.uniq.fasta', 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.uniq.parsi.fasta', 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()

def SNP_seq(seq1, seq2, total_length):
    SNP_total = 0
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
    return SNP_total

def find_close_relative(seq,uniq_seq,uniq_ID):
    SNP_set = dict()
    # look for closest relative from the closest timepoint
    allsample_order = sorted(uniq_seq,reverse=True)
    for samplenamesub in allsample_order:
        if samplenamesub != uniq_ID:
            meta_seq = uniq_seq[samplenamesub]
            SNP_total = SNP_seq(seq, meta_seq, total_length)
            SNP_set.setdefault(SNP_total,samplenamesub)
    SNP_set2 = sorted(SNP_set)
    print(SNP_set2,SNP_set)
    return [SNP_set2[0],SNP_set[SNP_set2[0]]]

def find_close_3_relative(close_relative_set,uniq_seq_count,pre_result,top_relative = 3):
    temp_result = close_relative_set.get(pre_result[1], [])
    closest_relative_ID = pre_result[1]
    if 'First' in closest_relative_ID:
        closest_relative_ID = find_first_uniqID(close_relative_set, closest_relative_ID)
    uniq_seq_count[closest_relative_ID] += 1
    print(temp_result,top_relative,pre_result)
    if pre_result[1] in cluster_set.get(donor_species,[]):
        # a pre set cluster
        return pre_result
    if temp_result != [] and top_relative > 0 and temp_result[0] <= SNP_cluster_cutoff:
        top_relative -= 1
        pre_result = find_close_3_relative(close_relative_set, uniq_seq_count, temp_result, top_relative)
    return pre_result

def find_first_uniqID(close_relative_set,closest_relative):
    for uniqID in close_relative_set:
        if close_relative_set[uniqID][1] == closest_relative:
            return uniqID

def outputfreq():
    timepoint = sorted(Sample_alignment_seq)
    cluster = True
    meta_seq_uniq = dict()
    uniq_seq = dict()
    # calculate uniq seq number
    for samplenamesub in timepoint:
        if Sample_alignment_freq[samplenamesub] != []:
            # not empty mapping
            meta_seq = Sample_alignment_seq[samplenamesub]
            SNP_short, closest_relative = [0, 'First']
            # uniq seq
            if meta_seq not in meta_seq_uniq:
                meta_seq_uniq.setdefault(meta_seq, [])
                uniq_seq.setdefault(samplenamesub, meta_seq)
    if len(uniq_seq) <= 10:
        cluster = False
    meta_seq_uniq = dict()
    uniq_seq = dict()
    output_freq = []
    close_relative_set = dict()
    close_relative_set_output = dict()
    uniq_seq_count = dict()
    # unique meta seq
    for samplenamesub in timepoint:
        if Sample_alignment_freq[samplenamesub] != []:
            # not empty mapping
            meta_seq = Sample_alignment_seq[samplenamesub]
            SNP_short, closest_relative = [0,'First']
            # uniq seq
            if meta_seq not in meta_seq_uniq:
                meta_seq_uniq.setdefault(meta_seq, [])
                uniq_seq.setdefault(samplenamesub,meta_seq)
                uniq_seq_count.setdefault(samplenamesub, 0)
            meta_seq_uniq[meta_seq].append(samplenamesub)
            uniq_ID = meta_seq_uniq[meta_seq][0]
            uniq_seq_count[uniq_ID] += 1
            if len(uniq_seq) == 1:
                # set up the first timepoint
                close_relative_set.setdefault(uniq_ID, [SNP_short,closest_relative])
                close_relative_set_output.setdefault(uniq_ID, [SNP_short, closest_relative])
                uniq_ID_secondstrain = uniq_ID.replace('.1', '.2')
                if uniq_ID_secondstrain in timepoint:
                    close_relative_set.setdefault(uniq_ID_secondstrain, [0, 'First2'])
                    close_relative_set_output.setdefault(uniq_ID_secondstrain, [0, 'First2'])
                uniq_ID_secondstrain = uniq_ID.replace('.1', '.3')
                if uniq_ID_secondstrain in timepoint:
                    close_relative_set.setdefault(uniq_ID_secondstrain, [0, 'First3'])
                    close_relative_set_output.setdefault(uniq_ID_secondstrain, [0, 'First3'])
            # find closest relative
            if uniq_ID not in close_relative_set:
                SNP_short, closest_relative = find_close_relative(meta_seq, uniq_seq,uniq_ID)
                close_relative_set.setdefault(uniq_ID, [SNP_short, closest_relative])
                if cluster:
                    SNP_short, closest_relative = find_close_3_relative(close_relative_set, uniq_seq_count,[SNP_short, closest_relative], top_relative)
                    closest_relative_ID = closest_relative
                    if 'First' in closest_relative:
                        closest_relative_ID = find_first_uniqID(close_relative_set,closest_relative)
                    SNP_short = SNP_seq(uniq_seq[closest_relative_ID], uniq_seq[uniq_ID], total_length)
                close_relative_set_output.setdefault(uniq_ID,[SNP_short,closest_relative])
            else:
                SNP_short, closest_relative = close_relative_set_output[uniq_ID]
            # output freq
            output_freq.append('%s\t%d\t%.3f\t%s\t%s\t\n'%(uniq_ID,
                                                   int(samplenamesub.replace('am','').split('.')[0]),
                                                   statistics.mean(Sample_alignment_freq[samplenamesub]),
                                                           closest_relative,SNP_short))
    outputtree(uniq_seq)
    vcf_file_filtered = open(vcf_file + '.freq', 'w')
    vcf_file_filtered.write('#seq_type\ttimepoint\tfreq\tclosest_relative\tSNP\n' + ''.join(output_freq))
    vcf_file_filtered.close()

################################################### Main ########################################################
# read MG size
MGsize = dict()
for lines in open(os.path.join(input_script + '/..','MGsize.txt')):
    lines_set = lines.split('\n')[0].split('\t')
    MGsize.setdefault(lines_set[0].split(fastq_name)[0],2*int(lines_set[-1]))

# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir + '/merge/rawdata','*.all.flt.snp.vcf'))
Coverage = []
Coverage.append('donor_species\ttimepoint\tavg_depth\tstdev_depth\tSNP_retrieved\ttotal_refSNP\t\n')
for vcf_file in all_vcf_file:
    try:
        vcf_file_filtered = open(vcf_file + '.freq', 'r')
    except FileNotFoundError:
        donor_species = os.path.split(vcf_file)[-1].split('.all.flt.snp.vcf')[0]
        if donor_species in lose_set:
            Sample_depth_cutoff = 1
        else:
            Sample_depth_cutoff = 2
        # always round 2 result
        ref_vcf_file = glob.glob(os.path.join(ref_merge_vcf_dir, '%s.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf' % (donor_species)))
        print(vcf_file,ref_vcf_file)
        Sample_name = []
        Sample_alignment = dict()
        Sample_alignment_seq = dict()
        Sample_alignment_freq = dict()
        Sample_alignment_cov = dict()
        Mapped_snp = 0
        ref_chr, ref_chr_order, total_length = load_ref_vcf(ref_vcf_file)
        for lines in open(os.path.join(input_script_merge_sub, '%s.sh' % (donor_species)), 'r'):
            if lines.startswith('bcftools mpileup '):
                # setup samples
                sample_set = lines.split('.fasta ')[1].split('\n')[0].split('  |')[0].split(' ')
                for samples in sample_set:
                    genomename = os.path.split(samples)[-1].split(fastq_name)[0]
                    Sample_name.append(genomename)
                    Sample_alignment_cov.setdefault(genomename,[])
                    Sample_alignment.setdefault(genomename,[genomename + '.1'])
                    Sample_alignment_seq.setdefault(genomename + '.1', create_SNP_seq(total_length))
                    Sample_alignment_freq.setdefault(genomename + '.1', [])
        print('running %s' % (donor_species))
        for lines in open(vcf_file, 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\n')[0].split('\t')
                CHR = lines_set[0]
                POS = int(lines_set[1])
                # a SNP confirmed in genome analysis
                CHR_POS = '%s__%s' % (CHR, POS)
                if CHR_POS in ref_chr:
                    REF_ref, ALT_ref = ref_chr[CHR_POS]
                    REF_meta, ALT_meta = lines_set[3:5]
                    # same ALTs
                    if [REF_ref, ALT_ref].sort() == [REF_meta, ALT_meta].sort():
                        # add to coverage
                        Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                        Mapped_snp += 1
                        # which ALT in SNP seq
                        SNP_pos = ref_chr_order[CHR_POS]
                        SNP_check_all(lines_set, SNP_pos)
        outputfreq()
        for genomename in Sample_alignment_cov:
            cov_genome = Sample_alignment_cov[genomename]
            cov_genome_notzero = [i for i in cov_genome if i > 0]
            if len(cov_genome) > 1:
                Coverage.append('%s\t%s\t%.1e\t%.1e\t%s\t%s\t\n'%(donor_species,
                                                              int(genomename.replace('am','').split('.')[0]),
                                                              statistics.mean(cov_genome),
                                                        stdev(cov_genome),len(cov_genome_notzero),
                                                                  total_length))
            elif len(cov_genome) == 1:
                Coverage.append('%s\t%s\t%.1e\t%.1e\t%s\t%s\t\n' % (donor_species,
                                                                    int(genomename.replace('am', '').split('.')[0]),
                                                                    cov_genome[0],
                                                                    0, len(cov_genome_notzero),
                                                                    total_length))

f1 = open(os.path.join(output_dir + '/merge', 'allcov.sum.txt'),'w')
f1.write(''.join(Coverage))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# after round 4 calculate dMRCA
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics

Round = 4
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round%s/merge_genome/vcf/'%(Round)
output_dir_merge2 = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round%s/merge_genome/vcf/'%(Round)
vcf_name = '.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf.singleevent.fasta'

def SNP_seq_diff(seq1, seq2):
    SNP_total = 0
    j = 0
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            j = i
    return SNP_total

# calculate dMRCA
all_vcf_file=glob.glob(os.path.join(output_dir_merge,'*%s'%(vcf_name)))+\
glob.glob(os.path.join(output_dir_merge2,'*%s'%(vcf_name)))
alloutput = []
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split('.all.flt.snp.vcf')[0]
    genome_num = 0
    SNP_to_ref_all = []
    record_seq = ''
    Ref_seq = ''
    for record in SeqIO.parse(vcf_file, 'fasta'):
        record_name = str(record.id)
        record_seq = str(record.seq)
        if 'reference' in record_name:
            Ref_seq = record_seq
        else:
            SNP_to_ref_all.append(SNP_seq_diff(Ref_seq, record_seq))
            genome_num += 1
    if genome_num > 0:
        alloutput.append('%s\t%s\t%.3f\t%.3f\n'%(donor_species,genome_num,statistics.mean(SNP_to_ref_all)/genome_num,
                                               statistics.stdev(SNP_to_ref_all)/genome_num))

f1 = open(os.path.join(input_script, 'alldmrca.txt'),'w')
f1.write('donor_species\tgenome_num\tdmrca\tdmrca_stdev\n' + ''.join(alloutput))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# round 4 merge all species together -> clonal population distance
import glob
import os

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/vcf_round4'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/'
genome_dir1 = '/scratch/users/anniz44/genomes/donor_species/species/'
genome_dir2 = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
allvcf_merge2sample = glob.glob(
    '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2/bwa/0/*filtered.*samplename.txt') + \
                      glob.glob(
                          '/scratch/users/anniz44/genomes/donor_species/species/vcf_merge_round2/bwa/0/*filtered.*samplename.txt')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4'
fastq_name = '_1.fastq'
fasta_name = '.all.spades.fasta'

os.system('rm -rf %s' % (input_script_sub))

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
    os.mkdir(output_dir + '/ref')
except IOError:
    pass

try:
    os.mkdir(input_script)
except IOError:
    pass

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# select representative fasta
species_set = dict()
sample_set = dict()
sample_output = []
for samplefiles in allvcf_merge2sample:
    for lines in open(samplefiles):
        sample = lines.split('\n')[0].split('\t')[0]
        donor_species = '_'.join(sample.split('_')[0:3])
        species = '_'.join(sample.split('_')[1:3])
        species_set.setdefault(species,set())
        samplenewname = os.path.split(samplefiles)[-1].replace('.all.flt.snp.vcf.filtered', '').replace(
            '.samplename.txt', '')
        sample_output.append('%s\t%s\t\n'%(sample,samplenewname))
        sample_set.setdefault(sample,samplenewname)
        break
    fastq = glob.glob(genome_dir1 + '*/fastq/%s%s' % (sample, fastq_name))+ \
            glob.glob(genome_dir2 + '*/fastq/%s%s' % (sample, fastq_name))+ \
            glob.glob(genome_dir1 + '*/fastq/%s%s' % (sample.replace('_sp_','_sp._'), fastq_name)) + \
            glob.glob(genome_dir2 + '*/fastq/%s%s' % (sample.replace('_sp_','_sp._'), fastq_name)) # sometimes sp. was shorten for sp
    if fastq != []:
        species_set[species].add(fastq[0])
    else:
        print('cannot find fastq for genome %s'%(sample))

output_file = open(os.path.join(input_script,'round4.sample.txt'),'w')
output_file.write('#genomename\tpopulation\t\n'+''.join(sample_output))
output_file.close()

# generate codes
ref_folder = output_dir + '/ref'
for species in species_set:
    fastq = species_set[species]
    if len(fastq) > 1:
        cmd = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
        species_folder_all = os.path.join(ref_folder, species + '_allspades')
        databasename = '%s/%s.all.spades.fasta' % (ref_folder, species)
        print(databasename)
        # select 0.5m reads
        for fastq_file in fastq:
            if '.all' + fastq_name not in fastq_file:
                donor_species_genome = fastq_file.split(fastq_name)[0]
                cmd += 'head -1000000 %s%s >> %s/%s.all%s\n' % (
                    donor_species_genome, fastq_name, ref_folder, species, fastq_name)
                cmd += 'tail -1000000 %s%s >> %s/%s.all%s\n' % (
                    donor_species_genome, fastq_name, ref_folder, species, fastq_name)
                cmd += 'head -1000000 %s%s >> %s/%s.all%s\n' % (
                    donor_species_genome, fastq_name.replace('1', '2'), ref_folder, species,
                    fastq_name.replace('1', '2'))
                cmd += 'tail -1000000 %s%s >> %s/%s.all%s\n' % (
                    donor_species_genome, fastq_name.replace('1', '2'), ref_folder, species,
                    fastq_name.replace('1', '2'))
        # spades
        cmd += 'spades.py --careful -1 %s/%s.all%s -2 %s/%s.all%s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' % \
               (ref_folder, species, fastq_name, ref_folder, species, fastq_name.replace('1', '2'),
                species_folder_all)
        cmd += 'mv %s/scaffolds.fasta %s\n' % (species_folder_all, databasename)
        # prodigal
        cmd += 'prodigal -a %s -d %s -q -i %s\n' % (
            databasename.replace(fasta_name, fasta_name.replace('.fasta', '.faa')),
            databasename.replace(fasta_name, fasta_name.replace('.fasta', '.fna')),
            databasename)
        # prokka not used
        cmd += '#prokka --kingdom Bacteria --outdir %s/prokka_%s  --locustag Bacter %s\n' % \
                      (ref_folder, species, databasename)
        cmd += '#mv %s/prokka_%s/*.gff %s/%s.all.gff\n' % (
            ref_folder, species, ref_folder, species)
        cmd += '#mv %s/prokka_%s/*.tsv %s/%s.all.tsv\n' % (
            ref_folder, species, ref_folder, species)
        cmd += '#mv %s/prokka_%s/*.faa %s/%s.all.faa\n' % (
            ref_folder, species, ref_folder, species)
        cmd += '#rm -rf %s/prokka_%s\n' % (ref_folder, species)
        # call SNPs
        cmd +=('rm -rf %s.fai\nbowtie2-build %s %s\n' % (databasename, databasename, databasename))
        sub_samples = []
        for files in fastq:
            if 'all' + fastq_name not in files:
                # not including the mixed sample
                donor_species_dir_file = os.path.split(files)[-1].split(fastq_name)[0]
                newfastqname = sample_set.get(donor_species_dir_file,
                                              donor_species_dir_file)
                print(donor_species_dir_file,newfastqname)
                tempbamoutput = os.path.join(output_dir + '/bwa/0', newfastqname)
                try:
                    f1 = open(tempbamoutput + '.sorted.bam')
                except IOError:
                    sub_samples.append(tempbamoutput + '.sorted.bam')
                    files2 = files.replace(fastq_name, fastq_name.replace('1', '2'))
                    cmd += 'bowtie2' + ' --threads %s --un-conc %s.unaligned -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                        min(40, 40), tempbamoutput, databasename, files, files2, 'samtools', min(40, 40),
                        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
                        tempbamoutput)
                    cmd += '#samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
                        tempbamoutput, tempbamoutput)
                    cmd += 'rm -rf %s.bam\n' % (tempbamoutput)
                    cmd += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
                    cmd += 'rm -rf %s.*.unaligned\n' % (tempbamoutput)
        # mileup
        tempbamoutput = os.path.join(output_dir + '/merge', species + '.all')
        try:
            f1 = open(tempbamoutput + '.sorted.bam')
        except IOError:
            if sub_samples != []:
                cmd += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s  | %s call  --threads %s -m > %s.raw.vcf\n' % (
                    'bcftools', min(40, 40), databasename,
                    ' '.join(sub_samples), 'bcftools', min(40, 40), tempbamoutput)
                cmd += '%s filter --threads %s -s  LowQual %s.raw.vcf > %s.flt.vcf \n' % (
                    'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
                cmd += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
                    'bcftools', tempbamoutput, tempbamoutput)
                cmd += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
        f1 = open(os.path.join(input_script_sub, '%s.sh' % species), 'w')
        f1.write('%s' % (''.join(cmd)))
        f1.close()

f1 = open(os.path.join(input_script, 'round4.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# round 4 SNP summary
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics

# setup path all clonal population
input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/vcf_round4'
input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/vcf_round4'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/tree_round3'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4'
output_notwanted = '/scratch/users/anniz44/genomes/donor_species/unwanted/outlier'
fastq_name = '_1.fastq'
genome_split = '_g'
deleting_file = []
Step = 2
Cov_dis = 50
rule1=[]
rule2=[]
rule3 = []
################################################### Function ########################################################
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def vcf_to_txt(lines,output_list,cluster_sub=[]):
    lines_set = lines.split('\n')[0].split('\t')
    if len(lines_set) >9:
        CHR = lines_set[0]
        POS = int(lines_set[1])
        temp_line = []
        temp_line.append(CHR)
        temp_line.append(str(POS))
        i = 9
        for Subdepth_all in lines_set[9:]:
            if (cluster_sub==[] and i not in deleting_set) or i in cluster_sub:
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                temp_line.append(str(total_sub_depth))
            i += 1
        output_list.append('\t'.join(temp_line)+'\n')
    else:
        print(lines)

def curate_REF(allels_set,lines_set_sub):
    SNP_count_genome_count = [0,0,0,0]
    for Subdepth_all in lines_set_sub:
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '')
        Subdepth_set = Subdepth.split(',')
        Subdepth_set_int = []
        for sub_depth in Subdepth_set:
            Subdepth_set_int.append(int(sub_depth))
        if sum(Subdepth_set_int) > 0:
            major_alt_frq = max(Subdepth_set_int)
            major_alt_frq_index = Subdepth_set_int.index(major_alt_frq)
            major_alt = allels_set[major_alt_frq_index]
            SNP_count_genome_count[major_alt_frq_index] += 1
    major_alt_frq = max(SNP_count_genome_count)
    major_alt_frq_index = SNP_count_genome_count.index(major_alt_frq)
    major_alt = allels_set[major_alt_frq_index]
    return [major_alt, major_alt_frq_index]

def outputvcf(output_name):
    vcf_file_filtered = open(vcf_file + '.%s.vcf' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list_vcf))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.snp.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.samplename.txt' % (output_name), 'w')
    vcf_file_filtered.write('\t'.join(Sample_name) + '\n')
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.POS.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_POS))
    vcf_file_filtered.close()

def outputcov(output_name,vcf_file_POS_candidate,cluster_sub=[]):
    if len(vcf_file_list) > 0:
        vcf_file_POS_candidate = '\n'.join(vcf_file_POS_candidate)
        vcf_file_POS_candidate_output = ('%s' % (vcf_file_POS_candidate))
        f1 = open(os.path.join(input_script,'grep.temp.txt'),'w')
        f1.write(vcf_file_POS_candidate_output)
        f1.close()
        os.system('grep -%s -f %s %s --no-group-separator > %s'% (
            Cov_dis, os.path.join(input_script,'grep.temp.txt'),
            vcf_file.split('.all.flt.snp.vcf')[0] + '.all.raw.vcf',
            vcf_file + '.%s.cov.temp' % (output_name)))
        os.system('cat %s | sort | uniq > %s' % (
            vcf_file + '.%s.cov.temp' % (output_name),
            vcf_file + '.%s.uniqcov.temp' % (output_name)))
        for lines in open(vcf_file + '.%s.uniqcov.temp' % (output_name), 'r'):
            if not lines.startswith("#"):
                vcf_to_txt(lines, cov_file_list,cluster_sub)
        os.system('rm -rf %s %s' % (vcf_file + '.%s.cov.temp' % (output_name),
                                    vcf_file + '.%s.uniqcov.temp' % (output_name)))
        vcf_file_filtered = open(vcf_file + '.%s.cov.txt' % (output_name), 'w')
        vcf_file_filtered.write(''.join(cov_file_list))
        vcf_file_filtered.close()

def outputtree(output_name):
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, SNP_alignment[genomename]))
            if '_cluster' in genomename:
                newgenomename = '%s_CL%s'%(genomename.split('_')[0],genomename.split('_cluster')[1].split('.')[0])
            else:
                newgenomename = '%s_CL0' % (genomename.split('_')[0])
            if '.cluster' in genomename:
                newgenomename = '%s.CL%s'%(newgenomename,genomename.split('.cluster')[1].split('.')[0])
            else:
                newgenomename = '%s.CL0' % (newgenomename)
            if 'reference' in genomename:
                newgenomename = genomename
            SNP_alignment_output_parsi.append('%s    %s\n' % (newgenomename, SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.%s.fasta' % (output_name), 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.parsi.fasta' %(output_name), 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()

def SNP_seq(seq1, seq2, POS_info,POS_info_CHR,POS_info_CHR_LEN,POS_info_output,G1,G2):
    SNP_total = 0
    j = 0
    POS_DIS = []
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            CHR = POS_info_CHR[i]
            POS = POS_info[i]
            LEN = POS_info_CHR_LEN[CHR]
            if CHR == POS_info_CHR[j]:  # same CHR
                DIS = abs(POS - POS_info[j])
                POS_DIS.append(DIS)  # POS diff
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, DIS, LEN))
            else:  # new CHR
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, 0, LEN))
            j = i
    return SNP_total

def SNP_distance_correct(distanceNEW, distanceREF, SNPREF):
    return int((distanceNEW/distanceREF)*SNPREF)

def find_neighbor(Cluster_SNP,neighbor,Cluster_SNP_set,cluster,Cluster_SNP_set_added):
    if neighbor != []:
        for record_name in neighbor:
            if record_name not in Cluster_SNP_set_added:
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    subneighbor = Cluster_SNP.get(record_name,[])
                    find_neighbor(Cluster_SNP, subneighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)

def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def dnORds(amino1, amino2):
    if amino1 == amino2:
        return 'S'
    else:
        return 'N'

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

def contig_to_gene(CHR, POS):
    all_genes = Mapping_loci.get(CHR,[])
    Reverse_chr = 0
    for a_gene in all_genes:
        POS1, POS2, GENE = a_gene
        if POS >= POS1 and POS <= POS2:
            Ref_seq_chr = Ref_seq.get(GENE, 'None')
            Gene_length = len(Ref_seq_chr)
            if GENE in Reverse:  # reversed
                POS_gene = Gene_length-(int(POS-POS1))
                Reverse_chr = 1
            else:
                POS_gene = int(POS-POS1)+1
            codon_start = POS_gene - 1 - int((POS_gene - 1) % 3)
            return [GENE,POS_gene,codon_start,Ref_seq_chr,Reverse_chr]
    return []

def SNP_check_all(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name,SNP_presence_cutoff,SNP_presence_sample_cutoff,cluster_sub=[]):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    Total_qualify = 0
    Total_qualify2 = 0
    Total_qualify_SNP = 0
    Total_qualify_SNP2 = 0
    Total_qualify_notSNP = 0
    Total_qualify_notSNP2 = 0
    Total_unqualify_alt_freq = 0
    Total_unqualify_alt_freq2 = 0
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    if cluster_sub!= []:
        lines_set_sub = [lines_set[i] for i in cluster_sub]
        Total_subsample = len(cluster_sub)
        if Total_subsample >= 15:
            SNP_presence_cutoff = 0.33  # for a large group of samples
        if Total_subsample <= 3:
            SNP_presence_cutoff = 1  # for a small group of samples
            SNP_presence_sample_cutoff = 2
    else:
        cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        if Step == 3:
            # re-calculate depth by a subset
            Depth4 = '%s,0,%s,0'%(sum([int(Subdepth_all.split(':')[-1].replace('\n', '').split(',')[0]) for Subdepth_all in lines_set_sub]),
                              sum([int(Subdepth_all.split(':')[-1].replace('\n', '').split(',')[1]) for Subdepth_all in
                                   lines_set_sub]))
        REF, REF_where = curate_REF(allels_set, lines_set_sub)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
                Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
                total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
                total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                        Allels_frq_sub[Allels[allels] * 2] += int(Subdepth_forward[num_allels])
                        Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_reverse[num_allels])
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    qualify_SNP = 0  # whether a qualified SNP
                    qualify_loci = 0
                    MLF = Major_ALT[1] / total_sub_depth
                    # rule 2
                    if donor_species in rule2:
                        if Major_ALT[1] >= 100 and MLF >= 0.7:
                            if Major_ALT[0] != REF:
                                # SNP
                                if total_sub_depth_forward - int(Subdepth_forward[REF_where]) >= 50 and \
                                total_sub_depth_reverse - int(Subdepth_reverse[REF_where]) >= 50:
                                    qualify_SNP = 1
                            else:
                                # REF
                                if int(Subdepth_forward[0]) >= 50 and \
                                int(Subdepth_reverse[0]) >= 50:
                                    qualify_loci = 1
                            temp_snp_line_pass += 'RULE'
                    # rule 3
                    if donor_species in rule3:
                        if Major_ALT[1] >= 30 and MLF >= 0.8:
                            if Major_ALT[0] != REF:
                                # SNP
                                if total_sub_depth_forward - int(Subdepth_forward[REF_where]) >= 15 and \
                                total_sub_depth_reverse - int(Subdepth_reverse[REF_where]) >= 15:
                                    qualify_SNP = 1
                            else:
                                # REF
                                if int(Subdepth_forward[0]) >= 15 and \
                                int(Subdepth_reverse[0]) >= 15:
                                    qualify_loci = 1
                            temp_snp_line_pass += 'RULE'
                    if total_sub_depth_forward >= Sample_depth_cutoff and \
                            total_sub_depth_reverse >= Sample_depth_cutoff:
                        # forward and reverse cutoff POS detected
                        if MLF >= Major_alt_freq_cutoff or qualify_loci == 1:
                        # major alt frequency cutoff
                            Total_qualify += 1
                            # check for qualified SNP
                            if Major_ALT[0] != REF:
                                SNP_seq[-1] = Major_ALT[0]  # unqualified SNP also include in alignment
                                # rule 1
                                if donor_species in rule1:
                                    if total_sub_depth_forward - int(Subdepth_forward[REF_where]) >= 3 and \
                                            total_sub_depth_reverse - int(Subdepth_reverse[REF_where]) >= 3 and \
                                            Major_ALT[1] >= 6 and MLF >= 0.95:
                                        qualify_SNP = 1
                                        temp_snp_line_pass += 'RULE'
                                # no specific rule
                                if total_sub_depth_forward - int(Subdepth_forward[REF_where]) >= SNP_depth_cutoff and \
                                        total_sub_depth_reverse - int(Subdepth_reverse[REF_where]) >= SNP_depth_cutoff and \
                                        Major_ALT[1] >= 2 * SNP_depth_cutoff:
                                    qualify_SNP = 1
                                    temp_snp_line_pass += 'PASS'
                            else:
                                Total_qualify_notSNP += 1
                            # tune cutoff
                            if Rough == 1:
                                if total_sub_depth_forward >= Sample_depth_cutoff2 and \
                                        total_sub_depth_reverse >= Sample_depth_cutoff2:
                                    # forward and reverse cutoff
                                    if MLF >= Major_alt_freq_cutoff2:
                                        # major alt frequency cutoff
                                        Total_qualify2 += 1
                                        if Major_ALT[0] != REF:
                                            if total_sub_depth_forward - int(Subdepth_forward[REF_where]) >= SNP_depth_cutoff2 and \
                                                    total_sub_depth_reverse - int(Subdepth_reverse[REF_where]) >= SNP_depth_cutoff2 and \
                                                    Major_ALT[1] >= 2 * SNP_depth_cutoff2:
                                                Total_qualify_SNP2 += 1
                                        else:
                                            Total_qualify_notSNP2 += 1
                                    else:
                                        Total_unqualify_alt_freq2 += 1
                        else:
                            # major alt frequency low
                            Total_unqualify_alt_freq += 1
                        # a qualified SNP
                        if qualify_SNP == 1:
                            Total_qualify_SNP += 1
                            SNP.add(genome_order)  # only take qualified SNP as valid SNP
                            SNP_seq[-1] = Major_ALT[0]
            sample_num += 1
        if Total_qualify / Total_subsample >= SNP_presence_cutoff and \
                Total_unqualify_alt_freq / Total_subsample <= Poor_MLF_freq_cutoff and\
                Total_qualify >= SNP_presence_sample_cutoff and \
                Total_qualify_SNP >= 1 and Total_qualify_SNP < Total_qualify and\
                Total_qualify_notSNP > 0:
            # -> qualified SNP
            # qualified samples cutoff + unqualified samples cutoff
            # at least 1 qualified SNP
            # calculate NS
            gene_info = contig_to_gene(CHR, POS)
            if gene_info!= []:
                Chr_gene, POS_gene,codon_start,Ref_seq_chr,Reverse_chr  = gene_info
                if Ref_seq_chr != 'None':
                    #  observed NS ratio calculated
                    temp_snp_line_NS= [Chr_gene,str(POS_gene),'']
                    if codon_start <= POS_gene - 1:
                        Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                        Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                        SNP_seq_chr = Ref_seq_chr
                        if len(Ref_seq_codon) == 3:
                            Ref_seq_aa = translate(Ref_seq_codon)[0]
                            temp_snp_line_AA += Ref_seq_aa
                            ALT_set = allels_set
                            ALT_set.remove(REF)
                            for ALT in ALT_set:
                                SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                SNP_seq_aa = translate(SNP_seq_codon)[0]
                                temp_snp_line_AA += SNP_seq_aa
                                temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                temp_snp_line_NS[-1]+=temp_NorS
            # output lines and output major alt
            if Rough == 1: #tune cutoff
                if Depth / Total_subsample >= median_coverage_cutoff2 and \
                        Total_qualify2 / Total_subsample >= SNP_presence_cutoff2 and \
                        Total_unqualify_alt_freq2 / Total_subsample <= Poor_MLF_freq_cutoff and \
                        Total_qualify2 >= SNP_presence_sample_cutoff2 \
                        and Total_qualify_SNP2 >= 1 and Total_qualify_SNP2 < Total_qualify2 \
                        and Total_qualify_notSNP2 > 0:
                    temp_snp_line_pass = 'PASS'
                else:
                    temp_snp_line_pass = 'NOT_PASS'
            else:
                if 'PASS' in temp_snp_line_pass:
                    temp_snp_line_pass = 'PASS'
                else:
                    temp_snp_line_pass = 'RULE'
            if CHR == CHR_old:
                # same CHR
                POS_DIS = abs(POS - POS_old)
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
            else:
                # diff CHR first SNP
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
            POS_old = POS
            CHR_old = CHR
            temp_snp_line.append(CHR)
            temp_snp_line.append(str(POS))
            temp_snp_line.append(REF)
            temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
            vcf_file_list.append('\t'.join(temp_snp_line)+ '\t' +'\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,'\t'.join(temp_snp_line_NS),temp_snp_line_AA))
            vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
            vcf_file_list_vcf.append('\t'.join(lines_set[0:9])+'\t'+'\t'.join(lines_set_sub)+'\n')
            if Rough != 1:
                i = 9
                j = 0
                SNP_alignment[reference_name] += REF
                for genomename in SNP_alignment:
                    if genomename != reference_name:
                        if i in cluster_sub:
                            SNP_alignment[genomename] += SNP_seq[j]
                            j += 1
                        i += 1
    return [CHR_old,POS_old]

################################################### Set up ########################################################
# set up output

try:
    os.mkdir(output_dir_merge + '/merge/tree')
except IOError:
    pass

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# set up steps
SNP_cluster = dict()
cluster_set = set()
if Step == 1:
    reference_set = ['']
    outputname_set = ['filteredrough']
elif Step == 2:
    reference_set = ['reference']
    outputname_set = ['filtered']
elif Step == 3:
    reference_set = []
    outputname_set = []
    # load cluster
    for lines in open(os.path.join(input_script, 'SNP_round2.sum'), 'r'):
        if not lines.startswith('donor_species'):
            lines_set = lines.split('\n')[0].split('\t')
            genomes = lines_set[1]
            cluster_type = lines_set[-3]
            cluster_total = int(lines_set[-2])
            if cluster_total > 1:
                SNP_cluster.setdefault(genomes, cluster_type)
                cluster_set.add(cluster_type)
    for cluster_type in cluster_set:
        reference_set.append('refer_%s'%(cluster_type))
        outputname_set.append('filtered.%s' % (cluster_type))
    cluster_set = list(cluster_set)

# set up cutoff super rough
if Step ==1 :
    Rough = 1  # tune cutoff
    median_coverage_cutoff = 4 # for group of samples -> used 3*2*0.66
    Major_alt_freq_cutoff = 0.7 # for SNP samples -> used
    SNP_depth_cutoff = 3 # for a SNP sample -> used
else:
    Rough = 0
    median_coverage_cutoff = 10  # avg coverage in all samples
    Major_alt_freq_cutoff = 0.9  # major alt freq in a sample
    SNP_depth_cutoff = 5  # both forward and reverse reads cutoff for SNP ALTs in a sample

# unchanged cutoff
Sample_depth_cutoff = 3  # both forward and reverse reads cutoff in a sample
SNP_presence_cutoff = 0.66  # percentage of samples passing the above criteria
SNP_presence_sample_cutoff = 3  # num of samples passing the above criteria
Poor_MLF_freq_cutoff = (1 - SNP_presence_cutoff2)*0.75 #the unqualifie samples should be mostly low cov but not two alleles (low major alt freq)

# set up strict cutoff
median_coverage_cutoff2 = 10 # avg coverage in all samples
Sample_depth_cutoff2 = 3 # both forward and reverse reads cutoff in a sample
Major_alt_freq_cutoff2 = 0.9 # major alt freq in a sample
SNP_presence_cutoff2 = 0.66 # percentage of samples passing the above criteria
SNP_presence_sample_cutoff2 = 3 # num of samples passing the above criteria
SNP_depth_cutoff2 = 5 # both forward and reverse reads cutoff for SNP ALTs in a sample
# cluster cutoff Step 3
SNP_total_cutoff_2 = 50
cluster_cutoff = 2
# coverage cutoff Step 1
total_coverage_cutoff = 0.7
genome_avg_coverage_cutoff = 10

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

################################################### Main ########################################################
# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir_merge + '/merge','*.all.flt.snp.vcf'))
for set_num in range(0,len(reference_set)):
    reference_name = reference_set[set_num]
    output_name = outputname_set[set_num]
    for vcf_file in all_vcf_file:
        SNP_presence_cutoff = SNP_presence_cutoff2 # for group of samples  -> used
        SNP_presence_sample_cutoff = SNP_presence_sample_cutoff2
        Poor_MLF_freq_cutoff = (1 - SNP_presence_cutoff) * 0.75
        print(vcf_file)
        Total = 0
        donor_species = os.path.split(vcf_file)[-1].split('.all.flt.snp.vcf')[0]
        SNP_tree_cmd = []
        SNP_tree_cmd2 = []
        vcf_file_list = []
        vcf_file_list_vcf = []
        Sample_name = []
        deleting_set = []
        vcf_file_POS = []
        vcf_file_POS_candidate = set()
        SNP_alignment = dict()
        SNP_alignment.setdefault(reference_name, '')
        cov_file_list = []
        Ref_seq = dict()
        Mapping = dict()
        Mapping_loci = dict()
        CHR_old = ''
        POS_old = 0
        SNP_cluster_donor_species = dict()
        for cluster_type in cluster_set:
            SNP_cluster_donor_species.setdefault(cluster_type,[])
        for lines in open(os.path.join(input_script_merge_sub, '%s.sh' % (donor_species)), 'r'):
            if lines.startswith('bcftools mpileup '):
                # setup samples
                sample_set = lines.split('.fasta ')[1].split('\n')[0].split('  |')[0].split(' ')
                samplenum = 9
                for samples in sample_set:
                    genomename = os.path.split(samples)[-1].split(fastq_name + '.sorted.bam')[0]
                    Sample_name.append(genomename.replace('.', ''))
                    if genomename in deleting_file:
                        deleting_set.append(samplenum)
                    else:
                        SNP_alignment.setdefault(genomename, '')
                        if SNP_cluster!= dict() and genomename in SNP_cluster:
                            SNP_cluster_donor_species[SNP_cluster[genomename]].append(samplenum)
                    samplenum += 1
        # subset samples
        Sample_subset = []
        if Step == 3:
            Sample_subset = SNP_cluster_donor_species.get(cluster_set[set_num], [])
            Sample_name = [Sample_name[i-9] for i in Sample_subset]
        if Step < 3 or (Sample_subset != [] and len(Sample_subset)>= cluster_cutoff):
            print('running %s' %(donor_species))
            # load database
            database_file = glob.glob(os.path.join(genome_root,
                                                   '%s/*all.spades.fna' % (donor_species)))
            if database_file != []:
                Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file[0])
            for lines in open(vcf_file, 'r'):
                if not lines.startswith("#"):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR = lines_set[0]
                    POS = int(lines_set[1])
                    Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                    if Total == 0:
                        Total = len(lines_set) - 9 - len(deleting_set)
                        if Total >= 15:
                            SNP_presence_cutoff = 0.33  # for a large group of genomes
                            Poor_MLF_freq_cutoff = 0.3
                        if Total <= 3:
                            SNP_presence_cutoff = 1  # for a small group of genomes
                            SNP_presence_sample_cutoff = 2
                    if Depth / Total >= median_coverage_cutoff:
                        # average depth in all samples cutoff
                        if "INDEL" not in lines_set[7] \
                                and (lines_set[6] != 'LowQual'):
                            if Step == 1:
                                CHR_old, POS_old =  SNP_check_all(lines_set, '',
                                                                  CHR_old,POS_old,reference_name,
                                                                  SNP_presence_cutoff,SNP_presence_sample_cutoff)
                            elif Step == 2:
                                CHR_old, POS_old = SNP_check_all(lines_set, '',
                                                                 CHR_old,POS_old,reference_name,
                                                                 SNP_presence_cutoff,SNP_presence_sample_cutoff)
                            else:
                                CHR_old, POS_old = SNP_check_all(lines_set, '',
                                                                 CHR_old, POS_old, reference_name,
                                                                 SNP_presence_cutoff, SNP_presence_sample_cutoff,
                                                                 Sample_subset)
            outputvcf(output_name)
            if Step != 1:
                # step 2 and 3 output fasta
                outputtree(output_name)
            if Step != 2:
                # output coverage
                try:
                    f1 = open(vcf_file + '.%s.cov.txt' % (output_name), 'r')
                except IOError:
                    outputcov(output_name, list(vcf_file_POS_candidate), Sample_subset)

# run parsi tree
if Step >1:
    parsi_files = glob.glob(os.path.join(output_dir_merge + '/merge', '*parsi.fasta'))
    for a_parsi_file in parsi_files:
        if 'RAxML_parsimonyTree' not in a_parsi_file:
            os.system('rm -rf %s %s' % (a_parsi_file + '.out.txt',
                                        a_parsi_file + '.out.tree'))
            SNP_tree_cmd3 = ('%s\n5\nV\n1\ny\n' % (a_parsi_file))
            f1 = open(os.path.join(input_script_sub, 'parsi.optionfile.txt'), 'w')
            f1.write(SNP_tree_cmd3)
            f1.close()
            os.system('rm -rf outfile outtree')
            os.system('dnapars < %s/parsi.optionfile.txt > %s/%s.parsi.output\n' % (
                input_script_sub, input_script_sub, os.path.split(a_parsi_file)[-1]))
            os.system('mv outfile %s' % (a_parsi_file + '.out.txt'))
            os.system('mv outtree %s' % (a_parsi_file + '.out.tree'))
    os.system('mv %s/merge/*.parsi* %s/merge/tree' % (
        output_dir_merge, output_dir_merge))

# run clustering and mafft tree
if Step == 2:
    # sum up SNP
    SNP_cutoff = []
    SNP_pair = []
    SNP_pair.append('G1\tG2\tSNP\t\n')
    POS_info_output = []
    POS_info_output.append('G1\tG2\tCHR\tPOS\tDIS\tLEN\t\n')
    Fasta_SNP = glob.glob(os.path.join(output_dir_merge + '/merge', '*.filtered.fasta'))
    tree_distance_output = []
    for fasta in Fasta_SNP:
        fasta_name = os.path.split(fasta)[-1]
        donor_species = fasta_name.split('.all.flt.snp.vcf.filtered.fasta')[0]
        POS_file = glob.glob(os.path.join(output_dir_merge + '/merge', donor_species + '.all.flt.snp.vcf.filtered.POS.txt'))[0]
        print(donor_species)
        tree_distance = dict()
        Seq_list = dict()
        Ref = ''
        max_cluster_diff = 0
        POS_info = []
        POS_info_CHR = []
        POS_info_CHR_LEN = dict()
        Cluster_SNP = dict()
        Cluster_SNP_set = dict()
        Cluster_SNP_set_added = set()
        # load SNP POS info
        for lines in open(POS_file, 'r'):
            CHR = lines.split('\t')[0]
            POS_info.append(int(lines.split('\t')[1]))
            POS_info_CHR.append(CHR)
            POS_info_CHR_LEN.setdefault(CHR, 0)
            POS_info_CHR_LEN[CHR] = max(POS_info_CHR_LEN[CHR], int(lines.split('\t')[1]))
        # load genome SNP fasta and calculate pair-wise SNPs
        for record in SeqIO.parse(fasta, 'fasta'):
            record_name = str(record.id)
            if 'reference' not in record_name:
                new_center = 1
                record_seq = str(record.seq)
                if Ref == '':
                    # set upt the first seq as ref
                    Ref = record_seq
                    REF_name = record_name
                    new_center = 0
                    SNP_total = 0
                    SNP_total_length = len(Ref)
                    #SNP_total_cutoff = max(SNP_total_cutoff_ratio * SNP_total_length, SNP_total_cutoff_2)
                    SNP_total_cutoff = SNP_total_cutoff_2
                else:
                    for record_before in Seq_list:
                        SNP_total = SNP_seq(Seq_list[record_before], record_seq, POS_info, POS_info_CHR,
                                            POS_info_CHR_LEN,
                                            POS_info_output, record_before, record_name)
                        SNP_pair.append('%s\t%s\t%s\t\n' % (record_before, record_name, SNP_total))
                        if SNP_total <= SNP_total_cutoff:
                            Cluster_SNP.setdefault(record_before, [])
                            Cluster_SNP[record_before].append(record_name)
                            max_cluster_diff = max(max_cluster_diff, SNP_total)
                    SNP_total = SNP_seq(Ref, record_seq, POS_info, POS_info_CHR, POS_info_CHR_LEN,
                                        POS_info_output, REF_name, record_name)
                    SNP_pair.append('%s\t%s\t%s\t\n' % (REF_name, record_name, SNP_total))
                    if SNP_total <= SNP_total_cutoff:
                        Cluster_SNP.setdefault(REF_name,[])
                        Cluster_SNP[REF_name].append(record_name)
                        max_cluster_diff = max(max_cluster_diff,SNP_total)
                tree_distance.setdefault(record_name, SNP_total)
                Seq_list.setdefault(record_name, record_seq)
        cluster = 0
        # cluster genomes by SNP distance
        for record_name in Cluster_SNP:
            neighbor = Cluster_SNP.get(record_name,[])
            if neighbor != [] and record_name not in Cluster_SNP_set_added:
                cluster += 1
                Cluster_SNP_set.setdefault(cluster,set())
                Cluster_SNP_set[cluster].add(record_name)
                Cluster_SNP_set_added.add(record_name)
                find_neighbor(Cluster_SNP, neighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)
        # output single genome
        for record_name in Seq_list:
            if record_name not in Cluster_SNP_set_added:
                cluster += 1
                Cluster_SNP_set.setdefault(cluster, set())
                Cluster_SNP_set[cluster].add(record_name)
                Cluster_SNP_set_added.add(record_name)
        Sub_cluster = len(Cluster_SNP_set)
        print(Cluster_SNP_set,fasta_name)
        for cluster in Cluster_SNP_set:
            tree_name_list = Cluster_SNP_set[cluster]
            tree_SNP_count = 'cluster%s' % (cluster)
            for tree_name in tree_name_list:
                tree_distance_output.append(
                    '%s\t%s\t%s\t%s\t%s\t\n' % (
                        donor_species, tree_name, tree_distance[tree_name], tree_SNP_count, Sub_cluster))
        SNP_cutoff.append('%s\t%s\t%s\t%s\t\n' % (donor_species, max_cluster_diff, SNP_total_cutoff, SNP_total_length))
    os.system('#mv %s %s' % (os.path.join(input_script, 'SNP_round2.sum'),
                            os.path.join(input_script, 'SNP_round2.old.sum')))
    f1 = open(os.path.join(input_script, 'SNP_round4.sum'), 'w')
    f1.write('donor_species\tGenome\tSNP_total\tcluster1\tsubcluster\t\n%s' % (''.join(tree_distance_output)))
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round4.sum.cutoff'), 'w')
    f1.write('donor_species\tmax_cluster_diff\tSNP_cutoff\tSNP_total_len\t\n%s' % (''.join(SNP_cutoff)))
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round4.POS.sum'), 'w')
    f1.write('%s' % (''.join(POS_info_output)))
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round4.allpair.sum'), 'w')
    f1.write('Genome1\tGenome2\tSNP_total\t\n%s' % (''.join(SNP_pair)))
    f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# extract species specific genes
import glob
import os

# setup path
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round*/*')+\
glob.glob('/scratch/users/anniz44/genomes/donor_species/jay/round*/*')
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome'
output_dir_merge2 = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round4/merge_genome'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome/selectedgenes'
summary_file = os.path.join(output_dir_merge,'all.gene.kegg.HS.anno.final.txt')
cluster_file = os.path.join(output_dir_merge + '/summary_genomeonly/','all.selected.gene.faa.High_select2.faa.uc')
#cluster_file = os.path.join(output_dir_merge + '/summary_genomeonly/','all.denovo.gene.faa.uc')
#summary_file = os.path.join(output_dir_merge,'all.gene.kegg.HS.anno.final.FP.txt')
#output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/merge_genome/notselectedgenes'
showdonor = False

def SNP_process(SNP_file,donor_species,genelist,all_SNP_output):
    for lines in open(SNP_file,'r'):
        if not lines.startswith('#'):
            lines = lines.split('\n')[0]
            lines_set = lines.split('\t')
            record_id = lines_set[0]
            record_id = '%s__C_%s_G_%s' % (donor_species,record_id.split('_')[1], record_id.split('_')[-1])
            if record_id in genelist:
                cluster, fun, donor_species_short, species, KO, KO1, KO2, KO3 = Genes[record_id]
                genus = species.split('_')[0]
                all_SNP_output.setdefault(genus,[])
                Genus_species.setdefault(genus,set())
                Genus_species[genus].add(species)
                all_SNP_output[genus].append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(donor_species_short,cluster,fun,
                                                                                 species,KO, KO1, KO2, KO3,
                                                                               record_id, '\t'.join(lines_set[1:])))
    return all_SNP_output

def major_minor(Major_ALT,ALT):
    if ALT == Major_ALT or ALT == '':
        return 'REF'
    else:
        return ALT

def SNP_process_to_matrix(output_list,genus,input_species = ''):
    Fun_set = []
    Fun_set_out = []
    Genome_type = dict()
    Genome_type_out = []
    All_out = []
    seq_len = 0
    # set up all genome types
    for lines in output_list:
        lines_set = lines.split('\t')
        donor_species = lines_set[0]
        if input_species == '' or input_species in donor_species:
            seq_len += 1
            genome_set = lines_set[15:]
            i = 0
            for genome_type in genome_set:
                i += 1
                genome_name = '%s__g%.4d' % (donor_species, i)
                Genome_type.setdefault(genome_name, [])
    for genome_name in Genome_type:
        Genome_type[genome_name] = [''] * seq_len
    # input all SNP types
    seq_len = 0
    for lines in output_list:
        lines_set = lines.split('\t')
        donor_species, cluster, fun, species, KO, KO1, KO2, KO3, geneID, POS, Major_ALT, Minor_ALT = lines_set[0:12]
        if showdonor:
            fun = '%s %s' % (fun, donor_species.split('_')[0])
        if input_species == '' or input_species in donor_species:
            seq_len += 1
            N_S = lines_set[14]
            if 'N' in N_S:
                N_S = 'N'
            else:
                N_S = 'S'
            genome_set = lines_set[15:]
            Fun_set.append('%s %s' % (fun, POS))
            Fun_set_out.append('%s %s\t%s\t%s\n' % (fun, POS, KO2, N_S))
            i = 0
            for genome_type in genome_set:
                i += 1
                genome_name = '%s__g%.4d' % (donor_species, i)
                Genome_type[genome_name][seq_len - 1] = major_minor(Major_ALT, genome_type)
    Genome_type_seq = dict()
    # unique genome types
    for genome_name in Genome_type:
        seq = '\t'.join(Genome_type[genome_name])
        Genome_type_seq.setdefault(seq, [])
        Genome_type_seq[seq].append(genome_name)
    # output Fun list
    temp_line = 'genome_name'
    for fun in Fun_set:
        temp_line += '\t%s' % (fun)
    temp_line += '\t\n'
    All_out.append(temp_line)
    # output one representative
    i = 0
    for seq in Genome_type_seq:
        i += 1
        rep_genome = Genome_type_seq[seq][0]
        genome_num = len(Genome_type_seq[seq])
        Genome_type_out.append('%s\t%s\n' % (rep_genome, genome_num))
        All_out.append('%s:%s\t%s\n' % (rep_genome.split('_')[0],i, seq))
    f1 = open(
        os.path.join(output_dir, 'details/all.vcf.frq.snp.%s.%s.matrix.txt' % (genus.replace(';', '_'), input_species
                                                                               )), 'w')
    f1.write(''.join(All_out))
    f1.close()
    f1 = open(os.path.join(output_dir, 'details/all.vcf.frq.snp.%s.%s.fun.txt' % (genus.replace(';', '_'), input_species
                                                                                  )), 'w')
    f1.write(''.join(Fun_set_out))
    f1.close()
    return 'finished'

def list_int_to_str(alist):
    i = 0
    for sub in alist:
        alist[i]=str(sub)
        i+=1
    return alist

def SNP_process_to_matrix_sum(output_list,genus,input_species = ''):
    Fun_set = dict()
    Fun_set_list = []
    Fun_set_out = []
    Genome_type = dict()
    Genome_type_S = dict()
    Genome_type_out = []
    All_out = []
    All_out_S = []
    seq_len = 0
    # set up all genome types
    for lines in output_list:
        lines_set = lines.split('\t')
        donor_species, cluster, fun, species, KO, KO1, KO2, KO3, geneID, POS, Major_ALT, Minor_ALT = lines_set[0:12]
        if input_species == '' or input_species in donor_species:
            if fun not in Fun_set_list:
                Fun_set_list.append(fun)
                Fun_set.setdefault(fun, [0, 0, KO1])
                seq_len += 1
            genome_set = lines_set[15:]
            i = 0
            for genome_type in genome_set:
                i += 1
                genome_name = '%s__g%.4d' % (donor_species, i)
                Genome_type.setdefault(genome_name, [])
                Genome_type_S.setdefault(genome_name, [])
    for genome_name in Genome_type:
        Genome_type[genome_name] = [0] * seq_len
        Genome_type_S[genome_name] = [0] * seq_len
    # input all SNP types
    any_S = 0
    for lines in output_list:
        lines_set = lines.split('\t')
        donor_species, cluster, fun, species, KO, KO1, KO2, KO3, geneID, POS, Major_ALT, Minor_ALT = lines_set[0:12]
        if input_species == '' or input_species in donor_species:
            seq_len = Fun_set_list.index(fun)
            N_S = lines_set[14]
            if 'N' in N_S:
                N_S = 'N'
                Fun_set[fun][0] += 1
            else:
                N_S = 'S'
                Fun_set[fun][1] += 1
                any_S += 1
            genome_set = lines_set[15:]
            i = 0
            for genome_type in genome_set:
                i += 1
                genome_name = '%s__g%.4d' % (donor_species, i)
                if genome_type != Major_ALT:
                    if N_S == 'N':
                        Genome_type[genome_name][seq_len - 1] += 1
                    else:
                        Genome_type_S[genome_name][seq_len - 1] += 1
    Genome_type_seq = dict()
    # unique genome types
    for genome_name in Genome_type:
        seq = '\t'.join(list_int_to_str(Genome_type[genome_name])) + '---' + '\t'.join(
            list_int_to_str(Genome_type_S[genome_name]))
        Genome_type_seq.setdefault(seq, [])
        Genome_type_seq[seq].append(genome_name)
    # output Fun list
    temp_line = 'genome_name'
    for fun in Fun_set:
        temp_line += '\t%s' % (fun)
        Fun_set_out.append('%s\t%s\t%s\t%s\n' % (fun, Fun_set[fun][-1], Fun_set[fun][0], Fun_set[fun][1]))
    temp_line += '\t\n'
    All_out.append(temp_line)
    All_out_S.append(temp_line)
    # output one representative
    for seq in Genome_type_seq:
        rep_genome = Genome_type_seq[seq][0]
        genome_num = len(Genome_type_seq[seq])
        Genome_type_out.append('%s\t%s\n' % (rep_genome, genome_num))
        All_out.append('%s\t%s\n' % (rep_genome, seq.split('---')[0]))
        All_out_S.append('%s\t%s\n' % (rep_genome, seq.split('---')[1]))
    f1 = open(
        os.path.join(output_dir, 'details/all.vcf.frq.snp.%s.%s.matrix.N.txt' % (genus.replace(';', '_'), input_species
                                                                                 )), 'w')
    f1.write(''.join(All_out))
    f1.close()
    if any_S > 0:
        f1 = open(
            os.path.join(output_dir, 'details/all.vcf.frq.snp.%s.%s.matrix.S.txt' % (genus.replace(';', '_'), input_species
                                                                                     )), 'w')
        f1.write(''.join(All_out_S))
        f1.close()
    f1 = open(os.path.join(output_dir, 'details/all.vcf.frq.snp.%s.%s.fun.txt' % (genus.replace(';', '_'), input_species
                                                                                  )), 'w')
    f1.write(''.join(Fun_set_out))
    f1.close()
    f1 = open(os.path.join(output_dir, 'details/all.vcf.frq.snp.%s.%s.genome.txt' % (genus.replace(';', '_'), input_species
                                                                                     )), 'w')
    f1.write(''.join(Genome_type_out))
    f1.close()
    return 'finished'

try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir + '/details')
except IOError:
    pass

Species_replace = dict()
Species_replace.setdefault('BA','Bifidobacterium_adolescentis')
Species_replace.setdefault('BL','Bifidobacterium_longum')
Species_replace.setdefault('PB','Parabacteroides_butyrate')
# load species specific genes
# both annotated and unannotated: not by selected, but filter by NS ratio
Clusters = dict()
for lines in open(summary_file):
    lines = lines.split('\n')[0]
    lines_set = lines.split('\t')
    gene,cluster,fun,species,donor_species,KO,KO1,KO2,KO3 = lines_set
    Clusters.setdefault(gene,[cluster,fun,donor_species,species,KO,KO1,KO2,KO3])

Anno = dict()
for lines in open(os.path.join(output_dir_merge,'all.cluster.anno.txt')):
    lines = lines.split('\n')[0]
    lines_set = lines.split('\t')
    cluster,fun,KO1,KO2,KO3 = lines_set
    Anno.setdefault(cluster,[fun,KO1,KO2,KO3])

Genes = dict()
Donor_species = dict()
for lines in open(cluster_file):
    lines = lines.split('\n')[0]
    lines_set = lines.split('\t')
    cluster = lines_set[1]
    gene = lines_set[8]
    if gene in Clusters:
        oldcluster, fun, donor_species, species, KO, KO1, KO2, KO3 = Clusters[gene]
        if donor_species in gene:
            if KO1 == 'None':
                if cluster in Anno and Anno[cluster][1] != '':
                    KO1, KO2, KO3 = Anno[cluster][1:]
            if cluster in Anno and Anno[cluster][0] != '':
                fun = Anno[cluster][0]
            if fun in ['NA', 'None', '']:
                if KO not in ['NA', 'None']:
                    fun = KO
                else:
                    fun = gene
            if KO == 'NA':
                KO = 'None'
            Genes.setdefault(gene, [cluster, fun, donor_species, species, KO, KO1, KO2, KO3])
            donor_species = gene.split('__')[0].split('_CL')[0]
            Donor_species.setdefault(donor_species,[])
            Donor_species[donor_species].append(gene)

# extract SNPs for species specific genes
all_SNP = glob.glob(os.path.join(output_dir_merge + '/vcf/genome_only/','*.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf.frq.snp'))+\
glob.glob(os.path.join(output_dir_merge2 + '/vcf/genome_only/','*.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.vcf.frq.snp'))

all_SNP_output = dict()
Genus_species = dict()
for SNP_file in all_SNP:
    donor_species = os.path.split(SNP_file)[-1].split('.all.flt.snp')[0]
    donor = donor_species.split('_')[0]
    donor_species_set = donor_species.split('_')
    try:
        donor_species = '%s_%s_%s' % (donor_species_set[0],
                                      donor_species_set[1][0:min(6, len(donor_species_set[1]))],
                                      donor_species_set[2][0:min(6, len(donor_species_set[2]))])
    except IndexError:
        donor_species = '%s_%s_%s' % (donor_species_set[0],
                                      donor_species_set[1],
                                      donor_species_set[1])
    if donor_species in Donor_species:
        genelist = list(Donor_species[donor_species])
        if 'cluster' in donor_species_set[-1]:
            try:
                donor_species += '_CL' + donor_species_set[-2].split('cluster')[1]
            except IndexError:
                donor_species += '_CL' + donor_species_set[-1].split('cluster')[1]
        all_SNP_output = SNP_process(SNP_file, donor_species, genelist, all_SNP_output)
    else:
        print('no genes for %s' % (donor_species))

for genus in all_SNP_output:
    f1 = open(os.path.join(output_dir,'all.vcf.frq.snp.%s.txt'%(genus.replace(';','_'))),'w')
    f1.write(
        '#donor_species\tcluster\tfun\tspecies\tKO\tKO1\tKO2\tKO3\tCHR\tPOS\tMajor_ALT\tMinor_ALT\tMajor_ALT_frq\tMinor_Alt_frq\tN_or_S\tgenotype_allgenomes\n')
    f1.write('\n'.join(all_SNP_output[genus]))
    f1.close()
    if ';' in genus:
        SNP_process_to_matrix(all_SNP_output[genus], genus, '')
    for species in Genus_species[genus]:
        if species in ['Bifido_adoles','Bifido_longum']:
            SNP_process_to_matrix(all_SNP_output[genus], genus, species)
        else:
            #SNP_process_to_matrix_sum(all_SNP_output[genus], genus, species)
            SNP_process_to_matrix(all_SNP_output[genus], genus, species)

################################################### END ########################################################
################################################### SET PATH ########################################################
# EMP -> gut microbiome abu
import glob
import os
import pandas as pd

metadata = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/EMP/EMP_HMP_meta.txt'
otu_table = '/scratch/users/anniz44/16S/public_16S/EMP/data/seqs/greengenes_cr/emp_cr_gg_13_8.qc_filtered.per_study_otu_tables/'
otu_taxon = '/scratch/users/anniz44/16S/public_16S/EMP/data/seqs/greengenes_cr/97_otu_taxonomy.txt'
species_list = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/EMP/species_order.new.txt'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/otu_table'

try:
    os.mkdir(output_dir)
except IOError:
    pass

def load_study(study_id,study_otu_table,Taxon):
    OTU_table = pd.read_csv(study_otu_table,sep='\t', header=0, index_col=0)
    df = pd.DataFrame(OTU_table)
    df_new = df.div(df.sum(axis=0, numeric_only=True), axis=1)
    alltaxonid = df_new.index.values.tolist()
    set_df = False
    for species in Taxon:
        Taxon_subset = list(Taxon[species])
        Taxon_subset = [i for i in Taxon_subset if i in alltaxonid]
        if Taxon_subset != []:
            df_new_subset = pd.DataFrame(df_new.loc[Taxon_subset].sum(axis=0, numeric_only=True))
            df_new_subset.columns = [species]
            if not set_df:
                df_species = df_new_subset
                set_df = True
            else:
                df_species = pd.merge(df_species, df_new_subset,
                                      how='outer', left_index=True, right_index=True)
        else:
            print(species)
    df_species.to_csv(os.path.join(output_dir, '%s.otu_table.txt' % (study_id)), sep='\t',
                      header=True)
    return 'finished'

def merge_study():
    allstudy_otu_table = glob.glob(os.path.join(output_dir, '*.otu_table.txt'))
    df_all = []
    for study_otu_table in allstudy_otu_table:
        OTU_table = pd.read_csv(study_otu_table, sep='\t', header=0, index_col=0)
        df = pd.DataFrame(OTU_table)
        df_all.append(df)
    df_all = pd.concat(df_all)
    df_all.to_csv(os.path.join(output_dir, 'all.otu_table.txt'), sep='\t',
                  header=True)
    return 'finished'


# import and filter healthy human gut studies
Study = set()
for lines in open(metadata,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    study_id = lines_set[1]
    sample_scientific_name = lines_set[4]
    DiseaseState = lines_set[-2]
    if sample_scientific_name == 'human gut metagenome' and DiseaseState == 'N':
        Study.add(study_id)

# import human gut species list
Taxon = dict()
for lines in open(species_list,'r'):
    if not lines.startswith('species_short'):
        lines_set = lines.split('\n')[0].split('\t')
        genus,species = lines_set[-1].split('_')
        if species not in ['sp','sp.']:
            taxon = '%s_%s'%(genus,species)
        else:
            taxon = '%s_' % (genus)
        Taxon.setdefault(taxon,set())

# map taxon id to species
for lines in open(otu_taxon,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    taxonid = lines_set[0]
    for taxon in Taxon:
        genus, species = taxon.split('_')
        if genus in lines and species in lines:
            Taxon[taxon].add(int(taxonid))

# use genus if species not found
Subtaxon = dict()
for taxon in Taxon:
    if Taxon[taxon] == set():
        genus, species = taxon.split('_')
        Subtaxon.setdefault(genus,set())

for lines in open(otu_taxon,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    taxonid = lines_set[0]
    for taxon in Subtaxon:
            if taxon in lines:
                Subtaxon[taxon].add(int(taxonid))

for taxon in Subtaxon:
    Taxon.setdefault(taxon,Subtaxon[taxon])

# summarize species abundance per sample
Study = list(Study)
for study_id in Study:
    study_id=study_id.replace(' ','')
    study_otu_table = glob.glob(os.path.join(otu_table,'emp_cr_gg_13_8.qc_filtered__study_id_%s__.biom.otu_table.txt'%(study_id)))
    if study_otu_table!= []:
        study_otu_table=study_otu_table[0]
        load_study(study_id,study_otu_table,Taxon)

merge_study()

################################################### END ########################################################
################################################### SET PATH ########################################################
# round 1 WGS am_Bacteroides_fra
import glob
import os

input_script_vcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_roundWGS'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
genome_dir = glob.glob('%s/round*/am_Bacteroides_fra*' % (genome_root))
database = '/scratch/users/anniz44/genomes/donor_species/selected_species/round1/am_Bacteroides_fragilis/am_Bacteroides_fragilis.all.spades1.fasta'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_roundWGS'
fastq_name = '_1.fastq'

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
    os.mkdir(input_script)
except IOError:
    pass

os.system('rm -rf %s'%(input_script_vcf))

try:
    os.mkdir(input_script_vcf)
except IOError:
    pass

def run_vcf_WGS(files,files2,database,tempbamoutput):
    # generate code
    cmds = 'bowtie2' + ' --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, files, files2, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    cmds += 'rm -rf %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    sample_set3 = ['%s.sorted.bam' % (tempbamoutput)]
    return [cmds, sample_set3]

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

def run_vcf_WGS_fast(files,files2,database,tempbamoutput):
    # generate code
    cmds = 'minimap2' + ' -ax sr -N 1000 -p 0.8 -t %s %s.mmi %s %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, files, files2, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    cmds += 'rm -rf %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    sample_set3 = ['%s.sorted.bam' % (tempbamoutput)]
    return [cmds, sample_set3]

# generate code bowtie2
sub_samples = []
tempbamoutputGenome = os.path.join(output_dir + '/merge', 'am_am_Bacteroides_fragilis' + '.all')
for folder in genome_dir:
    cmds2 = ''
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    # BN10 donors
    donor_species_fastqall = glob.glob(os.path.join(folder, 'fastq/*' + fastq_name)) + \
                             glob.glob(os.path.join(folder, '*' + fastq_name))
    if len(donor_species_fastqall) >= 3:
        for fastq_file in donor_species_fastqall:
            if '.all' + fastq_name not in fastq_file and '.all.spades' not in fastq_file:
                fastq_file_name = os.path.split(fastq_file)[-1]
                filename = fastq_file.split(fastq_name)[0]
                fastq_file2 = filename + fastq_name.replace('1', '2')
                # WGS for confirmation
                results = run_vcf_WGS(fastq_file, fastq_file2,
                                      database,
                                          os.path.join(output_dir + '/bwa/0',
                                                       fastq_file_name))
                cmds2 += results[0]
                sub_samples += results[1]
    f1 = open(os.path.join(input_script_vcf, '%s.vcf.sh' % (donor_species)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds2)))
    f1.close()

f1 = open(os.path.join(input_script, 'vcf_roundWGS.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
cmds = merge_sample(database, tempbamoutputGenome, sub_samples)
f1.write(cmds)
f1.close()

# generate code minimap2
input_script_vcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcf_roundWGS_fast'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_roundWGS_fast'


os.system('rm -rf %s'%(input_script_vcf))

try:
    os.mkdir(input_script_vcf)
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

sub_samples = []
tempbamoutputGenome = os.path.join(output_dir + '/merge', 'am_Bacteroides_fragilis' + '.all')
for folder in genome_dir:
    cmds2 = ''
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    # BN10 donors
    donor_species_fastqall = glob.glob(os.path.join(folder, 'fastq/*' + fastq_name)) + \
                             glob.glob(os.path.join(folder, '*' + fastq_name))
    if len(donor_species_fastqall) >= 3:
        for fastq_file in donor_species_fastqall:
            if '.all' + fastq_name not in fastq_file and '.all.spades' not in fastq_file:
                fastq_file_name = os.path.split(fastq_file)[-1]
                filename = fastq_file.split(fastq_name)[0]
                fastq_file2 = filename + fastq_name.replace('1', '2')
                # WGS for confirmation
                results = run_vcf_WGS_fast(fastq_file, fastq_file2,
                                      database,
                                          os.path.join(output_dir + '/bwa/0',
                                                       fastq_file_name))
                cmds2 += results[0]
                sub_samples += results[1]
    f1 = open(os.path.join(input_script_vcf, '%s.vcf.sh' % (donor_species)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds2)))
    f1.close()

f1 = open(os.path.join(input_script, 'vcf_roundWGS_fast.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
cmds = merge_sample(database, tempbamoutputGenome, sub_samples)
f1.write(cmds)
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
#changename
import glob
import os

input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round*/'
pangenome_dir = '/scratch/users/anniz44/genomes/pan-genome/roary'
fasta_name = '.fasta.corrected.fasta'

allgenome = glob.glob(os.path.join(genome_root,'*'))
Species_set = dict()
Species_set2 = dict()
for genome_folder in allgenome:
    if 'notused' not in genome_folder:
        genome_folder_name = os.path.split(genome_folder)[-1]
        genome_folder_name_set = genome_folder_name.split('_')
        species = '%s_%s'%(genome_folder_name_set[1],
                                              genome_folder_name_set[2])
        newgenome_folder_name = '%s%s%s' % (genome_folder_name_set[1][0:2],
                                              genome_folder_name_set[2][0:1].upper(),genome_folder_name_set[2][1:2])
        newgenome_folder_name = newgenome_folder_name.replace('BiAd','BA').replace('BiLo','BL')
        if newgenome_folder_name not in Species_set:
            Species_set.setdefault(newgenome_folder_name,[species])
            Species_set2.setdefault(species, newgenome_folder_name)
        elif species not in Species_set[newgenome_folder_name]:
            print('duplicated', species,Species_set[newgenome_folder_name])

for genome_folder in allgenome:
    if 'notused' not in genome_folder:
        genome_folder_root, genome_folder_name = os.path.split(genome_folder)
        genome_folder_name_set = genome_folder_name.split('_')
        species = '%s_%s' % (genome_folder_name_set[1],
                             genome_folder_name_set[2])
        newgenome_folder_name = Species_set2[species]
        os.system('mv %s %s' % (genome_folder, os.path.join(genome_folder_root,genome_folder_name.replace(species,newgenome_folder_name))))

allgenome = glob.glob(os.path.join(genome_root, '*'))
for genome_folder in allgenome:
    if 'notused' not in genome_folder:
        allfiles = glob.glob(os.path.join(genome_folder, '*')) + glob.glob(os.path.join(genome_folder + '/fastq', '*'))
        for afile in allfiles:
            file_root, file_name = os.path.split(afile)
            if file_name != 'fastq' and file_name != 'spades' and file_name != 'others':
                file_name_set = file_name.split('_')
                try:
                    species = '%s_%s' % (file_name_set[1],
                                         file_name_set[2].split('.')[0])
                    newgenome_folder_name = Species_set2[species]
                except KeyError:
                    try:
                        species = '%s_%s' % (file_name_set[0],
                                             file_name_set[1].split('.')[0])
                        newgenome_folder_name = Species_set2[species]
                    except KeyError:
                        try:
                            if 'aa_Parabacteroides' in file_name:
                                species = 'Parabacteroides_goldsteinii'
                                newgenome_folder_name = Species_set2[species]
                        except KeyError:
                            pass
                os.system('mv %s %s' % (
                    afile, os.path.join(file_root, file_name.replace(species, newgenome_folder_name))))

Output = []
for species in Species_set2:
    Output.append('%s\t%s\t\n'%(species,Species_set2[species]))

f1 = open(input_script + '/changespeciesname.txt','w')
f1.write(''.join(Output))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
#pangenome
import glob
import os

input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
input_script_pan = input_script + '/prokka'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round*/'
genome_root2 = '/scratch/users/anniz44/genomes/donor_species/jay/round*/'
pangenome_dir = '/scratch/users/anniz44/genomes/pan-genome/roary/species'
fasta_name = '.fasta.corrected.fasta'

try:
    os.mkdir(pangenome_dir)
except IOError:
    pass

os.system('rm -rf %s'%(input_script_pan))
try:
    os.mkdir(input_script_pan)
except IOError:
    pass

def prokka(output_dir, species, genome):
    genome_filename = os.path.split(genome)[-1].replace(fasta_name,'.gff')
    cmd = 'prokka --kingdom Bacteria --outdir %s/prokka_%s  --locustag Bacter %s\n' % \
           (output_dir, species, genome)
    cmd += 'mv %s/prokka_%s/*.gff %s/%s\n' % (
        output_dir, species, output_dir, genome_filename)
    cmd += 'rm -rf %s/prokka_%s\n' % (output_dir, species)
    return cmd

allgenome = glob.glob(os.path.join(genome_root,'*')) + glob.glob(os.path.join(genome_root2,'*'))

for genome_folder in allgenome:
    allfasta = glob.glob(os.path.join(genome_folder,'*%s'%(fasta_name)))
    if allfasta != []:
        species = os.path.split(genome_folder)[-1].split('_')[1]
        donor_species = os.path.split(genome_folder)[-1]
        species_output = os.path.join(pangenome_dir,species)
        cmd = ''
        try:
            os.mkdir(species_output)
        except IOError:
            pass
        for fasta in allfasta:
            fasta_filename = os.path.split(fasta)[-1]
            gff_file = glob.glob(os.path.join(species_output,fasta_filename.replace(fasta_name,'.gff')))
            if gff_file!= []:
                print('has %s'%(gff_file[0]))
            else:
                cmd += prokka(species_output, donor_species, fasta)
        if cmd!= '':
            f1 = open(os.path.join(input_script_pan, '%s.sh' % (donor_species)), 'a')
            f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s' % (''.join(cmd)))
            f1.close()

f1 = open(os.path.join(input_script, 'allprokka.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_pan, '*.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# before running pangenome.py for IBD dataset, subset 5 donors per cluster
import glob
import os
pangenome_dir = '/scratch/users/anniz44/genomes/pan-genome/roary/species'

allpan = glob.glob('%s/*IBD*'%(pangenome_dir))
for pan in allpan[1:]:
    print(pan)
    subpan = glob.glob('%s/*'%(pan))
    species = os.path.split(pan)[-1]
    print(subpan)
    for a_folder in subpan:
        allgenomes = glob.glob('%s/*.gff'%(a_folder))
        donor_set = set()
        for genome in allgenomes:
            donor_set.add(os.path.split(genome)[-1].split('_')[0])
        donor_set = list(donor_set)
        print(donor_set)
        total_folder = int(len(donor_set)/5)
        print(total_folder)
        for i in range(0,total_folder):
            newsubfolder = os.path.join(a_folder,'%s_%s'%(species,i))
            print(newsubfolder)
            try:
                os.mkdir(newsubfolder)
            except IOError:
                pass
            cmd = ''
            for donor in donor_set[i*5:(i+1)*5]:
                cmd +=('mv %s/%s_*.gff %s/\n'%(a_folder, donor,newsubfolder))
            os.system(cmd)
            if i == total_folder -1:
                cmd = ''
                for donor in donor_set[i:]:
                    cmd += ('mv %s/%s_*.gff %s/\n' % (a_folder, donor, newsubfolder))
                os.system(cmd)

################################################### END ########################################################
################################################### SET PATH ########################################################

# pipeline #newcluster.sh
# step 1 clonal population clustering
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round1/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/'
ref_dir = '/scratch/users/anniz44/genomes/donor_species/WGS_old/WGS/vcf_round1/co-assembly'
print('python mapping_WGS_cluster.py -i %s -s %s -o %s -fa %s -ref %s '%(genome_root, input_script,
                                                                               output_dir,'.fasta.corrected.fasta',ref_dir))
print('sh %s/allWGS_cluster.sh'%(input_script))
print('python vcf_process.py  -addqual False -i %s -s %s -o %s -rd 1'%(genome_root, input_script, output_dir))
# python SNPfilter_WGS_cluster_all.py
print('python SNPfilter_WGS_cluster.py -i %s/vcf_round1/merge/ -vcf .filtered.vcf -s %'%(output_dir,input_script))
print('python cluster_CP_mapping.py -i %s -s %s -clustering 1'%(output_dir,input_script))
# check files
print('python cluster_CP_mapping.py -i %s -s %s -clustering 2'%(output_dir,input_script))
#print('python cluster_CP.py -i %s -s %s -o %s'%(pangenome_dir,input_script,output_dir))
#print('python cluster_CP.py -i %s -s %s -o %s -clustering 2'%(pangenome_dir,input_script,output_dir))

# step 2 co-assembly and mapping
# Tag inside mapping_WGS.py sh mapping_WGS.sh
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round*/'
pangenome_dir = '/scratch/users/anniz44/genomes/donor_species/vcf_round1/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/'
ref_genome = '/scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/reference.genome.txt'

print('python mapping_WGS.py -i %s -s %s -o %s -fa %s -ref %s -cl %s'%(genome_root, input_script,
                                                                       output_dir,'.fasta.corrected.fasta','None', pangenome_dir))
print('sh %s/allWGS.sh'%(input_script))
# use clusters! -> vcfprocess.py
print('python vcf_process.py -i %s -s %s -o %s'%(genome_root, input_script, output_dir + '/WGS/'))

# step 3 remove recombination and bad-quality SNPs
vcf_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge'
vcf_format = '.filtered.vcf'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
# use clusters! SNPfilter_WGS_all.py
print('python SNPfilter_WGS.py -i %s -vcf %s -s %s'%(vcf_folder, vcf_format,input_script))
# check error
#'ls -lh *.err | grep -v 'K''
#'ls -lh *.all.raw.vcf | grep -v 'G' | grep -v 'M'  '
#'ls -lh *.all.flt.snp.vcf | grep -v 'K' | grep -v 'M' '
# -> move finished clonal population to finished -> run mapping_WGS again for unfinished ones
# annotate_allgenome.sh extract genome fasta and cluster ratio
# step 4 calculate dnds
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
cutoff_file = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/total_SNP_cutoff.txt'
#contig = 10000
contig = 5000
print('python dnds.py -s %s -o %s -cutoff %s -contig %s'%(input_script, output_dir_merge,cutoff_file,contig))

# step 5 parallel evolution
co_assembly_dir = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
cutoff_file = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/total_SNP_cutoff_lineage.txt'

print('python parallel_evolution.py -i %s -s %s -o %s -cutoff %s'%(co_assembly_dir,input_script, output_dir_merge,cutoff_file))
print('please run %s/%s'%(input_script,'allannotate.sh'))
print('please run %s/%s'%(input_script,'allannotate_all.sh'))

# step 6 dmrca
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
print('python dnds.py -s %s -o %s'%(input_script, output_dir_merge))

# sum annotation + runprokka.py + sumprokka.py + snpsum.py
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
print('python annotate_sum.py -o %s'%(output_dir_merge))

# step 7 core or flexible
#allgenome = '/scratch/users/anniz44/genomes/donor_species/*/round*/'
#input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
#output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
#print('python core_flexible.py -i \"%s\" -s %s -o %s'%(allgenome,input_script, output_dir_merge))
#print('please run %s/%s'%(input_script,'pangenome.sh'))
#output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
#print('python core_flexible_sum.py -o %s'%(output_dir_merge)) > not used

# sum core or flexible -> core_flexible2.py or gene core flexible in R cross-genus
# input R file all.genome.gene.faa.uc.species.sum from annotate_sum.py
# all.genome.gene.faa.uc usearch 0.9 cluster from /scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/*/*.noHM.fasta.faa
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
cutoff_file = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/total_SNP_cutoff.txt'
core_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.species.High_select2.all.gene.speciesnum.core_flexible.txt'
print('python dnds.py -s %s -o %s -cutoff %s -contig 9440 -core %s'%(input_script, output_dir_merge,cutoff_file,core_file))

# stepn co-assembly mapping to metagenomes
assembly_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/withPE/'
fastq1 = '/scratch/users/anniz44/Metagenomes/public_metagenomes/fecal_human'
output_dir = '/scratch/users/anniz44/genomes/donor_species/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/human_MG_ESSG/'

## ESSG
print('python mapping_meta.31marker.py -i %s -m %s -mfq .fasta -o %s -s %s'%(assembly_folder, fastq1, output_dir, input_script))
print('sh %s/allMGvcf.sh'%(input_script))

print('python SNPfilter_meta.31marker.py -i %s -mfq _1.fasta -o %s'%(assembly_folder,output_dir))

assembly_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/withPE/'
fastq1 = '/scratch/users/anniz44/Metagenomes/public_metagenomes/fecal_human'
output_dir = '/scratch/users/anniz44/genomes/donor_species/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/human_MG'
snp_dir = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge'
## all genomes map to MG
print('python mapping_meta.py -i %s -m %s -mfq .fasta -o %s -s %s -snp %s'%(assembly_folder, fastq1, output_dir, input_script,snp_dir))
print('sh %s/allMGvcf.sh'%(input_script))

print('python SNPsum_meta.py -o %s'%(output_dir))
################################################### END ########################################################
################################################### SET PATH ########################################################
# merge annotation and snpfiles -> mergeanno.py
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

allsnpfiles = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/*.all.parsi.fasta.sum.txt')
annotation_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.denovo.gene.faa.cluster.aa.all.eggnog.sum.species.sum'
changename = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.denovo.gene.faa.changename.txt'

genename = dict()
for lines in open(changename,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    donorspecies,nouse,oldgene,newgene = lines_set[0:4]
    genename.setdefault('%s\t%s'%(donorspecies,oldgene),newgene)

annotation_set = dict()
for lines in open(annotation_file,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    newgene,eggnog,anno,cog,cog1,cog2 = lines_set[4:10]
    annotation_set.setdefault(newgene,'%s\t%s\t%s\t%s'%(eggnog,anno[:40],cog1,cog2))

for snpfile in allsnpfiles:
    alloutput = []
    donorspecies = os.path.split(snpfile)[-1].split('.all.parsi.fasta.sum.txt')[0]
    donorspecies = donorspecies.replace('_PB_','_PaDi_')
    POS_all_0 = 0
    oldCHR = ''
    for lines in open(snpfile, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        oldgene = lines_set[5]
        if oldgene == 'Gene':
            lines = 'Newgene\tEggnog\tAnno\tCOG1\tCOG2\tPOSnew\t' + lines
        else:
            CHR = lines_set[0]
            POS = int(lines_set[1])
            if oldCHR == '':
                oldCHR = CHR
            if oldCHR!=CHR:
                try:
                    total_length = oldCHR.split('size')[1]
                except IndexError:
                    try:
                        total_length = oldCHR.split('length_')[1].split('_cov')[0]
                    except IndexError:
                        total_length = 0
                oldCHR = CHR
                POS_all_0 += int(total_length)
            newgene = genename.get('%s\t%s'%(donorspecies,oldgene),'None')
            annoall = annotation_set.get(newgene,'None\tNone\tNone\tNone')
            lines = '%s\t%s\t%s\t'%(newgene,annoall,POS+POS_all_0) + lines
        alloutput.append(lines)
    f1 = open(snpfile.replace('.all.parsi.fasta.sum.txt','.snpsum.txt'), 'w')
    f1.write(''.join(alloutput))
    f1.close()

os.system('mv /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/*snpsum.txt /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/snpsummary/')

################################################### END ########################################################
################################################### SET PATH ########################################################
# snpsum.py
# merge annotation and snpfiles IBD
################################################### END ########################################################
################################################### SET PATH ########################################################
# re cluster and remove rec
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

vcf_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/allvcfdetails/'
vcf_folderold = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/allvcfdetailsold/'
vcf_format = '.filtered.vcf'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
donor_species = '1_PB_IBD_0_clustercluster1'

#os.system('mv %s/%s*.donor.* %s/'%(vcf_folder, donor_species,vcf_folderold))

print('python SNPfilter_WGS.py -i %s -vcf %s -s %s -cluster %s'%(vcf_folder, vcf_format,input_script,donor_species))

input_script_sub = input_script + '/tree'
def run_parsi(a_parsi_file):
    if 'RAxML_parsimonyTree' not in a_parsi_file:
        donor_species = os.path.split(a_parsi_file)[-1]
        input_script_sub_temp = input_script_sub + '/' + donor_species
        try:
            os.mkdir(input_script_sub_temp)
        except IOError:
            pass
        os.system('rm -rf %s %s' % (a_parsi_file + '.out.txt',
                                    a_parsi_file + '.out.tree'))
        SNP_tree_cmd3 = ('%s\n5\nV\n1\ny\n' % (a_parsi_file))
        f1 = open(os.path.join(input_script_sub_temp, '%s.parsi.optionfile.txt'%(donor_species)), 'w')
        f1.write(SNP_tree_cmd3)
        f1.close()
        os.system('rm -rf outfile outtree')
        os.system('dnapars < %s/%s.parsi.optionfile.txt > %s/%s.parsi.output\n' % (
            input_script_sub_temp, donor_species,input_script_sub_temp, donor_species))
        os.system('mv outfile %s' % (a_parsi_file + '.out.txt'))
        os.system('mv outtree %s' % (a_parsi_file + '.out.tree'))

def outputtree_parsi(fastafile,output_file):
    SNP_alignment = dict()
    for record in SeqIO.parse(fastafile, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        SNP_alignment.setdefault(record_id,record_seq)
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        newgenomename = genomename
        if len(genomename) > 8:
            newgenomename = genomename[0:4] + '_' + genomename[-4:]
        if seq_len > 0:
            SNP_alignment_output_parsi.append('S%s    %s\n' % (newgenomename, SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(output_file, 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()

fastafile = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/allvcfdetails/1_PB_IBD_0_clustercluster1.all_subcluster0.donor.D14.raw.vcf.filtered.vcf.final.fasta'
output_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/allvcfdetails/tree/1_PB_IBD_0_clustercluster1.all_subcluster0.donor.D14.raw.vcf.filtered.vcf.final.parsi.fasta'
outputtree_parsi(fastafile,output_file)
run_parsi(output_file)
################################################### END ########################################################
################################################### SET PATH ########################################################
# check all trees -> checkalltree.py
import glob,os
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
input_script_sub = input_script + '/tree'
os.system('rm -f %s/*'%(input_script_sub))
parsi_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
tree_folder = '%s/tree/'%(parsi_folder)

def run_parsi(a_parsi_file):
    if 'RAxML_parsimonyTree' not in a_parsi_file:
        donor_species = os.path.split(a_parsi_file)[-1]
        input_script_sub_temp = input_script_sub + '/' + donor_species
        os.system('rm -r %s' % input_script_sub_temp)
        try:
            os.mkdir(input_script_sub_temp)
        except IOError:
            pass
        SNP_tree_cmd3 = ('%s\n5\nV\n1\ny\n' % (a_parsi_file))
        f1 = open(os.path.join(input_script_sub_temp, '%s.parsi.optionfile.txt'%(donor_species)), 'w')
        f1.write(SNP_tree_cmd3)
        f1.close()
        os.system('cd %s'%(input_script_sub_temp))
        os.system('rm -r outfile outtree')
        os.system('dnapars < %s/%s.parsi.optionfile.txt > %s/%s.parsi.output\n' % (
            input_script_sub_temp, donor_species,input_script_sub_temp, donor_species))
        os.system('mv outfile %s' % (a_parsi_file + '.out.txt'))
        os.system('mv outtree %s' % (a_parsi_file + '.out.tree'))

for fastafile in glob.glob('%s/*.all.parsi.fasta'%(parsi_folder)):
    output_file = '%s/%s.out.tree'%(tree_folder,os.path.split(fastafile)[-1])
    #print(fastafile)
    filesize = 0
    try:
        filesize = int(os.path.getsize(output_file))
    except FileNotFoundError:
        pass
    if filesize == 0:
        print('run tree %s'%(fastafile))
        run_parsi(fastafile)

os.system('mv %s/*.parsi.fasta.out* %s/'%(parsi_folder,tree_folder))
os.system('python treetopdf.py')
os.system('sh treetopdf.sh')
################################################### END ########################################################
################################################### SET PATH ########################################################
# merge MG results
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/withinHS'
alloutput = []
os.system('rm %s'%('%s/all.withinHS.snp.sum'%(outputdir)))
for files in glob.glob('%s/*withinHS.snp.sum'%(outputdir)):
    donor = os.path.split(files)[-1].split('.withinHS.snp.sum')[0]
    for lines in open(files, 'r'):
        alloutput.append('%s\t' % (donor) + lines)

f1 = open('%s/all.withinHS.snp.sum'%(outputdir),'w')
f1.write(''.join(alloutput))
f1.close()

outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/other'
alloutput = []
os.system('rm %s'%('%s/all.other.snp.sum'%(outputdir)))
for files in glob.glob('%s/*other.snp.sum'%(outputdir)):
    donor = os.path.split(files)[-1].split('.other.snp.sum')[0]
    for lines in open(files, 'r'):
        alloutput.append('%s\t' % (donor) + lines)

f1 = open('%s/all.other.snp.sum'%(outputdir),'w')
f1.write(''.join(alloutput))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# vcfprocess.py
# vcf processing per cluster
import glob
import os
input_script_vcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcfpro_round2'
output_dir = '/scratch/users/anniz44/genomes/donor_species/'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round*/'

os.system('rm -rf %s'%(input_script_vcf))
vcf_name = '.all*raw.vcf'
try:
    os.mkdir(input_script_vcf)
except IOError:
    pass

all_vcf_file=glob.glob(os.path.join(output_dir + '/vcf_round2/merge/','*%s'%(vcf_name)))
print(all_vcf_file)
for vcf_file in all_vcf_file:
    filesize = 0
    try:
        filesize = int(os.path.getsize(vcf_file + '.filtered.snpfreq.txt'))
        if filesize == 0:
            donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
            f1 = open(os.path.join(input_script_vcf, '%s.sh' % (donor_species)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n')
            f1.write('python %s/../vcf_process.py -i \"%s\" -s %s -o %s -cluster %s\n' % (
                input_script_vcf, genome_root, input_script_vcf, output_dir, donor_species))
            f1.close()
    except FileNotFoundError:
        donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
        f1 = open(os.path.join(input_script_vcf, '%s.sh'%(donor_species)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n')
        f1.write('python %s/../vcf_process.py -i \"%s\" -s %s -o %s -cluster %s\n'%(input_script_vcf,genome_root, input_script_vcf, output_dir,donor_species))
        f1.close()

f1 = open(os.path.join(input_script_vcf, '../allvcfprocessing.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run: sh %s/../allvcfprocessing.sh'%(input_script_vcf))

################################################### END ########################################################
################################################### SET PATH ########################################################
# remove rec per cluster
import glob
import os
vcf_format = '.filtered.vcf'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
vcf_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge'

input_script_vcf = os.path.join(input_script,'vcfremoverec')

rec_set = ['EsCo','PB','BL','BA','TuSa','BiPs','PaDi']
os.system('rm -rf %s'%(input_script_vcf))
try:
    os.mkdir(input_script_vcf)
except IOError:
    pass

all_vcf_file=glob.glob(os.path.join(vcf_folder,'*%s'%(vcf_format)))
Donor_species = set()
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
    if donor_species not in Donor_species:
        Donor_species.add(donor_species)
        parsi_output = os.path.join(os.path.split(vcf_file)[0], '%s.all.parsi.fasta' % (donor_species))
        filesize = 0
        try:
            filesize = int(os.path.getsize(parsi_output))
        except FileNotFoundError:
            pass
        if filesize == 0:
            if any(item in donor_species for item in rec_set):
                rec_cutoff = 50000
                contig_cutoff = 10000
            else:
                rec_cutoff = 5000
                contig_cutoff = 5000
            if donor_species == '1_PB_IBD_0_clustercluster1':
                # lots of big MGEs, small contigs are main contigs
                contig_cutoff = 5000
            f1 = open(os.path.join(input_script_vcf, '%s.sh'%(donor_species)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n')
            f1.write('python %s/SNPfilter_WGS_withdonor.py -i %s -vcf %s -cluster %s -s %s -rec %s -contig %s\n'%(input_script,vcf_folder, vcf_format, donor_species,input_script,rec_cutoff,contig_cutoff))
            f1.close()

f1 = open(os.path.join(input_script, 'allvcfremoverec.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    #f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.write('sh %s\n' % (sub_scripts))

f1.write('jobmit dnds.sh dnds.sh\n')
f1.close()
print('please run: sh %s/allvcfremoverec.sh'%(input_script))

################################################### END ########################################################
################################################### SET PATH ########################################################
# remove rec per cluster -> SNPfilter_WGS_all.py
import glob
import os
vcf_format = '.filtered.vcf'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
vcf_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge'
input_script_vcf = os.path.join(input_script,'vcfremoverec')
rec_cutoff = 5000
contig_cutoff = 10000

#rec_set = ['BL','BA','PaDi']
rec_set = []
os.system('rm -rf %s'%(input_script_vcf))
try:
    os.mkdir(input_script_vcf)
except IOError:
    pass

all_vcf_file=glob.glob(os.path.join(vcf_folder,'*%s'%(vcf_format)))
Donor_species = set()

for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
    if donor_species not in Donor_species:
        Donor_species.add(donor_species)
        parsi_output = os.path.join(os.path.split(vcf_file)[0], '%s.all.parsi.fasta' % (donor_species))
        filesize = 0
        try:
            filesize = int(os.path.getsize(parsi_output))
        except FileNotFoundError:
            pass
        if filesize == 0:
            f1 = open(os.path.join(input_script_vcf, '%s.sh'%(donor_species)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n')
            f1.write('python %s/SNPfilter_WGS.py -i %s -vcf %s -cluster %s -s %s -rec %s -contig %s\n'%(input_script,vcf_folder, vcf_format, donor_species,input_script,rec_cutoff,contig_cutoff))
            #f1.write('python %s/SNPfilter_WGS_cluster.py -i %s -vcf %s -cluster %s -s %s\n'%(input_script,vcf_folder, vcf_format, donor_species,input_script))
            f1.close()

f1 = open(os.path.join(input_script, 'allvcfremoverec.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    #f1.write('sh %s\n' % (sub_scripts))

f1.write('#jobmit dnds.sh dnds.sh\n')
f1.close()
print('please run: sh %s/allvcfremoverec.sh'%(input_script))

################################################### END ########################################################
################################################### SET PATH ########################################################

# check clonal population
from Bio import SeqIO
from Bio.Seq import Seq
import os
fasta = 'core_gene_alignment.aln'
os.system('prodigal -q -i %s -d %s.fna' % (fasta, fasta))
output = []
for record in SeqIO.parse(fasta, 'fasta'):
    record_id = str(record.id)
    record_seq = str(record.seq)
    output.append('>%s\n%s\n' % (record_id, record_seq))
    break

f1 = open(fasta + '.subset', 'w')
f1.write(''.join(output))
f1.close()
# check SNP position
#os.system('snp-sites -v %s.subset > %s.subset.vcf'%(output,output))
# check homologous in core genome
#os.system('blastn -subject %s.subset  -query %s.subset -out %s.subset.homologous.txt -outfmt 6 -max_target_seqs 100'%(output,output,output))
# check gene level SNP distribution
os.system('prodigal -q -i %s.subset -d %s.subset.fna' % (fasta, fasta))
#os.system('blastn -subject %s  -query %s -out %s.gene.divergence.txt -outfmt 6 -window_size 500 -perc_identity 80 -max_target_seqs 100'%(fasta,fasta,fasta))
os.system('%s -makeudb_usearch %s.subset.fna -output %s.subset.fna.udb' %
          ('usearch', fasta, fasta))
os.system('%s -ublast %s.fna -db %s.subset.fna.udb  -evalue 1e-2 -accel 0.5 -blast6out %s -threads 5 -strand plus -maxaccepts 500' %
          ('usearch', fasta, fasta, fasta + '.gene.divergence.txt'))

# check gene level SNP distribution of a suspicious subset
fasta = 'core_gene_alignment.aln'
subset = ['ao_BiPs_g0003','ao_BiPs_g0006']
output = []
for record in SeqIO.parse(fasta, 'fasta'):
    record_id = str(record.id)
    if record_id in subset:
        record_seq = str(record.seq)
        output.append('>%s\n%s\n' % (record_id, record_seq))

f1 = open(fasta + '.check.subset', 'w')
f1.write(''.join(output))
f1.close()
os.system('prodigal -q -i %s.check.subset -d %s.check.subset.fna' % (fasta, fasta))
os.system('%s -ublast %s.check.subset.fna -db %s.subset.fna.udb  -evalue 1e-2 -accel 0.5 -blast6out %s -threads 5 -strand plus -maxaccepts 500' %
          ('usearch', fasta, fasta, fasta + '.check.gene.divergence.txt'))
################################################### END ########################################################
################################################### SET PATH ########################################################
# set cutoff for recombination windows to calculate NS ratio
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics

input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge'
vcf_name = '.raw.vcf.filtered.snp.txt'
windowsize_set = [810000, 729000, 656100, 590490, 531441, 478296, 430466, 387419, 348677, 313809,
                  282428, 254185, 228766, 205889, 185300, 166770, 150093, 135083, 121574, 109416,
                  98474, 88626, 79763, 71786, 64607, 58146, 52331, 47097, 42387, 38148, 34333, 30899,
                  27809, 25028, 22525, 20272, 18244, 16419, 14777, 13299, 11969, 10772, 9694, 8724,
                  7851, 7065, 6358, 5722, 5149, 4634, 4170, 3753, 3377, 3039, 2735, 2461, 2214, 1992,
                  1792, 1612, 1450, 1305, 1174, 1056, 950, 855, 769, 692, 622, 559, 503, 452, 406, 365,
                  328, 295, 265, 238, 214, 192, 172, 154, 138, 124, 111, 99, 89, 80, 72, 64, 57, 51, 45,
                  40, 36, 32, 28, 25, 22, 19, 17, 15, 13, 11, 9, 8, 7, 6, 5, 4, 3, 2, 1]

try:
    os.mkdir(output_dir_merge + '/summary')
except IOError:
    pass

def windowing_concat(Seq_N):
    for windowsize in windowsize_set:
        N_S_set = [0, 0, 0]
        total_NS = ['']
        POS_old = 0
        for CHR in Seq_N:
            # calculate interval
            try:
                total_length = CHR.split('size')[1]
            except IndexError:
                total_length = CHR.split('length_')[1].split('_cov')[0]
            total_length = int(total_length)
            if POS_old + total_length >= windowsize:
                total_interval = int(total_length / windowsize) + 1
                total_NS += [''] * (total_interval)
            # windowing SNPs
            POS_set = Seq_N[CHR][0]
            NS_set = Seq_N[CHR][1]
            for i in range(0, len(POS_set)):
                POS = POS_set[i] + POS_old
                loci_POS = int(POS / windowsize)
                if total_NS[loci_POS] == '':
                    total_NS[loci_POS] = NS_set[i]
            POS_old += total_length
        N_S_set[0] += total_NS.count('N')
        N_S_set[1] += total_NS.count('S')
        try:
            N_S_set[2] = N_S_set[0] / N_S_set[1]
        except ZeroDivisionError:
            N_S_set[2] = 'N_only'
        Output.append('%s\t%s\t%s\t%s\t%s\t\n'%(donor_species,windowsize,
                                                N_S_set[0],N_S_set[1],N_S_set[2]))

def windowing(Seq_N):
    for windowsize in windowsize_set:
        N_S_set = [0, 0, 0]
        total_NS = ['']
        for CHR in Seq_N:
            # calculate interval
            try:
                total_length = CHR.split('size')[1]
            except IndexError:
                total_length = CHR.split('length_')[1].split('_cov')[0]
            total_length = int(total_length)
            if total_length >= windowsize:
                total_interval = int(total_length / windowsize) + 1
                total_NS += [''] * (total_interval)
            # windowing SNPs
            POS_set = Seq_N[CHR][0]
            NS_set = Seq_N[CHR][1]
            for i in range(0, len(POS_set)):
                POS = POS_set[i]
                loci_POS = int(POS / windowsize)
                if total_NS[loci_POS] == '':
                    total_NS[loci_POS] = NS_set[i]
        N_S_set[0] += total_NS.count('N')
        N_S_set[1] += total_NS.count('S')
        try:
            N_S_set[2] = N_S_set[0] / N_S_set[1]
        except ZeroDivisionError:
            N_S_set[2] = 'N_only'
        Output.append('%s\t%s\t%s\t%s\t%s\t\n'%(donor_species,windowsize,
                                                N_S_set[0],N_S_set[1],N_S_set[2]))

def readSNPfile(vcf_file):
    Seq_N = dict()
    for lines in open(vcf_file,'r'):
        lines_set = lines.replace('\n', '').replace('\r', '').split('\t')
        CHR = lines_set[0]
        POS = int(lines_set[1])
        N_S = lines_set[-2]
        if 'None' not in N_S and 'S' not in N_S:
            N_S = 'N'
        Seq_N.setdefault(CHR,[[],[]])
        Seq_N[CHR][0].append(POS)
        Seq_N[CHR][1].append(N_S)
    return Seq_N

# before remove rec
Output = []
all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
for vcf_file in all_vcf_file:
    if 'BaFr_clustercluster2' not in vcf_file:
        print(vcf_file)
        donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0]
        Seq_N = readSNPfile(vcf_file)
        windowing(Seq_N)

foutput = open(output_dir_merge + '/summary/all.donor.species.NSratio.txt', 'w')
foutput.write('donor_species\twindowsize\tN\tS\tNS_ratio\t\n')
foutput.write(''.join(Output))
foutput.close()

# after remove rec
vcf_name = '.raw.vcf.filtered.vcf.final.snp.txt'
Output = []
all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
for vcf_file in all_vcf_file:
    if 'BaFr_clustercluster2' not in vcf_file:
        print(vcf_file)
        donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0]
        Seq_N = readSNPfile(vcf_file)
        windowing(Seq_N)

foutput = open(output_dir_merge + '/summary/all.donor.species.NSratio.removerec.txt', 'w')
foutput.write('donor_species\twindowsize\tN\tS\tNS_ratio\t\n')
foutput.write(''.join(Output))
foutput.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# cluster ratio -> cluster_ratio.py
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

species_list = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/species_order.new.txt'
denovo_fasta = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.genome.gene.faa'
output_list = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/cluster_ratio.txt'

def split_line(lines):
    return lines.split('\n')[0].split('\t')

def find_species(record_id):
    species = record_id.split('_')[0]
    if species == '1':
        species = record_id.split('_')[1]
    return species

# load all species
All_species = dict()
for lines in open(species_list,'r'):
    if not lines.startswith('species_short'):
        lines_set = split_line(lines)
        All_species.setdefault(lines_set[2],[0,0])

# load denovo fasta
for record in SeqIO.parse(denovo_fasta, 'fasta'):
    record_id = str(record.id)
    All_species[find_species(record_id)][0]+=1

# load clustered denovo fasta
for record in SeqIO.parse(denovo_fasta + '.cluster.aa', 'fasta'):
    record_id = str(record.id)
    All_species[find_species(record_id)][1] += 1

# output cluster ratio
Output = []
Output.append('species\tseq_before_cluster\tseq_after_cluster\tcluster_ratio\t\n')
for species in All_species:
    before_cluster,aftercluster = All_species[species]
    cluster_ratio = 0
    if aftercluster > 0:
        cluster_ratio = before_cluster/aftercluster
    Output.append('%s\t%s\t%s\t%s\t\n'%(species,before_cluster,aftercluster,cluster_ratio))

foutput = open(output_list , 'w')
foutput.write(''.join(Output))
foutput.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# annotate all genomes -> annotate_allgenome.py, annotate_allgenome.sh
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

allreference_list = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/all.reference.list'
#allfasta = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/*/*.noHM.fasta.faa')
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/annotate_genome'

#os.system('head %s/*.all*.raw.vcf | grep \'##reference=\' > %s/all.reference.list'%(output_dir_merge,output_dir_merge))
try:
    os.mkdir(input_script_sub)
except IOError:
    pass

def annotation(all_filter_gene_fasta_file,pre_cluster = ''):
    all_filter_gene_fasta_file = all_filter_gene_fasta_file + '.cluster.aa'
    if pre_cluster!= '':
        os.system('#%s -makeudb_usearch %s -output %s.udb' %
                  ('usearch', pre_cluster, pre_cluster))
        os.system('%s -ublast %s -db %s.udb  -evalue 1e-2 -accel 0.5 -blast6out %s -threads 2'%
                  ('usearch', all_filter_gene_fasta_file,pre_cluster, all_filter_gene_fasta_file + '.ref.out.txt'))
    # run metacyc
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/mit_alm/database/metacyc/protseq.fsa'
    cmds = ("%s blastp --query %s --db %s.dmnd --out %s.metacyc.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %('diamond',all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    f1 = open(os.path.join(input_script_sub, 'metacyc.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    # run eggnog
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/eggnog/xaa.hmm'
    cmds = ('%s --tblout %s.eggnog.1.txt --cpu 40 -E %s %s %s\n') %('hmmsearch', all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.1.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xab.hmm'
    cmds = ('%s --tblout %s.eggnog.2.txt --cpu 40 -E %s %s %s\n') % (
        'hmmsearch',
        all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.2.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xac.hmm'
    cmds = ('%s --tblout %s.eggnog.3.txt --cpu 40 -E %s %s %s\n') % (
        'hmmsearch',
        all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.3.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmds))
    f1.close()
    # run kegg
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/kegg/kofam/profiles/prokaryote/prokaryote.hmm'
    cmds = ('%s --tblout %s.kegg.txt --cpu 40 -E %s %s %s\n') %('hmmsearch', all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'kegg.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    # run prokka
    cmdsprokka = 'py37\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n' + \
                 'prokka --kingdom Bacteria --force --outdir %s/prokka_%s  --protein %s --locustag Bacter %s/%s\n' % \
                 (output_dir_merge + '/summary', os.path.split(all_filter_gene_fasta_file)[-1],
                  all_filter_gene_fasta_file,
                  output_dir_merge + '/summary',
                  os.path.split(all_filter_gene_fasta_file)[-1].replace('.faa', '.fna'))
    f1 = open(os.path.join(input_script_sub, 'prokka.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmdsprokka))
    f1.close()
    # all scripts
    f1 = open(os.path.join(input_script, 'allannotategenome.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.sh')):
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()
    print('please run %s/%s'%(input_script,'allannotategenome.sh'))

allfasta = []
for lines in open(allreference_list, 'r'):
    if lines.startswith('##reference=file:'):
        # set database
        database_file = lines.split('##reference=file:')[1].split('\n')[0]
        allfasta.append(database_file)
        input_fasta = database_file + '.faa'
        input_fasta_dna = database_file + '.fna'
        try:
            f1 = open(input_fasta, 'r')
        except FileNotFoundError:
            os.system('prodigal -q -i %s -a %s' % (database_file, input_fasta))
        try:
            f1 = open(input_fasta_dna, 'r')
        except FileNotFoundError:
            os.system('prodigal -q -i %s -d %s' % (database_file, input_fasta_dna))

# merge all genes
output_gene = output_dir_merge + '/summary/all.genome.gene.faa'
output_fasta = []
try:
    f1 = open(output_gene,'r')
except IOError:
    change_name = set()
    for database in allfasta:
        database += '.faa'
        donor_species = os.path.split(database)[-1].split('.all')[0]
        donor_species_new = donor_species.replace('_clustercluster', '_CL')
        print(donor_species,donor_species_new)
        for record in SeqIO.parse(database, 'fasta'):
            record_id = str(record.id)
            temp_line = '%s\t%s\t%s\t' % (donor_species, donor_species_new, record_id)
            record_id = 'C_%s_G_%s' % (record_id.split('_')[1], record_id.split('_')[-1])
            output_fasta.append('>%s__%s\n%s\n' % (donor_species_new, record_id, str(record.seq)))
            temp_line += '%s__%s\t\n' % (donor_species_new, record_id)
            change_name.add(temp_line)
    f1 = open(output_gene, 'w')
    f1.write(''.join(output_fasta))
    f1.close()
    f1 = open(output_gene + '.changename.txt', 'w')
    f1.write(''.join(list(change_name)))
    f1.close()
    # run cluster
    cutoff = 0.9
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                   % ('usearch', output_fasta, cutoff, output_fasta,
                      output_fasta, 40))
    os.system(cmd_cluster)

annotation(output_gene)

output_gene = output_dir_merge + '/summary/all.genome.gene.fna'
output_fasta = []
try:
    f1 = open(output_gene,'r')
except IOError:
    change_name = set()
    for database in allfasta:
        database += '.fna'
        donor_species = os.path.split(database)[-1].split('.all')[0]
        donor_species_new = donor_species.replace('_clustercluster', '_CL')
        print(donor_species,donor_species_new)
        for record in SeqIO.parse(database, 'fasta'):
            record_id = str(record.id)
            temp_line = '%s\t%s\t%s\t' % (donor_species, donor_species_new, record_id)
            record_id = 'C_%s_G_%s' % (record_id.split('_')[1], record_id.split('_')[-1])
            output_fasta.append('>%s__%s\n%s\n' % (donor_species_new, record_id, str(record.seq)))
            temp_line += '%s__%s\t\n' % (donor_species_new, record_id)
            change_name.add(temp_line)
    f1 = open(output_gene, 'w')
    f1.write(''.join(output_fasta))
    f1.close()
    cutoff = 0.9
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                   % ('usearch', output_fasta, cutoff, output_fasta,
                      output_fasta, 40))
    os.system(cmd_cluster)

################################################### END ########################################################
################################################### SET PATH ########################################################
# core flexible of de novo genes, split HS annotation within lineage and across lineages -> core_flexible2.py
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'

allgenome = output_dir_merge + '/summary/all.genome.gene.faa.uc'
allgenome_denovo = output_dir_merge + '/summary/all.denovo.gene.faa.uc'
allgenome_HS = output_dir_merge + '/summary/all.selected.gene.faa.uc'

def cluster_uc(cluster_input):
    Clusters = dict()
    for lines in open(cluster_input, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        cluster = line_set[1]
        record_name = line_set[8].split(' ')[0]
        Clusters.setdefault(cluster, set())
        Clusters[cluster].add(record_name)
    return Clusters

# clusters based on all genes
Clusters = cluster_uc(allgenome)

def genus_list(species_list):
    return set([i.replace('BL','BiLo').replace('BA','BiAd')[0:2] for i in species_list])

def species_cluster(Clusters):
    Clusters_species = dict()
    Gene_cluster = dict()
    Output = []
    Output.append('record\tspecies_num\tallspecies\tgenus_num\t\n')
    for cluster in Clusters:
        allrecord = Clusters[cluster]
        Clusters_species.setdefault(cluster,set())
        for record in allrecord:
            record_new = record.split('_')[0]
            if record_new == '1':
                record_new = record.split('_')[1]
            record_new = record_new.replace('PB','PaDi').replace('Bfragilis','BaFr')
            Clusters_species[cluster].add(record_new)
        for record in allrecord:
            Gene_cluster.setdefault(record,Clusters_species[cluster])
    for record in Gene_cluster:
        allspecies = Gene_cluster[record]
        genus = genus_list(allspecies)
        Output.append('%s\t%s\t%s\t%s\t\n'%(record,len(allspecies),','.join(list(allspecies)),len(genus)
                                        ))
    f1 = open(allgenome + '.species.sum','w')
    f1.write(''.join(Output))
    f1.close()
    return [Clusters_species,Gene_cluster]

Clusters_species,Gene_cluster = species_cluster(Clusters)

def cluster_species_sum(fastaname,Gene_cluster):
    Clusters = cluster_uc(fastaname)
    Output = []
    Output.append('record\tspecies_num\tallspecies\tgenus_num\t\n')
    for cluster in Clusters:
        allrecord = Clusters[cluster]
        for record in allrecord:
            record_new = record.replace('.donor.'+record.split('.')[-1].split('__')[0],'')
            if '1_PaDi_IBD' in record_new:
                record_new = record_new.replace('1_PaDi_IBD','1_PB_IBD')
            allspecies = Gene_cluster[record_new]
            genus = genus_list(allspecies)
            Output.append('%s\t%s\t%s\t%s\t\n' % (record, len(allspecies), ','.join(list(allspecies)),len(genus)
                                                  ))
    f1 = open(fastaname + '.species.sum', 'w')
    f1.write(''.join(Output))
    f1.close()

cluster_species_sum(allgenome_denovo,Gene_cluster)
cluster_species_sum(allgenome_HS,Gene_cluster)
################################################### END ########################################################
################################################### SET PATH ########################################################
# core flexible of de novo genes eggnog ->core_flexible_eggnog.py
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details'
genus_cutoff = 2 # more than 2 genera as core
allgenome = output_dir_merge + '/summary/all.genome.gene.faa'
allgenome_denovo = output_dir_merge + '/summary/all.denovo.gene.faa'
allgenome_HS = output_dir_merge + '/summary/all.selected.gene.faa'

def genus_list(species_list):
    return set([i.replace('BL','BiLo').replace('BA','BiAd')[0:2] for i in species_list])

def function_species(fastaname):
    species_file = fastaname + '.uc.species.sum'
    eggnog_file = fastaname + '.cluster.aa.all.eggnog.sum'
    Species_fun = dict()
    for lines in open(species_file, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        record_name = line_set[0]
        species_num = line_set[1]
        species_name = line_set[2].split(',')
        tag_genus = 'False'
        if len(species_name) > 1:
            genus = genus_list(species_name)
            if len(genus) > genus_cutoff:
                tag_genus = 'True'
        Species_fun.setdefault(record_name, [species_num,tag_genus])
    Output = []
    for lines in open(eggnog_file, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        if lines.startswith('cluster'):
            Output.append('\t'.join(line_set) + '\tspecies_num\tcross_genus\n')
        else:
            if line_set[6] == '':
                line_set.pop(6)
            record_name = line_set[4]
            species_num, tag_genus = Species_fun[record_name]
            Output.append('\t'.join(line_set) + '\t%s\t%s\n'%(species_num, tag_genus))
    f1 = open(fastaname + '.cluster.aa.all.eggnog.sum.species.sum', 'w')
    f1.write(''.join(Output))
    f1.close()

function_species(allgenome_HS)
function_species(allgenome_denovo)
function_species(allgenome)

# split HS annotation within lineage and across lineages
def HS_lineage(filename,HS_lineagefasta):
    HS_lineage_set = set()
    for record in SeqIO.parse(HS_lineagefasta, 'fasta'):
        record_id = str(record.id)
        HS_lineage_set.add(record_id)
    Output_within = set()
    Output_across = set()
    for lines in open(filename,'r'):
        if lines.split('\t')[4] in HS_lineage_set:
            Output_within.add(lines)
        else:
            Output_across.add(lines)
    f1 = open(filename + '.within.sum', 'w')
    f1.write(''.join(list(Output_within)))
    f1.close()
    f1 = open(filename + '.across.sum', 'w')
    f1.write(''.join(list(Output_across)))
    f1.close()

HS_annotation_sum = output_dir_merge + '/summary/all.selected.gene.faa.cluster.aa.all.eggnog.sum.species.sum'
allgenome_HS_lineage = output_dir_merge + '/summary/all.selected.gene.faa'

HS_lineage(HS_annotation_sum,allgenome_HS_lineage)

################################################### END ########################################################
################################################### SET PATH ########################################################
# tree phylogenetic diversity
from Bio import Phylo
from Bio.Phylo import BaseTree
import numpy

def to_distance(tree):
    # phylogenetic diversity of a tree, sum of branch_length
    dis_all = 0
    for parent in tree.find_clades(terminal=False, order="level"):
        for child in parent.clades:
            if child.branch_length:
                dis_all += child.branch_length
    return dis_all

filename_tree = '1_BL_IBD_0_clustercluster1.donor.D77.all.parsi.fasta.out.tree'
tree = Phylo.read(filename_tree, "newick")
to_distance(tree)

################################################### END ########################################################
################################################### SET PATH ########################################################
# move MG results
import os,glob
result_dir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
snpsumdir = '/scratch/users/anniz44/genomes/donor_species/MG/snpsum/'
covsumdir = '/scratch/users/anniz44/genomes/donor_species/MG/covsum/'

try:
    os.mkdir(snpsumdir)
except IOError:
    pass

try:
    os.mkdir(covsumdir)
except IOError:
    pass

for snpfile in glob.glob('%s/*.snp.sum'%(result_dir)):
    try:
        outputdir = os.path.split(snpfile)[-1].split('.IN.')[1].split('.snp.sum')[0]
        outputdir = os.path.join('%s/%s'%(snpsumdir,outputdir))
        os.mkdir(outputdir)
    except IOError:
        pass
    os.system('mv %s %s/'%(snpfile,outputdir))


for snpfile in glob.glob('%s/*.coverage.sum'%(result_dir)) + glob.glob('%s/*.depth.MV.sum'%(result_dir)) + glob.glob('%s/*.depth.sum'%(result_dir)):
    try:
        outputdir = os.path.split(snpfile)[-1].split('.fasta.')[1].split('.raw.vcf')[0]
        outputdir = os.path.join('%s/%s'%(covsumdir,outputdir))
        os.mkdir(outputdir)
    except IOError:
        pass
    os.system('mv %s %s/'%(snpfile,outputdir))



for snpfile in glob.glob('%s/finished/*.zip'%(result_dir)):
    try:
        outputdir = os.path.split(snpfile)[-1].split('.fasta.')[1].split('.raw.vcf')[0]
        outputdir = os.path.join('%s/finished/%s'%(result_dir,outputdir))
        os.mkdir(outputdir)
    except IOError:
        pass
    os.system('mv %s %s/'%(snpfile,outputdir))
################################################### END ########################################################
################################################### SET PATH ########################################################
# PE truncation -> PE_trunc.py
# all SNPs truncation
import os,glob
snp_dir = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'

def loadHS(sum_file):
    HS_gene = dict()
    for lines in open(sum_file, 'r'):
        if not lines.startswith("#donor_species"):
            lines_set = lines.split('\n')[0].split('\t')
            if lines_set[-1] == 'True':
                lineage, gene = lines_set[0:2]
                HS_gene.setdefault(lineage,set())
                HS_gene[lineage].add(gene)
    return HS_gene

def load_snp(snpfile,alloutput):
    donor_species = os.path.split(snpfile)[-1].split('.raw.vcf.filtered.vcf.final.snp.txt')[0].replace('.all','')
    #if donor_species in HS_gene_within:
    for lines in open(snpfile,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            genename, POS,N_or_S,AAchange = lines_set[-4:]
            genenamenew, POS = lines_set[:2]
            if '*' in AAchange:
                if '*' == AAchange[0]:
                    Trunc_SNP = lines_set[2]
                else:
                    Trunc_SNP = lines_set[3]
                if donor_species in HS_gene_within and genename in HS_gene_within[donor_species]:
                    alloutput.add('%s\t%s\t%s\t%s\t%s\twithinHS\t%s\t%s\n'%(donor_species,genenamenew,POS,AAchange[0],AAchange[1],Trunc_SNP,genename))
                #elif genename in HS_gene_all[donor_species]:
                    #alloutput.add('%s\t%s\t%s\t%s\t%s\tacrossHS\t%s\t%s\n'%(donor_species,genenamenew,POS,AAchange[0],AAchange[1],Trunc_SNP,genename))
                else:
                    alloutput.add('%s\t%s\t%s\t%s\t%s\tothers\t%s\t%s\n' % (
                    donor_species, genenamenew, POS, AAchange[0], AAchange[1],Trunc_SNP,genename))
    else:
        print('no HS in %s'%(donor_species))
    return alloutput

# load HS genes
HS_gene_within = loadHS('%s/summary/all.species.txt'%(snp_dir))
#HS_gene_all = loadHS('%s/summary/all.species.txt.High_select2.txt'%(snp_dir))

alloutput=set()
snp_folder = glob.glob('%s/*.donor.*.raw.vcf.filtered.vcf.final.snp.txt' % (snp_dir))
print(snp_folder)
for snpfile in snp_folder:
    alloutput = load_snp(snpfile,alloutput)

f1 = open('%s/summary/all.species.txt.Truncated.sum'%(snp_dir), 'w')
f1.write(''.join(list(alloutput)))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# genomic diversity among strains carried different PE alleles of the same lineage
import os,glob
allmultisnp = '/scratch/users/anniz44/genomes/donor_species/MG/summary/all.withinHS.snp.genename.depthcheck.count.multigenotype.sum'
changename = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.selected.gene.faa.changename.txt'
genome_dir1 = '/scratch/users/anniz44/genomes/donor_species/selected_species/round1/'
genome_dir2 = '/scratch/users/mit_alm/IBD_Evo/BA/Assembly_for_gene_flow/'
genome_dir3 = '/scratch/users/mit_alm/IBD_Evo/BL/Assembly_for_gene_flow/'
genome_dir4 = '/scratch/users/mit_alm/IBD_Evo/PB/Assembly_for_gene_flow/'
snpfolder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/fastani'
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/fastani'
try:
    os.mkdir(outputdir)
except IOError:
    pass

try:
    os.mkdir(outputscript)
except IOError:
    pass

def loadsnpfile(snpfile, positionset):
    genomeset = dict()
    genomeorder = dict()
    for lines in open(snpfile,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        if lines.startswith('CHR'):
            i = 9
            for genome in lines_set[9:]:
                genomeorder.setdefault(i,genome)
                i += 1
        else:
            CHRPOS = '%s\t%s'%(lines_set[0],lines_set[1])
            if CHRPOS in positionset:
                major, minor = lines_set[2:4]
                genomeset.setdefault(CHRPOS,[[major,set()],[minor,set()]])
                i = 9
                for alleles in lines_set[9:]:
                    if alleles == major:
                        genomeset[CHRPOS][0][-1].add(genomeorder[i])
                    elif alleles == minor:
                        genomeset[CHRPOS][1][-1].add(genomeorder[i])
                    i += 1
    return genomeset

def findfasta(lineage,genomename):
    if lineage.startswith('1_'):
        fasta = glob.glob('%s/%s/scaffolds.fasta' % (genome_dir2, genomename))+\
        glob.glob('%s/%s/scaffolds.fasta' % (genome_dir3, genomename))+\
        glob.glob('%s/%s/scaffolds.fasta' % (genome_dir4, genomename))
    else:
        fasta = glob.glob('%s/%s/fasta/%s_final.scaffolds.fasta' % (genome_dir1, genomename.split('_g')[0], genomename))
    return fasta[0]

def run_fastani(genomeset,lineage):
    for CHRPOS in genomeset:
        genomeset1,genomeset2=genomeset[CHRPOS]
        genomeset1out = set()
        genomeset2out = set()
        CHRPOSset = CHRPOS.split('\t')
        for genome in genomeset1[-1]:
            genomeset1out.add(findfasta(lineage,genome))
        for genome in genomeset2[-1]:
            genomeset2out.add(findfasta(lineage,genome))
        outputname1 = '%s__%s__%s__%s'%(lineage,CHRPOSset[0],CHRPOSset[1],genomeset1[0])
        outputname2 = '%s__%s__%s__%s' % (lineage, CHRPOSset[0], CHRPOSset[1], genomeset2[0])
        f1 = open('%s/%s.list'%(outputdir,outputname1),'w')
        f1.write('\n'.join(list(genomeset1out))+'\n')
        f1.close()
        f1 = open('%s/%s.list' % (outputdir, outputname2), 'w')
        f1.write('\n'.join(list(genomeset2out))+'\n')
        f1.close()
        command = ('fastANI --rl %s --ql %s -o %s \n' %
                    ('%s/%s.list'%(outputdir,outputname1), '%s/%s.list'%(outputdir,outputname1),
                     '%s/%s.fastaniout' % (outputdir, outputname1)))
        command += ('fastANI --rl %s --ql %s -o %s \n' %
                   ('%s/%s.list' % (outputdir, outputname2), '%s/%s.list' % (outputdir, outputname2),
                    '%s/%s.fastaniout' % (outputdir, outputname2)))
        f1 = open('%s/%s.sh' % (outputscript, outputname1), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (command))
        f1.close()

newgenename = dict()
for lines in open(changename,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    newgenename.setdefault(lines_set[3],lines_set[2])

donorspecies = dict()
for lines in open(allmultisnp,'r'):
    if not lines.startswith('subgroupnew'):
        lines_set = lines.split('\t')
        lineage = lines_set[0].split('__')[0].replace('_CL','_clustercluster')
        donorspecies.setdefault(lineage,
                                set())
        donorspecies[lineage].add('%s\t%s'%('_'.join(newgenename[lines_set[0].split(' ')[0]].split('_')[:-1]),
                                            lines_set[0].split(' ')[1]))
for lineage in donorspecies:
    snpfile = '%s/%s.all.parsi.fasta.sum.txt'%(snpfolder,lineage)
    positionset = donorspecies[lineage]
    genomeset = loadsnpfile(snpfile, positionset)
    run_fastani(genomeset,lineage)

f1 = open(os.path.join(outputscript, '../allfastani.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(outputscript, '*.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/../allfastani.sh'%(outputscript))

################################################### END ########################################################
################################################### SET PATH ########################################################
# sum up fastani for genomic diversity among strains carried different PE alleles of the same lineage
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/fastani'

alloutput = []
for files in glob.glob('%s/*.fastaniout'%(outputdir)):
    donor_species_allele = os.path.split(files)[-1].split('.fastaniout')[0]
    for lines in open(files,'r'):
        alloutput.append(
            '%s\t%s'%(donor_species_allele,lines)
        )

f1 = open('%s/allfastaniout.sum'%(outputdir),'w')
f1.write('allele\tgenome1\tgenome2\tani\tgenemap\ttotalgenes\n')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# sum up ref genome gene length and non-ORF length -> cluster_length.py
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
import statistics

assemblyfolder1 = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/'
assemblyfolder2 = '/scratch/users/anniz44/genomes/donor_species/vcf_round1/co-assembly/'
assemblyfolder3 = '/scratch/users/anniz44/genomes/donor_species/WGS_old/WGS/vcf_round1/co-assembly/'
fasta = '.fasta.noHM.fasta'
outputfile = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/clonal_genelength.txt'

alloutput = []
for files in glob.glob('%s/*/*%s'%(assemblyfolder1,fasta))+\
             glob.glob('%s/*/*%s'%(assemblyfolder2,fasta))+\
             glob.glob('%s/*/AlOn*%s'%(assemblyfolder3,fasta))+\
             glob.glob('%s/*/BaFr*%s'%(assemblyfolder3,fasta)):
    donor_species = os.path.split(files)[-1].split(fasta)[0]
    species = donor_species.split('_')[0]
    print(files)
    if species!= 'prokka':
        total_length = 0
        total_genelength = []
        gene_num = 0
        for record in SeqIO.parse(files, 'fasta'):
            total_length += len(str(record.seq))
        filesize = 0
        try:
            filesize = int(os.path.getsize(files + '.fna'))
        except FileNotFoundError:
            pass
        if filesize == 0:
            os.system('prodigal -q -i %s -d %s -a %s' % (files, files + '.fna',files + '.faa'))
        for record in SeqIO.parse(files + '.fna', 'fasta'):
            total_genelength.append(len(str(record.seq)))
            gene_num += 1
        alloutput.append('%s\t%s\t%s\t%.1f\t%s\n'%(species,donor_species.replace('_PB_','_PaDi_'),
                                           total_length - sum(total_genelength),
                                           statistics.mean(total_genelength),gene_num))


f1 = open('%s'%(outputfile),'w')
f1.write('species\tdonor_species\tnonORFlength\tavgORFlength\ttotal_gene_num\n')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# sum up same SNPs across lineages -> sameSNP.py
import os,glob
SNPfiles = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/*.all.parsi.fasta.sum.txt')
HSfileout = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.species.txt.sameSNP.multipledonor.txt'
all_filter_gene_fasta_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.denovo.gene.faa'

cutoff = 0.7
cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.%s.cluster.aa -uc %s.%s.uc -threads %s\n'
                   % ('usearch', all_filter_gene_fasta_file, cutoff, all_filter_gene_fasta_file,cutoff,
                      all_filter_gene_fasta_file, cutoff,40))
#os.system(cmd_cluster)
cluster_file = all_filter_gene_fasta_file + '.%s.uc'%(cutoff)
#cluster_file = all_filter_gene_fasta_file + '.uc'
Clusters = dict()
for lines in open(cluster_file, 'r'):
    line_set = lines.split('\n')[0].split('\t')
    cluster = line_set[1]
    record_name = line_set[8]
    Clusters.setdefault(record_name, cluster)

HSsum = dict()
genetocluster = dict()
for aSNPfile in SNPfiles:
    donor_species = os.path.split(aSNPfile)[-1].split('.all.parsi.fasta.sum.txt')[0]
    CL = donor_species.split('.donor.')[0]
    donor = donor_species.split('.donor.')[1]
    donor_species_new = donor_species.replace('_clustercluster', '_CL').replace('PB_', 'PaDi_')
    for lines in open(aSNPfile,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\n')[0].split('\t')
            genename = lines_set[5]
            POS = lines_set[6]
            N_S = lines_set[4]
            if genename == 'None':
                genename = lines_set[0]
                POS = lines_set[1]
                newgenename = '%s\t%s\t%s\tNone\t%s' % (CL, genename, POS,N_S)
            else:
                genename2 = '%s__C_%s_G_%s' % (donor_species_new, genename.split('_')[1], genename.split('_')[-1])
                cluster = Clusters.get(genename2)
                newgenename = '%s\t%s' % (cluster,POS)
                genetocluster.setdefault(newgenename,set())
                genetocluster[newgenename].add('%s\t%s\t%s\t%s\t%s' % (CL, genename, POS,genename2,N_S))
            HSsum.setdefault(newgenename, set())
            HSsum[newgenename].add(donor_species)

alloutput = []
for newgenename in HSsum:
    donorcutoff = 3
    donornum = len(HSsum[newgenename])
    if donornum >= donorcutoff:
        if newgenename in genetocluster:
            for oldgenename in genetocluster[newgenename]:
                alloutput.append('%s\t%s\t%s\t%s\t\n' % (oldgenename,newgenename, donornum,
                                                         ';'.join(list(HSsum[newgenename]))
                                                         ))
        else:
            alloutput.append('%s\tNone\tNone\t%s\t%s\t\n'%(newgenename,donornum,
                                                         ';'.join(list(HSsum[newgenename]))
                                                           ))

f1 = open('%s'%(HSfileout),'w')
f1.write('Refgenome\tgenename\tPOS\tnewgenename\tN_or_S\tcluster\tPOS2\tnum_donor\tdonor_list\t\n')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# test MG filter
import os,glob
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf')
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/testMG.sh'
fastaformat = '.fasta'
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/\n'%(vcffile,fastaformat)

cmds += 'python SNPsum_meta.py  -o /scratch/users/anniz44/genomes/donor_species/\n'

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmds))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# re run MG filtering
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/refilterMG/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/refilterMG.sh'
jub_num = 200
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    i += 1
    if (i%jub_num) == 0:
        f1 = open('%s/refilterMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
        f1.write(''.join(cmds))
        f1.close()
        cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
        cmdsall += 'jobmit %s %s\n'%('%s/refilterMG.%s.sh'%(outputscript,int(i/jub_num)),
                                     'refilterMG.%s.sh' % (int(i / jub_num)))
    tempoutputfolder = '%s/refilterMG_%s'%(outputscript,int(i/jub_num))
    try:
        os.mkdir(tempoutputfolder)
    except IOError:
        pass
    cmds += 'unzip %s -d %s\n'%(vcffile,tempoutputfolder)
    cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n'%(tempoutputfolder,outputdir)
    cmds += 'rm -r %s/scratch/\n'%(tempoutputfolder)
    newvcffile = os.path.join(outputdir,os.path.split(vcffile)[-1].split('.zip')[0])
    if '.fasta' in newvcffile:
        fastaformat = '.fasta'
    else:
        fastaformat = '.1.fq'
    cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/\n'%(newvcffile,fastaformat)
    cmds += 'rm %s\n'%(newvcffile)

f1 = open('%s/refilterMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
f1.write(''.join(cmds))
f1.close()
cmdsall += 'jobmit %s %s\n'%('%s/refilterMG.%s.sh'%(outputscript,int(i/jub_num)),
                                     'refilterMG.%s.sh' % (int(i / jub_num)))

cmdsall += '#python SNPsum_meta.py  -o /scratch/users/anniz44/genomes/donor_species/\n'
f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# MG trunc filtering
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncMG/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncMG.sh'
jub_num = 200
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    if '.fasta' in vcffile:
        fastaformat = '.fasta'
    else:
        fastaformat = '.1.fq'
    samplename = os.path.split(vcffile)[-1].split(fastaformat)[0]
    donorspecies = os.path.split(vcffile)[-1].split(fastaformat + '.')[1].split('.raw.vcf')[0]
    try:
        f1 = open('%s/../truncsum/%s/%s.IN.%s.snp.sum'%(outputdir,donorspecies,samplename,donorspecies),'r')
    except IOError:
        try:
            f1 = open('%s/../truncsumHS/%s/%s.IN.%s.snp.sum' % (outputdir, donorspecies, samplename, donorspecies), 'r')
        except IOError:
            i += 1
            if (i%jub_num) == 0:
                f1 = open('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
                f1.write(''.join(cmds))
                f1.close()
                cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
                cmdsall += 'jobmit %s %s\n'%('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),
                                             'truncMG.%s.sh' % (int(i / jub_num)))
            tempoutputfolder = '%s/truncMG_%s'%(outputscript,int(i/jub_num))
            try:
                os.mkdir(tempoutputfolder)
            except IOError:
                pass
            cmds += 'unzip %s -d %s\n'%(vcffile,tempoutputfolder)
            cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n'%(tempoutputfolder,outputdir)
            cmds += 'rm -r %s/scratch/\n'%(tempoutputfolder)
            newvcffile = os.path.join(outputdir,os.path.split(vcffile)[-1].split('.zip')[0])
            cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta_trunc.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/\n'%(newvcffile,fastaformat)
            cmds += 'rm %s\n'%(newvcffile)

f1 = open('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
f1.write(''.join(cmds))
f1.close()
cmdsall += 'jobmit %s %s\n'%('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),
                                     'truncMG.%s.sh' % (int(i / jub_num)))

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# MG trunc sum
import os,glob
inputdir = glob.glob('/scratch/users/anniz44/genomes/donor_species/MG/truncsumHS/*/')+\
           glob.glob('/scratch/users/anniz44/genomes/donor_species/MG/truncsum/*/')
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/'

alloutput = []
for a_folder in inputdir:
    donorspecies = os.path.split(os.path.split(a_folder)[-2])[-1]
    for a_file in glob.glob('%s/*.snp.sum'%(a_folder)):
        samplename = os.path.split(a_file)[-1].split('.IN.')[0]
        for lines in open(a_file):
            if not lines.startswith('CHR'):
                alloutput.append('%s\t%s\t%s'%(donorspecies,samplename,lines))

f1 = open('%s/all.truncsum.txt'%(outputdir),'w')
f1.write('Lineage\tSample\tCHR\tPOS\tGene\tGenePOS\twithinHS\tMajor_ALT\tMinor_ALT\tMajor_ALT_freq\tMinor_ALT_freq\tRef_aa\tSNP_aa\n')
f1.write(''.join(alloutput))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# MG PE trunc all samples at trunc snp
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncMGHS/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncMGHS.sh'
snp_file = '/scratch/users/anniz44/genomes/donor_species/MG/summary/all.truncsum.filterfreq.trunclist.withHS.txt'
jub_num = 50
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

def load_snp(snpfile):
    SNP_list = dict()
    SNP_list2 = dict()
    if snpfile != 'None':
        for lines in open(snpfile,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            lineage,CHRPOS,Major,Minor,Gene,Genepos = lines_set[:6]
            SNP_list.setdefault(lineage,set())
            SNP_list[lineage].add(CHRPOS)
            SNP_list2.setdefault(CHRPOS, [Minor,Gene,Genepos])
    return[SNP_list,SNP_list2]

SNP_list,SNP_list2 = load_snp(snp_file)

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    if '.fasta' in vcffile:
        fastaformat = '.fasta'
    else:
        fastaformat = '.1.fq'
    samplename = os.path.split(vcffile)[-1].split(fastaformat)[0]
    donorspecies = os.path.split(vcffile)[-1].split(fastaformat + '.')[1].split('.raw.vcf')[0]
    if donorspecies in SNP_list:
        try:
            f1 = open('%s/../truncsumHSnew/%s/%s.IN.%s.snp.sum'%(outputdir,donorspecies,samplename,donorspecies),'r')
        except IOError:
            i += 1
            if (i%jub_num) == 0:
                f1 = open('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
                f1.write(''.join(cmds))
                f1.close()
                cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
                cmdsall += 'jobmit %s %s\n'%('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),
                                             'truncMG.%s.sh' % (int(i / jub_num)))
            tempoutputfolder = '%s/truncMG_%s'%(outputscript,int(i/jub_num))
            try:
                os.mkdir(tempoutputfolder)
            except IOError:
                pass
            cmds += 'unzip %s -d %s\n'%(vcffile,tempoutputfolder)
            cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n'%(tempoutputfolder,outputdir)
            cmds += 'rm -r %s/scratch/\n'%(tempoutputfolder)
            newvcffile = os.path.join(outputdir,os.path.split(vcffile)[-1].split('.zip')[0])
            cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta_trunc.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/ -snplist %s\n'%(newvcffile,fastaformat,snp_file)
            cmds += 'rm %s\n'%(newvcffile)

f1 = open('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
f1.write(''.join(cmds))
f1.close()

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# MG trunc all SNPs before trunc snp
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncallsnp/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncallsnp.sh'
snp_file = '/scratch/users/anniz44/genomes/donor_species/MG/summary/all.truncsum.txt'
jub_num = 100
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

def load_snp(snpfile):
    SNP_list = dict()
    SNP_list2 = dict()
    if snpfile != 'None':
        for lines in open(snpfile,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            lineage, sample, CHR, POS, Gene, Genepos = lines_set[:6]
            SNP_list.setdefault(lineage, set())
            SNP_list[lineage].add(sample)
            lineagesample = '%s\t%s' % (lineage, sample)
            SNP_list2.setdefault(lineagesample, [])
            SNP_list2[lineagesample].append(' '.join([CHR, POS, Gene, Genepos]))
    return[SNP_list,SNP_list2]

SNP_list,SNP_list2 = load_snp(snp_file)

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    if '.fasta' in vcffile:
        fastaformat = '.fasta'
    else:
        fastaformat = '.1.fq'
    samplename = os.path.split(vcffile)[-1].split(fastaformat)[0]
    donorspecies = os.path.split(vcffile)[-1].split(fastaformat + '.')[1].split('.raw.vcf')[0]
    if samplename in SNP_list.get(donorspecies):
        lineagesample = '%s\t%s' % (donorspecies, samplename)
        snplistall = ';'.join(SNP_list2[lineagesample])
        try:
            f1 = open('%s/../truncallsnps/%s/%s.IN.%s.snp.sum'%(outputdir,donorspecies,samplename,donorspecies),'r')
        except IOError:
            i += 1
            if (i%jub_num) == 0:
                f1 = open('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
                f1.write(''.join(cmds))
                f1.close()
                cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
                cmdsall += 'jobmit %s %s\n'%('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),
                                             'truncMG.%s.sh' % (int(i / jub_num)))
            tempoutputfolder = '%s/truncMG_%s'%(outputscript,int(i/jub_num))
            try:
                os.mkdir(tempoutputfolder)
            except IOError:
                pass
            cmds += 'unzip %s -d %s\n'%(vcffile,tempoutputfolder)
            cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n'%(tempoutputfolder,outputdir)
            cmds += 'rm -r %s/scratch/\n'%(tempoutputfolder)
            newvcffile = os.path.join(outputdir,os.path.split(vcffile)[-1].split('.zip')[0])
            cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta_trunc.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/ -trunclist \'%s\'\n'%(
                newvcffile,fastaformat,snplistall)
            cmds += 'rm %s\n'%(newvcffile)

f1 = open('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
f1.write(''.join(cmds))
f1.close()

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# MG trunc allSNPs sum (before the truncation sites)
import os,glob
inputdir = glob.glob('/scratch/users/anniz44/genomes/donor_species/MG/truncallsnps/*/')
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/'
snpfile = '/scratch/users/anniz44/genomes/donor_species/MG/summary/all.truncsum.withdetails.nosingle.txt'

lineagesample = set()
lineagesample2=set()
for lines in open(snpfile,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    lineage_sample = lines_set[1]
    gene = lines_set[7]
    lineagesample2.add('%s %s'%(lineage_sample,gene))
    lineagesample.add(lineage_sample)

alloutput = []
for a_folder in inputdir:
    donorspecies = os.path.split(os.path.split(a_folder)[-2])[-1]
    for a_file in glob.glob('%s/*.snp.sum'%(a_folder)):
        samplename = os.path.split(a_file)[-1].split('.IN.')[0]
        lineage_sample = '%s %s'%(donorspecies,samplename)
        if lineage_sample in lineagesample:
            for lines in open(a_file):
                if not lines.startswith('CHR'):
                    lines_set = lines.split('\n')[0].split('\t')
                    gene = lines_set[2]
                    if '%s %s'%(lineage_sample,gene) in lineagesample2:
                        alloutput.append('%s\t%s\t%s'%(donorspecies,samplename,lines))
        else:
            pass

f1 = open('%s/all.truncsum.allsnps.txt'%(outputdir),'w')
f1.write('Lineage\tSample\tCHR\tPOS\tGene\tGenePOS\twithinHS\tMajor_ALT\tMinor_ALT\tMajor_ALT_freq\tMinor_ALT_freq\tRef_aa\tSNP_aa\n')
f1.write(''.join(alloutput))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# MG trunc HS sum
import os,glob
inputdir = glob.glob('/scratch/users/anniz44/genomes/donor_species/MG/truncsumHSnew/*/')
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/'

alloutput = []
for a_folder in inputdir:
    donorspecies = os.path.split(os.path.split(a_folder)[-2])[-1]
    for a_file in glob.glob('%s/*.snp.sum'%(a_folder)):
        samplename = os.path.split(a_file)[-1].split('.IN.')[0]
        for lines in open(a_file):
            if not lines.startswith('CHR'):
                alloutput.append('%s\t%s\t%s'%(donorspecies,samplename,lines))

f1 = open('%s/all.truncsum.HS.txt'%(outputdir),'w')
f1.write('Lineage\tSample\tCHR\tPOS\tGene\tGenePOS\twithinHS\tMajor_ALT\tMinor_ALT\tMajor_ALT_freq\tMinor_ALT_freq\tRef_aa\tSNP_aa\n')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# MG indel filtering
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/indelMG/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/indelMG.sh'
jub_num = 200
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    if '.fasta' in vcffile:
        fastaformat = '.fasta'
    else:
        fastaformat = '.1.fq'
    samplename = os.path.split(vcffile)[-1].split(fastaformat)[0]
    donorspecies = os.path.split(vcffile)[-1].split(fastaformat + '.')[1].split('.raw.vcf')[0]
    try:
        f1 = open('%s/../indelsum/%s/%s.IN.%s.snp.sum'%(outputdir,donorspecies,samplename,donorspecies),'r')
    except IOError:
        i += 1
        if (i%jub_num) == 0:
            f1 = open('%s/indelMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
            f1.write(''.join(cmds))
            f1.close()
            cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
            cmdsall += 'jobmit %s %s\n'%('%s/indelMG.%s.sh'%(outputscript,int(i/jub_num)),
                                         'indelMG.%s.sh' % (int(i / jub_num)))
        tempoutputfolder = '%s/indelMG_%s'%(outputscript,int(i/jub_num))
        try:
            os.mkdir(tempoutputfolder)
        except IOError:
            pass
        cmds += 'unzip %s -d %s\n'%(vcffile,tempoutputfolder)
        cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n'%(tempoutputfolder,outputdir)
        cmds += 'rm -r %s/scratch/\n'%(tempoutputfolder)
        newvcffile = os.path.join(outputdir,os.path.split(vcffile)[-1].split('.zip')[0])
        cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta_indel.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/\n'%(newvcffile,fastaformat)
        cmds += 'rm %s\n'%(newvcffile)

f1 = open('%s/indelMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
f1.write(''.join(cmds))
f1.close()
cmdsall += 'jobmit %s %s\n'%('%s/indelMG.%s.sh'%(outputscript,int(i/jub_num)),
                                     'indelMG.%s.sh' % (int(i / jub_num)))

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# MG indel HS sum
import os,glob
inputdir = glob.glob('/scratch/users/anniz44/genomes/donor_species/MG/indelsum/*/')
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/'
HS_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.species.txt'

HS_set = dict()

for lines in open(HS_file):
    if 'True' in lines:
        lines_set=lines.split('\t')
        donorspecies,gene = lines_set[0:2]
        donorspecies=donorspecies.split('.donor')[0]
        HS_set.setdefault(donorspecies,set())
        HS_set[donorspecies].add(gene)

alloutput = []
for a_folder in inputdir:
    donorspecies = os.path.split(os.path.split(a_folder)[-2])[-1]
    donorspecies = donorspecies.replace('_PB_','_PaDi_')
    if donorspecies in HS_set:
        HS_set_gene = HS_set[donorspecies]
        for a_file in glob.glob('%s/*.snp.sum'%(a_folder)):
            samplename = os.path.split(a_file)[-1].split('.IN.')[0]
            for lines in open(a_file):
                if not lines.startswith('CHR'):
                    lines_set = lines.split('\t')
                    gene = lines_set[2]
                    if gene in HS_set_gene:
                        alloutput.append('%s\t%s\t%s'%(donorspecies,samplename,lines))
    else:
        print(donorspecies)

f1 = open('%s/all.indelsum.HS.txt'%(outputdir),'w')
f1.write('Lineage\tSample\tCHR\tPOS\tGene\tGenePOS\tMajor_ALT\tMinor_ALT\tMajor_ALT_freq\tMinor_ALT_freq\n')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# MG PE mapping length + all gene mapping results
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncPEsum/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncPEsum.sh'
outputsumdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/'
jub_num = 200
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

trunc_output  = outputsumdir + '/all.lineage.PESNP.sum'
f1 = open(trunc_output,'w')
f1.write('Lineage\tSample\tNo.PESNPs\n')
f1.close()

trunc_output_gene  = outputsumdir + '/all.lineage.genemapping.sum'
f1 = open(trunc_output_gene,'w')
f1.write('Lineage\tSample\tcluster_count_PE\tcluster_count_noPE\tgene_count_PE\tgene_count_noPE\n')
f1.close()

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    if '.fasta' in vcffile:
        fastaformat = '.fasta'
    else:
        fastaformat = '.1.fq'
    samplename = os.path.split(vcffile)[-1].split(fastaformat)[0]
    donorspecies = os.path.split(vcffile)[-1].split(fastaformat + '.')[1].split('.raw.vcf')[0]
    i += 1
    if (i % jub_num) == 0:
        f1 = open('%s/truncPEsum.%s.sh' % (outputscript, int(i / jub_num)), 'a')
        f1.write(''.join(cmds))
        f1.close()
        cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
        cmdsall += 'jobmit %s %s\n' % ('%s/truncPEsum.%s.sh' % (outputscript, int(i / jub_num)),
                                       'truncPEsum.%s.sh' % (int(i / jub_num)))
    tempoutputfolder = '%s/truncPEsum_%s' % (outputscript, int(i / jub_num))
    try:
        os.mkdir(tempoutputfolder)
    except IOError:
        pass
    cmds += 'unzip %s -d %s\n' % (vcffile, tempoutputfolder)
    cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n' % (tempoutputfolder, outputdir)
    cmds += 'rm -r %s/scratch/\n' % (tempoutputfolder)
    newvcffile = os.path.join(outputdir, os.path.split(vcffile)[-1].split('.zip')[0])
    cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta_genesum.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/\n' % (
    newvcffile, fastaformat)
    cmds += 'rm %s\n' % (newvcffile)

f1 = open('%s/truncPEsum.%s.sh'%(outputscript,int(i/jub_num)),'a')
f1.write(''.join(cmds))
f1.close()

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# MG snp median depth
import os,glob
outputfile = '/scratch/users/anniz44/genomes/donor_species/MG/summary/all.lineage.depth.sum'
alldepthfiles = glob.glob('/scratch/users/anniz44/genomes/donor_species/MG/covsum/*/*.raw.vcf.depth.sum')

alloutput = []
for files in alldepthfiles:
    sample = os.path.split(files)[-1].split('.fasta')[0].split('.1.fq')[0]
    lineage = os.path.split(os.path.split(files)[0])[-1]
    totalloci = 0
    mediandepth = 0
    depth_to_loci = dict()
    alldepth = list()
    for lines in open(files,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\n')[0].split('\t')
            depth = int(lines_set[1])
            if depth > 0:
                num_loci = int(lines_set[2])
                alldepth.append(depth)
                depth_to_loci.setdefault(depth,num_loci)
                totalloci += num_loci
    newtotalloci = 0
    for depth in sorted(alldepth):
        newtotalloci += depth_to_loci[depth]
        if newtotalloci >= totalloci/2:
            mediandepth = depth
            break
    alloutput.append('%s\t%s\t%.1f\t%s\n'%(sample,lineage,mediandepth,totalloci))

f1 = open('%s'%(outputfile),'w')
f1.write('sample\tlineage\tmediandepth\ttotalloci\n')
f1.write(''.join(alloutput))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# check all depth
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/MGdepth/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/MGdepth.sh'
subset = ['SRR2938117.fasta.BA_clustercluster6','SRR2912781.fasta.1_BL_IBD_1_clustercluster4',
          'SRR2938117.fasta.1_BL_IBD_1_clustercluster4','SRR2938079.fasta.1_BL_IBD_1_clustercluster4']
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    if any(vcfname in vcffile for vcfname in subset):
        i += 1
        if (i%2000) == 0:
            f1 = open('%s/MGdepth.%s.sh'%(outputscript,int(i/2000)),'a')
            f1.write(''.join(cmds))
            f1.close()
            cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
            cmdsall += 'jobmit %s %s\n'%('%s/MGdepth.%s.sh'%(outputscript,int(i/2000)),
                                         'MGdepth.%s.sh' % (int(i / 2000)))
        tempoutputfolder = '%s/MGdepth_%s'%(outputscript,int(i/2000))
        try:
            os.mkdir(tempoutputfolder)
        except IOError:
            pass
        cmds += 'unzip %s -d %s\n'%(vcffile,tempoutputfolder)
        cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n'%(tempoutputfolder,outputdir)
        cmds += 'rm -r %s/scratch/\n'%(tempoutputfolder)
        newvcffile = os.path.join(outputdir,os.path.split(vcffile)[-1].split('.zip')[0])
        if '.fasta' in newvcffile:
            fastaformat = '.fasta'
        else:
            fastaformat = '.1.fq'
        cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta.allloci.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/\n'%(newvcffile,fastaformat)
        cmds += 'rm %s\n'%(newvcffile)

f1 = open('%s/MGdepth.%s.sh'%(outputscript,int(i/2000)),'a')
f1.write(''.join(cmds))
f1.close()
cmdsall += 'jobmit %s %s\n'%('%s/MGdepth.%s.sh'%(outputscript,int(i/2000)),
                                     'MGdepth.%s.sh' % (int(i / 2000)))

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# simulation to loci info
import os
alloutput = []
for lines in open('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.genome.gene.fnatruncation.simulation.withmutfreq.txt','r'):
    GC = lines.split('\t')[0]
    truncation_loci = lines.split('\t')[-1].split('\n')[0]
    truncation_loci = truncation_loci.replace('[','').replace(']','')
    for loci in truncation_loci.split(','):
        alloutput.append('%s\t%s\n'%(GC,loci.replace(' ','')))

f1 = open('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.genome.gene.fnatruncation.simulation.withmutfreq.loci.txt','w')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# start codon near truncations not used
from Bio import SeqIO
from Bio.Seq import Seq
import glob
import os

fasta = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.genome.gene.fna'
truncfile = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.truncsum.withdetails.addtruncpercentage.uniqueGC.txt'
codon_check = 1 # search 10 codons
gene_trunc = dict()
for lines in open(truncfile,'r'):
    if not lines.startswith('gene_name'):
        lines_set = lines.split('\n')[0].split('\t')
        gene_name = lines_set[0]
        GC,gene_length,trunc_loci = lines_set[-3:]
        trunc_loci = float(trunc_loci)
        gene_trunc.setdefault(gene_name,[GC,int(gene_length),set()])
        gene_trunc[gene_name][-1].add(trunc_loci)

outputfasta = []
outputstartcodon = []
outputstartcodonsum = dict()
for record in SeqIO.parse(fasta, 'fasta'):
    record_id = str(record.id)
    record_seq = str(record.seq)
    if record_id in gene_trunc:
        outputfasta.append('>%s\n%s\n'%(record_id,record_seq))
        GC, gene_length, trunc_lociset = gene_trunc[record_id]
        for trunc_loci in trunc_lociset:
            trunc_loci2 = int(trunc_loci*gene_length/100)
            codon_loci = int(trunc_loci2/3)
            start_loci = int(max(1,codon_loci-codon_check)) # start with first codon
            stop_loci = min(int(gene_length/3),codon_loci + codon_check)
            record_seq_sub = [record_seq[i*3:(i*3+3)] for i in range(start_loci,stop_loci+1)]
            trunc_loci_interset = int(trunc_loci/5)
            outputstartcodonsum.setdefault(trunc_loci_interset,[0,0])
            if 'ATG' in record_seq_sub or 'CAT' in record_seq_sub:
                outputstartcodon.append('%s\t%s\t%s\tTrue\n' % (record_id, GC, trunc_loci))
                outputstartcodonsum[trunc_loci_interset][0] += 1
            else:
                outputstartcodon.append('%s\t%s\t%s\tFalse\n' % (record_id, GC, trunc_loci))
                outputstartcodonsum[trunc_loci_interset][1] += 1

outputsum = []
for trunc_loci_interset in outputstartcodonsum:
    outputsum.append('%s\t%s\t%s\n'%(trunc_loci_interset,
                                 outputstartcodonsum[trunc_loci_interset][0],
                                 outputstartcodonsum[trunc_loci_interset][1]))

f1 = open('%s.trunc.start.codon.sum.txt'%(fasta),'w')
f1.write('truncloci_per\twith_near_start_codon\twithout_near_start_codon\n')
f1.write(''.join(outputsum))
f1.close()

f1 = open('%s.trunc.fna'%(fasta),'w')
f1.write(''.join(outputfasta))
f1.close()

f1 = open('%s.trunc.start.codon.txt'%(fasta),'w')
f1.write('gene_name\tGC\ttruncloci_per\tnear_start_codon\n')
f1.write(''.join(outputstartcodon))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# check wrong vcf_round2
import os,glob
shfileall = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/*.sh')
filesize = 29

cmds = '#!/bin/bash\nsource ~/.bashrc\n'
for shfile in shfileall:
    try:
        filesize = int(os.path.getsize(shfile))
    except FileNotFoundError:
        pass
    if filesize == 29:
        donor_species = os.path.split(shfile)[-1].split('.all.flt.snp.vcf.sh')[0]
        print(donor_species)
        cmds += ('rm /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/%s.*all.parsi.fasta\n'%(
            donor_species
        ))
        ref_file = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/%s.*.ref.raw.vcf'%(donor_species))
        if ref_file!= []:
            cmds += ('bcftools filter -T /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/%s.all.flt.snp.vcf.CHRPOS %s > %s.filter\n'%(donor_species,ref_file[0],ref_file[0]))
        else:
            print(donor_species,'no ref')
        cmds += ('sh /scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcfremoverec/%s.sh\n'%(donor_species))

f1 = open('checkalltree2.sh','w')
f1.write(cmds)
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# check finished vcf_round2
# vcfround2finish.py
import os,glob
parsiall = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/*.all.parsi.fasta')

for parsi in parsiall:
    donor_species = os.path.split(parsi)[-1].split('_cluster')[0]
    clonal_file = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/%s.*'%(donor_species))
    if clonal_file!=[]:
        os.system('mv %s /scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/finished/'%(' '.join(clonal_file)))
################################################### END ########################################################
################################################### SET PATH ########################################################
#convert tree to tree.pdf -> treetopdf.py
import os,glob
alltree = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/tree/*.tree')
cmds = '#!/bin/bash\nsource ~/.bashrc\n'
cmds += 'cd /scratch/users/anniz44/bin/pro/FigTree_v1.4.4/\n'
for tree in alltree:
    try:
        f1 = open('%s.pdf'%(tree),'r')
    except FileNotFoundError:
        cmds += ('./figtree -graphic PDF %s %s.pdf\n'%(tree,tree))

f1 = open('/scratch/users/anniz44/scripts/1MG/donor_species/assembly/treetopdf.sh','w')
f1.write(cmds)
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
#compare cluster difference -> clustercheck.py
import os,glob
newcluster = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/*.genome.cluster.txt')
oldclusterfolder = '/scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/finished/'

def donor_num(clusterfile):
    donor_cluster = dict()
    for lines in open(clusterfile, 'r'):
        lines_set = lines.split('\t')
        donor = lines_set[1].split('_')[0]
        donor_cluster.setdefault(donor, set())
        donor_cluster[donor].add(lines_set[2])
    return donor_cluster

newredo = []
for cluster in newcluster:
    species = os.path.split(cluster)[-1].split('.')[0]
    oldcluster = glob.glob('%s/%s'%(oldclusterfolder,os.path.split(cluster)[-1]))[0]
    donor_clusterold = donor_num(oldcluster)
    donor_clusternew = donor_num(cluster)
    for donor in donor_clusternew:
        if len(donor_clusternew[donor])<len(donor_clusterold[donor]):
            print('something wrong with %s %s'%(species,donor))
        elif len(donor_clusternew[donor])>len(donor_clusterold[donor]):
            newredo.append('%s\t%s\t\n'%(species,donor))

f1 = open('/scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/needredo.txt','w')
f1.write(''.join(newredo))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
#check whether sim cutoff works
import os,glob
donor_species_cutoff = []
for lines in open('total_SNP_cutoff.txt','r'):
    donor_species_cutoff.append(lines.split('\t')[0])

output = []
for lines in open('all.species.txt','r'):
    if lines.split('\t')[0] in donor_species_cutoff and int(lines.split('\t')[4])<=2:
        output.append(lines.replace('False', 'FalseFalse').replace('True', 'FalseFalse'))
    else:
        output.append(lines)

f1 = open('all.species.cutoff.txt','w')
f1.write(''.join(output))
f1.close()
