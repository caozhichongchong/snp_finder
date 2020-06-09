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
donor_list =["cx","bk","av","ao","an","am","af","aa","8306RM","7652PC","6937AY","6587MW","6383DG","4117FE","3790QQ","3554DX","3321JW","3219KR"]
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

# round 1 run whole genomes fastq mapping to a reference genome/pan-genome
import glob
import os
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species/vcf_new'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/species/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/species/vcf_round1'
fastq_name = '_1.fastq'
#input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_new'
#input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/test'
#genome_dir = glob.glob('/scratch/users/mit_alm/IBD_Evo/Testing_Sets/*')
#output_dir = '/scratch/users/anniz44/genomes/donor_species/test/vcf_round1'
#fastq_name = '_1.fastq'

os.system('rm -rf %s'%(input_script_sub))

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
    os.mkdir(output_dir+'/bwa/0/old')
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


for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    if len(donor) == 2: # BN10 donors
        donor_species_genome = glob.glob(os.path.join(folder,'*.fasta'))
        sub_samples = []
        donor_species_fastq = glob.glob(os.path.join(folder, 'fastq/*' + fastq_name)) + glob.glob(os.path.join(folder, '*' + fastq_name))
        if donor_species_fastq != []:
            cmds = ''
            i = 0
            split_bash = 1
            if donor_species_genome!= []:
                donor_species_genome_size = dict()
                for fasta_file in donor_species_genome:
                    filesize = int(os.path.getsize(fasta_file))
                    donor_species_genome_size.setdefault(filesize,[])
                    donor_species_genome_size[filesize].append(fasta_file)
                print(donor_species_genome_size)
                max_filesize = sorted(donor_species_genome_size.keys())[-1]
                database = donor_species_genome_size[max_filesize][0]
                try:
                    f1 = open(database + '.1.bt2','r')
                except IOError:
                    os.system('bowtie2-build %s %s' %(database,database))
            else:
                donor_species_folder_all = os.path.join(folder, donor_species + '_allspades')
                for fastq_file in donor_species_fastq:
                    if '.all' + fastq_name not in fastq_file:
                        donor_species_genome = fastq_file.split(fastq_name)[0]
                        cmds += 'head -500000 %s%s >> %s/%s.all%s\n' % (
                        donor_species_genome,fastq_name, folder, donor_species,fastq_name)
                        cmds += 'tail -500000 %s%s >> %s/%s.all%s\n' % (
                        donor_species_genome, fastq_name,folder, donor_species,fastq_name)
                        cmds += 'head -500000 %s%s >> %s/%s.all%s\n' % (
                            donor_species_genome, fastq_name.replace('1','2'), folder, donor_species, fastq_name.replace('1','2'))
                        cmds += 'tail -500000 %s%s >> %s/%s.all%s\n' % (
                            donor_species_genome, fastq_name.replace('1','2'), folder, donor_species, fastq_name.replace('1','2'))
                cmds += 'spades.py --careful -1 %s/%s.all%s -2 %s/%s.all%s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' %\
                       (folder,donor_species,fastq_name,folder,donor_species,fastq_name.replace('1','2'), donor_species_folder_all)
                cmds += 'mv %s/scaffolds.fasta %s/%s.all.spades.fasta\n' %(donor_species_folder_all,folder,donor_species)
                database = '%s/%s.all.spades.fasta'%(folder,donor_species)
                cmds += 'bowtie2-build %s %s\n' %(database,database)
                split_bash = 0
            for files in donor_species_fastq:
                donor_species_dir_file = os.path.split(files)[-1]
                print(donor_species_dir_file)
                tempbamoutput = os.path.join(output_dir +'/bwa/0', donor_species_dir_file)
                try:
                    f1 = open(tempbamoutput + '.sorted.bam')
                except IOError:
                    files2 = files.replace(fastq_name, fastq_name.replace('1','2'))
                    cmds += 'bowtie2' + ' -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --threads %s --un-conc %s.unaligned -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                        min(40, 40),tempbamoutput, database, files, files2,'samtools', min(40, 40),
                        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
                        tempbamoutput)
                    sub_samples.append(tempbamoutput + '.sorted.bam')
                    cmds += 'samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
                        tempbamoutput, tempbamoutput)
                    cmds += 'rm -rf %s.*.unaligned\n' % (tempbamoutput)
                    cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
                    cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
                    if split_bash == 1:
                        i += 1
                        if i%5 == 0:
                            f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species,int(i/5))), 'a')
                            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                            f1.close()
                            cmds = ''
            f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species,int(i/5))), 'a')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            f1.close()

f1 = open(os.path.join(input_script, 'allnewvcf.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################

# run allnewvcf and finish
# round 1 clonal population selection vcf filering
import glob
import os
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species/vcf_new'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species'
genome_dir = '/scratch/users/anniz44/genomes/donor_species/species/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/species/vcf_round1/bwa/0'
output_notwanted = '/scratch/users/anniz44/genomes/donor_species/unwanted/outlier'
fastq_name = '_1.fastq'
genome_split = '_g'
#input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_new'
#input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/test'
#genome_dir = ('/scratch/users/mit_alm/IBD_Evo/Testing_Sets/')
#output_dir = '/scratch/users/anniz44/genomes/donor_species/test/vcf_round1/bwa/0'
#output_notwanted = '/scratch/users/mit_alm/IBD_Evo/Testing_Sets/unwanted'
#fastq_name = '_1.fastq'
#genome_split = '_'

try:
    os.mkdir(output_notwanted)
except IOError:
    pass


# set up cutoff
#total_coverage_cutoff = 0.7 #70% reads align to the ref genome
total_coverage_cutoff = 0.6
genome_avg_coverage_cutoff = 10


Donor_species = dict()
Coverage = []
Coverage.append('genome\treads_alignment_rate\tref_genome_length\taverage_coverage\twanted\t\n')

sh_files = glob.glob(os.path.join(input_script_sub,'*.sh'))
for sh_file in sh_files:
    Coverage1 = []
    Coverage2 = []
    os.system('grep \"bowtie2\" %s > %s/temp.sh'%(sh_file,input_script_sub))
    os.system('grep \"overall alignment rate\" %s.err > %s/temp.err' % (sh_file,input_script_sub))
    # load fastq name
    for lines in open('%s/temp.sh' % (input_script_sub),'r'):
        try:
            filename = lines.split('>')[1].split('\n')[0].split(fastq_name + '.bam')[0]
            filename = os.path.split(filename)[-1]
            Coverage1.append(filename)
        except IndexError:
            pass
    # load reads alignment rate
    for lines in open('%s/temp.err' % (input_script_sub),'r'):
        Coverage2.append(float(lines.split('%')[0]))
    i = 0
    for filename in Coverage1:
        # load average coverage of ref genome
        tempbamoutput = os.path.join(output_dir,filename + fastq_name )
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
        Donor_species.setdefault(donor_species,[[],0])
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
        temp_line = ('%s\t%s\t%s\t%s' % (filename, coverage_num, length_total, avg_coverage))
        if coverage_num < total_coverage_cutoff or avg_coverage < genome_avg_coverage_cutoff:
            # not qualified
            print(filename, coverage_num, avg_coverage)
            Donor_species[donor_species][0].append(filename)
            temp_line += ('\tunwanted\t\n')
        else:
            # qualified
            Donor_species[donor_species][1] += 1
            temp_line += ('\twanted\t\n')
        Coverage.append(temp_line)

cmd = ''
newcluster = 1
Donor_species_output = []
Donor_species_output.append('donor_species\tqualified\tnon-qualified\t\n')
for donor_species in Donor_species:
    non_qualified_list = Donor_species[donor_species][0]
    total_non_qualified = len(non_qualified_list)
    Donor_species_output.append('%s\t%s\t%s\t\n' % (
        donor_species,
        Donor_species[donor_species][1],
        total_non_qualified))
    if total_non_qualified >= 6:
        # new cluster
        new_cluster_file = '%s/%s_newcluster%s' %(genome_dir,donor_species,newcluster+1)
        try:
            os.mkdir(new_cluster_file)
        except IOError:
            pass
        print(donor_species,non_qualified_list,new_cluster_file)
        for filename in non_qualified_list:
            cmd += ('mv %s %s/\n' % (os.path.join(genome_dir, '*/%s.*' % (filename)), new_cluster_file))
            cmd += ('mv %s %s/\n' % (os.path.join(genome_dir, '*/%s' % (filename + fastq_name + '*' )), new_cluster_file))
            cmd += ('mv %s %s/\n' % (os.path.join(genome_dir, '*/%s' % (filename + fastq_name.replace('1','2') + '*' )), new_cluster_file))
            #cmd += ('rm -rf %s\n' % (os.path.join(output_dir, '*%s'+fastq_name+'*' % (filename))))
    else:
        for filename in non_qualified_list:
            cmd += ('mv %s %s/\n' % (os.path.join(genome_dir, '*/%s.*' % (filename)), output_notwanted))
            cmd += ('mv %s %s/\n' % (os.path.join(genome_dir, '*/%s'% (filename + fastq_name + '*' )), output_notwanted))
            cmd += ('mv %s %s/\n' % (os.path.join(genome_dir, '*/%s'% (filename + fastq_name.replace('1','2') + '*' )), output_notwanted))
            #cmd += ('rm -rf %s\n' % (os.path.join(output_dir, '*%s'+fastq_name+'*' % (filename))))

f1 = open(os.path.join(input_script, 'wholegenome.coverage.%s.sum'%(newcluster)), 'w')
f1.write(''.join(Coverage))
f1.close()

f1 = open(os.path.join(input_script, 'wholegenome.qualify.%s.sum.txt'%(newcluster)), 'w')
f1.write(''.join(Donor_species_output))
f1.close()

# check Donor_species_output first! maybe wrong reference
f1 = open(os.path.join(input_script, 'wholegenome.nonqualify.move.%s.sh'%(newcluster)), 'w')
f1.write(''.join(cmd))
f1.close()

# clean up
os.system('#rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round1/bwa/0/*.unaligned')
os.system('#rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round1/bwa/0/*raw.vcf')
os.system('#rm -rf /scratch/users/anniz44/genomes/donor_species/unwanted/outlier/*.fastq.*')
os.system('#rm -rf /scratch/users/anniz44/genomes/donor_species/test/vcf_round1/bwa/0/*.unaligned')
os.system('#rm -rf /scratch/users/anniz44/genomes/donor_species/test/vcf_round1/bwa/0/*raw.vcf')
os.system('rm -rf /scratch/users/anniz44/genomes/donor_species/species/vcf_round1/bwa/0/*raw.vcf')

################################################### END ########################################################
################################################### SET PATH ########################################################

# round 1 multiple whole genome samples mileup
import glob
import os
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species/vcf_merge_new'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/species/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/species/vcf_round1'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/species/vcf_merge'
fastq_name = '_1.fastq'

#input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/vcf_merge_new'
#input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/'
#genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')
#output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round1'
#output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge'
#fastq_name = '_1.fastq'

#input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_merge_new'
#input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/test'
#genome_dir = glob.glob('/scratch/users/mit_alm/IBD_Evo/Testing_Sets/*')
#output_dir = '/scratch/users/anniz44/genomes/donor_species/test/vcf_round1'
#output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/test/vcf_merge'
#fastq_name = '_1.fastq'

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(output_dir_merge)
except IOError:
    pass

try:
    os.mkdir(output_dir_merge+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir_merge+'/bwa/0')
except IOError:
    pass

i=0
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    if len(donor) == 2: # BN10 donors
        ref_genome = glob.glob(os.path.join(folder,'*.fasta.1.bt2'))
        tempbamoutput = os.path.join(output_dir_merge + '/bwa/0', donor_species + '.all')
        try:
            f1 = open(tempbamoutput+'.flt.snp.vcf','r')
        except IOError:
            if ref_genome != []:
                database = ref_genome[0].split('.1.bt2')[0]
                print(database)
                cmds = ''
                fastq_files = glob.glob(os.path.join(folder, '*'+fastq_name))
                if fastq_files == []:
                    fastq_files = glob.glob(os.path.join(folder, 'fastq/*' + fastq_name))
                sub_samples = []
                for fastq in fastq_files:
                    if 'all' + fastq_name not in fastq:
                        sub_samples.append(os.path.join(output_dir +'/bwa/0',os.path.split(fastq)[-1]+'.sorted.bam'))
                print(len(sub_samples),folder,fastq_files)
                cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -B -Ou -d3000 -f %s %s  | %s call --ploidy 1 --threads %s -m > %s.raw.vcf\n' % (
                    'bcftools', min(40, 40),database,
                    ' '.join(sub_samples), 'bcftools', min(40, 40), tempbamoutput)
                cmds += '%s filter --threads %s -s  LowQual %s.raw.vcf > %s.flt.vcf \n' % (
                    'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
                cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
                    'bcftools', tempbamoutput, tempbamoutput)
                cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
                f1 = open(os.path.join(input_script_sub, '%s.vcf.sh' % donor_species), 'w')
                f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                f1.close()
            else:
                print('database missing for %s'%(donor_species))


f1 = open(os.path.join(input_script, 'allnewmergevcf.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################

# round 1 SNP selection merge vcf filtering
# round 1 phylogenetic tree
import glob
import os
input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/vcf_merge_new'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/tree_round1'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round1'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge'
fastq_name = '_1.fastq'

input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species/vcf_merge_new'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species/tree_round1'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/species/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/species/vcf_round1'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/species/vcf_merge'
fastq_name = '_1.fastq'

#input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_merge_new'
#input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/tree_round1'
#input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/test'
#genome_dir = glob.glob('/scratch/users/mit_alm/IBD_Evo/Testing_Sets/*')
#output_dir = '/scratch/users/anniz44/genomes/donor_species/test/vcf_round1'
#output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/test/vcf_merge'
#fastq_name = '_1.fastq'

try:
    os.mkdir(input_script_sub)
except IOError:
    pass


# set up cutoff
median_coverage_cutoff = 10 # avg coverage in all samples
Sample_depth_cutoff = 5 # both forward and reverse reads cutoff in a sample
Major_alt_freq_cutoff = 0.9 # major alt freq in a sample
SNP_presence_cutoff = 0.66 # percentage of samples passing the above criteria
SNP_presence_sample_cutoff = 3 # num of samples passing the above criteria
SNP_depth_cutoff = 5 # both forward and reverse reads cutoff for SNP ALTs in a sample

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

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

def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = int(Subdepth[0]) + int(Subdepth[1])
    Subdepth_ALT = int(Subdepth[2]) + int(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return allels_set[1]
    else:
        return allels_set[0]

all_vcf_file=glob.glob(os.path.join(output_dir_merge + '/bwa/0','*.all.flt.snp.vcf'))
for vcf_file in all_vcf_file:
    Total = 0
    vcf_file_list = []
    vcf_file_POS = []
    vcf_file_POS.append('CHR\tPOS\tPOS_DIS\t\n')
    SNP_alignment = dict()
    SNP_alignment.setdefault('reference', '')
    donor_species = os.path.split(vcf_file)[-1].split('.all.flt.snp.vcf')[0]
    SNP_tree_cmd = []
    SNP_tree_cmd2 = []
    CHR_old = ''
    POS_old = 0
    try:
        f1 = open(vcf_file+'.filtered.fasta','r')
    except IOError:
        print(vcf_file)
        for lines in open(os.path.join(input_script_merge_sub,'%s.vcf.sh'%(donor_species)),'r'):
            if lines.startswith('bcftools mpileup '):
            # setup samples
                sample_set = lines.split('.fasta ')[1].split('\n')[0].split('  |')[0].split(' ')
                for samples in sample_set:
                    genomename = os.path.split(samples)[-1].split(fastq_name+'.sorted.bam')[0]
                    SNP_alignment.setdefault(genomename, '')
        for lines in open(vcf_file, 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\n')[0].split('\t')
                Total_qualify = 0
                Total_qualify_SNP = 0
                if Total == 0:
                    Total = len(lines_set) - 9
                Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                if Depth/Total >= median_coverage_cutoff:
                    # average depth in all samples cutoff
                    if  "INDEL" not in lines_set[7] \
                            and (lines_set[6] != 'LowQual'):
                        # SNP only
                        SNP = []
                        REF = lines_set[3]
                        allels_set = [REF]
                        if '.' not in lines_set[4]:
                            allels_set += lines_set[4].split(',')
                        Total_alleles = len(allels_set)
                        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
                        REF = curate_REF(allels_set, Depth4)
                        for Subdepth_all in lines_set[9:]:
                            SNP.append(REF) # set as reference
                            Allels_frq = [0, 0, 0, 0]
                            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                            Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
                            Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
                            total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
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
                            if total_sub_depth_forward >= Sample_depth_cutoff and \
                                    total_sub_depth_reverse >= Sample_depth_cutoff and \
                                    Major_ALT[1] / total_sub_depth >= Major_alt_freq_cutoff:
                                # forward and reverse of the SNP ALTs cutoff
                                # major alt frequency cutoff
                                    Total_qualify += 1
                                    if Major_ALT[0] != REF:
                                        SNP[-1] = Major_ALT[0]
                                        if total_sub_depth_forward - int(Subdepth_forward[0]) >= SNP_depth_cutoff and \
                                            total_sub_depth_reverse - int(Subdepth_reverse[0]) >= SNP_depth_cutoff and \
                                             Major_ALT[1] >= 2*SNP_depth_cutoff:
                                            Total_qualify_SNP += 1
                        if Total_qualify/Total >= SNP_presence_cutoff and Total_qualify >= SNP_presence_sample_cutoff and Total_qualify_SNP >= 1 and Total_qualify_SNP < Total_qualify :
                            # -> qualified SNP
                            # qualified samples cutoff
                            # at least 1 qualified SNP
                            # output lines and output major alt
                            CHR = lines_set[0]
                            POS = int(lines_set[1])
                            if CHR == CHR_old:
                                # same CHR
                                POS_DIS = abs(POS - POS_old)
                                vcf_file_POS.append('%s\t%s\t%s\t\n' % (CHR, POS, POS_DIS))
                                POS_old = POS
                            else:
                                # diff CHR first SNP
                                POS_old = 0
                                vcf_file_POS.append('%s\t%s\t%s\t\n' % (CHR, POS, POS))
                            CHR_old = CHR
                            vcf_file_list.append(lines)
                            i = 0
                            SNP_alignment['reference'] += REF
                            for genomename in SNP_alignment:
                                if genomename!='reference':
                                    SNP_alignment[genomename] += SNP[i]
                                    i+=1
        SNP_alignment_output = []
        for genomename in SNP_alignment:
            SNP_alignment_output.append('>%s\n%s\n'%(genomename,SNP_alignment[genomename]))
        vcf_file_filtered = open(vcf_file + '.filtered','w')
        vcf_file_filtered.write(''.join(vcf_file_list))
        vcf_file_filtered.close()
        vcf_file_filtered = open(vcf_file + '.filtered.POS', 'w')
        vcf_file_filtered.write(''.join(vcf_file_POS))
        vcf_file_filtered.close()
        vcf_file_filtered = open(vcf_file + '.filtered.fasta', 'w')
        vcf_file_filtered.write(''.join(SNP_alignment_output))
        vcf_file_filtered.close()
        fasta = vcf_file+'.filtered.fasta'
        SNP_tree_cmd.append('python /scratch/users/anniz44/scripts/maffttree/remove.duplicated.py -i %s\n' % (fasta))
        filesize = int(os.path.getsize(fasta))
        if filesize >= 10^7: #10Mb
            SNP_tree_cmd.append(
                'mafft --nuc --quiet --nofft --retree 2 --maxiterate 0 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                #'mafft --nuc --quiet --retree 2 --maxiterate 100 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                fasta, fasta)) # remove --adjustdirection
        else:
            SNP_tree_cmd.append(
                # 'mafft --nuc --quiet --nofft --retree 2 --maxiterate 0 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                'mafft --nuc --quiet --retree 2 --maxiterate 100 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                    fasta, fasta))  # remove --adjustdirection
        SNP_tree_cmd.append('FastTree -nt -quiet %s.mafft.align > %s.mafft.align.nwk\n' % (fasta, fasta))
        SNP_tree_cmd2.append(
            '#run_gubbins.py --threads 40  --tree_builder hybrid --use_time_stamp --prefix %s_gubbins --verbose %s.mafft.align\n' % (
                fasta, fasta))
        f1 = open(os.path.join(input_script_sub, '%s.tree.sh'%(donor_species)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n' + ''.join(SNP_tree_cmd)+ ''.join(SNP_tree_cmd2))
        f1.close()

f1 = open(os.path.join(input_script, 'alltree.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.tree.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

# fix round 1 move -> use if necessary
import glob
import os
#input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
#output_notwanted = '/scratch/users/anniz44/genomes/donor_species/unwanted/outlier'
#genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'

input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/test'
genome_root = '/scratch/users/mit_alm/IBD_Evo/Testing_Sets/'
genome_dir = glob.glob('/scratch/users/mit_alm/IBD_Evo/Testing_Sets/*')
output_notwanted = '/scratch/users/mit_alm/IBD_Evo/Testing_Sets/unwanted'

all_move = glob.glob(os.path.join(input_script,'*move*sh'))
all_move_scripts = []
for move_sh in all_move:
    move_sh_name = os.path.split(move_sh)[-1]
    if move_sh_name.startswith('SNP_round1'):
        for lines in open(move_sh,'r'):
            if lines.startswith('mv'):
                moved_files = os.path.split(lines.split(' ')[1])[-1]
                #donor_species = moved_files.split('_g')[0]
                donor_species = moved_files.split('_')[1] + '_' + moved_files.split('_')[0]
                genome_dir_sub = glob.glob('%s%s*'%(genome_root,moved_files.split('_')[1]))
                from_dir = lines.split(' ')[-1].split('\n')[0]
                all_move_scripts.append('mv %s %s//\n' %(os.path.join(from_dir,moved_files),genome_dir_sub))
            else:
                print(lines)

f1 = open(os.path.join(input_script, 'fix_round1_move.sh'), 'w')
f1.write('#!/bin/bash\n' + ''.join(all_move_scripts))
f1.close()

os.system('sh %s'%(os.path.join(input_script, 'fix_round1_move.sh')))
for move_sh in all_move:
    move_sh_name = os.path.split(move_sh)[-1]
    if move_sh_name.startswith('wholegenome.nonqualify'):
        os.system('sh %s' % (move_sh))

################################################### END ########################################################
################################################### SET PATH ########################################################

# alltree.sh finish
# round 1 SNP compare, SNP number and tree distance -> construct final cluster + remove unwanted
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/vcf_merge_new_finished'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round1'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge'
output_notwanted = '/scratch/users/anniz44/genomes/donor_species/unwanted/outlier'
cluster_cutoff = 4
fastq_name = '_1'
fastq_name_2 = '_2'

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species/vcf_merge_new_finished'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species'
genome_root = '/scratch/users/anniz44/genomes/donor_species/species/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/species/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/species/vcf_round1'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/species/vcf_merge'
output_notwanted = '/scratch/users/anniz44/genomes/donor_species/unwanted/outlier'
cluster_cutoff = 3
fastq_name = '_1'
fastq_name_2 = '_2'

#input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_merge_new'
#input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/test'
#genome_dir = glob.glob('/scratch/users/mit_alm/IBD_Evo/Testing_Sets/*')
#genome_root = '/scratch/users/mit_alm/IBD_Evo/Testing_Sets/'
#output_dir = '/scratch/users/anniz44/genomes/donor_species/test/vcf_round1'
#output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/test/vcf_merge'
#output_notwanted = '/scratch/users/mit_alm/IBD_Evo/Testing_Sets/unwanted'
#cluster_cutoff = 2
#fastq_name = '_1.'
#fastq_name_2 = '_2.'

# set up cutoff
SNP_total_cutoff_ratio = 0.2
SNP_total_cutoff_2 = 100

def SNP_seq(seq1, seq2, POS_info,POS_info_CHR,POS_info_CHR_LEN,POS_info_output,G1,G2):
    SNP_total = 0
    j = 0
    POS_DIS = []
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            if 'reference' not in G1 + G2:
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

Fasta_SNP = glob.glob(os.path.join(output_dir_merge + '/bwa/0','*.filtered.fasta'))
tree_distance_output = []
unwanted = []
second_strain = dict()
SNP_cutoff = []
SNP_pair = []
SNP_pair.append('G1\tG2\tSNP\t\n')

for fasta in Fasta_SNP:
    fasta_name = os.path.split(fasta)[-1]
    donor_species = fasta_name.split('.all.flt.snp.vcf.filtered.fasta')[0]
    gubbins_tree = glob.glob(os.path.join(output_dir_merge + '/bwa/0', fasta_name + '*.final_tree.tre'))
    if gubbins_tree == []:
        gubbins_tree = glob.glob(os.path.join(output_dir_merge + '/bwa/0',fasta_name + '*.align.nwk'))
    if gubbins_tree!= []:
        gubbins_tree = gubbins_tree[0]
        POS_file = \
        glob.glob(os.path.join(output_dir_merge + '/bwa/0', donor_species + '.all.flt.snp.vcf.filtered.POS'))[0]
        print(donor_species)
        Ref = ''
        distanceREF = 0
        SNP_total_length = 0
        tree_root_distance_list = dict()
        tree_root_distance_list.setdefault(0, [])
        tree_distance = dict()
        tree = Phylo.read(gubbins_tree, 'newick')
        Seq_list = dict()
        Cluster_SNP = dict()
        Cluster_SNP_set = dict()
        Cluster_SNP_set_added = set()
        max_cluster_diff = 0
        POS_info = []
        POS_info_CHR = []
        POS_info_CHR_LEN = dict()
        POS_info_output = []
        POS_info_output.append('G1\tG2\tCHR\tPOS\tDIS\tLEN\t\n')
        for lines in open(POS_file, 'r'):
            CHR = lines.split('\t')[0]
            POS_info.append(int(lines.split('\t')[1]))
            POS_info_CHR.append(CHR)
            POS_info_CHR_LEN.setdefault(CHR, 0)
            POS_info_CHR_LEN[CHR] = max(POS_info_CHR_LEN[CHR], int(lines.split('\t')[1]))
        for record in SeqIO.parse(fasta, 'fasta'):
            new_center = 1
            if Ref == '':
                #set upt the first seq as ref
                Ref = str(record.seq)
                REF_name = str(record.id)
                new_center = 0
                SNP_total = 0
                SNPREF = 0
                tree_distance.setdefault(str(record.id), SNP_total)
                SNP_total_length = len(Ref)
                SNP_total_cutoff = max(SNP_total_cutoff_ratio * SNP_total_length,SNP_total_cutoff_2)
            else:
                record_name = str(record.id)
                record_seq = str(record.seq)
                if SNPREF == 0 or distanceREF == 0:
                    distanceREF = tree.distance(record_name, REF_name)
                    SNP_total = SNP_seq(Ref, record_seq, POS_info, POS_info_CHR, POS_info_CHR_LEN,
                            POS_info_output, REF_name, record_name)
                    SNP_pair.append('%s\t%s\t%s\t\n'%(REF_name,record_name,SNP_total))
                    SNPREF = SNP_total
                    tree_distance.setdefault(record_name, SNP_total)
                    if SNP_total <= SNP_total_cutoff:
                        Cluster_SNP.setdefault(REF_name,[])
                        Cluster_SNP[REF_name].append(record_name)
                        max_cluster_diff = max(max_cluster_diff,SNP_total)
                else:
                    for record_before in Seq_list:
                        SNP_total = SNP_seq(Seq_list[record_before], record_seq, POS_info, POS_info_CHR, POS_info_CHR_LEN,
                                            POS_info_output, record_before, record_name)
                        #try:
                        #    SNP_total = SNP_distance_correct(tree.distance(record_name, REF_name), \
                        #                                      distanceREF, SNPREF)
                        #except ValueError:
                        #    # new seq not on the tree
                        #    print('%s not on the tree' % (str(record.id)))
                        #    SNP_total =  SNP_seq(Seq_list[record_before], record_seq, POS_info, POS_info_CHR, POS_info_CHR_LEN,
                        #                                             POS_info_output, record_before, record_name)
                        SNP_pair.append('%s\t%s\t%s\t\n' % (record_before, record_name, SNP_total))
                        if SNP_total <= SNP_total_cutoff:
                            Cluster_SNP.setdefault(record_before, [])
                            Cluster_SNP[record_before].append(record_name)
                            max_cluster_diff = max(max_cluster_diff, SNP_total)
                        if record_before == REF_name:
                            tree_distance.setdefault(record_name, SNP_total)
            Seq_list.setdefault(str(record.id), str(record.seq))
        cluster = 0
        SNP_cutoff.append('%s\t%s\t%s\t%s\t\n' % (donor_species,max_cluster_diff, SNP_total_cutoff, SNP_total_length))
        for record_name in Cluster_SNP:
            neighbor = Cluster_SNP.get(record_name,[])
            if neighbor != [] and record_name not in Cluster_SNP_set_added:
                cluster += 1
                Cluster_SNP_set.setdefault(cluster,set())
                Cluster_SNP_set[cluster].add(record_name)
                Cluster_SNP_set_added.add(record_name)
                find_neighbor(Cluster_SNP, neighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)
        SNP_perdistance = distanceREF/SNPREF
        for cluster in Cluster_SNP_set:
            tree_name_list = Cluster_SNP_set[cluster]
            if len(tree_name_list) >= cluster_cutoff:
                if cluster == 1:
                    for tree_name in tree_name_list:
                        donor_species_cluster = '%s' % (donor_species)
                        tree_distance_output.append(
                            '%s\t%s\t%s\tcluster%s\t%s\t\n' % (donor_species, tree_name, tree_distance[tree_name], cluster,SNP_perdistance))
                        second_strain.setdefault(donor_species_cluster, [])
                        second_strain[donor_species_cluster].append(tree_name)
                else:
                    donor_species_cluster = '%s_cluster%s' % (donor_species, cluster)
                    print(donor_species_cluster)
                    second_strain.setdefault(donor_species_cluster, [])
                    for tree_name in tree_name_list:
                        tree_distance_output.append(
                            '%s\t%s\t%s\tcluster%s\t%s\t\n' % (donor_species, tree_name, tree_distance[tree_name], cluster,SNP_perdistance))
                        second_strain[donor_species_cluster].append(tree_name)
            else:
                for tree_name in tree_name_list:
                    tree_distance_output.append(
                        '%s\t%s\t%s\tcluster_unwanted\t%s\t\n' % (donor_species, tree_name, tree_distance[tree_name],SNP_perdistance))
                    unwanted.append(tree_name)
        for record_name in Seq_list:
            if record_name not in Cluster_SNP_set_added:
                tree_distance_output.append(
                    '%s\t%s\t%s\tcluster_unwanted\t%s\t\n' % (
                    donor_species, record_name, tree_distance[record_name], SNP_perdistance))
                unwanted.append(record_name)
        f1 = open(POS_file + '.sum', 'w')
        f1.write('%s' % (''.join(POS_info_output)))
        f1.close()
    else:
        print('tree missing %s'%(donor_species))

os.system('mv %s %s'%(os.path.join(input_script, 'SNP_round1.allcompare.sum'),
                      os.path.join(input_script, 'SNP_round1.allcompare.old.sum')))

f1 = open(os.path.join(input_script, 'SNP_round1.allcompare.sum'), 'w')
f1.write('%s' % (''.join(tree_distance_output)))
f1.close()

f1 = open(os.path.join(input_script, 'SNP_round1.allcompare.sum.cutoff'), 'w')
f1.write('%s' % (''.join(SNP_cutoff)))
f1.close()

f1 = open(os.path.join(input_script, 'SNP_round1.allcompare.allpair.sum'), 'w')
f1.write('%s' % (''.join(SNP_pair)))
f1.close()


#os.system('cat %s | cut -f 2,4 |sort > %s'%(os.path.join(input_script, 'SNP_round1.allcompare.sum'),os.path.join(input_script, 'SNP_round1.allcompare.sum.brief')))
#os.system('cat %s | cut -f 2,4 |sort > %s'%(os.path.join(input_script, 'SNP_round1.sum'),os.path.join(input_script, 'SNP_round1.sum.brief')))
#os.system('diff %s %s > %s'%(os.path.join(input_script, 'SNP_round1.sum.brief'),
#                             os.path.join(input_script, 'SNP_round1.allcompare.sum.brief'),
#                             os.path.join(input_script, 'SNP_round1.sum.brief.diff')))

cmd_move=''
for files in unwanted:
    cmd_move += ('mv %s %s/\n'%(os.path.join(genome_root + '/*',files + '.*'),output_notwanted))
    cmd_move += ('mv %s %s/\n' % (os.path.join(genome_root + '/*', files + fastq_name + '*'), output_notwanted))
    cmd_move += ('mv %s %s/\n' % (os.path.join(genome_root + '/*', files + fastq_name_2 + '*'), output_notwanted))
    cmd_move += ('mv %s %s/\n' % (os.path.join(genome_root + '/*/fastq', files + fastq_name + '*'), output_notwanted))
    cmd_move += ('mv %s %s/\n' % (os.path.join(genome_root + '/*/fastq', files + fastq_name_2 + '*'), output_notwanted))

for donor_species in second_strain:
    folders_dir = os.path.join(genome_root,donor_species)
    try:
        os.mkdir(folders_dir)
    except IOError:
        pass
    for genome_files in second_strain[donor_species]:
        cmd_move+=('mv %s %s/\n'%(os.path.join(genome_root + '/*',genome_files + '.*'),folders_dir))
        cmd_move += ('mv %s %s/\n' % (os.path.join(genome_root + '/*', genome_files + fastq_name + '*'), folders_dir))
        cmd_move += ('mv %s %s/\n' % (os.path.join(genome_root + '/*', genome_files + fastq_name_2 + '*'), folders_dir))
        cmd_move += ('mv %s %s/\n' % (os.path.join(genome_root + '/*/fastq', genome_files + fastq_name + '*'), folders_dir))
        cmd_move += ('mv %s %s/\n' % (os.path.join(genome_root + '/*/fastq', genome_files + fastq_name_2 + '*'), folders_dir))
        cmd_move += ('mv %s %s/\n' % (os.path.join(output_notwanted, genome_files + '.*'), folders_dir))
        cmd_move += ('mv %s %s/\n' % (os.path.join(output_notwanted, genome_files + fastq_name + '*'), folders_dir))
        cmd_move += ('mv %s %s/\n' % (os.path.join(output_notwanted, genome_files + fastq_name_2 +'*'), folders_dir))

f1 = open(os.path.join(input_script, 'SNP_round1.move.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmd_move)))
f1.close()

# move first
# round 2
# not used: deleting genomes >= 6MB, apparent contamination
# first individual genome assembly
# am.E.coli and am.T needs to run

################################################### END ########################################################
################################################### SET PATH ########################################################

import glob
import os
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species/spades'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/species/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/species/vcf_round1'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/species/vcf_merge'
fastq_name = '_1.fastq'
output_notwanted = '/scratch/users/anniz44/genomes/donor_species/unwanted/outlier'

#input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/spades'
#input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/test'
#genome_dir = glob.glob('/scratch/users/mit_alm/IBD_Evo/Testing_Sets/*')
#genome_root = '/scratch/users/mit_alm/IBD_Evo/Testing_Sets/'
#output_dir = '/scratch/users/anniz44/genomes/donor_species/test/vcf_round1'
#output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/test/vcf_merge'
#output_notwanted = '/scratch/users/mit_alm/IBD_Evo/Testing_Sets/unwanted'
#fastq_name = '_1.fastq'
os.system('rm -rf %s'%(input_script_sub))
try:
    os.mkdir(input_script_sub)
except IOError:
    pass

for folder in genome_dir:
    if folder != output_notwanted:
        print(folder)
        fastq = glob.glob(os.path.join(folder,'*' + fastq_name)) + glob.glob(os.path.join(folder,'fastq/*' + fastq_name))
        print(fastq)
        if fastq != []:
            donor_species = os.path.split(folder)[-1]
            cmd = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
            donor_species_folder_all = os.path.join(folder,donor_species + '_allspades')
            # select 0.25m reads
            for fastq_file in fastq:
                os.system('rm -rf %s/%s.all*' %(folder, donor_species))
                os.system('rm -rf %s' % (donor_species_folder_all))
                if '.all' + fastq_name not in fastq_file:
                    donor_species_genome = fastq_file.split(fastq_name)[0]
                    #donor_species_genome_name = os.path.split(donor_species_genome)[-1]
                    #donor_species_folder = os.path.join(folder,donor_species_genome_name + '_spades')
                    #cmd += 'spades.py --careful -1 %s%s -2 %s%s -o %s --threads 40 --memory 500 --cov-cutoff 10\n' %(donor_species_genome,fastq_name,donor_species_genome,fastq_name,donor_species_folder)
                    #cmd += 'mv %s/scaffolds.fasta %s/%s.spades.fasta\n' %(donor_species_folder,folder,donor_species_genome_name)
                    cmd += 'head -500000 %s%s >> %s/%s.all%s\n' % (
                    donor_species_genome,fastq_name, folder, donor_species,fastq_name)
                    cmd += 'tail -500000 %s%s >> %s/%s.all%s\n' % (
                    donor_species_genome, fastq_name,folder, donor_species,fastq_name)
                    cmd += 'head -500000 %s%s >> %s/%s.all%s\n' % (
                        donor_species_genome, fastq_name.replace('1','2'), folder, donor_species, fastq_name.replace('1','2'))
                    cmd += 'tail -500000 %s%s >> %s/%s.all%s\n' % (
                        donor_species_genome, fastq_name.replace('1','2'), folder, donor_species, fastq_name.replace('1','2'))
            cmd += 'spades.py --careful -1 %s/%s.all%s -2 %s/%s.all%s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' %\
                   (folder,donor_species,fastq_name,folder,donor_species,fastq_name.replace('1','2'), donor_species_folder_all)
            cmd += 'mv %s/scaffolds.fasta %s/%s.all.spades.fasta\n' %(donor_species_folder_all,folder,donor_species)
            f1 = open(os.path.join(input_script_sub, '%s.sh' % donor_species), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmd)))
            f1.close()


f1 = open(os.path.join(input_script, 'allspades.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    #try:
    #    ferr = open(sub_scripts+'.err2','r')
    #except IOError:
    #    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################

# round 2 after spades annotate pangenome
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

input_script='/scratch/users/anniz44/scripts/1MG/donor_species'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/prokka'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*') + \
glob.glob('/scratch/users/anniz44/genomes/donor_species/species/*')
genome_root1 = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
genome_root2 = '/scratch/users/anniz44/genomes/donor_species/species/'
fasta_name = '.all.spades.fasta'

input_script='/scratch/users/anniz44/scripts/1MG/donor_species/test/'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/prokka'
genome_dir = glob.glob('/scratch/users/mit_alm/IBD_Evo/Testing_Sets/*')
genome_root = '/scratch/users/mit_alm/IBD_Evo/Testing_Sets/'
fasta_name = '.all.spades.fasta'

os.system('rm -rf %s'%(input_script_sub))
try:
    os.mkdir(input_script_sub)
except IOError:
    pass

i=0
#fasta_all = glob.glob(os.path.join(genome_root1,'*/*%s'%(fasta_name)))+\
#glob.glob(os.path.join(genome_root2,'*/*%s'%(fasta_name)))
fasta_all = glob.glob(os.path.join(genome_root,'*/*%s'%(fasta_name)))

for fasta_all_file in fasta_all:
    fasta_folder, fasta_file = os.path.split(fasta_all_file)
    donor_species = fasta_file.split(fasta_name)[0]
    cmdsprokka = ''
    cmdsprodigal = 'prodigal -a %s -d %s -q -i %s\n'%(
        fasta_all_file.replace(fasta_name,fasta_name.replace('.fasta','.faa')),
        fasta_all_file.replace(fasta_name, fasta_name.replace('.fasta', '.fna')),
        fasta_all_file)
    cmdsprokka += 'py37\nprokka --kingdom Bacteria --outdir %s/prokka_%s --genus Bacteroides --locustag Bacter %s\n' %\
                  (fasta_folder, donor_species,fasta_all_file)
    cmdsprokka += 'mv %s/prokka_%s/*.gff %s/%s.all.gff\n' %(fasta_folder,donor_species, fasta_folder,donor_species)
    cmdsprokka += 'mv %s/prokka_%s/*.tsv %s/%s.all.tsv\n' % (fasta_folder, donor_species, fasta_folder, donor_species)
    cmdsprokka += 'mv %s/prokka_%s/*.faa %s/%s.all.faa\n' % (fasta_folder, donor_species, fasta_folder, donor_species)
    cmdsprokka += 'rm -rf %s/prokka_%s\n' % (fasta_folder, donor_species)
    f1 = open(os.path.join(input_script_sub, '%s.%s.prokka.sh' % (donor_species,i)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmdsprokka)))
    f1.close()
    f1 = open(os.path.join(input_script_sub, '%s.%s.prodigal.sh' % (donor_species, i)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmdsprodigal)))
    f1.close()
    i += 1

f1 = open(os.path.join(input_script, 'allprokkanew.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.prokka.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'allprodigal.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.prodigal.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################

# round 2 run vcf -> need run allnewmergevcf.round2.sh
import glob
import os
input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/vcf_new_merge_round2'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/vcf_new_round2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round2'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2'
fastq_name = '_1.fastq'

#input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_new_round2'
#input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/test'
#genome_dir = glob.glob('/scratch/users/mit_alm/IBD_Evo/Testing_Sets/*')
#output_dir = '/scratch/users/anniz44/genomes/donor_species/test/vcf_round2'
#output_merge = '/scratch/users/anniz44/genomes/donor_species/test/vcf_merge_round2'
#input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_new_merge_round2'
#output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/test/vcf_merge_round2'
#fastq_name = '_1.fastq'

os.system('#rm -rf %s %s' %(input_script_merge_sub, input_script_sub))
os.system('#rm -rf %s %s' %(output_dir_merge, output_dir))

try:
    os.mkdir(input_script_merge_sub)
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
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(output_dir_merge)
except IOError:
    pass

try:
    os.mkdir(output_dir_merge+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir_merge+'/bwa/0')
except IOError:
    pass

# generate codes to run alignment and SNP calling
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    if len(donor) == 2: # BN10 donors
        database_files = glob.glob(os.path.join(folder,'%s.all.spades.fasta'%(donor_species)))
        if database_files != []:
            i = 0
            database = database_files[0]
            print(donor_species,database,database_files)
            try:
                f1 = open(database + '.1.bt2','r')
            except IOError:
                os.system('bowtie2-build %s %s\n' %(database,database))
            cmds = ''
            cmds_merge = ''
            sub_samples = []
            donor_species_fastq = glob.glob(os.path.join(folder,'*'+ fastq_name) ) + glob.glob(os.path.join(folder,'fastq/*'+ fastq_name) )
            for files in donor_species_fastq:
                if 'all' + fastq_name not in files:
                    # not including the mixed sample
                    donor_species_dir_file = os.path.split(files)[-1]
                    print(donor_species_dir_file)
                    tempbamoutput = os.path.join(output_dir +'/bwa/0', donor_species_dir_file)
                    try:
                        f1 = open(tempbamoutput + '.sorted.bam')
                    except IOError:
                        sub_samples.append(tempbamoutput + '.sorted.bam')
                        files2 = files.replace(fastq_name,fastq_name.replace('1','2'))
                        #files2 = files.replace('_1.fastq', '_2.fastq')
                        cmds += 'bowtie2' + ' -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --threads %s --un-conc %s.unaligned -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                            min(40, 40),tempbamoutput, database, files, files2,'samtools', min(40, 40),
                            tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
                            tempbamoutput)
                        cmds += '#samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
                            tempbamoutput, tempbamoutput)
                        cmds += 'rm -rf %s.bam\n' % (tempbamoutput)
                        cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
                        cmds += 'rm -rf %s.*.unaligned\n' % (tempbamoutput)
                        i += 1
                        if i%10 == 0:
                            f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species,int(i/10))), 'a')
                            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                            f1.close()
                            cmds = ''
            # mileup
            tempbamoutput = os.path.join(output_dir_merge + '/bwa/0', donor_species + '.all')
            try:
                f1 = open(tempbamoutput + '.sorted.bam')
            except IOError:
                if sub_samples != []:
                    cmds_merge += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -B -Ou -d3000 -f %s %s  | %s call --ploidy 1 --threads %s -m > %s.raw.vcf\n' % (
                        'bcftools', min(40, 40), database,
                        ' '.join(sub_samples), 'bcftools', min(40, 40), tempbamoutput)
                    cmds_merge += '%s filter --threads %s -s  LowQual %s.raw.vcf > %s.flt.vcf \n' % (
                        'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
                    cmds_merge += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
                        'bcftools', tempbamoutput, tempbamoutput)
                    cmds_merge += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
                    f2 = open(os.path.join(input_script_merge_sub, '%s.vcf.sh' % donor_species), 'w')
                    f2.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds_merge)))
                    f2.close()
            f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species,int(i/10))), 'a')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            f1.close()
        else:
            print('database missing',os.path.join(folder,'%s.all.spades.fasta'%(donor_species)))

f1 = open(os.path.join(input_script, 'allnewvcf.round2.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'allnewmergevcf.round2.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_merge_sub, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()


################################################### END ########################################################
################################################### SET PATH ########################################################
# need finish run allnewmergevcf.round2.sh
# round 2 SNP filtering setup parameter by matlab
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics

# set up path -> test
input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_new_round2'
input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_new_merge_round2'
#input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_new'
#input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_merge_new'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/tree_round2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/test'
genome_dir = glob.glob('/scratch/users/mit_alm/IBD_Evo/Testing_Sets/*')
genome_root = '/scratch/users/mit_alm/IBD_Evo/Testing_Sets/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/test/vcf_round2'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/test/vcf_merge_round2'
#output_dir = '/scratch/users/anniz44/genomes/donor_species/test/vcf_round1'
#output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/test/vcf_merge'
fastq_name = '_1.fastq'
deleting_file = ['H30_08','H30_16']
Cov_dis = 30
#Cov_dis = 50
rule1=['BA_H25_cluster4',
       'BA_H25_cluster3',
       'BA_H25','PD_H30']
rule2=['PD_H30']
rule3 = ['BA_H25_cluster4',
         'BA_H25']
# step1 tune cutoff -> matlab
#Step = 1
# step2 after tune cutoff cluster potential recombinations -> matlab check
Step = 2
# step3 final vcf filtering
#Step = 3

# set up path -> selected
input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/vcf_new_round2'
input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/vcf_new_merge_round2'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/tree_round2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round2'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2'
output_notwanted = '/scratch/users/anniz44/genomes/donor_species/unwanted/outlier'
fastq_name = '_1.fastq'
genome_split = '_g'
deleting_file = []
Cov_dis = 50
# other rules
# rule 1 SNP_depth_cutoff set as 3 and Major_alt_freq_cutoff as 0.95 for the SNP
rule1 = ['am_Turicibacter_sanguinis','am_Turicibacter_sanguinis_newcluster2','an_Parabacteroides_distasonis',
         'ao_Bacteroides_ovatus','av_Bifidobacterium_longum','av_Lactobacillus_ruminis',
         'cx_Bacteroides_ovatus','cx_Bifidobacterium_longum','cx_Escherichia_coli']
# rule 2 SNP_depth_cutoff set as 50 and Major_alt_freq_cutoff as 0.75
rule2 = ['cx_Bifidobacterium_longum','af_Turicibacter_sanguinis']
# rule 3 SNP_depth_cutoff set as 15 and Major_alt_freq_cutoff as 0.8
rule3 = []
# step1 tune cutoff -> matlab
#Step = 1
# step2 after tune cutoff cluster potential recombinations -> matlab check
#Step = 2
# step3 final vcf filtering
Step = 3

# set up path -> all_species
input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species/vcf_new_round2'
input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species/vcf_new_merge_round2'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species/tree_round2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/all_species'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/species/*')
genome_root = '/scratch/users/anniz44/genomes/donor_species/species/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/species/vcf_round2'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/species/vcf_merge_round2'
output_notwanted = '/scratch/users/anniz44/genomes/donor_species/unwanted/outlier'
fastq_name = '_1.fastq'
genome_split = '_g'
deleting_file = []
Cov_dis = 50
# other rules
# rule 1 SNP_depth_cutoff set as 3 and Major_alt_freq_cutoff as 0.95 for the SNP
rule1=['aa_Akkermansia_sp',
       'aa_Bifidobacterium_adolescentis_cluster2',
       'aa_Eubacterium_rectale',
       'aa_Faecalibacterium_prausnitzii_2_cluster2',
       'af_Bacteroides_fragilis',
       'af_Bacteroides_thetaiotaomicron',
       'af_Bacteroides_thetaiotaomicron_cluster2',
       'af_Blautia_wexlerae',
       'af_Blautia_wexlerae_cluster2',
       'af_Parabacteroides_merdae_cluster2',
       'am_Bacteroides_fragilis',
       'am_Bacteroides_fragilis_cluster2',
       'am_Bacteroides_salyersiae',
       'am_Bacteroides_vulgatus',
       'am_Bacteroides_vulgatus_cluster2',
       'am_Bifidobacterium_adolescentis',
       'am_Bifidobacterium_adolescentis_cluster2',
       'am_Bifidobacterium_pseudocatenulatum',
       'am_Collinsella_aerofaciens',
       'am_Parabacteroides_merdae',
       'am_Parasutterella_excrementihominis',
       'am_Parasutterella_excrementihominis_cluster2',
       'am_Ruthenibacterium_lactatiformans',
       'an_Bacteroides_xylanisolvens',
       'ao_Bifidobacterium_adolescentis',
       'ao_Bifidobacterium_adolescentis_cluster2',
       'ao_Bifidobacterium_pseudocatenulatum',
       'ao_Bifidobacterium_pseudocatenulatum_cluster2',
       'av_Bacteroides_vulgatus',
       'av_Bifidobacterium_bifidum',
       'av_Clostridium_beijerinckii',
       'bk_Bifidobacterium_adolescentis',
       'bq_Parabacteroides_distasonis',
       'cx_Bacteroides_thetaiotaomicron']
# rule 2 SNP_depth_cutoff set as 50 and Major_alt_freq_cutoff as 0.75
rule2=['aa_Bifidobacterium_adolescentis',
       'aa_Bifidobacterium_adolescentis_cluster2',
       'aa_Eubacterium_rectale',
       'am_Bifidobacterium_adolescentis',
       'am_Bifidobacterium_pseudocatenulatum',
       'am_Collinsella_aerofaciens',
       'am_Parabacteroides_merdae',
       'an_Enterococcus_hirae',
       'ao_Bifidobacterium_pseudocatenulatum',
       'bk_Bifidobacterium_adolescentis',
       'cx_Streptococcus_parasanguinis_cluster2']
# rule 3 SNP_depth_cutoff set as 15 and Major_alt_freq_cutoff as 0.8
rule3 = []
# step1 tune cutoff -> matlab
#Step = 1
# step2 after tune cutoff cluster potential recombinations -> matlab check
#Step = 2
# step3 final vcf filtering
Step = 3


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

def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = int(Subdepth[0]) + int(Subdepth[1])
    Subdepth_ALT = int(Subdepth[2]) + int(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return allels_set[1]
    else:
        return allels_set[0]

def outputvcf(output_name):
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
            SNP_alignment_output_parsi.append('S%s    %s\n' % (genomename[-8:], SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.%s.vcf' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list_vcf))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.fasta' % (output_name), 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.parsi.fasta' %(output_name), 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()
    if Step == 2:
        if SNP_alignment_output != []:
            # with SNP
            fasta = vcf_file + '.%s.fasta' % (output_name)
            SNP_tree_cmd.append('python /scratch/users/anniz44/scripts/maffttree/remove.duplicated.py -i %s\n' % (fasta))
            filesize = int(os.path.getsize(fasta))
            if filesize >= 10 ^ 7:  # 10Mb
                SNP_tree_cmd.append(
                    'mafft --nuc --quiet --nofft --retree 2 --maxiterate 0 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                        # 'mafft --nuc --quiet --retree 2 --maxiterate 100 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                        fasta, fasta))  # remove --adjustdirection
                SNP_tree_cmd2.append(
                    '#raxmlHPC -s %s -m GTRGAMMA -n %s.raxml -p 1000 -T 40\n' % (
                        fasta, fasta))
            else:
                SNP_tree_cmd.append(
                    # 'mafft --nuc --quiet --nofft --retree 2 --maxiterate 0 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                    'mafft --nuc --quiet --retree 2 --maxiterate 100 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                        fasta, fasta))  # remove --adjustdirection
                SNP_tree_cmd2.append(
                    '#raxmlHPC -s %s -m GTRGAMMA -n %s.raxml -p 10000 -T 40\n' % (
                        fasta, fasta))
            SNP_tree_cmd.append('FastTree -nt -quiet %s.mafft.align > %s.mafft.align.nwk\n' % (fasta, fasta))
            SNP_tree_cmd2.append(
                '#run_gubbins.py --threads 40  --tree_builder hybrid --use_time_stamp --prefix %s_gubbins --verbose %s.mafft.align\n' % (
                    fasta, fasta))
            SNP_tree_cmd2.append(
                'mv *%s.mafft* %s\n' % (
                    fasta,output_dir_merge + '/bwa/0/tree'))
            SNP_tree_cmd2.append('rm -rf %s.sorted*\n'%(fasta))
        f1 = open(os.path.join(input_script_sub, '%s.tree.all.sh' % (donor_species)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n' + ''.join(SNP_tree_cmd) + ''.join(SNP_tree_cmd2))
        f1.close()

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
        REF = curate_REF(allels_set, Depth4)  # as the major alt in the population
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
                                if total_sub_depth_forward - int(Subdepth_forward[0]) >= 50 and \
                                total_sub_depth_reverse - int(Subdepth_reverse[0]) >= 50:
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
                                if total_sub_depth_forward - int(Subdepth_forward[0]) >= 15 and \
                                total_sub_depth_reverse - int(Subdepth_reverse[0]) >= 15:
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
                                    if total_sub_depth_forward - int(Subdepth_forward[0]) >= 3 and \
                                            total_sub_depth_reverse - int(Subdepth_reverse[0]) >= 3 and \
                                            Major_ALT[1] >= 6 and MLF >= 0.95:
                                        qualify_SNP = 1
                                        temp_snp_line_pass += 'RULE'
                                # no specific rule
                                if total_sub_depth_forward - int(Subdepth_forward[0]) >= SNP_depth_cutoff and \
                                        total_sub_depth_reverse - int(Subdepth_reverse[0]) >= SNP_depth_cutoff and \
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
                                            if total_sub_depth_forward - int(Subdepth_forward[0]) >= SNP_depth_cutoff2 and \
                                                    total_sub_depth_reverse - int(Subdepth_reverse[0]) >= SNP_depth_cutoff2 and \
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
        tempbamoutput = os.path.join(output_dir + '/bwa/0', filename + fastq_name)
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
################################################### Set up ########################################################
# set up output

try:
    os.mkdir(output_dir_merge + '/bwa/0/tree')
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
total_coverage_cutoff = 0.6
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
all_vcf_file=glob.glob(os.path.join(output_dir_merge + '/bwa/0','*.all.flt.snp.vcf'))
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
        for lines in open(os.path.join(input_script_merge_sub, '%s.vcf.sh' % (donor_species)), 'r'):
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
                # step 1 and 3 output coverage, step 2 use step 1 coverage
                try:
                    f1 = open(vcf_file + '.%s.cov.txt' % (output_name), 'r')
                except IOError:
                    outputcov(output_name,list(vcf_file_POS_candidate),Sample_subset)

# run parsi tree
if Step >1:
    parsi_files = glob.glob(os.path.join(output_dir_merge + '/bwa/0', '*parsi.fasta'))
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
    os.system('mv %s/bwa/0/*.parsi* %s/bwa/0/tree' % (
        output_dir_merge, output_dir_merge))

# run clustering and mafft tree
if Step == 2:
    # sum up SNP
    SNP_cutoff = []
    SNP_pair = []
    SNP_pair.append('G1\tG2\tSNP\t\n')
    POS_info_output = []
    POS_info_output.append('G1\tG2\tCHR\tPOS\tDIS\tLEN\t\n')
    Fasta_SNP = glob.glob(os.path.join(output_dir_merge + '/bwa/0', '*.filtered.fasta'))
    tree_distance_output = []
    for fasta in Fasta_SNP:
        fasta_name = os.path.split(fasta)[-1]
        donor_species = fasta_name.split('.all.flt.snp.vcf.filtered.fasta')[0]
        POS_file = glob.glob(os.path.join(output_dir_merge + '/bwa/0', donor_species + '.all.flt.snp.vcf.filtered.POS.txt'))[0]
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
    f1 = open(os.path.join(input_script, 'SNP_round2.sum'), 'w')
    f1.write('donor_species\tGenome\tSNP_total\tcluster1\tsubcluster\t\n%s' % (''.join(tree_distance_output)))
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round2.sum.cutoff'), 'w')
    f1.write('donor_species\tmax_cluster_diff\tSNP_cutoff\tSNP_total_len\t\n%s' % (''.join(SNP_cutoff)))
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round2.POS.sum'), 'w')
    f1.write('%s' % (''.join(POS_info_output)))
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round2.allpair.sum'), 'w')
    f1.write('Genome1\tGenome2\tSNP_total\t\n%s' % (''.join(SNP_pair)))
    f1.close()
    # run mafft tree
    f1 = open(os.path.join(input_script, 'allround2alltree.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.tree.all.sh')):
        #f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
        f1.write('sh %s %s.out\n' % (sub_scripts, sub_scripts))
    f1.close()

# run coverage checking
if Step == 1:
    Donor_species = dict()
    Coverage = []
    Coverage.append('genome\treads_alignment_rate\tref_genome_length\taverage_coverage\tquality\t\n')
    sh_files = glob.glob(os.path.join(input_script_split_sub,'*.sh'))
    for sh_file in sh_files:
        coverage_check(sh_file, Coverage)
    f1 = open(os.path.join(input_script, 'SNP_round2.coverage.sum'), 'w')
    f1.write(''.join(Coverage))
    f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# after round 2 step 3 calculate NS ratio
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics

# set up path -> test
input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_new_round2'
input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/test/vcf_new_merge_round2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/test'
genome_dir = glob.glob('/scratch/users/mit_alm/IBD_Evo/Testing_Sets/*')
genome_root = '/scratch/users/mit_alm/IBD_Evo/Testing_Sets/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/test/vcf_round2'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/test/vcf_merge_round2'
output_dir_merge2 = ''
fastq_name = '_1.fastq'

# set up path -> selected
input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/vcf_new_round2'
input_script_merge_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/vcf_new_merge_round2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')+\
glob.glob('/scratch/users/anniz44/genomes/donor_species/species/*')
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round2'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2'
output_dir_merge2 = '/scratch/users/anniz44/genomes/donor_species/species/vcf_merge_round2'

fastq_name = '_1.fastq'

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
NSratioobserve_cutoff = 0.5
Min_SNP_highselect_cutoff = 1/1500
Max_SNP_highselect_cutoff = 0.02
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
    #if (REF in purines and ALT in purines):
    #    return [0,[REF,ALT]]
    #elif (REF in pyrimidines and ALT in pyrimidines):
    #    return [0,[REF,ALT]]
    #return [1,[REF,ALT]]

def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = int(Subdepth[0]) + int(Subdepth[1])
    Subdepth_ALT = int(Subdepth[2]) + int(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return allels_set[1]
    else:
        return allels_set[0]

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
    for Ref in Ref_NSratio:
        foutput_list.append('%s\t%s\t\n' % (Ref, Ref_NSratio[Ref]))
    foutput.write(''.join(foutput_list))
    foutput.close()
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
            and total_SNP_position >= 2 and \
            total_SNP_position / Gene_length >= Min_SNP_highselect_cutoff \
            and total_SNP_position / Gene_length <= Max_SNP_highselect_cutoff and\
        SNP_gene_temp.NSratio[0] >= SNP_gene_temp.NSratio[1] * NSratioobserve_cutoff:
        High_select = True
    new_line += '\t%s\n'%(High_select)
    return [new_line,High_select]

def freq_call(vcf_file,Ref_seq, Ref_NSratio,SNP_gene_species,SNP_gene_all,SNP_gene_all_highselect,Output2,donor_species):
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
        SNP_type.setdefault(Chr, [''] * Total)
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
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        REF = curate_REF(allels_set, Depth4)
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
            ALT_num = 1
            for ALT in ALT_set:
                ALT_frq = SNP_count_genome_count[0][ALT_num]
                SNP_pair = transitions(REF, ALT)
                # add to P
                SNP_gene_temp.addSNP_pair(SNP_pair, 2,
                                          ALT_frq, 1, Depth)
                SNP_gene_temp.addmutposition(Chr, position)
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
            total_genome_set = len(set(SNP_type[Chr])) - 1
            sumgene_line,High_select = sumgene(SNP_gene_temp,total_genome_set,donor_species,SNP_gene_species,Total,1)
            Output2.append(sumgene_line)
            if High_select:
                SNP_gene_all_highselect.addpredictSNP_pair(SNP_gene_temp.SNP_pair)
                for pair in SNP_gene_temp.SNP_pair_freq:
                    # use selected genes frequency * all genes NS ratio (codon NS sum of all genes)
                    SNP_gene_all_highselect.SNP_pair[pair][0] += SNP_gene_temp.SNP_pair[pair][0]
                    SNP_gene_all_highselect.SNP_pair[pair][1] += SNP_gene_temp.SNP_pair[pair][1]
                    SNP_gene_all_highselect.SNP_pair_freq[pair] += SNP_gene_temp.SNP_pair_freq[pair]
                SNP_gene_all_highselect.NSratio[0] += SNP_gene_temp.NSratio[0]
                SNP_gene_all_highselect.NSratio[1] += SNP_gene_temp.NSratio[1]
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

# set up sum all species
SNP_gene_all = SNP_gene() # all denovo mutation
SNP_gene_all.init('allspecies')
SNP_gene_all_highselect = SNP_gene() # all mutations of highly selected genes
SNP_gene_all_highselect.init('allspecies_highselect')

# process each vcf file
Output2 = []
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    if len(donor) == 2:
        print(donor_species)
        database_file = glob.glob(os.path.join(folder,
                                               '*.all.spades.fna'))
        vcf_file = glob.glob(os.path.join(output_dir_merge,'bwa/0/%s.all.flt.snp.vcf.filtered.vcf'%(donor_species)))+ \
                   glob.glob(os.path.join(output_dir_merge2, 'bwa/0/%s.all.flt.snp.vcf.filtered.vcf' % (donor_species)))
        vcf_file_cluster = glob.glob(os.path.join(output_dir_merge, 'bwa/0/%s.all.flt.snp.vcf.filtered.cluster*.vcf' % (donor_species)))+ \
                           glob.glob(os.path.join(output_dir_merge2,
                                                  'bwa/0/%s.all.flt.snp.vcf.filtered.cluster*.vcf' % (donor_species)))
        vcf_file2 = glob.glob(os.path.join(output_dir_merge, 'bwa/0/%s.all.flt.snp.vcf' % (donor_species)))+ \
                    glob.glob(os.path.join(output_dir_merge2, 'bwa/0/%s.all.flt.snp.vcf' % (donor_species)))
        if database_file!=[] and (vcf_file!=[] or vcf_file_cluster!=[]) and vcf_file2!=[]:
            if int(os.path.getsize(vcf_file[0])) > 0:
                # not empty snps
                print('running %s'%donor_species)
                Ref_seq, Ref_NSratio,Mapping,Mapping_loci,Reverse = loaddatabase(database_file[0])
                SNP_gene_species = SNP_gene()  # all mutations of highly selected genes
                SNP_gene_species.init(donor_species)
                if vcf_file_cluster!= []:
                    print('running subcluster for %s' % donor_species)
                    for vcf_file_sub_cluster in vcf_file_cluster:
                        Total = freq_call(vcf_file_sub_cluster, Ref_seq, Ref_NSratio, SNP_gene_species,SNP_gene_all, SNP_gene_all_highselect,Output2, donor_species)
                else:
                    Total = freq_call(vcf_file[0],Ref_seq, Ref_NSratio,SNP_gene_species,SNP_gene_all,SNP_gene_all_highselect,Output2,donor_species)
                if sum(SNP_gene_species.NSratio) > 0:
                    sumgene_line, High_select = sumgene(SNP_gene_species, 1, donor_species, SNP_gene_species, Total, 0)
                    Output2.append(sumgene_line)
                    print(SNP_gene_species.expectNSratio, SNP_gene_species.dNdS,SNP_gene_species.NSratio,SNP_gene_species.NSratiosum)
        else:
            print('missing files %s %s %s %s'%(donor_species,database_file,vcf_file,vcf_file2))

# sum all species dNdS
sumgene_line,High_select = sumgene(SNP_gene_all,1,'allspecies',SNP_gene_all,'None',0)
Output2.append(sumgene_line)
print(SNP_gene_all.expectNSratio,SNP_gene_all.dNdS,SNP_gene_all.NSratio,SNP_gene_all.NSratiosum)
sumgene_line,High_select = sumgene(SNP_gene_all_highselect,1,'allspecies',SNP_gene_all,'None',1)
Output2.append(sumgene_line)
print(SNP_gene_all_highselect.expectNSratio,SNP_gene_all_highselect.dNdS,SNP_gene_all_highselect.NSratio,SNP_gene_all_highselect.NSratiosum)

# output
foutput = open(output_dir_merge + '/summary/all.donor.species.sum.txt', 'w')
foutput.write('#donor_species\tgene\tNo.genome\tgene_length\ttotal_SNP_genomeset\tNo.SNP\tNo.SNP_position\tN\tS\tOther\tobserved_ratio\texpected_ratio\tdNdS\t' +\
            'A-T_freq\tA-T_N:S\tA-C_freq\tA-C_N:S\tG-C_freq\tG-C_N:S\tG-T_freq\tG-T_N:S\tA-G_freq\tA-G_N:S\tG-A_freq\tG-A_N:S\tHigh_selected\n')
foutput.write(''.join(Output2))
foutput.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# after calculate NS ratio, extract all genes
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics


# set up path -> selected
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/annotate'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')+\
glob.glob('/scratch/users/anniz44/genomes/donor_species/species/*')
genome_root1 = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
genome_root2 = '/scratch/users/anniz44/genomes/donor_species/species/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round2'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2'
fastq_name = '_1.fastq'

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# function
def sum_donor_species(input_summary,selected = 1):
    Donor_species = dict()
    for lines in open(input_summary,'r'):
        lines_set = lines.replace('\n','').split('\t')
        if selected == 0 or lines_set[-1] == 'True':
            # highly selected genes or all genes
            gene_name = lines_set[1]
            donor_species = lines_set[0]
            Donor_species.setdefault(donor_species,[])
            Donor_species[donor_species].append(gene_name)
    return Donor_species

def extract_donor_species(donor_species,input_fasta,gene_name_list,output_fasta):
    gene_name_extract = []
    donor_species_set = donor_species.split('_')
    record_id_set = set()
    donor_species = '%s_%s_%s'%(donor_species_set[0],
                                donor_species_set[1][0:min(6,len(donor_species_set[1]))],
                                donor_species_set[2][0:min(6,len(donor_species_set[2]))])
    if 'cluster' in donor_species_set[-1]:
        donor_species += '_CL' + donor_species_set[3].split('cluster')[1]
    gene_name_list = set(gene_name_list)
    for record in SeqIO.parse(input_fasta, 'fasta'):
        record_id = str(record.id)
        if record_id in gene_name_list and record_id not in record_id_set:
            gene_name_extract.append(record_id)
            record_id = 'C_%s_G_%s'%(record_id.split('_')[1], record_id.split('_')[-1])
            output_fasta.append('>%s__%s\n%s\n'%(donor_species, record_id, str(record.seq)))
            record_id_set.add(str(record.id))
    if len(gene_name_extract) < len(gene_name_list):
        print('missing genes in files %s extracted %s genes of %s genes'%(input_fasta,len(gene_name_extract),len(gene_name_list)))
        print('missing genes ' + ' '.join([gene_name for gene_name in gene_name_list if gene_name not in gene_name_extract]))
    return output_fasta

def annotation(all_filter_gene_fasta_file):
    # run prokka
    cmdsprokka = 'py37\nprokka --kingdom Bacteria --outdir %s/prokka_%s --genus Bacteroides --protein %s --locustag Bacter %s/%s\n' % \
                 (output_dir_merge + '/summary', os.path.split(all_filter_gene_fasta_file)[-1],
                  all_filter_gene_fasta_file,
                  output_dir_merge + '/summary',
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
    f1 = open(os.path.join(input_script, 'allannotate.sh'), 'w')
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

# extract sequences for all genes
input_summary = output_dir_merge + '/summary/all.donor.species.sum.txt'
output_gene = output_dir_merge + '/summary/all.selected.gene.all.faa'
output_gene_dna = output_dir_merge + '/summary/all.selected.gene.all.fna'
output_fasta = []
output_fasta_dna = []
Donor_species = sum_donor_species(input_summary,0)
for donor_species in Donor_species:
    print(donor_species)
    input_fasta = glob.glob(os.path.join(genome_root1,'%s/%s.all.spades.faa'%(donor_species,donor_species))) + \
                      glob.glob(os.path.join(genome_root2, '%s/%s.all.spades.faa' % (donor_species, donor_species)))
    input_fasta_dna = glob.glob(os.path.join(genome_root1, '%s/%s.all.spades.fna' % (donor_species, donor_species))) + \
                      glob.glob(os.path.join(genome_root2, '%s/%s.all.spades.fna' % (donor_species, donor_species)))
    gene_name_list = Donor_species[donor_species]
    if gene_name_list != []:
        if input_fasta != [] and input_fasta_dna != []:
            output_fasta = extract_donor_species(donor_species,input_fasta[0], gene_name_list, output_fasta)
            output_fasta_dna = extract_donor_species(donor_species,input_fasta_dna[0], gene_name_list, output_fasta_dna)
        else:
            print('missing files for %s' % (donor_species))
    else:
        print('missing genes for %s' % (donor_species))

f1 = open(output_gene, 'w')
f1.write(''.join(output_fasta))
f1.close()
f1 = open(output_gene_dna, 'w')
f1.write(''.join(output_fasta_dna))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# find highly selected genes in one species across populations and donors
import os
import glob
import copy
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
all_fasta = os.path.join(output_dir + '/summary', 'all.selected.gene.all.faa')
input_summary = output_dir + '/summary/all.donor.species.sum.txt'

# functions
# clustering
def cluster_uc(cluster_input):
    Clusters = dict()
    High_select2 = set()
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
                High_select2.add(record_name2)
                High_select2.add(record_name)
    if High_select2 != set():
        for record_name in High_select2:
            High_select2_output.append('%s\tTrue\t\n'%(record_name))
    return [Clusters,High_select2_output,High_select2]

# correct highly selected genes dNdS
def format_freq(freq):
    freq = freq.split(':')
    return [int(freq[0]),int(freq[1])]

def init_highselect(line_set):
    allspecies_highselect = line_set
    for i in [3,5,6,7,8,-13,-11,-9,-7,-5,-3]:
        allspecies_highselect[i] = 0
    for i in [-12,-10,-8,-6,-4,-2]:
        allspecies_highselect[i] = format_freq(line_set[i])
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
    # NS ratio
    allspecies_highselect[10] = allspecies_highselect[7]/allspecies_highselect[8]
    # expected NS ratio
    tempNS = [0,0]
    i = -12
    for freq in allspecies:
        tempNS[0] += freq * allspecies_highselect[i][0]
        tempNS[1] += freq * allspecies_highselect[i][1]
        i += 2
    allspecies_highselect[11] = tempNS[0] / tempNS[1]
    # dNdS
    allspecies_highselect[12] = allspecies_highselect[10] / allspecies_highselect[11]
    return allspecies_highselect


def add_new_selectiong(input_summary,High_select2):
    newoutput = []
    # use selected codon NS ratio (SNP pair) * all genes freq
    for lines in open(input_summary, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        Selected = line_set[-1]
        if Selected == 'False':
            donor_species = line_set[0]
            record_id = line_set[1]
            # all species
            if record_id == 'allspecies':
                # all genes freq
                newoutput.append(lines)
                allspecies = [int(line_set[-13]), int(line_set[-11]), int(line_set[-9]),
                              int(line_set[-7]), int(line_set[-5]), int(line_set[-3])]
            # all highly selected genes
            elif record_id == 'allspecies_highselect':
                # selected codon NS ratio (SNP pair)
                allspecies_highselect = init_highselect(line_set)
    for lines in open(input_summary, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        Selected = line_set[-1]
        if Selected == 'False':
            donor_species = line_set[0]
            record_id = line_set[1]
            # species or genes
            if record_id not in ['allspecies', 'allspecies_highselect']:
                donor_species_set = donor_species.split('_')
                donor_species = '%s_%s_%s' % (donor_species_set[0],
                                              donor_species_set[1][0:min(6, len(donor_species_set[1]))],
                                              donor_species_set[2][0:min(6, len(donor_species_set[2]))])
                record_id = '%s__C_%s_G_%s' % (donor_species, record_id.split('_')[1], record_id.split('_')[-1])
                # highly selected genes with NS ratio > 1
                if record_id in High_select2 and (line_set[10] == 'observe_N_only' or float(line_set[10]) >= 1):
                    newoutput.append('\t'.join(line_set[:-1]) + '\tTrue\n')
                    allspecies_highselect = add_highselect(line_set, allspecies_highselect)
                # other genes and species
                else:
                    newoutput.append(lines)
        # original highly selected
        elif Selected == 'True':
            newoutput.append(lines)
            allspecies_highselect = add_highselect(line_set, allspecies_highselect)
    # output new_allspecies_highselect
    allspecies_highselect = calculate_NS(allspecies, allspecies_highselect)
    print(allspecies_highselect, allspecies)
    for i in range(0, len(allspecies_highselect)):
        if i in [14, 16, 18, 20, 22, 24]:
            allspecies_highselect[i] = '%s:%s' % (allspecies_highselect[i][0],
                                                  allspecies_highselect[i][1])
        else:
            allspecies_highselect[i] = str(allspecies_highselect[i])
    print(allspecies_highselect, allspecies)
    newoutput.append('\t'.join(allspecies_highselect) + '\n')
    # output new summary
    f1 = open(input_summary + '.High_select2.txt', 'w')
    f1.write(''.join(newoutput))
    f1.close()


# clustering
Clusters_gene, High_select2_output,High_select2 = cluster_uc(all_fasta + '.uc')
f1 = open(all_fasta + '.High_select2.txt', 'w')
f1.write(''.join(High_select2_output))
f1.close()
# correcting
add_new_selectiong(input_summary,High_select2)


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
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/annotate'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')+\
glob.glob('/scratch/users/anniz44/genomes/donor_species/species/*')
genome_root1 = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
genome_root2 = '/scratch/users/anniz44/genomes/donor_species/species/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round2'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2'
fastq_name = '_1.fastq'
input_summary = output_dir_merge + '/summary/all.donor.species.sum.txt.High_select2.txt'
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
        if selected == 0 or lines_set[-1] == 'True':
            # highly selected genes or all genes
            gene_name = lines_set[1]
            donor_species = lines_set[0]
            Donor_species.setdefault(donor_species,[])
            Donor_species[donor_species].append(gene_name)
    return Donor_species

def extract_donor_species(donor_species,input_fasta,gene_name_list,output_fasta):
    gene_name_extract = []
    donor_species_set = donor_species.split('_')
    record_id_set = set()
    donor_species = '%s_%s_%s'%(donor_species_set[0],
                                donor_species_set[1][0:min(6,len(donor_species_set[1]))],
                                donor_species_set[2][0:min(6,len(donor_species_set[2]))])
    if 'cluster' in donor_species_set[-1]:
        donor_species += '_CL' + donor_species_set[3].split('cluster')[1]
    gene_name_list = set(gene_name_list)
    for record in SeqIO.parse(input_fasta, 'fasta'):
        record_id = str(record.id)
        if record_id in gene_name_list and record_id not in record_id_set:
            gene_name_extract.append(record_id)
            record_id = 'C_%s_G_%s'%(record_id.split('_')[1], record_id.split('_')[-1])
            output_fasta.append('>%s__%s\n%s\n'%(donor_species, record_id, str(record.seq)))
            record_id_set.add(str(record.id))
    if len(gene_name_extract) < len(gene_name_list):
        print('missing genes in files %s extracted %s genes of %s genes'%(input_fasta,len(gene_name_extract),len(gene_name_list)))
        print('missing genes ' + ' '.join([gene_name for gene_name in gene_name_list if gene_name not in gene_name_extract]))
    return output_fasta

def annotation(all_filter_gene_fasta_file):
    # run prokka
    cmdsprokka = 'py37\nprokka --kingdom Bacteria --outdir %s/prokka_%s --genus Bacteroides --protein %s --locustag Bacter %s/%s\n' % \
                 (output_dir_merge + '/summary', os.path.split(all_filter_gene_fasta_file)[-1],
                  all_filter_gene_fasta_file,
                  output_dir_merge + '/summary',
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
    f1 = open(os.path.join(input_script, 'allannotate.sh'), 'w')
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

# extract sequences
output_fasta = []
output_fasta_dna = []
Donor_species = sum_donor_species(input_summary)
for donor_species in Donor_species:
    print(donor_species)
    input_fasta = glob.glob(os.path.join(genome_root1,'%s/%s.all.spades.faa'%(donor_species,donor_species))) + \
                      glob.glob(os.path.join(genome_root2, '%s/%s.all.spades.faa' % (donor_species, donor_species)))
    input_fasta_dna = glob.glob(os.path.join(genome_root1, '%s/%s.all.spades.fna' % (donor_species, donor_species))) + \
                      glob.glob(os.path.join(genome_root2, '%s/%s.all.spades.fna' % (donor_species, donor_species)))
    gene_name_list = Donor_species[donor_species]
    if gene_name_list != []:
        if input_fasta != [] and input_fasta_dna != []:
            output_fasta = extract_donor_species(donor_species,input_fasta[0], gene_name_list, output_fasta)
            output_fasta_dna = extract_donor_species(donor_species,input_fasta_dna[0], gene_name_list, output_fasta_dna)
        else:
            print('missing files for %s' % (donor_species))
    else:
        print('missing genes for %s' % (donor_species))

f1 = open(output_gene, 'w')
f1.write(''.join(output_fasta))
f1.close()
f1 = open(output_gene_dna, 'w')
f1.write(''.join(output_fasta_dna))
f1.close()

# run clustering
annotation(output_gene)
#annotation_dna(output_gene_dna)

################################################### END ########################################################
################################################### SET PATH ########################################################
# gene summary for each species
import os
import glob
import copy
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
all_fasta = os.path.join(output_dir+'/summary', 'all.selected.gene.faa')
database_metacyc = '/scratch/users/mit_alm/database/metacyc/genes.col'
database_eggnog = '/scratch/users/mit_alm/database/eggnog/2_annotations.tsv'
database_kegg = '/scratch/users/anniz44/scripts/database/Kegg/ko_formated'

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
    Blast_output = dict()
    DB_name = dict()
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
    for database in [database_eggnog]:
        for lines in open(database, 'r'):
            if not lines.startswith('#'):
                lines_set = lines.split('\n')[0].split('\t')
                db_name = lines_set[1]
                function_name = ''
                annotation_fun = lines_set[3]
                if db_name in DB_name:
                    DB_name[db_name][0] = function_name
                    DB_name[db_name][1] = annotation_fun
    best_hit(Blast_output, 1)
    return [Blast_output,DB_name]

def annotate_kegg(blast_search):
    Blast_output = dict()
    DB_name = dict()
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
    best_hit(Blast_output, 1)
    return [Blast_output,DB_name]

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
        record_name2 = line_set[9]
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

def cluster_genes(All_annotation,Clusters):
    for gene_name in Clusters:
        cluster = Clusters.get(gene_name,'')
        donor = gene_name.split('_')[0]
        donor_species = gene_name.split('__')[0]
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
            new_species = species.split('_cluster')[0].split('_newcluster')[0]
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

all_fasta_folder, all_fasta_file = os.path.split(all_fasta)
database_prokka = glob.glob(os.path.join(all_fasta_folder, 'prokka_%s/*.tsv' % (all_fasta_file)))[0]
database_prokka_fasta = glob.glob(os.path.join(all_fasta_folder, 'prokka_%s/*.faa' % (all_fasta_file)))[0]
# set up gene annotation and clustering
Blast_output1, DB_name1 = annotate_kegg(all_fasta + '.cluster.aa.kegg.txt')
Blast_output2, DB_name2 = annotate_metacyc(all_fasta + '.cluster.aa.metacyc.txt')
Blast_output3, DB_name3 = annotate_eggnog(glob.glob(all_fasta + '.cluster.aa.eggnog.*.txt'))
Blast_output4, DB_name4 = annotate_prokka(all_fasta)
# sum up
Clusters_gene = cluster_uc(all_fasta + '.uc')
sum_annotation(Blast_output1, DB_name1, Blast_output2, DB_name2, Blast_output3, DB_name3, Blast_output4, DB_name4,
               Clusters_gene)


################################################### END ########################################################
################################################### SET PATH ########################################################

# phylogeny finished
input_sum = 'all.donor.species.sum.txt.High_select2.short.txt'
Edge_table = []
Edge_table.append('node1\tnode2\n') # as class
Node_table = []
Node_table.append('node\tphylo\tanno\n')

edges = set()
nodes = set()
phylo_set = [19,20,21,22,23,2]#phylum to species
nodes.add('Bacteria')
for lines in open(input_sum,'r'):
    if not lines.startswith('population'):
        if not lines.startswith('allspecies'):
            lines_set = lines.split('\n')[0].split('\t')
            phylo_anno = lines_set[20]
            for i in phylo_set:
                phylo = lines_set[i]
                phylo_short = phylo
                phylo_before = lines_set[i - 1]
                if i != 2:
                    phylo_before = lines_set[i - 1]
                else:
                    phylo_before = lines_set[23]
                    phylo_short = lines_set[5]
                nodes.add('%s\t%s\t%s\n' % (phylo, phylo_anno,phylo_short))
                if phylo == 'NA':
                    if phylo_before!='NA':
                        phylo= '%s_%s'%(phylo_before,i)
                if phylo_before == 'NA':
                    phylo_before = '%s_%s'%(lines_set[i - 2],i)
                phylo_edge = '%s\t%s\n' % (phylo_before, phylo)
                edges.add(phylo_edge)

for phylo_edge in edges:
    Edge_table.append(phylo_edge)

f1 = open(input_sum + '.edge.txt','w')
f1.write(''.join(Edge_table))
f1.close()

for phylo in nodes:
    Node_table.append(phylo)

f1 = open(input_sum + '.node.txt','w')
f1.write(''.join(Node_table))
f1.close()

# fix annotation
import os
import glob
import copy
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species'
all_fasta_cluster = os.path.join(output_dir+'/summary', 'all.selected.gene.faa.uc')
all_fasta_cluster_old = os.path.join(output_dir+'/summary_highlyselected1', 'all.selected.gene.faa.uc')

def load_cluster(cluster_file):
    Cluster = dict()
    for lines in open(cluster_file,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        cluster = lines_set[1]
        genename = lines_set[8]
        Cluster.setdefault(genename,cluster)
    return Cluster

def compare_cluster(oldcluster,newcluster,filename):
    clustermap = []
    clustermap.append('oldcluster\tnewcluster\t\n')
    for genename in oldcluster:
        if genename in newcluster:
            clustermap.append('%s\t%s\t\n'%(oldcluster[genename],newcluster[genename]))
    f1 = open(filename,'w')
    f1.write(''.join(clustermap))
    f1.close()

Clusterold = load_cluster(all_fasta_cluster_old)
Clusternew = load_cluster(all_fasta_cluster)
compare_cluster(Clusterold,Clusternew,all_fasta_cluster+'.mapping.txt')

################################################### END ########################################################
# clean up
os.system('#rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/*/prokka_*')
os.system('#rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge/bwa/0/*.all.raw.vcf')
os.system('rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge/bwa/0/*.sorted*')
os.system('mv /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge/bwa/0/*.mafft.align* /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge/bwa/0/tree')
os.system('mv /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge/bwa/0/*_gubbins* /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge/bwa/0/tree')
os.system('#rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2/bwa/0/*.all.raw.vcf')
os.system('rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2/bwa/0/*.sorted*')
os.system('mkdir /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2/bwa/0/tree')
os.system('mv /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2/bwa/0/*.mafft.align* /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2/bwa/0/tree')
os.system('mv /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2/bwa/0/*_gubbins* /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2/bwa/0/tree')
os.system('#rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/*/*_spades')

################################################### END ########################################################
################################################### SET PATH ########################################################

# metagenomes for highly selected genes
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
import statistics

genome_dir = '/scratch/users/anniz44/genomes/donor_species'
metagenome_summary = '/scratch/users/anniz44/Metagenomes/donor_species/summary_BN10'
#metagenome_summary = '/scratch/users/anniz44/Metagenomes/donor_species_test/summary'
metagenome_snp = '/scratch/users/anniz44/Metagenomes/donor_species/bwa/0'
output_dir_snps = '/scratch/users/anniz44/genomes/donor_species/filtering_SNPs'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/MG'

Wanted = []
Gene_info = dict()
Gene_pos = dict()
for lines in open(output_dir_snps + '/summary/all.gene.list.filtered.loci.txt'):
    gene_ID = lines.split('\t')[1]
    Wanted.append(gene_ID)
    # gene ID, POS, REF, ALT, N or S
    Gene_pos.setdefault(gene_ID,[])
    pos = lines.split('\t')[2]
    Gene_pos[gene_ID].append(pos)
    # donor_species, REF ALT N_or_S
    Gene_info.setdefault('%s_%s'%(gene_ID,pos),[lines.split('\t')[0],lines.split('\t')[3],lines.split('\t')[4],lines.split('\t')[5].split('\n')[0]])

Wanted.append ('all')

Cluster_fun = dict()
for lines in open(output_dir_snps + '/fasta/cluster.gene.function.select.txt'):
    gene_ID = lines.split('\t')[1]
    Cluster_fun.setdefault(gene_ID,lines)

# filter metagenomes results
all_snp_sum = glob.glob(os.path.join(metagenome_snp,'*.raw.vcf'))
# non_fixed_ALT, non_fixed_REF, non_fixed_ALT_other
def ALT_freq_mut_type(REF,ALT,Other,Total):
    MAX = max(REF, ALT, Other)
    if REF == MAX:
        return 'REF'
    elif ALT == MAX:
        return 'ALT'
    elif Other == MAX:
        return 'Other'

Allels_order = ['A','T','G','C']
def ALT_freq_frq(Allels_count, REF_genome, ALT_genome):
    REF_metagenome = [REF_genome,0]
    ALT_metagneome = [ALT_genome,0]
    Other_metagenome = ['',0]
    total_count = 0
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT = Allels_order[alleles]
        total_count += ALT_frq
        if ALT_frq > 1:
            # remove only 1 read mapping to the snps
            if ALT == REF_genome:
                REF_metagenome[1] = ALT_frq
            elif ALT == ALT_genome:
                ALT_metagneome[1] = ALT_frq
            else:
                Other_metagenome[0] += ALT
                Other_metagenome[1] += ALT_frq
    total_ALT_frq = ALT_metagneome[1] + Other_metagenome[1]
    if total_ALT_frq == 0:
        return '%s\t%s\t\t\t\t\tno_mutation\t' %(REF_genome,REF_metagenome[1]/total_count)
    else:
        mutation_type = ALT_freq_mut_type(REF_metagenome[1],ALT_metagneome[1],Other_metagenome[1],total_count)
        return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (REF_genome,REF_metagenome[1]/total_count,
                                             ALT_genome,ALT_metagneome[1]/total_count,
                                             Other_metagenome[0],Other_metagenome[1]/total_count,mutation_type)

def ALT_freq(REF_meta,ALT_meta,Subdepth_set,REF_genome, ALT_genome):
    REF_metagenome = [REF_genome,0]
    ALT_metagneome = [ALT_genome,0]
    Other_metagenome = ['',0]
    total_count = 0
    i = 0
    for alleles in [REF_meta] + ALT_meta:
        if alleles != '.':
            depth_allele = int(Subdepth_set[i])
            if depth_allele > 1:
                # remove only 1 read mapping to the snps
                total_count += depth_allele
                if alleles == REF_genome:
                    REF_metagenome[1] = depth_allele
                elif alleles == ALT_genome:
                    ALT_metagneome[1] = depth_allele
                else:
                    Other_metagenome[0] += alleles
                    Other_metagenome[1] += depth_allele
        i+=1
    total_ALT_frq = ALT_metagneome[1] + Other_metagenome[1]
    if total_count == 0:
        return ''
    elif total_ALT_frq == 0:
        return '%s\t%s\t\t\t\t\tno_mutation\t%s\t' %(REF_genome,REF_metagenome[1]/total_count,total_count)
    else:
        mutation_type = ALT_freq_mut_type(REF_metagenome[1],ALT_metagneome[1],Other_metagenome[1],total_count)
        return '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (REF_genome,REF_metagenome[1]/total_count,
                                             ALT_genome,ALT_metagneome[1]/total_count,
                                             Other_metagenome[0],Other_metagenome[1]/total_count,mutation_type,total_count)

all_output_select = []
all_output_select.append(
        '#donor_species\tmetagenome\tdonor\ttime\tN_or_S\tREF_genome\tREF_frq\tALT_genome\tALT_frq\tOther_ALT\tOther_ALT_frq\tmutation_type_norm\ttotal_count\tCHR\tPOS\n')

output_select_file = open(os.path.join(metagenome_summary,'all.filtered.selected.frq.snp'),'w')
output_select_file.write(''.join(all_output_select))
output_select_file.close()
all_output_select = []
i = 0
Mut_gene = dict()
for files in all_snp_sum:
    files_dir, file_name = os.path.split(files)
    print(file_name)
    i += 1
    if i%20 == 0:
        print('%s output results for %s samples' % (datetime.now(), i))
        output_select_file = open(os.path.join(metagenome_summary, 'all.filtered.selected.frq.snp'), 'a')
        output_select_file.write(''.join(all_output_select))
        output_select_file.close()
        all_output_select = []
    for lines in open(files,'r'):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            gene_ID = lines_set[0]
            pos = lines_set[1]
            if gene_ID in Gene_pos and pos in Gene_pos.get(gene_ID,[]):
                    #print('%s find a qualified gene in sample %s' % (datetime.now(), file_name))
                    donor_species, REF, ALT, N_or_S = Gene_info['%s_%s'%(gene_ID,pos)]
                    Subdepth = lines_set[9].split(':')[-1].replace('\n', '')
                    Subdepth_set = Subdepth.split(',')
                    REF_meta = lines_set[3]
                    ALT_meta = lines_set[4].split(',')
                    temp_result = ALT_freq(REF_meta,ALT_meta,Subdepth_set,REF,ALT)
                    #print('%s finishe processing gene %s' % (datetime.now(), gene_ID))
                    if temp_result!= '':
                        donor = file_name[0:2]
                        temp_line = '%s\t%s\t%s\t%s\t%s\t'%(donor_species,file_name.split('.raw.vcf')[0],donor,file_name[2:6],N_or_S)+ temp_result + '\t'.join(lines_set[0:2])
                        all_output_select.append(temp_line)
                        if 'no_mutation' not in temp_result:
                            Mut_gene.setdefault(donor,set())
                            Mut_gene[donor].add('%s_%s'%(gene_ID,pos))

output_select_file = open(os.path.join(metagenome_summary, 'all.filtered.selected.frq.snp'), 'a')
output_select_file.write(''.join(all_output_select))
output_select_file.close()

all_output_select = []
all_output_select.append(
        '#donor_species\tmetagenome\tdonor\ttime\tN_or_S\tREF_genome\tREF_frq\tALT_genome\tALT_frq\tOther_ALT\tOther_ALT_frq\tmutation_type_norm\ttotal_count\tCHR\tPOS\n')

for lines in open(os.path.join(metagenome_summary, 'all.filtered.selected.frq.snp'), 'r'):
    if not lines.startswith('#'):
        lines_set = lines.split('\n')[0].split('\t')
        gene_ID = lines_set[12]
        pos = lines_set[13]
        donor = lines_set[2]
        if '%s_%s'%(gene_ID,pos) in Mut_gene.get(donor,[]):
            all_output_select.append(lines)

output_select_file = open(os.path.join(metagenome_summary, 'all.filtered.selected.frq.snp.snp.only'), 'w')
output_select_file.write(''.join(all_output_select))
output_select_file.close()

