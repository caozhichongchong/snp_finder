# start
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
