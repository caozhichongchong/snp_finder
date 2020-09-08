#---------------------------------------- pan-genome ----------------------------------------
# extract genomes for pan-genomes
input_script = '/scratch/users/anniz44/scripts/1MG/pan-genome'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/pan-genome/roary'
output_dir = '/scratch/users/anniz44/genomes/pan-genome'
#os.system('rm -rf %s' %(output_dir))
try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(output_dir + '/gubbins')
except IOError:
    pass

list_file = glob.glob(os.path.join(input_script,'*.all'))
Filelist = os.path.join(input_script,'allfilelist.txt')
Set = dict()
Species = dict()
Species_correct = dict()
for files in list_file:
    for lines in open(files,'r'):
        genome = lines.split('\t')[0]
        species = lines.split('\t')[-5].replace(' ','_')
        print(species)
        try:
            if len(species.split('_')[1])==1:
                Species_correct.setdefault(species, species.split('_')[0]+'_'+species.split('_')[2])
                print(species, species.split('_')[0]+'_'+species.split('_')[2])
                species = species.split('_')[0]+'_'+species.split('_')[2]
            if species.split('_')[1].startswith('sp'):
                Species_correct.setdefault(species, species.split('_')[0])
                print(species, species.split('_')[0])
                species = species.split('_')[0]
        except IndexError:
            pass
        Set.setdefault(genome,species)
        Species.setdefault(species,0)
        if not genome.startswith('GCA_'):
            Species[species] = 1

for species in dict(Species):
    if Species[species] == 0:
        print(species)
        try:
            newspecies = species.split('_')[0]+'_'+species.split('_')[1]
            if Species.get(newspecies,'None')== 1:
                    Species_correct.setdefault(species, newspecies)
                    Species_correct[species]=newspecies
            else:
                newspecies = species.split('_')[0]
                Species[newspecies]=1
                Species_correct.setdefault(species, newspecies)
                Species_correct[species] = newspecies
                print(species,newspecies)
        except IndexError:
            pass

i = 0
cmds2 = ''
for species in Species:
    if Species[species] == 1:
        try:
            cmds = ''
            output_subdir = os.path.join(output_dir, species)
            cmds += 'roary -p 40 -o roary_%s -e -n --mafft -f %s/roary_%s -i 80 -cd 90 %s/*.gff\n' % (
            species, output_dir,species, output_subdir)
            cmds2 += 'run_gubbins.py --threads 40  --tree_builder hybrid --use_time_stamp --prefix %s/gubbins/gubbins_%s_core_gene_alignment --verbose %s/roary_%s/core_gene_alignment.aln\n' % (
                output_dir, species, output_dir, species)
            #f1 = open(os.path.join(input_script_sub, '%s.roary.sh'%(i)), 'w')
            #f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            #f1.write('py37\n%s' % (''.join(cmds2)))
            #f1.close()
            i+=1
            os.mkdir(output_subdir)
        except OSError:
            pass

f1 = open(os.path.join(output_dir, '/gubbins/gubbins.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s' % (''.join(cmds2)))
f1.close()

for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.roary.sh')):
    print('jobmit %s %s' % (sub_scripts,os.path.split(sub_scripts)[-1]))


f1 = open(os.path.join(input_script,'species.correct.txt'),'w')
Species_correct_list=[]
for species in Species_correct:
    Species_correct_list.append('%s\t%s\t\n'%(species,Species_correct[species]))

f1.write('%s'%(''.join(Species_correct_list)))
f1.close()


cmds = ''
for lines in open(Filelist,'r'):
    if 'GCA' in lines:
        genome_dir = '_'.join(lines.split('_')[0:2])
    else:
        genome_dir = lines.split('_final.scaffolds.fasta')[0].split('.fasta')[0]
    genome_name = os.path.split(genome_dir)[-1]
    if genome_name in Set:
        species_name = Set[genome_name]
        if Species[species_name] == 0:
            species_name = Species_correct.get(species_name,'None')
        if species_name!= 'None':
            cmds += 'cp %s %s\n' %(genome_dir+'_*.fasta',os.path.join(output_dir,species_name + '/'+genome_name + '.fasta'))
            cmds += 'cp %s %s\n' % (genome_dir + '_*.gff', os.path.join(output_dir, species_name + '/'+genome_name+'.gff'))
            cmds += 'cp %s %s\n' % (genome_dir + '_*.fna', os.path.join(output_dir, species_name + '/'+genome_name + '.fasta'))
            cmds += 'cp %s %s\n' % (genome_dir + '_*.add', os.path.join(output_dir, species_name + '/'+genome_name+'.faa.add'))
        else:
            print(genome_name,Set[genome_name])

f1 = open(os.path.join(input_script,'copy.sh'),'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(''.join(cmds)))
f1.close()

#flexible gneome selection
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

input_script = '/scratch/users/anniz44/scripts/1MG/pan-genome'
input_dir1 = '/scratch/users/anniz44/genomes/pan-genome/roary/roary_*'
input_dir2 = '/scratch/users/anniz44/genomes/donor_species/roary/roary_*'
input_dir3 = '/scratch/users/anniz44/genomes/Jay/roary/roary_*'
output_dir = '/scratch/users/anniz44/genomes/pan-genome/allpangenome/subcluster'

try:
    os.mkdir(output_dir)
except IOError:
    pass

output_file = open(os.path.join(output_dir,  ('../sumsub/all.flexible.genome.annotation')), 'w')
output_file.write('species\tgeneID\tfunction\tannotation\tgene_copy\tgene_length\n')
output_file.close()
Annotation = dict()

roary_folder = glob.glob(input_dir1)
for roary_folder1 in roary_folder:
    species = roary_folder1.split('roary_')[1]
    species2 = species.split('_')[0]
    pan_genome_annotation = os.path.join(roary_folder1,'gene_presence_absence.csv')
    pan_genome_fa = os.path.join(roary_folder1, 'pan_genome_reference.fa')
    accessory_tab = os.path.join(roary_folder1, 'accessory.tab')
    accessory_gene = dict()
    accessory_gene_length = dict()
    output_annotation = []
    output_fasta = []
    print(species, pan_genome_fa,accessory_tab)
    for lines in open(accessory_tab,'r'):
        if '/gene=' in lines:
            accessory_gene.setdefault(lines.split('/gene=')[1].split('\n')[0],'')
    for record in SeqIO.parse(pan_genome_fa, 'fasta'):
        geneID = str(record.description).split(str(record.id)+' ')[1]
        if geneID in accessory_gene:
            recordID = str(record.id)
            if accessory_gene[geneID]!='':
                print(species,geneID,accessory_gene[geneID],recordID)
            accessory_gene[geneID]=recordID
            gene_length = len(str(record.seq))
            accessory_gene_length.setdefault(geneID, gene_length)
            #output_fasta.append('>%s\t%s\n%s\n' % (fun_name, str(record.id), str(record.seq)))
    #output_file = open(os.path.join(output_dir, '%s_%s' % (species, 'flexible.genome.fa')), 'a')
    #output_file.write(''.join(output_fasta))
    #output_file.close()
    for lines in open(pan_genome_annotation, 'r'):
        lines_set = lines.replace('\n','').replace('\r','').split(',')
        geneID = lines_set[0].replace('\"', '')
        if geneID in accessory_gene:
            annotation = geneID + '\t' + lines_set[2].replace('\"', '')
            annotation1 = geneID.split('_')[0] + '\t' + lines_set[2].replace('\"', '')
            recordID = accessory_gene[geneID]
            gene_length = accessory_gene_length.get(geneID, 1000)
            copy_population = str(lines_set[4].replace('\"',''))
            output_annotation.append('%s\t%s\t%s\t%s\t%s\t\n' % (species, recordID, annotation,copy_population,gene_length))
            if 'group' not in geneID:
                Annotation.setdefault(annotation1, set())
                Annotation[annotation1].add(species2)
    output_file = open(os.path.join(output_dir, ('../sumsub/all.flexible.genome.annotation')), 'a')
    output_file.write(''.join(output_annotation))
    output_file.close()

roary_folder = glob.glob(input_dir2)
for roary_folder1 in roary_folder:
    if '_newcluster' not in roary_folder1:
        # for donor species
        species = '_'.join(roary_folder1.split('roary_')[1].split('_cluster')[0].split('_newcluster')[0].split('_')[1:])
        species2 = species.split('_')[0]
        pan_genome_annotation = os.path.join(roary_folder1, 'gene_presence_absence.csv')
        pan_genome_fa = os.path.join(roary_folder1, 'pan_genome_reference.fa')
        accessory_tab = os.path.join(roary_folder1, 'accessory.tab')
        accessory_gene = dict()
        accessory_gene_length = dict()
        output_annotation = []
        output_fasta = []
        print(species, pan_genome_fa, accessory_tab)
        for lines in open(accessory_tab, 'r'):
            if '/gene=' in lines:
                accessory_gene.setdefault(lines.split('/gene=')[1].split('\n')[0], '')
        for record in SeqIO.parse(pan_genome_fa, 'fasta'):
            geneID = str(record.description).split(str(record.id) + ' ')[1]
            if geneID in accessory_gene:
                recordID = str(record.id)
                if accessory_gene[geneID] != '':
                    print(species, geneID, accessory_gene[geneID], recordID)
                accessory_gene[geneID] = recordID
                gene_length = len(str(record.seq))
                accessory_gene_length.setdefault(geneID, gene_length)
                # output_fasta.append('>%s\t%s\n%s\n' % (fun_name, str(record.id), str(record.seq)))
        # output_file = open(os.path.join(output_dir, '%s_%s' % (species, 'flexible.genome.fa')), 'a')
        # output_file.write(''.join(output_fasta))
        # output_file.close()
        for lines in open(pan_genome_annotation, 'r'):
            lines_set = lines.replace('\n','').replace('\r','').split(',')
            geneID = lines_set[0].replace('\"', '')
            if geneID in accessory_gene:
                annotation = geneID + '\t' + lines_set[2].replace('\"', '')
                annotation1 = geneID.split('_')[0] + '\t' + lines_set[2].replace('\"', '')
                recordID = accessory_gene[geneID]
                gene_length = accessory_gene_length.get(geneID,1000)
                copy_population = str(lines_set[4].replace('\"',''))
                output_annotation.append(
                    '%s\t%s\t%s\t%s\t%s\t\n' % (species, recordID, annotation, copy_population, gene_length))
                if 'group' not in geneID:
                    Annotation.setdefault(annotation1, set())
                    Annotation[annotation1].add(species2)
        output_file = open(os.path.join(output_dir, ('../sumsub/all.flexible.genome.annotation')), 'a')
        output_file.write(''.join(output_annotation))
        output_file.close()

os.system('cp %s %s'%(os.path.join(output_dir, ('../sumsub/all.flexible.genome.annotation')),
                      os.path.join(output_dir, ('../sumsub/all.flexible.genome.noJay.annotation'))))

roary_folder = glob.glob(input_dir3)
for roary_folder1 in roary_folder:
    if '_newcluster' not in roary_folder1:
        # for donor species
        species = 'Jay_Bacteroides_' + (roary_folder1.split('roary_')[1].split('_cluster')[0].split('_newcluster')[0])
        species2 = 'Bacteroides'
        pan_genome_annotation = os.path.join(roary_folder1, 'gene_presence_absence.csv')
        pan_genome_fa = os.path.join(roary_folder1, 'pan_genome_reference.fa')
        accessory_tab = os.path.join(roary_folder1, 'accessory.tab')
        accessory_gene = dict()
        accessory_gene_length = dict()
        output_annotation = []
        output_fasta = []
        print(species, pan_genome_fa, accessory_tab)
        for lines in open(accessory_tab, 'r'):
            if '/gene=' in lines:
                accessory_gene.setdefault(lines.split('/gene=')[1].split('\n')[0], '')
        for record in SeqIO.parse(pan_genome_fa, 'fasta'):
            geneID = str(record.description).split(str(record.id) + ' ')[1]
            if geneID in accessory_gene:
                recordID = str(record.id)
                if accessory_gene[geneID] != '':
                    print(species, geneID, accessory_gene[geneID], recordID)
                accessory_gene[geneID] = recordID
                gene_length = len(str(record.seq))
                accessory_gene_length.setdefault(geneID, gene_length)
                # output_fasta.append('>%s\t%s\n%s\n' % (fun_name, str(record.id), str(record.seq)))
        # output_file = open(os.path.join(output_dir, '%s_%s' % (species, 'flexible.genome.fa')), 'a')
        # output_file.write(''.join(output_fasta))
        # output_file.close()
        for lines in open(pan_genome_annotation, 'r'):
            lines_set = lines.replace('\n','').replace('\r','').split(',')
            geneID = lines_set[0].replace('\"', '')
            if geneID in accessory_gene:
                annotation = geneID + '\t' + lines_set[2].replace('\"', '')
                annotation1 = geneID.split('_')[0] + '\t' + lines_set[2].replace('\"', '')
                recordID = accessory_gene[geneID]
                gene_length = accessory_gene_length.get(geneID, 1000)
                copy_population = str(lines_set[4].replace('\"',''))
                output_annotation.append(
                    '%s\t%s\t%s\t%s\t%s\t\n' % (species, recordID, annotation, copy_population, gene_length))
                if 'group' not in geneID:
                    Annotation.setdefault(annotation1, set())
                    Annotation[annotation1].add(species2)
        output_file = open(os.path.join(output_dir, ('../sumsub/all.flexible.genome.annotation')), 'a')
        output_file.write(''.join(output_annotation))
        output_file.close()

output_file = open(os.path.join(output_dir, ('../sumsub/all.flexible.genome.multispecies')), 'w')
output_annotation = []
for annotation1 in Annotation:
    number_species = len(Annotation[annotation1])
    if number_species > 1:
        output_annotation.append('%s\t%s\t%s\n'%(annotation1,number_species,';'.join(Annotation[annotation1])))

output_file.write(''.join(output_annotation))
output_file.close()

#roary_folder = glob.glob(input_dir3)
#for roary_folder1 in roary_folder:
#    species = roary_folder1.split('roary_')[1]
#    pan_genome_fa = os.path.join(roary_folder1,'pan_genome_reference.fa')
#    accessory_tab = os.path.join(roary_folder1,'accessory.tab')
#    accessory_gene = []
#    output_fasta = []
#    print(species, pan_genome_fa,accessory_tab)
#    os.system('mv %s %s/'%(os.path.join(output_dir,'%s_%s'%(species,'flexible.genome.fa')),os.path.join(output_dir,'Jay_Bacteroides_%s_%s'%(species,'flexible.genome.fa'))))

# cluster genes
import os
import glob
output_dir_sub = '/scratch/users/anniz44/genomes/pan-genome/allpangenome/subcluster'
output_dir = '/scratch/users/anniz44/genomes/pan-genome/allpangenome'
input_script = '/scratch/users/anniz44/scripts/1MG/pan-genome'
cutoff = 0.7

try:
    os.mkdir(output_dir_sub)
except IOError:
    pass

cmd = ''
all_filter_gene_fasta_file = glob.glob(os.path.join(output_dir, 'Jay_Bacteroides*flexible.genome.fa'))
for files in all_filter_gene_fasta_file:
    if 'Jay_Bacteroides_all_flexible.genome.fa' not in files:
        outputfile = os.path.join(output_dir_sub,os.path.split(files)[-1])
        try:
            f1 = open(outputfile +'.uc','r')
        except IOError:
            cmd += ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                                      % ('usearch', files, cutoff, outputfile, outputfile, 40))

cmd +=('cat %s/Jay_Bacteroides_*.fa.cluster.aa > %s/Jay_Bacteroides_all_flexible.genome.fa\n' %(output_dir_sub,output_dir))
files = '%s/Jay_Bacteroides_all_flexible.genome.fa'%(output_dir)
cmd += ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                                      % ('usearch', files, cutoff, files, files, 40))

f1 = open(os.path.join(input_script, 'all.pan-genome.usearch.sh'),'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n'+cmd)
f1.close()


# filter out functions carried by multiple species
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

input_script = '/scratch/users/anniz44/scripts/1MG/pan-genome'
output_dir = '/scratch/users/anniz44/genomes/pan-genome/allpangenome/sumsub'

# use cluster information
Candidate_geneID = set()
Annotation = dict()
annotation_file = os.path.join(output_dir,'all.flexible.genome.annotation')

for lines in open(annotation_file):
        lines_set = lines.split('\n')[0].split('\t')
        species = lines_set[0]
        if 'Jay' in species:
            species = 'Bacteroides'
        else:
            species = species.split('_')[0]
        geneID = lines_set[1]
        annotation = lines_set[3]
        Annotation.setdefault(geneID,[species,annotation])

cluster_file = glob.glob(os.path.join(output_dir,'*.uc'))
Cluster = dict()
for a_cluster_file in cluster_file:
    for lines in open(a_cluster_file):
        lines_set = lines.split('\n')[0].split('\t')
        if len(lines_set)>11:
            geneID = lines_set[-1].split('\n')[0]
            geneID2 = lines_set[9]
            if geneID in Annotation:
                species = Annotation[geneID][0]
                Cluster.setdefault(geneID,set())
                Cluster[geneID].add(species)
                if geneID2 in Annotation:
                    species = Annotation[geneID2][0]
                    Cluster[geneID].add(species)

output_file = open(os.path.join(output_dir, 'all.flexible.genome.multispecies.bycluster'), 'w')
output_annotation = []
for geneID in Cluster:
    number_species = len(Cluster[geneID])
    if number_species > 1:
        Candidate_geneID.add(geneID)
        output_annotation.append('%s\t%s\t%s\t%s\n'%(geneID,Annotation[geneID][1], number_species,';'.join(Cluster[geneID])))

output_file.write(''.join(output_annotation))
output_file.close()

# use annotation
Candidate_function = set()
for lines in open(os.path.join(output_dir,'all.flexible.genome.multispecies'),'r'):
    function_name =lines.split('\t')[0]
    Candidate_function.add(function_name)

# output multiple species
for fasta_files in glob.glob(os.path.join(output_dir,'*cluster.aa')) + glob.glob(os.path.join(output_dir,'*.fasta')):
    multi_species_fasta = []
    for record in SeqIO.parse(fasta_files, 'fasta'):
        recordID = str(record.description).replace(' ','\t')
        function_name = recordID.split('_')[0]
        geneID = '%s_%s' %(recordID.split('_')[-2], str(record.id).split('_')[-1])
        if function_name in Candidate_function or geneID in Candidate_geneID:
            multi_species_fasta.append('>%s\n%s\n' % (recordID, str(record.seq)))
    output_file = open(fasta_files + '.multispecies.fasta', 'w')
    output_file.write(''.join(multi_species_fasta))
    output_file.close()


# Jay pangenome
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

input_script = '/scratch/users/anniz44/scripts/1MG/pan-genome'
input_dir = '/scratch/users/anniz44/genomes/Jay/data'
output_dir = '/scratch/users/anniz44/genomes/Jay'

try:
    os.mkdir(output_dir + '/species')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/roary')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/gubbins')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/prokka')
except IOError:
    pass

try:
    os.mkdir(input_script + '/Jay_roary')
except IOError:
    pass

try:
    os.mkdir(input_script + '/Jay_prokka')
except IOError:
    pass

i=0
cmdsgubbins = ''
fasta_all = glob.glob(os.path.join(output_dir + '/species/*/','*.fasta'))
donor_species_set = []
for fasta_all_file in fasta_all:
    genome = os.path.split(fasta_all_file)[-1].split('.fasta')[0]
    donor_species = fasta_all_file.split('_')[-1].split('.')[0]
    species_dir = os.path.join(output_dir +'/species',donor_species)
    roary_dir = output_dir + '/roary/roary_' + donor_species
    try:
        os.mkdir(roary_dir)
    except IOError:
        pass
    try:
        os.mkdir(species_dir)
    except IOError:
        pass
    cmdsprokka = ''
    cmdsprokka += '#!/bin/bash\nsource ~/.bashrc\npy37\nprokka --kingdom Bacteria --outdir %s/prokka/prokka_%s --genus Bacteroides --locustag Bacter %s\n' %(output_dir, genome,fasta_all_file)
    cmdsprokka += '#mv %s/prokka/prokka_%s/*.gff %s/%s.gff\n' %(output_dir,genome, species_dir,genome)
    cmdsprokka += '#mv %s/prokka/prokka_%s/*.tsv %s/%s.tsv\n' % (output_dir, genome, species_dir, genome)
    cmdsprokka += '#mv %s %s/\n' % (fasta_all_file, species_dir)
    cmdsprokka += '#mv %s %s/\n' % (fasta_all_file.replace('.fasta','.faa'), species_dir)
    cmdsprokka += '#mv %s/%s*.add %s\n' % (input_dir,genome, species_dir)
    if donor_species not in donor_species_set:
        cmdsroary = ''
        cmdsroary += 'roary -p 40 -o roary_%s -e -n --mafft -f %s -i 80 -cd 95 %s/*.gff\n' % (
            donor_species, roary_dir, species_dir)
        cmdsgubbins += 'run_gubbins.py --threads 40  --tree_builder hybrid --use_time_stamp --prefix %s/gubbins/gubbins_%s_core_gene_alignment --verbose %s/core_gene_alignment.aln\n' % (
            output_dir, donor_species, roary_dir)
        f1 = open(os.path.join(input_script + '/Jay_roary', '%s.roary.sh' % (i)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmdsroary)))
        f1.close()
        donor_species_set.append(donor_species)
    f1 = open(os.path.join(input_script + '/Jay_prokka', '%s.prokka.sh' % (i)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmdsprokka)))
    f1.close()
    i += 1

f1 = open(os.path.join(input_script, 'Jaygubbins.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s' % (''.join(cmdsgubbins)))
f1.close()

f1 = open(os.path.join(input_script, 'Jayroarynew.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script + '/Jay_roary', '*.roary.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'Jayprokkanew.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script + '/Jay_prokka', '*.prokka.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

# metagenomes running
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
fasta_dir = '/scratch/users/anniz44/genomes/pan-genome/allpangenome'
fasta_file = os.path.join(fasta_dir,'all_flexible.genome.fa.cluster.aa')
fasta_file_set = []
fasta_file_set_length = []
for record in SeqIO.parse(fasta_file, 'fasta'):
    record_id = str(record.id)
    fasta_file_set.append('>%s\n%s\n' % (str(record.description).replace('\t','_'), str(record.seq)))
    fasta_file_set_length.append('%s\t%s\n'%(str(record.description).replace('\t','_'), len(str(record.seq))))

output_file = open(fasta_file + '.nom.fasta','w')
output_file.write(''.join(fasta_file_set))
output_file.close()
output_file = open(fasta_file + '.length','w')
output_file.write(''.join(fasta_file_set_length))
output_file.close()
os.system('mv %s %s/.old' %(fasta_file,fasta_file))
os.system('mv %s.nom.fasta %s' %(fasta_file,fasta_file))

import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
fasta_dir = '/scratch/users/anniz44/genomes/pan-genome/allpangenome'
fasta_file = os.path.join(fasta_dir,'Jay_Bacteroides_all_flexible.genome.fa.cluster.aa')
fasta_file_set = []
fasta_file_set_length = []
for record in SeqIO.parse(fasta_file, 'fasta'):
    record_id = str(record.id)
    fasta_file_set.append('>%s\n%s\n' % (str(record.description).replace('\t','_'), str(record.seq)))
    fasta_file_set_length.append('%s\t%s\n'%(str(record.description).replace('\t','_'), len(str(record.seq))))

output_file = open(fasta_file + '.nom.fasta','w')
output_file.write(''.join(fasta_file_set))
output_file.close()
output_file = open(fasta_file + '.length','w')
output_file.write(''.join(fasta_file_set_length))
output_file.close()
os.system('mv %s %s/.old' %(fasta_file,fasta_file))
os.system('mv %s.nom.fasta %s' %(fasta_file,fasta_file))

# metagenomes enriched genes
#from matplotlib import style
import os
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
#import skbio
#from scipy.spatial import distance
#from sklearn import datasets
#from sklearn.decomposition import PCA
#from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
#import ecopy as ep
#from sklearn import manifold

genemultispecies = '/scratch/users/anniz44/genomes/pan-genome/allpangenome/all_flexible.genome.fa.cluster.aa.multispecies.fasta'
geneenrich = '/scratch/users/anniz44/genomes/pan-genome/allpangenome/subcluster/all_flexible.genome.fa.cluster.aa.enriched.fasta.multispecies.fasta'

genemultispecies_set = []
for record in SeqIO.parse(genemultispecies, 'fasta'):
    record_id = str(record.id)
    genemultispecies_set.append(record_id)

geneenrich_set = []
for record in SeqIO.parse(geneenrich, 'fasta'):
    record_id = str(record.id)
    geneenrich_set.append(record_id)

summary_filename = glob.glob(os.path.join(input_dir,'*.gene.summary.cell.txt'))[0]

candidate_set = set(geneenrich_set).intersection(set(genemultispecies_set))

colnames_set = [0]
summary_filename_set = []
for lines in open(summary_filename,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    candidate_set2 = set(lines_set).intersection(set(candidate_set))
    if lines.startswith('SampleID'):
        for geneID in candidate_set2:
            colnames_set.append(lines_set.index(geneID))
    print('finished finding a subset')
    summary_filename_set.append('\t'.join(list(np.array(lines_set)[colnames_set]))+'\n')

output_file = open(summary_filename + '.multispecies.humanenrich','w')
output_file.write(''.join(summary_filename_set))
output_file.close()

import os
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

input_dir = '/scratch/users/anniz44/Metagenomes/pangenome/summary_Jay'
metadata = '/scratch/users/anniz44/scripts/1MG/all_MGD_GMD_metagenome.metadata.txt'

summary_filename = glob.glob(os.path.join(input_dir,'*.gene.summary.cell.txt'))[0]

metadata_dict = dict()
metadata_type = dict()
for lines in open(metadata,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    sampleID = lines_set[0]
    metainfo = lines_set[-2]
    metadata_dict.setdefault(sampleID, metainfo)
    metadata_type.setdefault(metainfo,[])


OTU_table = pd.read_csv(summary_filename,
sep='\t',header=0,index_col=0)
df = pd.DataFrame(OTU_table)
samples = df.index
i = 0
for a_sample in samples:
    metadata_type[metadata_dict[a_sample]].append(i)
    i+=1

all_avg = []
colnames = []
for enviroment in metadata_type:
    colnames.append(enviroment+'_abu')
    colnames.append(enviroment + '_pre')
    sub_set=metadata_type[enviroment]
    sub_set_avg = df.iloc[sub_set].sum(axis=0)/len(sub_set) # for a subset of samples, sum gene
    sub_set_avg_pre = df.iloc[sub_set].fillna(0).astype(bool).sum(axis=0) / len(sub_set)  # for a subset of samples, count gene presence
    all_avg.append(sub_set_avg)
    all_avg.append(sub_set_avg_pre)
    #sub_set_avg.to_csv(os.path.join(input_dir, 'all_flexible.genome.fa.cluster.aa.gene.summary.cell.%s.txt')%(enviroment), sep='\t',
    #          header=True)

all_avg_result = pd.concat(all_avg,axis = 1,names = colnames)
all_avg_result.columns = colnames
all_avg_result.to_csv(summary_filename + '.sumenv.txt', sep='\t',
              header=True)

# R first then upload file */*enrich
# pick up enriched genes
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq

input_sum = '/scratch/users/anniz44/Metagenomes/pangenome/summary_Jay'
fasta_dir = '/scratch/users/anniz44/genomes/pan-genome/allpangenome'
enrich_file1 = glob.glob(os.path.join(input_sum,'*enrich.1e6.txt'))[0]
enrich_file1_set = []
enrich_file2 = glob.glob(os.path.join(input_sum,'*.enrich.1e2.human.txt'))[0]
enrich_file2_set = []
#fasta_file = os.path.join(fasta_dir,'all_flexible.genome.fa.cluster.aa.multispecies.fasta')
fasta_file = os.path.join(fasta_dir,'Jay_Bacteroides_all_flexible.genome.fa.cluster.aa.multispecies.fasta')
for lines in open(enrich_file1,'r'):
    enrich_file1_set.append(lines.split('\t')[0])

for lines in open(enrich_file2,'r'):
    enrich_file2_set.append(lines.split('\t')[0])

output_file1_set = []
output_file2_set = []
for record in SeqIO.parse(fasta_file, 'fasta'):
    record_id = str(record.id)
    if record_id in enrich_file1_set:
        output_file1_set.append('>%s\t%s\n%s\n' %(record_id,str(record.description),str(record.seq)))
    if record_id in enrich_file2_set:
        output_file2_set.append('>%s\t%s\n%s\n' %(record_id,str(record.description),str(record.seq)))

output_file = open(fasta_file + '.enriched.fasta','w')
output_file.write(''.join(output_file1_set))
output_file.close()

output_file = open(fasta_file + '.enriched.human.fasta','w')
output_file.write(''.join(output_file2_set))
output_file.close()

#fasta_file = os.path.join(fasta_dir,'subcluster/all_flexible.genome.fa.cluster.aa')

output_file1_set = []
output_file2_set = []
for lines in open(fasta_file + '.length', 'r'):
    record_id = str(lines).split('\t')[0]
    if record_id in enrich_file1_set:
        output_file1_set.append(lines)
    if record_id in enrich_file2_set:
        output_file2_set.append(lines)

#fasta_file = os.path.join(fasta_dir,'all_flexible.genome.fa.cluster.aa.multispecies.fasta')
output_file = open(fasta_file + '.enriched.length','w')
output_file.write(''.join(output_file1_set))
output_file.close()

output_file = open(fasta_file + '.enriched.human.length','w')
output_file.write(''.join(output_file2_set))
output_file.close()

# gene selected percentage
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq

# including Jay pangenome
input_script = '/scratch/users/anniz44/scripts/1MG/pan-genome'
input_dir1 = '/scratch/users/anniz44/genomes/pan-genome/roary/roary_*'
input_dir2 = '/scratch/users/anniz44/genomes/donor_species/roary/roary_*'
input_dir3 = '/scratch/users/anniz44/genomes/Jay/roary/roary_*'
output_dir = '/scratch/users/anniz44/genomes/pan-genome/allpangenome'
mutation_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2/summary'
mutation_dir2 = '/scratch/users/anniz44/genomes/donor_species/species/vcf_merge_round2'
mutation_fasta = os.path.join(mutation_dir,'all.selected.gene.faa')
all_fasta = os.path.join(output_dir,'all_flexible.Jay.cluster.multispecies.fasta')
annotation_file = os.path.join(output_dir,'all.flexible.genome.annotation')
#enrichment_file = os.path.join(output_dir,'all_flexible.Jay.cluster.multispecies.enriched.human.fasta')
enrichment_file = os.path.join(output_dir,'all_flexible.Jay.cluster.multispecies.enriched.fasta')
mutation_file = os.path.join(mutation_dir,'all.donor.species.sum.txt.High_select2.txt')
mutation_mapping_file = os.path.join(output_dir,'all_flexible.Jay.cluster.multispecies.fasta.selectedgene.usearch.txt')
pangenome_info = os.path.join(output_dir,'sumsub/all.mutation.sum') # generate previously by mapping genomes to pan-genome genes

# map enriched genes to mutation genes
cutoff = 50
cutoff2 = 80
cmds = ("#diamond blastx --query %s --db %s.dmnd --out %s --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_fasta,mutation_fasta,mutation_mapping_file,cutoff,cutoff2))
os.system(cmds)

# calculate pangenome size and copy
donor_species_pan = dict()
donor_species_mut = dict()
for lines in open(pangenome_info,'r'):
    if not lines.startswith('#'):
        lines_set = lines.split('\n')[0].split('\t')
        donor_species = lines_set[0]
        pan_genome_length = float(lines_set[1])
        genome_depth = float(lines_set[2])
        donor_species_pan.setdefault(donor_species,[pan_genome_length,genome_depth])
        donor_species2 = '_'.join(lines_set[0].split('_cluster')[0].split('_newcluster')[0].replace('_sp.','_sp').split('_')[1:])
        if donor_species2 not in donor_species_mut:
            donor_species_mut.setdefault(donor_species2, [float(pan_genome_length) * float(genome_depth),
                                                          0])
        else:
            donor_species_mut[donor_species2][0] += float(pan_genome_length) * float(genome_depth)

# load mutation rate
donor_species_mut_gene = dict()
donor_species_mut_list = []
for lines in open(mutation_file,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    geneID = lines_set[1]
    donor_species = lines_set[0]
    total_SNP = lines_set[5]
    if not donor_species.startswith('allspecies'):
        if donor_species == geneID:
            donor_species2 = '_'.join(lines_set[0].split('_cluster')[0].split('_newcluster')[0].replace('_sp.','_sp').split('_')[1:])
            pan_genome_length,genome_depth = donor_species_pan.get(donor_species,[1,0])
            if donor_species2 not in donor_species_mut:
                donor_species_mut.setdefault(donor_species2,[float(pan_genome_length)*float(genome_depth),
                                                             float(total_SNP)])
            else:
                donor_species_mut[donor_species2][1] += float(total_SNP)
            donor_species_mut_list.append('%s\t%s\t%s\t%s\t\n'%(donor_species,pan_genome_length,genome_depth,total_SNP))
        elif lines_set[-1]=='True':
            donor_species_set = donor_species.split('_')
            donor_species3 = '%s_%s_%s' % (donor_species_set[0],
                                          donor_species_set[1][0:min(6, len(donor_species_set[1]))],
                                          donor_species_set[2][0:min(6, len(donor_species_set[2]))])
            record_id = '%s__C_%s_G_%s' % (donor_species3,geneID.split('_')[1], geneID.split('_')[-1])
            donor_species_mut_gene.setdefault(record_id,float(total_SNP))

output_file = open(os.path.join(mutation_dir,'all.mutation.sum'),'w')
output_file.write('#donor_species\tpan_genome_length\tavg_pangenome_copy\ttotal_SNP\n'+''.join(donor_species_mut_list))
output_file.close()

# load all genes
Cluster = set()
for lines in open(all_fasta,'r'):
    if lines.startswith('>'):
        lines_set = lines.split('\n')[0].split('\t')
        geneID = '_'.join(lines_set[0].split('_')[-2:])
        Cluster.add(geneID)

# load enriched genes
Enrichment = set()
for lines in open(enrichment_file,'r'):
    if lines.startswith('>'):
        lines_set = lines.split('\n')[0].split('\t')
        geneID = '_'.join(lines_set[0].split('_')[-2:])
        if geneID not in Cluster:
            print('wrong geneID not in cluster',geneID)
        else:
            Enrichment.add(geneID)

# load annotation files
Annotation = dict()
i = 0
for lines in open(annotation_file,'r'):
    if not lines.startswith('species'):
        try:
            lines_set = lines.split('\n')[0].split('\t')
            geneID = lines_set[1]
            species = lines_set[0]
            gene_copy_length = float(lines_set[4]) * float(lines_set[5])
            Annotation.setdefault(geneID, [species,gene_copy_length])
        except ValueError:
            print('wrong line',lines_set)
            i+=1

# load mutation mapping
mutation_mapping = dict()
for lines in open(mutation_mapping_file,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    geneID = '_'.join(lines_set[0].split('_')[-2:])
    geneID2 = lines_set[1]
    mutation_mapping.setdefault(geneID, geneID2)
    if geneID not in Annotation:
        print('wrong geneID not in annotation %s'%(geneID))

def annotate_gene(geneID,Annotation,donor_species_mut,Enriched_set,mutation_mapping,gene_SNP=0):
    species = ''
    gene_copy_length = 0
    pan_genome_length_depth = 0
    total_SNP = 0
    if geneID in Annotation:
        species, gene_copy_length = Annotation[geneID]
        if species in donor_species_mut:
            pan_genome_length_depth, total_SNP = donor_species_mut[species]
        else:
            for allspecies in donor_species_mut:
                if species in allspecies:
                    pan_genome_length_depth, total_SNP = donor_species_mut[allspecies]
                    break
        print(species,pan_genome_length_depth, total_SNP)
    if geneID in mutation_mapping:
        gene_SNP = donor_species_mut_gene.get(mutation_mapping[geneID],0)
    if total_SNP > 0:
        Enriched_set[0] += gene_SNP/total_SNP
        Enriched_set[1] += gene_copy_length / pan_genome_length_depth
    return [species,gene_copy_length,gene_SNP,pan_genome_length_depth,total_SNP]

Enriched = [0,0]
Enriched_list = []
Not_enriched = [0,0]
Not_enriched_list = []
for geneID in Cluster:
    if geneID in Enrichment:
        if geneID in Annotation:
            species, gene_copy_length, gene_SNP_geneID, pan_genome_length_depth, total_SNP = \
                annotate_gene(geneID,Annotation,donor_species_mut,Enriched,mutation_mapping)
            Enriched_list.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n' % (
                geneID, geneID, species, gene_copy_length, gene_SNP_geneID, pan_genome_length_depth, total_SNP))
    else:
        if geneID in Annotation:
            species, gene_copy_length, gene_SNP_geneID, pan_genome_length_depth, total_SNP = \
                annotate_gene(geneID,Annotation,donor_species_mut,Not_enriched,mutation_mapping)
            Not_enriched_list.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n' % (
                geneID, geneID, species, gene_copy_length, gene_SNP_geneID, pan_genome_length_depth, total_SNP))

output_file = open(os.path.join(mutation_dir,'all.mutation.enriched.nonhuman.sum'),'w')
output_file.write('#geneID\tgeneID2\tspecies\tgene_copy_length\tgene_SNP\tpan_genome_length_depth\ttotal_SNP\t\n'+''.join(Enriched_list))
output_file.close()

output_file = open(os.path.join(mutation_dir,'all.mutation.notenriched.nonhuman.sum'),'w')
output_file.write('#geneID\tgeneID2\tspecies\tgene_copy_length\tgene_SNP\tpan_genome_length_depth\ttotal_SNP\t\n'+''.join(Not_enriched_list))
output_file.close()

output_file = open(os.path.join(mutation_dir,'all.mutation.enriched.ratio.nonhuman.sum'),'w')
output_file.write('#enriched_SNP_ratio\tenriched_length_ratio\tnon_enriched_SNP_ratio\tnon_enriched_length_ratio\tmutation_rate_enriched\tmutation_rate_nonenriched\tenriched_to_nonenriched\t\n'+
                  '%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n'%(Enriched[0],Enriched[1],Not_enriched[0],Not_enriched[1],
                                        Enriched[0]/Enriched[1],Not_enriched[0]/Not_enriched[1],
                                        Enriched[0] / Enriched[1]/Not_enriched[0]* Not_enriched[1]))
output_file.close()
# clean up
import os
os.system('#rm -rf /scratch/users/anniz44/genomes/Jay/data')
os.system('#rm -rf /scratch/users/anniz44/genomes/Jay/prokka')
os.system('rm -rf /scratch/users/anniz44/Metagenomes/donor_species/bwa/0/*.flt.vcf')
os.system('rm -rf /scratch/users/anniz44/Metagenomes/donor_species/bwa/0/*.sorted.bam.*')
os.system('rm -rf  roary_*/pan_genome_sequences')
os.system('rm -rf  roary_*/_*')
os.system('rm -rf  roary_*/*.Rtab')
os.system('rm -rf  roary_*/fixed_input_files')
os.system('rm -rf  /scratch/users/anniz44/Metagenomes/pangenome_all/usearch')
os.system('rm -rf  /scratch/users/anniz44/Metagenomes/pangenome_all/search_output/0/*.txt.filter.aa')
os.system('rm -rf  /scratch/users/anniz44/Metagenomes/pangenome_all/search_output/0/*.blast.txt')
os.system('rm -rf  /scratch/users/anniz44/Metagenomes/pangenome/usearch')
os.system('rm -rf  /scratch/users/anniz44/Metagenomes/pangenome/search_output/0/*.txt.filter.aa')
os.system('rm -rf  /scratch/users/anniz44/Metagenomes/pangenome/search_output/0/*.blast.txt')
# after traits_finder sum_meta
os.system('rm -rf  /scratch/users/anniz44/Metagenomes/pangenome_all/search_output')
os.system('rm -rf  /scratch/users/anniz44/Metagenomes/pangenome/search_output')

