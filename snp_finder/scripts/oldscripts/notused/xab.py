# start
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
