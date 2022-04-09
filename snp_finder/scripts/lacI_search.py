import os,glob
input_fasta = glob.glob('/scratch/users/mit_alm/IBD_Evo/BA/Assembly_for_gene_flow/*/scaffolds.fasta')+\
              glob.glob('/scratch/users/mit_alm/IBD_Evo/BL/Assembly_for_gene_flow/*/scaffolds.fasta')
print(len(input_fasta))
output_faa = '/scratch/users/anniz44/genomes/Jay/species/'
input_faa = glob.glob('/scratch/users/anniz44/genomes/GHM/newgenomes/GMC_Genome_Library_01312020/*.add')+\
              glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round1/*_BA/fasta/*.add')+\
              glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round1/*_BL/fasta/*.add')+\
              glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round1/*_BiPs/fasta/*.add')
output_folder = '/scratch/users/anniz44/genomes/HTH/diamond/'
output_scripts = '/scratch/users/anniz44/genomes/HTH/scripts/'
database = '/scratch/users/anniz44/genomes/HTH/TF_sequences_modified.fasta'
print(len(input_faa))
def run_diamond(input_faa,database,output_file):
    cmds = 'diamond blastp --query %s --db %s.dmnd --out %s --outfmt 6 --max-target-seqs 50 --evalue 0.01 --threads 40\n'%(input_faa,database,output_file)
    return cmds

def run_prodigal(input_fasta,output_aa,database,output_file):
    cmds = 'prodigal -q -a %s -i %s\n'%(output_aa,input_fasta)
    cmds += run_diamond(output_aa,database,output_file)
    return cmds

# process jay's data
cmdsset = '#!/bin/bash\nsource ~/.bashrc\n'
cmdsset += 'export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n'
i = 1
cmds = ''
for fasta in input_fasta:
    filename = os.path.split(os.path.split(fasta)[0])[-1]
    species = os.path.split(os.path.split(os.path.split(os.path.split(fasta)[0])[0])[0])[-1]
    output_aa = '%s/%s/%s.faa'%(output_faa,species,filename)
    output_file = '%s/%s.txt'%(output_folder,filename)
    cmds += run_prodigal(fasta,output_aa,database,output_file)
    if i%100 == 0:
        print('process %s fasta'%(i))
        f1 = open('%s/%s.fasta.sh'%(output_scripts,i),'w')
        f1.write(cmdsset)
        f1.write(cmds)
        f1.close()
        cmds = ''
    i+=1

# process GMBC + BN10
i = 1
cmds = ''
for fasta in input_faa:
    filename = os.path.split(fasta)[-1].split('_prokka')[0]
    output_file = '%s/%s.txt'%(output_folder,filename)
    cmds += run_diamond(fasta,database,output_file)
    if i%100 == 0:
        print('process %s faa' % (i))
        f1 = open('%s/%s.faa.sh'%(output_scripts,i),'w')
        f1.write(cmdsset)
        f1.write(cmds)
        f1.close()
        cmds = ''
    i += 1


