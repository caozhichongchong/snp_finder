import os,glob
from Bio import SeqIO

input_bs_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS/binding_results_ccpA.txt'
motif_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/BL_clustercluster33/motifs.meme'
script_file = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/fimo/'
script_file2 = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/predictaa/'

output_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS/'
os.system('#rm %s'%(script_file))
try:
    os.mkdir(script_file)
except IOError:
    pass
os.system('rm %s'%(script_file2))
try:
    os.mkdir(script_file2)
except IOError:
    pass
CL = dict()
for lines in open(input_bs_file,'r'):
    if not lines.startswith('AA_POS_ref'):
        lines_set = lines.split('\t')
        species = lines_set[4].split('_')[0]
        donor = lines_set[5]
        CL.setdefault(species,set())
        CL[species].add(donor)

for species in CL:
    for donor in CL[species]:
        output_file = '%s/%s_%s/' % (output_folder, species, donor)
        try:
            os.mkdir(output_file)
        except IOError:
            pass
        if donor.startswith('D') or donor.startswith('H'):
            input_fasta=glob.glob('/scratch/users/mit_alm/IBD_Evo/%s/Assembly_for_gene_flow/%s_*/scaffolds.fasta'%(species,donor))
            for fasta in input_fasta:
                cmds = '#!/bin/bash\nsource ~/.bashrc\n'
                genomename = os.path.split(os.path.split(fasta)[0])[-1]
                cmds += 'fimo -o %s/%s %s %s\n'%(output_file,genomename,motif_file,fasta)
                f1 = open('%s/%s.sh'%(script_file,genomename), 'w')
                f1.write(cmds)
                f1.close()
                try:
                    os.mkdir('%s/%s/' % (output_file, genomename))
                except IOError:
                    pass
                cmds = '#!/bin/bash\nsource ~/.bashrc\n'
                cmds += 'prodigal -q -i %s -a %s/%s/%s.faa' % (
                fasta, output_file, genomename,genomename)
                f1 = open('%s/%s.sh' % (script_file2, genomename), 'w')
                f1.write(cmds)
                f1.close()
        else:
            input_fasta = glob.glob(
                '/scratch/users/anniz44/genomes/donor_species/selected_species/round1/%s_%s/fasta/*_final.scaffolds.fasta' % (donor,species))
            for fasta in input_fasta:
                cmds = '#!/bin/bash\nsource ~/.bashrc\n'
                genomename = os.path.split(fasta)[-1].split('_final.scaffolds.fasta')[0]
                cmds += 'fimo -o %s/%s %s %s\n' % (output_file, genomename, motif_file, fasta)
                f1 = open('%s/%s.sh' % (script_file, genomename), 'w')
                f1.write(cmds)
                f1.close()
                try:
                    os.mkdir('%s/%s/' % (output_file, genomename))
                except IOError:
                    pass
                cmds = '#!/bin/bash\nsource ~/.bashrc\n'
                cmds += 'prodigal -q -i %s -a %s/%s/%s.faa' % (
                    fasta, output_file, genomename,genomename)
                f1 = open('%s/%s.sh' % (script_file2, genomename), 'w')
                f1.write(cmds)
                f1.close()

f1 = open(os.path.join(script_file, '../allfimo.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(script_file, '*.sh')):
    f1.write('jobmit %s %s small1\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/../allfimo.sh'%(script_file))

f1 = open(os.path.join(script_file, '../allpredictaa.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(script_file2, '*.sh')):
    f1.write('jobmit %s %s small1\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/../allpredictaa.sh'%(script_file2))