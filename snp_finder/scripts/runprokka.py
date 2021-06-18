# prokka of all ref genomes
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
all_genomes = '/scratch/users/anniz44/genomes/donor_species//vcf_round2/co-assembly/*/*.noHM.fasta'

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
                os.system('%s -q -i %s -d %s.fna -a %s.faa' % ('prodigal', genome, genome,genome))
            shortname(genome,True)
            #shortname(genome +'.faa')
            cmds = '%s --kingdom Bacteria --force --outdir %s --locustag Bacter %s\n' % \
                         ('prokka', output_dir,
                          genome + '.prokka.fasta')
    return cmds

cmds = ''
for genome in glob.glob(all_genomes):
    cmds += runprokka(genome)

f1 = open(os.path.join(input_script, 'allprokka.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
f1.write('py37\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n')
f1.write(cmds)
f1.close()
