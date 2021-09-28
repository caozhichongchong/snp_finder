import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'

input_coassembly = '%s/co-assembly'%(output_folder)
input_bs_file = '%s/binding_results_ccpA.txt'%(output_folder)
motif_file = '%s/motifs.meme'%(output_folder)
script_file = '%s/fimo/'%(input_script)
script_file2 = '%s/predictaa/'%(input_script)
length_halfBS = 7

os.system('rm -r %s'%(script_file))
try:
    os.mkdir(script_file)
except IOError:
    pass
os.system('rm -r %s'%(script_file2))
try:
    os.mkdir(script_file2)
except IOError:
    pass

def find_BS(input_fasta):
    allBS = []
    allBS.append('BS\tcontig\tlocus\n')
    for record in SeqIO.parse(input_fasta, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        for i in range(0, (len(record_seq) - length_halfBS * 2)):
            seq1 = record_seq[i: (i + length_halfBS)]
            seq2 = record_seq[(i + length_halfBS):(i + 2 * length_halfBS)]
            if seq1 == str(Seq(seq2).reverse_complement()) and seq1 != 'NNNNNNN':
                allBS.append('%s%s\t%s\t%s\n' % (seq1, seq2,
                                                              record_id, i))
    f1 = open('%s/%s.BS.txt'%(output_folder, genomename), 'w')
    f1.write(''.join(allBS))
    f1.close()

for fasta in glob.glob('%s/*.fasta'%(input_coassembly)):
    genomename = os.path.split(fasta)[-1].split('.all.spades2.fasta')[0]
    cmds = '#!/bin/bash\nsource ~/.bashrc\n'
    cmds += 'fimo -o %s/%s %s %s\n' % (output_folder, genomename, motif_file, fasta)
    cmds += 'mv %s/%s/fimo.tsv %s/%s.fimo.tsv\n'%(output_folder, genomename,output_folder, genomename)
    cmds += 'rm -r %s/%s\n' % (output_folder, genomename)
    f1 = open('%s/%s.sh' % (script_file, genomename), 'w')
    f1.write(cmds)
    f1.close()
    cmds = '#!/bin/bash\nsource ~/.bashrc\n'
    cmds += 'prodigal -q -i %s -a %s.faa\n' % (
        fasta, fasta)
    f1 = open('%s/%s.sh' % (script_file2, genomename), 'w')
    f1.write(cmds)
    f1.close()
    try:
        f1 = open('%s/%s.BS.txt' % (output_folder, genomename), 'r')
    except FileNotFoundError:
        find_BS(fasta)

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