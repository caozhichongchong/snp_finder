import glob,os
from Bio import SeqIO

input_results = glob.glob('/scratch/users/anniz44/genomes/donor_species/SNP_curate/xq_data/merge/*.allfreqsum.txt')
input_reference = glob.glob('/scratch/users/xy43/mapper_testing/reference_genomes/*.fa')
allfastalength = dict()
for fasta in input_reference:
    totallength = 0
    genome = os.path.split(fasta)[-1].split('.fa')[0]
    for record in SeqIO.parse(fasta, 'fasta'):
        totallength += len(str(record.seq))
    allfastalength.setdefault(genome,totallength)

allresults = []
allresults.append('test\tmethod\tgenome\tsumdepth\tmean_depth\tdepth_0.1\tmedian_depth\tdepth_0.9\tgenome_coverage\tgenome_length\n')
for files in input_results:
    filesname = os.path.split(files)[-1].split('.')[0]
    method = os.path.split(files)[-1].split('.')[1]
    for lines in open(files):
        if not lines.startswith('genome'):
            allresults.append('%s\t%s\t%s\t%s\n'%(filesname,method,lines.split('\n')[0],allfastalength[lines.split('\t')[0]]))

f1 = open('/scratch/users/anniz44/genomes/donor_species/SNP_curate/xq_data/merge/allfreqsum.txt', 'w')
f1.write(''.join(allresults))
f1.close()
