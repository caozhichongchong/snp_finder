import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

allsnpfiles = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/*.all.parsi.fasta.sum.txt')
annotation_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.denovo.gene.faa.cluster.aa.all.eggnog.sum.species.sum'
changename = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.denovo.gene.faa.changename.txt'

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
