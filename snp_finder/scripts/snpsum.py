################################################### SET PATH ########################################################
# snpsum.py
# merge annotation and snpfiles IBD
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

allsnpfiles = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/*.all.parsi.fasta.sum.txt')
annotation_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.denovo.gene.faa.cluster.aa.all.eggnog.sum'
annotation_file_prokka = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.denovo.gene.faa.prokka.txt'
annotation_file_kegg = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.denovo.gene.faa.cluster.aa.all.kegg.sum'
changename = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.denovo.gene.faa.changename.txt'

genename = dict()
for lines in open(changename,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    donorspecies,nouse,oldgene,newgene = lines_set[0:4]
    genename.setdefault('%s\t%s'%(donorspecies,oldgene),newgene)

annotation_set = dict()
for lines in open(annotation_file,'r'):
    if not lines.startswith('cluster'):
        lines_set = lines.split('\n')[0].split('\t')
        newgene,eggnog,nouse,anno,cog,cog1,cog2 = lines_set[4:11]
        annotation_set.setdefault(newgene,'%s\t%s\t%s\t%s'%(eggnog,anno[:40],cog1,cog2))

annotation_set2 = dict()
for lines in open(annotation_file_prokka,'r'):
    if not lines.startswith('cluster'):
        lines_set = lines.split('\n')[0].split('\t')
        newgene,gene_length,short_genename,gene,COG,product = lines_set
        annotation_set2.setdefault(newgene,'%s\t%s\t%s\t%s'%(gene_length,gene,COG,product))

annotation_set3 = dict()
for lines in open(annotation_file_kegg,'r'):
    if not lines.startswith('cluster'):
        lines_set = lines.split('\n')[0].split('\t')
        try:
            newgene, KO, KO1, KO2, KO3 = lines_set[4:9]
        except ValueError:
            newgene, KO = lines_set[4:6]
            KO1 = ''
            KO2 = ''
            KO3 = ''
        annotation_set3.setdefault(newgene,'%s\t%s\t%s\t%s'%(KO, KO1, KO2, KO3))

for snpfile in allsnpfiles:
    alloutput = []
    donorspecies = os.path.split(snpfile)[-1].split('.all.parsi.fasta.sum.txt')[0]
    donorspecies = donorspecies.replace('_PB_','_PaDi_')
    POS_all_0 = 0
    oldCHR = ''
    vcf_file = '%s/details/%s.raw.vcf.filtered.vcf.final.removerec.snp.txt'%(os.path.split(snpfile)[0],donorspecies.replace('.donor','.all.donor').replace('_PaDi_','_PB_'))
    gene_snp = dict()
    for lines in open(vcf_file,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        gene, genepos, N_S, SNP = lines_set[-4:]
        if SNP == 'None' or SNP == '':
            gene_snp.setdefault('%s_%s' % (gene, genepos), 'None')
        else:
            gene_snp.setdefault('%s_%s'%(gene, genepos),'%s%s%s'%(SNP[0],genepos,SNP[1]))
    for lines in open(snpfile, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        oldgene = lines_set[5]
        if oldgene == 'Gene':
            lines = 'Newgene\tSNP\tGene_length\tProkka_gene\tProkka_COG\tProkka_product\tEggnog\tAnno\tCOG1\tCOG2\tKO\tKO1\tKO2\tKO3\tPOSnew\t' + lines
        else:
            CHR = lines_set[0]
            POS = int(lines_set[1])
            gene,genepos = lines_set[5:7]
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
            annoall = annotation_set2.get(newgene, 'None\tNone\tNone\tNone')
            annoall += '\t' + annotation_set.get(newgene, 'None\tNone\tNone\tNone')
            annoall += '\t' + annotation_set3.get(newgene, 'None\tNone\tNone\tNone')
            lines = '%s\t%s\t%s\t%s\t'%(newgene,gene_snp['%s_%s'%(gene,genepos)],annoall,POS+POS_all_0) + lines
        alloutput.append(lines)
    f1 = open(snpfile.replace('.all.parsi.fasta.sum.txt','.snpsum.txt'), 'w')
    f1.write(''.join(alloutput))
    f1.close()

os.system('mkdir  /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/snpsummary/')
os.system('mv /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/*snpsum.txt /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/snpsummary/')
