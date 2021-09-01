import os,glob
allsumfile = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/*.sum.txt')

allout = []
for sumfile in allsumfile:
    total_snp = 0
    snp_type = dict()
    donor_species = os.path.split(sumfile)[-1].split('.all.parsi.fasta.sum.txt ')[0]
    for lines in open(sumfile,'r'):
        if not lines.startswith('CHR'):
            snptype = lines.split('\t')[8]
            snp_type.setdefault(snptype,[0,0,0])
            total_snp += 1
            if lines.split('\t')[4] == 'N':
                snp_type[snptype][0] += 1
            elif lines.split('\t')[4] == 'S':
                snp_type[snptype][1] += 1
            else:
                snp_type[snptype][2] += 1
    for snptype in snp_type:
        allout.append('%s\t%s\t%s\t%s\t%s\t%s\n'%(donor_species,snptype,snp_type[snptype][0],
                                          snp_type[snptype][1],snp_type[snptype][2],total_snp))

f1 = open('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/all.genotype.sum.txt','w')
f1.write('donor_species\tsnptype\tN\tS\tothers\ttotal_SNP\n')
f1.write(''.join(allout))
f1.close()
