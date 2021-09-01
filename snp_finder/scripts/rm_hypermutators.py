import os,glob
input_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
sum_file = '%s/all.genotype.sum.txt'%(input_folder)
allsumfile = glob.glob('%s/*.sum.txt'%(input_folder))

total_SNP_cutoff = 24
ratio_cutoff = 0.2
try:
    os.mkdir(input_folder+'/withhypermutator/')
except IOError:
    pass

snp_type = dict()
for lines in open(sum_file,'r'):
    if not lines.startswith('donor_species'):
        donor_species,snptype,N,S,others,total_SNP =lines.split('\n')[0].split('\t')
        total_genotype = int(N)+int(S)+int(others)
        total_SNP = int(total_SNP)
        if total_SNP > total_SNP_cutoff and total_genotype/total_SNP > ratio_cutoff:
            snp_type.setdefault(donor_species,set())
            snp_type[donor_species].add(snptype)

for sumfile in allsumfile:
    donor_species = os.path.split(sumfile)[-1].split('.all.parsi.fasta.sum.txt ')[0]
    if donor_species in snp_type:
        remove_snptype = snp_type[donor_species]
        print(donor_species,remove_snptype)
        newoutput = []
        for lines in open(sumfile,'r'):
            if not lines.startswith('CHR'):
                snptype = lines.split('\t')[8]
                if snptype not in remove_snptype:
                    newoutput.append(lines)
            else:
                newoutput.append(lines)
        os.system('mv %s %s/withhypermutator/' % (sumfile, input_folder))
        f1 = open(sumfile, 'w')
        f1.write(''.join(newoutput))
        f1.close()
