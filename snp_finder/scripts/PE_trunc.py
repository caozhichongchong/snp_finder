import os,glob
snp_dir = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/patients/'

def loadHS(sum_file):
    HS_gene = dict()
    for lines in open(sum_file, 'r'):
        if not lines.startswith("#donor_species"):
            lines_set = lines.split('\n')[0].split('\t')
            if lines_set[-1] == 'True':
                lineage, gene = lines_set[0:2]
                HS_gene.setdefault(lineage,set())
                HS_gene[lineage].add(gene)
    return HS_gene

def load_snp(snpfile,alloutput):
    donor_species = os.path.split(snpfile)[-1].split('.raw.vcf.filtered.vcf.final.snp.txt')[0].replace('.all','')
    #if donor_species in HS_gene_within:
    for lines in open(snpfile,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            genename, POS,N_or_S,AAchange = lines_set[-4:]
            genenamenew, POS = lines_set[:2]
            if '*' in AAchange:
                if '*' == AAchange[0]:
                    Trunc_SNP = lines_set[2]
                else:
                    Trunc_SNP = lines_set[3]
                if donor_species in HS_gene_within and genename in HS_gene_within[donor_species]:
                    alloutput.add('%s\t%s\t%s\t%s\t%s\twithinHS\t%s\t%s\n'%(donor_species,genenamenew,POS,AAchange[0],AAchange[1],Trunc_SNP,genename))
                #elif genename in HS_gene_all[donor_species]:
                    #alloutput.add('%s\t%s\t%s\t%s\t%s\tacrossHS\t%s\t%s\n'%(donor_species,genenamenew,POS,AAchange[0],AAchange[1],Trunc_SNP,genename))
                else:
                    alloutput.add('%s\t%s\t%s\t%s\t%s\tothers\t%s\t%s\n' % (
                    donor_species, genenamenew, POS, AAchange[0], AAchange[1],Trunc_SNP,genename))
    else:
        print('no HS in %s'%(donor_species))
    return alloutput

# load HS genes
HS_gene_within = loadHS('%s/summary/all.species.txt'%(snp_dir))
#HS_gene_all = loadHS('%s/summary/all.species.txt.High_select2.txt'%(snp_dir))

alloutput=set()
snp_folder = glob.glob('%s/*.donor.*.raw.vcf.filtered.vcf.final.snp.txt' % (snp_dir))
print(snp_folder)
for snpfile in snp_folder:
    alloutput = load_snp(snpfile,alloutput)

f1 = open('%s/summary/all.species.txt.Truncated.sum'%(snp_dir), 'w')
f1.write(''.join(list(alloutput)))
f1.close()
