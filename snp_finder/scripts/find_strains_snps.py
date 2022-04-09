import os,glob

metadata = '/scratch/users/anniz44/genomes/donor_species/selected_species/file.change.name.new.txt.addtime.addclonal.txt.newname.txt'
vcf_dir = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/final'
snpmeta = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/final/alltargetsnp.txt'
def load_meta():
    genomenewname = dict()
    for lines in open(metadata,'r'):
        lines_set = lines.split('\t')
        oldname = lines_set[1]
        newname = lines_set[2]
        genomenewname.setdefault(newname,oldname)
    return genomenewname

def load_snp():
    snpstrains = dict()
    for lines in open(snpmeta,'r'):
        if not lines.startswith('SNP'):
            lines_set = lines.split('\n')[0].split('\t')
            lineage = lines_set[1].replace('CL','clustercluster').split('__')[0] + '.all.parsi.fasta.linktrunc.sum.txt'
            wanted_strain = [int(i) for i in lines_set[-1].split(';')]
            unwanted_strain = [int(i) for i in lines_set[-2].split(';')]
            lines_set_short = '\t'.join(lines_set[0:4]) + '\t' + lines_set[-1]
            snpstrains[lineage] = snpstrains.setdefault(lineage, [])
            snpstrains[lineage].append([wanted_strain,unwanted_strain,lines_set_short])
    return snpstrains

def findstrains(snpstrains):
    alloutput = []
    alloutput.append('SNP\tGenename\tDonor\tSNP_region\tGenome_SNP_set\tMT\tWT\n')
    for lineage in snpstrains:
        for lines in open('%s/%s'%(vcf_dir,lineage), 'r'):
            if lines.startswith('CHR'):
                lines_set = lines.split('\n')[0].split('\t')[9:]
                for wanted_strain,unwanted_strain,lines_set_short in snpstrains[lineage]:
                    wanted_strain_set = [genomenewname.get(lines_set[i-1],lines_set[i-1]) for i in wanted_strain]
                    unwanted_strain_set = [genomenewname.get(lines_set[i - 1], lines_set[i - 1]) for i in unwanted_strain]
                    alloutput.append('%s\t%s\t%s\n'%(lines_set_short,';'.join(wanted_strain_set),';'.join(unwanted_strain_set)))
            break
    f1 = open('%s/alltargetsnp.strains.txt'%(vcf_dir),'w')
    f1.write(''.join(alloutput))
    f1.close()

genomenewname = load_meta()
snpstrains = load_snp()
print(snpstrains)
findstrains(snpstrains)
