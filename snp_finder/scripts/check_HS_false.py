import os,glob
input_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
sum_file = '%s/summary/all.species.txt'%(input_folder)
Min_SNP_highselect_cutoff = 1/3000

HS_false = dict()
for lines in open(sum_file,'r'):
    if 'NODE_' in lines and 'False' in lines:
        lines_set = lines.split('\t')
        if int(lines_set[4])>=2 and int(lines_set[9]) == 0 and int(lines_set[5])/int(lines_set[3])>=Min_SNP_highselect_cutoff:
            HS_false.setdefault(lines_set[0],set())
            HS_false[lines_set[0]].add(lines_set[1])

alloutput = []
for donor_species in HS_false:
    print(donor_species)
    vcf_file = '%s/moredetails/%s.raw.vcf.filtered.vcf.noremoverec.snpfreq.txt'%(
        input_folder,donor_species.replace('.donor','.all.donor'))
    allgene = HS_false[donor_species]
    oldline = ''
    oldline2 = ''
    nextline = False
    for lines in open(vcf_file,'r'):
        lines_set = lines.split('\t')
        if lines_set[7] in allgene:
            alloutput.append('%s\t%s\t%s' % (donor_species, lines_set[7], oldline2))
            alloutput.append('%s\t%s\t%s' % (donor_species, lines_set[7], oldline))
            alloutput.append('%s\t%s\t%s' % (donor_species, lines_set[7], lines))
            nextline = True
            targetgene = lines_set[7]
        elif nextline == True and targetgene !='':
            alloutput.append('%s\t%s\t%s' % (donor_species, targetgene, lines))
            nextline = False
            targetgene = ''
        oldline2 = oldline
        oldline = lines

f1 = open('%s/summary/all.HS.false.txt'%(input_folder), 'w')
f1.write(''.join(alloutput))
f1.close()
