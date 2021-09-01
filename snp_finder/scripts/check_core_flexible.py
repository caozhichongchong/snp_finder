import os,glob
input_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
sum_file = '%s/summary/all.species.txt'%(input_folder)
core_flexible = '%s/summary/oldspecies/all.denovo.gene.faa.uc.species.sum'%(input_folder)

core_flexible_set = dict()
for lines in open(core_flexible,'r'):
    lines_set = lines.split('\t')
    core_flexible_set.setdefault(lines_set[0],lines_set[3].split('\n')[0])

newoutput = []
for lines in open(sum_file,'r'):
    lines = lines.split('\n')[0]
    if lines.startswith('#donor_species'):
        newoutput.append('%s\tgenus_num\n' % (lines))
    else:
        lines_set = lines.split('\t')
        donor_species = lines_set[0]
        donor_species_new = donor_species.replace('_clustercluster', '_CL').replace('_PB_', '_PaDi_')
        record_id = lines_set[1]
        if 'other' not in record_id and record_id!='0' and 'allspecies' not in record_id:
            #print(record_id)
            record_id = 'C_%s_G_%s' % (record_id.split('_')[1], record_id.split('_')[-1])
            newrecord_id = '%s__%s' % (donor_species_new, record_id)
            if newrecord_id in core_flexible_set:
                newoutput.append('%s\t%s\n'%(lines,core_flexible_set[newrecord_id]))
            else:
                print('%s not in cluster'%(newrecord_id))

f1 = open('%s/summary/all.species.genusnum.txt'%(input_folder), 'w')
f1.write(''.join(newoutput))
f1.close()
