import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details'
genus_cutoff = 3 # >= 3 genera as core
allgenome = output_dir_merge + '/summary/all.genome.gene.faa'
allgenome_denovo = output_dir_merge + '/summary/all.denovo.gene.faa'
allgenome_HS = output_dir_merge + '/summary/all.selected.gene.faa'

def genus_list(species_list):
    return set([i.replace('BL','BiLo').replace('BA','BiAd')[0:2] for i in species_list])

def function_species(fastaname):
    species_file = fastaname + '.uc.species.sum'
    eggnog_file = fastaname + '.cluster.aa.all.eggnog.sum'
    kegg_file = fastaname + '.cluster.aa.all.kegg.sum'
    Species_fun = dict()
    for lines in open(species_file, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        record_name = line_set[0]
        species_num = line_set[1]
        species_name = line_set[2].split(',')
        tag_genus = 'False'
        if len(species_name) > 1:
            genus = genus_list(species_name)
            if len(genus) >= genus_cutoff:
                tag_genus = 'True'
        Species_fun.setdefault(record_name, [species_num,tag_genus])
    Output = []
    for lines in open(eggnog_file, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        if lines.startswith('cluster'):
            Output.append('\t'.join(line_set) + '\tspecies_num\tcross_genus\n')
        else:
            if line_set[6] == '':
                line_set.pop(6)
            record_name = line_set[4]
            if record_name in Species_fun:
                species_num, tag_genus = Species_fun[record_name]
            else:
                species_num, tag_genus = [1,'False']
            Output.append('\t'.join(line_set) + '\t%s\t%s\n'%(species_num, tag_genus))
    f1 = open(fastaname + '.cluster.aa.all.eggnog.sum.species.sum', 'w')
    f1.write(''.join(Output))
    f1.close()
    Output = []
    for lines in open(kegg_file, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        if lines.startswith('cluster'):
            Output.append('\t'.join(line_set) + '\tspecies_num\tcross_genus\n')
        else:
            if line_set[6] == '':
                line_set.pop(6)
            record_name = line_set[4]
            if record_name in Species_fun:
                species_num, tag_genus = Species_fun[record_name]
            else:
                species_num, tag_genus = [1,'False']
            Output.append('\t'.join(line_set) + '\t%s\t%s\n'%(species_num, tag_genus))
    f1 = open(fastaname + '.cluster.aa.all.kegg.sum.species.sum', 'w')
    f1.write(''.join(Output))
    f1.close()

function_species(allgenome_HS)
function_species(allgenome_denovo)
function_species(allgenome)

# split HS annotation within lineage and across lineages
def HS_lineage(filename,HS_lineagefasta):
    HS_lineage_set = set()
    for record in SeqIO.parse(HS_lineagefasta, 'fasta'):
        record_id = str(record.id)
        HS_lineage_set.add(record_id)
    Output_within = set()
    Output_across = set()
    for lines in open(filename,'r'):
        if lines.split('\t')[4] in HS_lineage_set:
            Output_within.add(lines)
        else:
            Output_across.add(lines)
    f1 = open(filename + '.within.sum', 'w')
    f1.write(''.join(list(Output_within)))
    f1.close()
    f1 = open(filename + '.across.sum', 'w')
    f1.write(''.join(list(Output_across)))
    f1.close()

#HS_annotation_sum = output_dir_merge + '/summary/all.selected.gene.faa.cluster.aa.all.eggnog.sum.species.sum'
#allgenome_HS_lineage = output_dir_merge + '/summary/all.selected.gene.faa'

#HS_lineage(HS_annotation_sum,allgenome_HS_lineage)
