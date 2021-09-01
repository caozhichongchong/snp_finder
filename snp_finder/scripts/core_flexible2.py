import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
cluster_cutoff = 0.9
allgenomeref = output_dir_merge + '/summary/all.genome.gene.faa'
allgenome = output_dir_merge + '/summary/all.genome.gene.faa.uc'
allgenome_denovo = output_dir_merge + '/summary/all.denovo.gene.faa.uc'
allgenome_HS = output_dir_merge + '/summary/all.selected.gene.faa.uc'

def cluster_uc(cluster_input):
    # run cluster
    try:
        f1 = open(cluster_input,'r')
    except IOError:
        cutoff = cluster_cutoff
        cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -uc %s -threads %s\n'
                       % ('usearch', allgenomeref, cutoff,
                          cluster_input, 40))
        os.system(cmd_cluster)
    Clusters = dict()
    for lines in open(cluster_input, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        cluster = line_set[1]
        record_name = line_set[8].split(' ')[0]
        Clusters.setdefault(cluster, set())
        Clusters[cluster].add(record_name)
    return Clusters

# clusters based on all genes
Clusters = cluster_uc(allgenome)

def genus_list(species_list):
    return set([i.replace('BL','BiLo').replace('BA','BiAd')[0:2] for i in species_list])

def species_cluster(Clusters):
    Clusters_species = dict()
    Gene_cluster = dict()
    Output = []
    Output.append('record\tspecies_num\tallspecies\tgenus_num\t\n')
    for cluster in Clusters:
        allrecord = Clusters[cluster]
        Clusters_species.setdefault(cluster,set())
        for record in allrecord:
            species = record.split('_')[0]
            if species == '1':
                species = record.split('_')[1]
            species = species.replace('PB','PaDi').replace('Bfragilis','BaFr')
            Clusters_species[cluster].add(species)
        for record in allrecord:
            Gene_cluster.setdefault(record,Clusters_species[cluster])
    for record in Gene_cluster:
        allspecies = Gene_cluster[record]
        genus = genus_list(allspecies)
        Output.append('%s\t%s\t%s\t%s\t\n'%(record,len(allspecies),','.join(list(allspecies)),len(genus)
                                        ))
    #f1 = open(allgenome + '.species.sum','w')
    #f1.write(''.join(Output))
    #f1.close()
    return [Clusters_species,Gene_cluster]

Clusters_species,Gene_cluster = species_cluster(Clusters)

def cluster_species_sum(fastaname,Gene_cluster):
    Clusters = cluster_uc(fastaname)
    Output = []
    Output.append('record\tspecies_num\tallspecies\tgenus_num\t\n')
    for cluster in Clusters:
        allrecord = Clusters[cluster]
        for record in allrecord:
            if record in Gene_cluster:
                record_new = record
            else:
                record_new = record.replace('.donor.'+record.split('.')[-1].split('__')[0],'')
                if record_new not in Gene_cluster:
                    record_new = record.split('_')[0] + '__' + record.split('__')[1]
            if '1_PaDi_IBD' in record_new:
                record_new = record_new.replace('1_PaDi_IBD','1_PB_IBD')
            if record_new in Gene_cluster:
                allspecies = Gene_cluster[record_new]
                genus = genus_list(allspecies)
                Output.append('%s\t%s\t%s\t%s\t\n' % (
                    record, len(allspecies), ','.join(list(allspecies)),len(genus)))
            else:
                print('missing genome uc for %s %s'%(record,record_new))
    f1 = open(fastaname + '.species.sum', 'w')
    f1.write(''.join(Output))
    f1.close()

cluster_species_sum(allgenome_denovo,Gene_cluster)
cluster_species_sum(allgenome_HS,Gene_cluster)
