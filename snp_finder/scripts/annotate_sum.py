# start
# gene annotation summary
import os
import glob
import copy
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
# optional output setup
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='snp_output/',
                      metavar='snp_output/')
################################################## Definition ########################################################
args = parser.parse_args()

output_dir_merge = args.o
database_metacyc = '/scratch/users/mit_alm/database/metacyc/genes.col'
database_eggnog = '/scratch/users/mit_alm/database/eggnog/2_annotations.tsv'
database_eggnog2 = '/scratch/users/mit_alm/database/eggnog/COG_catogary.txt'
database_kegg = '/scratch/users/anniz44/scripts/database/Kegg/ko_formated'
all_fasta_set = glob.glob(os.path.join(output_dir_merge + '/summary', '*.cluster.aa'))

# function annotation
def best_hit(Blast_output,small = 0):
    for gene_name in Blast_output:
        db_name,identity = Blast_output[gene_name]
        if len(identity) > 2:
            identity2 = copy.deepcopy(identity)
            if small == 1:
                identity2.sort()
            else:
                identity2.sort(reverse=True)
            top1 = identity2[0]
            top2 = identity2[1]
            Blast_output[gene_name]=[db_name[identity.index(top1)],
                                     db_name[identity.index(top2)]]
        else:
            Blast_output[gene_name]=db_name
    return Blast_output

def annotate_metacyc(blast_search):
    Blast_output = dict()
    DB_name = dict()
    for lines in open(blast_search, 'r'):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            gene_name = lines_set[0]
            db_name = lines_set[1].split('|')[-1].split('-MONOMER')[0]
            identity = float(lines_set[2])
            Blast_output.setdefault(gene_name, [[],[]])
            Blast_output[gene_name][0].append(db_name)
            Blast_output[gene_name][1].append(identity)
            DB_name.setdefault(db_name, ['',''])
    for database in [database_metacyc]:
        for lines in open(database, 'r', encoding="utf8", errors='ignore'):
            if not lines.startswith('#'):
                lines_set = lines.split('\n')[0].split('\t')
                db_name = lines_set[0].split('-MONOMER')[0]
                function_name = lines_set[1]
                annotation_fun = lines_set[2]
                if db_name in DB_name:
                    DB_name[db_name][0] = function_name
                    DB_name[db_name][1] = annotation_fun
    best_hit(Blast_output, 0)
    return [Blast_output,DB_name]

def annotate_eggnog(blast_search):
    Anno_eggnog = dict()
    for lines in open(database_eggnog2,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        Anno_eggnog.setdefault(lines_set[0],lines.split('\n')[0])
    Blast_output = dict()
    DB_name = dict()
    DB_name_2 = dict()
    for blast_search_file in blast_search:
        for lines in open(blast_search_file, 'r'):
            if not lines.startswith('#'):
                db_name = ''
                identity = 0
                lines_set = lines.split('\n')[0].split(' ')
                gene_name = lines_set[0]
                for sub_line in lines_set[1:]:
                    if sub_line !='' and sub_line !='-':
                        if db_name == '':
                            db_name = sub_line.split('.')[0]
                        elif identity == 0:
                            identity = float(sub_line)
                            break
                Blast_output.setdefault(gene_name, [[], []])
                Blast_output[gene_name][0].append(db_name)
                Blast_output[gene_name][1].append(identity)
                DB_name.setdefault(db_name, ['', ''])
                DB_name_2.setdefault(db_name, '')
    for database in [database_eggnog]:
        for lines in open(database, 'r'):
            if not lines.startswith('#'):
                lines_set = lines.split('\n')[0].split('\t')
                db_name = lines_set[1]
                function_name = ''
                annotation_catogary = ''
                for COG in lines_set[2]:
                    annotation_catogary += '\t' + Anno_eggnog[COG]
                annotation_fun = lines_set[3]
                if db_name in DB_name:
                    DB_name[db_name][0] = function_name
                    DB_name[db_name][1] = annotation_fun
                    DB_name_2[db_name] = '%s\t%s%s'%(function_name,annotation_fun,annotation_catogary)
    best_hit(Blast_output, 1)
    return [Blast_output,DB_name,DB_name_2]

def annotate_kegg(blast_search):
    Blast_output = dict()
    DB_name = dict()
    DB_name2 = dict()
    for lines in open(blast_search, 'r'):
        if not lines.startswith('#'):
            db_name = ''
            identity = 0
            lines_set = lines.split('\n')[0].split(' ')
            gene_name = lines_set[0]
            for sub_line in lines_set[1:]:
                if sub_line != '' and sub_line != '-':
                    if db_name == '':
                        db_name = sub_line
                    elif identity == 0:
                        identity = float(sub_line)
                        break
            Blast_output.setdefault(gene_name, [[], []])
            Blast_output[gene_name][0].append(db_name)
            Blast_output[gene_name][1].append(identity)
            DB_name.setdefault(db_name, ['', ''])
            DB_name2.setdefault(db_name, '')
    for database in [database_kegg]:
        for lines in open(database, 'r'):
            if not lines.startswith('#'):
                lines_set = lines.split('\n')[0].split('\t')
                db_name = lines_set[0]
                function_name = lines_set[-6]
                annotation_fun = lines_set[-5]
                if db_name in DB_name:
                    DB_name[db_name][0] = function_name
                    DB_name[db_name][1] = annotation_fun
                    if 'Drug resistance' in lines_set[5] or ('Cellular community - eukaryotes' not in lines_set[5] and 'Overview' not in lines_set[5] \
                            and 'Human Diseases' not in lines_set[4] and 'Organismal Systems' not in lines_set[4]):
                        DB_name2[db_name] = '\t'.join(lines_set[4:])
    best_hit(Blast_output, 1)
    return [Blast_output,DB_name,DB_name2]

def annotate_prokka(fasta_input):
    cutoff = 80
    cutoff2 = 80
    cmds = ("diamond makedb --in %s -d %s.dmnd\ndiamond blastp --query %s --db %s.dmnd --out %s.prokka.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                % (database_prokka_fasta,database_prokka_fasta,
                   fasta_input, database_prokka_fasta, fasta_input, cutoff, cutoff2))
    os.system(cmds)
    blast_search = '%s.prokka.txt' %(fasta_input)
    Blast_output = dict()
    DB_name = dict()
    for lines in open(blast_search, 'r'):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            gene_name = lines_set[0]
            db_name = lines_set[1]
            identity = float(lines_set[2])
            Blast_output.setdefault(gene_name, [[], []])
            Blast_output[gene_name][0].append(db_name)
            Blast_output[gene_name][1].append(identity)
            DB_name.setdefault(db_name, ['', ''])
    for database in [database_prokka]:
        for lines in open(database, 'r'):
            if not lines.startswith('#'):
                lines_set = lines.split('\n')[0].split('\t')
                db_name = lines_set[0]
                function_name = lines_set[3].split('_')[0]
                annotation_fun = lines_set[6]
                if db_name in DB_name:
                    DB_name[db_name][0] = function_name
                    DB_name[db_name][1] = annotation_fun
    best_hit(Blast_output, 1)
    return [Blast_output,DB_name]

def cluster_uc(cluster_input):
    Clusters = dict()
    for lines in open(cluster_input, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        cluster = line_set[1]
        record_name = line_set[8]
        Clusters.setdefault(record_name, cluster)
    return Clusters

def annotate_gene(Blast_output1,DB_name1,All_annotation):
    for gene_name in Blast_output1:
        All_annotation.setdefault(gene_name, [[], set(), set(),[], [], []])
        for db_name in Blast_output1[gene_name]:
            if db_name in DB_name1:
                fun, anno = DB_name1[db_name]
                All_annotation[gene_name][1].add(fun)
                All_annotation[gene_name][2].add(anno)
    return All_annotation

def sum_kegg_eggnog(Blast_outputkegg,DB_namekegg,Blast_outputegg,DB_nameegg,Clusters):
    Output = []
    for gene_name in Blast_outputkegg:
        cluster = Clusters.get(gene_name, '')
        donor = gene_name.split('_')[0]
        donor_species = gene_name.split('__')[0].split('_CL')[0].replace('BA_BA', 'Bifido_adoles').replace('BL_BL',
                                                                                                           'Bifido_longum').replace(
            'PB_PB', 'Paraba_butyra').replace('BA_cluste', 'Bifido_adoles').replace('BL_cluste',
                                                                                    'Bifido_longum').replace(
            'PB_cluste', 'Paraba_butyra')
        species = donor_species.replace(donor + '_', '')
        for db_name in Blast_outputkegg[gene_name]:
            if db_name in DB_namekegg:
                Output.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(cluster,donor,species,donor_species,gene_name,db_name,DB_namekegg[db_name]))
    f1 = open(all_fasta + '.kegg.sum','w')
    f1.write('cluster\tdonor\tspecies\tdonor_species\tgene_name\tKO\tBRITE_KO1\tBRITE_KO2\tBRITE_KO3\n' + ''.join(Output))
    f1.close()
    Output = []
    for gene_name in Blast_outputegg:
        cluster = Clusters.get(gene_name, '')
        donor = gene_name.split('_')[0]
        donor_species = gene_name.split('__')[0].split('_CL')[0].replace('BA_BA', 'Bifido_adoles').replace('BL_BL',
                                                                                                           'Bifido_longum').replace(
            'PB_PB', 'Paraba_butyra').replace('BA_cluste', 'Bifido_adoles').replace('BL_cluste',
                                                                                    'Bifido_longum').replace(
            'PB_cluste', 'Paraba_butyra')
        species = donor_species.replace(donor + '_', '')
        for db_name in Blast_outputegg[gene_name]:
            if db_name in DB_nameegg:
                Output.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                cluster, donor, species, donor_species, gene_name, db_name, DB_nameegg[db_name]))
    f1 = open(all_fasta + '.eggnog.sum', 'w')
    f1.write(
        'cluster\tdonor\tspecies\tdonor_species\tgene_name\tEGGNOG\tannotation\tCOG\tCOG1\tCOG2\n' + ''.join(Output))
    f1.close()

def cluster_genes(All_annotation,Clusters):
    for gene_name in Clusters:
        cluster = Clusters.get(gene_name,'')
        donor = gene_name.split('_')[0]
        donor_species = gene_name.split('__')[0].split('_CL')[0].replace('BA_BA','Bifido_adoles').replace('BL_BL','Bifido_longum').replace('PB_PB','Paraba_butyra').replace('BA_cluste','Bifido_adoles').replace('BL_cluste','Bifido_longum').replace('PB_cluste','Paraba_butyra')
        species = donor_species.replace(donor + '_', '')
        All_annotation.setdefault(gene_name,[[], set(), set(),[], [], []])
        All_annotation[gene_name][3].append(donor)
        All_annotation[gene_name][4].append(species)
        All_annotation[gene_name][5].append(donor_species)
        if cluster!= []:
            All_annotation[gene_name][0].append(cluster)
        else:
            print('missing cluster for gene %s'%(gene_name))
    return All_annotation

def cluster_species(All_annotation):
    Clusters = dict()
    Species = dict()
    Donorspecies = dict()
    Donor = dict()
    output_cluster = []
    Species_species_count = dict()
    output_sum = []
    output_cluster.append('cluster\tTag_species\tTag_donor\tdonor\tspecies\tdonor_species\tfun\tanno\t\n')
    output_sum.append('gene_name\tcluster\tdonor\tspecies\tdonor_species\tfun\tanno\t\n')
    for gene_name in All_annotation:
        clusters, fun, anno, donor, species, donor_species = All_annotation[gene_name]
        All_annotation[gene_name][0] = ';'.join(list(set(clusters)))
        for cluster in clusters:
            Clusters.setdefault(cluster, [[], [], [], set(), set()])
            Clusters[cluster][0].append(donor[0])
            Clusters[cluster][1].append(species[0])
            Clusters[cluster][2].append(donor_species[0])
            All_annotation[gene_name][1] = set()
            All_annotation[gene_name][2] = set()
            for a_fun in fun:
                if a_fun != '':
                    All_annotation[gene_name][1].add(a_fun)
                    Clusters[cluster][3].add(a_fun)
            for a_anno in anno:
                if a_anno != '':
                    All_annotation[gene_name][2].add(a_anno)
                    Clusters[cluster][4].add(a_anno)
        cluster, fun, anno, donor, species, donor_species = All_annotation[gene_name]
        output_sum.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n' % (gene_name,
                                                              cluster, ';'.join(donor), ';'.join(species),
                                                              ';'.join(donor_species), ';'.join(fun), ';'.join(anno)
                                                              ))
    Clusters_sum = dict()
    for cluster in Clusters:
        Tag_species = 'unique_species'
        Tag_donor = 'unique_donor'
        species_set = set()
        genus_set = set()
        alldonor, allspecies, alldonorspecies, fun, anno = Clusters[cluster]
        for donor in alldonor:
            Donor.setdefault(donor, [])
            Donor[donor].append(cluster)
        for donorspecies in alldonorspecies:
            Donorspecies.setdefault(donorspecies, [])
            Donorspecies[donorspecies].append(cluster)
        for species in allspecies:
            new_species = species.split('_cluste')[0].split('_newcluster')[0].split('_CL')[0]
            species_set.add(new_species)
            Species.setdefault(new_species, [])
            Species[new_species].append(cluster)
            new_genus = new_species.split('_')[0]
            genus_set.add(new_genus)
            Species_species_count.setdefault(new_species, [])
        if len(allspecies) > 1:
            if len(alldonor) > 1:
                Tag_donor = 'multi_donor'
            if len(species_set) > 1:
                Tag_species = 'multi_species'
            if len(genus_set) > 1:
                Tag_species = 'multi_genera'
        Clusters_sum.setdefault(cluster, [Tag_species, Tag_donor, ';'.join(list(fun)), ';'.join(list(anno))])
        Tag_species, Tag_donor, fun, anno = Clusters_sum[cluster]
        output_cluster.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n' % (cluster, Tag_species, Tag_donor,
                                                         ';'.join(alldonor), ';'.join(list(species_set)),';'.join(alldonorspecies),
                                                          fun, anno))
        for new_species in species_set:
            for new_species2 in species_set:
                if new_species < new_species2:
                    Species_species_count[new_species].append(new_species2)
                    Species_species_count[new_species2].append(new_species)
    Species_sum = []
    Species_sum.append('#species\tgenus\tcluster\tTag_species\tTag_donor\t\n')
    for new_species in Species:
        new_genus = new_species.split('_')[0]
        clusters = list(set(Species[new_species]))
        for cluster in clusters:
            Species_sum.append('%s\t%s\t%s\t%s\t%s\t\n' % (
                new_species, new_genus, cluster, Clusters_sum[cluster][0], Clusters_sum[cluster][1]))
    Donor_sum = []
    Donor_sum.append('#donor\tcluster\tTag_species\tTag_donor\t\n')
    for donor in Donor:
        clusters = list(set(Donor[donor]))
        for cluster in clusters:
            Donor_sum.append(
                '%s\t%s\t%s\t%s\t\n' % (donor, cluster, Clusters_sum[cluster][0], Clusters_sum[cluster][1]))
    Donorspecies_sum = []
    Donorspecies_sum.append('#donorspecies\tcluster\tTag_species\tTag_donor\t\n')
    for donorspecies in Donorspecies:
        clusters = list(set(Donorspecies[donorspecies]))
        for cluster in clusters:
            Donorspecies_sum.append(
                '%s\t%s\t%s\t%s\t\n' % (donorspecies, cluster, Clusters_sum[cluster][0], Clusters_sum[cluster][1]))
    f1 = open(all_fasta + '.gene.sum', 'w')
    f1.write(''.join(output_sum))
    f1.close()
    f1 = open(all_fasta + '.cluster.sum', 'w')
    f1.write(''.join(output_cluster))
    f1.close()
    f1 = open(all_fasta + '.species.sum', 'w')
    f1.write(''.join(Species_sum))
    f1.close()
    f1 = open(all_fasta + '.donor.sum', 'w')
    f1.write(''.join(Donor_sum))
    f1.close()
    f1 = open(all_fasta + '.donor.species.sum', 'w')
    f1.write(''.join(Donorspecies_sum))
    f1.close()
    Output = []
    Output.append('#species1\tspecies2\tgenus1\tgenus2\tnum_cluster_shared\t\n')
    for new_species in Species_species_count:
        all_species = Species_species_count[new_species]
        all_species_unique = set(all_species)
        new_genus1 = new_species.split('_')[0]
        for new_species2 in all_species_unique:
            new_genus2 = new_species2.split('_')[0]
            Output.append(
                '%s\t%s\t%s\t%s\t%s\t\n' % (new_species, new_species2, new_genus1,new_genus2, all_species.count(new_species2)))
    f1 = open(all_fasta + '.speciesTospecies.sum', 'w')
    f1.write(''.join(Output))
    f1.close()

def sum_annotation(Blast_output1,DB_name1,Blast_output2,DB_name2,Blast_output3,DB_name3,Blast_output4,DB_name4,Clusters_gene):
    All_annotation = dict()
    All_annotation = annotate_gene(Blast_output1, DB_name1, All_annotation)
    All_annotation = annotate_gene(Blast_output2, DB_name2, All_annotation)
    All_annotation = annotate_gene(Blast_output3, DB_name3, All_annotation)
    All_annotation = annotate_gene(Blast_output4, DB_name4, All_annotation)
    All_annotation = cluster_genes(All_annotation, Clusters_gene)
    cluster_species(All_annotation)

for all_fasta in all_fasta_set:
    all_fasta_folder, all_fasta_file = os.path.split(all_fasta)
    try:
        database_prokka = glob.glob(os.path.join(all_fasta_folder, 'prokka_%s/*.tsv' % (all_fasta_file)))[0]
        database_prokka_fasta = glob.glob(os.path.join(all_fasta_folder, 'prokka_%s/*.faa' % (all_fasta_file)))[0]
        Blast_output4, DB_name4 = annotate_prokka(all_fasta)
    except IndexError:
        Blast_output4 = dict()
        DB_name4 = dict()
    # set up gene annotation and clustering
    Blast_output1, DB_name1, DB_name1_2 = annotate_kegg(all_fasta + '.kegg.txt')
    Blast_output2, DB_name2 = annotate_metacyc(all_fasta + '.metacyc.txt')
    Blast_output3, DB_name3, DB_name3_2 = annotate_eggnog(glob.glob(all_fasta + '.eggnog.*.txt'))
    # sum up
    Clusters_gene = cluster_uc(all_fasta.replace('.cluster.aa','.uc'))
    sum_kegg_eggnog(Blast_output1,DB_name1_2,Blast_output3,DB_name3_2,Clusters_gene)
    sum_annotation(Blast_output1, DB_name1, Blast_output2, DB_name2, Blast_output3, DB_name3, Blast_output4, DB_name4,
                   Clusters_gene)
################################################### END ########################################################
