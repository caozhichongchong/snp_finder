import os
import glob
import copy
from Bio.Seq import Seq

database_metacyc = '/scratch/users/mit_alm/database/metacyc/genes.col'
database_eggnog = '/scratch/users/mit_alm/database/eggnog/2_annotations.tsv'
database_eggnog2 = '/scratch/users/mit_alm/database/eggnog/COG_catogary.txt'
database_kegg = '/scratch/users/anniz44/scripts/database/Kegg/ko_formated'
#output_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS/Motif_by_DNArec/'
output_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS/'
all_fasta_set = glob.glob(os.path.join(output_file, '*.cluster.aa'))

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
            # check the top 4 hits, keep COG annotation
            try:
                if not db_name[identity.index(top1)].startswith('COG'):
                    top1 = identity2[1]
                    if not db_name[identity.index(top1)].startswith('COG'):
                        top1 = identity2[2]
                        if not db_name[identity.index(top1)].startswith('COG'):
                            top1 = identity2[3]
                            if not db_name[identity.index(top1)].startswith('COG'):
                                top1 = identity2[0]
            except IndexError:
                pass
            Blast_output[gene_name]=[db_name[identity.index(top1)]#,
                                     #db_name[identity.index(top2)]
             ]
        else:
            Blast_output[gene_name]=db_name
    return Blast_output

def annotate_eggnog(blast_search):
    Anno_eggnog = dict()
    for lines in open(database_eggnog2,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        Anno_eggnog.setdefault(lines_set[0],lines.split('\n')[0])
    Blast_output = dict()
    DB_name = dict()
    DB_name_2 = dict()
    try:
        f1 = open(all_fasta + '.eggnog.all.txt','r')
    except IOError:
        os.system('cat %s > %s'%(' '.join(blast_search),all_fasta + '.eggnog.all.txt'))
    blast_search_file = all_fasta + '.eggnog.all.txt'
    for lines in open(blast_search_file, 'r'):
        if not lines.startswith('#'):
            db_name = ''
            identity = 0
            lines_set = lines.split('\n')[0].split(' ')
            gene_name = lines_set[0]
            for sub_line in lines_set[1:]:
                if sub_line != '' and sub_line != '-':
                    if db_name == '':
                        db_name = sub_line.split('.')[0]
                    elif identity == 0:
                        identity = float(sub_line)
                        break
            Blast_output.setdefault(gene_name, [[], []])
            Blast_output[gene_name][0].append(db_name)
            Blast_output[gene_name][1].append(identity)
            DB_name.setdefault(db_name, ['', ''])
            DB_name_2.setdefault(db_name, [])
    for database in [database_eggnog]:
        for lines in open(database, 'r'):
            if not lines.startswith('#'):
                lines_set = lines.split('\n')[0].split('\t')
                db_name = lines_set[1]
                function_name = ''
                annotation_catogary = []
                for COG in lines_set[2]:
                    annotation_catogary.append(Anno_eggnog[COG])
                annotation_fun = lines_set[3]
                if db_name in DB_name:
                    DB_name[db_name][0] = function_name
                    DB_name[db_name][1] = annotation_fun
                    for annotation_catogary_sub in annotation_catogary:
                        DB_name_2[db_name].append('%s%s\t%s'%(function_name,annotation_fun,annotation_catogary_sub))
    best_hit(Blast_output, 1)
    return [Blast_output,DB_name,DB_name_2]

def cluster_uc(cluster_input):
    Clusters = dict()
    for lines in open(cluster_input, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        cluster = line_set[1]
        record_name = line_set[8].split(' ')[0]
        Clusters.setdefault(record_name, cluster)
    return Clusters

def sum_kegg_eggnog(Blast_outputegg,DB_nameegg,Clusters):
    Clusters_cluster = dict()
    for gene in Clusters:
        cluster = Clusters[gene]
        Clusters_cluster.setdefault(cluster,set())
        Clusters_cluster[cluster].add(gene)
    Output = set()
    BS_anno = dict()
    for gene_name in Blast_outputegg:
        cluster = Clusters.get(gene_name, '')
        all_gene = Clusters_cluster.get(cluster,'None')
        if all_gene != 'None':
            for gene_name_sub in all_gene:
                donor = gene_name_sub.split('_')[0]
                donor_species = gene_name_sub.split('__')[0].split('_CL')[0].replace('BA_BA', 'Bifido_adoles').replace('BL_BL',
                                                                                                                   'Bifido_longum').replace(
                    'PB_PB', 'Paraba_butyra').replace('BA_cluste', 'Bifido_adoles').replace('BL_cluste',
                                                                                            'Bifido_longum').replace(
                    'PB_cluste', 'Paraba_butyra')
                species = donor_species.replace(donor + '_', '')
                for db_name in Blast_outputegg[gene_name]:
                    if db_name in DB_nameegg:
                        for funegg in DB_nameegg[db_name]:
                            Output.add('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                            cluster, donor, species, donor_species, gene_name_sub, db_name, funegg))
                            BS_anno.setdefault(donor,'%s\t%s\t%s\t'%(gene_name_sub,db_name,funegg))
    f1 = open(all_fasta + '.all.eggnog.sum', 'w')
    f1.write(
        'cluster\tdonor\tspecies\tdonor_species\tgene_name\tEGGNOG\tannotation\tCOG\tCOG1\tCOG2\n' + ''.join(
            list(Output)))
    f1.close()
    # filter for BS with changes on BS itself
    BS_set = set()
    for lines in open(BS_file,'r'):
        if not lines.startswith('BS'):
            lines_set = lines.split('\t')
            BS = lines_set[0]
            POS_on_BS = int(lines_set[4])
            if POS_on_BS >= -20 and POS_on_BS <= 20:
                # BS with changes on BS itself
                BS_set.add(BS)
    BS_file_out = []
    for lines in open(BS_file,'r'):
        if lines.startswith('BS'):
            lines = 'gene_name\tEGGNOG\tannotation\tCOG\tCOG1\tCOG2\t' + lines
            BS_file_out.append(lines)
        elif lines != '':
            seq = lines.split('\t')[0]
            if seq in BS_set:
                # BS with changes on BS itself
                seq_rc = str(Seq(seq).reverse_complement())
                if seq_rc in BS_anno:
                    seq = seq_rc
                lines = BS_anno.get(seq,'\t\t\t\t\t\t') + lines
                BS_file_out.append(lines)
    f1 = open(BS_file + '.addano.txt', 'w')
    f1.write(''.join(BS_file_out))
    f1.close()

for all_fasta in all_fasta_set:
    all_fasta_folder, all_fasta_file = os.path.split(all_fasta)
    BS_file = all_fasta.replace('.BS.faa.cluster.aa','.BSsum.filtered.txt')
    print(BS_file)
    # sum up
    Clusters_gene = cluster_uc(all_fasta.replace('.cluster.aa', '.uc'))
    Blast_output3, DB_name3, DB_name3_2 = annotate_eggnog(glob.glob(all_fasta + '.eggnog.*.txt'))
    sum_kegg_eggnog(Blast_output3, DB_name3_2, Clusters_gene)
