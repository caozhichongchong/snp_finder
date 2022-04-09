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
optional.add_argument("-i",
                      help="a folder to store all kegg hmmer results",
                      type=str, default='kegg_output/',
                      metavar='kegg_output/')
################################################## Definition ########################################################
args = parser.parse_args()

kegg_folder = args.i
database_kegg = '/scratch/users/anniz44/scripts/database/Kegg/ko_formated'
allkegg_results= glob.glob(os.path.join(kegg_folder, '*.txt'))

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
            #top2 = identity2[1]
            Blast_output[gene_name]=[db_name[identity.index(top1)]#,
                                     #db_name[identity.index(top2)]
             ]
        else:
            Blast_output[gene_name]=db_name
    return Blast_output

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

def sum_kegg(Blast_outputkegg,DB_namekegg):
    Output = set()
    for gene_name in Blast_outputkegg:
        for db_name in Blast_outputkegg[gene_name]:
            if db_name in DB_namekegg:
                Output.add('%s\t%s\t%s\n' % (gene_name, db_name, DB_namekegg[db_name]))
    f1 = open(kegg_result + '.kegg.sum','w')
    f1.write('gene_name\tKO\tBRITE_KO1\tBRITE_KO2\tBRITE_KO3\n' + ''.join(list(Output)))
    f1.close()

for kegg_result in allkegg_results:
    kegg_result_folder, kegg_result_file = os.path.split(kegg_result)
    # set up gene annotation and clustering
    Blast_output1, DB_name1, DB_name1_2 = annotate_kegg(kegg_result)
    # sum up
    sum_kegg(Blast_output1,DB_name1_2)
################################################### END ########################################################
