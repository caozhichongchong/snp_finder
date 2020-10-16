import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
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

mutation_dir = args.o + '/summary'
mutation_select_fasta = os.path.join(mutation_dir,'all.selected.gene.faa.High_select2.faa')
mutation_fasta = os.path.join(mutation_dir,'all.denovo.gene.faa')
output_dir = args.o + '/pangenome'
allresult = glob.glob(os.path.join(output_dir,'*.denovo.txt'))
core_cutoff = 0.8
core_cutoff_species_num = 15
# load diamond results
def load_result(file,speciesname,genomename,Gene,Gene_copy,alloutput):
    resultset = []
    tempspeciesset = dict()
    for lines in open(file,'r'):
        lines_set = lines.split('\t')
        geneID = lines_set[1]
        resultset.append(geneID)
        tempspeciesset.setdefault(geneID,set())
        tempspeciesset[geneID].add(speciesname)
        alloutput.append('%s\t%s\t%s\t\n'%(geneID,genomename,speciesname))
        Gene_copy[geneID].append(speciesname)
    for geneID in tempspeciesset:
        Gene[geneID]+=list(tempspeciesset[geneID])
    return [Gene,Gene_copy,alloutput]

# selected genes
# load selected genes
Gene_select = []
for record in SeqIO.parse(mutation_select_fasta, 'fasta'):
    geneID = str(record.id)
    Gene_select.append(geneID)

# load all genes
Gene = dict()
Gene_copy = dict()
Gene_summary = dict()
Gene_length = dict()
for record in SeqIO.parse(mutation_fasta, 'fasta'):
    geneID = str(record.id)
    gene_length = len(str(record.seq))
    Gene.setdefault(geneID,[])
    Gene_copy.setdefault(geneID, [])
    Gene_summary.setdefault(geneID, set())
    Gene_length.setdefault(geneID,gene_length)

# check all genes
for geneID in Gene_select:
    if geneID not in Gene:
        print('missing geneID %s in %s'%(geneID,mutation_fasta))

# summarize all genomes for a species
Species = dict()
Species_set = set()
try:
    f1 = open(mutation_fasta + '.allpangenome.txt','r')
    Gene_temp = dict()
    for lines in f1:
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            geneID,genomename,speciesname = lines_set[0:3]
            Species.setdefault(speciesname, set())
            Species[speciesname].add(genomename)
            Gene_temp.setdefault(geneID, set())
            Gene_temp[geneID].add(genomename)
            Gene_copy[geneID].append(speciesname)
            Species_set.add(speciesname)
    for geneID in Gene_temp:
        for genomename in Gene_temp[geneID]:
            genomename_set = genomename.split('_')
            speciesname = genomename_set[1].replace('.', '').replace('Parabacteroides', 'PB')
            if genomename.startswith('Bfrag') or speciesname in ['Bact', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7',
                                                                 'T8']:
                speciesname = 'BaFr'
            Gene[geneID].append(speciesname)
except IOError:
    alloutput = []
    alloutput.append('#geneID\tgenome\tspecies\t\n')
    for resultfile in allresult:
        genomename = os.path.split(resultfile)[-1].split('.faa.selected.txt')[0]
        genomename_set = genomename.split('_')
        speciesname = genomename_set[1].replace('.','').replace('Parabacteroides','PB')
        if genomename.startswith('Bfrag') or speciesname in ['Bact','T1','T2','T3','T4','T5','T6','T7','T8']:
            speciesname = 'BaFr'
        Species.setdefault(speciesname,set())
        Species[speciesname].add(genomename)
        Gene,Gene_copy,alloutput = load_result(resultfile,speciesname,genomename,Gene,Gene_copy,alloutput)
        Species_set.add(speciesname)
    f1 = open(mutation_fasta + '.allpangenome.txt','w')
    f1.write(''.join(alloutput))
    f1.close()
print(Species_set)

# summarize gene distribution in a species
#allsummary = []
allsummary2 = set()
geneID_list = 'Species\tTotalgenome'
for geneID in Gene:
    geneID_list += '\t%s'%(geneID)

geneID_list += '\n'
#allsummary.append(geneID_list)
for speciesname in Species:
    Totalgenome = len(Species[speciesname])
    geneID_list = '%s\t%s'%(speciesname,Totalgenome)
    for geneID in Gene:
        Totalgenome_geneID = Gene[geneID].count(speciesname)
        Totalcopy_geneID = Gene_copy[geneID].count(speciesname)
        Percentage = Totalgenome_geneID/Totalgenome
        Percentage_copy = Totalcopy_geneID/Totalgenome
        geneID_list += '\t%.3f' % (Percentage)
        if speciesname in geneID:
            # same species
            if Percentage >= core_cutoff:
                Gene_summary[geneID].add('origin_core')
            elif Percentage > 0:
                    Gene_summary[geneID].add('origin_flexible')
        else:
            if Percentage >= core_cutoff:
                Gene_summary[geneID].add('other_core')
            elif Percentage > 0:
                    Gene_summary[geneID].add('other_flexible')
        if geneID in Gene_select:
            tag = 'selected'
        else:
            tag = 'nonselected'
        if Percentage > 0:
            allsummary2.add('%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n'%(geneID,speciesname,Gene_length[geneID],Totalgenome_geneID,
                                                          Percentage,Percentage_copy,tag))
    geneID_list += '\n'
    #allsummary.append(geneID_list)

Gene_summary_output=[]
Gene_summary_output.append('#geneID\tspecies_num\tselected\tmultispecies\tspecies_list\tcore_or_flexible\t')
for geneID in Gene:
    allspecies = Gene[geneID]
    allspecies_unique = list(set(allspecies))
    Totalspecies = len(allspecies_unique)
    tempoutput = ('%s\t%s\t'%(geneID,Totalspecies))
    if geneID in Gene_select:
        tempoutput += 'selected\t'
    else:
        tempoutput += 'nonselected\t'
    if Totalspecies > 1:
        tempoutput += ('multispecies\t')
    else:
        tempoutput += ('singlespecies\t')
    tempoutput += ';'.join(allspecies_unique) + '\t'
    if 'origin_flexible' in Gene_summary[geneID]:
        tag = 'species_flexible'
    elif 'origin_core' in Gene_summary[geneID]:
        if Totalspecies >= core_cutoff_species_num and 'other_core' in Gene_summary[geneID]:
            # potential core gene
            tag = 'all_core'
        else:
            tag = 'origin_core'
    else:
        tag = 'other'
    tempoutput += tag + '\t'
    for tag in Gene_summary[geneID]:
        tempoutput += '%s\t'%(tag)
    Gene_summary_output.append(tempoutput)

# output results
f1 = open(mutation_fasta + '.allpangenome.sum.txt','w')
f1.write('\n'.join(Gene_summary_output))
f1.close()

f1 = open(mutation_fasta + '.allpangenome.sum.species.txt','w')
f1.write('#geneID\tspecies\tgene_length\tgenome_percentage\tavg_copy_num\tselected\t\n' + ''.join(list(allsummary2)))
f1.close()
