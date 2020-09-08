# start
# find highly selected genes in one species across populations and donors
import os
import glob
import copy
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of folders of WGS of each species",
                      type=str, default='.',
                      metavar='input/')
required.add_argument("-fq",
                      help="file extension of WGS fastq #1 files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')
# optional output setup

optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='snp_output/',
                      metavar='snp_output/')
# optional search parameters
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 1)",
                      metavar="1 or more", action='store', default=1, type=int)
optional.add_argument('-job',
                      help="Optional: command to submit jobs",
                      metavar="nohup or customized",
                      action='store', default='jobmit', type=str)
# requirement for software calling
optional.add_argument('-bw', '--bowtie',
                          help="Optional: complete path to bowtie if not in PATH",
                          metavar="/usr/local/bin/bowtie",
                          action='store', default='bowtie', type=str)
optional.add_argument('-sp', '--spades',
                          help="Optional: complete path to spades if not in PATH",
                          metavar="/usr/local/bin/spades",
                          action='store', default='spades', type=str)
optional.add_argument('-pro', '--prodigal',
                      help="Optional: complete path to prodigal if not in PATH, None for no prodigal (default)",
                      metavar="/usr/local/bin/prodigal",
                      action='store', default='None', type=str)
optional.add_argument('-bcf', '--bcftools',
                      help="Optional: complete path to bcftools if not in PATH",
                      metavar="/usr/local/bin/bcftools",
                      action='store', default='bcftools', type=str)
optional.add_argument('-sam', '--samtools',
                      help="Optional: complete path to bwa if not in PATH",
                      metavar="/usr/local/bin/samtools",
                      action='store', default='samtools', type=str)
optional.add_argument('-mini', '--minimap2',
                      help="Optional: complete path to minimap2 if not in PATH",
                      metavar="/usr/local/bin/minimap2",
                      action='store', default='minimap2', type=str)
optional.add_argument('--u','--usearch',
                        help="Optional: cluster genes with SNPs",
                        metavar="usearch",
                        action='store', default='usearch', type=str)
################################################## Definition ########################################################
args = parser.parse_args()

input_script_sub = args.s +'/annotate'
input_script = args.s
genome_root = args.i + '/round*'
output_dir_merge = args.o +'/vcf_round%s/merge_genome/'%(Round)
all_fasta = os.path.join(output_dir_merge + '/summary', 'all.denovo.gene.faa')
all_fasta_HS = os.path.join(output_dir_merge + '/summary', 'all.selected.gene.faa')
input_summary = output_dir_merge + '/summary/all.donor.species.dnds.txt'

# functions
# clustering
def cluster_uc(cluster_input):
    Clusters = dict()
    High_select2 = set()
    High_select2_output = []
    for lines in open(cluster_input, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        cluster = line_set[1]
        record_name = line_set[8]
        record_name2 = line_set[9]
        Clusters.setdefault(record_name, cluster)
        if record_name2!= '*':
            donor = record_name.split('_')[0]
            donor_species = record_name.split('__')[0]
            species = donor_species.replace(donor + '_', '')
            donor2 = record_name2.split('_')[0]
            donor_species2 = record_name2.split('__')[0]
            species2 = donor_species2.replace(donor2 + '_', '')
            if species == species2:
                High_select2.add(record_name2)
                High_select2.add(record_name)
    if High_select2 != set():
        for record_name in High_select2:
            High_select2_output.append('%s\tTrue\t\n'%(record_name))
    return [Clusters,High_select2_output,High_select2]

# correct highly selected genes dNdS
def format_freq(freq):
    freq = freq.split(':')
    return [int(freq[0]),int(freq[1])]

def init_highselect(line_set):
    allspecies_highselect = line_set
    for i in [3,5,6,7,8,-13,-11,-9,-7,-5,-3]:
        allspecies_highselect[i] = 0
    for i in [-12,-10,-8,-6,-4,-2]:
        allspecies_highselect[i] = format_freq(line_set[i])
    return allspecies_highselect

def add_freq(freq_new,freq_old):
    temp_result = format_freq(freq_new)
    freq_old[0] += temp_result[0]
    freq_old[1] += temp_result[1]

def add_highselect(line_set,allspecies_highselect):
    for i in [3,5,6,7,8,-13,-11,-9,-7,-5,-3]:
        allspecies_highselect[i] += int(line_set[i])
    for i in [-12,-10,-8,-6,-4,-2]:
        add_freq(line_set[i], allspecies_highselect[i])
    return allspecies_highselect

def calculate_NS(allspecies,allspecies_highselect):
    # NS ratio
    allspecies_highselect[10] = allspecies_highselect[7]/allspecies_highselect[8]
    # expected NS ratio
    tempNS = [0,0]
    i = -12
    for freq in allspecies:
        tempNS[0] += freq * allspecies_highselect[i][0]
        tempNS[1] += freq * allspecies_highselect[i][1]
        i += 2
    allspecies_highselect[11] = tempNS[0] / tempNS[1]
    # dNdS
    allspecies_highselect[12] = allspecies_highselect[10] / allspecies_highselect[11]
    return allspecies_highselect

def add_new_selection(input_summary,High_select2):
    newoutput = []
    # use selected codon NS ratio (SNP pair) * all genes freq
    for lines in open(input_summary, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        Selected = line_set[-1]
        if Selected == 'False':
            donor_species = line_set[0]
            record_id = line_set[1]
            # all species
            if record_id == 'allspecies':
                # all genes freq
                newoutput.append(lines)
                allspecies = [int(line_set[-13]), int(line_set[-11]), int(line_set[-9]),
                              int(line_set[-7]), int(line_set[-5]), int(line_set[-3])]
            # all highly selected genes
            elif record_id == 'allspecies_highselect':
                # selected codon NS ratio (SNP pair)
                allspecies_highselect = init_highselect(line_set)
    for lines in open(input_summary, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        Selected = line_set[-1]
        if Selected == 'False':
            donor_species = line_set[0]
            record_id = line_set[1]
            # species or genes
            if record_id not in ['allspecies', 'allspecies_highselect']:
                donor_species_set = donor_species.split('_')
                try:
                    donor_species = '%s_%s_%s' % (donor_species_set[0],
                                                  donor_species_set[1][0:min(6, len(donor_species_set[1]))],
                                                  donor_species_set[2][0:min(6, len(donor_species_set[2]))])
                except IndexError:
                    donor_species = '%s_%s_%s' % (donor_species_set[0],
                                                  donor_species_set[1],
                                                  donor_species_set[1])
                record_id = '%s__C_%s_G_%s' % (donor_species, record_id.split('_')[1], record_id.split('_')[-1])
                # highly selected genes
                if record_id in High_select2:
                    if (line_set[10] == 'observe_N_only' or float(line_set[10]) >= 1):
                        newoutput.append('\t'.join(line_set[:-1]) + '\tTrue\n')
                        allspecies_highselect = add_highselect(line_set, allspecies_highselect)
                    else:
                        High_select2.remove(record_id)
                # other genes and species
                else:
                    newoutput.append(lines)
        # original highly selected
        elif Selected == 'True':
            newoutput.append(lines)
            allspecies_highselect = add_highselect(line_set, allspecies_highselect)
    # output new_allspecies_highselect
    allspecies_highselect = calculate_NS(allspecies, allspecies_highselect)
    print(allspecies_highselect, allspecies)
    for i in range(0, len(allspecies_highselect)):
        if i in [14, 16, 18, 20, 22, 24]:
            allspecies_highselect[i] = '%s:%s' % (allspecies_highselect[i][0],
                                                  allspecies_highselect[i][1])
        else:
            allspecies_highselect[i] = str(allspecies_highselect[i])
    print(allspecies_highselect, allspecies)
    newoutput.append('\t'.join(allspecies_highselect) + '\n')
    # output new summary
    f1 = open(input_summary + '.High_select2.txt', 'w')
    f1.write(''.join(newoutput))
    f1.close()
    return High_select2

def sum_gene(input_summary,High_select2, selected = 1):
    genelist = []
    genelist += High_select2
    for lines in open(input_summary,'r'):
        lines_set = lines.replace('\n','').split('\t')
        if not lines.startswith('#'):
            if selected == 0 or lines_set[-1] == 'True':
                # highly selected genes or all genes
                gene_name = lines_set[1]
                genelist.append(gene_name)
    print(len(genelist))
    genelist = list(set(genelist))
    print(len(genelist))
    Output = []
    output_set = []
    for record in SeqIO.parse(all_fasta_HS, 'fasta'):
        Output.append('>%s\n%s\n' % (str(record.id), str(record.seq)))
        output_set.append(str(record.id))
    for record in SeqIO.parse(all_fasta, 'fasta'):
        if str(record.id) in genelist and str(record.id) not in output_set:
            Output.append('>%s\n%s\n' % (str(record.id), str(record.seq)))
            output_set.append(str(record.id))
    f1 = open(all_fasta_HS + '.High_select2.faa' , 'w')
    f1.write(''.join(Output))
    f1.close()

# clustering
Clusters_gene, High_select2_output,High_select2 = cluster_uc(all_fasta + '.uc')
f1 = open(all_fasta + '.High_select2.txt', 'w')
f1.write(''.join(High_select2_output))
f1.close()
# correcting
High_select2 = add_new_selection(input_summary,High_select2)

# run clustering
sum_gene(input_summary,High_select2,1)

################################################### END ########################################################
