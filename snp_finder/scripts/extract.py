# start
# after calculate NS ratio, extract genes and run annotation
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
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
ref_filename = '.all.spades*.fasta'
input_summary = output_dir_merge + '/summary/all.donor.species.dnds.txt'
output_gene = output_dir_merge + '/summary/all.selected.gene.faa'
output_gene_dna = output_dir_merge + '/summary/all.selected.gene.fna'

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# function
def sum_donor_species(input_summary,selected = 1):
    Donor_species = dict()
    for lines in open(input_summary,'r'):
        lines_set = lines.replace('\n','').split('\t')
        if not lines.startswith('#'):
            if selected == 0 or lines_set[-1] == 'True':
                # highly selected genes or all genes
                gene_name = lines_set[1]
                donor_species = lines_set[0]
                Donor_species.setdefault(donor_species,[])
                Donor_species[donor_species].append(gene_name)
    return Donor_species

def sum_vcf(vcf,Donor_species):
    donor_species = os.path.split(vcf)[-1].split('.all.flt.snp.vcf')[0]
    Donor_species.setdefault(donor_species, [])
    for lines in open(vcf,'r'):
        lines_set = lines.replace('\n','').split('\t')
        if not lines.startswith('#'):
            if '*' in lines_set[-1]:
                # highly selected genes or all genes
                gene_name = lines_set[-4]
                Donor_species[donor_species].append(gene_name)
    return Donor_species

def extract_donor_species(donor_species,input_fasta,gene_name_list,output_fasta):
    gene_name_extract = []
    donor_species_set = donor_species.split('_')
    record_id_set = set()
    try:
        donor_species = '%s_%s_%s' % (donor_species_set[0],
                                      donor_species_set[1][0:min(6, len(donor_species_set[1]))],
                                      donor_species_set[2][0:min(6, len(donor_species_set[2]))])
    except IndexError:
        donor_species = '%s_%s_%s' % (donor_species_set[0],
                                      donor_species_set[1],
                                      donor_species_set[1])
    if 'cluster' in donor_species_set[-1]:
        try:
            donor_species += '_CL' + donor_species_set[-2].split('cluster')[1]
        except IndexError:
            donor_species += '_CL' + donor_species_set[-1].split('cluster')[1]
    gene_name_list = set(gene_name_list)
    for record in SeqIO.parse(input_fasta, 'fasta'):
        record_id = str(record.id)
        temp_line = '%s\t%s\t%s\t'%('_'.join(donor_species_set),donor_species,record_id)
        if record_id in gene_name_list and record_id not in record_id_set:
            gene_name_extract.append(record_id)
            record_id = 'C_%s_G_%s'%(record_id.split('_')[1], record_id.split('_')[-1])
            output_fasta.append('>%s__%s\n%s\n'%(donor_species, record_id, str(record.seq)))
            record_id_set.add(str(record.id))
            temp_line += '%s__%s\t\n'%(donor_species,record_id)
            change_name.add(temp_line)
    if len(gene_name_extract) < len(gene_name_list):
        print('missing genes in files %s extracted %s genes of %s genes'%(input_fasta,len(gene_name_extract),len(gene_name_list)))
        print('missing genes ' + ' '.join([gene_name for gene_name in gene_name_list if gene_name not in gene_name_extract]))
    return output_fasta

# extract highly selected sequences
output_fasta = []
output_fasta_dna = []
change_name = set()
Donor_species = sum_donor_species(input_summary)
for donor_species in Donor_species:
    print(donor_species)
    database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
    database = database[0]
    ref_dir, ref_name = os.path.split(database)
    input_fasta = database.replace('.fasta', '.faa')
    input_fasta_dna = database.replace('.fasta', '.faa')
    gene_name_list = Donor_species[donor_species]
    try:
        f1 = open(input_fasta, 'r')
    except FileNotFoundError:
        os.system(args.pro + ' -q -i %s -a %s' % (database, input_fasta))
    try:
        f1 = open(input_fasta_dna, 'r')
    except FileNotFoundError:
        os.system(args.pro + ' -q -i %s -d %s' % (database, input_fasta_dna))
    if gene_name_list != []:
        if input_fasta != [] and input_fasta_dna != []:
            output_fasta = extract_donor_species(donor_species,input_fasta, gene_name_list, output_fasta)
            output_fasta_dna = extract_donor_species(donor_species,input_fasta_dna, gene_name_list, output_fasta_dna)
        else:
            print('missing files for %s' % (donor_species))
    else:
        print('missing genes for %s' % (donor_species))

f1 = open(output_gene, 'w')
f1.write(''.join(output_fasta))
f1.close()
f1 = open(output_gene + '.changename.txt', 'w')
f1.write(''.join(list(change_name)))
f1.close()
f1 = open(output_gene_dna, 'w')
f1.write(''.join(output_fasta_dna))
f1.close()

# extract all sequences
output_gene = output_dir_merge + '/summary_jay/all.denovo.gene.faa'
output_gene_dna = output_dir_merge + '/summary_jay/all.denovo.gene.fna'
output_fasta = []
output_fasta_dna = []
Donor_species = sum_donor_species(input_summary,0)
change_name = set()
for donor_species in Donor_species:
    print(donor_species)
    database = glob.glob('%s/%s/%s%s' % (genome_root, donor_species, donor_species, ref_filename))
    database = database[0]
    ref_dir, ref_name = os.path.split(database)
    input_fasta = database.replace('.fasta', '.faa')
    input_fasta_dna = database.replace('.fasta', '.faa')
    gene_name_list = Donor_species[donor_species]
    try:
        f1 = open(input_fasta, 'r')
    except FileNotFoundError:
        os.system('prodigal -q -i %s -a %s' % (database, input_fasta))
    try:
        f1 = open(input_fasta_dna, 'r')
    except FileNotFoundError:
        os.system('prodigal -q -i %s -d %s' % (database, input_fasta_dna))
    if gene_name_list != []:
        if input_fasta != [] and input_fasta_dna != []:
            output_fasta = extract_donor_species(donor_species,input_fasta, gene_name_list, output_fasta)
            output_fasta_dna = extract_donor_species(donor_species,input_fasta_dna, gene_name_list, output_fasta_dna)
        else:
            print('missing files for %s' % (donor_species))
    else:
        print('missing genes for %s' % (donor_species))

f1 = open(output_gene, 'w')
f1.write(''.join(output_fasta))
f1.close()
f1 = open(output_gene + '.changename.txt', 'w')
f1.write(''.join(list(change_name)))
f1.close()
f1 = open(output_gene_dna, 'w')
f1.write(''.join(output_fasta_dna))
f1.close()

# run cluster
cutoff = 0.7
cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                   % (args.u, output_gene, cutoff, output_gene,
                      output_gene, args.t))
os.system(cmd_cluster)

################################################### END ########################################################
