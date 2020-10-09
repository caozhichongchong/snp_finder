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
    High_select2_3SNP = dict()
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
                if species not in SNP_cutoff:
                    High_select2.add(record_name2)
                    High_select2.add(record_name)
                else:
                    High_select2_3SNP.setdefault(record_name,set())
                    High_select2_3SNP[record_name].add(record_name2)
                    High_select2_3SNP.setdefault(record_name2, set())
                    High_select2_3SNP[record_name2].add(record_name)
    if High_select2_3SNP!= dict():
        for record in High_select2_3SNP:
            if len(High_select2_3SNP[record]) >= 3:
                High_select2.append(record)
    return [Clusters,High_select2_output,High_select2]

# correct highly selected genes dNdS
def format_freq(freq):
    freq = freq.split(':')
    return [int(freq[0]),int(freq[1])]

def init_highselect_notHS(line_set):
    allspecies_highselect = line_set
    for i in [3,5,6,7,8,-13,-11,-9,-7,-5,-3]:
        allspecies_highselect[i] = int(allspecies_highselect[i])
    for i in [-12,-10,-8,-6,-4,-2]:
        allspecies_highselect[i] = format_freq(line_set[i])
    return allspecies_highselect

def init_highselect_empty(line_set):
    allspecies_highselect = line_set
    for i in [3]:
        allspecies_highselect[i] = int(allspecies_highselect[i])
    for i in [5,6,7,8,9,10,11,12,-13,-11,-9,-7,-5,-3]:
        allspecies_highselect[i] = 0
    for i in [-12,-10,-8,-6,-4,-2]:
        allspecies_highselect[i] = format_freq('0:0')
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
    # expect NS
    # expected NS ratio
    tempNS = [0, 0]
    i = -12
    for freq in allspecies:
        tempNS[0] += freq * allspecies_highselect[i][0]
        tempNS[1] += freq * allspecies_highselect[i][1]
        i += 2
    if tempNS[1] > 0:
        allspecies_highselect[11] = tempNS[0] / tempNS[1] # expect
    else:
        allspecies_highselect[11] = 'expect_N_only'
    # NS ratio
    if allspecies_highselect[8] > 0:
        allspecies_highselect[10] = allspecies_highselect[7]/allspecies_highselect[8] # observe
        if allspecies_highselect[11] != 'expect_N_only':
            # dNdS
            allspecies_highselect[12] = allspecies_highselect[10] / allspecies_highselect[11]# dnds
    elif allspecies_highselect[7] > 0:
        allspecies_highselect[10] = 'observe_N_only'
        allspecies_highselect[12] = 'observe_N_only'
    return allspecies_highselect

def output_highlight(allspecies_highselect,allspecies,newoutput):
    allspecies_highselect = calculate_NS(allspecies, allspecies_highselect)
    for i in range(0, len(allspecies_highselect)):
        if i in [14, 16, 18, 20, 22, 24]:
            allspecies_highselect[i] = '%s:%s' % (allspecies_highselect[i][0],
                                                  allspecies_highselect[i][1])
        else:
            allspecies_highselect[i] = str(allspecies_highselect[i])
    newoutput.append('\t'.join(allspecies_highselect) + '\n')
    return newoutput

def genus_species(donor_species):
    genus = donor_species.split('_')[1].replace('BA', 'Bifidobacterium').replace('BL',
                                                                                 'Bifidobacterium').replace(
        'PB', 'Parabacteroides')
    if any(item in donor_species for item in ['BA','BL','PB']):
        species = donor_species.split('_')[1].replace('BA', 'Bifidobacterium_adolescentis').replace('BL',
                                                                                                      'Bifidobacterium_longum').replace(
            'PB', 'Parabacteroides_butyrate')
    else:
        species = (donor_species.split('_')[1] + '_' +
                   donor_species.split('_')[2])
    return [genus,species]

def add_new_selection(input_summary,High_select2):
    Donor_species = dict()
    Donor_species_notHS = dict()
    Genus = dict()
    Genus_notHS = dict()
    Species = dict()
    Species_notHS = dict()
    newoutput = []
    # use selected codon NS ratio (SNP pair) * all genes freq of no single event
    for lines in open(input_summary_nosingleevent, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        Selected = line_set[-1]
        if Selected == 'False':
            donor_species = line_set[0]
            record_id = line_set[1]
            # all species
            if record_id == 'allspecies':
                # all genes freq
                allspecies = [int(line_set[-13]), int(line_set[-11]), int(line_set[-9]),
                              int(line_set[-7]), int(line_set[-5]), int(line_set[-3])]
                allspecies_highselect = init_highselect_empty(line_set)
                allspecies_highselect[1] = 'allspecies_highselect'
                newoutput.append(lines)
            elif record_id in ['allspecies_flexible', 'allspecies_core']:
                newoutput.append(lines)
            elif donor_species == record_id:
                # donor_species not high select set up
                line_set = lines.split('\n')[0].split('\t')
                donor_species_highselect = init_highselect_notHS(line_set)
                Donor_species_notHS.setdefault(donor_species,
                                               donor_species_highselect)
                line_set = lines.split('\n')[0].split('\t')
                donor_species_highselect = init_highselect_empty(line_set)
                donor_species_highselect[1] = donor_species + '_highselect'
                Donor_species.setdefault(donor_species,
                                         donor_species_highselect)
                genus, species = genus_species(donor_species)
                line_set = lines.split('\n')[0].split('\t')
                if genus not in Genus_notHS:
                    genus_highselect = init_highselect_notHS(line_set)
                    genus_highselect[0] = genus
                    genus_highselect[1] = genus
                    Genus_notHS.setdefault(genus, genus_highselect)
                    line_set = lines.split('\n')[0].split('\t')
                    genus_highselect = init_highselect_empty(line_set)
                    genus_highselect[0] = genus
                    genus_highselect[1] = genus + '_highselect'
                    Genus.setdefault(genus, genus_highselect)
                else:
                    genus_highselect = Genus_notHS[genus]
                    Genus_notHS[genus] = add_highselect(line_set, genus_highselect)
                line_set = lines.split('\n')[0].split('\t')
                if species not in Species_notHS:
                    species_highselect = init_highselect_notHS(line_set)
                    species_highselect[0] = species
                    species_highselect[1] = species
                    Species_notHS.setdefault(species, species_highselect)
                    line_set = lines.split('\n')[0].split('\t')
                    species_highselect = init_highselect_empty(line_set)
                    species_highselect[0] = species
                    species_highselect[1] = species + '_highselect'
                    Species.setdefault(species, species_highselect)
                else:
                    species_highselect = Species_notHS[species]
                    Species_notHS[species] = add_highselect(line_set, species_highselect)
    # calculate NS for new HS genes
    for lines in open(input_summary, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        Selected = line_set[-1]
        donor_species = line_set[0]
        record_id = line_set[1]
        # species or genes
        if Selected == 'False':
            if record_id not in ['0', 'gene', 'allspecies', 'allspecies_highselect',
                                 'allspecies_core', 'allspecies_flexible',
                                 'allspecies_highselect_noNS']:
                genus, species = genus_species(donor_species)
                donor_species_set = donor_species.split('_')
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
                record_id = '%s__C_%s_G_%s' % (donor_species, record_id.split('_')[1], record_id.split('_')[-1])
                # highly selected genes
                if record_id in High_select2:
                    if (line_set[10] == 'observe_N_only' or float(line_set[10]) > NSratioobserve_cutoff):
                        newoutput.append('\t'.join(line_set[:-1]) + '\tTrue\n')
                        allspecies_highselect = add_highselect(line_set, allspecies_highselect)
                        Donor_species[line_set[0]] = add_highselect(line_set, Donor_species[line_set[0]])
                        Genus[genus] = add_highselect(line_set, Genus[genus])
                        Species[species] = add_highselect(line_set, Species[species])
                    else:
                        newoutput.append('\t'.join(line_set[:-1]) + '\tFalse\n')
                        High_select2.remove(record_id)
                # other genes and species
                elif line_set[0] not in line_set[1]:
                    newoutput.append(lines)
                else:
                    print('not output %s' % (line_set[1]))
        # original highly selected
        elif Selected == 'True':
            newoutput.append(lines)
            allspecies_highselect = add_highselect(line_set, allspecies_highselect)
            Donor_species[line_set[0]] = add_highselect(line_set, Donor_species[line_set[0]])
            Genus[genus] = add_highselect(line_set, Genus[genus])
            Species[species] = add_highselect(line_set, Species[species])
        else:
            newoutput.append(lines)
    # output new highselect
    newoutput = output_highlight(allspecies_highselect, allspecies, newoutput)
    for donor_species in Donor_species:
        line_set = Donor_species_notHS[donor_species]
        line_set_sub = [int(line_set[-13]), int(line_set[-11]), int(line_set[-9]),
                        int(line_set[-7]), int(line_set[-5]), int(line_set[-3])]
        newoutput = output_highlight(Donor_species[donor_species],
                                     line_set_sub, newoutput)
        newoutput = output_highlight(Donor_species_notHS[donor_species],
                                     line_set_sub, newoutput)
    for genus in Genus:
        line_set = Genus_notHS[genus]
        line_set_sub = [int(line_set[-13]), int(line_set[-11]), int(line_set[-9]),
                        int(line_set[-7]), int(line_set[-5]), int(line_set[-3])]
        newoutput = output_highlight(Genus[genus],
                                     line_set_sub, newoutput)
        newoutput = output_highlight(Genus_notHS[genus],
                                     line_set_sub, newoutput)
    for species in Species:
        line_set = Species_notHS[species]
        line_set_sub = [int(line_set[-13]), int(line_set[-11]), int(line_set[-9]),
                        int(line_set[-7]), int(line_set[-5]), int(line_set[-3])]
        newoutput = output_highlight(Species[species],
                                     line_set_sub, newoutput)
        newoutput = output_highlight(Species_notHS[species],
                                     line_set_sub, newoutput)
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

def annotation(all_filter_gene_fasta_file,pre_cluster = ''):
    # run cluster
    cutoff = 0.7
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                   % ('usearch', all_filter_gene_fasta_file, cutoff, all_filter_gene_fasta_file,
                      all_filter_gene_fasta_file, 40))
    os.system(cmd_cluster)
    all_filter_gene_fasta_file = all_filter_gene_fasta_file + '.cluster.aa'
    if pre_cluster!= '':
        os.system('#%s -makeudb_usearch %s -output %s.udb' %
                  ('usearch', pre_cluster, pre_cluster))
        os.system('%s -ublast %s -db %s.udb  -evalue 1e-2 -accel 0.5 -blast6out %s -threads 2'%
                  ('usearch', all_filter_gene_fasta_file,pre_cluster, all_filter_gene_fasta_file + '.ref.out.txt'))
    # run prokka
    cmdsprokka = 'py37\nprokka --kingdom Bacteria --outdir %s/prokka_%s  --protein %s --locustag Bacter %s/%s\n' % \
                 (output_dir_merge + '/summary', os.path.split(all_filter_gene_fasta_file)[-1],
                  all_filter_gene_fasta_file,
                  output_dir_merge + '/summary',
                  os.path.split(all_filter_gene_fasta_file)[-1].replace('.faa', '.fna'))
    f1 = open(os.path.join(input_script_sub, 'prokka.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmdsprokka))
    f1.close()
    # run metacyc
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/mit_alm/database/metacyc/protseq.fsa'
    cmds = ("diamond blastp --query %s --db %s.dmnd --out %s.metacyc.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    f1 = open(os.path.join(input_script_sub, 'metacyc.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    # run eggnog
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/eggnog/xaa.hmm'
    cmds = ('hmmsearch --tblout %s.eggnog.1.txt --cpu 40 -E %s %s %s\n') %(all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.1.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xab.hmm'
    cmds = ('hmmsearch --tblout %s.eggnog.2.txt --cpu 40 -E %s %s %s\n') % (
    all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.2.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xac.hmm'
    cmds = ('hmmsearch --tblout %s.eggnog.3.txt --cpu 40 -E %s %s %s\n') % (
    all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.3.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmds))
    f1.close()
    # run kegg
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/kegg/kofam/profiles/prokaryote/prokaryote.hmm'
    cmds = ('hmmsearch --tblout %s.kegg.txt --cpu 40 -E %s %s %s\n') %(all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'kegg.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    # run customed database
    cutoff = 80
    cutoff2 = 80
    cmds = ''
    database = '/scratch/users/anniz44/scripts/database/SARG.db.fasta'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.SARG.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 50
    cutoff2 = 50
    database = '/scratch/users/anniz44/scripts/database/AHR.aa.db'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.AHR.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 60
    cutoff2 = 80
    database = '/scratch/users/anniz44/scripts/database/Butyrate.pro.aa'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.buty.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/anniz44/scripts/database/IntI1_database.fasta'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.int.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/anniz44/scripts/database/SRB.AA'
    cmds += ("diamond blastp --query %s --db %s.dmnd --out %s.SRB.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    cutoff = 0.01
    database = '/scratch/users/anniz44/scripts/database/NR.hmm'
    cmds += ('hmmsearch --tblout %s.NR.txt --cpu 40 -E %s %s %s\n') %(all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'customed.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
    f1.close()
    # all scripts
    f1 = open(os.path.join(input_script, 'allannotate.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.sh')):
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
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
annotation(all_fasta_HS + '.High_select2.faa',pre_cluster)

################################################### END ########################################################
