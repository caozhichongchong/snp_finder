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
                      help="path of all co-assemblies",
                      type=str, default='.',
                      metavar='input/')
optional.add_argument("-meta",
                      help="a metadata of taxonomy of each lineage",
                      type=str, default='None',
                      metavar='metadata.txt')
optional.add_argument("-cutoff",
                      help="a file of cutoff of how many SNPs on a gene for each clonal population to call parallel evolution",
                      type=str, default='None',
                      metavar='total_SNP_cutoff.txt')
optional.add_argument("-sig",
                      help="a file of lineages with PE significance (PEsig2)",
                      type=str, default='None',
                      metavar='total_SNP_cutoff.txt')
# optional output setup
optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='snp_output/',
                      metavar='snp_output/')
optional.add_argument("-trunc",
                      help="extract truncated genes",
                      type=str, default='True',
                      metavar='False or True')
# optional search parameters
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 1)",
                      metavar="1 or more", action='store', default=1, type=int)
optional.add_argument('-job',
                      help="Optional: command to submit jobs",
                      metavar="nohup or customized",
                      action='store', default='jobmit', type=str)
# requirement for software calling
optional.add_argument('-pro',
                      help="Optional: complete path to prodigal if not in PATH",
                      metavar="/usr/local/bin/prodigal",
                      action='store', default='prodigal', type=str)
optional.add_argument('--u',
                        help="Optional: complete path to usearch if not in PATH",
                        metavar="/usr/local/bin/usearch",
                        action='store', default='usearch', type=str)
optional.add_argument('--dm',
                        help="Optional: complete path to diamond if not in PATH",
                        metavar="/usr/local/bin/diamond",
                        action='store', default='diamond', type=str)
optional.add_argument('--hmm',
                        help="Optional: complete path to hmmsearch if not in PATH",
                        metavar="/usr/local/bin/hmmsearch",
                        action='store', default='hmmsearch', type=str)
optional.add_argument('--prokka',
                        help="Optional: complete path to prokka if not in PATH",
                        metavar="/usr/local/bin/prokka",
                        action='store', default='prokka', type=str)


################################################## Definition ########################################################
args = parser.parse_args()

input_script_sub = args.s +'/annotate'
input_script_sub_all = args.s +'/annotate_all'
input_script_sub_trunc = args.s +'/annotate_trunc'
input_script = args.s
genome_root = args.i
output_dir_merge = args.o
ref_filename = '.noHM.fasta'
input_summary = output_dir_merge + '/summary/all.species.txt'
output_gene = output_dir_merge + '/summary/all.selected.gene.faa'
output_gene_dna = output_dir_merge + '/summary/all.selected.gene.fna'
compute_parallel = False # whether to find PE across lineages
annotate_customized = False
cluster_cutoff = 0.9
parallel_across_donor = True # across donor or across lineage
include_across_donor_PE = False # including PE across donors when computing dnds
try:
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(input_script_sub_all)
except IOError:
    pass

if args.trunc!= 'False':
    try:
        os.mkdir(input_script_sub_trunc)
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
                if '_other' not in gene_name and 'NODE' in gene_name:
                    donor_species = lines_set[0]
                    Donor_species.setdefault(donor_species,[])
                    Donor_species[donor_species].append(gene_name)
    return Donor_species

def sum_vcf(allvcf):
    Donor_species = dict()
    for vcf in allvcf:
        donor_species = os.path.split(vcf)[-1].split('.raw.vcf.filtered.vcf.final.snp.txt')[0].replace('.all','')
        Donor_species.setdefault(donor_species, [])
        for lines in open(vcf,'r'):
            lines_set = lines.replace('\n','').split('\t')
            if not lines.startswith('#'):
                if '*' in lines_set[-1]:
                    # truncated genes
                    gene_name = lines_set[-4]
                    Donor_species[donor_species].append(gene_name)
    return Donor_species

def extract_donor_species(donor_species,input_fasta,gene_name_list,output_fasta):
    gene_name_extract = []
    record_id_set = set()
    donor_species_new = donor_species.replace('_clustercluster','_CL').replace('_PB_','_PaDi_')
    gene_name_list = set(gene_name_list)
    for record in SeqIO.parse(input_fasta, 'fasta'):
        record_id = str(record.id)
        temp_line = '%s\t%s\t%s\t'%(donor_species,donor_species_new,record_id)
        if record_id in gene_name_list and record_id not in record_id_set:
            gene_name_extract.append(record_id)
            record_id = 'C_%s_G_%s' % (record_id.split('_')[1], record_id.split('_')[-1])
            output_fasta.append('>%s__%s\n%s\n' % (donor_species_new, record_id, str(record.seq)))
            record_id_set.add(str(record.id))
            temp_line += '%s__%s\t\n' % (donor_species_new, record_id)
            change_name.add(temp_line)
    if len(gene_name_extract) < len(gene_name_list):
        print('missing genes in files %s extracted %s genes of %s genes'%(input_fasta,len(gene_name_extract),len(gene_name_list)))
        print('missing genes ' + ' '.join([gene_name for gene_name in gene_name_list if gene_name not in gene_name_extract]))
    return output_fasta

def get_all_cluster(High_select2_SNP,record,allrecord):
    for record2 in High_select2_SNP[record]:
        if record2 not in allrecord:
            allrecord.add(record2)
            get_all_cluster(High_select2_SNP, record2, allrecord)
    return allrecord

def get_all_donor(allrecord):
    donorset = set()
    for record in allrecord:
        donorset.add(record.replace('new','').split('donor.')[1].split('__')[0])
    return donorset

def get_all_lineage(allrecord):
    donorset = set()
    for record in allrecord:
        donorset.add(record.split('__')[0])
    return donorset

# clustering
def cluster_uc(cluster_input):
    Clusters = dict()
    High_select2 = set()
    High_select2_SNP = dict()
    High_select2_output = []
    for lines in open(cluster_input, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        cluster = line_set[1]
        record_name = line_set[8]
        record_name2 = line_set[9]
        Clusters.setdefault(record_name, cluster)
        if compute_parallel:
            if record_name2!= '*':
                donor_species = record_name.split('__')[0]
                species = donor_species.split('_')[0]
                donor_species2 = record_name2.split('__')[0]
                species2 = donor_species2.split('_')[0]
                if species == species2:
                    High_select2_SNP.setdefault(record_name, set())
                    High_select2_SNP[record_name].add(record_name2)
                    High_select2_SNP.setdefault(record_name2, set())
                    High_select2_SNP[record_name2].add(record_name)
    if compute_parallel and High_select2_SNP!= dict():
        for record in High_select2_SNP:
            donor_species = record.split('__')[0]
            species = donor_species.split('_')[0]
            allrecord = set()
            allrecord.add(record)
            allrecord = get_all_cluster(High_select2_SNP,record,allrecord)
            if parallel_across_donor:
                donorset = get_all_donor(allrecord)
            else:
                donorset = get_all_lineage(allrecord)
            if len(donorset) >= Donor_species_cutoff.get(species,2):
                # at least 2 or more lineages
                High_select2.add(record)
                print(record,Donor_species_cutoff.get(species,2),species,donorset)
    return [Clusters,High_select2_output,High_select2]

# correct highly selected genes dNdS
def format_freq(freq):
    freq = freq.split(':')
    return [int(freq[0]),int(freq[1])]

def init_highselect_notHS(line_set):
    allspecies_highselect = line_set
    for i in [2,3,5,6,7,8,-13,-11,-9,-7,-5,-3]:
        allspecies_highselect[i] = int(allspecies_highselect[i])
    for i in [-12,-10,-8,-6,-4,-2]:
        allspecies_highselect[i] = format_freq(line_set[i])
    return allspecies_highselect

def init_highselect_empty(line_set):
    allspecies_highselect = line_set
    for i in [3]:
        allspecies_highselect[i] = int(allspecies_highselect[i])
    for i in [2,5,6,7,8,9,10,11,12,-13,-11,-9,-7,-5,-3]:
        allspecies_highselect[i] = 0
    for i in [-12,-10,-8,-6,-4,-2]:
        allspecies_highselect[i] = format_freq('0:0')
    return allspecies_highselect

def add_freq(freq_new,freq_old):
    temp_result = format_freq(freq_new)
    freq_old[0] += temp_result[0]
    freq_old[1] += temp_result[1]

def add_highselect(line_set,allspecies_highselect):
    for i in [2,3,5,6,7,8,-13,-11,-9,-7,-5,-3]:
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
    species = donor_species.split('_')[0]
    if species == '1':
        species = donor_species.split('_')[1]
    genus = species
    if species in Genus_species:
        genus = Genus_species[species]
    return [genus,species]

def add_new_selection(input_summary,High_select2):
    print(len(High_select2))
    Donor_species = dict()
    Donor_species_notHS = dict()
    Genus = dict()
    Genus_notHS = dict()
    Species = dict()
    Species_notHS = dict()
    newoutput = []
    # use selected codon NS ratio (SNP pair) * all genes freq
    for lines in open(input_summary, 'r'):
        if not lines.startswith('#'):
            try:
                line_set = lines.split('\n')[0].split('\t')
                Selected = line_set[-1]
                if Selected == 'False' and len(line_set) > 20:
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
                        genus, species = genus_species(donor_species)
                        if species != donor_species:
                            # avoid double recording SNPs
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
            except ValueError:
                print(lines)
    # calculate NS for new HS genes
    for lines in open(input_summary, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        Selected = line_set[-1]
        donor_species = line_set[0]
        record_id = line_set[1]
        # species or genes
        if Selected == 'False' and len(line_set) > 20:
            if record_id not in ['0', 'gene', 'allspecies', 'allspecies_highselect',
                                 'allspecies_core', 'allspecies_flexible',
                                 'allspecies_highselect_noNS']:
                genus, species = genus_species(donor_species)
                donor_species_new = donor_species.replace('_clustercluster', '_CL').replace('_PB_','_PaDi_')
                if 'NODE' in record_id:
                    record_id = '%s__C_%s_G_%s' % (donor_species_new, record_id.split('_')[1], record_id.split('_')[-1])
                # highly selected genes
                if record_id in High_select2:
                    if line_set[10] == 'observe_None':
                        print(line_set)
                        newoutput.append('\t'.join(line_set[:-1]) + '\tFalse\n')
                        High_select2.remove(record_id)
                    elif (line_set[10] == 'observe_N_only' or float(line_set[10]) >= 0):
                        newoutput.append('\t'.join(line_set[:-1]) + '\tTrueacross\n')
                        if include_across_donor_PE:
                            if len(PE_lineage) == 0 or donor_species in PE_lineage:
                                allspecies_highselect = add_highselect(line_set, allspecies_highselect)
                            try:
                                if len(PE_lineage)==0 or donor_species in PE_lineage:
                                    # significant
                                    if line_set[0] in Donor_species:
                                        Donor_species[line_set[0]] = add_highselect(line_set, Donor_species[line_set[0]])
                                    Genus[genus] = add_highselect(line_set, Genus[genus])
                                    Species[species] = add_highselect(line_set, Species[species])
                            except KeyError:
                                print('no big contig in ',line_set[0])
                    else:
                        newoutput.append('\t'.join(line_set[:-1]) + '\tFalse\n')
                        High_select2.remove(record_id)
                # other genes and species
                elif line_set[0] not in line_set[1]:
                    newoutput.append(lines)
                else:
                    pass
                    #print('not output %s' % ('\t'.join(line_set[0:5])))
        # original highly selected
        elif Selected == 'True' and (len(PE_lineage)==0 or donor_species in PE_lineage):
            newoutput.append(lines)
            allspecies_highselect = add_highselect(line_set, allspecies_highselect)
            if line_set[0] in Donor_species:
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
            if selected == 0 or lines_set[-1] == 'Trueacross':
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
    cutoff = cluster_cutoff
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                   % (args.u, all_filter_gene_fasta_file, cutoff, all_filter_gene_fasta_file,
                      all_filter_gene_fasta_file, 40))
    os.system(cmd_cluster)
    all_filter_gene_fasta_file = all_filter_gene_fasta_file + '.cluster.aa'
    if pre_cluster!= '':
        os.system('#%s -makeudb_usearch %s -output %s.udb' %
                  (args.u, pre_cluster, pre_cluster))
        os.system('%s -ublast %s -db %s.udb  -evalue 1e-2 -accel 0.5 -blast6out %s -threads 2'%
                  (args.u, all_filter_gene_fasta_file,pre_cluster, all_filter_gene_fasta_file + '.ref.out.txt'))
    # run prokka
    cmdsprokka = 'py37\n'+\
    '%s --kingdom Bacteria --force --outdir %s/prokka_%s  --protein %s --locustag Bacter %s/%s\n' % \
                 (args.prokka,output_dir_merge + '/summary', os.path.split(all_filter_gene_fasta_file)[-1],
                  all_filter_gene_fasta_file,
                  output_dir_merge + '/summary',
                  os.path.split(all_filter_gene_fasta_file)[-1].replace('.faa', '.fna'))
    f1 = open(os.path.join(input_script_sub, 'prokka.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmdsprokka))
    f1.close()
    # run metacyc
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/mit_alm/database/metacyc/protseq.fsa'
    cmds = ('py37\n'+\
            "%s blastp --query %s --db %s.dmnd --out %s.metacyc.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    f1 = open(os.path.join(input_script_sub, 'metacyc.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    # run eggnog
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/eggnog/xaa.hmm'
    cmds = ('%s --tblout %s.eggnog.1.txt --cpu 40 -E %s %s %s\n') %(args.hmm, all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.1.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xab.hmm'
    cmds = ('%s --tblout %s.eggnog.2.txt --cpu 40 -E %s %s %s\n') % (
        args.hmm,
        all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.2.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xac.hmm'
    cmds = ('%s --tblout %s.eggnog.3.txt --cpu 40 -E %s %s %s\n') % (
        args.hmm,
        all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.3.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmds))
    f1.close()
    # run kegg
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/kegg/kofam/profiles/prokaryote/prokaryote.hmm'
    cmds = ('%s --tblout %s.kegg.txt --cpu 40 -E %s %s %s\n') %(args.hmm, all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'kegg.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    if annotate_customized:
        # run customed database
        cutoff = 80
        cutoff2 = 80
        cmds = ''
        database = '/scratch/users/anniz44/scripts/database/SARG.db.fasta'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.SARG.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 50
        cutoff2 = 50
        database = '/scratch/users/anniz44/scripts/database/AHR.aa.db'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.AHR.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 60
        cutoff2 = 80
        database = '/scratch/users/anniz44/scripts/database/Butyrate.pro.aa'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.buty.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 50
        cutoff2 = 80
        database = '/scratch/users/anniz44/scripts/database/IntI1_database.fasta'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.int.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 50
        cutoff2 = 80
        database = '/scratch/users/anniz44/scripts/database/SRB.AA'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.SRB.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 0.01
        database = '/scratch/users/anniz44/scripts/database/NR.hmm'
        cmds += ('%s --tblout %s.NR.txt --cpu 40 -E %s %s %s\n') %(args.hmm, all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
        f1 = open(os.path.join(input_script_sub, 'customed.sh'), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
        f1.close()
    # all scripts
    f1 = open(os.path.join(input_script, 'allannotate.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.sh')):
        if 'jobmit' in args.job:
            f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
        else:
            f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()
    print('please run %s/%s'%(input_script,'allannotate.sh'))

def annotation_all(all_filter_gene_fasta_file,pre_cluster = ''):
    # run cluster
    all_filter_gene_fasta_file = all_filter_gene_fasta_file + '.cluster.aa'
    if pre_cluster!= '':
        os.system('#%s -makeudb_usearch %s -output %s.udb' %
                  (args.u, pre_cluster, pre_cluster))
        os.system('%s -ublast %s -db %s.udb  -evalue 1e-2 -accel 0.5 -blast6out %s -threads 2'%
                  (args.u, all_filter_gene_fasta_file,pre_cluster, all_filter_gene_fasta_file + '.ref.out.txt'))
    # run prokka
    cmdsprokka = 'py37\n'+\
            '%s --kingdom Bacteria --force --outdir %s/prokka_%s  --protein %s --locustag Bacter %s/%s\n' % \
                 (args.prokka,output_dir_merge + '/summary', os.path.split(all_filter_gene_fasta_file)[-1],
                  all_filter_gene_fasta_file,
                  output_dir_merge + '/summary',
                  os.path.split(all_filter_gene_fasta_file)[-1].replace('.faa', '.fna'))
    f1 = open(os.path.join(input_script_sub_all, 'prokka.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmdsprokka))
    f1.close()
    # run metacyc
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/mit_alm/database/metacyc/protseq.fsa'
    cmds = ('py37\n'+\
            "%s blastp --query %s --db %s.dmnd --out %s.metacyc.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    f1 = open(os.path.join(input_script_sub_all, 'metacyc.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    # run eggnog
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/eggnog/xaa.hmm'
    cmds = ('%s --tblout %s.eggnog.1.txt --cpu 40 -E %s %s %s\n') %(args.hmm, all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub_all, 'eggnog.1.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xab.hmm'
    cmds = ('%s --tblout %s.eggnog.2.txt --cpu 40 -E %s %s %s\n') % (
        args.hmm,
        all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub_all, 'eggnog.2.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xac.hmm'
    cmds = ('%s --tblout %s.eggnog.3.txt --cpu 40 -E %s %s %s\n') % (
        args.hmm,
        all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub_all, 'eggnog.3.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmds))
    f1.close()
    # run kegg
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/kegg/kofam/profiles/prokaryote/prokaryote.hmm'
    cmds = ('%s --tblout %s.kegg.txt --cpu 40 -E %s %s %s\n') %(args.hmm, all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub_all, 'kegg.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    # run customed database
    if annotate_customized:
        cutoff = 80
        cutoff2 = 80
        cmds = ''
        database = '/scratch/users/anniz44/scripts/database/SARG.db.fasta'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.SARG.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 50
        cutoff2 = 50
        database = '/scratch/users/anniz44/scripts/database/AHR.aa.db'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.AHR.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 60
        cutoff2 = 80
        database = '/scratch/users/anniz44/scripts/database/Butyrate.pro.aa'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.buty.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 50
        cutoff2 = 80
        database = '/scratch/users/anniz44/scripts/database/IntI1_database.fasta'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.int.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 50
        cutoff2 = 80
        database = '/scratch/users/anniz44/scripts/database/SRB.AA'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.SRB.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 0.01
        database = '/scratch/users/anniz44/scripts/database/NR.hmm'
        cmds += ('%s --tblout %s.NR.txt --cpu 40 -E %s %s %s\n') %(args.hmm, all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
        f1 = open(os.path.join(input_script_sub_all, 'customed.sh'), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
        f1.close()
    # all scripts
    f1 = open(os.path.join(input_script, 'allannotate_all.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script_sub_all, '*.sh')):
        if 'jobmit' in args.job:
            f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
        else:
            f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()
    print('please run %s/%s'%(input_script,'allannotate_all.sh'))

def annotation_trunc(all_filter_gene_fasta_file,pre_cluster = ''):
    # run cluster
    cutoff = cluster_cutoff
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                   % (args.u, all_filter_gene_fasta_file, cutoff, all_filter_gene_fasta_file,
                      all_filter_gene_fasta_file, 40))
    os.system(cmd_cluster)
    all_filter_gene_fasta_file = all_filter_gene_fasta_file + '.cluster.aa'
    if pre_cluster!= '':
        os.system('#%s -makeudb_usearch %s -output %s.udb' %
                  (args.u, pre_cluster, pre_cluster))
        os.system('%s -ublast %s -db %s.udb  -evalue 1e-2 -accel 0.5 -blast6out %s -threads 2'%
                  (args.u, all_filter_gene_fasta_file,pre_cluster, all_filter_gene_fasta_file + '.ref.out.txt'))
    # run prokka
    cmdsprokka = 'py37\n'+\
                 '%s --kingdom Bacteria --force --outdir %s/prokka_%s  --protein %s --locustag Bacter %s/%s\n' % \
                 (args.prokka,output_dir_merge + '/summary', os.path.split(all_filter_gene_fasta_file)[-1],
                  all_filter_gene_fasta_file,
                  output_dir_merge + '/summary',
                  os.path.split(all_filter_gene_fasta_file)[-1].replace('.faa', '.fna'))
    f1 = open(os.path.join(input_script_sub_trunc, 'prokka.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmdsprokka))
    f1.close()
    # run metacyc
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/mit_alm/database/metacyc/protseq.fsa'
    cmds = ('py37\n'+\
            "%s blastp --query %s --db %s.dmnd --out %s.metacyc.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    f1 = open(os.path.join(input_script_sub_trunc, 'metacyc.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    # run eggnog
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/eggnog/xaa.hmm'
    cmds = ('%s --tblout %s.eggnog.1.txt --cpu 40 -E %s %s %s\n') %(args.hmm, all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub_trunc, 'eggnog.1.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xab.hmm'
    cmds = ('%s --tblout %s.eggnog.2.txt --cpu 40 -E %s %s %s\n') % (
        args.hmm,
        all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub_trunc, 'eggnog.2.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xac.hmm'
    cmds = ('%s --tblout %s.eggnog.3.txt --cpu 40 -E %s %s %s\n') % (
        args.hmm,
        all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub_trunc, 'eggnog.3.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmds))
    f1.close()
    # run kegg
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/kegg/kofam/profiles/prokaryote/prokaryote.hmm'
    cmds = ('%s --tblout %s.kegg.txt --cpu 40 -E %s %s %s\n') %(args.hmm, all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub_trunc, 'kegg.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n#export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    # run customed database
    if annotate_customized:
        cutoff = 80
        cutoff2 = 80
        cmds = ''
        database = '/scratch/users/anniz44/scripts/database/SARG.db.fasta'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.SARG.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 50
        cutoff2 = 50
        database = '/scratch/users/anniz44/scripts/database/AHR.aa.db'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.AHR.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 60
        cutoff2 = 80
        database = '/scratch/users/anniz44/scripts/database/Butyrate.pro.aa'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.buty.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 50
        cutoff2 = 80
        database = '/scratch/users/anniz44/scripts/database/IntI1_database.fasta'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.int.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 50
        cutoff2 = 80
        database = '/scratch/users/anniz44/scripts/database/SRB.AA'
        cmds += ("%s blastp --query %s --db %s.dmnd --out %s.SRB.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
                %(args.dm,all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
        cutoff = 0.01
        database = '/scratch/users/anniz44/scripts/database/NR.hmm'
        cmds += ('%s --tblout %s.NR.txt --cpu 40 -E %s %s %s\n') %(args.hmm, all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
        f1 = open(os.path.join(input_script_sub_trunc, 'customed.sh'), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s'%(cmds))
        f1.close()
    # all scripts
    f1 = open(os.path.join(input_script, 'allannotate_trunc.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script_sub_trunc, '*.sh')):
        if 'jobmit' in args.job:
            f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
        else:
            f1.write('nohup sh %s > %s.out &\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()
    print('please run %s/%s'%(input_script,'allannotate_trunc.sh'))

def find_database(donor_species):
    try:
        database = glob.glob('%s/%s/%s*%s' % (genome_root, donor_species.split('.donor')[0],
                                            donor_species.split('.donor')[0], ref_filename))[0]
    except IndexError:
        try:
            database = glob.glob('%s/%s/%s*%s' % (
            genome_root.replace('vcf_round2', 'vcf_round1'), donor_species.split('_')[0],donor_species.split('_')[0], ref_filename))[0]
        except IndexError:
            database = glob.glob(
                '/scratch/users/anniz44/genomes/donor_species/WGS_old/WGS/vcf_round1/co-assembly/%s/%s*%s' % (
                donor_species.split('.donor')[0],donor_species.split('.donor')[0], ref_filename))[0]
    return database

def load_siglineages(lineagefile):
    PE_lineage = []
    for lines in open(lineagefile,'r'):
        lines_set = lines.replace('\n', '').split('\t')
        if not lines.startswith('X.donor_species'):
            if lines_set[24] == 'True' :
                PE_lineage.append(lines_set[0])
    return PE_lineage


################################################## Definition ########################################################
# load significant lineages
PE_lineage = []
if args.sig!='None':
    PE_lineage = load_siglineages(args.sig)

# extract highly selected sequences
print('extract high selected genes')
try:
    f1 = open(output_gene,'r')
except IOError:
    output_fasta = []
    output_fasta_dna = []
    change_name = set()
    Donor_species = sum_donor_species(input_summary)
    for donor_species in Donor_species:
        gene_name_list = Donor_species[donor_species]
        if gene_name_list != []:
            print('processing %s' % donor_species)
            database = find_database(donor_species)
            input_fasta = database + '.faa'
            input_fasta_dna = database + '.fna'
            try:
                f1 = open(input_fasta, 'r')
            except FileNotFoundError:
                os.system(args.pro + ' -q -i %s -a %s' % (database, input_fasta))
            try:
                f1 = open(input_fasta_dna, 'r')
            except FileNotFoundError:
                os.system(args.pro + ' -q -i %s -d %s' % (database, input_fasta_dna))
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
print('extract genes with mutations')
output_gene = output_dir_merge + '/summary/all.denovo.gene.faa'
output_gene_dna = output_dir_merge + '/summary/all.denovo.gene.fna'
output_fasta = []
output_fasta_dna = []
try:
    f1 = open(output_gene,'r')
except IOError:
    Donor_species = sum_donor_species(input_summary,0)
    change_name = set()
    for donor_species in Donor_species:
        gene_name_list = Donor_species[donor_species]
        if gene_name_list != []:
            print('processing %s' % donor_species)
            database = find_database(donor_species)
            input_fasta = database + '.faa'
            input_fasta_dna = database + '.fna'
            try:
                f1 = open(input_fasta, 'r')
            except FileNotFoundError:
                os.system('prodigal -q -i %s -a %s' % (database, input_fasta))
            try:
                f1 = open(input_fasta_dna, 'r')
            except FileNotFoundError:
                os.system('prodigal -q -i %s -d %s' % (database, input_fasta_dna))
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
    cutoff = cluster_cutoff
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                       % (args.u, output_gene, cutoff, output_gene,
                          output_gene, args.t))
    os.system(cmd_cluster)

# extract all truncated genes
if args.trunc!= 'False':
    print('extract genes with truncations')
    output_gene = output_dir_merge + '/summary/all.trunc.gene.faa'
    output_gene_dna = output_dir_merge + '/summary/all.trunc.gene.fna'
    output_fasta = []
    output_fasta_dna = []
    Donor_species = dict()
    try:
        f1 = open(output_gene,'r')
    except IOError:
        output_fasta = []
        output_fasta_dna = []
        change_name = set()
        allvcf = glob.glob(output_dir_merge + '/*donor*.raw.vcf.filtered.vcf.final.snp.txt')
        Donor_species = sum_vcf(allvcf)
        for donor_species in Donor_species:
            gene_name_list = Donor_species[donor_species]
            if gene_name_list != []:
                print('processing %s' % donor_species)
                database = find_database(donor_species)
                input_fasta = database + '.faa'
                input_fasta_dna = database + '.fna'
                try:
                    f1 = open(input_fasta, 'r')
                except FileNotFoundError:
                    os.system(args.pro + ' -q -i %s -a %s' % (database, input_fasta))
                try:
                    f1 = open(input_fasta_dna, 'r')
                except FileNotFoundError:
                    os.system(args.pro + ' -q -i %s -d %s' % (database, input_fasta_dna))
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

# find highly selected genes in one species across populations and donors
all_fasta = os.path.join(output_dir_merge + '/summary', 'all.denovo.gene.faa')
all_fasta_HS = os.path.join(output_dir_merge + '/summary', 'all.selected.gene.faa')
all_fasta_trunc = os.path.join(output_dir_merge + '/summary', 'all.trunc.gene.faa')

Genus_species = dict()
if args.meta!= 'None':
    for lines in open(args.meta,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        Genus_species.setdefault(lines_set[-1],lines_set[0].split('_')[0])

Donor_species_cutoff = dict()
if args.cutoff != 'None':
    for lines in open(args.cutoff):
        lines_set = lines.replace('\n', '').replace('\r', '').split('\t')
        species = lines_set[0].split('_')[0]
        cutoff = int(lines_set[1])
        Donor_species_cutoff.setdefault(species,cutoff)

# clustering
Clusters_gene, High_select2_output,High_select2 = cluster_uc(all_fasta + '.uc')

# correcting
High_select2 = add_new_selection(input_summary,High_select2)

# run clustering
if False:
    sum_gene(input_summary,High_select2,1)
    annotation(all_fasta_HS,'')
    #annotation(all_fasta_HS + '.High_select2.faa','')
    annotation_all(all_fasta,'')
    if args.trunc!= 'False':
        annotation_trunc(all_fasta_trunc,'')
################################################### END ########################################################
