# start
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
# set up path
import argparse
from datetime import datetime
import difflib
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of folders of all genomes",
                      type=str, default='.',
                      metavar='input/')
# optional output setup
optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
# optional search parameters
optional.add_argument('-clustering',
                      help="Optional: 1 for calculating pair-wise SNPs; 2 for clustering genomes based on SNP distance",
                      metavar="1 or 2", action='store', default=1, type=int, choices=[1,2])
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 1)",
                      metavar="1 or more", action='store', default=40, type=int)
optional.add_argument('-job',
                      help="Optional: command to submit jobs",
                      metavar="nohup or customized",
                      action='store', default='jobmit', type=str)
# requirement for software calling
optional.add_argument('-roary',
                          help="Optional: complete path to roary if not in PATH",
                          metavar="/usr/local/bin/roary",
                          action='store', default='roary', type=str)
optional.add_argument('-snp',
                          help="Optional: complete path to snp-sites if not in PATH",
                          metavar="/usr/local/bin/snp-sites",
                          action='store', default='snp-sites', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
genome_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/round1/'
input_script_pan = args.s + '/cluster_CP'
input_script = args.s
vcf_dir = args.i + '/vcf_round1/merge/'
output_dir = args.i + '/vcf_round1/clonal_population'
vcf_final_dir = args.i + '/vcf_round2/merge'
bam_final_dir = args.i + '/vcf_round2/bwa'
empty_case = ('-', 'N')
max_empty = 0.9 # each strain should cover 10% of the genome sequences with SNPs
SNPs_subsample = True # sabsample 0.1% SNPs
empty_species_cutoff = {
'BaOv': 0.4,
    'BL':0.3,
    'EsCo':0.3,
    'FaPr':0.6,
    'BA':0.4
}
core_gene_cutoff = 1000
cluster_cutoff2 = 2 # cutoff for clustering 2 SNP/1kbp
spetial_species = {
    'AlOn':10,
    'BA':3,
    'BaUn':0.5,
    'BaVu':3.0,
    'BL':3,
    'BlWe':5,
    'CoAe':10,
    'EgSp':5,
    'EsCo':5,
    'FaSp':15,
    'FiSp':30,
    'PaDi':3,
    'PaMe':1,
    'StSa':1,
    'TuSa':1.5
}
changename = {
    'bk_BaCa_g0001':'bj_BaCa_g0003',
    'aa_CiAm_g0001':'ao_CiAm_g0015',
    'aa_LaRh_g0001':'av_LaRh_g0010',
'aa_LaRh_g0002':'av_LaRh_g0011',
    'cx_CoSp_g0001':'af_CoSp_g0010',
    'aa_FaPr_g0001':'cx_FaPr_g0010',
'aa_FaPr_g0009':'cx_FaPr_g0011',
'aa_FaPr_g0011':'cx_FaPr_g0012',
'bq_RuTo_g0001':'aa_RuTo_g0010',
    'am_BaUn_g0001':'cx_BaUn_g0010'
}
changedonor = {
"D28":"af",
"D82":"bj",
"D134":"D109",#BL
"P57":"P56",#BL
"P60":"P50",#BL, BA
"P45":"P72",#BA
    "P38":"P06",#BA
    "H23":"H25",#BA
    "P73":"HX",#BA
}
try:
    os.mkdir(args.i + '/vcf_round2/')
except IOError:
    pass

try:
    os.mkdir(vcf_final_dir)
except IOError:
    pass

try:
    os.mkdir(bam_final_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(input_script_pan)
except IOError:
    pass

################################################## Function ########################################################

def SNP_seq(seq1, seq2):
    SNP_total = 0
    total_length = 0
    total_length_all = min(len(seq1),len(seq2))
    if SNPs_subsample and total_length_all >= 50000:
        #subsample 0.1%
        for i in range(0, total_length_all,100):
            if seq1[i] not in empty_case and seq2[i] not in empty_case:
                # excluding empty/missing allele
                total_length += 1
                if seq1[i] != seq2[i]:
                    # a SNP
                    SNP_total += 1
        return [SNP_total,total_length]
    else:
        for i in range(0, total_length_all):
            if seq1[i] not in empty_case and seq2[i] not in empty_case:
                # excluding empty/missing allele
                total_length += 1
                if seq1[i] != seq2[i]:
                    # a SNP excluding empty/missing allele
                    SNP_total += 1
        return [SNP_total,total_length]

def SNP_seq_slow(seq1, seq2):
    s = difflib.SequenceMatcher(None, seq1, seq2).get_matching_blocks()
    SNP_total = 0
    total_length = min(len(seq1),len(seq2))
    sameloci = [range(loci[0],loci[0]+loci[2]) for loci in s if loci[0] == loci[1]]
    diffloci = [loci for loci in range(0, total_length) if loci not in sameloci]
    for loci in diffloci:
        if seq1[loci] not in empty_case and seq2[loci] not in empty_case:
            # a SNP excluding empty/missing allele
            SNP_total += 1
    return SNP_total

def find_neighbor(Cluster_SNP,neighbor,Cluster_SNP_set,cluster,Cluster_SNP_set_added):
    if neighbor != []:
        for record_name in neighbor:
            if record_name not in Cluster_SNP_set_added:
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    subneighbor = Cluster_SNP.get(record_name,[])
                    find_neighbor(Cluster_SNP, subneighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)

def tree_distance(record_before, record_name, tree,SNP_tree_distance):
    try:
        tree_dis = tree.distance(record_before,record_name)
        return [tree_dis/SNP_tree_distance,tree_dis]
    except ValueError:
        if len(record_before) > 8:
            record_before = 'S' + record_before[0:4] + '_' + record_before[-4:]
        else:
            record_before = 'S' + record_before
        if len(record_name) > 8:
            record_name = 'S' + record_name[0:4] + '_' + record_name[-4:]
        else:
            record_name = 'S' + record_name
        try:
            tree_dis = tree.distance(record_before, record_name)
            return [tree_dis / SNP_tree_distance, tree_dis]
        except ValueError:
            print('wrong tree name',record_before,record_name)
            return [0,0]

def find_donor_species(recordname):
    recordnameset = recordname.split('_')
    return '%s_%s'%(recordnameset[0],recordnameset[1])

def changegenomename(recordname,changeanygenome):
    needfix = False
    recordnamenew = recordname
    if recordname in changename:
        # change genome name
        recordnamenew = changename[recordname]
        needfix = True
    #if species not in recordname:
        # add species name
    #    recordnamenew = '%s_%s_%s'%(recordname.split('_')[0],species,recordname.split('_')[1])
    #    needfix = True
    if needfix:
        if recordname not in changeanygenome:
            changeanygenome.add(recordname)
            if recordname.startswith('P') or recordname.startswith('H') or recordname.startswith('D'):
                if recordnamenew.startswith('P') or recordnamenew.startswith('H') or recordnamenew.startswith('D'):
                    olddonor = recordname.split('_')[0]
                    newdonor = recordnamenew.split('_')[0]
                    old_dir = '/scratch/users/mit_alm/IBD_Evo/%s/%s/%s' % (species,olddonor,recordname)
                    new_dir = '/scratch/users/mit_alm/IBD_Evo/%s/%s/%s' % (species,newdonor,recordnamenew)
                    os.system('mv %s %s'%(old_dir,new_dir))
                    print('mv %s %s' % (old_dir, new_dir))
                else:
                    print(recordname,recordnamenew)
            else:
                os.system('mv %s/%s/fastq/%s_1.fastq %s/%s/fastq/%s_1.fastq'%(
                    genome_dir,find_donor_species(recordname),recordname,
                    genome_dir, find_donor_species(recordnamenew), recordnamenew))
                os.system('mv %s/%s/fastq/%s_2.fastq %s/%s/fastq/%s_2.fastq' % (
                    genome_dir, find_donor_species(recordname), recordname,
                    genome_dir, find_donor_species(recordnamenew), recordnamenew))
                os.system('mv %s/%s/fasta/%s_final.scaffolds.fasta %s/%s/fasta/%s_final.scaffolds.fasta' % (
                    genome_dir, find_donor_species(recordname), recordname,
                    genome_dir, find_donor_species(recordnamenew), recordnamenew))
                os.system('mv %s/%s/fasta/%s_prokka.gff %s/%s/fasta/%s_prokka.gff' % (
                    genome_dir, find_donor_species(recordname), recordname,
                    genome_dir, find_donor_species(recordnamenew), recordnamenew))
                os.system('mv %s/%s/fasta/%s_prokka.faa %s/%s/fasta/%s_prokka.faa' % (
                    genome_dir, find_donor_species(recordname), recordname,
                    genome_dir, find_donor_species(recordnamenew), recordnamenew))
                os.system('mv %s/%s/fasta/%s_prokka.faa.add %s/%s/fasta/%s_prokka.faa.add' % (
                    genome_dir, find_donor_species(recordname), recordname,
                    genome_dir, find_donor_species(recordnamenew), recordnamenew))
            if 'WGS_old' not in database_file:
                os.system('mv %s/../bwa/%s_1.fastq.sorted.bam %s/%s_1.fastq.sorted.bam'%(
                    output_dir,recordname,bam_final_dir,recordnamenew))
    return [recordnamenew,changeanygenome]
################################################## Main ########################################################
if args.clustering == 1:
    # calculating pair-wise SNPs
    allsnpfasta = glob.glob(vcf_dir + '/*.all.parsi.fasta.normal.fasta')
    for snpfasta in allsnpfasta:
        species = os.path.split(snpfasta)[-1].split('.all.parsi.fasta.normal.fasta')[0]
        SNP_result_file = os.path.join(output_dir, '%s.SNP.pair.sum.txt' % (species))
        try:
            f1 = open(SNP_result_file, 'r')
        except FileNotFoundError:
            print('start analyzing SNP distribution in %s' % (species))
            # load tree
            try:
                tree = Phylo.read(glob.glob(vcf_dir + '/tree/%s.all.parsi.fasta.out.*tree' % (species))[0], "newick")
            except FileNotFoundError:
                print('no tree for %s in %s/tree/%s.all.parsi.fasta.out.tree' % (species, vcf_dir, species))
                tree = dict()
            # load SNP fasta and calculate pair-wise SNPs
            Ref = ''
            Seq_list = dict()
            SNP_pair = []
            Bad_sample = set()
            #SNP_length = []
            SNP_pair.append('Genome1\tGenome2\tSNPs\tCore_gene_length\n')
            for record in SeqIO.parse(snpfasta, 'fasta'):
                record_name = str(record.id)
                record_seq = str(record.seq)
                total_length = len(record_seq)
                count_empty = record_seq.count('-')
                # SNP_length.append('%s\t%s\t%s\t%s\n'%(record_name,
                #                                          (count_empty/total_length),
                #                      count_empty,
                #                      total_length))
                if count_empty < empty_species_cutoff.get(species, max_empty) * total_length:
                    # less than min_coverage% empty
                    if Ref == '':
                        # set up the first seq as ref
                        Ref = record_seq
                        REF_name = record_name
                        SNP_tree_distance = 1 / len(Ref)
                    else:
                        # print('start analysis one seq', datetime.now())
                        for record_before in Seq_list:
                            SNP_total, total_length = SNP_seq(Seq_list[record_before], record_seq)
                            # SNP_total_tree, tree_distance_2record = tree_distance(record_before, record_name, tree,
                            #                                                      SNP_tree_distance)
                            SNP_pair.append('%s\t%s\t%s\t%s\n' % (
                                record_before, record_name, SNP_total,
                                total_length))
                        print('finish %s' % (record_name))
                    Seq_list.setdefault(record_name, record_seq)
                else:
                    Bad_sample.add(record_name)
                    print(record_name, count_empty / total_length)
            if len(SNP_pair) == 1:
                # no SNPs in core genes
                print('no SNPs found in %s' % (species))
                for record in SeqIO.parse(snpfasta, 'fasta'):
                    record_name = str(record.id)
                    if record_name not in Bad_sample:
                        for record_before in Seq_list:
                            SNP_pair.append(
                                '%s\t%s\t%s\t%s\n' % (record_before, record_name, 0, 0))
                        Seq_list.setdefault(record_name, '')
            f1 = open(SNP_result_file, 'w')
            f1.write(''.join(SNP_pair))
            f1.close()
            #f1 = open('all.genome.empty.length.txt', 'a')
            #f1.write(''.join(SNP_length))
            #f1.close()
            print('finish', datetime.now())

if args.clustering == 2:
    # save all reference genomes
    allref = []
    # clustering genomes based on SNP distance
    all_results = glob.glob(os.path.join(output_dir, 'PaDi.SNP.pair.sum.txt'))
    for SNP_result_file in all_results:
        changeanygenome = set()
        species = os.path.split(SNP_result_file)[-1].split('.SNP.pair.sum.txt')[0]
        try:
            f1 = open(os.path.join(output_dir, '%s.genome.cluster.txt') % (species), 'r')
        except FileNotFoundError:
            # check ref file
            rawvcf = glob.glob('%s/../merge/%s.all*.raw.vcf' % (output_dir, species))
            try:
                for lines in open(rawvcf[0], 'r'):
                    if lines.startswith('##reference=file:'):
                        # set database
                        database_file = lines.split('##reference=file:')[1].split('\n')[0]
                        break
            except IndexError:
                database_file = 'WGS_old'
            # load clustering file
            Cluster_SNP = dict()
            Cluster_SNP_set_added = set()
            Cluster_SNP_set = dict()
            all_record = set()
            Cluster_output = []
            cluster_cutoff = spetial_species.get(species,cluster_cutoff2)
            print(species,cluster_cutoff)
            for lines in open(SNP_result_file, 'r'):
                if not lines.startswith('Genome1'):
                    lines_set = lines.split('\n')[0].split('\t')
                    record_before, record_new, SNP, tree_distance_2record, SNP_curated, total_length = lines_set[0:6]
                    record_before,changeanygenome = changegenomename(record_before,changeanygenome)
                    record_new,changeanygenome = changegenomename(record_new,changeanygenome)
                    Cluster_SNP.setdefault(record_before, [])
                    Cluster_SNP.setdefault(record_new, [])
                    all_record.add(record_before)
                    all_record.add(record_new)
                    if float(SNP_curated)/float(total_length)*1000 <= cluster_cutoff:
                        Cluster_SNP[record_before].append(record_new)
                        Cluster_SNP[record_new].append(record_before)
            cluster = 0
            for record_name in Cluster_SNP:
                neighbor = Cluster_SNP.get(record_name, [])
                if neighbor != [] and record_name not in Cluster_SNP_set_added:
                    cluster += 1
                    Cluster_SNP_set.setdefault(cluster, set())
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    find_neighbor(Cluster_SNP, neighbor, Cluster_SNP_set, cluster, Cluster_SNP_set_added)
            # output single genome
            for record_name in all_record:
                if record_name not in Cluster_SNP_set_added:
                    cluster += 1
                    Cluster_SNP_set.setdefault(cluster, set())
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
            Sub_cluster = len(Cluster_SNP_set)
            for cluster in Cluster_SNP_set:
                record_list = Cluster_SNP_set[cluster]
                record_cluster = 'cluster%s' % (cluster)
                for record_name in record_list:
                    Cluster_output.append(
                        '%s\t%s\t%s\t%s\t\n' % (
                            species, record_name, record_cluster, Sub_cluster))
            if 'WGS_old' not in database_file:
                # no need to re-do co-assembly
                allbam = glob.glob('%s/../bwa/*%s*sorted.bam'%(output_dir,species))
                allref.append('%s\t%s\n' % (species,database_file))
                os.system('mv %s %s/' % (' '.join(allbam), bam_final_dir))
                if len(changeanygenome) == 0: # no changed genome
                    # move to a finished dir
                    fltvcf = glob.glob('%s/../merge/%s.all*.flt.snp.vcf' % (output_dir, species))
                    os.system('mv %s %s/'%(' '.join(rawvcf),vcf_final_dir))
                    os.system('mv %s %s/' % (' '.join(fltvcf), vcf_final_dir))
            f1 = open(os.path.join(output_dir, '%s.genome.cluster.txt') % (species), 'w')
            f1.write('species\tgenome\tcluster\ttotal_cluster\t\n%s' % (''.join(Cluster_output)))
            f1.close()
    f1 = open(os.path.join(output_dir, 'reference.genome.txt'), 'a')
    f1.write(''.join(allref))
    f1.close()