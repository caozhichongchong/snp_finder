# start
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
# set up path
import argparse
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
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='pan-genome/',
                      metavar='pan-genome/')
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
input_script_pan = args.s + '/pangenome'
input_script = args.s
genome_root = args.i
genome_dir = glob.glob('%s/*'%(args.i))
pangenome_dir = args.o + '/pangenome'
output_dir = args.o + '/clonal_population'

core_gene_cutoff = 1000
cluster_cutoff2 = 0 # cutoff for clustering 0 SNPs
spetial_species = {
    '1_BA_IBD_0':1000,
    '1_BA_IBD_1':1000,
    '1_BA_IBD_2':0,
    '1_BA_IBD_3':400,
    '1_BA_IBD_4':0,
'1_BA_IBD_5':0,
'1_BL_IBD_0':2,
'1_BL_IBD_1':0,
'1_BL_IBD_2':0,
'1_BL_IBD_3':0,
#'1_BL_IBD_4':0,
'1_PB_IBD_0':0,
#'1_PB_IBD_1':0,
#'1_PB_IBD_3':100,
#'AkMu': 100,
#'AkSp': 110,
    #'AlOn':0,
    'BA':0,
    #'BaFr':0,
#'BaOv':0,
#?'BaSa':100,
#'BaSt':500,
    #'BaTh':0,
#?'BaVu':0,
  #'BaXy':400,
   # 'BiBi':400,
    'BiPs':0,
    'BL':0,
    'BlWe':0,
    #'CiAm':60,
    #'ClBe':50,
    #'ClIn':15,
    'CoAe':100,
    #'EnDu':500,
    #'EnHi':15,
    #'EnMu':25,
    'EsCo_95':0,
    #'FaPr':10000,
    #'LaRu':40,
    #'PaDi':0,
    'PaEx':110,
    'PaGo':20,
    'PaMe':220,
    #'PhFa':6,
    #'PsSp':200,
    'RuLa':100,
    #'StPa':0,
    'TuSa':50,
    '2_BA_IBD_0':12,
'2_BA_IBD_3':4000, # need check
    '2_BA_IBD_4':60,
'2_BL_IBD_0':1000,
'2_BL_IBD_1':120,
'2_BL_IBD_2':60,
'2_BL_IBD_3':1000,
    '2_PB_IBD_1':600,
    '2_PB_IBD_4':5000,
    '3_BL_IBD_0':2000,
    '3_BL_IBD_1':3000,
    '3_BL_IBD_2':200,
    '3_PB_IBD_0':60,
    '3_PB_IBD_2':200,
    '3_PB_IBD_3':500
}
try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(input_script_pan)
except IOError:
    pass

################################################## Function ########################################################
def checknum(summary_file):
    for lines in open(summary_file,'r'):
        if lines.startswith('Core'):
            lines_set = lines.split('\n')[0].split('\t')
            if int(lines_set[-1]) >= core_gene_cutoff:
                return 'pass'
            else:
                return 're-run'

def SNP_seq(seq1, seq2, empty_keep = True):
    SNP_total = 0
    j = 0
    total_length = len(seq1)
    for i in range(0, total_length):
        try:
            if seq1[i] != seq2[i]:
                if empty_keep:
                    # a SNP considering empty/missing allele
                    SNP_total += 1
                elif seq1[i] != '-' and seq2[i] !='-':
                    # a SNP excluding empty/missing allele
                        SNP_total += 1
        except IndexError:
            print('wrong length', len(seq1),len(seq2),i)
    return SNP_total

def pangenome(input_dir,species,output_dir,cd_cutoff):
    cmds = '%s -p %s -o roary_%s_%s -e -n --mafft -f %s/roary_%s_%s -i %s -cd %s %s/%s/*.gff\n' % (
        args.roary,args.t, species,cd_cutoff, output_dir, species, cd_cutoff, 90, cd_cutoff, input_dir,species)
    cmds += '%s %s/roary_%s/core_gene_alignment.aln > %s/roary_%s/core_gene_alignment.snp_sites.aln\n' %(
        args.snp, output_dir, species, output_dir, species)
    return cmds

def find_neighbor(Cluster_SNP,neighbor,Cluster_SNP_set,cluster,Cluster_SNP_set_added):
    if neighbor != []:
        for record_name in neighbor:
            if record_name not in Cluster_SNP_set_added:
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    subneighbor = Cluster_SNP.get(record_name,[])
                    find_neighbor(Cluster_SNP, subneighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)

def correct_empty(fasta):
    newout = []
    for record in SeqIO.parse(fasta, 'fasta'):
        record_seq = str(record.seq)
        record_seq = record_seq.replace('-','N')
        newout.append('>%s\n%s\n'%(str(record.id),record_seq))
    f1 = open(fasta + '.noempty.aln', 'w')
    f1.write(''.join(newout))
    f1.close()
    print('finished correction empty alleles in %s'%(fasta))

def parsi_tree(core_gene_alignment_snp,tree_output):
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    print('output parsi fasta')
    for record in SeqIO.parse(core_gene_alignment_snp, 'fasta'):
        record_seq = str(record.seq)
        genome = str(record.id)
        seq_len = len(record_seq)
        newgenomename = genome
        if len(genome) > 8:
            newgenomename = genome[0:4] + '_' + genome[-4:]
        if seq_len > 0:
            if seq_len > 1000:
                keep_SNP = 200
                seq_len = keep_SNP
                interval_SNP = int(seq_len / keep_SNP)
                subset = list(range(1, int(seq_len / interval_SNP))) * interval_SNP
                seq_len = len(subset)
                seq_subset = ''.join([record_seq[i] for i in subset])
                SNP_alignment_output_parsi.append('S%s    %s\n' % (newgenomename,
                                                                   seq_subset))
            else:
                SNP_alignment_output_parsi.append('S%s    %s\n' % (newgenomename, record_seq))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    parsi_output = open(core_gene_alignment_snp + '.parsi.fasta', 'w')
    parsi_output.write(temp_line + ''.join(SNP_alignment_output_parsi))
    parsi_output.close()
    print('compute parsi tree')
    core_gene_alignment_snp = core_gene_alignment_snp + '.parsi.fasta'
    SNP_tree_cmd3 = ('%s\n5\nV\n1\ny\n' % (core_gene_alignment_snp))
    f1 = open(os.path.join(input_script_pan, 'parsi.optionfile.txt'), 'w')
    f1.write(SNP_tree_cmd3)
    f1.close()
    os.system('rm -rf outfile outtree')
    os.system('dnapars < %s/parsi.optionfile.txt > %s/%s.parsi.output\n' % (
        input_script_pan, input_script_pan, os.path.split(core_gene_alignment_snp)[-1]))
    os.system('mv outfile %s' % (tree_output + '.txt'))
    os.system('mv outtree %s' % (tree_output))

def tree_distance(record_before, record_name, tree,SNP_tree_distance):
    if len(record_before) > 8:
        record_before = 'S' + record_before[0:4] + '_' + record_before[-4:]
    else:
        record_before = 'S' + record_before
    if len(record_name) > 8:
        record_name = 'S' + record_name[0:4] + '_' + record_name[-4:]
    else:
        record_name = 'S' + record_name
    try:
        return [tree.distance(record_before,record_name)/SNP_tree_distance,tree.distance(record_before,record_name)]
    except ValueError:
        print('wrong tree name',record_before,record_name)
        return [0,0]

################################################## Main ########################################################
if args.clustering == 1:
    # calculating pair-wise SNPs
    roary_folder = glob.glob(pangenome_dir + '/roary_*')
    for roary_folder1 in roary_folder:
        species = os.path.split(roary_folder1)[-1].split('roary_')[1]
        try:
            SNP_result_file = os.path.join(output_dir, '%s.SNP.pair.sum.txt' % (species))
            f1 = open(SNP_result_file,'r')
        except FileNotFoundError:
            print('start analyzing SNP distribution in %s'%(species))
            core_gene_alignment = os.path.join(roary_folder1,'core_gene_alignment.aln')
            try:
                f1 = open(core_gene_alignment,'r')
                try:
                    core_gene_alignment_snp = os.path.join(roary_folder1,'core_gene_alignment.snp_sites.aln')
                    f1 = open(core_gene_alignment_snp, 'r')
                except FileNotFoundError:
                    cmd = '%s -c %s > %s\n' % (
                        args.snp, core_gene_alignment, core_gene_alignment_snp)
                    os.system(cmd)
                if checknum(os.path.join(roary_folder1,'summary_statistics.txt')) == 're-run':
                    print('%s has < 1000 core genes. Sorry but we have to re-run the roary with core cutoff of 95'%(species))
                    cmds = pangenome(genome_root, species, roary_folder1 + '_cd95', 95)
                    print('please run:\n%s'%cmds)
                    os.system('rm -rf %s'%(roary_folder1))
                else:
                    # run tree
                    tree_output = '%s/%s.tree' % (output_dir, species)
                    try:
                        f1 = open(tree_output, 'r')
                    except FileNotFoundError:
                        parsi_tree(core_gene_alignment_snp, tree_output)
                    try:
                        tree = Phylo.read(tree_output, "newick")
                    except FileNotFoundError:
                        tree = dict()
                    # load core genome alignment SNP fasta and calculate pair-wise SNPs
                    Ref = ''
                    Seq_list = dict()
                    SNP_pair = []
                    SNP_pair.append('Genome1\tGenome2\tSNPs\ttree_distance\tSNPs_curated_bytree\tCore_gene_length\t\n')
                    Total_length = 0
                    # load total length
                    for record in SeqIO.parse(core_gene_alignment, 'fasta'):
                        Total_length = len(str(record.seq))
                        break
                    for record in SeqIO.parse(core_gene_alignment_snp, 'fasta'):
                        record_name = str(record.id)
                        record_seq = str(record.seq)
                        if Ref == '':
                            # set up the first seq as ref
                            Ref = record_seq
                            REF_name = record_name
                            SNP_tree_distance = 1 / len(Ref)
                        else:
                            for record_before in Seq_list:
                                if 'IBD' in species:
                                    SNP_total = 0
                                else:
                                    SNP_total = SNP_seq(Seq_list[record_before], record_seq, False)
                                SNP_total_tree, tree_distance_2record = tree_distance(record_before, record_name, tree, SNP_tree_distance)
                                SNP_pair.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (record_before, record_name, SNP_total,tree_distance_2record, SNP_total_tree,Total_length))
                        Seq_list.setdefault(record_name, record_seq)
                    if len(SNP_pair) == 1:
                        # no SNPs in core genes
                        print('no SNPs found in %s'%(species))
                        for record in SeqIO.parse(core_gene_alignment, 'fasta'):
                            record_name = str(record.id)
                            for record_before in Seq_list:
                                SNP_pair.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (record_before, record_name, 0,0,0, Total_length))
                            Seq_list.setdefault(record_name, '')
                    f1 = open(SNP_result_file, 'w')
                    f1.write(''.join(SNP_pair))
                    f1.close()
            except FileNotFoundError:
                print('core alignment for %s does not exist'%(species))

if args.clustering == 2:
    #clustering genomes based on SNP distance
    all_results = glob.glob(os.path.join(output_dir, '*.SNP.pair.sum.txt'))
    for SNP_result_file in all_results:
        species = os.path.split(SNP_result_file)[-1].split('.SNP.pair.sum.txt')[0]
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
                Cluster_SNP.setdefault(record_before, [])
                Cluster_SNP.setdefault(record_new, [])
                all_record.add(record_before)
                all_record.add(record_new)
                if int(float(SNP_curated)) <= cluster_cutoff:
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
        f1 = open(os.path.join(output_dir, '%s.genome.cluster.txt') % (species), 'w')
        f1.write('species\tgenome\tcluster\ttotal_cluster\t\n%s' % (''.join(Cluster_output)))
        f1.close()