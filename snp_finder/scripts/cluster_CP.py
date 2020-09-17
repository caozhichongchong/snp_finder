# start
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
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
optional.add_argument('-roary', '--roary',
                          help="Optional: complete path to roary if not in PATH",
                          metavar="/usr/local/bin/roary",
                          action='store', default='roary', type=str)


################################################## Definition ########################################################
args = parser.parse_args()

input_script_pan = args.s + '/pangenome'
input_script = args.s
genome_root = args.i
genome_dir = glob.glob('%s/*'%(args.i))
pangenome_dir = args.o + '/pangenome'
output_dir = args.o + '/clonal_population'

core_gene_cutoff = 1000

try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(input_script_pan)
except IOError:
    pass


################################################## Function ########################################################
def pangenome(input_dir,species,output_dir,input_cd):
    cmds = '%s -p %s -o roary_%s -e -n --mafft -f %s/roary_%s -i %s -cd %s %s/*.gff\n' % (
        args.roary,args.t, species, output_dir, species, args.id, input_cd, input_dir)
    cmds += '%s %s/roary_%s/core_gene_alignment.aln > %s/roary_%s/core_gene_alignment.snp_sites.aln\n' %(
        args.snp, output_dir, species, output_dir, species
    )
    return cmds

def checknum(summary_file):
    for lines in open(summary_file,'r'):
        if lines.startswith('Core'):
            lines_set = lines.split('\n')[0].split('\t')
            if int(lines_set[-1]) >= core_gene_cutoff:
                return 'pass'
            else:
                return 're-run'

def SNP_seq(seq1, seq2):
    SNP_total = 0
    j = 0
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
    return SNP_total

################################################## Main ########################################################
roary_folder = glob.glob(pangenome_dir + '/roary_*')
for roary_folder1 in roary_folder:
    species = os.path.split(roary_folder1)[-1].split('roary_')[1]
    core_gene_alignment = os.path.join(roary_folder1,'core_gene_alignment.snp_sites.aln')
    core_gene_alignment_snp = os.path.join(roary_folder1,'core_gene_alignment.snp_sites.aln')
    if checknum(os.path.join(roary_folder1,'summary_statistics.txt')) == 're-run':
        print('%s has < 1000 core genes. Sorry but we have to re-run the roary with core cutoff of 96'%(species))
        cmds = pangenome(genome_root, species, roary_folder1 + '_cd95', 95)
        print('please run:\n%s'%cmds)
    else:
        # load core genome alignment SNP fasta and calculate pair-wise SNPs
        Ref = ''
        Seq_list = dict()
        SNP_pair = []
        SNP_pair.append('Genome1\tGenome2\tSNPs\tCore_gene_length\t\n')
        Total_length = 0
        # load total length
        for record in SeqIO.parse(core_gene_alignment, 'fasta'):
            Total_length = len(str(record.seq))
            break
        for record in SeqIO.parse(core_gene_alignment_snp, 'fasta'):
            record_name = str(record.id)
            record_seq = str(record.seq)
            if Ref == '':
                # set upt the first seq as ref
                Ref = record_seq
                REF_name = record_name
            else:
                for record_before in Seq_list:
                    SNP_total = SNP_seq(Seq_list[record_before], record_seq)
                    SNP_pair.append('%s\t%s\t%s\t%s\t\n' % (record_before, record_name, SNP_total, Total_length))
            Seq_list.setdefault(record_name, record_seq)
        f1 = open(os.path.join(output_dir, '%s.SNP.pair.sum.txt' % (species)), 'a')
        f1.write(''.join(SNP_pair))
        f1.close()