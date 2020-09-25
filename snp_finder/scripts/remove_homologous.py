import glob
import os
import statistics
from Bio import SeqIO
from Bio.Seq import Seq
# set up path
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of genome",
                      type=str, default='.',
                      metavar='input/genome.fasta')
optional.add_argument('-blastn',
                      help="Optional: complete path to blastn if not in PATH",
                      metavar="/usr/local/bin/blastn",
                      action='store', default='blastn', type=str)

############################################ Functions ##############################################
args = parser.parse_args()
# homologous cutoff
length_cutoff = 500 # 500 bp for homologous region
identity_cutoff = 95 # 95% identity for homologous region
CHR_length_cutoff = 2000 # minimum contig lengths for reference genome

def length_CHR(CHR):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        total_length = CHR.split('length_')[1].split('_cov')[0]
    return int(total_length)

def remove_homologous(genome):
    try:
        f1 = open('%s.noHM.fasta' %(genome),'r')
    except FileNotFoundError:
        print('searching homologous regions for genome %s'%(genome))
        try:
            f1 = open('%s.homologous.txt'%(genome),'r')
        except FileNotFoundError:
            os.system(os.path.join(os.path.split('args.blastn')[0], 'makeblastdb') + ' -in %s -dbtype nucl'%(genome))
            os.system(args.blastn + ' -db %s -query %s -word_size %s -perc_identity %s -out %s.homologous.txt -outfmt 6 -max_target_seqs 100' % (
                genome,genome,length_cutoff,identity_cutoff,genome))
        print('removing homologous regions for genome %s' % (genome))
        HM_region = []
        for lines in open('%s.homologous.txt'%(genome),'r'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR1,CHR2 = lines_set[0:2]
            L1, L2 = [length_CHR(CHR1), length_CHR(CHR2)]
            # contig length cutoff
            if L1 < CHR_length_cutoff:
                HM_region.append(CHR1)
            if L2 < CHR_length_cutoff:
                HM_region.append(CHR2)
            if CHR1!=CHR2:
                # homologous region cutoff
                if L1 < L2:
                    HM_region.append(CHR1)
                else:
                    HM_region.append(CHR2)
        Newgenome = []
        for record in SeqIO.parse(genome, 'fasta'):
            record_id = str(record.id)
            if record_id not in HM_region:
                Newgenome.append('>%s\n%s\n'%(record_id, str(record.seq)))
        f1 = open('%s.noHM.fasta' %(genome),'w')
        f1.write(''.join(Newgenome))
        f1.close()
    return '%s.noHM.fasta' %(genome)
############################################ Main ##############################################
remove_homologous(args.i)
