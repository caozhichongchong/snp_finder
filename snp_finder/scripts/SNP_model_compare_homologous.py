################################################## END ########################################################
################################################### SET PATH ########################################################
# compare bowtie call SNPs VS mapper call SNPs
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import argparse
from datetime import datetime
import pandas as pd

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path to all vcf files",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/SNP_model_newnew/',
                      metavar='input/')
required.add_argument("-duplication",
                      help="reference regions with duplication",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/test_data/new',
                      metavar='database')
# optional input genome
optional.add_argument("-cluster",
                      help="a cluster to run, default is all clusters",
                      type=str, default='',
                      metavar='cluster1')


################################################## Definition ########################################################
args = parser.parse_args()
# set up path
mapper_checkFP = True
kmer_size = 20 + 10 # kmer + length of indels
# function
def load_vcf_ref(vcf_file):
    vcf_input = set()
    for lines in open(vcf_file):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = lines_set[1]
            CHRPOS = '%s\t%s\t'%(CHR,POS)
            vcf_input.add(CHRPOS)
    return vcf_input

def load_vcf(vcf_file):
    vcf_input = set()
    for lines in open(vcf_file):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = lines_set[1]
            vcf_input.add('%s\t%s\t'%(CHR,POS))
    return vcf_input

def compare_vcf(vcf,duplication_result,print_vcf_set=False):
    summary_file_output = []
    print(datetime.now(), 'started processing', vcf)
    vcf_input1 = load_vcf(vcf)
    # duplication region
    vcf_input1_duplication = [x for x in vcf_input1 if x in duplication_result]
    vcf_inputref_duplication = [x for x in vcf_inputref if x in duplication_result]
    FN_diff_duplication = compare_vcf_detail(vcf, vcf_inputref_duplication,vcf_input1_duplication , False)  # vcf diff in ref not in 1, FN
    FP_diff_duplication = compare_vcf_detail(vcf, vcf_input1_duplication, vcf_inputref_duplication, False)  # vcf diff in 1 not in ref, FP
    # nonduplication region
    vcf_input1_noduplication = [x for x in vcf_input1 if x not in duplication_result]
    vcf_inputref_noduplication = [x for x in vcf_inputref if x not in duplication_result]
    print('scaffold20|size54746\t31899\t' in vcf_inputref_noduplication)
    print('scaffold20|size54746\t31920\t' in vcf_inputref_noduplication)
    print('scaffold34|size40049\t463\t' in vcf_inputref_noduplication)
    print('scaffold30|size44032\t223\t' in vcf_inputref_noduplication)
    FN_diff = compare_vcf_detail(vcf, vcf_inputref_noduplication, vcf_input1_noduplication,False) # vcf diff in ref not in 1, FN
    FP_diff = compare_vcf_detail(vcf, vcf_input1_noduplication, vcf_inputref_noduplication,print_vcf_set) # vcf diff in 1 not in ref, FP
    summary_file_output.append('%s\t%s\t%s\t%s\t%s\t%s\n'%(
        os.path.split(vcf)[-1],
        FP_diff,
        FN_diff,FP_diff_duplication,FN_diff_duplication,total_SNP_ref
    ))
    f1 = open(summary_file, 'a')
    f1.write(''.join(summary_file_output))
    f1.close()
    print(datetime.now(), 'finished processing', vcf)

def compare_vcf_detail(vcf, vcf_input1, vcf_input2, print_vcf_set):
    # CHRPOS in 1 not in 2
    vcf_set = [CHRPOS for CHRPOS in vcf_input1 if CHRPOS not in vcf_input2]
    if print_vcf_set:
        # print max 100 FPs
        temp_output = os.path.join(output_dir, 'grep.temp.txt')
        f1 = open(temp_output, 'w')
        f1.write('\n'.join(vcf_set))
        f1.close()
        os.system('grep -T -f %s %s --no-group-separator > %s' % (
            temp_output,
            vcf,
            vcf + '.temp'))
        os.system('sort -k3 -n %s | sort -k2 > %s' %
                  (vcf + '.temp', vcf + '.' + 'FP')
                  )
        os.system('rm -rf %s' % (vcf + '.temp'))
        os.system('rm -rf %s' % (temp_output))
    return len(vcf_set)

def load_duplication(duplication_file):
    duplication_result = pd.read_csv(duplication_file, sep='\t')
    duplication_result['CHRPOS']=[('%s\t%s\t')%(x,y) for x,y in zip(duplication_result['CHR'],
                                                                    duplication_result['POS'])]
    duplication_list = list(duplication_result['CHRPOS'])
    for i in range(-kmer_size,kmer_size + 1):
        duplication_list += [('%s\t%s\t')%(x,y + i) for x,y in zip(duplication_result['CHR'],
                                                                    duplication_result['POS'])]
    return set(duplication_list)

# ref
allresultdir = glob.glob('%s/SNP_model*/'%(args.i))
allduplication_set = dict()
for resultdir in allresultdir:
    resultdirfilename = os.path.split(resultdir)[-1]
    output_dir = '%s/merge/'%(resultdir)
    ref_dir = '%s/data/'%(resultdir)
    vcf_ref = glob.glob(ref_dir + '/%s*.snp.txt'%(args.cluster))
    print(vcf_ref)
    # summarize results
    summary_file = output_dir + '/model.sum.duplicate.%s.txt'%(resultdirfilename)
    f1=open(summary_file,'w')
    f1.write('sample\tFP_nonduplication\tFN_nonduplication\tFP_duplication\tFN_duplication\ttotal_refSNP\n')
    f1.close()
    vcf_ref.sort()
    for vcf_ref_file in vcf_ref:
        vcf_ref_file_name = os.path.split(vcf_ref_file)[-1].split('.snp.txt')[0]
        print('process vcfs for reference %s'%(vcf_ref_file_name))
        database_name = vcf_ref_file_name.split('.fasta')[0]
        # load duplication
        if database_name not in allduplication_set:
            duplication_file = '%s/%s.duplicate.txt'%(args.duplication,database_name)
            print('loading %s' % (duplication_file))
            duplication_result = load_duplication(duplication_file)
            print('scaffold20|size54746\t31899\t' in  duplication_result)
            print('scaffold20|size54746\t31920\t' in duplication_result)
            print('scaffold34|size40049\t463\t' in duplication_result)
            print('scaffold30|size44032\t223\t' in duplication_result)
            allduplication_set.setdefault(database_name,duplication_result)
        else:
            duplication_result = allduplication_set[database_name]
        # load ref
        vcf_inputref = load_vcf_ref(vcf_ref_file)
        total_SNP_ref = len(vcf_inputref)
        # try:
        #     # bowtie
        #     compare_vcf(glob.glob(output_dir + vcf_ref_file_name + '.bowtie.flt.snp.vcf.final.vcf')[0],duplication_result)
        # except IndexError:
        #     pass
        # try:
        #     # minimap2
        #     compare_vcf(glob.glob(output_dir + vcf_ref_file_name + '.minimap.flt.snp.vcf.final.vcf')[0],duplication_result)
        # except IndexError:
        #     pass
        # try:
        #     # bwa
        #     compare_vcf(glob.glob(output_dir + vcf_ref_file_name + '.bwa.flt.snp.vcf.final.vcf')[0],duplication_result)
        # except IndexError:
        #     pass
        try:
            # mapper
            compare_vcf(glob.glob(output_dir + vcf_ref_file_name + '.mapper1.vcf.final.vcf')[0],duplication_result, mapper_checkFP)
        except IndexError:
            pass
        # try:
        #     # mapper sam
        #     compare_vcf(glob.glob(output_dir + vcf_ref_file_name + '.mapper1sam.flt.snp.vcf.final.vcf')[0],duplication_result, mapper_checkFP)
        # except IndexError:
        #     pass
        # try:
        #     # mapper sam
        #     compare_vcf(glob.glob(output_dir + vcf_ref_file_name + '.kmermapper1.vcf.final.vcf')[0],duplication_result, mapper_checkFP)
        # except IndexError:
        #     pass