################################################## END ########################################################
################################################### SET PATH ########################################################
# compare bowtie call SNPs VS mapper call SNPs
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path to all vcf files",
                      type=str, default='.',
                      metavar='input/')
required.add_argument("-ref",
                      help="path to ref files",
                      type=str, default='.',
                      metavar='ref/')
################################################## Definition ########################################################
args = parser.parse_args()
# set up path
output_dir = args.i
ref_dir = args.ref
# ref
vcf_ref = glob.glob(ref_dir + '/*.snp.txt')
# uncorrected SNPs
SNP_uncorrected = dict()
SNP_uncorrected.setdefault('am_BaFr_gS1T203',['NODE_41_length_7188_cov_24.9891_ID_6233\t189\t\n',
                                                                    'NODE_27_length_32950_cov_21.0268_ID_5538\t3\t\n'])
# function
def load_vcf_ref(vcf_file):
    vcf_input = []
    vcf_ref = dict()
    for lines in open(vcf_file):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = lines_set[1]
            CHRPOS = '%s\t%s\t\n'%(CHR,POS)
            Gene = lines_set[4].replace('N','Gene').replace('S','Gene')
            vcf_input.append(CHRPOS)
            vcf_ref.setdefault(CHRPOS,Gene)
    return [vcf_input,vcf_ref]

def load_vcf(vcf_file):
    vcf_input = []
    for lines in open(vcf_file):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = lines_set[1]
            vcf_input.append('%s\t%s\t\n'%(CHR,POS))
    return vcf_input


def compare_vcf(vcf_input1, vcf_input2, vcf_file1, vcf_file2, vcf_file3,output_file):
    vcf_1_diff = dict()
    vcf_1_diff.setdefault('all',[])
    vcf_1_diff.setdefault('Gene', 0)
    vcf_1_diff.setdefault('Other', 0)
    # CHRPOS in 1 not in 2
    for CHRPOS in vcf_input1:
        if CHRPOS not in vcf_input2:
            vcf_1_diff['all'].append(CHRPOS)
            if CHRPOS in ref_gene:
                vcf_1_diff[ref_gene[CHRPOS]] += 1
    if len(vcf_1_diff['all']) > 0 and len(vcf_1_diff['all']) <= 1000:
        temp_output = os.path.join(output_dir, 'grep.temp.txt')
        f1 = open(temp_output, 'w')
        f1.write(''.join(vcf_1_diff['all']))
        f1.close()
        os.system('grep -T -f %s %s %s %s --no-group-separator > %s' % (
            temp_output,
            vcf_file1, vcf_file2,vcf_file3,
            output_file + '.temp'))
        os.system('sort -k3 -n %s | sort -k2 > %s' %
                  (output_file + '.temp', output_file + '.' + 'all')
                  )
        os.system('rm -rf %s' % (output_file + '.temp'))
    else:
        os.system('rm -rf %s' % (output_file + '.' + 'all'))
    return vcf_1_diff

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def compare_vcf_all(vcf_1,vcf_1_all,vcf_2):
    vcf_input1 = load_vcf(vcf_1)
    FN_diff = compare_vcf(vcf_inputref, vcf_input1, vcf_ref_file, vcf_1_all, vcf_2,vcf_1 + '.FN') # vcf diff in ref not in 1, FN
    FP_diff = compare_vcf(vcf_input1, vcf_inputref, vcf_1, vcf_ref_file, vcf_2, vcf_1 + '.FP') # vcf diff in 1 not in ref, FP
    summary_file_output.append('%s\t%s\t%s\t%s\n'%(
        os.path.split(vcf_1)[-1],
        len(FP_diff['all']),
        len(FN_diff['all']),total_SNP_ref
    ))

def compare_vcf_nogrep(vcf_input1, vcf_input2):
    vcf_1_diff = dict()
    vcf_1_diff.setdefault('all',[])
    vcf_1_diff.setdefault('Gene', 0)
    vcf_1_diff.setdefault('Other', 0)
    # CHRPOS in 1 not in 2
    for CHRPOS in vcf_input1:
        if CHRPOS not in vcf_input2:
            vcf_1_diff['all'].append(CHRPOS)
            if CHRPOS in ref_gene:
                vcf_1_diff[ref_gene[CHRPOS]] += 1
    return vcf_1_diff

def compare_vcf_bowtie(vcf_1):
    vcf_input1 = load_vcf(vcf_1)
    FN_diff = compare_vcf_nogrep(vcf_inputref, vcf_input1) # vcf diff in ref not in 1, FN
    FP_diff = compare_vcf_nogrep(vcf_input1, vcf_inputref) # vcf diff in 1 not in ref, FP
    summary_file_output.append('%s\t%s\t%s\t%s\n'%(
        os.path.split(vcf_1)[-1],
        len(FP_diff['all']),
        len(FN_diff['all']),total_SNP_ref
    ))
# WGS
summary_file = output_dir + '/model.sum.txt'
summary_file_output = []
summary_file_output.append('sample\tFP\tFN\ttotal_refSNP\n')
for vcf_ref_file in vcf_ref:
    vcf_ref_file_name = os.path.split(vcf_ref_file)[-1].split('.snp.txt')[0]
    print('process vcfs for reference %s'%(vcf_ref_file_name))
    try:
        # load ref
        vcf_inputref, ref_gene = load_vcf_ref(vcf_ref_file)
        total_SNP_ref = len(vcf_inputref)
        # add uncorrected SNPs
        genomename = vcf_ref_file_name.split('.')[0]
        if genomename in SNP_uncorrected:
            vcf_inputref += SNP_uncorrected[genomename]
        # bowtie
        vcf_1 = glob.glob(output_dir + vcf_ref_file_name + '.bowtie.flt.snp.vcf.final.vcf')[0]
        #vcf_1_all = vcf_1.replace('.bowtie.flt.snp.vcf.final.vcf','.bowtie.raw.vcf')
        compare_vcf_bowtie(vcf_1)
        # mapper
        vcf_2 = glob.glob(output_dir + vcf_ref_file_name + '.mapper1.vcf.final.vcf')[0]
        vcf_2_all = vcf_2.replace('.mapper1.vcf.final.vcf', '.mapper1.vcf')
        compare_vcf_all(vcf_2,vcf_2_all,vcf_1)
    except IndexError:
        pass

f1=open(summary_file,'w')
f1.write(''.join(summary_file_output))
f1.close()
os.system('rm -rf %s'%(os.path.join(output_dir, 'grep.temp.txt')))
