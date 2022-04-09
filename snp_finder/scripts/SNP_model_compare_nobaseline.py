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
end_cutoff = 5
min_qual_for_call = 20 #Remove sample*candidate that has lower than this quality

# function
def load_vcf_ref(vcf_file):
    vcf_input = []
    vcf_refset = dict()
    for lines in open(vcf_file):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = lines_set[1]
            CHRPOS = '%s\t%s'%(CHR,POS)
            Gene = lines_set[4].replace('N','Gene').replace('S','Gene')
            vcf_input.append(CHRPOS)
            vcf_refset.setdefault(CHRPOS,Gene)
    return [vcf_input,vcf_refset]

def load_vcf(vcf_file):
    vcf_input = []
    for lines in open(vcf_file):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = lines_set[1]
            vcf_input.append('%s\t%s'%(CHR,POS))
    return vcf_input


def FP_end(CHRPOS):
    CHR,POS,nouse = CHRPOS.split('\t')
    if CHR in CHR_length:
        total_length = CHR_length[CHR]
    else:
        try:
            total_length = CHR.split('size')[1]
        except IndexError:
            try:
                total_length = CHR.split('length_')[1].split('_cov')[0]
            except IndexError:
                return False
        total_length = int(total_length)
        CHR_length.setdefault(CHR,total_length)
    if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
        return True
    else:
        return False

def compare_vcf(vcf_input1, vcf_input2, vcf_file1, vcf_file2, vcf_file3,output_file,checkend = False):
    vcf_1_diff = dict()
    vcf_1_diff.setdefault('all',[])
    vcf_1_diff.setdefault('Gene', 0)
    vcf_1_diff.setdefault('Other', 0)
    # CHRPOS in 1 not in 2
    for CHRPOS in vcf_input1:
        if CHRPOS not in vcf_input2:
            # for FP, check end + CHRPOS in CHRPOS_set of bowtie corrected ones
            if not checkend or (not FP_end(CHRPOS) and CHRPOS in CHRPOS_set):
                # check SNP at the ends
                vcf_1_diff['all'].append(CHRPOS)
                if CHRPOS in ref_gene:
                    vcf_1_diff[ref_gene[CHRPOS]] += 1
    if False:
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

def compare_vcf_all(vcf_1, comparetobowtie = False):
    vcf_input1 = load_vcf(vcf_1)
    FN_diff = compare_vcf_nogrep(vcf_inputref, vcf_input1) # vcf diff in ref not in 1, FN
    FP_diff = compare_vcf_nogrep(vcf_input1, vcf_inputref,comparetobowtie) # vcf diff in 1 not in ref, FP
    summary_file_output.append('%s\t%s\t%s\t%s\n'%(
        os.path.split(vcf_1)[-1],
        FP_diff,
        FN_diff,total_SNP_ref
    ))

def compare_vcf_nogrep(vcf_input1, vcf_input2,FP = False):
    vcf_1_diff = 0
    # CHRPOS in 1 not in 2
    for CHRPOS in vcf_input1:
        if CHRPOS not in vcf_input2:
            # for FP, CHRPOS in CHRPOS_set of bowtie corrected ones
            if not FP or (CHRPOS_set!=[] and CHRPOS in CHRPOS_set):
                # check SNP at the ends
                vcf_1_diff += 1
    return vcf_1_diff

def compare_vcf_bowtie(vcf_1, comparetobowtie = False):
    vcf_input1 = load_vcf(vcf_1)
    FN_diff = compare_vcf_nogrep(vcf_inputref, vcf_input1) # vcf diff in ref not in 1, FN
    FP_diff = compare_vcf_nogrep(vcf_input1, vcf_inputref,comparetobowtie) # vcf diff in 1 not in ref, FP
    summary_file_output.append('%s\t%s\t%s\t%s\n'%(
        os.path.split(vcf_1)[-1],
        FP_diff,
        FN_diff,total_SNP_ref
    ))

def load_CHRPOS(vcf_0):
    CHRPOS_set = set()
    i = 0
    for lines in open(vcf_0):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            SNP_quality = float(lines_set[5])
            if SNP_quality >= min_qual_for_call:
                CHR = lines_set[0]
                POS = lines_set[1]
                CHRPOS = '%s\t%s'%(CHR,POS)
                CHRPOS_set.add(CHRPOS)
                i+=1
                if i%1000000 == 0:
                    print(datetime.now(),'input %s lines'%(i))
    return CHRPOS_set

# summarize results
allCHRPOS = dict()
summary_file = output_dir + '/model.sum.txt'
f1=open(summary_file,'w')
f1.write('sample\tFP\tFN\ttotal_refSNP\n')
f1.close()
vcf_ref.sort()
print(vcf_ref)
for vcf_ref_file in vcf_ref:
    vcf_ref_file_name = os.path.split(vcf_ref_file)[-1].split('.snp.txt')[0]
    print('process vcfs for reference %s'%(vcf_ref_file_name))
    # load all POS mapped by bowtie to POS with no SNPs
    SNP_num = vcf_ref_file_name.split('.SNP.fasta')[0].split('.')[-1]
    samplename = vcf_ref_file_name.split('.%s.SNP.fasta' % (SNP_num))[0]
    if samplename not in allCHRPOS:
        vcf_0 = glob.glob(
            output_dir + vcf_ref_file_name.replace('.%s.SNP.fasta' % (SNP_num),
                                                   '.0.SNP.fasta') + '.bowtie.raw.vcf')[0]
        print(datetime.now(), 'load bowtie CHRPOS', vcf_0)
        CHRPOS_set = load_CHRPOS(vcf_0)
        allCHRPOS = dict()
        allCHRPOS.setdefault(samplename, CHRPOS_set)
    CHRPOS_set = allCHRPOS.get(samplename,[])
    comparetobowtie = len(CHRPOS_set)!=0
    # load ref
    vcf_inputref, ref_gene = load_vcf_ref(vcf_ref_file)
    total_SNP_ref = len(vcf_inputref)
    # add uncorrected SNPs
    genomename = vcf_ref_file_name.split('.')[0]
    # store CHR length
    CHR_length = dict()
    try:
        # bowtie
        summary_file_output = []
        vcf_1 = glob.glob(output_dir + vcf_ref_file_name + '.bowtie.flt.snp.vcf.final.vcf')[0]
        print(datetime.now(), 'started processing', vcf_1)
        compare_vcf_bowtie(vcf_1, comparetobowtie)
        f1 = open(summary_file, 'a')
        f1.write(''.join(summary_file_output))
        f1.close()
        print(datetime.now(), 'finished processing', vcf_1)
    except IndexError:
        pass
    try:
        # minimap2
        summary_file_output = []
        vcf_2 = glob.glob(output_dir + vcf_ref_file_name + '.minimap.flt.snp.vcf.final.vcf')[0]
        print(datetime.now(), 'started processing', vcf_2)
        compare_vcf_bowtie(vcf_2, comparetobowtie)
        f1 = open(summary_file, 'a')
        f1.write(''.join(summary_file_output))
        f1.close()
        print(datetime.now(), 'finished processing', vcf_2)
    except IndexError:
        pass
    try:
        # bwa
        summary_file_output = []
        vcf_3 = glob.glob(output_dir + vcf_ref_file_name + '.bwa.flt.snp.vcf.final.vcf')[0]
        print(datetime.now(), 'started processing', vcf_3)
        compare_vcf_bowtie(vcf_3, comparetobowtie)
        f1 = open(summary_file, 'a')
        f1.write(''.join(summary_file_output))
        f1.close()
        print(datetime.now(), 'finished processing', vcf_3)
    except IndexError:
        pass
    try:
        # mapper
        summary_file_output = []
        vcf_4 = glob.glob(output_dir + vcf_ref_file_name + '.mapper1.vcf.final.vcf')[0]
        #vcf_4_all = vcf_4.replace('.mapper1.vcf.final.vcf', '.mapper1.vcf')
        print(datetime.now(), 'started processing', vcf_4)
        compare_vcf_all(vcf_4,comparetobowtie)
        f1 = open(summary_file, 'a')
        f1.write(''.join(summary_file_output))
        f1.close()
        print(datetime.now(), 'finished processing', vcf_4)
    except IndexError:
        pass

#os.system('rm -rf %s'%(os.path.join(output_dir, 'grep.temp.txt')))
