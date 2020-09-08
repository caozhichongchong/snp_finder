# start
# compare genome call SNPs VS WGS call SNPs
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import random
import argparse

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="subsample * 20", type=int,
                    default=0,metavar='0 to 11')
################################################## Definition ########################################################
args = parser.parse_args()
i = int(args.i)
# set up path
output_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round4/merge/'
ref_dir = '/scratch/users/anniz44/genomes/donor_species/jay/vcf_round4/merge_genome/'
# genome
vcf_ref = glob.glob(ref_dir + '/*.all.flt.snp.vcf.filtered.snp.txt')
# set up cutoff
end_cutoff = 70 # 10 bp at the ends of a contig, separate

# function
def load_vcf_ref(vcf_file,end_count = False):
    end_set = ['end', 'notend']
    vcf_count = dict()
    vcf_count.setdefault('end', 0)
    vcf_count.setdefault('notend', 0)
    vcf_input = []
    vcf_ref = dict()
    for lines in open(vcf_file):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = lines_set[1]
            CHRPOS = '%s\t%s\t\n'%(CHR,POS)
            Gene = lines_set[-2]
            if Gene == 'None':
                Gene = Gene.replace('None','Other')
            else:
                Gene = 'Gene'
            vcf_input.append(CHRPOS)
            vcf_ref.setdefault(CHRPOS,Gene)
            if end_count:
                vcf_count[contig_end(CHR,POS)] += 1
    if end_count:
        summary_file_output.append('%s\t%s\t%s\t0\t0\n'%(os.path.split(vcf_file)[-1],vcf_count['end'],
                                                               vcf_count['notend']))
        print(vcf_count)
    return [vcf_input,vcf_ref]

def load_vcf(vcf_file,end_count = False):
    end_set = ['end', 'notend']
    vcf_count = dict()
    vcf_count.setdefault('end', 0)
    vcf_count.setdefault('notend', 0)
    vcf_input = []
    for lines in open(vcf_file):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = lines_set[1]
            vcf_input.append('%s\t%s\t\n'%(CHR,POS))
            if end_count:
                vcf_count[contig_end(CHR,POS)] += 1
    if end_count:
        summary_file_output.append('%s\t%s\t%s\t0\t0\n'%(os.path.split(vcf_file)[-1],vcf_count['end'],
                                                               vcf_count['notend']))
        print(vcf_count)
    return vcf_input

def contig_end(CHR,POS):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        total_length = CHR.split('length_')[1].split('_cov')[0]
    total_length = int(total_length)
    if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
        return 'end'
    else:
        return 'notend'

def compare_vcf(vcf_input1, vcf_input2, vcf_file1, vcf_file2, output_file):
    #end_set = ['end','notend']
    end_set = ['notend']
    vcf_1_diff = dict()
    vcf_1_diff.setdefault('end',[])
    vcf_1_diff.setdefault('notend', [])
    vcf_1_diff.setdefault('Gene', 0)
    vcf_1_diff.setdefault('Other', 0)
    # CHRPOS in 1 not in 2
    for CHRPOS in vcf_input1:
        if CHRPOS not in vcf_input2:
            CHR,POS = CHRPOS.split('\t')[0:2]
            contig_end_tag = contig_end(CHR,POS)
            vcf_1_diff[contig_end_tag].append(CHRPOS)
            if contig_end_tag == 'notend' and CHRPOS in ref_gene:
                vcf_1_diff[ref_gene[CHRPOS]] += 1
    for contig_end_tag in end_set:
        if len(vcf_1_diff[contig_end_tag]) > 0:
            print(len(vcf_1_diff[contig_end_tag]))
            temp_output = os.path.join(output_dir, 'grep.temp.%s.txt' % (contig_end_tag))
            print(temp_output)
            f1 = open(temp_output, 'w')
            f1.write(''.join(vcf_1_diff[contig_end_tag]))
            f1.close()
            #os.system('grep -T -f %s %s %s --no-group-separator > %s' % (
            #    temp_output,
            #    vcf_file1,vcf_file2,
            #    output_file + '.temp'))
            #os.system('sort -k3 -n %s | sort -k2 > %s' %
            #          (output_file + '.temp', output_file + '.' + contig_end_tag)
            #          )
            #os.system('rm -rf %s' % (output_file+ '.temp'))
        else:
            os.system('rm -rf %s'%(output_file + '.' + contig_end_tag))
    return vcf_1_diff

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def compare_vcf_all(vcf_1,vcf_2 = ''):
    vcf_input1 = load_vcf(vcf_1)
    vcf_1_all = glob.glob(vcf_1.split('.flt.snp.vcf')[0] + '.raw.vcf')
    if vcf_1_all != []:
        vcf_1_all = vcf_1_all[0]
    else:
        vcf_1_all = vcf_1
    if vcf_2 != '':
        vcf_input1 = intersection(vcf_input1, load_vcf(vcf_2))
        vcf_input1 = list(set(vcf_input1))
        vcf_1 = vcf_1 + '.all'
    FN_diff = compare_vcf(vcf_inputref, vcf_input1, vcf_ref_file, vcf_1_all, vcf_1 + '.FN') # vcf diff in ref not in 1, FN
    FP_diff = compare_vcf(vcf_input1, vcf_inputref, vcf_1, vcf_ref_file, vcf_1 + '.FP') # vcf diff in 1 not in ref, FP
    summary_file_output.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(
        os.path.split(vcf_1)[-1],
        len(FN_diff['end']),len(FN_diff['notend']),
        len(FP_diff['end']), len(FP_diff['notend']),
                                 FN_diff['Gene'],FN_diff['Other']
    ))

def depthcheck(vcf_1,vcf_2,vcf_3):
    vcf_input2 = load_vcf(vcf_2)
    vcf_input3 = load_vcf(vcf_3)
    Length = dict()
    Total = 0
    for lines in open(vcf_1):
        if not lines.startswith('#'):
            lines_set = lines.split('\n')[0].split('\t')
            Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
            if Total == 0:
                Total = len(lines_set) - 9
            if Depth / Total <= 100:
                CHR = lines_set[0]
                if CHR not in Length:
                    try:
                        total_length = CHR.split('size')[1]
                    except IndexError:
                        total_length = CHR.split('length_')[1].split('_cov')[0]
                    total_length = int(total_length)
                    Length.setdefault(CHR,total_length)
                total_length = Length[CHR]
                if total_length >= 5000:
                    POS = lines_set[1]
                    CHRPOS = '%s\t%s\t\n' % (CHR, POS)
                    lines_set_sub = lines_set[9:]
                    for Subdepth_all in lines_set_sub:
                        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                        total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                        # Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
                        # Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
                        # forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
                        # reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
                        if total_sub_depth > 0:
                            Depth_set.setdefault(total_sub_depth, [0, 0, 0, 0])
                            if CHRPOS in vcf_inputref:
                                Depth_set[total_sub_depth][1] += 1
                                if CHRPOS in vcf_input3:
                                    Depth_set[total_sub_depth][3] += 1
                            else:
                                Depth_set[total_sub_depth][0] += 1
                            if CHRPOS in vcf_input2:
                                Depth_set[total_sub_depth][2] += 1
    return 'done'

# WGS VS Genome
summary_file = output_dir + '/diff.sum.nofilter.txt'
summary_file_output = []
summary_file_output.append('sample\tFN_end\tFN_notend\tFP_end\tFP_notend\tFN_notend_gene\tFN_notend_other\n')
for vcf_ref_file in vcf_ref:
    vcf_ref_file_name = os.path.split(vcf_ref_file)[-1].split('.all.')[0]
    try:
        vcf_inputref, ref_gene = load_vcf_ref(vcf_ref_file,True)
        # WGS filtered
        vcf_1 = glob.glob(output_dir + vcf_ref_file_name + '*.all.fq.flt.snp.vcf.filtered.snp.txt')[0]
        compare_vcf_all(vcf_1)
        # WGS not filtered
        vcf_1 = glob.glob(output_dir + vcf_ref_file_name + '*.all.fq.flt.snp.vcf')[0]
        compare_vcf_all(vcf_1)
    except IndexError:
        pass

f1=open(summary_file,'w')
f1.write(''.join(summary_file_output))
f1.close()
os.system('rm -rf %s'%(os.path.join(output_dir, 'grep.temp.*.txt')))

# Depth check of Genome mutations
summary_file = output_dir + '/SNP.depth.sum.nofilter%s.txt' %(i)
f1=open(summary_file,'w')
f1.write('Depth\tGe_noSNP\tGe_SNP\tWGS_filterSNP\tWGS_Ge_SNP\n')
f1.close()
for vcf_ref_file in vcf_ref[i*20:(i+1)*20]:
    vcf_ref_file_name = os.path.split(vcf_ref_file)[-1].split('.all.')[0]
    Depth_set = dict()
    summary_file_output = []
    try:
        # WGS
        vcf_1 = glob.glob(output_dir + vcf_ref_file_name + '*.all.fq.raw.vcf')[0]
        vcf_2 = glob.glob(output_dir + vcf_ref_file_name + '*.all.fq.flt.snp.vcf.filtered.snp.txt')[0]
        vcf_3 = glob.glob(output_dir + vcf_ref_file_name + '*.all.fq.flt.snp.vcf')[0]
        vcf_inputref = load_vcf(vcf_ref_file,True)
        depthcheck(vcf_1, vcf_2, vcf_3)
        summary_file_output = []
        for Depth in Depth_set:
            summary_file_output.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (vcf_ref_file_name,Depth,
                                                                 Depth_set[Depth][0], Depth_set[Depth][1],
                                                                 Depth_set[Depth][2], Depth_set[Depth][3]))

        f1 = open(summary_file, 'a')
        f1.write(''.join(summary_file_output))
        f1.close()
        print(Depth_set.get(3, 'None'), Depth_set.get(10, 'None'))
    except IndexError:
        pass

################################################### END ########################################################
