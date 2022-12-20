import glob
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from pandas import Series, DataFrame,date_range
import warnings
from scipy.stats import poisson
from datetime import datetime
warnings.filterwarnings('ignore')

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-vcf",
                      help="path of folders of all vcfs",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/moredetails/moredetails/',
                      metavar='input/')
required.add_argument("-MW",
                      help="path of folders of all vcfs",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/MV_cov',
                      metavar='input/')

required.add_argument("-o",
                      help="path of folders of all vcfs",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/mergevcf/',
                      metavar='input/')
required.add_argument("-ref",
                      help="path of folders of all vcfs",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/',
                      metavar='input/')

################################################## Definition ########################################################
args = parser.parse_args()
min_maf = .9 # min major allele frequency
min_depth = 6 # for each genome and for each SNP
min_quality_good_alignment = 100 # quality for good alignment
max_bad_alignment_samples = .4 # % of samples with bad alignment
min_quality_alignment = 50 # quality for alignment
min_sample_alignment = 0.95 # at least X% samples have alignment >= min_quality_alignment
min_CI_coverage_for_genome = 0.05 # at least 5% CI for minimum coverage
depth_range = [0.01,0.99]  #quantile depth for POSs to call SNPs
depth_distance = 2 # remove up and down X*1000 bp regions with low depth
depth_distance_diff = 10 # neighbour average depth has X fold difference

exclusion_list = ['D25_BA_27','D25_BA_26','D25_BA_24','D25_BA_18','D25_BA_17','D25_BA_16','D25_BA_03',
                  'H26_BA_07','H26_BA_13','H26_BA_14','H26_BA_06','H26_BA_15',
                  'H26_BA_21','H26_BA_34','H26_BA_37','H26_BA_39','H26_BA_27',
                  'H26_BA_20','H26_BA_26','H26_BA_04','H26_BA_02',
                  'H29_BA_30','H29_BA_21','H29_BA_19','H29_BA_16','H29_BA_06',
                  'P76_PB_10','P76_PB_30','P76_PB_27','P76_PB_23',
                  'P76_PB_03','P76_PB_31','P76_PB_11','P76_PB_08','P76_PB_02',
                  'H42_PB_19','H42_PB_24','H42_PB_17','H42_PB_15','H42_PB_13','H42_PB_09',
                  'H42_PB_29','H42_PB_26','H42_PB_23','H42_PB_22','H42_PB_21','H42_PB_20','H42_PB_02',
                  'H09_BL_24','H09_BL_19','H09_BL_17','H09_BL_01',
                  'P45_PB_27','P45_PB_18','P45_PB_10','P45_PB_09','P45_PB_25','P45_PB_24','P45_PB_11','P45_PB_02',
                  'H23_BA_11','H23_BA_02','H23_BA_01',
                  'D503_PB_04','D503_PB_29','D503_PB_20',
                  'D503_PB_17','D503_PB_28','D503_PB_22',
                  'D503_PB_18','D503_PB_16','D503_PB_10','D503_PB_06','D503_PB_03',
                  'P50_PB_12','P50_PB_17','P50_PB_10','P50_PB_06',
                  'P50_PB_21','P50_PB_23','P50_PB_19','P50_PB_18','P50_PB_09','P50_PB_02',
                  'D503_BA_22','D503_BA_12','D503_BA_19','D503_BA_08',
                  'H30_BA_22','H30_BA_19','H30_BA_16','H30_BA_13',
                  'H30_BA_12','H30_BA_11','H30_BA_09','H30_BA_07',
                  'H30_BA_05','D92_BL_28','D92_BL_26',
                  'HX_PB_06','HX_PB_19','HX_PB_08','HX_PB_13','HX_PB_22','HX_PB_21',
                  'HX_PB_20','HX_PB_29','HX_PB_25','HX_PB_18','HX_PB_03','HX_PB_02',
                  'P44_PB_13','P44_PB_10','P44_PB_17','P44_PB_11','P44_PB_09','P44_PB_08','P44_PB_07','P44_PB_02',
                  'H13_BA_29','H13_BA_28','H13_BA_32','H13_BA_31','H13_BA_30','H13_BA_27','H13_BA_26','H13_BA_25','H13_BA_24',
                  'Bfrag_32','S1_T7_07','S1_T8_04','S1_T5_03','S1_T7_06','S1_T6_13','S1_T8_14','S1_T7_04','S1_T8_13','S1_T7_12',
                  'S1_T6_14','S1_T5_06','S1_T6_02','S1_T6_11','S1_T8_10','S1_T4_05','S1_T8_02','S1_T1_01','am_RuLa_g0010',
                  'D466_BA_19','D466_BA_14','D466_BA_07','D466_BA_17','D466_BA_20','D466_BA_12','D466_BA_02','H32_BA_27',
                  'H32_BA_26','H32_BA_23','H32_BA_21','H32_BA_17','H32_BA_02','D544_BA_14','D05_BA_14','D05_BA_18','D05_BA_12','D05_BA_11',
                  'D05_BA_22','D05_BA_05','D05_BA_08','D05_BA_04','D549_BL_25','D549_BL_26','P70_BL_12','P70_BL_20',
                  'P70_BL_17','P70_BL_09','P70_BL_13','P70_BL_05','P02_BA_02','P02_BA_11','P02_BA_08','P02_BA_07','P02_BA_01',
                   'P70_BL_10','P70_BL_08','P70_BL_02','H17_BA_15','H17_BA_24','H17_BA_31','H17_BA_22','H17_BA_11','H17_BA_08',
                  'D465_BA_21', 'D465_BA_18', 'D465_BA_12', 'D465_BA_04', 'D465_BA_03', 'D465_BA_01',
                  'D465_BA_19','D14_BA_18','D14_BA_22','D14_BA_15','H37_BA_26','H37_BA_14','H37_BA_06','H37_BA_03',
                  'am_BL_g0059','am_BL_g0056','am_BL_g0052','am_BL_g0032','am_BL_g0041','am_BL_g0036','am_BL_g0046','am_BL_g0019',
                  'am_BL_g0040','am_BL_g0016','am_BL_g0008','am_BL_g0039','am_BL_g0035','am_BL_g0034','am_BL_g0024','am_BL_g0022',
                  'am_BL_g0014','am_BL_g0004','am_BL_g0044','am_BL_g0031','am_BL_g0015','am_BL_g0009','am_BL_g0003','am_BL_g0002',
                  'HX_BA_14','HX_BA_09','HX_BA_01','H02_BA_05','H02_BA_18','D442_BA_08',
                  ]

exclusion_list2 = ['HX_PB_16','af_TuSa_g0054','af_TuSa_g0031','af_TuSa_g0073','af_TuSa_g0058','af_TuSa_g0053','af_TuSa_g0023','af_TuSa_g0001',
                   'bj_EsCo_g0023','P48_PB_08','P48_PB_07','P48_PB_04','P46_PB_07','P46_PB_13','P79_PB_08',
                   'af_BlWe_g0010','af_BlWe_g0007','af_BlWe_g0001','bf_BaVu_g0001','H01_BL_30','H01_BL_23','H02_BA_06','H02_BA_03',
                   'ao_BL_g0007','ao_BL_g0001','ao_BL_g0006','H26_BL_03','D442_BL_06','D442_BL_01','D442_BL_19','D442_BL_03',
                   'D549_BL_19','D549_BL_22','D549_BL_08','D549_BL_21','D549_BL_10','D549_BL_07','D549_BL_20','D549_BL_16','D549_BL_09',
                   'D84_BL_03','D84_BL_01','P67_BL_25','P67_BL_12','P67_BL_11','P67_BL_13','H25_BL_05','H25_BL_02','H25_BL_08','H25_BL_01',
                   'H11_BL_26','H11_BL_27','H11_BL_18','H11_BL_17','H11_BL_16','H11_BL_15','H11_BL_14','H11_BL_12','H11_BL_05','H11_BL_04',
                   'H11_BL_03','H11_BL_01','D92_BL_30','D92_BL_23','D92_BL_17','D92_BL_18','D92_BL_11',
                   'H17_BL_03','H17_BL_08','P70_BL_23','P70_BL_30','P70_BL_07',
                   'P70_BL_22','H02_BL_24','H02_BL_20','H02_BL_19','H02_BL_13','H02_BL_07','H02_BL_06','H02_BL_05','H02_BL_04','H02_BL_03',
                   'H02_BL_02','P76_BL_30','P76_BL_26','P76_BL_15','P76_BL_14','P76_BL_11','P76_BL_09','P76_BL_03','P76_BL_01','D114_BA_27',
                   'D114_BA_26','D114_BA_23','D114_BA_21','D114_BA_19','D114_BA_16','D114_BA_11','D114_BA_13','D114_BA_10','D114_BA_06',
                   'D114_BA_04','H17_BA_17','H17_BA_14','H17_BA_30','H17_BA_29','H17_BA_21','H17_BA_16','H17_BA_35','H17_BA_34','H17_BA_33',
                   'H17_BA_32','H17_BA_28','H17_BA_23','H17_BA_20','H17_BA_18','H17_BA_09','H28_BA_09','H28_BA_22','an_BA_g0006','P02_BA_19',
                   'P02_BA_09','P02_BA_22','P02_BA_16','P02_BA_14','P02_BA_12','P02_BA_04','D466_BA_22','D466_BA_21','D466_BA_18','D466_BA_11',
                   'ao_BA_g0046','ao_BA_g0033','ao_BA_g0029','ao_BA_g0031','ao_BA_g0024','ao_BA_g0008','ao_BA_g0006','D115_BA_18','D115_BA_17',
                   'D115_BA_07','D115_BA_04','D115_BA_02','D115_BA_01','bf_BA_g0014','bf_BA_g0006','D465_BA_05','D465_BA_02',
                   'H39_BA_16','H39_BA_05','H39_BA_20','H39_BA_19','H39_BA_22','H39_BA_15','H39_BA_26','H39_BA_04','D448_BA_04','D448_BA_12',
                   'D448_BA_11','D448_BA_10','D448_BA_02','H30_BA_08','H30_BA_24','H30_BA_20','H30_BA_18','H30_BA_15','H30_BA_04',
                   'H29_BA_26','H29_BA_05','H29_BA_28','H29_BA_20','H29_BA_23',
                   'H29_BA_13','H29_BA_12','H29_BA_10','H29_BA_02','H21_BA_05','H21_BA_26','H21_BA_29','H21_BA_08','H21_BA_27','H21_BA_25',
                   'H21_BA_04','H21_BA_22','H21_BA_03','D98_BA_11','D98_BA_07','D98_BA_06','D98_BA_19','D98_BA_13','D98_BA_08','D98_BA_17',
                   'D98_BA_12','D98_BA_05','D98_BA_02','H20_BA_21','H20_BA_29','H20_BA_27','H20_BA_24','H20_BA_23','H20_BA_20',
                   'H20_BA_17','H20_BA_14','H20_BA_11','H20_BA_08','H20_BA_03','H20_BA_02','H20_BA_01','H13_BA_15','H13_BA_14',
                   'D544_BA_20','D544_BA_22','D544_BA_23','D544_BA_18','D05_BA_15','D01_BA_18','D01_BA_13','D01_BA_08','D14_BA_13','D14_BA_02',
                   'D14_BA_16','D14_BA_29','D14_BA_26','D14_BA_25','D14_BA_24','D14_BA_23','D14_BA_21','D14_BA_17','D14_BA_09','D14_BA_08',
                   'D14_BA_01','H37_BA_24','H37_BA_21','H37_BA_16','H37_BA_11','H32_BA_21','H32_BA_17','am_EsCo_g0046',
                   'am_BL_g0017','am_BL_g0051','H09_BL_25','H09_BL_28','H09_BL_26','H09_BL_27','H09_BL_23',
                   'am_BL_g0050','am_BL_g0005','am_BL_g0053','am_BL_g0012','am_BL_g0011','am_BL_g0010','am_BL_g0057','am_BL_g0025','am_BL_g0021',
                   'am_BL_g0049','am_BL_g0023','am_BL_g0001','D82_BL_17','D82_BL_07''D442_BA_05','D442_BA_03',]

exclusion_list3 = ['H13_BA_17','H21_BA_12','H21_BA_21','H21_BA_02','D98_BA_04', 'D465_BA_16', 'D465_BA_15', 'D465_BA_13', 'D465_BA_09','bf_BA_g0028','P70_BL_26','H28_BA_07','H30_BA_02','H30_BA_21','H20_BA_15','D05_BA_17','D05_BA_09','D05_BA_06','D05_BA_03','H37_BA_22',
'H32_BA_29','H32_BA_28','H32_BA_19','H32_BA_15','H32_BA_12','H32_BA_04','H29_BA_15','H29_BA_25','P45_BA_29','H17_BA_19','H17_BA_27','D14_BA_10',
                   'af_PsSp_g0019','af_PsSp_g0012','af_PsSp_g0014','af_PsSp_g0007','af_PsSp_g0016','af_PsSp_g0002','aa_PaGo_g0011',
                   'aa_PaGo_g0008','aa_PaGo_g0003','aa_PaGo_g0002','an_PaDi_g0010','an_PaDi_g0009','an_PaDi_g0004','an_PaDi_g0003',
                   'an_PaDi_g0002','D25_PB_19','D25_PB_15','D14_PB_24','D14_PB_19','D14_PB_17','D14_PB_11','HX_PB_29',
                   'H40_PB_12','H40_PB_11','H40_PB_05','H40_PB_03','H30_PB_19','H30_PB_09','H15_PB_19','H15_PB_17','H15_PB_08',
                   'P46_PB_21','P42_PB_22','P42_PB_20','P42_PB_19','P42_PB_18','P42_PB_11','P42_PB_03','H32_PB_10','H32_PB_12',
                   'H32_PB_02','P67_PB_32','P67_PB_18','P67_PB_27','P67_PB_21','P67_PB_17','P67_PB_13','P67_PB_10','P67_PB_006',
                   'P67_PB_35','P67_PB_31','P67_PB_25','P67_PB_14','P67_PB_12','P67_PB_11','P67_PB_09','P67_PB_005','P67_PB_001','P02_PB_14',
                   'P02_PB_08','P02_PB_06','P02_PB_005','P02_PB_17','P02_PB_11','P02_PB_003','P02_PB_001','av_LaRu_g0021','bj_EsCo_g0025',
                   'bj_EsCo_g0020','bj_EsCo_g0008','av_BiBi_g0016','av_BiBi_g0011','av_BiBi_g0007','av_BiBi_g0002','bj_BaXy_g0007','bj_BaXy_g0009',
                   'bj_BaXy_g0006','ao_BL_g0001','ao_BL_g0010','D01_BL_01','P67_BL_13','P41_BL_01','D82_BL_16','D82_BL_11','D82_BL_14',
                   'D82_BL_13','D82_BL_12','D82_BL_10','D82_BL_08','D82_BL_01','H16_BL_16','H16_BL_14','H16_BL_10','D92_BL_14','D92_BL_16',
                   'D92_BL_05','D92_BL_04','D542_BL_25','D542_BL_17','D542_BL_16','D542_BL_14','D542_BL_12','D542_BL_07','D542_BL_06',
                   'H02_BL_15','H02_BL_11','H18_BL_06','am_BL_g0030','H09_BL_14','H09_BL_16','H09_BL_10','H09_BL_09','H09_BL_07','H09_BL_04','H09_BL_02','D549_BL_14',
                   'P02_BA_15', 'P02_BA_13','P02_BA_05']
exclusion_list4 = ['H17_BA_13','H12_BA_13','D92_BL_15','H29_BA_09','H32_BA_09','H32_BA_30','H32_BA_14','H32_BA_11','H32_BA_20','H29_BA_24','H29_BA_17','am_BL_g0020','am_BL_g0047','am_BL_g0037',
                  'H32_BA_25','H32_BA_22','H32_BA_08','H30_BA_01','P70_BL_28','P70_BL_16','P70_BL_04','P70_BL_14','P70_BL_11','P70_BL_03','H37_BA_17','H37_BA_25','H37_BA_04']
def load_ref(ref_genome):
    Mapping_loci = dict()
    Reverse = []
    Ref_seq = dict()
    print('process reference database %s' % (ref_genome))
    for record in SeqIO.parse(ref_genome, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0:  # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
        Ref_seq.setdefault(record_id, record_seq)
    return [Ref_seq,Mapping_loci,Reverse]

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

def allele_set_convert(REF,ALT,AD,QUA):
    ALTset = [REF]
    if ALT != '.':
        ALTset += ALT.split(',')
    AD = AD.split(':')
    ADF = AD[-3].split(',')
    ADR = AD[-2].split(',')
    total_depth = [int(x) for x in ADF] + [int(x) for x in ADR]
    total_depth = sum(total_depth)
    allele_set = [0]*8
    QUA = float(QUA)
    main_allele = '-'
    i = 0
    for allele in ALTset:
            depth_ALT = int(ADF[i]) + int(ADR[i])
            allele_set[Allels[allele] * 2] = ADF[i]
            allele_set[Allels[allele] * 2 + 1] = ADR[i]
            if QUA >= min_quality_good_alignment and depth_ALT >= min_depth and depth_ALT/total_depth >= min_maf:
                # good depth and good major allele frequency
                main_allele = allele
            i += 1
    allele_set = [QUA] + allele_set
    allele_set = ','.join([str(x) for x in allele_set])
    return [main_allele,allele_set]

def merge_vcf(vcfset,low_quality_genomes):
    VCF_result = dict()
    sampleset = []
    for vcf in vcfset:
        samplename = os.path.basename(vcf).split('_1.fastq')[0].split('.')[1]
        if samplename not in low_quality_genomes:
            # remove low quality genomes
            sampleset.append(samplename)
    sample_length = len(sampleset)
    sampleID = 0
    for vcf in vcfset:
        samplename = os.path.basename(vcf).split('_1.fastq')[0].split('.')[1]
        if samplename not in low_quality_genomes:
            # remove low quality genomes
            for lines in open(vcf,'r'):
                if not lines.startswith('#') and 'INDE' not in lines:
                    lines_set = lines.split('\n')[0].split('\t')
                    contig, POS, nothing, REF, ALT, QUA = lines_set[:6]
                    contigPOS = '%s\t%s\t%s'%(contig,POS,REF)
                    main_allele, allele_set = allele_set_convert(REF, ALT, lines_set[-1],QUA)
                    VCF_result.setdefault(contigPOS,[['-']*sample_length,['0,0,0,0,0,0,0,0,0']*sample_length,0]) # main alleles, allele set, num of samples with alignment
                    VCF_result[contigPOS][0][sampleID] = main_allele
                    VCF_result[contigPOS][1][sampleID] = allele_set
                    if float(QUA) >= min_quality_alignment:
                        VCF_result[contigPOS][2] += 1
            sampleID += 1
    return [sampleset,VCF_result]

def contig_to_gene(CHR, POS,Mapping_loci,Reverse,Ref_seq):
    all_genes = Mapping_loci.get(CHR,[])
    Reverse_chr = 0
    for a_gene in all_genes:
        POS1, POS2, GENE = a_gene
        if POS >= POS1 and POS <= POS2:
            Ref_seq_chr = Ref_seq.get(GENE, 'None')
            Gene_length = len(Ref_seq_chr)
            if GENE in Reverse:  # reversed
                POS_gene = Gene_length-(int(POS-POS1))
                Reverse_chr = 1
            else:
                POS_gene = int(POS-POS1)+1
            codon_start = POS_gene - 1 - int((POS_gene - 1) % 3)
            return [GENE,POS_gene,codon_start,Ref_seq_chr,Reverse_chr]
    return []

def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def dnORds(amino1, amino2):
    if amino1 == amino2:
        return 'S'
    else:
        return 'N'

def causeSNP(seq,position,ALT,Reverse_chr,checkREF = True):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    if checkREF:
        if seq[position - 1] != ALT:
            print('WRONG REF %s %s %s'%(position, ALT, seq[position - 1]))
    else:
        seq = list(seq)
        seq[position - 1]=ALT
    return ''.join(seq)

def vcf_output(sampleset,VCF_result,Ref_seq, Mapping_loci, Reverse,bad_CHRPOS_list):
    alloutput = ['CHR\tPOS\tMajor_allele\tMinor_allele\tAlleles\tGenome_set\tGene\tGene_POS\tN_or_S\tAA_change\t%s\n'%(
        '\t'.join(sampleset)
    )]
    sample_length = len(sampleset)
    for contigPOS in VCF_result:
        contig, POS, REF = contigPOS.split('\t')
        bad_list = bad_CHRPOS_list.get(contig, [])
        if int(POS) not in bad_list:
            # good quality POS
            alleles, allallele_set, num_sample_with_alignment = VCF_result[contigPOS]
            ALT_set = dict()
            for allele in alleles:
                ALT_set.setdefault(allele,alleles.count(allele))
            if '-' in ALT_set:
                del ALT_set['-']
            if num_sample_with_alignment >= sample_length * min_sample_alignment and len(ALT_set) > 1 and alleles.count('-') < max_bad_alignment_samples*sample_length:
                # no samples have no alignment, a SNP, low quality alignment samples < max_bad_alignment_samples
                sample_lengthnew = (sample_length - alleles.count('-'))
                # find minor allele
                minor_allele_set = []
                Major = ''
                Minor = ''
                for allele in ALT_set:
                    if ALT_set[allele] >= 0.5 * sample_lengthnew and Major == '':
                        # major allele, set major allele as the first one passing 50% if empty
                        Major = allele
                    if ALT_set[allele] <= 0.5 * sample_lengthnew and allele != Major:
                        Minor = allele
                        # minor allele
                        for i in range(0, sample_length):
                            if alleles[i] == allele:
                                minor_allele_set.append(sampleset[i])
                    if Major != '' and Minor != '':
                        break
                # SNPs + good quality samples >= 1-max_bad_alignment_samples
                contig, POS, REF = contigPOS.split('\t')
                if REF in ALT_set:
                    del ALT_set[REF]
                POS = int(POS)
                # SNPs on a gene
                temp_snp_line_NS = ['None', 'None', 'None']
                temp_snp_line_AA = ''
                gene_info = contig_to_gene(contig,POS,Mapping_loci,Reverse,Ref_seq)
                if gene_info != []:
                    Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
                    if Ref_seq_chr != 'None':
                        #  observed NS
                        temp_snp_line_NS = [Chr_gene, str(POS_gene), '']
                        if codon_start <= POS_gene - 1:
                            Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr,True)
                            Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                            SNP_seq_chr = Ref_seq_chr
                            if len(Ref_seq_codon) == 3:
                                Ref_seq_aa = translate(Ref_seq_codon)[0]
                                temp_snp_line_AA += Ref_seq_aa
                                for ALT in ALT_set:
                                    if ALT != REF:
                                        SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr,False)
                                        SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                        SNP_seq_aa = translate(SNP_seq_codon)[0]
                                        temp_snp_line_AA += SNP_seq_aa
                                        temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                        temp_snp_line_NS[-1] += temp_NorS
                        if 'N' in temp_snp_line_NS[-1]:
                            temp_snp_line_NS[-1] = 'N'
                alloutput.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(contig,POS,Major,Minor,''.join(alleles),';'.join(minor_allele_set),
                                                       '\t'.join(temp_snp_line_NS),temp_snp_line_AA,'\t'.join(allallele_set)
                                                       ))
    return alloutput

def genome_screening(allSNPscov):
    exclusion_list = []
    total_genome_length = allSNPscov.shape[0]
    genome_coverage_info = dict()
    coverage_genomeset = []
    allgenomes = list(allSNPscov.columns)[3:]
    allgenomes.sort()
    for genome in allgenomes:
        average_genome_depth = np.mean(allSNPscov.loc[:,genome])
        total_genome_coverage = len([x for x in allSNPscov.loc[:,genome] if x >= min_depth])
        coverage_genomeset.append(total_genome_coverage/total_genome_length)
        genome_coverage_info.setdefault(genome,[average_genome_depth,total_genome_coverage/total_genome_length])
    coverage_genomeset_cutoff = poisson.ppf(min_CI_coverage_for_genome,total_genome_coverage*np.mean(coverage_genomeset))/total_genome_coverage
    for genome in genome_coverage_info:
        average_genome_depth,genome_coverage = genome_coverage_info[genome]
        if average_genome_depth < min_depth or genome_coverage < coverage_genomeset_cutoff:
            # low depth or low coverage
            # exclude genomes
            exclusion_list.append(genome)
    return exclusion_list

def depth_screening(allSNPscov,low_quality_genomes):
        allSNPscov = allSNPscov.loc[:,~allSNPscov.columns.isin(low_quality_genomes)]
        allSNPscov = allSNPscov[allSNPscov['POS'] > 100]  # not considering the depth of the first POS of each contig
        bad_CHRPOS_list = dict()
        if allSNPscov.shape[1] > 3 + 3:
            # at least 3 samples
            # avg depth of non zero
            allSNPscov['max_depth'] = 0
            allrows = allSNPscov.shape[0]
            allSNPscov.index = range(0, allrows)
            for i in allSNPscov.index:
                non_zero_depth = [x for x in allSNPscov.iloc[i, 3:] if x > 0]
                if len(non_zero_depth) > 1:
                    allSNPscov.loc[i, 'max_depth'] = max(non_zero_depth)
                    allSNPscov.loc[i, 'avg_depth'] = np.mean(non_zero_depth)
                else:
                    allSNPscov.loc[i, 'avg_depth'] = 0
                    # depth screening
            avg_depth_set = [x for x in allSNPscov['avg_depth'] if x > 0]
            max_depth_set = [x for x in allSNPscov['max_depth'] if x > 0]
            depth_range_lineage = np.quantile(avg_depth_set, depth_range)
            depth_range_lineage[0] = max(depth_range_lineage[0], min_depth)  # at least min_depth
            depth_range_max_lineage = np.quantile(max_depth_set, depth_range)
            depth_range_max_lineage[0] = max(depth_range_max_lineage[0], min_depth)  # at least min_depth
            invalid_CHRPOS = allSNPscov.loc[(allSNPscov['avg_depth'] <= depth_range_lineage[0]) | (
                    allSNPscov['avg_depth'] >= depth_range_lineage[1]) | (
                                                    allSNPscov['max_depth'] >= depth_range_max_lineage[1]),
                             :]  # for example, for plasmids
            # add regions with a drop of depth
            invalid_CHRPOS_list = list(invalid_CHRPOS.index)
            for i in allSNPscov.index[depth_distance:-depth_distance]:
                depth_neighbour = allSNPscov.loc[(i - depth_distance):(i + depth_distance), 'avg_depth']
                if max(depth_neighbour) / max(min(depth_neighbour), 0.1) >= depth_distance_diff:
                    # regions with a drop of depth
                    invalid_CHRPOS_list += [i]
            for i in list(set(invalid_CHRPOS_list)):
                CHR = allSNPscov.loc[i, 'CHR']
                POS = allSNPscov.loc[i, 'POS']
                bad_CHRPOS_list.setdefault(CHR,[])
                for i in range(POS - 1000*depth_distance, POS + 1000*depth_distance):
                    bad_CHRPOS_list[CHR].append(i)
        return bad_CHRPOS_list


################################################## Main ########################################################
# path to lineage vcf and ref
allvcf = glob.glob('%s/*.raw.vcf.filter'%(args.vcf))
allvcf.sort()
Lineage_genome = dict()
Lineage_ref = dict()
for vcf in allvcf:
    lineage = os.path.basename(vcf).split('.')[0]
    donor = os.path.basename(vcf).split('.')[1].split('_')[0]
    samplename = os.path.basename(vcf).split('_1.fastq')[0].split('.')[1]
    ref_genome = glob.glob('%s/%s/%s.all.spades2.fasta.noHM.fasta.fna'%(args.ref,lineage,lineage))
    if ref_genome != []:
        Lineage_ref.setdefault(lineage,ref_genome[0])
    if samplename in exclusion_list4:
        lineage_donor = '%s.donor.%s_lineage5' % (lineage, donor)
    elif samplename in exclusion_list3:
        lineage_donor = '%s.donor.%s_lineage4' % (lineage, donor)
    elif samplename in exclusion_list2:
        lineage_donor = '%s.donor.%s_lineage3' % (lineage, donor)
    elif samplename in exclusion_list:
        lineage_donor = '%s.donor.%s_lineage2' % (lineage, donor)
    else:
        lineage_donor = '%s.donor.%s'%(lineage,donor)
    Lineage_genome.setdefault(lineage_donor,[])
    Lineage_genome[lineage_donor].append(vcf)

# load ref
for lineage in Lineage_ref:
    ref_genome = Lineage_ref[lineage]
    Lineage_ref[lineage] = load_ref(ref_genome)

# load vcf
for lineage_donor in Lineage_genome:
    print(datetime.now(),'processing %s'%(lineage_donor))
    # load MW coverage
    try:
        try:
            f1 = open('%s/%s.vcf' % (args.o, lineage_donor), 'r')
        except FileNotFoundError:
            allSNPscov = pd.read_csv(os.path.join(args.MW,'%s.raw.vcf.filtered.cov.MW.txt') % (
                lineage_donor.replace('.donor', '.all.donor')).split('_lineage')[0]
                    ,sep='\t')
            # remove low quality genomes
            low_quality_genomes = genome_screening(allSNPscov)
            print(datetime.now(),'excluding genomes with low quality',low_quality_genomes)
            # remove low quality depth POS
            bad_CHRPOS_list = depth_screening(allSNPscov,low_quality_genomes)
            print(datetime.now(), 'excluding POSs with low depth')
            # load vcfs
            sampleset, VCF_result = merge_vcf(Lineage_genome[lineage_donor],low_quality_genomes)
            print(datetime.now(), 'merging vcfs')
            if len(sampleset) >= 3:
                # at least 3 genomes
                # compare alleles
                lineage = lineage_donor.split('.donor.')[0]
                print(lineage_donor,lineage)
                if lineage in Lineage_ref:
                    Ref_seq, Mapping_loci, Reverse = Lineage_ref[lineage]
                    alloutput = vcf_output(sampleset,VCF_result,Ref_seq, Mapping_loci, Reverse,bad_CHRPOS_list)
                    print(datetime.now(), 'outputing vcfs')
                    if len(alloutput) > 1:
                        f1 = open('%s/%s.vcf'%(args.o,lineage_donor),'w')
                        f1.write(''.join(alloutput))
                        f1.close()
                else:
                    print('no reference for %s'%(lineage))
    except FileNotFoundError:
        print('MW file not found for %s'%(lineage_donor))


################################################### END ########################################################
