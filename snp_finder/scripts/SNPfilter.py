import glob
import os
from datetime import datetime
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path to all vcf files",
                      type=str, default='.',
                      metavar='input/')
# optional input genome
optional.add_argument("-cluster",
                      help="a cluster to run, default is all clusters",
                      type=str, default='',
                      metavar='cluster1')
################################################## Definition ########################################################
args = parser.parse_args()
################################################### Set up ########################################################
# Set up cutoff
min_qual_for_call = 40 #Remove sample*candidate that has lower than this quality
min_maf_for_call = .9 #Remove sample*candidate
min_cov = 3 # at least 3 reads mapped to POS
if 'human' in args.i:
    min_cov = 10
bowtievcfsuffix = 'flt.snp.vcf'#'.flt.snp.vcf'
mappervcfsuffix = 'mapper1.vcf'
compare_to_baseline = False
Middle_only = False
if 'noindel' in args.i:
    compare_to_baseline = True
################################################### Function ########################################################
def filter_snp(depthall,depthsnp):
    MAF = depthsnp/depthall
    if (MAF == 1 and depthsnp >= min_cov) or \
            (MAF >= min_maf_for_call and depthsnp >= min_cov * 2):
        return True
    return False

def load_baselinesnp(bowtievcf):
    baseline_chrpos = set()
    for lines in open(bowtievcf + '.final.txt', 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        baseline_chrpos.add('%s\t%s'%(lines_set[0],lines_set[1]))
    return baseline_chrpos

def load_bowtie(bowtievcf):
    try:
        f1 = open(bowtievcf + '.final.txt', 'r')
    except IOError:
        baseline_chrpos = set()
        if compare_to_baseline:
            # load baseline
            numsnp = os.path.split(bowtievcf)[-1].split('.SNP.fasta')[0].split('.')[-1]
            tool = os.path.split(bowtievcf)[-1].split('.SNP.fasta.')[-1].split('.flt.snp.vcf')[0]
            if numsnp != '0':
                baselinegenome = bowtievcf.replace('%s.SNP.fasta'%(numsnp),'0.SNP.fasta')
                print('%s start loading baseline %s' % (datetime.now(), baselinegenome))
                if baselinegenome not in baseline_set:
                    baseline_set.setdefault(baselinegenome + tool,load_baselinesnp(baselinegenome))
                baseline_chrpos = baseline_set[baselinegenome + tool]
                print('%s finished loading baseline %s' % (datetime.now(), baselinegenome))
        # load depth of qualified indel to support snps on POSs of that indel
        indel_file = bowtievcf.replace('.flt.snp.vcf','.flt.indel.vcf')
        indel_depth = dict()
        print('%s start loading indel depth %s' % (datetime.now(), indel_file))
        for lines in open(indel_file, 'r'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS, USELESS, REF, ALT, QUAL = lines_set[:6]
            if float(QUAL) >= min_qual_for_call:
                POS = int(POS)
                DPset = [float(x) for x in lines_set[9].split(':')[-1].replace('\n', '').split(',')]
                DP = sum(DPset)
                for POSsub in range(POS, POS + len(REF)):
                    indel_depth.setdefault('%s\t%s'%(CHR,POSsub),DP)
        print('%s finished loading indel depth %s' % (datetime.now(), indel_file))
        # grep all snps
        vcfoutput = []
        snpoutput = ['CHR\tPOS\tREF\tALT\tDPall\tDPsnp\n']
        # load each line of snps
        snpline = 0
        for lines in open(bowtievcf, 'r'):
            snpline += 1
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS, USELESS, REF, ALT, QUAL = lines_set[:6]
            ALT = ALT.split(',')
            if REF != 'N' and float(QUAL) >= min_qual_for_call:
                if baseline_chrpos== set() or '%s\t%s'%(CHR,POS) not in baseline_chrpos:
                    withsnp = False
                    DPset = [float(x) for x in lines_set[9].split(':')[-1].replace('\n', '').split(',')]
                    DP = sum(DPset) + indel_depth.get('%s\t%s'%(CHR,POS),0)
                    for i in range(0, len(ALT)):
                        # each potential ALT
                        depthsnp = DPset[i+1]# skip REF
                        if filter_snp(DP, depthsnp):
                            # a qualified snp
                            withsnp = True
                            snpoutput.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (CHR, POS, REF, ALT[i], DP, depthsnp))
                    if withsnp:
                        vcfoutput.append(lines)
            if snpline % 1000000 == 0:
                print(datetime.now(), 'processed %s lines' % (snpline))
        # output qualified snps
        f1 = open(bowtievcf + '.final.vcf','w')
        f1.write(''.join(vcfoutput))
        f1.close()
        f1 = open(bowtievcf + '.final.txt', 'w')
        f1.write(''.join(snpoutput))
        f1.close()

def load_mapper(mappervcf):
    try:
        f1 = open(mappervcf + '.final2.txt', 'r')
    except IOError:
        # grep all snps
        os.system('grep \';\' %s > %s.snp' % (mappervcf, mappervcf))
        vcfoutput = ['']
        snpoutput = ['CHR\tPOS\tREF\tALT\tDPall\tDPsnp\n']
        # load each line of snps
        snpline = 0
        baseline_chrpos = set()
        if compare_to_baseline:
            # load baseline
            numsnp = os.path.split(mappervcf)[-1].split('.SNP.fasta')[0].split('.')[-1]
            if numsnp != '0':
                baselinegenome = mappervcf.replace('%s.SNP.fasta' % (numsnp), '0.SNP.fasta')
                tool = os.path.split(mappervcf)[-1].split('.SNP.fasta.')[-1].split('.flt.snp.vcf')[0]
                print('%s start loading baseline %s' % (datetime.now(), baselinegenome))
                if baselinegenome not in baseline_set:
                    baseline_set.setdefault(baselinegenome + tool, load_baselinesnp(baselinegenome))
                baseline_chrpos = baseline_set[baselinegenome + tool]
                print('%s finished loading baseline %s' % (datetime.now(), baselinegenome))
        for lines in open(mappervcf + '.snp', 'r'):
            snpline += 1
            if not lines.startswith("CHR"):
                lines_set = lines.split('\n')[0].split('\t')
                CHR,POS,REF,ALT,DP,middleDP,endDP = lines_set[:7]
                ALT = ALT.split(',')
                if REF not in ['N','-']:
                    if baseline_chrpos == set() or '%s\t%s' % (CHR, POS) not in baseline_chrpos:
                        withsnp = False
                        DP = float(DP)
                        middleDP = [x.split(',') for x in middleDP.split(';')]
                        endDP = [x.split(',') for x in endDP.split(';')]
                        if Middle_only:
                            DP = sum([float(x[0]) + float(x[1]) for x in middleDP])
                        if DP > 0:
                            for i in range(0,len(ALT)):
                                # each potential ALT
                                if ALT[i] != '-':
                                    # not indel
                                    # using both middle and end depth
                                    if Middle_only:
                                        depthsnp = sum([float(x) for x in middleDP[i + 1]]) # skip REF
                                    else:
                                        depthsnp = sum([float(x) for x in middleDP[i + 1]]) + sum([float(x) for x in endDP[i + 1]]) # skip REF
                                    if filter_snp(DP, depthsnp):
                                        # a qualified snp
                                        withsnp = True
                                        snpoutput.append('%s\t%s\t%s\t%s\t%s\t%s\n'%(CHR,POS,REF,ALT[i],DP,depthsnp))
                            if withsnp:
                                vcfoutput.append(lines)
            if snpline % 1000000 == 0:
                print(datetime.now(), 'processed %s lines' % (snpline))
        # output qualified snps
        f1 = open(mappervcf + '.final.vcf','w')
        f1.write(''.join(vcfoutput))
        f1.close()
        f1 = open(mappervcf + '.final.txt', 'w')
        f1.write(''.join(snpoutput))
        f1.close()

# find all snp vcfs
allvcf_bowtie = glob.glob('%s/SNP_model*/merge/%s*%s'%(args.i,args.cluster,bowtievcfsuffix))
allvcf_bowtie.sort()
allvcf_mapper = glob.glob('%s/SNP_model*/merge/%s*%s'%(args.i,args.cluster,mappervcfsuffix))
allvcf_mapper.sort()
# process bowtie vcfs
baseline_set = dict()
for bowtievcf in allvcf_bowtie:
    print('%s start processing %s' % (datetime.now(), bowtievcf))
    load_bowtie(bowtievcf)
    print('%s finished processing %s' % (datetime.now(), bowtievcf))
# process mapper vcfs
for mappervcf in allvcf_mapper:
    print('%s start processing %s' % (datetime.now(), mappervcf))
    load_mapper(mappervcf)
    print('%s finished processing %s' % (datetime.now(), mappervcf))