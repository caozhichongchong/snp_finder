# start
# merge allele sum of metagenomes for clonal populations
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
from statistics import stdev
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
# optional output setup
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='snp_output/',
                      metavar='snp_output/')

################################################## Definition ########################################################
args = parser.parse_args()
# setup path all clonal population
output_dir = args.o + '/MG/snpsum/'
sumoutput_dir = args.o + '/MG/summary/'
try:
    os.mkdir(sumoutput_dir)
except IOError:
    pass

try:
    os.mkdir(sumoutput_dir + '/withinHS')
except IOError:
    pass

try:
    os.mkdir(sumoutput_dir +'/acrossHS')
except IOError:
    pass

try:
    os.mkdir(sumoutput_dir +'/acrossdonor')
except IOError:
    pass

try:
    os.mkdir(sumoutput_dir +'/other')
except IOError:
    pass
# setup cutoff
#min_sum_allele_freq = 10
min_separate_allele_freq = 5
max_major_allele_freq = 0.95
#max_sum_allele_freq = 1000
################################################### Function ########################################################
__metaclass__ = type

class SNP_lineage:
    # create a class to store SNP_lineage
    'a class to store SNP_lineage'
    def init(self, lineage):
        self.lineage = lineage
        self.position = dict()
        self.metasample = []
    def addSNP(self,CHR_POS,Major_ALT,Minor_ALT,withinHS,allHS,acrossdonor,genename):
        self.position.setdefault(CHR_POS,[Major_ALT,Minor_ALT,withinHS,allHS,acrossdonor,genename,[],[]])
    def addmetasample(self,samplename):
        if samplename not in self.metasample:
            self.metasample.append(samplename)
    def addmeta(self,CHR_POS,samplename,Major_ALT_freq,Minor_ALT_freq):
        loci =self.metasample.index(samplename)
        if loci > len(self.position[CHR_POS][-2]) - 1:
            self.position[CHR_POS][-2].append(Major_ALT_freq)
            self.position[CHR_POS][-1].append(Minor_ALT_freq)
        else:
            self.position[CHR_POS][-2][loci] = Major_ALT_freq
            self.position[CHR_POS][-1][loci] = Minor_ALT_freq
            print('conflict samples for %s'%(samplename))
            print(len(self.position[CHR_POS][-2]),len(self.metasample))

# set up functions
def init_lineage(sumfiles,temp_SNP_lineage):
    for lines in open(sumfiles,'r'):
        if not lines.startswith("CHR\tPOS"):
            lines_set = lines.split('\t')
            # add SNP info
            CHR, POS, Major_ALT, Minor_ALT = lines_set[0:4]
            CHR_POS = '%s\t%s' % (CHR, POS)
            genename = lines_set[9]
            withinHS = lines_set[7]
            allHS = lines_set[8]
            acrossdonor = lines_set[9]
            temp_SNP_lineage.addSNP(CHR_POS,Major_ALT,Minor_ALT,withinHS,allHS,acrossdonor,genename)
            # add allele info

def remove_small_freq(freq):
    if freq >= min_separate_allele_freq:
        return freq
    else:
        return 0

def remove_small_freq2(freq,freq_sum):
    if freq == freq_sum:
        return freq_sum
    else:
        return 0

def compute_freq(freq1,freq2,freq3):
    freq1 = int(freq1)
    freq2 = int(freq2)
    freq3 = int(freq3)
    freq1 = remove_small_freq(freq1)
    freq2 = remove_small_freq(freq2)
    freq3 = remove_small_freq(freq3)
    freq_sum = freq1 + freq2 + freq3
    if freq_sum > 0:
        MAF = max(freq1, freq2, freq3)/ freq_sum
        if MAF == 1 or MAF <= max_major_allele_freq:
            # fixed alleles, or max major_ALT <= 0.95
            return [freq1/freq_sum,freq2/freq_sum,freq3/freq_sum,compare_freq(freq1, freq2, freq3),freq_sum]
        else:
            # MLF > 0.95, mostly likely to be fixed
            freq_sum = max(freq1, freq2, freq3)
            freq1 = remove_small_freq2(freq1, freq_sum)
            freq2 = remove_small_freq2(freq2, freq_sum)
            freq3 = remove_small_freq2(freq3, freq_sum)
            return [freq1 / freq_sum, freq2 / freq_sum, freq3 / freq_sum, compare_freq(freq1, freq2, freq3), freq_sum]
    else:
        return [0,0,0,'None\tNone',freq_sum]


def compare_freq(freq1,freq2,freq3):
    freq_max = max(freq1,freq2,freq3)
    freq_sum = freq1 + freq2 + freq3
    if freq_max > 0:
        if freq1 == freq_max:
            Tag1 = 'Major'
        elif freq2 == freq_max:
            Tag1 = 'Minor'
        else:
            Tag1 = 'Other'
        if freq_max == freq_sum:
            Tag2 = 'fixed'
        else:
            Tag2 = 'Poly'
        return  '%s\t%s'%(Tag1,Tag2)
    else:
        return 'None\tNone'


def output_sum_sub(newlines,output_list,Major_freq,Minor_freq,other_freq):
    if Major_freq > 0:
        output_list.append('%s\t%s\t%s\t\n'%(newlines,Major_freq,'Major'))
    if Minor_freq > 0:
        output_list.append('%s\t%s\t%s\t\n' % (newlines, Minor_freq, 'Minor'))
    if other_freq  > 0:
        output_list.append('%s\t%s\t%s\t\n' % (newlines, other_freq, 'Other'))
    return output_list

def output_sum(sumfiles,output_list_within,output_list_across,output_list_other,output_list_acrossdonor):
    samplename = os.path.split(sumfiles)[-1].split('.IN.')[0]
    for lines in open(sumfiles,'r'):
        if not lines.startswith("CHR\tPOS"):
            lines_set = lines.split('\t')
            Major_freq,Minor_freq,other_freq = lines_set[4:7]
            Major_freq, Minor_freq, other_freq, Tag, sum_freq = compute_freq(Major_freq,Minor_freq,other_freq)
            within_HS, HS_all,acrossdonor = lines_set[7:10]
            newlines = '\t'.join(lines_set[0:4]) + '\t%s\t%s\t%s\t%s'%(lines_set[10],samplename,sum_freq,Tag)
            if within_HS == 'True':
                output_list_within = output_sum_sub(newlines, output_list_within, Major_freq, Minor_freq, other_freq)
            elif HS_all == 'True':
                output_list_across = output_sum_sub(newlines, output_list_across, Major_freq, Minor_freq, other_freq)
            elif acrossdonor == 'True':
                output_list_acrossdonor= output_sum_sub(newlines, output_list_acrossdonor, Major_freq, Minor_freq, other_freq)
            output_list_other = output_sum_sub(newlines, output_list_other, Major_freq, Minor_freq, other_freq)


################################################### Main ########################################################
allsum = glob.glob('%s/*/*.IN.*.snp.sum'%(output_dir))
# split sum into sum per lineage
alllineage = dict()
for sumfiles in allsum:
    lineage = os.path.split(sumfiles)[-1].split('.IN.')[1].split('.snp.sum')[0]
    alllineage.setdefault(lineage,set())
    alllineage[lineage].add(sumfiles)

# process allele freq sum for each lineage
for lineage in alllineage:
    output_list_other = []
    output_list_within = []
    output_list_across = []
    output_list_acrossdonor = []
    allsumfiles = alllineage[lineage]
    for sumfiles in allsumfiles:
        output_sum(sumfiles, output_list_within,output_list_across,output_list_other,output_list_acrossdonor)
    if len(output_list_within) > 0:
        f1 = open(os.path.join(sumoutput_dir,'withinHS/%s.withinHS.snp.sum'%(lineage)),'w')
        f1.write('CHR\tPOS\tMajor_ALT\tMinor_ALT\tgenename\tsamplename\tSum_freq\tMax_ALT\tFixed\tFreq\tFreq_tag\t\n')
        f1.write(''.join(output_list_within))
        f1.close()
    if len(output_list_across) > 0:
        f1 = open(os.path.join(sumoutput_dir,'acrossHS/%s.acrossHS.snp.sum'%(lineage)),'w')
        f1.write('CHR\tPOS\tMajor_ALT\tMinor_ALT\tgenename\tsamplename\tSum_freq\tMax_ALT\tFixed\tFreq\tFreq_tag\t\n')
        f1.write(''.join(output_list_across))
        f1.close()
    if len(output_list_other) > 0:
        f1 = open(os.path.join(sumoutput_dir,'other/%s.other.snp.sum'%(lineage)),'w')
        f1.write('CHR\tPOS\tMajor_ALT\tMinor_ALT\tgenename\tsamplename\tSum_freq\tMax_ALT\tFixed\tFreq\tFreq_tag\t\n')
        f1.write(''.join(output_list_other))
        f1.close()
    if len(output_list_acrossdonor) > 0:
        f1 = open(os.path.join(sumoutput_dir,'acrossdonor/%s.other.snp.sum'%(lineage)),'w')
        f1.write('CHR\tPOS\tMajor_ALT\tMinor_ALT\tgenename\tsamplename\tSum_freq\tMax_ALT\tFixed\tFreq\tFreq_tag\t\n')
        f1.write(''.join(output_list_acrossdonor))
        f1.close()

# merge MG results
outputdir = '%s/withinHS'%(sumoutput_dir)
alloutput = []
os.system('rm %s'%('%s/all.withinHS.snp.sum'%(outputdir)))
for files in glob.glob('%s/*withinHS.snp.sum'%(outputdir)):
    donor = os.path.split(files)[-1].split('.withinHS.snp.sum')[0]
    for lines in open(files, 'r'):
        alloutput.append('%s\t' % (donor) + lines)

f1 = open('%s/all.withinHS.snp.sum'%(outputdir),'w')
f1.write(''.join(alloutput))
f1.close()

outputdir = '%s/other'%(sumoutput_dir)
alloutput = []
os.system('rm %s'%('%s/all.other.snp.sum'%(outputdir)))
for files in glob.glob('%s/*other.snp.sum'%(outputdir)):
    donor = os.path.split(files)[-1].split('.other.snp.sum')[0]
    for lines in open(files, 'r'):
        alloutput.append('%s\t' % (donor) + lines)

f1 = open('%s/all.other.snp.sum'%(outputdir),'w')
f1.write(''.join(alloutput))
f1.close()

outputdir = '%s/acrossHS'%(sumoutput_dir)
alloutput = []
os.system('rm %s'%('%s/all.acrossHS.snp.sum'%(outputdir)))
for files in glob.glob('%s/*acrossHS.snp.sum'%(outputdir)):
    donor = os.path.split(files)[-1].split('.acrossHS.snp.sum')[0]
    for lines in open(files, 'r'):
        alloutput.append('%s\t' % (donor) + lines)

f1 = open('%s/all.acrossHS.snp.sum'%(outputdir),'w')
f1.write(''.join(alloutput))
f1.close()

os.system('mv %s/*/all*.sum %s/'%(sumoutput_dir,sumoutput_dir))
