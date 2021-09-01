cmds = 'python genome_sum.py -s . -o /scratch/users/anniz44/genomes/donor_species/SNP_curate/bwa'
cmds = 'python genome_curate.py -s . -i /scratch/users/anniz44/genomes/donor_species/SNP_curate/test_data -o /scratch/users/anniz44/genomes/donor_species/SNP_curate/bwa'
cmds = 'python SNP_model.py -s /scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/ -i /scratch/users/anniz44/genomes/donor_species/SNP_curate/test_data -o /scratch/users/anniz44/genomes/donor_species/SNP_curate/'
cmds = 'sh allsnpmodel.sh'
cmds = 'python SNPfilter_WGS_singlesample.py -i /scratch/users/anniz44/genomes/donor_species/SNP_curate/SNP_model/merge -vcf .flt.snp.vcf -s . -smp /scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/SNP_model/'
cmds = 'python SNP_model_compare.py -i /scratch/users/anniz44/genomes/donor_species/SNP_curate/SNP_model/merge/ -ref /scratch/users/anniz44/genomes/donor_species/SNP_curate/SNP_model/data/'
# indel only
cmds = 'python SNP_model_indelonly.py -s /scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/ -i /scratch/users/anniz44/genomes/donor_species/SNP_curate/test_data -o /scratch/users/anniz44/genomes/donor_species/SNP_curate/'
cmds = 'sh allsnpmodel.sh'
cmds = 'python SNPfilter_WGS_singlesample.py -i /scratch/users/anniz44/genomes/donor_species/SNP_curate/indel_model/merge -vcf .flt.snp.vcf -s . -smp /scratch/users/anniz44/scripts/1MG/donor_species/snp_curate/indel_model/'
#cmds = 'python SNP_model_compare.py -i /scratch/users/anniz44/genomes/donor_species/SNP_curate/indel_model/merge/ -ref /scratch/users/anniz44/genomes/donor_species/SNP_curate/indel_model/data/'

# crispr
cmds = 'python spacer_parsimony.py'
################################################### END ########################################################
# step 1 whole genomes fastq mapping to the assembly genome
import glob
import os
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/vcf_new'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')
genome_dir = ['/scratch/users/anniz44/genomes/donor_species/selected_species/af_Bacteroides_ovatus']
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_currate1'
fastq_name = '_1.fastq'

os.system('rm -rf %s'%(input_script_sub))

try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa/0')
except IOError:
    pass


try:
    os.mkdir(input_script)
except IOError:
    pass

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# generate mapping cmd
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    if len(donor) == 2: # BN10 donors
        # pick up one genome
        donor_species_fastq = []
        donor_species_genome_all = glob.glob(os.path.join(folder, '*.spades.fasta')) + \
                                   glob.glob(os.path.join(folder, '*.fasta'))
        total_genome = len(donor_species_genome_all)
        # pick up one genome
        #i = 0
        #while donor_species_fastq == [] and i < total_genome:
        #    donor_species_genome = donor_species_genome_all[i]
        #    if 'all.spades.fasta' not in donor_species_genome:
        #        genome_name = os.path.split(donor_species_genome)[-1].split('spades.fasta')[0].split('.fasta')[0]
        #        donor_species_fastq = glob.glob(os.path.join(folder, 'fastq/%s%s'%(genome_name,fastq_name)))
        #    i += 1
        # all genomes
        i = 0
        for donor_species_genome in donor_species_genome_all:
            if 'all.spades.fasta' not in donor_species_genome:
                genome_name = os.path.split(donor_species_genome)[-1].split('spades.fasta')[0].split('.fasta')[0]
                donor_species_fastq = glob.glob(os.path.join(folder, 'fastq/%s%s'%(genome_name,fastq_name)))
                if donor_species_fastq != []:
                    cmds = ''
                    database = donor_species_genome
                    try:
                        f1 = open(database + '.1.bt2', 'r')
                    except IOError:
                        os.system('bowtie2-build %s %s' % (database, database))
                    for files in donor_species_fastq:
                        donor_species_dir_file = os.path.split(files)[-1]
                        tempbamoutput = os.path.join(output_dir +'/bwa/0', donor_species_dir_file)
                        try:
                            f1 = open(tempbamoutput + '.sorted.bam')
                        except IOError:
                            files2 = files.replace(fastq_name, fastq_name.replace('1','2'))
                            cmds += 'bowtie2' + ' -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --threads %s --un-conc %s.unaligned -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                                min(40, 40),tempbamoutput, database, files, files2,'samtools', min(40, 40),
                                tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
                                tempbamoutput)
                            cmds += 'samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
                                tempbamoutput, tempbamoutput)
                            cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s.sorted.bam  | %s call  --threads %s -m > %s.raw.vcf\n' % (
                                'bcftools', min(40, 40), database,
                                tempbamoutput, 'bcftools', min(40, 40), tempbamoutput)
                            cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
                                'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
                            cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
                                'bcftools', tempbamoutput, tempbamoutput)
                            cmds += 'rm -rf %s.*.unaligned\n' % (tempbamoutput)
                            cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
                            f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species,i)), 'a')
                            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                            f1.close()
                            cmds = ''
                            i += 1
                else:
                    print('no reference for genome',donor_species_genome)

f1 = open(os.path.join(input_script, 'allnewvcf.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
# step 2 check SNPs
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
# set up path
input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/vcf_new'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_currate1'
fastq_name = '_1.fastq'
genome_split = '_g'

# set up cutoff
# good mapping
Major_alt_freq_cutoff = 0.8 # major alt freq in a sample
Sample_depth_cutoff = 5 # both forward and reverse reads cutoff in a sample
# good coverage
total_coverage_cutoff = 0.8 # at least X reads map to its original genome
genome_avg_coverage_cutoff = 10 # genome average coverage cutoff
# reasonable curation
Major_alt_freq_cutoff2 = 0.6 # major alt freq in a sample
Sample_depth_cutoff2 = 3 # both forward and reverse reads cutoff in a sample

Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def coverage_check(sh_file,Coverage):
    Coverage1 = []
    Coverage2 = []
    os.system('rm -rf %s/temp.sh %s/temp.err'%(input_script_split_sub,input_script_split_sub))
    os.system('grep \"bowtie2\" %s > %s/temp.sh' % (sh_file, input_script_split_sub))
    os.system('grep \"overall alignment rate\" %s.err > %s/temp.err' % (sh_file, input_script_split_sub))
    # load fastq name
    for lines in open('%s/temp.sh' % (input_script_split_sub), 'r'):
        try:
            filename = lines.split('>')[1].split('\n')[0].split(fastq_name + '.bam')[0]
            filename = os.path.split(filename)[-1]
            Coverage1.append(filename)
        except IndexError:
            pass
    # load reads alignment rate
    for lines in open('%s/temp.err' % (input_script_split_sub), 'r'):
        Coverage2.append(float(lines.split('%')[0]))
    i = 0
    for filename in Coverage1:
        # load average coverage of ref genome
        tempbamoutput = os.path.join(output_dir + '/bwa/0', filename + fastq_name)
        coverage_file = glob.glob(tempbamoutput + '.sorted.bam.avgcov')
        if coverage_file == []:
            cmds = 'samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
                tempbamoutput, tempbamoutput)
            os.system(cmds)
            coverage_file = tempbamoutput + '.sorted.bam.avgcov'
        else:
            coverage_file = coverage_file[0]
        length_total = 0
        coverage_total = 0
        donor_species = filename.split(genome_split)[0]
        Donor_species.setdefault(donor_species, [[], 0])
        try:
            for lines in open(coverage_file, 'r'):
                length_contig = float(lines.split('\t')[1])
                length_total += length_contig
                coverage_contig = float(lines.split('\t')[2])
                coverage_total += coverage_contig * length_contig
            coverage_num = Coverage2[i] / 100
            avg_coverage = coverage_total / length_total
            i += 1
        except IOError:
            coverage_num = total_coverage_cutoff - 0.1
            avg_coverage = genome_avg_coverage_cutoff - 1
        temp_line = ('%s\t%.2f\t%.1f\t%.1f' % (filename, coverage_num, length_total, avg_coverage))
        if coverage_num < total_coverage_cutoff or avg_coverage < genome_avg_coverage_cutoff:
            # not qualified
            print(filename, coverage_num, avg_coverage)
            Donor_species[donor_species][0].append(filename)
            temp_line += ('\tnot_qualified\t\n')
        else:
            # qualified
            Donor_species[donor_species][1] += 1
            temp_line += ('\tqualified\t\n')
        Coverage.append(temp_line)

def SNP_check(lines,donor_species,vcf_file_list):
    # CHR, POS, REF, ALT, good assembly, qualified mapping
    lines_set = lines.split('\n')[0].split('\t')
    report_line = ''
    temp_report_line = ['T','T']
    need_curation = 'F'
    REF = lines_set[3]
    allels_set = [REF]
    if '.' not in lines_set[4]:
        allels_set += lines_set[4].split(',')
    Total_alleles = len(allels_set)
    for Subdepth_all in lines_set[9:]:
            Allels_frq = [0, 0, 0, 0]
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
            total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
            Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
            Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
            total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
            total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
            for num_allels in range(0, Total_alleles):
                allels = allels_set[num_allels]
                Subdepth_alleles = int(Subdepth[num_allels])
                if allels in Allels:
                    Allels_frq[Allels[allels]] += Subdepth_alleles
                else:
                    pass
            # find major alt and calculate frequency
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            MLF = Major_ALT[1] / total_sub_depth
            if REF != Major_ALT[0]:
                # wrong assembly
                temp_report_line[0] = 'F'
            if total_sub_depth_forward < Sample_depth_cutoff or \
                    total_sub_depth_reverse < Sample_depth_cutoff or \
                    MLF < Major_alt_freq_cutoff:
                # unqualified mapping
                temp_report_line[1] = 'F'
            if total_sub_depth_forward >= Sample_depth_cutoff2 and \
                    total_sub_depth_reverse >= Sample_depth_cutoff2 and \
                    MLF >= Major_alt_freq_cutoff2:
                # can be curated
                need_curation = 'T'
    if temp_report_line != ['T','T']:
        report_line += '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t\n' %(donor_species,lines_set[0],lines_set[1],
                                                    REF,','.join([ALT for ALT in allels_set if ALT != REF]),
                                                    temp_report_line[0],temp_report_line[1],need_curation,
                                                     Major_ALT[0], MLF,
                                                     total_sub_depth_forward, total_sub_depth_reverse)
        vcf_file_list.append(donor_species + '\t' + lines)
    return report_line

# check major alt
all_vcf_file=glob.glob(os.path.join(output_dir + '/bwa/0','*.flt.snp.vcf'))
vcf_file_report = []
vcf_file_list = []
vcf_file_report.append('donor_species\tCHR\tPOS\tREF\tALT\tAssembly\tMapping_quality\tNeed_curation\tMajor_ALT\tMajor_ALT_frq\tDepth_F\tDepth_R\n')
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split(fastq_name)[0]
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            vcf_file_report.append(SNP_check(lines,donor_species,vcf_file_list))

f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.sum'), 'w')
f1.write(''.join(vcf_file_report))
f1.close()

f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.vcf'), 'w')
f1.write(''.join(vcf_file_list))
f1.close()

# check coverage
Donor_species = dict()
Coverage = []
Coverage.append('genome\treads_alignment_rate\tref_genome_length\taverage_coverage\tquality\t\n')
sh_files = glob.glob(os.path.join(input_script_split_sub, '*.sh'))
for sh_file in sh_files:
    coverage_check(sh_file, Coverage)

f1 = open(os.path.join(input_script, 'SNP_currate1.coverage.sum'), 'w')
f1.write(''.join(Coverage))
f1.close()

################################################### END ########################################################
# step 3 curate genome
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
# set up path
input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/vcf_new'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_currate1'
fastq_name = '_1.fastq'
genome_split = '_g'
species_select = ''
# set up cutoff
end_cutoff = 10 # don't curate the first and last end_cutoff bp of a contig

# function

def loadreport(vcf_file_report,coverage_report):
    assembly_quality = dict()
    assembly_quality_chr = dict()
    report_curate = dict()
    coverage = dict()
    remove_list = []
    try:
        for lines in open(coverage_report,'r'):
            lines_set = lines.split('\t')
            genome = lines_set[0]
            quality = lines_set[4]
            coverage.setdefault(genome,quality)
    except FileNotFoundError:
        pass
    for lines in open(vcf_file_report,'r'):
        if species_select == '' or species_select in lines:
            lines_set = lines.split('\t')
            genome = lines_set[0]
            if coverage.get(genome,'') == 'qualified' or coverage == {}:
                assembly_quality.setdefault(genome,[])
                CHR = lines_set[1]
                POS = lines_set[2]
                REF = lines_set[3]
                ALT_major = lines_set[8]
                Assembly = lines_set[5]
                Mapping_quality = lines_set[6]
                Need_curation = lines_set[7]
                genome_chr = '%s_%s'%(genome,CHR)
                report_curate['%s_%s_%s' % (genome, CHR, POS)] = lines
                if Need_curation == 'T':
                    # can be curated
                    assembly_quality_chr.setdefault(genome_chr,set())
                    assembly_quality_chr[genome_chr].add('__'.join([POS, REF, ALT_major]))
                    assembly_quality[genome].append(CHR)
                elif Assembly == 'F' or Mapping_quality == 'F':
                    # wrong assembly or bad mapping quality
                    assembly_quality_chr.setdefault(genome_chr, set())
                    assembly_quality_chr[genome_chr].add('__'.join([POS, REF, 'N'])) # use ambiguous characters instead
                    assembly_quality[genome].append(CHR)
                else:
                    print(lines)
            else:
                remove_list.append(genome)
    return [assembly_quality,assembly_quality_chr,report_curate]

def checkREF(seq,position,REF):
    return seq[position - 1].upper() == REF.upper()

def causeSNP(seq,position,ALT):
    seq = list(seq)
    seq[position - 1] = ALT
    return ''.join(seq)

def correct_genome(database,assembly_quality_genome,assembly_quality_chr,report_curate):
    Output = []
    Output_report = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        total_length = len(record_seq)
        if record_id in assembly_quality_genome:
            genome_chr = '%s_%s' % (genome, record_id)
            allposition = assembly_quality_chr.get(genome_chr,set())
            for a_position in allposition:
                POS, REF, ALT = a_position.split('__')
                POS = int(POS)
                if POS > end_cutoff and POS < total_length-end_cutoff + 1:
                    # don't curate the first and last end_cutoff bp of a contig
                    if not checkREF(record_seq,int(POS),REF):
                        # wrong POS that needs to be checked
                        print(genome_chr, allposition)
                        print('wrong SNP POS %s %s for %s %s in %s'%(POS,REF,record_seq[POS - 1],
                                                                  record_id,database))
                    record_seq = causeSNP(record_seq, int(POS), ALT)
                    Output_report.append(report_curate['%s_%s_%s' % (genome, record_id, POS)])
        Output.append('>%s\n%s\n'%(record_id,record_seq))
    f1 = open(database + '.corrected.fasta','w')
    f1.write(''.join(Output))
    f1.close()
    return Output_report

# load report
assembly_quality,assembly_quality_chr,report_curate = loadreport(os.path.join(input_script, 'SNP_currate1.assembly.sum'),
           os.path.join(input_script, 'SNP_currate1.coverage.sum'))

# output curated loci
f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.sum.curated'), 'w')
f1.write('donor_species\tCHR\tPOS\tREF\tALT\tAssembly\tMapping_quality\tNeed_curation\tMajor_ALT\tMajor_ALT_frq\tDepth_F\tDepth_R\n')
f1.close()
# correct genomes
for genome in assembly_quality:
    assembly_quality_genome = assembly_quality.get(genome, [])
    if assembly_quality_genome!=[]:
        donor_species = '_'.join(genome.split(genome_split)[:-1])
        donor = donor_species.split('_')[0]
        database = glob.glob(os.path.join(genome_root,'%s*/%s.fasta'%(donor_species,genome)))
        if database!= []:
            database = database[0]
            correct_genome(database, assembly_quality_genome,assembly_quality_chr,report_curate)
        else:
            print('missing input fasta for %s'%(genome))
    else:
        os.system('cp %s %s.corrected.fasta'%(genome,genome))

################################################### END ########################################################
# step 4 all spades map to all curated genomes
import glob
import os

# set up path
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_currate2'
input_script_sub_merge = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_currate2_merge'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
#genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')
genome_dir = ['/scratch/users/anniz44/genomes/donor_species/selected_species/af_Bacteroides_ovatus']
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_currate2'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_currate2_merge'
fastq_name = '.fasta.corrected.fasta'
genome_split = '_g'

try:
    os.mkdir(output_dir_merge)
except IOError:
    pass

try:
    os.mkdir(output_dir_merge+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir_merge+'/bwa/0')
except IOError:
    pass


try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa/0')
except IOError:
    pass

os.system('rm -rf %s %s'%(input_script_sub,input_script_sub_merge))
try:
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(input_script_sub_merge)
except IOError:
    pass

i=0
for folder in genome_dir:
    donor_species = os.path.split(folder)[-1]
    donor = donor_species.split('_')[0]
    if len(donor) == 2: # BN10 donors
        database = glob.glob(os.path.join(folder,'*.all.spades.fasta'))
        sub_samples = []
        sub_samples2 = []
        genomes = glob.glob(os.path.join(folder, '*corrected.fasta'))
        if database!=[] and genomes!=[]:
            database = database[0]
            database2 = database.replace('.all.spades.fasta','.all.spades.fna')
            try:
                f1 = open(database + '.mmi', 'r')
            except IOError:
                os.system('minimap2 -d %s.mmi %s'%(database,database))
            try:
                f1 = open(database2 + '.mmi', 'r')
            except IOError:
                os.system('minimap2 -d %s.mmi %s'%(database2,database2))
            cmds = ''
            for files in genomes:
                donor_species_dir_file = os.path.split(files)[-1]
                print(donor_species_dir_file)
                tempbamoutput = os.path.join(output_dir +'/bwa/0', donor_species_dir_file)
                try:
                    f1 = open(tempbamoutput + '.sorted.bam')
                except IOError:
                    cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                        min(40, 40), database, files, 'samtools', min(40, 40),
                        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
                        tempbamoutput)
                    # fna for genes
                    #cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.fna.bam\n%s sort -@ %s %s.fna.bam -o %s.fna.sorted.bam\n%s index -@ %s %s.fna.sorted.bam\n' % (
                    #    min(40, 40), database2, files, 'samtools', min(40, 40),
                    #    tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
                    #    tempbamoutput)
                    #cmds += 'rm -rf %s.bam %s.fna.bam\n' % (tempbamoutput,tempbamoutput)
                    sub_samples.append(tempbamoutput + '.sorted.bam')
                    #sub_samples2.append(tempbamoutput + '.fna.sorted.bam')
                    i += 1
                    if i % 5 == 0:
                        f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species, int(i / 5))), 'a')
                        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
                        f1.close()
                        cmds = ''
            f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (donor_species,int(i/5))), 'a')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            f1.close()
            # run all mileup
            tempbamoutput = os.path.join(output_dir_merge +'/bwa/0', donor_species)
            cmds = '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call  --threads %s -m > %s.raw.vcf\n' % (
                'bcftools', min(40, 40), database,
                ' '.join(sub_samples), 'bcftools', min(40, 40), tempbamoutput)
            cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
                'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
            cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
                'bcftools', tempbamoutput, tempbamoutput)
            cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
            # fna for genes
            #cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call  --threads %s -m > %s.fna.raw.vcf\n' % (
            #    'bcftools', min(40, 40), database2,
            #    ' '.join(sub_samples2), 'bcftools', min(40, 40), tempbamoutput)
            #cmds += '%s filter --threads %s -s LowQual %s.fna.raw.vcf > %s.fna.flt.vcf \n' % (
            #    'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
            #cmds += '%s view -H -v snps %s.fna.flt.vcf > %s.fna.flt.snp.vcf \n' % (
            #    'bcftools', tempbamoutput, tempbamoutput)
            #cmds += 'rm -rf %s.fna.flt.vcf\n' % (tempbamoutput)
            f1 = open(os.path.join(input_script_sub_merge, '%s.vcf.sh' % donor_species), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            f1.close()
        else:
            print('data missing for %s, database %s genomes %s' % (donor_species, database, genomes))

f1 = open(os.path.join(input_script, 'allnewvcf.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'allnewmergevcf.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub_merge, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
# step 5 filter vcf
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics

# set up path
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_currate2_tree'
input_script_sub_merge = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_currate2_merge'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/'
#genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/*')
genome_dir = ['/scratch/users/anniz44/genomes/donor_species/selected_species/af_Bacteroides_ovatus']
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_currate2'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_currate2_merge'
fastq_name = '.fasta.corrected.fasta'
genome_split = '_g'
Step = 1
deleting_file = []
Cov_dis = 20

# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def vcf_to_txt(lines,output_list,cluster_sub=[]):
    lines_set = lines.split('\n')[0].split('\t')
    if len(lines_set) >9:
        CHR = lines_set[0]
        POS = int(lines_set[1])
        temp_line = []
        temp_line.append(CHR)
        temp_line.append(str(POS))
        i = 9
        for Subdepth_all in lines_set[9:]:
            if (cluster_sub==[] and i not in deleting_set) or i in cluster_sub:
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                temp_line.append(str(total_sub_depth))
            i += 1
        output_list.append('\t'.join(temp_line)+'\n')
    else:
        print(lines)

def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = int(Subdepth[0]) + int(Subdepth[1])
    Subdepth_ALT = int(Subdepth[2]) + int(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return [allels_set[1],1]
    else:
        return [allels_set[0],0]

def outputvcf(output_name):
    vcf_file_filtered = open(vcf_file + '.%s.snp.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.samplename.txt' % (output_name), 'w')
    vcf_file_filtered.write('\t'.join(Sample_name) + '\n')
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.POS.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_POS))
    vcf_file_filtered.close()

def outputcov(output_name,vcf_file_POS_candidate,cluster_sub=[]):
    if len(vcf_file_list) > 0:
        vcf_file_POS_candidate = '\n'.join(vcf_file_POS_candidate)
        vcf_file_POS_candidate_output = ('%s' % (vcf_file_POS_candidate))
        f1 = open(os.path.join(input_script,'grep.temp.txt'),'w')
        f1.write(vcf_file_POS_candidate_output)
        f1.close()
        if '.fna.flt.snp.vcf' in vcf_file:
            cov_file = vcf_file.split('.fna.flt.snp.vcf')[0] + '.fna.raw.vcf'
        else:
            cov_file = vcf_file.split('.flt.snp.vcf')[0] + '.raw.vcf'
        os.system('grep -%s -f %s %s --no-group-separator > %s'% (
            Cov_dis, os.path.join(input_script,'grep.temp.txt'),
            cov_file,
            vcf_file + '.%s.cov.temp' % (output_name)))
        os.system('cat %s | sort | uniq > %s' % (
            vcf_file + '.%s.cov.temp' % (output_name),
            vcf_file + '.%s.uniqcov.temp' % (output_name)))
        for lines in open(vcf_file + '.%s.uniqcov.temp' % (output_name), 'r'):
            if not lines.startswith("#"):
                vcf_to_txt(lines, cov_file_list,cluster_sub)
        os.system('rm -rf %s %s %s' % (vcf_file + '.%s.cov.temp' % (output_name),
                                    vcf_file + '.%s.uniqcov.temp' % (output_name),
                                       os.path.join(input_script, 'grep.temp.txt')))
        vcf_file_filtered = open(vcf_file + '.%s.cov.txt' % (output_name), 'w')
        vcf_file_filtered.write(''.join(cov_file_list))
        vcf_file_filtered.close()

def outputtree(output_name):
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, SNP_alignment[genomename]))
            SNP_alignment_output_parsi.append('S%s    %s\n' % (genomename[-8:], SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.%s.vcf' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list_vcf))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.fasta' % (output_name), 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.parsi.fasta' %(output_name), 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()
    if Step == 2:
        if SNP_alignment_output != []:
            # with SNP
            fasta = vcf_file + '.%s.fasta' % (output_name)
            SNP_tree_cmd.append('python /scratch/users/anniz44/scripts/maffttree/remove.duplicated.py -i %s\n' % (fasta))
            filesize = int(os.path.getsize(fasta))
            if filesize >= 10 ^ 7:  # 10Mb
                SNP_tree_cmd.append(
                    'mafft --nuc --quiet --nofft --retree 2 --maxiterate 0 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                        # 'mafft --nuc --quiet --retree 2 --maxiterate 100 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                        fasta, fasta))  # remove --adjustdirection
                SNP_tree_cmd2.append(
                    '#raxmlHPC -s %s -m GTRGAMMA -n %s.raxml -p 1000 -T 40\n' % (
                        fasta, fasta))
            else:
                SNP_tree_cmd.append(
                    # 'mafft --nuc --quiet --nofft --retree 2 --maxiterate 0 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                    'mafft --nuc --quiet --retree 2 --maxiterate 100 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                        fasta, fasta))  # remove --adjustdirection
                SNP_tree_cmd2.append(
                    '#raxmlHPC -s %s -m GTRGAMMA -n %s.raxml -p 10000 -T 40\n' % (
                        fasta, fasta))
            SNP_tree_cmd.append('FastTree -nt -quiet %s.mafft.align > %s.mafft.align.nwk\n' % (fasta, fasta))
            SNP_tree_cmd2.append(
                '#run_gubbins.py --threads 40  --tree_builder hybrid --use_time_stamp --prefix %s_gubbins --verbose %s.mafft.align\n' % (
                    fasta, fasta))
            SNP_tree_cmd2.append(
                'mv *%s.mafft* %s\n' % (
                    fasta,output_dir_merge + '/bwa/0/tree'))
            SNP_tree_cmd2.append('rm -rf %s.sorted*\n'%(fasta))
        f1 = open(os.path.join(input_script_sub, '%s.tree.all.sh' % (donor_species)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n' + ''.join(SNP_tree_cmd) + ''.join(SNP_tree_cmd2))
        f1.close()

def SNP_seq(seq1, seq2, POS_info,POS_info_CHR,POS_info_CHR_LEN,POS_info_output,G1,G2):
    SNP_total = 0
    j = 0
    POS_DIS = []
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            CHR = POS_info_CHR[i]
            POS = POS_info[i]
            LEN = POS_info_CHR_LEN[CHR]
            if CHR == POS_info_CHR[j]:  # same CHR
                DIS = abs(POS - POS_info[j])
                POS_DIS.append(DIS)  # POS diff
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, DIS, LEN))
            else:  # new CHR
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, 0, LEN))
            j = i
    return SNP_total

def SNP_distance_correct(distanceNEW, distanceREF, SNPREF):
    return int((distanceNEW/distanceREF)*SNPREF)

def find_neighbor(Cluster_SNP,neighbor,Cluster_SNP_set,cluster,Cluster_SNP_set_added):
    if neighbor != []:
        for record_name in neighbor:
            if record_name not in Cluster_SNP_set_added:
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    subneighbor = Cluster_SNP.get(record_name,[])
                    find_neighbor(Cluster_SNP, subneighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)

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

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

def contig_to_gene(CHR, POS):
    all_genes = Mapping_loci.get(CHR,[])
    Reverse_chr = 0
    if all_genes == []:
        # database is a gene database
        codon_start = POS - 1 - int((POS - 1) % 3)
        Ref_seq_chr = Ref_seq.get(CHR, 'None')
        return [CHR, POS, codon_start, Ref_seq_chr, Reverse_chr]
    else:
        # database is a contig database
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

def SNP_check_all(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name,SNP_presence_cutoff,SNP_presence_sample_cutoff,cluster_sub=[]):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    Total_qualify = 0
    Total_qualify_SNP = 0
    Total_qualify_notSNP = 0
    Total_unqualify_alt_freq = 0
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    if cluster_sub!= []:
        lines_set_sub = [lines_set[i] for i in cluster_sub]
        Total_subsample = len(cluster_sub)
        if Total_subsample >= 15:
            SNP_presence_cutoff = 0.33  # for a large group of samples
        if Total_subsample <= 3:
            SNP_presence_cutoff = 1  # for a small group of samples
            SNP_presence_sample_cutoff = 2
    else:
        cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        if Step == 3:
            # re-calculate depth by a subset
            Depth4 = '%s,0,%s,0'%(sum([int(Subdepth_all.split(':')[-1].replace('\n', '').split(',')[0]) for Subdepth_all in lines_set_sub]),
                              sum([int(Subdepth_all.split(':')[-1].replace('\n', '').split(',')[1]) for Subdepth_all in
                                   lines_set_sub]))
        REF, REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    qualify_loci = 0
                    MLF = Major_ALT[1] / total_sub_depth
                    if MLF >= Major_alt_freq_cutoff or qualify_loci == 1:
                        # major alt frequency cutoff
                        Total_qualify += 1
                        # check for qualified SNP
                        if Major_ALT[0] != REF:
                            SNP_seq[-1] = Major_ALT[0]  # unqualified SNP also include in alignment
                            # no specific rule
                            # a qualified SNP
                            temp_snp_line_pass += 'PASS'
                            Total_qualify_SNP += 1
                            SNP.add(genome_order)  # only take qualified SNP as valid SNP
                            SNP_seq[-1] = Major_ALT[0]
                        else:
                            Total_qualify_notSNP += 1
                    else:
                        # major alt frequency low
                        Total_unqualify_alt_freq += 1
            sample_num += 1
        if Total_qualify / Total_subsample >= SNP_presence_cutoff and \
                Total_unqualify_alt_freq <= Poor_MLF_freq_cutoff and\
                Total_qualify >= SNP_presence_sample_cutoff and \
                Total_qualify_SNP >= 1 and Total_qualify_SNP < Total_qualify and\
                Total_qualify_notSNP > 0:
            # -> qualified SNP
            # qualified samples cutoff + unqualified samples cutoff
            # at least 1 qualified SNP
            # calculate NS
            gene_info = contig_to_gene(CHR, POS)
            if gene_info!= []:
                Chr_gene, POS_gene,codon_start,Ref_seq_chr,Reverse_chr  = gene_info
                if Ref_seq_chr != 'None':
                    #  observed NS ratio calculated
                    temp_snp_line_NS= [Chr_gene,str(POS_gene),'']
                    if codon_start <= POS_gene - 1:
                        Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                        Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                        SNP_seq_chr = Ref_seq_chr
                        if len(Ref_seq_codon) == 3:
                            Ref_seq_aa = translate(Ref_seq_codon)[0]
                            temp_snp_line_AA += Ref_seq_aa
                            ALT_set = allels_set
                            ALT_set.remove(REF)
                            for ALT in ALT_set:
                                SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                SNP_seq_aa = translate(SNP_seq_codon)[0]
                                temp_snp_line_AA += SNP_seq_aa
                                temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                temp_snp_line_NS[-1]+=temp_NorS
            # output lines and output major alt
            if Rough == 1: #tune cutoff
                if Total_qualify / Total_subsample >= SNP_presence_cutoff2 and \
                        Total_unqualify_alt_freq <= Poor_MLF_freq_cutoff2 and \
                        Total_qualify >= SNP_presence_sample_cutoff2 \
                        and Total_qualify_SNP >= 1 and Total_qualify_SNP <= Total_qualify - no_SNP_cutoff \
                        and Total_qualify_notSNP >= no_SNP_cutoff:
                    temp_snp_line_pass = 'PASS'
                else:
                    temp_snp_line_pass = 'NOT_PASS'
            else:
                if 'PASS' in temp_snp_line_pass:
                    temp_snp_line_pass = 'PASS'
                else:
                    temp_snp_line_pass = 'RULE'
            if CHR == CHR_old:
                # same CHR
                POS_DIS = abs(POS - POS_old)
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
            else:
                # diff CHR first SNP
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
            POS_old = POS
            CHR_old = CHR
            temp_snp_line.append(CHR)
            temp_snp_line.append(str(POS))
            temp_snp_line.append(REF)
            temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
            vcf_file_list.append('\t'.join(temp_snp_line)+ '\t' +'\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,'\t'.join(temp_snp_line_NS),temp_snp_line_AA))
            vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
            vcf_file_list_vcf.append('\t'.join(lines_set[0:9])+'\t'+'\t'.join(lines_set_sub)+'\n')
            if Rough != 1:
                i = 9
                j = 0
                SNP_alignment[reference_name] += REF
                for genomename in SNP_alignment:
                    if genomename != reference_name:
                        if i in cluster_sub:
                            SNP_alignment[genomename] += SNP_seq[j]
                            j += 1
                        i += 1
    return [CHR_old,POS_old]

# set up output

try:
    os.mkdir(output_dir_merge + '/bwa/0/tree')
except IOError:
    pass

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# set up steps
SNP_cluster = dict()
cluster_set = set()
if Step == 1:
    reference_set = ['']
    outputname_set = ['filteredrough']
elif Step == 2:
    reference_set = ['reference']
    outputname_set = ['filtered']
elif Step == 3:
    reference_set = []
    outputname_set = []
    # load cluster
    for lines in open(os.path.join(input_script, 'SNP_round2.sum'), 'r'):
        if not lines.startswith('donor_species'):
            lines_set = lines.split('\n')[0].split('\t')
            genomes = lines_set[1]
            cluster_type = lines_set[-3]
            cluster_total = int(lines_set[-2])
            if cluster_total > 1:
                SNP_cluster.setdefault(genomes, cluster_type)
                cluster_set.add(cluster_type)
    for cluster_type in cluster_set:
        reference_set.append('refer_%s'%(cluster_type))
        outputname_set.append('filtered.%s' % (cluster_type))
    cluster_set = list(cluster_set)

# set up cutoff super rough
if Step ==1 :
    Rough = 1  # tune cutoff
    SNP_presence_cutoff = 0.33 # avg presence in all samples
    Poor_MLF_freq_cutoff = 0  # not set, the unqualifie samples should be mostly low cov but not two alleles, homologous genes (low major alt freq)
else:
    Rough = 0
    SNP_presence_cutoff = 0.66  # avg presence in all samples
    Poor_MLF_freq_cutoff = 1  # all unqualifie samples should be mostly low cov but not two alleles, homologous genes (low major alt freq)

# unchanged cutoff
SNP_presence_sample_cutoff = 3  # num of samples passing the above criteria
Major_alt_freq_cutoff = 0.9 # major alt freq in a genome, do not allow multiple homolougs genes
no_SNP_cutoff = 1

# set up strict cutoff
SNP_presence_cutoff2 = 0.66 # avg coverage in all samples
SNP_presence_sample_cutoff2 = 3  # num of samples passing the above criteria
Poor_MLF_freq_cutoff2 = 1 # all unqualifie samples should be mostly low cov but not two alleles (low major alt freq)

# cluster cutoff Step 3
SNP_total_cutoff_2 = 50
cluster_cutoff = 2

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir_merge + '/bwa/0','*.flt.snp.vcf'))
for set_num in range(0,len(reference_set)):
    reference_name = reference_set[set_num]
    output_name = outputname_set[set_num]
    for vcf_file in all_vcf_file:
        SNP_presence_cutoff = SNP_presence_cutoff2 # for group of samples  -> used
        SNP_presence_sample_cutoff = SNP_presence_sample_cutoff2
        no_SNP_cutoff = 1
        print(vcf_file)
        Total = 0
        donor_species = os.path.split(vcf_file)[-1].split('.fna.flt.snp.vcf')[0].split('.flt.snp.vcf')[0]
        Sample_name = []
        deleting_set = []
        Ref_seq = dict()
        Mapping = dict()
        Mapping_loci = dict()
        SNP_cluster_donor_species = dict()
        for cluster_type in cluster_set:
            SNP_cluster_donor_species.setdefault(cluster_type,[])
        for lines in open(os.path.join(input_script_sub_merge, '%s.vcf.sh' % (donor_species)), 'r'):
            if lines.startswith('bcftools mpileup '):
                # setup samples
                sample_set = lines.split('all.spades')[1].split('\n')[0].split(' |')[0].split(' ')
                samplenum = 9
                for samples in sample_set[1:]:
                    genomename = os.path.split(samples)[-1].split(fastq_name)[0]
                    Sample_name.append(genomename.replace('.', ''))
                    if genomename in deleting_file:
                        deleting_set.append(samplenum)
                    else:
                        if SNP_cluster!= dict() and genomename in SNP_cluster:
                            SNP_cluster_donor_species[SNP_cluster[genomename]].append(samplenum)
                    samplenum += 1
                break
        # subset samples
        Sample_subset = []
        if Step == 3:
            Sample_subset = SNP_cluster_donor_species.get(cluster_set[set_num], [])
            Sample_name = [Sample_name[i-9] for i in Sample_subset]
        print(Sample_name)
        if Step < 3 or (Sample_subset != [] and len(Sample_subset)>= cluster_cutoff):
            print('running %s' %(donor_species))
            # load database
            database_file = glob.glob(os.path.join(genome_root,
                                                   '%s/*all.spades.fna' % (donor_species)))
            if database_file != []:
                Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file[0])
            vcf_file_all = [vcf_file,vcf_file.replace('.fna.flt.snp.vcf','.flt.snp.vcf')]
            for vcf_file in vcf_file_all:
                print(vcf_file)
                SNP_tree_cmd = []
                SNP_tree_cmd2 = []
                vcf_file_list = []
                vcf_file_list_vcf = []
                vcf_file_POS = []
                vcf_file_POS_candidate = set()
                SNP_alignment = dict()
                SNP_alignment.setdefault(reference_name, '')
                cov_file_list = []
                CHR_old = ''
                POS_old = 0
                for genomename in Sample_name:
                    SNP_alignment.setdefault(genomename, '')
                for lines in open(vcf_file, 'r'):
                    if not lines.startswith("#"):
                        lines_set = lines.split('\n')[0].split('\t')
                        CHR = lines_set[0]
                        POS = int(lines_set[1])
                        Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                        if Total == 0:
                            Total = len(lines_set) - 9 - len(deleting_set)
                            if Total >= 15:
                                SNP_presence_cutoff = 0.33  # for a large group of genomes
                            elif Total in [3,4]:
                                SNP_presence_cutoff = 1  # for a small group of genomes
                                SNP_presence_sample_cutoff = 2
                            elif Total in [1,2]:
                                # compare to reference
                                no_SNP_cutoff = 0
                                SNP_presence_cutoff = 1
                                SNP_presence_sample_cutoff = 1
                        if Depth / Total >= SNP_presence_cutoff:
                            # average depth in all samples cutoff
                            if "INDEL" not in lines_set[7] \
                                    and (lines_set[6] != 'LowQual'):
                                if Step == 1:
                                    CHR_old, POS_old =  SNP_check_all(lines_set, '',
                                                                      CHR_old,POS_old,reference_name,
                                                                      SNP_presence_cutoff,SNP_presence_sample_cutoff)
                                elif Step == 2:
                                    CHR_old, POS_old = SNP_check_all(lines_set, '',
                                                                     CHR_old,POS_old,reference_name,
                                                                     SNP_presence_cutoff,SNP_presence_sample_cutoff)
                                else:
                                    CHR_old, POS_old = SNP_check_all(lines_set, '',
                                                                     CHR_old, POS_old, reference_name,
                                                                     SNP_presence_cutoff, SNP_presence_sample_cutoff,
                                                                     Sample_subset)
                outputvcf(output_name)
                if Step != 1:
                    # step 2 and 3 output fasta
                    outputtree(output_name)
                if Step != 2:
                    # step 1 and 3 output coverage, step 2 use step 1 coverage
                    try:
                        f1 = open(vcf_file + '.%s.cov.txt' % (output_name), 'r')
                    except IOError:
                        outputcov(output_name,list(vcf_file_POS_candidate),Sample_subset)

# run clustering and mafft tree
if Step == 2:
    # sum up SNP
    SNP_cutoff = []
    SNP_pair = []
    SNP_pair.append('G1\tG2\tSNP\t\n')
    POS_info_output = []
    POS_info_output.append('G1\tG2\tCHR\tPOS\tDIS\tLEN\t\n')
    Fasta_SNP = glob.glob(os.path.join(output_dir_merge + '/bwa/0', '*.filtered.fasta'))
    tree_distance_output = []
    for fasta in Fasta_SNP:
        fasta_name = os.path.split(fasta)[-1]
        donor_species = fasta_name.split('.flt.snp.vcf.filtered.fasta')[0]
        POS_file = glob.glob(os.path.join(output_dir_merge + '/bwa/0', donor_species + '.flt.snp.vcf.filtered.POS.txt'))[0]
        print(donor_species)
        tree_distance = dict()
        Seq_list = dict()
        Ref = ''
        max_cluster_diff = 0
        POS_info = []
        POS_info_CHR = []
        POS_info_CHR_LEN = dict()
        Cluster_SNP = dict()
        Cluster_SNP_set = dict()
        Cluster_SNP_set_added = set()
        # load SNP POS info
        for lines in open(POS_file, 'r'):
            CHR = lines.split('\t')[0]
            POS_info.append(int(lines.split('\t')[1]))
            POS_info_CHR.append(CHR)
            POS_info_CHR_LEN.setdefault(CHR, 0)
            POS_info_CHR_LEN[CHR] = max(POS_info_CHR_LEN[CHR], int(lines.split('\t')[1]))
        # load genome SNP fasta and calculate pair-wise SNPs
        for record in SeqIO.parse(fasta, 'fasta'):
            record_name = str(record.id)
            if 'reference' not in record_name:
                new_center = 1
                record_seq = str(record.seq)
                if Ref == '':
                    # set upt the first seq as ref
                    Ref = record_seq
                    REF_name = record_name
                    new_center = 0
                    SNP_total = 0
                    SNP_total_length = len(Ref)
                    #SNP_total_cutoff = max(SNP_total_cutoff_ratio * SNP_total_length, SNP_total_cutoff_2)
                    SNP_total_cutoff = SNP_total_cutoff_2
                else:
                    for record_before in Seq_list:
                        SNP_total = SNP_seq(Seq_list[record_before], record_seq, POS_info, POS_info_CHR,
                                            POS_info_CHR_LEN,
                                            POS_info_output, record_before, record_name)
                        SNP_pair.append('%s\t%s\t%s\t\n' % (record_before, record_name, SNP_total))
                        if SNP_total <= SNP_total_cutoff:
                            Cluster_SNP.setdefault(record_before, [])
                            Cluster_SNP[record_before].append(record_name)
                            max_cluster_diff = max(max_cluster_diff, SNP_total)
                    SNP_total = SNP_seq(Ref, record_seq, POS_info, POS_info_CHR, POS_info_CHR_LEN,
                                        POS_info_output, REF_name, record_name)
                    SNP_pair.append('%s\t%s\t%s\t\n' % (REF_name, record_name, SNP_total))
                    if SNP_total <= SNP_total_cutoff:
                        Cluster_SNP.setdefault(REF_name,[])
                        Cluster_SNP[REF_name].append(record_name)
                        max_cluster_diff = max(max_cluster_diff,SNP_total)
                tree_distance.setdefault(record_name, SNP_total)
                Seq_list.setdefault(record_name, record_seq)
        cluster = 0
        # cluster genomes by SNP distance
        for record_name in Cluster_SNP:
            neighbor = Cluster_SNP.get(record_name,[])
            if neighbor != [] and record_name not in Cluster_SNP_set_added:
                cluster += 1
                Cluster_SNP_set.setdefault(cluster,set())
                Cluster_SNP_set[cluster].add(record_name)
                Cluster_SNP_set_added.add(record_name)
                find_neighbor(Cluster_SNP, neighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)
        # output single genome
        for record_name in Seq_list:
            if record_name not in Cluster_SNP_set_added:
                cluster += 1
                Cluster_SNP_set.setdefault(cluster, set())
                Cluster_SNP_set[cluster].add(record_name)
                Cluster_SNP_set_added.add(record_name)
        Sub_cluster = len(Cluster_SNP_set)
        print(Cluster_SNP_set,fasta_name)
        for cluster in Cluster_SNP_set:
            tree_name_list = Cluster_SNP_set[cluster]
            tree_SNP_count = 'cluster%s' % (cluster)
            for tree_name in tree_name_list:
                tree_distance_output.append(
                    '%s\t%s\t%s\t%s\t%s\t\n' % (
                        donor_species, tree_name, tree_distance[tree_name], tree_SNP_count, Sub_cluster))
        SNP_cutoff.append('%s\t%s\t%s\t%s\t\n' % (donor_species, max_cluster_diff, SNP_total_cutoff, SNP_total_length))
    os.system('#mv %s %s' % (os.path.join(input_script, 'SNP_round2.sum'),
                            os.path.join(input_script, 'SNP_round2.old.sum')))
    f1 = open(os.path.join(input_script, 'SNP_round2.sum'), 'w')
    f1.write('donor_species\tGenome\tSNP_total\tcluster1\tsubcluster\t\n%s' % (''.join(tree_distance_output)))
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round2.sum.cutoff'), 'w')
    f1.write('donor_species\tmax_cluster_diff\tSNP_cutoff\tSNP_total_len\t\n%s' % (''.join(SNP_cutoff)))
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round2.POS.sum'), 'w')
    f1.write('%s' % (''.join(POS_info_output)))
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_round2.allpair.sum'), 'w')
    f1.write('Genome1\tGenome2\tSNP_total\t\n%s' % (''.join(SNP_pair)))
    f1.close()
    # run mafft tree
    f1 = open(os.path.join(input_script, 'allround2alltree.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.tree.all.sh')):
        #f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
        f1.write('sh %s %s.out\n' % (sub_scripts, sub_scripts))
    f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# compare genome call SNPs VS WGS call SNPs
# Question: why genome mapping shows SNPs that have no alternative alleles in raw.vcf of WGS mapping?
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics

# set up path
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate'
# bowtie WGS
vcf_1 = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_currate2/bowtie2.all.flt.snp.vcf'
# minimap genome
vcf_2 = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_currate2/test2.all.flt.snp.vcf'
# bowtie WGS
vcf_1 = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_merge_round2/bwa/0/af_Bacteroides_ovatus.all.flt.snp.vcf.filteredrough.snp.txt'
# minimap genome
vcf_2 = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_currate2_merge/bwa/0/af_Bacteroides_ovatus.flt.snp.vcf.filteredrough.snp.txt'

# function
def load_vcf(vcf_file):
    vcf_input = []
    for lines in open(vcf_file):
        lines_set = lines.split('\n')[0].split('\t')
        CHR = lines_set[0]
        POS = lines_set[1]
        vcf_input.append('%s\t%s\t\n'%(CHR,POS))
    return vcf_input

def compare_vcf(vcf_input1, vcf_input2, vcf_file1,vcf_file2):
    vcf_1_diff = []
    temp_output = os.path.join(input_script,'grep.temp.txt')
    for CHRPOS in vcf_input1:
        if CHRPOS not in vcf_input2:
            vcf_1_diff.append(CHRPOS)
    if len(vcf_1_diff) > 0:
        f1 = open(temp_output, 'w')
        f1.write(''.join(vcf_1_diff))
        f1.close()
        os.system('grep -T -f %s %s %s --no-group-separator > %s' % (
            temp_output,
            vcf_file1.split('.flt.snp.vcf')[0] + '.raw.vcf',
            vcf_file2.split('.flt.snp.vcf')[0] + '.raw.vcf',
            vcf_file1 + '.diff'))
        os.system('sort -k3 -n %s | sort -k2 > %s'%
                  (vcf_file1 + '.diff',vcf_file1 + '.diff.sort')
                  )
        os.system('rm -rf %s'%(vcf_file1 + '.diff'))

vcf_input1 = load_vcf(vcf_1)
vcf_input2 = load_vcf(vcf_2)
compare_vcf(vcf_input1, vcf_input2, vcf_1,vcf_2) # vcf diff in 1 not in 2
compare_vcf(vcf_input2, vcf_input1, vcf_2,vcf_1) # vcf diff in 2 not in 1

################################################### END ########################################################
################################################### SET PATH ########################################################
# step 1 modelling SNPs and test 2 methods
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import random

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_test_step1'
input_script_sub2 = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_test_step2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test'
fastq_name = '_1.fastq'

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# Set up N or S
N_S_set = dict()
N_S_set['N']=0
N_S_set['S']=1
purines=['A','G']
pyrimidines=['C','T']
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
# Set up mutation rate
mut_rate = 1/1000/100 # mutation rate 1 per 10000 bp
mut_time = 1 # generate mut_time strains with SNPs -> a population with mut_time + 1 strains
mut_set = [10,25,50,100,250,500,1000,2500]
if len(mut_set)!= 0 :
    mut_time=len(mut_set)

cause_SNP = False
minor_allele_freq_cutoff = 0.2
pop_size = 30 #not used
fastq_length = 150
SNP_times = 1000 # snp reads depth as 1500, both forward and reverse

try:
    os.mkdir(output_dir + '/data')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa/0')
except IOError:
    pass

try:
    os.mkdir(input_script)
except IOError:
    pass

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(input_script_sub2)
except IOError:
    pass

# function
def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

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

def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        if ALT_frq > 0:
            ALT_set.setdefault(ALT_frq, set())
            ALT_set[ALT_frq].add(alleles)
            ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def transitions(REF,ALT):
    if REF in pyrimidines:
        REF = complement[REF]
        ALT = complement[ALT]
    return '%s-%s'%(REF,ALT)

def contig_to_gene(CHR, POS):
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

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        Ref_seq.setdefault(record_id, record_seq)
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

def simulatefastqsub(seq,loci1,loci2):
    loci1=int(loci1)
    loci2=int(loci2)+1
    temp1 = '@SL-NXA:HNMMVBGX2170612:HNMMVBGX2:4:23612:%s:%s 1:N:0:CAGAGAGGATGTCAAT' % (loci1, loci2)
    temp2 = '@SL-NXA:HNMMVBGX2170612:HNMMVBGX2:4:23612:%s:%s 2:N:0:CAGAGAGGATGTCAAT' % (loci1, loci2)
    seq1 = seq[loci1:loci2]
    seq2 = Seq(seq[loci1:loci2])
    seq2 = str(seq2.reverse_complement())
    seq1 = '%s\n%s\n%s\n%s\n' %(temp1,seq1,'+',
                                'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF')
    seq2 = '%s\n%s\n%s\n%s\n' %(temp2,seq2,'+',
                                'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF')
    return(seq1,seq2)

def simulatefastq(seq,position,total_length):
    if position <= fastq_length/2:
        return simulatefastqsub(seq, 0, fastq_length - 1)
    elif position >= total_length - fastq_length/2:
        return simulatefastqsub(seq, total_length-fastq_length, total_length-1)
    else:
        return simulatefastqsub(seq, position - fastq_length / 2, position + fastq_length / 2 -1)

def modelSNP(seq,Chr,num_mut_chr):
    total_length = len(seq)
    position_set = random.choices(range(0, total_length-1), k=num_mut_chr)
    seq = list(seq)
    seq1 = ''
    seq2 = ''
    SNP_output = []
    for position in position_set:
        gene_info = contig_to_gene(Chr, position)
        REF = seq[position]
        temp_ALT = ['A', 'T', 'G', 'C']
        try:
            temp_ALT.remove(REF)
        except ValueError:
            pass
        ALT = random.choices(temp_ALT, k=1)[0]
        seq[position] = ALT
        fastqsimulated = simulatefastq(''.join(seq), position, total_length)
        seq1 += fastqsimulated[0]
        seq2 += fastqsimulated[1]
        temp_line = [Chr,str(position+1),REF,ALT,'Other','None']
        if gene_info != []:
            # a gene
            Chr_gene, position_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
            codon_start = position_gene - 1 - int((position_gene - 1) % 3)
            if codon_start <= position_gene - 1:
                Ref_seq_chr = Ref_seq[Chr_gene]
                SNP_seq_chr = Ref_seq_chr
                Ref_seq_chr = causeSNP(Ref_seq_chr, position_gene, REF,Reverse_chr)
                Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                if len(Ref_seq_codon) == 3:
                    Ref_seq_aa = translate(Ref_seq_codon)[0]
                    SNP_seq_chr = causeSNP(SNP_seq_chr, position_gene, ALT, Reverse_chr)
                    SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                    SNP_seq_aa = translate(SNP_seq_codon)[0]
                    temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                    temp_line[-1]=''.join([Ref_seq_aa,SNP_seq_aa])
                    temp_line[-2]=temp_NorS
        SNP_output.append('\t'.join(temp_line)+'\n')
    return [''.join(seq),seq1, seq2, SNP_output]

def run_vcf(files,files2,genome_file,database,temp_folder,copy_num,copy_num2,copy_num3):
    # generate code
    # for fastq files
    tempbamoutput = files
    # pair end
    #cmds = 'bowtie2' + ' -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
    #    min(40, 40), database, files, files2, 'samtools', min(40, 40),
    #    tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
    #    tempbamoutput)
    # single end
    cmds = 'bowtie2' + ' -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --threads %s -x %s %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, files, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    sample_set1 = ['%s.sorted.bam'%(tempbamoutput)]
    while copy_num > 1:
        newbam = '%s.%s.sorted.bam'%(tempbamoutput,copy_num)
        cmds += 'cp %s.sorted.bam %s\n'%(tempbamoutput, newbam)
        copy_num -= 1
        sample_set1.append(newbam)
    # for genome files
    tempbamoutput = genome_file
    cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, genome_file, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    sample_set2 = ['%s.sorted.bam' % (tempbamoutput)]
    while copy_num2 > 1:
        newbam = '%s.%s.sorted.bam' % (tempbamoutput, copy_num2)
        cmds += 'cp %s.sorted.bam %s\n' % (tempbamoutput, newbam)
        copy_num2 -= 1
        sample_set2.append(newbam)
    # for curated genome
    genome_file = genome_file + '.corrected.fasta'
    tempbamoutput = genome_file
    cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, genome_file, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    sample_set3 = ['%s.sorted.bam' % (tempbamoutput)]
    while copy_num3 > 1:
        newbam = '%s.%s.sorted.bam' % (tempbamoutput, copy_num3)
        cmds += 'cp %s.sorted.bam %s\n' % (tempbamoutput, newbam)
        copy_num3 -= 1
        sample_set3.append(newbam)
    return [cmds,sample_set1,sample_set2,sample_set3]

def run_vcf_curate(files,files2,database,tempbamoutput):
    # generate code
    cmds = 'rm -rf %s.fai\nbowtie2-build %s %s\n' % (database,database, database)
    cmds += 'bowtie2' + ' -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, files, files2, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    cmds += 'samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
        tempbamoutput, tempbamoutput)
    cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s.sorted.bam  | %s call  --threads %s -m > %s.raw.vcf\n' % (
        'bcftools', min(40, 40), database,
        tempbamoutput, 'bcftools', min(40, 40), tempbamoutput)
    cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
        'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
    cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
        'bcftools', tempbamoutput, tempbamoutput)
    cmds += 'rm -rf %s.*.unaligned\n' % (tempbamoutput)
    cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
    return cmds

def merge_sample(database,tempbamoutputWGS,tempbamoutputGenome,samplesetWGS,samplesetGenome,samplesetGenomecorrect):
    # WGS
    cmds = '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
        'bcftools', min(40, 40), database,
        ' '.join(samplesetWGS), 'bcftools', min(40, 40), tempbamoutputWGS)
    cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
        'bcftools', min(40, 40), tempbamoutputWGS, tempbamoutputWGS)
    cmds += '%s view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
        'bcftools', tempbamoutputWGS, tempbamoutputWGS)
    cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutputWGS)
    # genomes
    cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
        'bcftools', min(40, 40), database,
        ' '.join(samplesetGenome), 'bcftools', min(40, 40), tempbamoutputGenome)
    cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
        'bcftools', min(40, 40), tempbamoutputGenome, tempbamoutputGenome)
    cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
        'bcftools', tempbamoutputGenome, tempbamoutputGenome)
    cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutputGenome)
    # corrected genomes
    tempbamoutputGenome = tempbamoutputGenome + '.corrected'
    cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
        'bcftools', min(40, 40), database,
        ' '.join(samplesetGenomecorrect), 'bcftools', min(40, 40), tempbamoutputGenome)
    cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
        'bcftools', min(40, 40), tempbamoutputGenome, tempbamoutputGenome)
    cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
        'bcftools', tempbamoutputGenome, tempbamoutputGenome)
    cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutputGenome)
    return cmds

def modelSNPall(Input_seq,Input_id,num_mut):
    Output = []
    Output_SNP = []
    seq1all = ''
    seq2all = ''
    chr_set = random.sample(Input_id*num_mut, k=num_mut)
    unique_chr_set = list(set(chr_set))
    for chr in Input_id:
        if chr in unique_chr_set:
            # mutated
            num_mut_chr = chr_set.count(chr)
            newseq, seq1, seq2, newoutput = modelSNP(Input_seq[chr], chr,num_mut_chr)
            Output_SNP += newoutput
            seq1all += seq1
            seq2all += seq2
        else:
            # not mutated
            newseq = Input_seq[chr]
        Output.append('>%s\n%s\n' % (chr,
                                     newseq))
    # output mutated genome
    output_fasta = os.path.join(output_dir, 'test.%s.SNP.fasta' % (num_mut))
    f1 = open(output_fasta, 'w')
    f1.write(''.join(Output))
    f1.close()
    f1 = open(output_fasta + '.snp.txt', 'w')
    f1.write(''.join(Output_SNP))
    f1.close()
    # output mutated fastq
    f1 = open(files + '.temp', 'w')
    f1.write(''.join(seq1all * SNP_times))
    f1.write(''.join(seq2all * SNP_times))
    f1.close()
    f1 = open(files2 + '.temp', 'w')
    f1.write(''.join(seq1all * SNP_times))
    f1.write(''.join(seq2all * SNP_times))
    f1.close()
    os.system('cat %s %s > %s' % (os.path.join(output_dir, 'data/test.fasta1.fq'),
                                  files + '.temp', files))
    os.system('cat %s %s > %s' % (os.path.join(output_dir, 'data/test.fasta2.fq'),
                                  files2 + '.temp', files2))
    os.system('rm -rf %s %s'%(files + '.temp',files2 + '.temp'))
    return 'done'

# load database
database = os.path.join(output_dir,'test.fasta')
database_file = os.path.join(output_dir,'test.fna')
os.system('#prodigal -q -i %s -d %s'%(database,database_file))
Ref_seq, Mapping,Mapping_loci,Reverse = loaddatabase(database_file)
Input_seq = dict()
Input_id = []
for record in SeqIO.parse(database, 'fasta'):
    seq = str(record.seq)
    seq_length = len(seq)
    if seq_length >= 8000:
        Input_seq.setdefault(str(record.id), seq)
        Input_id.append(str(record.id))

# cause SNP
# call SNPs for database as ref bam
fastq_file = os.path.join(output_dir, 'data/test.fasta')
files = fastq_file + '1.fq'
files2 = fastq_file + '2.fq'
genome_file = os.path.join(output_dir, 'data/test.genome.fasta')
temp_folder = os.path.join(output_dir, 'data/test.genome')
# call SNPs by WGS and curated genome
results = run_vcf(files, files2, genome_file, database, temp_folder, 1, 1, 1)
cmds = results[0]
refWGS = results[1]
refGenome = results[2]
refGenomecorrect = results[2]  # not corrected genome, use ref genome
f1 = open(os.path.join(input_script_sub2, 'test.0.vcf.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
f1.close()
cmdsmerge = ''
# cause SNP
while mut_time > 0:
    num_mut = mut_set[mut_time - 1]
    # set up simulated genome and fastq
    output_fasta = os.path.join(output_dir, 'test.%s.SNP.fasta' % (num_mut))
    fastq_file = os.path.join(output_dir, 'data/test.%s.SNP.fasta' % (num_mut))
    files = fastq_file + '1.fq'
    files2 = fastq_file + '2.fq'
    # set up merge output
    tempbamoutputWGS = os.path.join(output_dir, 'bwa/0/test.%s.all.SNP.fasta1.fq'%(num_mut))
    tempbamoutputGenome = os.path.join(output_dir, 'bwa/0/test.%s.all.SNP.fasta.genome.fasta'%(num_mut))
    # set up merge files
    samplesetWGS = []
    samplesetGenome = []
    samplesetGenomecorrect = []
    # cause SNP
    if cause_SNP:
        # simulate fastq files for mutated strains
        modelSNPall(Input_seq, Input_id, num_mut)
    mutated_genome = os.path.join(output_dir, 'test.%s.SNP.fasta'%(num_mut))
    genome_file = os.path.join(output_dir, 'data/test.%s.SNP.fasta.genome.fasta' % (num_mut))
    genome_file_fastq = os.path.join(output_dir, 'data/test.%s.SNP.fasta.genome.fasta' % (num_mut))
    temp_folder = os.path.join(output_dir, 'data/test.%s.SNP.fasta.genome' % (num_mut))
    cmds = ''
    #cmds = 'art_illumina -ss HS25 -sam -i %s -p -l 150 -f 50 -m 200 -s 10 -o %s\n' \
    #       % (mutated_genome, genome_file_fastq)
    #cmds += 'spades.py --careful -1 %s -2 %s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' % \
    #        (genome_file_fastq + '1.fq', genome_file_fastq + '2.fq', temp_folder)
    #cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_folder, genome_file)
    # mapping assembly to simulated WGS to curate genome
    cmds += run_vcf_curate(genome_file_fastq + '1.fq' , genome_file_fastq + '2.fq', genome_file, genome_file + '.curate')
    f1 = open(os.path.join(input_script_sub, 'test.%s.vcf.sh' % (mut_time)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
    f1.close()
    # call SNPs by WGS and curated genome
    results = run_vcf(files, files2, genome_file, database, temp_folder,1,1,1)
    cmds = results[0]
    samplesetWGS += results[1]
    samplesetGenome += results[2]
    samplesetGenomecorrect += results[3]
    f1 = open(os.path.join(input_script_sub2, 'test.%s.vcf.sh' % (mut_time)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
    f1.close()
    samplesetWGS += refWGS
    samplesetGenome += refGenome
    samplesetGenomecorrect += refGenomecorrect
    cmdsmerge += merge_sample(database, tempbamoutputWGS, tempbamoutputGenome, samplesetWGS, samplesetGenome,
                        samplesetGenomecorrect)
    mut_time -= 1

f1 = open(os.path.join(input_script, 'alltestvcfstep1.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'alltestvcfstep2.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub2, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'alltestvcfstep3.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmdsmerge)))
f1.write('#rm -rf %s/*.bam*\n' % (os.path.join(output_dir, 'data')))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# run alltestvcfstep1.sh
# step 2 check SNPs
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
# set up path

input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_test_step1/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_test_step1'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/data'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/data')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/data'
fastq_name = '.curate.flt.snp.vcf'
genome_split = '_g'


# set up cutoff
# good mapping
Major_alt_freq_cutoff = 0.7 # major alt freq in a sample
Sample_depth_cutoff = 3 # both forward and reverse reads cutoff in a sample

# good coverage
total_coverage_cutoff = 0.8 # at least X reads map to its original genome
genome_avg_coverage_cutoff = 10 # genome average coverage cutoff

Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def coverage_check(sh_file,Coverage):
    Coverage1 = []
    Coverage2 = []
    os.system('rm -rf %s/temp.sh %s/temp.err'%(input_script_split_sub,input_script_split_sub))
    os.system('grep \"bowtie2\" %s > %s/temp.sh' % (sh_file, input_script_split_sub))
    os.system('grep \"overall alignment rate\" %s.err > %s/temp.err' % (sh_file, input_script_split_sub))
    # load fastq name
    for lines in open('%s/temp.sh' % (input_script_split_sub), 'r'):
        try:
            filename = lines.split('>')[1].split('\n')[0].split(fastq_name + '.bam')[0]
            filename = os.path.split(filename)[-1]
            Coverage1.append(filename)
        except IndexError:
            pass
    # load reads alignment rate
    for lines in open('%s/temp.err' % (input_script_split_sub), 'r'):
        Coverage2.append(float(lines.split('%')[0]))
    i = 0
    for filename in Coverage1:
        # load average coverage of ref genome
        tempbamoutput = os.path.join(output_dir, filename + fastq_name)
        coverage_file = glob.glob(tempbamoutput + '.sorted.bam.avgcov')
        if coverage_file == []:
            cmds = 'samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
                tempbamoutput, tempbamoutput)
            os.system(cmds)
            coverage_file = tempbamoutput + '.sorted.bam.avgcov'
        else:
            coverage_file = coverage_file[0]
        length_total = 0
        coverage_total = 0
        donor_species = filename.split(genome_split)[0]
        Donor_species.setdefault(donor_species, [[], 0])
        try:
            for lines in open(coverage_file, 'r'):
                length_contig = float(lines.split('\t')[1])
                length_total += length_contig
                coverage_contig = float(lines.split('\t')[2])
                coverage_total += coverage_contig * length_contig
            coverage_num = Coverage2[i] / 100
            avg_coverage = coverage_total / length_total
            i += 1
        except IOError:
            coverage_num = total_coverage_cutoff - 0.1
            avg_coverage = genome_avg_coverage_cutoff - 1
        temp_line = ('%s\t%.2f\t%.1f\t%.1f' % (filename, coverage_num, length_total, avg_coverage))
        if coverage_num < total_coverage_cutoff or avg_coverage < genome_avg_coverage_cutoff:
            # not qualified
            print(filename, coverage_num, avg_coverage)
            Donor_species[donor_species][0].append(filename)
            temp_line += ('\tnot_qualified\t\n')
        else:
            # qualified
            Donor_species[donor_species][1] += 1
            temp_line += ('\tqualified\t\n')
        Coverage.append(temp_line)

def SNP_check(lines,donor_species,vcf_file_list):
    # CHR, POS, REF, ALT, good assembly, qualified mapping
    lines_set = lines.split('\n')[0].split('\t')
    report_line = ''
    temp_report_line = ['T','T']
    need_curation = 'F'
    REF = lines_set[3]
    allels_set = [REF]
    if '.' not in lines_set[4]:
        allels_set += lines_set[4].split(',')
    Total_alleles = len(allels_set)
    for Subdepth_all in lines_set[9:]:
            Allels_frq = [0, 0, 0, 0]
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
            total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
            Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
            Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
            total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
            total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
            for num_allels in range(0, Total_alleles):
                allels = allels_set[num_allels]
                Subdepth_alleles = int(Subdepth[num_allels])
                if allels in Allels:
                    Allels_frq[Allels[allels]] += Subdepth_alleles
                else:
                    pass
            # find major alt and calculate frequency
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            MLF = Major_ALT[1] / total_sub_depth
            if REF != Major_ALT[0]:
                # wrong assembly
                temp_report_line[0] = 'F' # F: not good assembly
            if MLF < Major_alt_freq_cutoff or total_sub_depth_forward < Sample_depth_cutoff or \
                    total_sub_depth_reverse < Sample_depth_cutoff:
                temp_report_line[1] = 'F'  # F: bad mapping
            if temp_report_line[1] == 'T' and temp_report_line[0] == 'F':
                # not good assembly but good mapping quality can be curated
                need_curation = 'T'  # T: need curation
    if temp_report_line != ['T','T']:
        report_line += '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t\n' %(donor_species,lines_set[0],lines_set[1],
                                                    REF,','.join([ALT for ALT in allels_set if ALT != REF]),
                                                    temp_report_line[0],temp_report_line[1],need_curation,
                                                     Major_ALT[0], MLF,
                                                     total_sub_depth_forward, total_sub_depth_reverse)
        vcf_file_list.append(donor_species + '\t' + lines)
    return report_line

# check major alt
all_vcf_file=glob.glob(os.path.join(output_dir,'*%s'%(fastq_name)))
vcf_file_report = []
vcf_file_list = []
vcf_file_report.append('donor_species\tCHR\tPOS\tREF\tALT\tAssembly\tMapping_quality\tNeed_curation\tMajor_ALT\tMajor_ALT_frq\tDepth_F\tDepth_R\n')
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split(fastq_name)[0]
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            vcf_file_report.append(SNP_check(lines,donor_species,vcf_file_list))

f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.sum'), 'w')
f1.write(''.join(vcf_file_report))
f1.close()

f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.vcf'), 'w')
f1.write(''.join(vcf_file_list))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# step 3 curate genome
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
# set up path

input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_test_step1/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_test_step1'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/data'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/data')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test'
fastq_name = '.curate.flt.snp.vcf'
genome_split = '_g'

# set up cutoff
end_cutoff = 10 # don't curate the first and last end_cutoff bp of a contig
# function

def loadreport(vcf_file_report,coverage_report):
    assembly_quality = dict()
    assembly_quality_chr = dict()
    report_curate = dict()
    coverage = dict()
    remove_list = []
    try:
        for lines in open(coverage_report,'r'):
            lines_set = lines.split('\t')
            genome = lines_set[0]
            quality = lines_set[4]
            coverage.setdefault(genome,quality)
    except FileNotFoundError:
        pass
    for lines in open(vcf_file_report,'r'):
        lines_set = lines.split('\t')
        genome = lines_set[0]
        if coverage.get(genome,'') == 'qualified' or coverage == {}:
            assembly_quality.setdefault(genome,[])
            CHR = lines_set[1]
            POS = lines_set[2]
            REF = lines_set[3]
            ALT_major = lines_set[8]
            Assembly = lines_set[5]
            Mapping_quality = lines_set[6]
            Need_curation = lines_set[7]
            genome_chr = '%s_%s'%(genome,CHR)
            report_curate['%s_%s_%s' % (genome, CHR, POS)] = lines
            if Need_curation == 'T':
                # can be curated
                assembly_quality_chr.setdefault(genome_chr,[])
                assembly_quality_chr[genome_chr].append([POS, REF, ALT_major])
                assembly_quality[genome].append(CHR)
            elif Assembly == 'F' or Mapping_quality == 'F':
                # wrong assembly or bad mapping quality
                assembly_quality_chr.setdefault(genome_chr, [])
                assembly_quality_chr[genome_chr].append([POS, REF, 'N']) # use ambiguous characters instead
                assembly_quality[genome].append(CHR)
            else:
                print(lines)
        else:
            remove_list.append(genome)
    return [assembly_quality,assembly_quality_chr,report_curate]

def checkREF(seq,position,REF):
    return seq[position - 1] == REF

def causeSNP(seq,position,ALT):
    seq = list(seq)
    seq[position - 1] = ALT
    return ''.join(seq)

def correct_genome(database,assembly_quality_genome,assembly_quality_chr,report_curate):
    Output = []
    Output_report = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        total_length = len(record_seq)
        if record_id in assembly_quality_genome:
            genome_chr = '%s_%s' % (genome, record_id)
            allposition = assembly_quality_chr.get(genome_chr,[])
            for a_position in allposition:
                POS, REF, ALT = a_position
                POS = int(POS)
                if POS > end_cutoff and POS < total_length-end_cutoff + 1:
                    # don't curate the first and last end_cutoff bp of a contig
                    if checkREF(record_seq,int(POS),REF):
                        record_seq = causeSNP(record_seq, int(POS),ALT)
                        Output_report.append(report_curate['%s_%s_%s' % (genome, record_id,POS)])
                    else:
                        print(genome_chr,allposition)
                        print('wrong SNP POS %s %s for %s in %s'%(POS,REF,record_id,database))
        Output.append('>%s\n%s\n'%(record_id,record_seq))
    f1 = open(database + '.corrected.fasta','w')
    f1.write(''.join(Output))
    f1.close()
    f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.sum.curated'), 'a')
    f1.write(''.join(Output_report))
    f1.close()

# load report
assembly_quality,assembly_quality_chr,report_curate = loadreport(os.path.join(input_script, 'SNP_currate1.assembly.sum'),
           os.path.join(input_script, 'SNP_currate1.coverage.sum'))

# output curated loci
f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.sum.curated'), 'w')
f1.write('donor_species\tCHR\tPOS\tREF\tALT\tAssembly\tMapping_quality\tNeed_curation\tMajor_ALT\tMajor_ALT_frq\tDepth_F\tDepth_R\n')
f1.close()

# correct genomes
for genome in assembly_quality:
    assembly_quality_genome = assembly_quality.get(genome, [])
    if assembly_quality_genome!=[]:
        database = glob.glob(os.path.join(genome_root,genome.split(fastq_name)[0]))
        if database!= []:
            database = database[0]
            correct_genome(database, assembly_quality_genome,assembly_quality_chr,report_curate)
        else:
            print('missing input fasta for %s'%(genome))
    else:
        os.system('cp %s %s.corrected.fasta'%(genome,genome))

################################################### END ########################################################
################################################### SET PATH ########################################################
# run alltestvcfstep2.sh and alltestvcfstep3.sh
# step 4 SNP filtering for genome
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics

# set up path

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_test_tree'
input_script_sub_merge = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_test_step2'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test'
genome_dir = ['/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/data']
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/'
fastq_name = 'all.SNP.fasta.genome.fasta.corrected'
ref_name = 'test.fasta'
ref_fna = 'test.fna'
genome_split = '_g'
Step = 2
deleting_file = []
Cov_dis = 20

# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def vcf_to_txt(lines,output_list,cluster_sub=[]):
    lines_set = lines.split('\n')[0].split('\t')
    if len(lines_set) >9:
        CHR = lines_set[0]
        POS = int(lines_set[1])
        temp_line = []
        temp_line.append(CHR)
        temp_line.append(str(POS))
        i = 9
        for Subdepth_all in lines_set[9:]:
            if (cluster_sub==[] and i not in deleting_set) or i in cluster_sub:
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                temp_line.append(str(total_sub_depth))
            i += 1
        output_list.append('\t'.join(temp_line)+'\n')
    else:
        print(lines)

def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = int(Subdepth[0]) + int(Subdepth[1])
    Subdepth_ALT = int(Subdepth[2]) + int(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return [allels_set[1],1]
    else:
        return [allels_set[0],0]

def outputvcf(output_name):
    vcf_file_filtered = open(vcf_file + '.%s.snp.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.samplename.txt' % (output_name), 'w')
    vcf_file_filtered.write('\t'.join(Sample_name) + '\n')
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.POS.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_POS))
    vcf_file_filtered.close()

def outputcov(output_name,vcf_file_POS_candidate,cluster_sub=[]):
    if len(vcf_file_list) > 0:
        vcf_file_POS_candidate = '\n'.join(vcf_file_POS_candidate)
        vcf_file_POS_candidate_output = ('%s' % (vcf_file_POS_candidate))
        f1 = open(os.path.join(input_script,'grep.temp.txt'),'w')
        f1.write(vcf_file_POS_candidate_output)
        f1.close()
        if '.fna.flt.snp.vcf' in vcf_file:
            cov_file = vcf_file.split('.fna.flt.snp.vcf')[0] + '.fna.raw.vcf'
        else:
            cov_file = vcf_file.split('.flt.snp.vcf')[0] + '.raw.vcf'
        os.system('grep -%s -f %s %s --no-group-separator > %s'% (
            Cov_dis, os.path.join(input_script,'grep.temp.txt'),
            cov_file,
            vcf_file + '.%s.cov.temp' % (output_name)))
        os.system('cat %s | sort | uniq > %s' % (
            vcf_file + '.%s.cov.temp' % (output_name),
            vcf_file + '.%s.uniqcov.temp' % (output_name)))
        for lines in open(vcf_file + '.%s.uniqcov.temp' % (output_name), 'r'):
            if not lines.startswith("#"):
                vcf_to_txt(lines, cov_file_list,cluster_sub)
        os.system('rm -rf %s %s %s' % (vcf_file + '.%s.cov.temp' % (output_name),
                                    vcf_file + '.%s.uniqcov.temp' % (output_name),
                                       os.path.join(input_script, 'grep.temp.txt')))
        vcf_file_filtered = open(vcf_file + '.%s.cov.txt' % (output_name), 'w')
        vcf_file_filtered.write(''.join(cov_file_list))
        vcf_file_filtered.close()

def outputtree(output_name):
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, SNP_alignment[genomename]))
            SNP_alignment_output_parsi.append('S%s    %s\n' % (genomename[-8:], SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.%s.vcf' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list_vcf))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.fasta' % (output_name), 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.parsi.fasta' %(output_name), 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()
    if Step == 2:
        if SNP_alignment_output != []:
            # with SNP
            fasta = vcf_file + '.%s.fasta' % (output_name)
            SNP_tree_cmd.append('python /scratch/users/anniz44/scripts/maffttree/remove.duplicated.py -i %s\n' % (fasta))
            filesize = int(os.path.getsize(fasta))
            if filesize >= 10 ^ 7:  # 10Mb
                SNP_tree_cmd.append(
                    'mafft --nuc --quiet --nofft --retree 2 --maxiterate 0 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                        # 'mafft --nuc --quiet --retree 2 --maxiterate 100 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                        fasta, fasta))  # remove --adjustdirection
                SNP_tree_cmd2.append(
                    '#raxmlHPC -s %s -m GTRGAMMA -n %s.raxml -p 1000 -T 40\n' % (
                        fasta, fasta))
            else:
                SNP_tree_cmd.append(
                    # 'mafft --nuc --quiet --nofft --retree 2 --maxiterate 0 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                    'mafft --nuc --quiet --retree 2 --maxiterate 100 --thread 40 %s.sorted.dereplicated.id.fasta > %s.mafft.align\n' % (
                        fasta, fasta))  # remove --adjustdirection
                SNP_tree_cmd2.append(
                    '#raxmlHPC -s %s -m GTRGAMMA -n %s.raxml -p 10000 -T 40\n' % (
                        fasta, fasta))
            SNP_tree_cmd.append('FastTree -nt -quiet %s.mafft.align > %s.mafft.align.nwk\n' % (fasta, fasta))
            SNP_tree_cmd2.append(
                '#run_gubbins.py --threads 40  --tree_builder hybrid --use_time_stamp --prefix %s_gubbins --verbose %s.mafft.align\n' % (
                    fasta, fasta))
            SNP_tree_cmd2.append(
                'mv *%s.mafft* %s\n' % (
                    fasta,output_dir_merge + '/bwa/0/tree'))
            SNP_tree_cmd2.append('rm -rf %s.sorted*\n'%(fasta))
        f1 = open(os.path.join(input_script_sub, '%s.tree.all.sh' % (donor_species)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n' + ''.join(SNP_tree_cmd) + ''.join(SNP_tree_cmd2))
        f1.close()

def SNP_seq(seq1, seq2, POS_info,POS_info_CHR,POS_info_CHR_LEN,POS_info_output,G1,G2):
    SNP_total = 0
    j = 0
    POS_DIS = []
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            CHR = POS_info_CHR[i]
            POS = POS_info[i]
            LEN = POS_info_CHR_LEN[CHR]
            if CHR == POS_info_CHR[j]:  # same CHR
                DIS = abs(POS - POS_info[j])
                POS_DIS.append(DIS)  # POS diff
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, DIS, LEN))
            else:  # new CHR
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, 0, LEN))
            j = i
    return SNP_total

def SNP_distance_correct(distanceNEW, distanceREF, SNPREF):
    return int((distanceNEW/distanceREF)*SNPREF)

def find_neighbor(Cluster_SNP,neighbor,Cluster_SNP_set,cluster,Cluster_SNP_set_added):
    if neighbor != []:
        for record_name in neighbor:
            if record_name not in Cluster_SNP_set_added:
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    subneighbor = Cluster_SNP.get(record_name,[])
                    find_neighbor(Cluster_SNP, subneighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)

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

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

def contig_to_gene(CHR, POS):
    all_genes = Mapping_loci.get(CHR,[])
    Reverse_chr = 0
    if all_genes == []:
        # database is a gene database
        codon_start = POS - 1 - int((POS - 1) % 3)
        Ref_seq_chr = Ref_seq.get(CHR, 'None')
        return [CHR, POS, codon_start, Ref_seq_chr, Reverse_chr]
    else:
        # database is a contig database
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

def SNP_check_all(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name,SNP_presence_cutoff,SNP_presence_sample_cutoff,no_SNP_cutoff,cluster_sub=[]):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    Total_qualify = 0
    Total_qualify_SNP = 0
    Total_qualify_notSNP = 0
    Total_unqualify_alt_freq = 0
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    REF_where=0
    if cluster_sub!= []:
        lines_set_sub = [lines_set[i] for i in cluster_sub]
        Total_subsample = len(cluster_sub)
        if Total_subsample >= 15:
            SNP_presence_cutoff = 0.33  # for a large group of samples
        elif Total_subsample in [3,4]:
            SNP_presence_cutoff = 1  # for a small group of samples
            SNP_presence_sample_cutoff = 2
        elif Total_subsample in [1,2]:
            SNP_presence_cutoff = 1  # for a small group of samples
            SNP_presence_sample_cutoff = 1
            no_SNP_cutoff = 0
    else:
        cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        if Step == 3:
            # re-calculate depth by a subset
            Depth4 = '%s,0,%s,0'%(sum([int(Subdepth_all.split(':')[-1].replace('\n', '').split(',')[0]) for Subdepth_all in lines_set_sub]),
                              sum([int(Subdepth_all.split(':')[-1].replace('\n', '').split(',')[1]) for Subdepth_all in
                                   lines_set_sub]))
        if Total_subsample > 2:
            REF,REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    qualify_loci = 0
                    MLF = Major_ALT[1] / total_sub_depth
                    if MLF >= Major_alt_freq_cutoff or qualify_loci == 1:
                        # major alt frequency cutoff
                        Total_qualify += 1
                        # check for qualified SNP
                        if Major_ALT[0] != REF:
                            SNP_seq[-1] = Major_ALT[0]  # unqualified SNP also include in alignment
                            # no specific rule
                            # a qualified SNP
                            temp_snp_line_pass += 'PASS'
                            Total_qualify_SNP += 1
                            SNP.add(genome_order)  # only take qualified SNP as valid SNP
                            SNP_seq[-1] = Major_ALT[0]
                        else:
                            Total_qualify_notSNP += 1
                    else:
                        # major alt frequency low
                        Total_unqualify_alt_freq += 1
            sample_num += 1
        if Total_qualify / Total_subsample >= SNP_presence_cutoff and \
                Total_unqualify_alt_freq <= Poor_MLF_freq_cutoff and\
                Total_qualify >= SNP_presence_sample_cutoff and \
                Total_qualify_SNP >= 1 and Total_qualify_SNP <= Total_qualify - no_SNP_cutoff and\
                Total_qualify_notSNP >= no_SNP_cutoff:
            # -> qualified SNP
            # qualified samples cutoff + unqualified samples cutoff
            # at least 1 qualified SNP
            # calculate NS
            gene_info = contig_to_gene(CHR, POS)
            if gene_info!= []:
                Chr_gene, POS_gene,codon_start,Ref_seq_chr,Reverse_chr  = gene_info
                if Ref_seq_chr != 'None':
                    #  observed NS ratio calculated
                    temp_snp_line_NS= [Chr_gene,str(POS_gene),'']
                    if codon_start <= POS_gene - 1:
                        Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                        Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                        SNP_seq_chr = Ref_seq_chr
                        if len(Ref_seq_codon) == 3:
                            Ref_seq_aa = translate(Ref_seq_codon)[0]
                            temp_snp_line_AA += Ref_seq_aa
                            ALT_set = allels_set
                            ALT_set.remove(REF)
                            for ALT in ALT_set:
                                SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                SNP_seq_aa = translate(SNP_seq_codon)[0]
                                temp_snp_line_AA += SNP_seq_aa
                                temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                temp_snp_line_NS[-1]+=temp_NorS
            # output lines and output major alt
            if Rough == 1: #tune cutoff
                if Total_qualify / Total_subsample >= SNP_presence_cutoff2 and \
                        Total_unqualify_alt_freq <= Poor_MLF_freq_cutoff2 and \
                        Total_qualify >= SNP_presence_sample_cutoff2 \
                        and Total_qualify_SNP >= 1 and Total_qualify_SNP <= Total_qualify - no_SNP_cutoff \
                        and Total_qualify_notSNP >= no_SNP_cutoff:
                    temp_snp_line_pass = 'PASS'
                else:
                    temp_snp_line_pass = 'NOT_PASS'
            else:
                if 'PASS' in temp_snp_line_pass:
                    temp_snp_line_pass = 'PASS'
                else:
                    temp_snp_line_pass = 'RULE'
            if CHR == CHR_old:
                # same CHR
                POS_DIS = abs(POS - POS_old)
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
            else:
                # diff CHR first SNP
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
            POS_old = POS
            CHR_old = CHR
            temp_snp_line.append(CHR)
            temp_snp_line.append(str(POS))
            temp_snp_line.append(REF)
            temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
            vcf_file_list.append('\t'.join(temp_snp_line)+ '\t' +'\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,'\t'.join(temp_snp_line_NS),temp_snp_line_AA))
            vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
            vcf_file_list_vcf.append('\t'.join(lines_set[0:9])+'\t'+'\t'.join(lines_set_sub)+'\n')
            if Rough != 1:
                i = 9
                j = 0
                SNP_alignment[reference_name] += REF
                for genomename in SNP_alignment:
                    if genomename != reference_name:
                        if i in cluster_sub:
                            SNP_alignment[genomename] += SNP_seq[j]
                            j += 1
                        i += 1
    return [CHR_old,POS_old]

# set up output

try:
    os.mkdir(output_dir_merge + '/bwa/0/tree')
except IOError:
    pass

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# set up steps
SNP_cluster = dict()
cluster_set = set()
if Step == 1:
    reference_set = ['']
    outputname_set = ['filteredrough']
elif Step == 2:
    reference_set = ['reference']
    outputname_set = ['filtered']
elif Step == 3:
    reference_set = []
    outputname_set = []
    # load cluster
    for lines in open(os.path.join(input_script, 'SNP_round2.sum'), 'r'):
        if not lines.startswith('donor_species'):
            lines_set = lines.split('\n')[0].split('\t')
            genomes = lines_set[1]
            cluster_type = lines_set[-3]
            cluster_total = int(lines_set[-2])
            if cluster_total > 1:
                SNP_cluster.setdefault(genomes, cluster_type)
                cluster_set.add(cluster_type)
    for cluster_type in cluster_set:
        reference_set.append('refer_%s'%(cluster_type))
        outputname_set.append('filtered.%s' % (cluster_type))
    cluster_set = list(cluster_set)

# set up cutoff super rough
if Step ==1 :
    Rough = 1  # tune cutoff
    SNP_presence_cutoff = 0.33 # avg presence in all samples
else:
    Rough = 0
    SNP_presence_cutoff = 0.66  # avg presence in all samples

# unchanged cutoff
SNP_presence_sample_cutoff = 3  # num of samples passing the above criteria
Major_alt_freq_cutoff = 0.9 # major alt freq in a genome, do not allow multiple homolougs genes
no_SNP_cutoff = 1
Poor_MLF_freq_cutoff = 1 # no sample should have homologous genes (low major alt freq)
# set up strict cutoff
SNP_presence_cutoff2 = 0.66 # avg coverage in all samples
SNP_presence_sample_cutoff2 = 3  # num of samples passing the above criteria
Poor_MLF_freq_cutoff2 = 1 # no sample should have homologous genes (low major alt freq)

# cluster cutoff Step 3
SNP_total_cutoff_2 = 50
cluster_cutoff = 2

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir_merge + '/bwa/0','*%s*.flt.snp.vcf'%(fastq_name)))
for set_num in range(0,len(reference_set)):
    reference_name = reference_set[set_num]
    output_name = outputname_set[set_num]
    for vcf_file in all_vcf_file:
        SNP_presence_cutoff = SNP_presence_cutoff2 # for group of samples  -> used
        SNP_presence_sample_cutoff = SNP_presence_sample_cutoff2
        no_SNP_cutoff = 1
        print(vcf_file)
        Total = 0
        donor_species = os.path.split(vcf_file)[-1].split('.fna.flt.snp.vcf')[0].split('.flt.snp.vcf')[0]
        Sample_name = []
        deleting_set = []
        Ref_seq = dict()
        Mapping = dict()
        Mapping_loci = dict()
        SNP_cluster_donor_species = dict()
        for cluster_type in cluster_set:
            SNP_cluster_donor_species.setdefault(cluster_type,[])
        for lines in open(vcf_file.replace('.flt.snp.vcf','.raw.vcf'), 'r'):
            if lines.startswith('##bcftoolsCommand=mpileup '):
                # setup samples
                sample_set = lines.split(ref_name + ' ')[1].split('\n')[0].split(' |')[0].split(' ')
                samplenum = 9
                for samples in sample_set[1:]:
                    genomename = os.path.split(samples)[-1].split(fastq_name)[0]
                    Sample_name.append(genomename.replace('.', ''))
                    if genomename in deleting_file:
                        deleting_set.append(samplenum)
                    else:
                        if SNP_cluster!= dict() and genomename in SNP_cluster:
                            SNP_cluster_donor_species[SNP_cluster[genomename]].append(samplenum)
                    samplenum += 1
                break
        # subset samples
        Sample_subset = []
        if Step == 3:
            Sample_subset = SNP_cluster_donor_species.get(cluster_set[set_num], [])
            Sample_name = [Sample_name[i-9] for i in Sample_subset]
        print(Sample_name)
        if Step < 3 or (Sample_subset != [] and len(Sample_subset)>= cluster_cutoff):
            print('running %s' %(donor_species))
            # load database
            database_file = glob.glob(os.path.join(genome_root,
                                                   '%s' % (ref_fna)))
            if database_file != []:
                Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file[0])
            vcf_file_all = [vcf_file,vcf_file.replace('.fna.flt.snp.vcf','.flt.snp.vcf')]
            for vcf_file in vcf_file_all:
                print(vcf_file)
                SNP_tree_cmd = []
                SNP_tree_cmd2 = []
                vcf_file_list = []
                vcf_file_list_vcf = []
                vcf_file_POS = []
                vcf_file_POS_candidate = set()
                SNP_alignment = dict()
                SNP_alignment.setdefault(reference_name, '')
                cov_file_list = []
                CHR_old = ''
                POS_old = 0
                for genomename in Sample_name:
                    SNP_alignment.setdefault(genomename, '')
                for lines in open(vcf_file, 'r'):
                    if not lines.startswith("#"):
                        lines_set = lines.split('\n')[0].split('\t')
                        CHR = lines_set[0]
                        POS = int(lines_set[1])
                        Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                        if Total == 0:
                            Total = len(lines_set) - 9 - len(deleting_set)
                            if Total >= 15:
                                SNP_presence_cutoff = 0.33  # for a large group of genomes
                            elif Total in [3,4]:
                                SNP_presence_cutoff = 1  # for a small group of genomes
                                SNP_presence_sample_cutoff = 2
                            elif Total in [1,2]:
                                SNP_presence_cutoff = 1  # for only 1 or 2 samples, compare to ref
                                SNP_presence_sample_cutoff = 1
                                no_SNP_cutoff = 0
                        if Depth / Total >= SNP_presence_cutoff:
                            # average depth in all samples cutoff
                            if "INDEL" not in lines_set[7] \
                                    and (lines_set[6] != 'LowQual'):
                                if Step == 1:
                                    CHR_old, POS_old =  SNP_check_all(lines_set, '',
                                                                      CHR_old,POS_old,reference_name,
                                                                      SNP_presence_cutoff,SNP_presence_sample_cutoff,no_SNP_cutoff)
                                elif Step == 2:
                                    CHR_old, POS_old = SNP_check_all(lines_set, '',
                                                                     CHR_old,POS_old,reference_name,
                                                                     SNP_presence_cutoff,SNP_presence_sample_cutoff,no_SNP_cutoff)
                                else:
                                    CHR_old, POS_old = SNP_check_all(lines_set, '',
                                                                     CHR_old, POS_old, reference_name,
                                                                     SNP_presence_cutoff, SNP_presence_sample_cutoff,no_SNP_cutoff,
                                                                     Sample_subset)
                outputvcf(output_name)
                if Step != 1:
                    # step 2 and 3 output fasta
                    outputtree(output_name)
                if Step != 2:
                    # step 1 and 3 output coverage, step 2 use step 1 coverage
                    try:
                        f1 = open(vcf_file + '.%s.cov.txt' % (output_name), 'r')
                    except IOError:
                        outputcov(output_name,list(vcf_file_POS_candidate),Sample_subset)

################################################### END ########################################################
################################################### SET PATH ########################################################
# step 5 filter SNPs for WGS
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics

# set up path

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_test_tree'
input_script_sub_merge = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_test'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test'
genome_dir = ['/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/data']
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/'
fastq_name = 'all.SNP.fasta1.fq'
genome_split = '_g'
ref_name = 'test.fasta'
ref_fna = 'test.fna'
deleting_file = []
Step = 1
#Step = 2 # sometimes need to change MLF when model add not enough reads
Cov_dis = 20
rule1=[]
rule2=[]
rule3 = []

################################################### Function ########################################################
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def vcf_to_txt(lines,output_list,cluster_sub=[]):
    lines_set = lines.split('\n')[0].split('\t')
    if len(lines_set) >9:
        CHR = lines_set[0]
        POS = int(lines_set[1])
        temp_line = []
        temp_line.append(CHR)
        temp_line.append(str(POS))
        i = 9
        for Subdepth_all in lines_set[9:]:
            if (cluster_sub==[] and i not in deleting_set) or i in cluster_sub:
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                temp_line.append(str(total_sub_depth))
            i += 1
        output_list.append('\t'.join(temp_line)+'\n')
    else:
        print(lines)

def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = int(Subdepth[0]) + int(Subdepth[1])
    Subdepth_ALT = int(Subdepth[2]) + int(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return [allels_set[1],1]
    else:
        return [allels_set[0],0]

def outputvcf(output_name):
    vcf_file_filtered = open(vcf_file + '.%s.snp.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.samplename.txt' % (output_name), 'w')
    vcf_file_filtered.write('\t'.join(Sample_name) + '\n')
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.POS.txt' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_POS))
    vcf_file_filtered.close()

def outputcov(output_name,vcf_file_POS_candidate,cluster_sub=[]):
    if len(vcf_file_list) > 0:
        vcf_file_POS_candidate = '\n'.join(vcf_file_POS_candidate)
        vcf_file_POS_candidate_output = ('%s' % (vcf_file_POS_candidate))
        f1 = open(os.path.join(input_script,'grep.temp.txt'),'w')
        f1.write(vcf_file_POS_candidate_output)
        f1.close()
        os.system('grep -%s -f %s %s --no-group-separator > %s'% (
            Cov_dis, os.path.join(input_script,'grep.temp.txt'),
            vcf_file.split('.flt.snp.vcf')[0] + '.raw.vcf',
            vcf_file + '.%s.cov.temp' % (output_name)))
        os.system('cat %s | sort | uniq > %s' % (
            vcf_file + '.%s.cov.temp' % (output_name),
            vcf_file + '.%s.uniqcov.temp' % (output_name)))
        for lines in open(vcf_file + '.%s.uniqcov.temp' % (output_name), 'r'):
            if not lines.startswith("#"):
                vcf_to_txt(lines, cov_file_list,cluster_sub)
        os.system('rm -rf %s %s' % (vcf_file + '.%s.cov.temp' % (output_name),
                                    vcf_file + '.%s.uniqcov.temp' % (output_name)))
        vcf_file_filtered = open(vcf_file + '.%s.cov.txt' % (output_name), 'w')
        vcf_file_filtered.write(''.join(cov_file_list))
        vcf_file_filtered.close()

def outputtree(output_name):
    SNP_alignment_output = []
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        if seq_len > 0:
            SNP_alignment_output.append('>%s\n%s\n' % (genomename, SNP_alignment[genomename]))
            if '_cluster' in genomename:
                newgenomename = '%s_CL%s'%(genomename.split('_')[0],genomename.split('_cluster')[1].split('.')[0])
            else:
                newgenomename = '%s_CL0' % (genomename.split('_')[0])
            if '.cluster' in genomename:
                newgenomename = '%s.CL%s'%(newgenomename,genomename.split('.cluster')[1].split('.')[0])
            else:
                newgenomename = '%s.CL0' % (newgenomename)
            if 'reference' in genomename:
                newgenomename = genomename
            SNP_alignment_output_parsi.append('%s    %s\n' % (newgenomename, SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(vcf_file + '.%s.vcf' % (output_name), 'w')
    vcf_file_filtered.write(''.join(vcf_file_list_vcf))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.fasta' % (output_name), 'w')
    vcf_file_filtered.write(''.join(SNP_alignment_output))
    vcf_file_filtered.close()
    vcf_file_filtered = open(vcf_file + '.%s.parsi.fasta' %(output_name), 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()

def SNP_seq(seq1, seq2, POS_info,POS_info_CHR,POS_info_CHR_LEN,POS_info_output,G1,G2):
    SNP_total = 0
    j = 0
    POS_DIS = []
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            CHR = POS_info_CHR[i]
            POS = POS_info[i]
            LEN = POS_info_CHR_LEN[CHR]
            if CHR == POS_info_CHR[j]:  # same CHR
                DIS = abs(POS - POS_info[j])
                POS_DIS.append(DIS)  # POS diff
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, DIS, LEN))
            else:  # new CHR
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, 0, LEN))
            j = i
    return SNP_total

def SNP_distance_correct(distanceNEW, distanceREF, SNPREF):
    return int((distanceNEW/distanceREF)*SNPREF)

def find_neighbor(Cluster_SNP,neighbor,Cluster_SNP_set,cluster,Cluster_SNP_set_added):
    if neighbor != []:
        for record_name in neighbor:
            if record_name not in Cluster_SNP_set_added:
                    Cluster_SNP_set[cluster].add(record_name)
                    Cluster_SNP_set_added.add(record_name)
                    subneighbor = Cluster_SNP.get(record_name,[])
                    find_neighbor(Cluster_SNP, subneighbor, Cluster_SNP_set, cluster,Cluster_SNP_set_added)

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

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

def contig_to_gene(CHR, POS):
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

def SNP_check_all_singleend(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name,SNP_presence_cutoff,SNP_presence_sample_cutoff,no_SNP_cutoff,cluster_sub=[]):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    Total_qualify = 0
    Total_qualify2 = 0
    Total_qualify_SNP = 0
    Total_qualify_SNP2 = 0
    Total_qualify_notSNP = 0
    Total_qualify_notSNP2 = 0
    Total_unqualify_alt_freq = 0
    Total_unqualify_alt_freq2 = 0
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    REF_where= 0
    if cluster_sub!= []:
        lines_set_sub = [lines_set[i] for i in cluster_sub]
        Total_subsample = len(cluster_sub)
        if Total_subsample >= 15:
            SNP_presence_cutoff = 0.33  # for a large group of samples
        elif Total_subsample in [3,4]:
            SNP_presence_cutoff = 1  # for a small group of samples
            SNP_presence_sample_cutoff = 2
        elif Total_subsample in [1,2]:
            SNP_presence_cutoff = 1  # for a small group of samples
            SNP_presence_sample_cutoff = 1
            no_SNP_cutoff = 0
    else:
        cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        if Step == 3:
            # re-calculate depth by a subset
            Depth4 = '%s,0,%s,0'%(sum([int(Subdepth_all.split(':')[-1].replace('\n', '').split(',')[0]) for Subdepth_all in lines_set_sub]),
                              sum([int(Subdepth_all.split(':')[-1].replace('\n', '').split(',')[1]) for Subdepth_all in
                                   lines_set_sub]))
        if Total_subsample > 2:
            REF,REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
                Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
                total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
                total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                        Allels_frq_sub[Allels[allels] * 2] += int(Subdepth_forward[num_allels])
                        Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_reverse[num_allels])
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    qualify_SNP = 0  # whether a qualified SNP
                    qualify_loci = 0
                    MLF = Major_ALT[1] / total_sub_depth
                    if total_sub_depth_forward >= Sample_depth_cutoff or \
                            total_sub_depth_reverse >= Sample_depth_cutoff:
                        # forward and reverse cutoff POS detected
                        if MLF >= Major_alt_freq_cutoff or qualify_loci == 1:
                        # major alt frequency cutoff
                            Total_qualify += 1
                            # check for qualified SNP
                            if Major_ALT[0] != REF:
                                SNP_seq[-1] = Major_ALT[0]  # unqualified SNP also include in alignment
                                # no specific rule
                                if (total_sub_depth_forward - int(Subdepth_forward[REF_where]) >= SNP_depth_cutoff or \
                                        total_sub_depth_reverse - int(Subdepth_reverse[REF_where]) >= SNP_depth_cutoff) and \
                                        Major_ALT[1] >= 2 * SNP_depth_cutoff:
                                    qualify_SNP = 1
                                    temp_snp_line_pass += 'PASS'
                            else:
                                Total_qualify_notSNP += 1
                            # tune cutoff
                            if Rough == 1:
                                if total_sub_depth_forward >= Sample_depth_cutoff2 and \
                                        total_sub_depth_reverse >= Sample_depth_cutoff2:
                                    # forward and reverse cutoff
                                    if MLF >= Major_alt_freq_cutoff2:
                                        # major alt frequency cutoff
                                        Total_qualify2 += 1
                                        if Major_ALT[0] != REF:
                                            if total_sub_depth_forward - int(Subdepth_forward[REF_where]) >= SNP_depth_cutoff2 and \
                                                    total_sub_depth_reverse - int(Subdepth_reverse[REF_where]) >= SNP_depth_cutoff2 and \
                                                    Major_ALT[1] >= 2 * SNP_depth_cutoff2:
                                                Total_qualify_SNP2 += 1
                                        else:
                                            Total_qualify_notSNP2 += 1
                                    else:
                                        Total_unqualify_alt_freq2 += 1
                        else:
                            # major alt frequency low
                            Total_unqualify_alt_freq += 1
                        # a qualified SNP
                        if qualify_SNP == 1:
                            Total_qualify_SNP += 1
                            SNP.add(genome_order)  # only take qualified SNP as valid SNP
                            SNP_seq[-1] = Major_ALT[0]
            sample_num += 1
        if Total_qualify / Total_subsample >= SNP_presence_cutoff and \
            Total_unqualify_alt_freq / Total_subsample <= Poor_MLF_freq_cutoff and\
            Total_qualify >= SNP_presence_sample_cutoff and \
            Total_qualify_SNP >= 1 and Total_qualify_SNP <= Total_qualify - no_SNP_cutoff and\
            Total_qualify_notSNP >= no_SNP_cutoff:
            # -> qualified SNP
            # qualified samples cutoff + unqualified samples cutoff
            # at least 1 qualified SNP
            # calculate NS
            gene_info = contig_to_gene(CHR, POS)
            if gene_info!= []:
                Chr_gene, POS_gene,codon_start,Ref_seq_chr,Reverse_chr  = gene_info
                if Ref_seq_chr != 'None':
                    #  observed NS ratio calculated
                    temp_snp_line_NS= [Chr_gene,str(POS_gene),'']
                    if codon_start <= POS_gene - 1:
                        Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                        Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                        SNP_seq_chr = Ref_seq_chr
                        if len(Ref_seq_codon) == 3:
                            Ref_seq_aa = translate(Ref_seq_codon)[0]
                            temp_snp_line_AA += Ref_seq_aa
                            ALT_set = allels_set
                            ALT_set.remove(REF)
                            for ALT in ALT_set:
                                SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                SNP_seq_aa = translate(SNP_seq_codon)[0]
                                temp_snp_line_AA += SNP_seq_aa
                                temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                temp_snp_line_NS[-1]+=temp_NorS
            # output lines and output major alt
            if Rough == 1: #tune cutoff
                pass
            else:
                if 'PASS' in temp_snp_line_pass:
                    temp_snp_line_pass = 'PASS'
                else:
                    temp_snp_line_pass = 'RULE'
            if CHR == CHR_old:
                # same CHR
                POS_DIS = abs(POS - POS_old)
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
            else:
                # diff CHR first SNP
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
            POS_old = POS
            CHR_old = CHR
            temp_snp_line.append(CHR)
            temp_snp_line.append(str(POS))
            temp_snp_line.append(REF)
            temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
            vcf_file_list.append('\t'.join(temp_snp_line)+ '\t' +'\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,'\t'.join(temp_snp_line_NS),temp_snp_line_AA))
            vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
            vcf_file_list_vcf.append('\t'.join(lines_set[0:9])+'\t'+'\t'.join(lines_set_sub)+'\n')
            if Rough != 1:
                i = 9
                j = 0
                SNP_alignment[reference_name] += REF
                for genomename in SNP_alignment:
                    if genomename != reference_name:
                        if i in cluster_sub:
                            SNP_alignment[genomename] += SNP_seq[j]
                            j += 1
                        i += 1
    return [CHR_old,POS_old]

def SNP_check_all(lines_set,temp_snp_line_pass,CHR_old,POS_old,reference_name,SNP_presence_cutoff,SNP_presence_sample_cutoff,no_SNP_cutoff,cluster_sub=[]):
    temp_snp_line = []
    temp_snp_line_frq = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    Total_qualify = 0
    Total_qualify2 = 0
    Total_qualify_SNP = 0
    Total_qualify_SNP2 = 0
    Total_qualify_notSNP = 0
    Total_qualify_notSNP2 = 0
    Total_unqualify_alt_freq = 0
    Total_unqualify_alt_freq2 = 0
    SNP = set()
    SNP_seq = []
    REF = lines_set[3]
    allels_set = [REF]
    Total_subsample = Total
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    REF_where= 0
    if cluster_sub!= []:
        lines_set_sub = [lines_set[i] for i in cluster_sub]
        Total_subsample = len(cluster_sub)
        if Total_subsample >= 15:
            SNP_presence_cutoff = 0.33  # for a large group of samples
        elif Total_subsample in [3,4]:
            SNP_presence_cutoff = 1  # for a small group of samples
            SNP_presence_sample_cutoff = 2
        elif Total_subsample in [1,2]:
            SNP_presence_cutoff = 1  # for a small group of samples
            SNP_presence_sample_cutoff = 1
            no_SNP_cutoff = 0
    else:
        cluster_sub = list(range(9, len(lines_set)))
    if Total_subsample > 0:
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        genome_order = 0
        Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
        if Step == 3:
            # re-calculate depth by a subset
            Depth4 = '%s,0,%s,0'%(sum([int(Subdepth_all.split(':')[-1].replace('\n', '').split(',')[0]) for Subdepth_all in lines_set_sub]),
                              sum([int(Subdepth_all.split(':')[-1].replace('\n', '').split(',')[1]) for Subdepth_all in
                                   lines_set_sub]))
        if Total_subsample > 2:
            REF,REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
        sample_num = 9
        for Subdepth_all in lines_set_sub:
            if sample_num not in deleting_set:
                genome_order += 1
                Allels_frq = [0, 0, 0, 0]
                Allels_frq_sub = [0, 0, 0, 0, 0, 0, 0, 0]
                Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
                Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
                Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
                total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
                total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
                for num_allels in range(0, Total_alleles):
                    allels = allels_set[num_allels]
                    Subdepth_alleles = int(Subdepth[num_allels])
                    if allels in Allels:
                        Allels_frq[Allels[allels]] += Subdepth_alleles
                        Allels_frq_sub[Allels[allels] * 2] += int(Subdepth_forward[num_allels])
                        Allels_frq_sub[Allels[allels] * 2 + 1] += int(Subdepth_reverse[num_allels])
                    else:
                        pass
                # find major alt and calculate frequency
                Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
                temp_snp_line_frq.append(';'.join(str(frq_sub) for frq_sub in Allels_frq_sub))
                SNP_seq.append(REF)  # set as reference
                if total_sub_depth > 0:
                    qualify_SNP = 0  # whether a qualified SNP
                    qualify_loci = 0
                    MLF = Major_ALT[1] / total_sub_depth
                    if total_sub_depth_forward >= Sample_depth_cutoff or \
                            total_sub_depth_reverse >= Sample_depth_cutoff:
                        # forward and reverse cutoff POS detected
                        if MLF >= Major_alt_freq_cutoff or qualify_loci == 1:
                        # major alt frequency cutoff
                            Total_qualify += 1
                            # check for qualified SNP
                            if Major_ALT[0] != REF:
                                SNP_seq[-1] = Major_ALT[0]  # unqualified SNP also include in alignment
                                # no specific rule
                                if (total_sub_depth_forward - int(Subdepth_forward[REF_where]) >= SNP_depth_cutoff or \
                                        total_sub_depth_reverse - int(Subdepth_reverse[REF_where]) >= SNP_depth_cutoff) and \
                                        Major_ALT[1] >= 2 * SNP_depth_cutoff:
                                    qualify_SNP = 1
                                    temp_snp_line_pass += 'PASS'
                            else:
                                Total_qualify_notSNP += 1
                            # tune cutoff
                            if Rough == 1:
                                if total_sub_depth_forward >= Sample_depth_cutoff2 and \
                                        total_sub_depth_reverse >= Sample_depth_cutoff2:
                                    # forward and reverse cutoff
                                    if MLF >= Major_alt_freq_cutoff2:
                                        # major alt frequency cutoff
                                        Total_qualify2 += 1
                                        if Major_ALT[0] != REF:
                                            if total_sub_depth_forward - int(Subdepth_forward[REF_where]) >= SNP_depth_cutoff2 and \
                                                    total_sub_depth_reverse - int(Subdepth_reverse[REF_where]) >= SNP_depth_cutoff2 and \
                                                    Major_ALT[1] >= 2 * SNP_depth_cutoff2:
                                                Total_qualify_SNP2 += 1
                                        else:
                                            Total_qualify_notSNP2 += 1
                                    else:
                                        Total_unqualify_alt_freq2 += 1
                        else:
                            # major alt frequency low
                            Total_unqualify_alt_freq += 1
                        # a qualified SNP
                        if qualify_SNP == 1:
                            Total_qualify_SNP += 1
                            SNP.add(genome_order)  # only take qualified SNP as valid SNP
                            SNP_seq[-1] = Major_ALT[0]
            sample_num += 1
        if Total_qualify / Total_subsample >= SNP_presence_cutoff and \
            Total_unqualify_alt_freq / Total_subsample <= Poor_MLF_freq_cutoff and\
            Total_qualify >= SNP_presence_sample_cutoff and \
            Total_qualify_SNP >= 1 and Total_qualify_SNP <= Total_qualify - no_SNP_cutoff and\
            Total_qualify_notSNP >= no_SNP_cutoff:
            # -> qualified SNP
            # qualified samples cutoff + unqualified samples cutoff
            # at least 1 qualified SNP
            # calculate NS
            gene_info = contig_to_gene(CHR, POS)
            if gene_info!= []:
                Chr_gene, POS_gene,codon_start,Ref_seq_chr,Reverse_chr  = gene_info
                if Ref_seq_chr != 'None':
                    #  observed NS ratio calculated
                    temp_snp_line_NS= [Chr_gene,str(POS_gene),'']
                    if codon_start <= POS_gene - 1:
                        Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                        Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                        SNP_seq_chr = Ref_seq_chr
                        if len(Ref_seq_codon) == 3:
                            Ref_seq_aa = translate(Ref_seq_codon)[0]
                            temp_snp_line_AA += Ref_seq_aa
                            ALT_set = allels_set
                            ALT_set.remove(REF)
                            for ALT in ALT_set:
                                SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                                SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                                SNP_seq_aa = translate(SNP_seq_codon)[0]
                                temp_snp_line_AA += SNP_seq_aa
                                temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                                temp_snp_line_NS[-1]+=temp_NorS
            # output lines and output major alt
            if Rough == 1: #tune cutoff
                pass
            else:
                if 'PASS' in temp_snp_line_pass:
                    temp_snp_line_pass = 'PASS'
                else:
                    temp_snp_line_pass = 'RULE'
            if CHR == CHR_old:
                # same CHR
                POS_DIS = abs(POS - POS_old)
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, POS_DIS))
            else:
                # diff CHR first SNP
                vcf_file_POS.append('%s\t%s\t%s\n' % (CHR, POS, 0))
            POS_old = POS
            CHR_old = CHR
            temp_snp_line.append(CHR)
            temp_snp_line.append(str(POS))
            temp_snp_line.append(REF)
            temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
            vcf_file_list.append('\t'.join(temp_snp_line)+ '\t' +'\t'.join(temp_snp_line_frq) + '\t\"%s\"\t%s\t%s\t%s\n' % (
                ';'.join(str(genome_order) for genome_order in SNP), temp_snp_line_pass,'\t'.join(temp_snp_line_NS),temp_snp_line_AA))
            vcf_file_POS_candidate.add('%s\t%s\t' % (CHR, POS))
            vcf_file_list_vcf.append('\t'.join(lines_set[0:9])+'\t'+'\t'.join(lines_set_sub)+'\n')
            if Rough != 1:
                i = 9
                j = 0
                SNP_alignment[reference_name] += REF
                for genomename in SNP_alignment:
                    if genomename != reference_name:
                        if i in cluster_sub:
                            SNP_alignment[genomename] += SNP_seq[j]
                            j += 1
                        i += 1
    return [CHR_old,POS_old]

################################################### Set up ########################################################
# set up output
try:
    os.mkdir(output_dir_merge + '/bwa/0/tree')
except IOError:
    pass

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# set up steps
SNP_cluster = dict()
cluster_set = set()
if Step == 1:
    reference_set = ['']
    outputname_set = ['filteredrough']
elif Step == 2:
    reference_set = ['reference']
    outputname_set = ['filtered']
elif Step == 3:
    reference_set = []
    outputname_set = []
    # load cluster
    for lines in open(os.path.join(input_script, 'SNP_round2.sum'), 'r'):
        if not lines.startswith('donor_species'):
            lines_set = lines.split('\n')[0].split('\t')
            genomes = lines_set[1]
            cluster_type = lines_set[-3]
            cluster_total = int(lines_set[-2])
            if cluster_total > 1:
                SNP_cluster.setdefault(genomes, cluster_type)
                cluster_set.add(cluster_type)
    for cluster_type in cluster_set:
        reference_set.append('refer_%s'%(cluster_type))
        outputname_set.append('filtered.%s' % (cluster_type))
    cluster_set = list(cluster_set)

# set up cutoff super rough
if Step == 1:
    Rough = 1  # tune cutoff
    median_coverage_cutoff = 4 # for group of samples -> used 3*2*0.66
    Major_alt_freq_cutoff = 0.7 # for SNP samples -> used
    SNP_depth_cutoff = 3 # for a SNP sample -> used
else:
    Rough = 0
    median_coverage_cutoff = 10  # avg coverage in all samples
    Major_alt_freq_cutoff = 0.9  # major alt freq in a sample
    SNP_depth_cutoff = 5  # both forward and reverse reads cutoff for SNP ALTs in a sample

# unchanged cutoff
Sample_depth_cutoff = 3  # both forward and reverse reads cutoff in a sample
SNP_presence_cutoff = 0.66  # percentage of samples passing the above criteria
SNP_presence_sample_cutoff = 3  # num of samples passing the above criteria
no_SNP_cutoff = 1
# set up strict cutoff
median_coverage_cutoff2 = 10 # avg coverage in all samples
Sample_depth_cutoff2 = 3 # both forward and reverse reads cutoff in a sample
Major_alt_freq_cutoff2 = 0.9 # major alt freq in a sample
SNP_presence_cutoff2 = 0.66 # percentage of samples passing the above criteria
SNP_presence_sample_cutoff2 = 3 # num of samples passing the above criteria
SNP_depth_cutoff2 = 5 # both forward and reverse reads cutoff for SNP ALTs in a sample
Poor_MLF_freq_cutoff = (1 - SNP_presence_cutoff2)*0.75 #the unqualifie samples should be mostly low cov but not two alleles (low major alt freq)

# cluster cutoff Step 3
SNP_total_cutoff_2 = 50
cluster_cutoff = 2
# coverage cutoff Step 1
total_coverage_cutoff = 0.7
genome_avg_coverage_cutoff = 10

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

################################################### Main ########################################################
# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir_merge + '/bwa/0','*%s*.flt.snp.vcf'%(fastq_name)))
for set_num in range(0,len(reference_set)):
    reference_name = reference_set[set_num]
    output_name = outputname_set[set_num]
    for vcf_file in all_vcf_file:
        SNP_presence_cutoff = SNP_presence_cutoff2 # for group of samples  -> used
        SNP_presence_sample_cutoff = SNP_presence_sample_cutoff2
        Poor_MLF_freq_cutoff = (1 - SNP_presence_cutoff) * 0.75
        no_SNP_cutoff = 1
        print(vcf_file)
        Total = 0
        donor_species = os.path.split(vcf_file)[-1].split('.flt.snp.vcf')[0]
        SNP_tree_cmd = []
        SNP_tree_cmd2 = []
        vcf_file_list = []
        vcf_file_list_vcf = []
        Sample_name = []
        deleting_set = []
        vcf_file_POS = []
        vcf_file_POS_candidate = set()
        SNP_alignment = dict()
        SNP_alignment.setdefault(reference_name, '')
        cov_file_list = []
        Ref_seq = dict()
        Mapping = dict()
        Mapping_loci = dict()
        CHR_old = ''
        POS_old = 0
        SNP_cluster_donor_species = dict()
        for cluster_type in cluster_set:
            SNP_cluster_donor_species.setdefault(cluster_type,[])
        for lines in open(vcf_file.replace('.flt.snp.vcf','.raw.vcf'), 'r'):
            if lines.startswith('##bcftoolsCommand=mpileup '):
                # setup samples
                sample_set = lines.split(ref_name + ' ')[1].split('\n')[0].split('  |')[0].split(' ')
                samplenum = 9
                for samples in sample_set:
                    genomename = os.path.split(samples)[-1].split(fastq_name + '.sorted.bam')[0]
                    Sample_name.append(genomename.replace('.', ''))
                    if genomename in deleting_file:
                        deleting_set.append(samplenum)
                    else:
                        SNP_alignment.setdefault(genomename, '')
                        if SNP_cluster!= dict() and genomename in SNP_cluster:
                            SNP_cluster_donor_species[SNP_cluster[genomename]].append(samplenum)
                    samplenum += 1
        # subset samples
        Sample_subset = []
        if Step == 3:
            Sample_subset = SNP_cluster_donor_species.get(cluster_set[set_num], [])
            Sample_name = [Sample_name[i-9] for i in Sample_subset]
        if Step < 3 or (Sample_subset != [] and len(Sample_subset)>= cluster_cutoff):
            print('running %s' %(donor_species))
            # load database
            database_file = glob.glob(os.path.join(genome_root,
                                                   '%s' % (ref_fna)))
            if database_file != []:
                Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file[0])
            for lines in open(vcf_file, 'r'):
                if not lines.startswith("#"):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR = lines_set[0]
                    POS = int(lines_set[1])
                    Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                    if Total == 0:
                        Total = len(lines_set) - 9 - len(deleting_set)
                        if Total >= 15:
                            SNP_presence_cutoff = 0.33  # for a large group of genomes
                            Poor_MLF_freq_cutoff = 0.3
                        elif Total in [3,4]:
                            SNP_presence_cutoff = 1  # for a small group of genomes
                            SNP_presence_sample_cutoff = 2
                        elif Total in [1,2]:
                            SNP_presence_cutoff = 1
                            SNP_presence_sample_cutoff = 1
                            no_SNP_cutoff = 0
                    if Depth / Total >= median_coverage_cutoff:
                        # average depth in all samples cutoff
                        if "INDEL" not in lines_set[7] \
                                and (lines_set[6] != 'LowQual'):
                            if Step == 1:
                                CHR_old, POS_old =  SNP_check_all_singleend(lines_set, '',
                                                                  CHR_old,POS_old,reference_name,
                                                                  SNP_presence_cutoff,SNP_presence_sample_cutoff,no_SNP_cutoff)
                            elif Step == 2:
                                CHR_old, POS_old = SNP_check_all_singleend(lines_set, '',
                                                                 CHR_old,POS_old,reference_name,
                                                                 SNP_presence_cutoff,SNP_presence_sample_cutoff,no_SNP_cutoff)
                            else:
                                CHR_old, POS_old = SNP_check_all_singleend(lines_set, '',
                                                                 CHR_old, POS_old, reference_name,
                                                                 SNP_presence_cutoff, SNP_presence_sample_cutoff,no_SNP_cutoff,
                                                                 Sample_subset)
            outputvcf(output_name)
            outputtree(output_name)
            if Step != 2:
                # output coverage
                try:
                    f1 = open(vcf_file + '.%s.cov.txt' % (output_name), 'r')
                except IOError:
                    outputcov(output_name, list(vcf_file_POS_candidate), Sample_subset)

################################################### END ########################################################
################################################### SET PATH ########################################################
# compare genome call SNPs VS WGS call SNPs
# Question: why genome mapping shows SNPs that have no alternative alleles in raw.vcf of WGS mapping?
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics

# set up path
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/bwa/0/'
ref_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/'
# ref
vcf_ref = glob.glob(ref_dir + '/test.*.SNP.fasta.snp.txt')
# set up cutoff
end_cutoff = 70 # 10 bp at the ends of a contig, separate

# function
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
                vcf_count[contig_end('%s\t%s\t\n'%(CHR,POS))] += 1
    if end_count:
        summary_file_output.append('%s\t%s\t%s\t0\t0\n'%(vcf_file,vcf_count['end'],
                                                               vcf_count['notend']))
        print(vcf_count)
    return vcf_input

def contig_end(CHRPOS):
    CHR, POS = CHRPOS.split('\t')[0:2]
    total_length = CHR.split('size')[1]
    total_length = int(total_length)
    if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
        return 'end'
    else:
        return 'notend'

def compare_vcf(vcf_input1, vcf_input2, vcf_file1, vcf_file2, output_file):
    end_set = ['end','notend']
    vcf_1_diff = dict()
    vcf_1_diff.setdefault('end',[])
    vcf_1_diff.setdefault('notend', [])
    # CHRPOS in 1 not in 2
    for CHRPOS in vcf_input1:
        if CHRPOS not in vcf_input2:
            contig_end_tag = contig_end(CHRPOS)
            vcf_1_diff[contig_end_tag].append(CHRPOS)
    for contig_end_tag in end_set:
        if len(vcf_1_diff[contig_end_tag]) > 0:
            print(len(vcf_1_diff[contig_end_tag]))
            temp_output = os.path.join(output_dir, 'grep.temp.%s.txt' % (contig_end_tag))
            print(temp_output)
            f1 = open(temp_output, 'w')
            f1.write(''.join(vcf_1_diff[contig_end_tag]))
            f1.close()
            os.system('grep -T -f %s %s %s --no-group-separator > %s' % (
                temp_output,
                vcf_file1,vcf_file2,
                output_file + '.temp'))
            os.system('sort -k3 -n %s | sort -k2 > %s' %
                      (output_file + '.temp', output_file + '.' + contig_end_tag)
                      )
            os.system('rm -rf %s' % (output_file+ '.temp'))
        else:
            os.system('rm -rf %s'%(output_file + '.' + contig_end_tag))
    return vcf_1_diff

def compare_vcf_all(vcf_1):
    vcf_input1 = load_vcf(vcf_1)
    vcf_1_all = vcf_1.split('.flt.snp.vcf')[0] + '.raw.vcf'
    FN_diff = compare_vcf(vcf_inputref, vcf_input1, vcf_ref_file, vcf_1_all, vcf_1 + '.FN') # vcf diff in ref not in 1, FN
    FP_diff = compare_vcf(vcf_input1, vcf_inputref, vcf_1_all, vcf_ref_file, vcf_1 + '.FP') # vcf diff in 1 not in ref, FP
    summary_file_output.append('%s\t%s\t%s\t%s\t%s\n'%(
        vcf_1,
        len(FN_diff['end']),len(FN_diff['notend']),
        len(FP_diff['end']), len(FP_diff['notend'])
    ))

# WGS
summary_file = output_dir + '/model.sum.loosecorrected.txt'
summary_file_output = []
for vcf_ref_file in vcf_ref:
    vcf_ref_file_name = os.path.split(vcf_ref_file)[-1].split('.snp.txt')[0].replace('.SNP.fasta','.all.SNP.fasta')
    vcf_1 = glob.glob(output_dir + vcf_ref_file_name + '1.fq.flt.snp.vcf.filtered.vcf')[0]
    vcf_2 = glob.glob(output_dir + vcf_ref_file_name + '1.fq.flt.snp.vcf.filteredrough.vcf')[0]
    # genome
    vcf_3 = glob.glob(output_dir + vcf_ref_file_name + '.genome.fasta.flt.snp.vcf.filtered.vcf')[0]
    vcf_4 = glob.glob(output_dir + vcf_ref_file_name + '.genome.fasta.corrected.flt.snp.vcf.filtered.vcf')[0]
    vcf_inputref = load_vcf(vcf_ref_file,True)
    compare_vcf_all(vcf_1)
    compare_vcf_all(vcf_2)
    compare_vcf_all(vcf_3)
    compare_vcf_all(vcf_4)

f1=open(summary_file,'w')
f1.write(''.join(summary_file_output))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# step 1 comparing genomes VS corrected genomes
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import random

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/genome_test1'
input_script_sub2 = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/genome_test2'
input_script_sub3 = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/genome_test3'
input_script_sub4 = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/genome_test4'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round2/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test'
fastq_name = '_1.fastq'
genome_name = '.fasta.corrected.fasta'

try:
    os.mkdir(output_dir + '/genome_test')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/genome_test/data')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/genome_test/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/genome_test/bwa/0')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/genome_test/bwa/1')
except IOError:
    pass

try:
    os.mkdir(output_dir + '/genome_test/bwa/2')
except IOError:
    pass

try:
    os.mkdir(input_script)
except IOError:
    pass

os.system('#rm -rf %s %s %s %s'%(input_script_sub,input_script_sub2,input_script_sub3,input_script_sub4))

try:
    os.mkdir(input_script_sub)
except IOError:
    pass

try:
    os.mkdir(input_script_sub2)
except IOError:
    pass

try:
    os.mkdir(input_script_sub3)
except IOError:
    pass

try:
    os.mkdir(input_script_sub4)
except IOError:
    pass

# function
def run_vcf_curate(files,files2,database,tempbamoutput):
    # generate code
    cmds = 'rm -rf %s.fai\nbowtie2-build %s %s\n' % (database,database, database)
    cmds += 'bowtie2' + ' -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, files, files2, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    cmds += 'samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
        tempbamoutput, tempbamoutput)
    cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s.sorted.bam  | %s call  --threads %s -m > %s.raw.vcf\n' % (
        'bcftools', min(40, 40), database,
        tempbamoutput, 'bcftools', min(40, 40), tempbamoutput)
    cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
        'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
    cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
        'bcftools', tempbamoutput, tempbamoutput)
    cmds += 'rm -rf %s.*.unaligned\n' % (tempbamoutput)
    cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
    return cmds

def run_vcf(genome_file,database,tempbamoutput):
    # generate code
    # for genome files
    cmds = 'minimap2 -d %s.mmi %s\n'%(database,database)
    cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, genome_file, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s.sorted.bam | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
        'bcftools', min(40, 40), database,
        tempbamoutput, 'bcftools', min(40, 40), tempbamoutput)
    cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
        'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
    cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
        'bcftools', tempbamoutput, tempbamoutput)
    cmds += 'rm -rf %s.flt.vcf %s.bam %s.sorted.bam.bai %s.sorted.bam\n' % (tempbamoutput,tempbamoutput,tempbamoutput,tempbamoutput)
    return cmds

def run_vcf_WGS(files,files2,database,tempbamoutput):
    cmds = 'rm -rf %s.fai\nbowtie2-build %s %s\n' % (database, database, database)
    cmds += 'bowtie2' + ' -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, files, files2, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s.sorted.bam | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
        'bcftools', min(40, 40), database,
        tempbamoutput, 'bcftools', min(40, 40), tempbamoutput)
    cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
        'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
    cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
        'bcftools', tempbamoutput, tempbamoutput)
    cmds += 'rm -rf %s.flt.vcf %s.bam %s.bam.bai %s.sorted.bam\n' % (
    tempbamoutput, tempbamoutput, tempbamoutput, tempbamoutput)
    return cmds

# subsample genomes
donor_species_set = set()
for folder in genome_dir:
    allfasta = glob.glob(os.path.join(folder, '*%s' % genome_name))
    donor_species = '_'.join(os.path.split(folder)[-1].split('_')[0:3])
    if allfasta != [] and donor_species not in donor_species_set:
        correctedfasta = allfasta[0]
        fastafile = correctedfasta.replace(genome_name,'.fasta')
        genomefile = os.path.split(fastafile)[-1].split('.fasta')[0]
        os.system('#cp %s %s/'%(fastafile,
                               output_dir + '/genome_test/'))
        donor_species_set.add(donor_species)
        fastq_file = glob.glob(os.path.join(folder, 'fastq/%s%s' % (genomefile,fastq_name)))[0]
        fastq_file2 = glob.glob(os.path.join(folder, 'fastq/%s%s' % (genomefile, fastq_name.replace('1','2'))))[0]
        tempbamoutput = os.path.join(output_dir + '/genome_test/bwa/1',genomefile)
        cmds = run_vcf_WGS(fastq_file, fastq_file2, fastafile, tempbamoutput)
        tempbamoutput = os.path.join(output_dir + '/genome_test/bwa/1', genomefile + '.corrected')
        cmds +=run_vcf_WGS(fastq_file, fastq_file2, correctedfasta, tempbamoutput)
        f1 = open(os.path.join(input_script_sub3, 'test.%s.vcf.sh' % (donor_species)), 'a')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
        f1.close()

allgenomes = glob.glob(output_dir + '/genome_test/*.fasta')
for genome in allgenomes:
    # simulate fastq files for each ref genome
    genome_name = os.path.split(genome)[-1].split('.fasta')[0]
    donor_species = '_'.join(genome_name.split('_')[0:3])
    genome_fastq = output_dir + '/genome_test/data/' + genome_name
    temp_folder = output_dir + '/genome_test/data/' + genome_name + '_genome'
    genome_file = output_dir + '/genome_test/data/' + genome_name + '.genome.fasta'
    tempbamoutput = output_dir + '/genome_test/bwa/0/' + genome_name + '.genome.fasta'
    print(donor_species,genome_name,genome_fastq,temp_folder,genome_file)
    cmds = 'art_illumina -ss HS25 -sam -i %s -p -l 150 -f 50 -m 200 -s 10 -o %s\n' \
           % (genome, genome_fastq)
    cmds += 'spades.py --careful -1 %s -2 %s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' % \
            (genome_fastq + '1.fq', genome_fastq + '2.fq', temp_folder)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_folder, genome_file)
    cmds += 'rm -rf %s *.aln *.sam\n' %(temp_folder)
    # mapping assembly to simulated WGS to curate genome
    cmds += run_vcf_curate(genome_fastq + '1.fq', genome_fastq + '2.fq', genome_file, genome_file + '.curate')
    cmds2 = run_vcf(genome_file, genome, tempbamoutput)
    genome_file = genome_file + '.corrected.fasta'
    tempbamoutput = tempbamoutput + '.corrected.fasta'
    cmds2 += run_vcf(genome_file, genome, tempbamoutput)
    tempbamoutput = output_dir + '/genome_test/bwa/2/' + genome_name + '1.fq'
    cmds3 = run_vcf_WGS(genome_fastq + '1.fq', genome_fastq + '2.fq',genome,tempbamoutput)
    f1 = open(os.path.join(input_script_sub, 'test.%s.vcf.sh' % (donor_species)), 'a')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
    f1.close()
    f1 = open(os.path.join(input_script_sub2, 'test.%s.vcf.sh' % (donor_species)), 'a')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds2)))
    f1.close()
    f1 = open(os.path.join(input_script_sub4, 'test.%s.vcf.sh' % (donor_species)), 'a')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds3)))
    f1.close()

f1 = open(os.path.join(input_script, 'genometest1.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'genometest2.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub2, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'genometest3.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub3, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

f1 = open(os.path.join(input_script, 'genometest4.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub4, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# run genometest1.sh
# step 2 check SNPs
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
# set up path

input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/genome_test1/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/genome_test/data'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/genome_test/data')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/genome_test/data'
fastq_name = '.curate.flt.snp.vcf'
genome_split = '_g'

# set up cutoff
# good mapping
Major_alt_freq_cutoff = 0.7 # major alt freq in a sample
Sample_depth_cutoff = 3 # both forward and reverse reads cutoff in a sample

# good coverage
total_coverage_cutoff = 0.8 # at least X reads map to its original genome
genome_avg_coverage_cutoff = 10 # genome average coverage cutoff

# curation cutoff
end_cutoff = 10 # don't curate the first and last end_cutoff bp of a contig

Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def coverage_check(sh_file,Coverage):
    Coverage1 = []
    Coverage2 = []
    os.system('rm -rf %s/temp.sh %s/temp.err'%(input_script_split_sub,input_script_split_sub))
    os.system('grep \"bowtie2\" %s > %s/temp.sh' % (sh_file, input_script_split_sub))
    os.system('grep \"overall alignment rate\" %s.err > %s/temp.err' % (sh_file, input_script_split_sub))
    # load fastq name
    for lines in open('%s/temp.sh' % (input_script_split_sub), 'r'):
        try:
            filename = lines.split('>')[1].split('\n')[0].split(fastq_name + '.bam')[0]
            filename = os.path.split(filename)[-1]
            Coverage1.append(filename)
        except IndexError:
            pass
    # load reads alignment rate
    for lines in open('%s/temp.err' % (input_script_split_sub), 'r'):
        Coverage2.append(float(lines.split('%')[0]))
    i = 0
    for filename in Coverage1:
        # load average coverage of ref genome
        tempbamoutput = os.path.join(output_dir, filename + fastq_name)
        coverage_file = glob.glob(tempbamoutput + '.sorted.bam.avgcov')
        if coverage_file == []:
            cmds = 'samtools depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
                tempbamoutput, tempbamoutput)
            os.system(cmds)
            coverage_file = tempbamoutput + '.sorted.bam.avgcov'
        else:
            coverage_file = coverage_file[0]
        length_total = 0
        coverage_total = 0
        donor_species = filename.split(genome_split)[0]
        Donor_species.setdefault(donor_species, [[], 0])
        try:
            for lines in open(coverage_file, 'r'):
                length_contig = float(lines.split('\t')[1])
                length_total += length_contig
                coverage_contig = float(lines.split('\t')[2])
                coverage_total += coverage_contig * length_contig
            coverage_num = Coverage2[i] / 100
            avg_coverage = coverage_total / length_total
            i += 1
        except IOError:
            coverage_num = total_coverage_cutoff - 0.1
            avg_coverage = genome_avg_coverage_cutoff - 1
        temp_line = ('%s\t%.2f\t%.1f\t%.1f' % (filename, coverage_num, length_total, avg_coverage))
        if coverage_num < total_coverage_cutoff or avg_coverage < genome_avg_coverage_cutoff:
            # not qualified
            print(filename, coverage_num, avg_coverage)
            Donor_species[donor_species][0].append(filename)
            temp_line += ('\tnot_qualified\t\n')
        else:
            # qualified
            Donor_species[donor_species][1] += 1
            temp_line += ('\tqualified\t\n')
        Coverage.append(temp_line)

def SNP_check(lines,donor_species,vcf_file_list):
    # CHR, POS, REF, ALT, good assembly, qualified mapping
    lines_set = lines.split('\n')[0].split('\t')
    report_line = ''
    temp_report_line = ['T','T']
    need_curation = 'F'
    REF = lines_set[3]
    allels_set = [REF]
    if '.' not in lines_set[4]:
        allels_set += lines_set[4].split(',')
    Total_alleles = len(allels_set)
    for Subdepth_all in lines_set[9:]:
            Allels_frq = [0, 0, 0, 0]
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
            total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
            Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
            Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
            total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
            total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
            for num_allels in range(0, Total_alleles):
                allels = allels_set[num_allels]
                Subdepth_alleles = int(Subdepth[num_allels])
                if allels in Allels:
                    Allels_frq[Allels[allels]] += Subdepth_alleles
                else:
                    pass
            # find major alt and calculate frequency
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            MLF = Major_ALT[1] / total_sub_depth
            if REF != Major_ALT[0]:
                # wrong assembly
                temp_report_line[0] = 'F' # F: not good assembly
            if MLF < Major_alt_freq_cutoff or total_sub_depth_forward < Sample_depth_cutoff or \
                    total_sub_depth_reverse < Sample_depth_cutoff:
                temp_report_line[1] = 'F'  # F: bad mapping
            if temp_report_line[1] == 'T' and temp_report_line[0] == 'F':
                # not good assembly but good mapping quality can be curated
                need_curation = 'T'  # T: need curation
    if temp_report_line != ['T','T']:
        report_line += '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t\n' %(donor_species,lines_set[0],lines_set[1],
                                                    REF,','.join([ALT for ALT in allels_set if ALT != REF]),
                                                    temp_report_line[0],temp_report_line[1],need_curation,
                                                     Major_ALT[0], MLF,
                                                     total_sub_depth_forward, total_sub_depth_reverse)
        vcf_file_list.append(donor_species + '\t' + lines)
    return report_line

def loadreport(vcf_file_report,coverage_report):
    assembly_quality = dict()
    assembly_quality_chr = dict()
    report_curate = dict()
    coverage = dict()
    remove_list = []
    try:
        for lines in open(coverage_report,'r'):
            lines_set = lines.split('\t')
            genome = lines_set[0]
            quality = lines_set[4]
            coverage.setdefault(genome,quality)
    except FileNotFoundError:
        pass
    for lines in open(vcf_file_report,'r'):
        lines_set = lines.split('\t')
        genome = lines_set[0]
        if coverage.get(genome,'') == 'qualified' or coverage == {}:
            assembly_quality.setdefault(genome,[])
            CHR = lines_set[1]
            POS = lines_set[2]
            REF = lines_set[3]
            ALT_major = lines_set[8]
            Assembly = lines_set[5]
            Mapping_quality = lines_set[6]
            Need_curation = lines_set[7]
            genome_chr = '%s_%s'%(genome,CHR)
            report_curate['%s_%s_%s' % (genome, CHR, POS)] = lines
            if Need_curation == 'T':
                # can be curated
                assembly_quality_chr.setdefault(genome_chr,[])
                assembly_quality_chr[genome_chr].append([POS, REF, ALT_major])
                assembly_quality[genome].append(CHR)
            elif Assembly == 'F' or Mapping_quality == 'F':
                # wrong assembly or bad mapping quality
                assembly_quality_chr.setdefault(genome_chr, [])
                assembly_quality_chr[genome_chr].append([POS, REF, 'N']) # use ambiguous characters instead
                assembly_quality[genome].append(CHR)
            else:
                print(lines)
        else:
            remove_list.append(genome)
    return [assembly_quality,assembly_quality_chr,report_curate]

def checkREF(seq,position,REF):
    return seq[position - 1] == REF

def causeSNP(seq,position,ALT):
    seq = list(seq)
    seq[position - 1] = ALT
    return ''.join(seq)

def correct_genome(database,assembly_quality_genome,assembly_quality_chr,report_curate):
    Output = []
    Output_report = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        total_length = len(record_seq)
        if record_id in assembly_quality_genome:
            genome_chr = '%s_%s' % (genome, record_id)
            allposition = assembly_quality_chr.get(genome_chr,[])
            for a_position in allposition:
                POS, REF, ALT = a_position
                POS = int(POS)
                if POS > end_cutoff and POS < total_length-end_cutoff + 1:
                    # don't curate the first and last end_cutoff bp of a contig
                    if checkREF(record_seq,int(POS),REF):
                        record_seq = causeSNP(record_seq, int(POS),ALT)
                        Output_report.append(report_curate['%s_%s_%s' % (genome, record_id,POS)])
                    else:
                        print(genome_chr,allposition)
                        print('wrong SNP POS %s %s for %s in %s'%(POS,REF,record_id,database))
        Output.append('>%s\n%s\n'%(record_id,record_seq))
    f1 = open(database + '.corrected.fasta','w')
    f1.write(''.join(Output))
    f1.close()
    f1 = open(os.path.join(input_script, 'genometest1.assembly.sum.curated'), 'a')
    f1.write(''.join(Output_report))
    f1.close()

# check major alt
all_vcf_file=glob.glob(os.path.join(output_dir,'*%s'%(fastq_name)))
vcf_file_report = []
vcf_file_list = []
vcf_file_report.append('donor_species\tCHR\tPOS\tREF\tALT\tAssembly\tMapping_quality\tNeed_curation\tMajor_ALT\tMajor_ALT_frq\tDepth_F\tDepth_R\n')
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split(fastq_name)[0]
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            vcf_file_report.append(SNP_check(lines,donor_species,vcf_file_list))

f1 = open(os.path.join(input_script, 'genometest1.assembly.sum'), 'w')
f1.write(''.join(vcf_file_report))
f1.close()

f1 = open(os.path.join(input_script, 'genometest1.assembly.vcf'), 'w')
f1.write(''.join(vcf_file_list))
f1.close()

# load report
assembly_quality,assembly_quality_chr,report_curate = loadreport(os.path.join(input_script, 'genometest1.assembly.sum'),
           os.path.join(input_script, 'genometest1.coverage.sum'))

# output curated loci
f1 = open(os.path.join(input_script, 'genometest1.assembly.sum.curated'), 'w')
f1.write('donor_species\tCHR\tPOS\tREF\tALT\tAssembly\tMapping_quality\tNeed_curation\tMajor_ALT\tMajor_ALT_frq\tDepth_F\tDepth_R\n')
f1.close()

# correct genomes
for genome in assembly_quality:
    assembly_quality_genome = assembly_quality.get(genome, [])
    if assembly_quality_genome!=[]:
        database = glob.glob(os.path.join(genome_root,genome.split(fastq_name)[0]))
        if database!= []:
            database = database[0]
            correct_genome(database, assembly_quality_genome,assembly_quality_chr,report_curate)
        else:
            print('missing input fasta for %s'%(genome))
    else:
        os.system('cp %s %s.corrected.fasta'%(genome,genome))

################################################### END ########################################################
################################################### SET PATH ########################################################
# run genometest2.sh
# step 3 check SNPs of WGS and genomes
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
# set up path

input_script_split_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/genome_test2/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/genome_test/data'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/genome_test/data')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/genome_test/'
fastq_name = '.flt.snp.vcf'
genome_split = '_g'

# set up cutoff
# SNP cutoff for WGS
Major_alt_freq_cutoff = 0.7 # major alt freq in a sample
Sample_depth_cutoff = 3 # both forward and reverse reads cutoff in a sample
# SNP cutoff for genome
Major_alt_freq_cutoff2 = 0.9 # major alt freq in a genome, do not allow multiple homolougs genes

# curation cutoff
end_cutoff = 70 # don't curate the first and last end_cutoff bp of a contig

Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def contig_end(CHR,POS):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        total_length = CHR.split('length_')[1].split('_cov')[0]
    total_length = int(total_length)
    if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
        return True
    else:
        return False

def SNP_check(lines):
    # CHR, POS, REF, ALT, good assembly, qualified mapping
    lines_set = lines.split('\n')[0].split('\t')
    CHR, POS = lines_set[0:2]
    if not contig_end(CHR,POS):
        REF = lines_set[3]
        allels_set = [REF]
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        for Subdepth_all in lines_set[9:]:
            Allels_frq = [0, 0, 0, 0]
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
            total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
            Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
            Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
            total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
            total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
            for num_allels in range(0, Total_alleles):
                allels = allels_set[num_allels]
                Subdepth_alleles = int(Subdepth[num_allels])
                if allels in Allels:
                    Allels_frq[Allels[allels]] += Subdepth_alleles
                else:
                    pass
            # find major alt and calculate frequency
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            MLF = Major_ALT[1] / total_sub_depth
            if MLF >= Major_alt_freq_cutoff and total_sub_depth_forward >= Sample_depth_cutoff and \
                    total_sub_depth_reverse >= Sample_depth_cutoff and REF != Major_ALT[0]:
                return lines
    return ''

def SNP_check_genome(lines):
    # CHR, POS, REF, ALT, good assembly, qualified mapping
    lines_set = lines.split('\n')[0].split('\t')
    CHR, POS = lines_set[0:2]
    if not contig_end(CHR,POS):
        REF = lines_set[3]
        allels_set = [REF]
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        for Subdepth_all in lines_set[9:]:
            Allels_frq = [0, 0, 0, 0]
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
            total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
            Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
            Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
            for num_allels in range(0, Total_alleles):
                allels = allels_set[num_allels]
                Subdepth_alleles = int(Subdepth[num_allels])
                if allels in Allels:
                    Allels_frq[Allels[allels]] += Subdepth_alleles
                else:
                    pass
            # find major alt and calculate frequency
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            MLF = Major_ALT[1] / total_sub_depth
            if MLF >= Major_alt_freq_cutoff2 and REF != Major_ALT[0]:
                return lines
    return ''

# check SNP for genome
all_vcf_file=glob.glob(os.path.join(output_dir + '/bwa/0/','*%s'%(fastq_name)))
snp_report = dict()
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split(fastq_name)[0]
    snp_report.setdefault(donor_species,0)
    vcf_file_filtered = []
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            output_line = SNP_check_genome(lines)
            if output_line != '':
                vcf_file_filtered.append(output_line)
                snp_report[donor_species] += 1
    f1 = open(vcf_file + '.filtered.vcf', 'w')
    f1.write(''.join(vcf_file_filtered))
    f1.close()

snp_report_sum = ['genome\tcorrected\tuncorrected\t\n']
for donor_species in snp_report:
    if 'corrected' in donor_species:
        snp_report_sum.append('%s\t%s\t%s\t\n'%(donor_species,
                                            snp_report[donor_species],
                                                snp_report[donor_species.split('.corrected')[0]]))

f1 = open(os.path.join(input_script, 'genometest1.genome.sum.txt'), 'w')
f1.write(''.join(snp_report_sum))
f1.close()

# check SNP for WGS
all_vcf_file=glob.glob(os.path.join(output_dir + '/bwa/1/','*%s'%(fastq_name)))
snp_report = dict()
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split(fastq_name)[0]
    snp_report.setdefault(donor_species,0)
    vcf_file_filtered = []
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            output_line = SNP_check(lines)
            if output_line != '':
                vcf_file_filtered.append(output_line)
                snp_report[donor_species] += 1
    f1 = open(vcf_file + '.filtered.vcf', 'w')
    f1.write(''.join(vcf_file_filtered))
    f1.close()

snp_report_sum = ['genome\tcorrected\tuncorrected\t\n']
for donor_species in snp_report:
    if 'corrected' in donor_species:
        snp_report_sum.append('%s\t%s\t%s\t\n'%(donor_species,
                                            snp_report[donor_species],
                                                snp_report[donor_species.split('.corrected')[0]]))

f1 = open(os.path.join(input_script, 'genometest1.WGS.sum.txt'), 'w')
f1.write(''.join(snp_report_sum))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# step 4 modelling SNPs and test 2 methods
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import random

input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/SNP_test_step1'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate'
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test'
genome_name = '.fasta.corrected.fasta'
fastq_name = '_1.fastq'

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# Set up N or S
N_S_set = dict()
N_S_set['N']=0
N_S_set['S']=1
purines=['A','G']
pyrimidines=['C','T']
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
# Set up mutation rate
mut_rate = 1/1000/100 # mutation rate 1 per 10000 bp
mut_time = 1 # generate mut_time strains with SNPs -> a population with mut_time + 1 strains
mut_set = [10,25,50,100,250,500,1000,2500]
if len(mut_set)!= 0 :
    mut_time=len(mut_set)

cause_SNP = True
minor_allele_freq_cutoff = 0.2
pop_size = 30 #not used
fastq_length = 150
SNP_times = 1000 # snp reads depth as 1500, both forward and reverse

try:
    os.mkdir(output_dir+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa/0')
except IOError:
    pass

try:
    os.mkdir(input_script)
except IOError:
    pass

os.system('rm -rf %s' %(input_script_sub))
try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# function
def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

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

def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        if ALT_frq > 0:
            ALT_set.setdefault(ALT_frq, set())
            ALT_set[ALT_frq].add(alleles)
            ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def transitions(REF,ALT):
    if REF in pyrimidines:
        REF = complement[REF]
        ALT = complement[ALT]
    return '%s-%s'%(REF,ALT)

def contig_to_gene(CHR, POS):
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

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        Ref_seq.setdefault(record_id, record_seq)
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

def modelSNP(seq,Chr,num_mut_chr):
    total_length = len(seq)
    position_set = random.choices(range(0, total_length-1), k=num_mut_chr)
    seq = list(seq)
    SNP_output = []
    for position in position_set:
        gene_info = contig_to_gene(Chr, position)
        REF = seq[position]
        temp_ALT = ['A', 'T', 'G', 'C']
        try:
            temp_ALT.remove(REF)
        except ValueError:
            pass
        ALT = random.choices(temp_ALT, k=1)[0]
        seq[position] = ALT
        temp_line = [Chr,str(position+1),REF,ALT,'Other','None']
        if gene_info != []:
            # a gene
            Chr_gene, position_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
            codon_start = position_gene - 1 - int((position_gene - 1) % 3)
            if codon_start <= position_gene - 1:
                Ref_seq_chr = Ref_seq[Chr_gene]
                SNP_seq_chr = Ref_seq_chr
                Ref_seq_chr = causeSNP(Ref_seq_chr, position_gene, REF,Reverse_chr)
                Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                if len(Ref_seq_codon) == 3:
                    Ref_seq_aa = translate(Ref_seq_codon)[0]
                    SNP_seq_chr = causeSNP(SNP_seq_chr, position_gene, ALT, Reverse_chr)
                    SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                    SNP_seq_aa = translate(SNP_seq_codon)[0]
                    temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                    temp_line[-1]=''.join([Ref_seq_aa,SNP_seq_aa])
                    temp_line[-2]=temp_NorS
        SNP_output.append('\t'.join(temp_line)+'\n')
    return [''.join(seq), SNP_output]

def run_vcf_WGS(files,files2,database,tempbamoutput):
    cmds = 'rm -rf %s.fai\nbowtie2-build %s %s\n' % (database, database, database)
    cmds += 'bowtie2' + ' -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, files, files2, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s.sorted.bam | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
        'bcftools', min(40, 40), database,
        tempbamoutput, 'bcftools', min(40, 40), tempbamoutput)
    cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
        'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
    cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
        'bcftools', tempbamoutput, tempbamoutput)
    cmds += 'rm -rf %s.flt.vcf %s.bam %s.sorted.bam.bai %s.sorted.bam\n' % (
    tempbamoutput, tempbamoutput, tempbamoutput, tempbamoutput)
    return cmds

def run_vcf(genome_file,database,tempbamoutput):
    # generate code
    # for genome files
    cmds = 'minimap2 -d %s.mmi %s\n'%(database,database)
    cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
        min(40, 40), database, genome_file, 'samtools', min(40, 40),
        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
        tempbamoutput)
    cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s.sorted.bam | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
        'bcftools', min(40, 40), database,
        tempbamoutput, 'bcftools', min(40, 40), tempbamoutput)
    cmds += '%s filter --threads %s -s LowQual %s.raw.vcf > %s.flt.vcf \n' % (
        'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
    cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
        'bcftools', tempbamoutput, tempbamoutput)
    cmds += 'rm -rf %s.flt.vcf %s.bam %s.sorted.bam.bai %s.sorted.bam\n' % (tempbamoutput,tempbamoutput,tempbamoutput,tempbamoutput)
    return cmds

def modelSNPall(Input_seq,Input_id,num_mut):
    Output = []
    Output_SNP = []
    chr_set = random.sample(Input_id*num_mut, k=num_mut)
    unique_chr_set = list(set(chr_set))
    for chr in Input_id:
        if chr in unique_chr_set:
            # mutated
            num_mut_chr = chr_set.count(chr)
            newseq, newoutput = modelSNP(Input_seq[chr], chr,num_mut_chr)
            Output_SNP += newoutput
        else:
            # not mutated
            newseq = Input_seq[chr]
        Output.append('>%s\n%s\n' % (chr,
                                     newseq))
    # output mutated genome
    output_fasta = os.path.join(output_dir, '%s.%s.SNP.fasta' % (genome_filename,num_mut))
    f1 = open(output_fasta, 'w')
    f1.write(''.join(Output))
    f1.close()
    f1 = open(output_fasta + '.snp.txt', 'w')
    f1.write(''.join(Output_SNP))
    f1.close()
    return 'done'

allgenomes = glob.glob(os.path.join(output_dir,'*%s'%(genome_name)))
for genome in allgenomes:
    # load genome
    database_file = genome + '.fna'
    genome_filename = os.path.split(genome)[-1].split(genome_name)[0]
    os.system('#prodigal -q -i %s -d %s'%(genome,database_file))
    Ref_seq, Mapping,Mapping_loci,Reverse = loaddatabase(database_file)
    Input_seq = dict()
    Input_id = []
    for record in SeqIO.parse(genome, 'fasta'):
        seq = str(record.seq)
        seq_length = len(seq)
        if seq_length >= 8000:
            Input_seq.setdefault(str(record.id), seq)
            Input_id.append(str(record.id))
    # cause SNP
    fastq_file = genome.split(genome_name)[0]
    files = fastq_file + fastq_name
    files2 = fastq_file + fastq_name.replace('1','2')
    if len(mut_set) != 0:
        mut_time = len(mut_set)
    while mut_time > 0:
        num_mut = mut_set[mut_time - 1]
        # cause SNP
        if cause_SNP:
            # simulate fastq files for mutated strains
            modelSNPall(Input_seq, Input_id, num_mut)
        mutated_genome = os.path.join(output_dir, '%s.%s.SNP.fasta'%(genome_filename,num_mut))
        tempbamoutput = os.path.join(output_dir, 'bwa/0/%s.%s.SNP.fasta.fq'%(genome_filename,num_mut))
        cmds = run_vcf_WGS(files, files2, mutated_genome, tempbamoutput)
        tempbamoutput = os.path.join(output_dir, 'bwa/0/%s.%s.SNP.fasta.genome' % (genome_filename,num_mut))
        cmds += run_vcf(genome, mutated_genome, tempbamoutput)
        f1 = open(os.path.join(input_script_sub, '%s.%s.vcf.sh' % (genome_filename,mut_time)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
        f1.close()
        mut_time -= 1

f1 = open(os.path.join(input_script, 'alltestvcfstep1.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))

f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# run alltestvcfstep1.sh
# step 3 check SNPs of WGS and genomes
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
# set up path

input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/SNP_currate/'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/SNP_test/'
genome_dir = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/SNP_test/'
fastq_name = '.flt.snp.vcf'

# set up cutoff
# SNP cutoff for WGS
Major_alt_freq_cutoff = 0.7 # major alt freq in a sample
Sample_depth_cutoff = 1 # both forward and reverse reads cutoff in a sample
# SNP cutoff for genome
Major_alt_freq_cutoff2 = 0.9 # major alt freq in a genome, do not allow multiple homolougs genes

# curation cutoff
end_cutoff = 70 # don't curate the first and last end_cutoff bp of a contig

Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
# set up functions
def ALT_freq(Allels_count):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
    return [Major_ALT,Minor_ALT]

def contig_end(CHR,POS):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        total_length = CHR.split('length_')[1].split('_cov')[0]
    total_length = int(total_length)
    if int(POS) <= end_cutoff or int(POS) >= total_length - end_cutoff + 1:
        return True
    else:
        return False

def SNP_check(lines):
    # CHR, POS, REF, ALT, good assembly, qualified mapping
    lines_set = lines.split('\n')[0].split('\t')
    CHR, POS = lines_set[0:2]
    if not contig_end(CHR,POS):
        REF = lines_set[3]
        allels_set = [REF]
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        for Subdepth_all in lines_set[9:]:
            Allels_frq = [0, 0, 0, 0]
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
            total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
            Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
            Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
            total_sub_depth_forward = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_forward)
            total_sub_depth_reverse = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth_reverse)
            for num_allels in range(0, Total_alleles):
                allels = allels_set[num_allels]
                Subdepth_alleles = int(Subdepth[num_allels])
                if allels in Allels:
                    Allels_frq[Allels[allels]] += Subdepth_alleles
                else:
                    pass
            # find major alt and calculate frequency
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            MLF = Major_ALT[1] / total_sub_depth
            if MLF >= Major_alt_freq_cutoff and (total_sub_depth_forward >= Sample_depth_cutoff or \
                    total_sub_depth_reverse >= Sample_depth_cutoff) and REF != Major_ALT[0]:
                return lines
    return ''

def SNP_check_genome(lines):
    # CHR, POS, REF, ALT, good assembly, qualified mapping
    lines_set = lines.split('\n')[0].split('\t')
    CHR, POS = lines_set[0:2]
    if not contig_end(CHR,POS):
        REF = lines_set[3]
        allels_set = [REF]
        if '.' not in lines_set[4]:
            allels_set += lines_set[4].split(',')
        Total_alleles = len(allels_set)
        for Subdepth_all in lines_set[9:]:
            Allels_frq = [0, 0, 0, 0]
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
            total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
            Subdepth_forward = Subdepth_all.split(':')[-3].split(',')
            Subdepth_reverse = Subdepth_all.split(':')[-2].split(',')
            for num_allels in range(0, Total_alleles):
                allels = allels_set[num_allels]
                Subdepth_alleles = int(Subdepth[num_allels])
                if allels in Allels:
                    Allels_frq[Allels[allels]] += Subdepth_alleles
                else:
                    pass
            # find major alt and calculate frequency
            Major_ALT, Minor_ALT = ALT_freq(Allels_frq)
            MLF = Major_ALT[1] / total_sub_depth
            if MLF >= Major_alt_freq_cutoff2 and REF != Major_ALT[0]:
                return lines
    return ''

# check SNP for genome
all_vcf_file=glob.glob(os.path.join(output_dir + '/bwa/0/','*genome%s'%(fastq_name)))
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split(fastq_name)[0]
    vcf_file_filtered = []
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            output_line = SNP_check_genome(lines)
            if output_line != '':
                vcf_file_filtered.append(output_line)
    f1 = open(vcf_file + '.filtered.vcf', 'w')
    f1.write(''.join(vcf_file_filtered))
    f1.close()

# check SNP for WGS
all_vcf_file=glob.glob(os.path.join(output_dir + '/bwa/0/','*fq%s'%(fastq_name)))
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split(fastq_name)[0]
    vcf_file_filtered = []
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            output_line = SNP_check(lines)
            if output_line != '':
                vcf_file_filtered.append(output_line)
    f1 = open(vcf_file + '.filtered.vcf', 'w')
    f1.write(''.join(vcf_file_filtered))
    f1.close()

################################################## END ########################################################
################################################### SET PATH ########################################################
# compare genome call SNPs VS WGS call SNPs -> SNP_model_compare.py
