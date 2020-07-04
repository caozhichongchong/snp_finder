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
                            cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -B -Ou -d3000 -f %s %s.sorted.bam  | %s call --ploidy 1 --threads %s -m > %s.raw.vcf\n' % (
                                'bcftools', min(40, 40), database,
                                tempbamoutput, 'bcftools', min(40, 40), tempbamoutput)
                            cmds += '%s filter --threads %s -i \'QUAL>30\' %s.raw.vcf > %s.flt.vcf \n' % (
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

# function

def loadreport(vcf_file_report,coverage_report):
    assembly_quality = dict()
    assembly_quality_chr = dict()
    coverage = dict()
    remove_list = []
    for lines in open(coverage_report,'r'):
        lines_set = lines.split('\t')
        genome = lines_set[0]
        quality = lines_set[4]
        coverage.setdefault(genome,quality)
    for lines in open(vcf_file_report,'r'):
        lines_set = lines.split('\t')
        genome = lines_set[0]
        if coverage.get(genome,'') == 'qualified':
            assembly_quality.setdefault(genome,[])
            CHR = lines_set[1]
            POS = lines_set[2]
            REF = lines_set[3]
            ALT = lines_set[4]
            Assembly = lines_set[5]
            Mapping_quality = lines_set[6]
            Need_curation = lines_set[7]
            genome_chr = '%s_%s'%(genome,CHR)
            if Need_curation == 'T':
                assembly_quality_chr.setdefault(genome_chr,[])
                assembly_quality_chr[genome_chr].append([POS, REF, ALT])
                assembly_quality[genome].append(CHR)
            elif Assembly == 'F' or Mapping_quality == 'F':
                assembly_quality_chr.setdefault(genome_chr, [])
                assembly_quality_chr[genome_chr].append([POS, REF, 'N'])
                assembly_quality[genome].append(CHR)
            else:
                print(lines)
        else:
            remove_list.append(genome)
    return [assembly_quality,assembly_quality_chr]

def checkREF(seq,position,REF):
    return seq[position - 1] == REF

def causeSNP(seq,position,ALT):
    seq = list(seq)
    seq[position - 1] = ALT
    return ''.join(seq)

def correct_genome(database,assembly_quality_genome,assembly_quality_chr):
    Output = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        if record_id in assembly_quality_genome:
            genome_chr = '%s_%s' % (genome, record_id)
            allposition = assembly_quality_chr.get(genome_chr,[])
            for a_position in allposition:
                POS, REF, ALT = a_position
                if checkREF(record_seq,int(POS),REF):
                    record_seq = causeSNP(record_seq, int(POS),ALT)
                else:
                    print(genome_chr,allposition)
                    print('wrong SNP POS %s %s for %s in %s'%(POS,REF,record_id,database))
        Output.append('>%s\n%s\n'%(record_id,record_seq))
    f1 = open(database+'.corrected.fasta','w')
    f1.write(''.join(Output))
    f1.close()

# load report
assembly_quality,assembly_quality_chr = loadreport(os.path.join(input_script, 'SNP_currate1.assembly.sum'),
           os.path.join(input_script, 'SNP_currate1.coverage.sum'))

# correct genomes
for genome in assembly_quality:
    assembly_quality_genome = assembly_quality.get(genome, [])
    if assembly_quality_genome!=[]:
        donor_species = '_'.join(genome.split(genome_split)[:-1])
        donor = donor_species.split('_')[0]
        database = glob.glob(os.path.join(genome_root,'%s*/%s.fasta'%(donor_species,genome)))
    if database!= []:
        database = database[0]
        correct_genome(database, assembly_quality_genome,assembly_quality_chr)
    else:
        print('missing input fasta for %s'%(genome))

################################################### END ########################################################
# step 4 all spades map to all genomes
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
                    cmds += 'minimap2' + ' -ax asm5 -N 1 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                        min(40, 40), database, files, 'samtools', min(40, 40),
                        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
                        tempbamoutput)
                    # fna for genes
                    cmds += 'minimap2' + ' -ax asm5 -N 1 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.fna.bam\n%s sort -@ %s %s.fna.bam -o %s.fna.sorted.bam\n%s index -@ %s %s.fna.sorted.bam\n' % (
                        min(40, 40), database2, files, 'samtools', min(40, 40),
                        tempbamoutput, 'samtools', min(40, 40), tempbamoutput, tempbamoutput, 'samtools', min(40, 40),
                        tempbamoutput)
                    cmds += 'rm -rf %s.bam %s.fna.bam\n' % (tempbamoutput,tempbamoutput)
                    sub_samples.append(tempbamoutput + '.sorted.bam')
                    sub_samples2.append(tempbamoutput + '.fna.sorted.bam')
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
            cmds = '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -B -Ou -d3000 -f %s %s | %s call --ploidy 1 --threads %s -m > %s.raw.vcf\n' % (
                'bcftools', min(40, 40), database,
                ' '.join(sub_samples), 'bcftools', min(40, 40), tempbamoutput)
            cmds += '%s filter --threads %s -i \'QUAL>30\' %s.raw.vcf > %s.flt.vcf \n' % (
                'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
            cmds += '%s view -H -v snps %s.flt.vcf > %s.flt.snp.vcf \n' % (
                'bcftools', tempbamoutput, tempbamoutput)
            cmds += 'rm -rf %s.flt.vcf\n' % (tempbamoutput)
            cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -B -Ou -d3000 -f %s %s | %s call --ploidy 1 --threads %s -m > %s.fna.raw.vcf\n' % (
                'bcftools', min(40, 40), database2,
                ' '.join(sub_samples2), 'bcftools', min(40, 40), tempbamoutput)
            cmds += '%s filter --threads %s -i \'QUAL>30\' %s.fna.raw.vcf > %s.fna.flt.vcf \n' % (
                'bcftools', min(40, 40), tempbamoutput, tempbamoutput)
            cmds += '%s view -H -v snps %s.fna.flt.vcf > %s.fna.flt.snp.vcf \n' % (
                'bcftools', tempbamoutput, tempbamoutput)
            cmds += 'rm -rf %s.fna.flt.vcf\n' % (tempbamoutput)
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
        return allels_set[1]
    else:
        return allels_set[0]

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
        REF = curate_REF(allels_set, Depth4)  # as the major alt in the population
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
                        and Total_qualify_SNP >= 1 and Total_qualify_SNP < Total_qualify \
                        and Total_qualify_notSNP > 0:
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
all_vcf_file=glob.glob(os.path.join(output_dir_merge + '/bwa/0','*.fna.flt.snp.vcf'))
for set_num in range(0,len(reference_set)):
    reference_name = reference_set[set_num]
    output_name = outputname_set[set_num]
    for vcf_file in all_vcf_file:
        SNP_presence_cutoff = SNP_presence_cutoff2 # for group of samples  -> used
        SNP_presence_sample_cutoff = SNP_presence_sample_cutoff2
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
                            if Total <= 3:
                                SNP_presence_cutoff = 1  # for a small group of genomes
                                SNP_presence_sample_cutoff = 2
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
                if Total <= 3:
                    SNP_presence_cutoff = 1  # for a small group of genomes
                    SNP_presence_sample_cutoff = 2
            if Depth / Total >= SNP_presence_cutoff:
                # average depth in all samples cutoff
                if "INDEL" not in lines_set[7] \
                        and (lines_set[6] != 'LowQual'):
                    if Step == 1:
                        CHR_old, POS_old = SNP_check_all(lines_set, '',
                                                         CHR_old, POS_old, reference_name,
                                                         SNP_presence_cutoff, SNP_presence_sample_cutoff)
                    elif Step == 2:
                        CHR_old, POS_old = SNP_check_all(lines_set, '',
                                                         CHR_old, POS_old, reference_name,
                                                         SNP_presence_cutoff, SNP_presence_sample_cutoff)
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
            outputcov(output_name, list(vcf_file_POS_candidate), Sample_subset)

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

# Question: why genome mapping shows SNPs that have no alternative alleles in raw.vcf of WGS mapping?
################################################### CLEAN UP ########################################################
os.system('rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/SNP_currate2')