# step 1 modelling SNPs and test 2 methods
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import random
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of folders of WGS of each species",
                      type=str, default='.',
                      metavar='input/')
required.add_argument("-fa",
                      help="file extension of fasta files",
                      type=str, default='.fasta',
                      metavar='.corrected.fasta')
required.add_argument("-fq",
                      help="file extension of fastq files",
                      type=str, default='_1.fastq',
                      metavar='_1.fastq')
# optional output setup
optional.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='scripts/',
                      metavar='scripts/')
optional.add_argument("-o",
                      help="a folder to store all output",
                      type=str, default='snp_output/',
                      metavar='snp_output/')
optional.add_argument("-rmhm",
                      help="remove homologous regions reported in genome_removeduplicate.py output",
                      type=str, default='False',
                      metavar='.duplicate_50kmer.txt')
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 40)",
                      metavar="1 or more", action='store', default=40, type=int)
optional.add_argument("-indel",
                      help="whether to insert indels",
                      type=str, default='True',
                      metavar='False or True')

################################################## Definition ########################################################
args = parser.parse_args()
input_script = args.s
genome_root = args.i
output_dir = args.o + '/SNP_model'
genome_name = args.fa
fastq_name=args.fq
fastq_name2=args.fq.replace('1','2')
cause_SNP = True
mapping_file = True
time_evaluation = False
homologous_region_length = 50
if time_evaluation:
    input_script_sub = '%s/SNP_model'%(input_script)
else:
    input_script_sub = '%s/SNP_model_onetime'%(input_script)
if args.indel == 'False':
    output_dir += '_noindel'
    input_script_sub += '_noindel'

latest_mapper = glob.glob('%s/mapper-1*.jar'%(args.s))[0]
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
mut_set = [0,5,10,50,100,1000,5000,10000,50000,100000,200000]
#mut_set = [0,5,10,50,100,1000,5000,10000,50000,100000,200000,300000,400000,500000]
indel_time = 100 # how many indels in a genome
if 'human' in args.i:
    mut_set = [int(x) for x in [0,1e2,1e3,1e4,1e5,1e6,1e7]]
if args.indel == 'False':
    indel_time = 0 # how many indels in a genome
#indel_orf = [-10,-7,-4, 4, 7, 10]
indel_nonorf = list(range(2,17))
indel_nonorf.extend(list(range(-17,-2)))
indel_orf = indel_nonorf
try:
    os.mkdir(output_dir)
except IOError:
    pass

try:
    os.mkdir(output_dir + '/data')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/bwa')
except IOError:
    pass

try:
    os.mkdir(output_dir+'/merge')
except IOError:
    pass

os.system('rm -r %s'%(input_script_sub))
try:
    os.mkdir(input_script_sub)
except IOError:
    pass

# function
def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq[position - 1]=ALT
    return seq

def translate(seq):
    seq = Seq(''.join(seq))
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
    return []

def loaddatabase(database_aa,database):
    # load database seq
    Length = []
    reference_database = os.path.split(database_aa)[-1]
    print('reference database_aa set as %s' % (reference_database))
    Ref_seq = dict()
    Input_seq = dict()
    Input_id = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        seq_length = len(record_seq)
        Input_seq.setdefault(record_id, record_seq)
        Input_id.append(record_id)
        Length.append(seq_length)
    for record in SeqIO.parse(database_aa, 'fasta'):
        record_id = str(record.id)
        record_seq = list(str(record.seq))
        Ref_seq.setdefault(record_id, record_seq)
    return [Ref_seq,Length,Input_seq,Input_id]

def load_homologous_regions(homologous_file):
    homologous_regions = dict()
    for lines in open(homologous_file, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        CHR, POS, No_match = lines_set
        homologous_regions.setdefault(CHR,set())
        POS = int(POS)
        for i in range(POS, POS + homologous_region_length + 1):
            homologous_regions[CHR].add(i)
    return homologous_regions

def modelindel(seq,Chr,indel_set):
    SNP_output = []
    indel_set.sort()
    for position in indel_set:
        total_length = len(seq)
        if position < total_length:
            REF = seq[position]
            gene_info = contig_to_gene(Chr, position)
            temp_ALT = ['A', 'T', 'G', 'C']
            indel_size = random.choices(indel_nonorf, k=1)[0]
            if False:
                if gene_info != []:
                    # a gene
                    # indel = + 3*n
                    indel_size = random.choices(indel_orf, k=1)[0]
                else:
                    # not a gene
                    indel_size = random.choices(indel_nonorf, k=1)[0]
            if indel_size > 0:# insertion on ref
                ALT = random.choices(temp_ALT, k=indel_size)
                seq = seq[:position + 1] + ALT + seq[position + 1:]
                # including the REF, add insertion after the REF
                temp_line = [Chr, str(position + 1), REF, REF + ''.join(ALT)]
            else:# deletion on ref
                REF = ''.join(seq[(position+indel_size):position])
                del seq[(position+indel_size):position]
                temp_line = [Chr, str(position + 1 +indel_size), REF, '-'*(-indel_size)]
            SNP_output.append('\t'.join(temp_line) + '\n')
        else:
            print('position %s out of the reference %s'%(position,total_length))
    return [seq, SNP_output]

def modelSNP(seq,Chr,num_mut_chr,num_indel_chr):
    total_length = len(seq)
    # indel modelling
    indel_output = []
    seq = list(seq)
    homologous_regions_chr = homologous_regions.get(Chr, set())
    if num_indel_chr > 0:
        candidate_position = [i for i in range(0, total_length) if seq[i] not in ['-','N'] and i not in homologous_regions_chr]
        indel_set = random.sample(candidate_position, k=num_indel_chr)
        seq, indel_output = modelindel(seq,Chr,indel_set)
    total_length = len(seq)
    # SNP modelling
    # if add indels, need re-adjust homologous_regions_chr!
    candidate_position = [i for i in range(0, total_length) if seq[i] not in ['-', 'N']]
    num_mut_chr = min(num_mut_chr,len(candidate_position))
    position_set = random.sample(candidate_position, k=num_mut_chr)
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
        if False and gene_info != []:
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
    return [''.join(seq),SNP_output,indel_output]

def run_vcf_WGS(files,files2,database,tempbamoutput):
    # generate code
    cmds = 'bowtie2-build %s %s\n'%(database,database)
    try:
        f1 = open('%s.sorted.bam' % (tempbamoutput),'r')
    except IOError:
        cmds += 'bowtie2 --threads %s -x %s -1 %s -2 %s |samtools view -@ %s -S -b >%s.bam\nsamtools sort -@ %s %s.bam -o %s.sorted.bam\nsamtools index -@ %s %s.sorted.bam\n' % (
            min(40, args.t), database, files, files2,  min(40, args.t),
            tempbamoutput,  min(40, args.t), tempbamoutput, tempbamoutput, min(40, args.t),
            tempbamoutput)
        cmds += 'rm -r %s.sam %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def merge_sample(database,vcfoutput,allsam):
    cmds = ''
    cmds += 'bcftools mpileup --ff UNMAP,QCFAIL,SECONDARY --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -Ou -d3000 -A -f  %s %s | bcftools call -Ov -A -M -c --threads %s > %s.raw.vcf\n' % (
             min(40, args.t), database,
            ' '.join(allsam), min(40, args.t), vcfoutput)
    if not time_evaluation:
        cmds += 'bcftools view -H -v snps -i \'QUAL>=20 && MIN(DP)>=3\' %s.raw.vcf > %s.flt.snp.vcf \n' % (
                vcfoutput, vcfoutput)
        cmds += 'bcftools view -H -v indels -i \'QUAL>=20 && MIN(DP)>=3\' %s.raw.vcf > %s.flt.indel.vcf \n' % (
                vcfoutput, vcfoutput)
    cmds += 'rm %s\n'%(' '.join(allsam))
    if '.0.SNP.fasta.bowtie' not in vcfoutput:
        cmds += 'rm %s.raw.vcf\n' % (vcfoutput)
    return cmds

def run_minimap(files,files2,database,tempbamoutput):
    cmds = 'minimap2 -d %s.mmi %s \n' % (database, database)
    cmds += 'minimap2 -ax sr -N 10 -p 0.9 -t %s %s.mmi %s %s >%s.sam\npy39\nsamtools view -@ %s -S -b %s.sam >%s.bam\nsamtools sort -@ %s %s.bam -o %s.sorted.bam\nsamtools index -@ %s %s.sorted.bam\n' % (
        min(40, args.t), database, files, files2, tempbamoutput, min(40, args.t),tempbamoutput,
            tempbamoutput,  min(40, args.t), tempbamoutput, tempbamoutput, min(40, args.t),
            tempbamoutput)
    cmds += 'rm -r %s.sam %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def run_bwa(files,files2,database,tempbamoutput):
    cmds = 'bwa index %s\n'%(database)
    cmds += 'bwa mem %s %s %s > %s.sam\npy39\nsamtools view -@ %s -S -b %s.sam >%s.bam\nsamtools sort -@ %s %s.bam -o %s.sorted.bam\nsamtools index -@ %s %s.sorted.bam\n' % (
         database, files, files2, tempbamoutput, min(40, args.t),tempbamoutput,
            tempbamoutput,  min(40, args.t), tempbamoutput, tempbamoutput, min(40, args.t),
            tempbamoutput)
    cmds += 'rm -r %s.sam %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def modelSNPall(Input_seq, Input_id, Length,num_mut,database_name):
    Output = []
    Output_SNP = []
    Output_indel = []
    chr_set = random.choices(Input_id, k=num_mut,weights=Length)
    chr_set_indel = random.choices(Input_id, k=indel_time, weights=Length)
    unique_chr_set = list(set(chr_set))
    unique_chr_set_indel = list(set(chr_set_indel))
    for chr in Input_id:
        if chr in unique_chr_set:
            # mutated
            num_mut_chr = chr_set.count(chr)
            num_indel_chr = 0
            if chr in unique_chr_set_indel:
                num_indel_chr = chr_set_indel.count(chr)
            newseq, newoutput, newoutputindel = modelSNP(Input_seq[chr], chr,num_mut_chr,num_indel_chr)
            Output_SNP += newoutput
            Output_indel += newoutputindel
        elif chr in unique_chr_set_indel:
            # not mutated, add indels
            num_indel_chr = chr_set_indel.count(chr)
            newseq, newoutput, newoutputindel = modelSNP(Input_seq[chr], chr, 0, num_indel_chr)
            Output_indel += newoutputindel
        else:
            # not mutated, no indels
            newseq = Input_seq[chr]
        Output.append('>%s\n%s\n' % (chr,
                                     newseq))
    # output mutated genome
    output_fasta = os.path.join(output_dir, 'data/%s.%s.SNP.fasta' % (database_name,num_mut))
    f1 = open(output_fasta, 'w')
    f1.write(''.join(Output))
    f1.close()
    f1 = open(output_fasta + '.snp.txt', 'w')
    f1.write(''.join(Output_SNP))
    f1.close()
    f1 = open(output_fasta + '.indel.txt', 'w')
    f1.write(''.join(Output_indel))
    f1.close()
    print('done %s mutations in %s'%(num_mut,database_name))
    return output_fasta

def run_mapper(files,files2,database,tempbamoutput):
    max_penalty = 0.1 # default
    if '200000' in database:
        max_penalty = 0.2
    if 'human' in args.i:
        cmds = '/usr/bin/time -v java -Xms850g -Xmx850g -jar %s  --max-penalty %s --num-threads 40 --reference %s --queries %s  --queries %s --out-vcf %s.vcf\n' % (
            latest_mapper, max_penalty, database, files, files2, tempbamoutput)
    else:
        cmds = '/usr/bin/time -v java -Xms10g -Xmx10g -jar %s --max-penalty %s --num-threads 40 --reference %s --queries %s  --queries %s --out-vcf %s.vcf\n' % (
            latest_mapper, max_penalty, database, files, files2, tempbamoutput)
    return cmds

# load homologous regions
homologous_regions = dict()
if args.rmhm != 'False':
    homologous_regions = load_homologous_regions(args.rmhm)

# load database
allgenome = glob.glob('%s/*%s'%(genome_root,genome_name))
for database in allgenome:
    database_name = os.path.split(database)[-1]
    database_file = '%s.fna'%(database)
    try:
        open(database_file,'r')
    except IOError:
        os.system('prodigal -q -i %s -d %s'%(database,database_file))
    Ref_seq, Length, Input_seq, Input_id = loaddatabase(database_file,database)
    # find fastq
    fastq_file = '%s/%s'%(genome_root,database_name.replace(genome_name,fastq_name))
    fastq_file2 = '%s/%s'%(genome_root,database_name.replace(genome_name,fastq_name2))
    # cause SNP
    if len(mut_set) != 0:
        mut_time = len(mut_set)
    while mut_time > 0:
        num_mut = mut_set[mut_time - 1]
        # cause SNP
        if num_mut == 0:
            # corrected genome and control
            mutated_genome = os.path.join(output_dir, 'data/%s.%s.SNP.fasta' % (database_name, num_mut))
            try:
                ftry = open(mutated_genome,'r')
            except IOError:
                os.system('cp %s %s'%(database,mutated_genome))
                f1 = open(mutated_genome + '.snp.txt', 'w')
                f1.close()
        elif cause_SNP:
            # simulate fastq files for mutated strains
            mutated_genome = modelSNPall(Input_seq, Input_id, Length,num_mut,database_name)
        else:
            mutated_genome = os.path.join(output_dir, 'data/%s.%s.SNP.fasta' % (database_name, num_mut))
        mutated_genome_filename = os.path.split(mutated_genome)[-1]
        if mapping_file:
            # call SNPs by time bowtie2
            results = run_vcf_WGS(fastq_file, fastq_file2,
                                  mutated_genome,
                                  os.path.join(output_dir + '/bwa',
                                               mutated_genome_filename + '.bowtie'))
            outputvcf = os.path.join(output_dir + '/merge',
                                               mutated_genome_filename + '.bowtie')
            cmds = results[0]
            cmds += merge_sample(mutated_genome, outputvcf, [results[1]])
            f1 = open(os.path.join(input_script_sub, '%s.bowtie.vcf.sh' % (mutated_genome_filename)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\npy39\n%s' % (''.join(cmds)))
            f1.close()
            # call SNPs by time bwa
            results = run_bwa(fastq_file, fastq_file2,
                                  mutated_genome,
                                  os.path.join(output_dir + '/bwa',
                                               mutated_genome_filename + '.bwa'))
            outputvcf = os.path.join(output_dir + '/merge',
                                     mutated_genome_filename + '.bwa')
            cmds = results[0]
            cmds += merge_sample(mutated_genome, outputvcf, [results[1]])
            f1 = open(os.path.join(input_script_sub, '%s.bwa.vcf.sh' % (mutated_genome_filename)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s' % (''.join(cmds)))
            f1.close()
            # call SNPs by time minimap
            results = run_minimap(fastq_file, fastq_file2,
                                  mutated_genome,
                                  os.path.join(output_dir + '/bwa',
                                               mutated_genome_filename + '.minimap'))
            outputvcf = os.path.join(output_dir + '/merge',
                                     mutated_genome_filename + '.minimap')
            cmds = results[0]
            cmds += merge_sample(mutated_genome, outputvcf, [results[1]])
            f1 = open(os.path.join(input_script_sub, '%s.minimap.vcf.sh' % (mutated_genome_filename)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            f1.close()
            # call SNPs by our mapper
            # run bowtie and mapper on the same node
            cmds = run_mapper(fastq_file, fastq_file2, mutated_genome, os.path.join(output_dir + '/merge',
                                                                                    mutated_genome_filename + '.mapper1'))
            cmds += '/usr/bin/time -v #sh %s\n'%(os.path.join(input_script_sub, '%s.bowtie.vcf.sh' % (mutated_genome_filename)))
            cmds += '/usr/bin/time -v #sh %s\n' % (os.path.join(input_script_sub, '%s.minimap.vcf.sh' % (mutated_genome_filename)))
            cmds += '/usr/bin/time -v #sh %s\n' % (os.path.join(input_script_sub, '%s.bwa.vcf.sh' % (mutated_genome_filename)))
            f1 = open(os.path.join(input_script_sub, '%s.mapper1.vcf.sh' % (mutated_genome_filename)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            f1.close()


        mut_time -= 1

f1 = open(os.path.join(input_script, 'allsnpmodel.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
total_test = 1
if time_evaluation:
    total_test = 10 # run each pipeline 10 times
for m in range(0,total_test):
    for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.mapper1.vcf.sh')):
        os.system('cp %s %s%s'%(sub_scripts,sub_scripts,m))
        if 'human' in args.i:
            f1.write('jobmit %s%s %s%s big\n' % (sub_scripts,m,os.path.split(sub_scripts)[-1],m))
        else:
            f1.write('jobmit %s%s %s%s small\n' % (sub_scripts, m, os.path.split(sub_scripts)[-1], m))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
