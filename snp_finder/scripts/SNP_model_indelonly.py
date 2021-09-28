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
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 40)",
                      metavar="1 or more", action='store', default=40, type=int)

################################################## Definition ########################################################
args = parser.parse_args()
input_script = args.s
genome_root = args.i
output_dir = args.o + '/indel_model'
genome_name = args.fa
fastq_name=args.fq
fastq_name2=args.fq.replace('1','2')
input_script_sub = '%s/indel_model'%(input_script)

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
mut_set = range(2,10)
indel_time = 2 # how many indels in a genome
cause_SNP = False
mapping_file = True
indel_orf = [-10,-7,-4, 4, 7, 10]
indel_nonorf = list(range(2,11))
indel_nonorf.extend(list(range(-10,-1)))
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

def loaddatabase(database_aa,database):
    # load database seq
    Length = []
    Mapping_loci = dict()
    reference_database = os.path.split(database_aa)[-1]
    print('reference database_aa set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    Input_seq = dict()
    Input_id = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        seq_length = len(record_seq)
        if seq_length >= 5000:
            Input_seq.setdefault(record_id, record_seq)
            Input_id.append(record_id)
            Length.append(seq_length)
    for record in SeqIO.parse(database_aa, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        Ref_seq.setdefault(record_id, record_seq)
        if float(description[3]) == -1.0:  # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Length,Mapping_loci,Reverse,Input_seq,Input_id]

def modelindel(seq,Chr,indel_set):
    SNP_output = []
    indel_set.sort()
    record_indel = dict()
    for position in indel_set:
        total_length = len(seq)
        if position < total_length:
            REF = seq[position]
            gene_info = contig_to_gene(Chr, position)
            temp_ALT = ['A', 'T', 'G', 'C']
            if gene_info != []:
                # a gene
                # indel = + 3*n
                indel_size = random.choices(indel_orf, k=1)[0]
            else:
                # not a gene
                indel_size = random.choices(indel_nonorf, k=1)[0]
            if indel_size > 0:# insertion on ref
                ALT = random.choices(temp_ALT, k=indel_size)
                seq = seq[:position] + ALT + seq[position+1:]
                temp_line = [Chr, str(position + 1), REF, ''.join(ALT)]
                record_indel.setdefault(position, [REF,''.join(ALT),indel_size])
            else:# deletion on ref
                REF_after = ''.join(seq[position:(position-indel_size)])
                REF = ''.join(seq[(position+indel_size):position])
                del seq[(position+indel_size):position]
                temp_line = [Chr, str(position + 1), REF, '-'*(-indel_size)]
                record_indel.setdefault(position, [REF_after, '-'*(-indel_size),indel_size])
            SNP_output.append('\t'.join(temp_line) + '\n')
        else:
            print('position %s out of the reference %s'%(position,total_length))
    return [seq, SNP_output]

def modelSNP(seq,Chr,num_indel_chr):
    total_length = len(seq)
    # indel modelling
    indel_output = []
    seq = list(seq)
    if num_indel_chr > 0:
        candidate_position = [i for i in range(0, total_length) if seq[i] not in ['-', 'N']]
        indel_set = random.sample(candidate_position, k=num_indel_chr)
        seq, indel_output = modelindel(seq, Chr, indel_set)
    return [''.join(seq),indel_output]

def run_vcf_WGS(files,files2,database,tempbamoutput):
    # generate code
    cmds = 'time bowtie2-build %s %s\n'%(database,database)
    try:
        f1 = open('%s.sorted.bam' % (tempbamoutput),'r')
    except IOError:
        cmds += 'time bowtie2 --threads %s -x %s -1 %s -2 %s |time samtools view -@ %s -S -b >%s.bam\ntime samtools sort -@ %s %s.bam -o %s.sorted.bam\ntime samtools index -@ %s %s.sorted.bam\n' % (
            min(40, args.t), database, files, files2,  min(40, args.t),
            tempbamoutput,  min(40, args.t), tempbamoutput, tempbamoutput, min(40, args.t),
            tempbamoutput)
        cmds += 'rm -rf %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def merge_sample(database,vcfoutput,allsam):
    cmds = ''
    try:
        f1 = open('%s.raw.vcf' % (vcfoutput), 'r')
    except FileNotFoundError:
        cmds += 'time bcftools mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | time bcftools call -c -Ov --threads %s > %s.raw.vcf\n' % (
             min(40, args.t), database,
            ' '.join(allsam), min(40, args.t), vcfoutput)
    try:
        f1 = open('%s.flt.snp.vcf' % (vcfoutput))
    except FileNotFoundError:
        cmds += 'time bcftools view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
            vcfoutput, vcfoutput)
    return cmds

def run_minimap(files,files2,database,tempbamoutput):
    #os.system('time minimap2 -d %s.mmi %s \n' % (database, database))
    cmds = 'time minimap2 -ax sr -N 1 -p 0.99 -t %s %s.mmi %s %s >%s.sam\nsource deactivate py37\ntime samtools view -@ %s -S -b %s.sam >%s.bam\ntime samtools sort -@ %s %s.bam -o %s.sorted.bam\ntime samtools index -@ %s %s.sorted.bam\n' % (
        min(40, args.t), database, files, files2, tempbamoutput, min(40, args.t),tempbamoutput,
            tempbamoutput,  min(40, args.t), tempbamoutput, tempbamoutput, min(40, args.t),
            tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def modelSNPall(Input_seq, Input_id, Length,num_mut,database_name):
    Output = []
    Output_indel = []
    for chr in Input_id:
        # change indel
        newseq, newoutputindel = modelSNP(Input_seq[chr], chr, indel_time)
        Output_indel += newoutputindel
        Output.append('>%s\n%s\n' % (chr,
                                     newseq))
    # output mutated genome
    output_fasta = os.path.join(output_dir, 'data/%s.%s.SNP.fasta' % (database_name,num_mut))
    f1 = open(output_fasta, 'w')
    f1.write(''.join(Output))
    f1.close()
    f1 = open(output_fasta + '.indel.txt', 'w')
    f1.write(''.join(Output_indel))
    f1.close()
    print('done %s mutations in %s'%(num_mut,database_name))
    return output_fasta

def run_mapper(files,database,tempbamoutput):
    cmds = 'time java -jar %s/mapper1.5.jar --reference %s --queries %s --out-vcf %s.vcf\n' % (args.s,database, files, tempbamoutput)
    return cmds

# load database
allgenome = glob.glob('%s/indeltest*%s'%(genome_root,genome_name))
for database in allgenome:
    database_name = os.path.split(database)[-1]
    database_file = '%s.fna'%(database)
    try:
        open(database_file,'r')
    except IOError:
        os.system('prodigal -q -i %s -d %s'%(database,database_file))
    Ref_seq, Length, Mapping_loci, Reverse, Input_seq, Input_id = loaddatabase(database_file,database)
    # find fastq
    fastq_file = '%s/%s'%(genome_root,database_name.replace(genome_name,fastq_name))
    fastq_file2 = '%s/%s'%(genome_root,database_name.replace(genome_name,fastq_name2))
    # cause SNP
    if len(mut_set) != 0:
        mut_time = len(mut_set)
    while mut_time > 0:
        num_mut = mut_set[mut_time - 1]
        # cause SNP
        if cause_SNP:
            # simulate fastq files for mutated strains
            mutated_genome = modelSNPall(Input_seq, Input_id, Length,num_mut,database_name)
        else:
            mutated_genome = os.path.join(output_dir, 'data/%s.%s.SNP.fasta' % (database_name, num_mut))
        mutated_genome_filename = os.path.split(mutated_genome)[-1]
        if mapping_file:
            # call SNPs by our mapper
            cmds = run_mapper(fastq_file,mutated_genome,os.path.join(output_dir + '/merge',
                                               mutated_genome_filename + '.mapper1'))
            f1 = open(os.path.join(input_script_sub, '%s.mapper1.vcf.sh' % (mutated_genome_filename)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (''.join(cmds)))
            f1.close()

        mut_time -= 1

f1 = open(os.path.join(input_script, 'allsnpmodel.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.vcf.sh')):
    f1.write('jobmit %s %s small1\n' % (sub_scripts,os.path.split(sub_scripts)[-1]))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
