# start
# step 3 curate genome
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics
# set up path
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of folders of WGS of each species",
                      type=str, default='.',
                      metavar='input/')
required.add_argument("-fq",
                      help="file extension of WGS fastq #1 files",
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
# optional search parameters
optional.add_argument('-t',
                      help="Optional: set the thread number assigned for running XXX (default 1)",
                      metavar="1 or more", action='store', default=1, type=int)
optional.add_argument('-job',
                      help="Optional: command to submit jobs",
                      metavar="nohup or customized",
                      action='store', default='jobmit', type=str)
# requirement for software calling
optional.add_argument('-bw', '--bowtie',
                          help="Optional: complete path to bowtie if not in PATH",
                          metavar="/usr/local/bin/bowtie",
                          action='store', default='bowtie', type=str)
optional.add_argument('-sp', '--spades',
                          help="Optional: complete path to spades if not in PATH",
                          metavar="/usr/local/bin/spades",
                          action='store', default='spades', type=str)
optional.add_argument('-pro', '--prodigal',
                      help="Optional: complete path to prodigal if not in PATH, None for no prodigal (default)",
                      metavar="/usr/local/bin/prodigal",
                      action='store', default='None', type=str)
optional.add_argument('-bcf', '--bcftools',
                      help="Optional: complete path to bcftools if not in PATH",
                      metavar="/usr/local/bin/bcftools",
                      action='store', default='bcftools', type=str)
optional.add_argument('-sam', '--samtools',
                      help="Optional: complete path to bwa if not in PATH",
                      metavar="/usr/local/bin/samtools",
                      action='store', default='samtools', type=str)
optional.add_argument('-mini', '--minimap2',
                      help="Optional: complete path to minimap2 if not in PATH",
                      metavar="/usr/local/bin/minimap2",
                      action='store', default='minimap2', type=str)

################################################## Definition ########################################################
args = parser.parse_args()

input_script = args.i
genome_root = args.i + '/*'
genome_dir = glob.glob(args.i + '/*')
output_dir = '/scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round1/SNP_curate'
fastq_name = args.o + '/vcf_round1/SNP_curate'
samplename = args.fq
genome_name = '.fasta'
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

# correct genomes
cmds = ''
Output_report = []
for genome in assembly_quality:
    assembly_quality_genome = assembly_quality.get(genome, [])
    if assembly_quality_genome != []:
        database = glob.glob(os.path.join(genome_root, genome.split(fastq_name)[0].replace(samplename, genome_name)))
        if database != []:
            database = database[0]
            Output_report += correct_genome(database, assembly_quality_genome, assembly_quality_chr, report_curate)
            # clean up old genomes
            cmds += 'rm -rf %s\n' % (database)
        else:
            print('missing input fasta for %s' % (genome))
    else:
        os.system('cp %s %s.corrected.fasta'%(genome,genome))

f1 = open(os.path.join(input_script, 'cleanoldgenome.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n' + 'rm -rf %s\n'%(output_dir))
f1.write(''.join(cmds))
f1.close()

f1 = open(os.path.join(input_script, 'SNP_currate1.assembly.sum.curated'), 'w')
f1.write(''.join(Output_report))
f1.close()

################################################### END ########################################################
