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
                      help="all vcf files",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/human/SNP_model_noindel/SNP_compare/SNP_compare_output/',
                      metavar='SNP_compare_output')

################################################## Definition ########################################################
args = parser.parse_args()

alloutput = ['tool\tseq\tcontig\tloci\tcoverage\tSNPs\tNs\n']
allmapperfiles = glob.glob('%s/*.mapper1.vcf'%(args.i))
for mapperresult in allmapperfiles:
    try:
        seqfile = mapperresult.replace('.mapper1.vcf','.fq')
        i = 0
        for lines in open(seqfile,'r'):
            i += 1
            if i ==2:
                seq = lines.split('\n')[0]
        bowtieresult = seqfile.replace('.fq','.bowtie.flt.snp.vcf')
        # parse bowtie result
        bowtiesum = [set(),set(),0,0,0]# contig loci coverage no.SNPs no.Ns
        for lines in open(bowtieresult,'r'):
            lines_set = lines.split('\t')
            # coverage
            bowtiesum[2] += 1
            contig = lines_set[0]
            loci = lines_set[1]
            if contig not in bowtiesum[0]:
                bowtiesum[1].add(loci)
                bowtiesum[0].add(contig)
            if lines_set[4]!= '.' and lines_set[3] != 'N':
                #a SNP or indel
                if len(lines_set[4]) > 1:
                    bowtiesum[3] += abs(len(lines_set[4]) - len(lines_set[3]))
                else:
                    bowtiesum[3] += 1
            if lines_set[3] == 'N':
                # an N ref
                bowtiesum[4] += 1
        # parse mapper result
        mappersum = [set(), set(), 0, 0, 0]  # contig loci coverage no.SNPs no.Ns
        for lines in open(mapperresult, 'r'):
            if not lines.startswith('CHR'):
                lines_set = lines.split('\t')
                # coverage
                mappersum[2] += 1
                contig = lines_set[0]
                loci = lines_set[1]
                if contig not in mappersum[0]:
                    mappersum[1].add(loci)
                    mappersum[0].add(contig)
                if lines_set[3] != '' and lines_set[2] != 'N':
                    # a SNP
                    mappersum[3] += 1
                if lines_set[2] == 'N':
                    # an N ref
                    mappersum[4] += 1
        alloutput.append('bowtie\t%s\t%s\t%s\t%s\t%s\t%s\n'%(seq,';'.join(bowtiesum[0]),';'.join(bowtiesum[1]),bowtiesum[2],bowtiesum[3],bowtiesum[4]
                                                     ))
        alloutput.append('mapper\t%s\t%s\t%s\t%s\t%s\t%s\n' % (seq,
        ';'.join(mappersum[0]), ';'.join(mappersum[1]), mappersum[2], mappersum[3], mappersum[4]
        ))
    except IOError:
        print(mapperresult,seqfile,bowtieresult)

f1 = open('%s/../allsnpcompare.sum.txt'%(args.i),'w')
f1.write(''.join(alloutput))
f1.close()
