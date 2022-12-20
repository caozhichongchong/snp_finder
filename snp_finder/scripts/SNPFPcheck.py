import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime

FP_folder = '/scratch/users/anniz44/genomes/donor_species/SNP_curate/SNP_model_newnew/SNP_model_new/merge/'
BWA_folder = '/scratch/users/anniz44/genomes/donor_species/SNP_curate/SNP_model_newnew/SNP_model_new/bwa/'
genome_folder = '/scratch/users/anniz44/genomes/donor_species/SNP_curate/SNP_model_newnew/SNP_model_new/data'

allFP = glob.glob('%s/*FP'%(FP_folder))
allsum = ['genome\tCHR\tPOS\tlength\n']
for FPfile in allFP:
    genomename = os.path.split(FPfile)[-1].split('.mapper1')[0]
    if '.0.SNP' not in genomename:
        print(datetime.now(), 'processing %s'%(genomename))
        genome = os.path.join(genome_folder,genomename)
        genome_length = dict()
        for record in SeqIO.parse(genome, 'fasta'):
            genome_length.setdefault(str(record.id),len(str(record.seq)))
        print(datetime.now(), 'finish loading CHR length %s' % (genomename))
        genome_set = genomename.split('.')
        genome_set_all = '.'.join(genome_set[:2])
        genome_set[2] = '0'
        genome_set_0 = '.'.join(genome_set)
        FP_control = set()
        for lines in open('%s/%s.mapper1.vcf.final.vcf.FP'%(FP_folder,genome_set_0),'r'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR, POS = lines_set[:2]
            FP_control.add(CHR)
        for lines in open(FPfile,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            CHR,POS = lines_set[:2]
            os.system('echo %s >> %s/allFPreadsam.txt'%(genomename,FP_folder))
            if CHR not in FP_control:
                read = lines_set[-1].split(',')[0]
                read = read.replace('[','').replace(']','')
                allsum.append('%s\t%s\t%s\t%s\n'%(genomename,CHR,POS,genome_length[CHR]))
                os.system('grep \"%s$(printf \'\\t\')%s$(printf \'\\t\')\" %s/%s*FP >> %s/allFPreadsam.txt' % (
                    CHR, POS, FP_folder, genome_set_all, FP_folder
                ))
                os.system('grep \"%s$(printf \'\\t\')%s$(printf \'\\t\')\" %s/%s*mapper1.vcf >> %s/allFPreadsam.txt' % (
                    CHR,POS, FP_folder, genome_set_all, FP_folder
                ))
                os.system('grep %s %s/%s*mapper1.sam >> %s/allFPreadsam.txt'%(
                    read, FP_folder,genome_set_all, FP_folder
                    ))
                os.system('grep %s %s/%s.*.sam >> %s/allFPreadsam.txt' % (
                    read, BWA_folder, genomename, FP_folder
                    ))
                print(datetime.now(), 'finish greping FP read %s %s %s' % (CHR, POS, read))

f1 = open('%s/allFPsum.txt'%(FP_folder),'w')
f1.write(''.join(allsum))
f1.close()
