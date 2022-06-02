import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import random
import argparse
import subprocess
import sys
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="all vcf files",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/SNP_model/SNP_model/SNP_compare/SNP_compare_output/',
                      metavar='SNP_compare_output')

################################################## Definition ########################################################
args = parser.parse_args()
snp_penalty = 1
new_indel_penalty = 2
extend_indel_penalty = 0.5
nomatch_penalty = 0.1
ambiguous_penalty = 0.1

def find_ref(reference,contig,locus,length,tempoutput):
    out = subprocess.check_output('source ~/.bashrc\nsqe %s \'%s\' %s %s'%(reference,contig,locus,length), shell=True).decode(sys.stdout.encoding)
    f1 = open(tempoutput,'w')
    f1.write('>1\n%s'%(out))
    f1.close()
    return [out,len(out)]

def call_blast(reffile,seqfile,outputfile):
    shline = '#!/bin/bash\nsource ~/.bashrc\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n'
    shline += 'makeblastdb -in %s -dbtype nucl\n'%(reffile)
    shline += 'blastn -subject %s -query %s.fa -out %s -outfmt 6\n' % (reffile,seqfile,outputfile)
    f1 = open(outputfile + '.sh', 'w')
    f1.write(shline)
    f1.close()
    os.system('sh %s'%(outputfile + '.sh'))

def process_blast(blastresult,lenseq,ref):
    try:
        for lines in open(blastresult, 'r'):
            linesset = lines.split('\t')
            ID, Length,Mismatch,Gaps,qstart,qend,dstart,dend = linesset[2:10]
            dstart = int(dstart)
            dend = int(dend)
            Gaps = int(Gaps)
            Mismatch = int(Mismatch)
            totalN = ref[min(dend,dstart)-1:max(dend,dstart)].count('N')
            # snp penalty + nomatch penalty + ambiguous penalty
            total_panelty = Mismatch * snp_penalty + (lenseq - int(Length)) * nomatch_penalty + totalN * ambiguous_penalty
            if Gaps > 0:
                # indel penalty
                total_panelty += new_indel_penalty + (Gaps-1) * extend_indel_penalty
            return ['%s\t%s\t%s\t%s\t%s\t%s\t%s'%(totalN, ID, Length,Mismatch,Gaps,qstart,qend),total_panelty]
    except IndexError:
        rint('no lines in %s' % (blastresult))
    return ['0\t0\t0\t0\t0\t0',0]

def loadfq(seqfile):
    i = 0
    for lines in open(seqfile, 'r'):
        i += 1
        if i == 2:
            seq = lines.strip()
            writefq(seqfile, seq)
            return seq
    return ''

def writefq(seqfile,seq):
    f1 = open(seqfile + '.fa','w')
    f1.write('>1\n%s\n'%(seq))
    f1.close()

def process_vcfresult(vcffile,reference,seqfile,seq,tool):
    CHR = ''
    POS = []
    reffile = vcffile + '.ref.fasta'
    blastoutput = vcffile + '.blast.txt'
    try:
        for lines in open(vcffile,'r'):
            if not lines.startswith('CHR'):
                lines_set = lines.split('\t')
                CHR = lines_set[0]
                POS.append(int(lines_set[1]))
        ref,lenref = find_ref(reference, CHR, POS[0], POS[-1] - POS[0] + 1, reffile)
        call_blast(reffile, seqfile, blastoutput)
        blastoutputresult, total_panelty = process_blast(blastoutput,len(seq),ref)
        return [total_panelty,genomename,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (tool, genomename, CHR, POS[0], total_panelty, lenref, blastoutputresult,seq)]
    except IndexError:
        print('no lines in %s'%(vcffile))
    return [100,genomename,'%s\t%s\tNone\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t%s\n'%(tool,genomename,seq)]

alloutput = ['tool\tgenomename\tcontig\tloci\ttotal_penalty\tlenref\tNo.N_ref\tID\tLength\tMismatch\tGaps\tQstart\tQend\tseq\n']
allfq = glob.glob('%s/*.fq'%(args.i))
alloutputtool = dict()
for seqfile in allfq:
    seq = loadfq(seqfile)
    if seq != '':
        genomename = os.path.split(seqfile)[-1].split('.seq')[0]
        reference = '%s/../../data/%s'%(args.i,genomename)
        # process bowtie results
        vcffile = seqfile.replace('.fq','.bowtie.flt.snp.vcf')
        total_panelty, genomename, process_vcfresultline = process_vcfresult(vcffile,reference,seqfile,seq,'bowtie')
        alloutput.append(process_vcfresultline)
        genomename_seq = '%s\t%s'%(genomename,seq)
        alloutputtool.setdefault(genomename_seq,[total_panelty,0,0,0])
        # process minimap results
        vcffile = seqfile.replace('.fq', '.minimap.flt.snp.vcf')
        total_panelty, genomename, process_vcfresultline = process_vcfresult(vcffile, reference, seqfile, seq, 'minimap')
        alloutput.append(process_vcfresultline)
        alloutputtool[genomename_seq][1]=total_panelty
        # process bwa results
        if 'human' not in args.i:
            vcffile = seqfile.replace('.fq', '.bwa.flt.snp.vcf')
            total_panelty, genomename, process_vcfresultline = process_vcfresult(vcffile, reference, seqfile, seq, 'bwa')
            alloutput.append(process_vcfresultline)
            alloutputtool[genomename_seq][2] = total_panelty
        # process mapper results
        vcffile = seqfile.replace('.fq', '.mapper1.vcf')
        total_panelty, genomename, process_vcfresultline = process_vcfresult(vcffile, reference, seqfile, seq, 'mapper1')
        alloutput.append(process_vcfresultline)
        alloutputtool[genomename_seq][3] = total_panelty
    else:
        print('warning: empty seqfile %s'%(seqfile))

f1 = open('%s/../allsnpcompare.sum.txt'%(args.i),'w')
f1.write(''.join(alloutput))
f1.close()

alloutput = ['genomename\tseq\tbowtie\tminimap\tbwa\tmapper\n']
for genomename_seq in alloutputtool:
    total_panelty_set = [str(x) for x in alloutputtool[genomename_seq]]
    alloutput.append('%s\t%s\n'%(genomename_seq,'\t'.join(total_panelty_set)))

f1 = open('%s/../allsnpcompare.sumtool.txt'%(args.i),'w')
f1.write(''.join(alloutput))
f1.close()