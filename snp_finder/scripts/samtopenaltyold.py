# step 1 modelling SNPs and test 2 methods
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import re
from multiprocessing import Process
from multiprocessing import Manager
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-sam",
                      help="sam file",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/SNP_model/SNP_model/bwaresults/',
                      metavar='mapper.sam')
required.add_argument("-ref",
                      help="reference genome",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/SNP_model/SNP_model/data/',
                      metavar='reference.fasta')
# optional output setup
optional.add_argument("-match",
                      help="penalty for matches",
                      type=float, default=0,
                      metavar='0')
optional.add_argument("-mismatch",
                      help="penalty for mismatches",
                      type=float, default=1,
                      metavar='-1')
optional.add_argument("-gapopen",
                      help="penalty for opening gaps",
                      type=float, default=2,
                      metavar='-2')
optional.add_argument("-gapextend",
                      help="penalty for extending gaps",
                      type=float, default=0.5,
                      metavar='-0.5')
optional.add_argument("-n",
                      help="penalty for ambiguous matches",
                      type=float, default=0.1,
                      metavar='0.9')
optional.add_argument("-t",
                      help="num of threads",
                      type=float, default=10,
                      metavar='10')
args = parser.parse_args()
################################################## Definition ########################################################
def mismatch_check(readset,refset,penalty):
    totallen = min(len(readset),len(refset))
    mismatch = 0
    nDNA = 0
    for i in range(0,totallen):
        if refset[i] == 'N' or readset[i]=='N':
            nDNA += 1
        elif readset[i] != refset[i]:
            mismatch += 1
    #print('total pennlty=',penalty+mismatch*args.mismatch + nDNA*args.n,'mismatch = ',mismatch,'nDNA=',nDNA)
    final_penalty = (penalty + mismatch*args.mismatch + nDNA*args.n + (totallen-mismatch-nDNA)*args.match)
    #print('final penalty=',final_penalty,totallen,(totallen-mismatch-nDNA)*args.match)
    return final_penalty

def ref_read_chopbyCIGAR(ReadID,read,ref,CIGAR,soft_clip,allpenalty):
    CIGARlenset = re.split('S|D|M|I', CIGAR)[:-1]
    i = 0
    kread = 0
    kref = 0
    readset = []
    refset = []
    penalty = 0
    #print('start',read,ref,CIGAR,CIGARlenset)
    CIGARlocus = 0
    for CIGARlen in CIGARlenset:
        i += len(CIGARlen)
        CIGARtype = CIGAR[i]
        #print('process CIGAR',CIGARtype,CIGARlen)
        i += 1
        if CIGARtype == 'M':
            # matches
            readset += list(read[kread:(kread + int(CIGARlen))])
            refset += list(ref[kref:(kref + int(CIGARlen))])
            kread += int(CIGARlen)
            kref += int(CIGARlen)
            #print('M',penalty,kread,kref,''.join(readset),''.join(refset))
        elif CIGARtype == 'S':
            # clipping in read treat as ambiguous matches
            if soft_clip:
                kread += int(CIGARlen)
                penalty += int(CIGARlen)*args.n
                #print('S', penalty,kread, kref)
            else:
                # matches
                readset += list(read[kread:(kread + int(CIGARlen))])
                refset += list(ref[kref:(kref + int(CIGARlen))])
                kread += int(CIGARlen)
                kref += int(CIGARlen)
                #print('S->M', penalty, kread, kref, ''.join(readset), ''.join(refset))
        elif CIGARtype == 'I':
            # insertion in read
            kread += int(CIGARlen)
            penalty += args.gapopen + (int(CIGARlen) - 1) * args.gapextend
            #print('I', penalty,kread, kref)
        elif CIGARtype == 'D':
            # insertion in ref
            kref += int(CIGARlen)
            penalty += args.gapopen + (int(CIGARlen) - 1) * args.gapextend
            #print('D', penalty,kread, kref)
        CIGARlocus += 1
    allpenalty.setdefault(ReadID,mismatch_check(readset,refset,penalty))

def load_ref(database):
    Ref = dict()
    for record in SeqIO.parse(database, 'fasta'):
        Ref.setdefault(str(record.id),str(record.seq))
    return Ref

def parse_samfile(samfile,Ref):
    manager = Manager()
    allpenalty = manager.dict()
    alllist = []
    procs = []
    i = 0
    for lines in open(samfile,'r'):
        if not lines.startswith('@'):
            lines_set = lines.split('\t')
            CHR, POS = lines_set[2:4]
            if CHR!='*':
                # read matches somewhere
                ReadID = lines_set[0]
                CIGAR = lines_set[5]
                read = lines_set[9]
                reflen = len(Ref[CHR])
                readlen = len(read)
                POS = int(POS)
                POSend = POS + 2 * readlen
                soft_clip = False
                if POS <= readlen or POSend >= reflen:
                    # possible soft clip when it's the end of the read
                    soft_clip = True
                ref = Ref[CHR][(POS - 1):min(POSend, reflen)]
                i += 1
                alllist.append([ReadID,read, ref, CIGAR,soft_clip])
                if i % args.t == 0:
                    # start parallel
                    for seqlist in alllist:
                        proc = Process(target=ref_read_chopbyCIGAR, args=(seqlist[0],seqlist[1],seqlist[2],seqlist[3],seqlist[4],allpenalty,))
                        procs.append(proc)
                        proc.start()
                    # finish parallel
                    for proc in procs:
                        proc.join()
                        alllist = []
                        procs = []
                if i%100000 == 0:
                    print('process %s lines'%(i))
    # process the last batch
    # start parallel
    #print(alllist)
    for seqlist in alllist:
        proc = Process(target=ref_read_chopbyCIGAR, args=(seqlist[0],seqlist[1],seqlist[2],seqlist[3],seqlist[4],allpenalty,))
        procs.append(proc)
        proc.start()
    # finish parallel
    for proc in procs:
        proc.join()
    return allpenalty

################################################### MAIN ########################################################
Toolset = {
    'bowtie':0,
    'bwa':1,
    'minimap':2,
    'mapper1':3
}
# load database
allref = glob.glob('%s/am_PaDi_g0001.fasta.1000.SNP.fasta'%(args.ref))
for refgenome in allref:
    Ref = load_ref(refgenome)
    # set output
    genomename = os.path.split(refgenome)[-1]
    print('process ref %s'%(genomename))
    # load all sam files
    allsamfile = glob.glob('%s/%s*.sam'%(args.sam,genomename))
    allsamfile.sort()
    allpenaltyall = dict()
    for samfile in allsamfile:
        print('process',samfile)
        samfilename = os.path.split(samfile)[-1].split('.')[-2]
        toolorder = Toolset[samfilename]
        allpenalty = parse_samfile(samfile,Ref)
        for ReadID in allpenalty:
            allpenaltyall.setdefault(ReadID,[0,0,0,0])
            allpenaltyall[ReadID][toolorder] = allpenalty[ReadID]
    print(allpenaltyall)
    alloutput = ['readID\tbowtie\tbwa\tminimap\tmapper\n']
    for ReadID in allpenaltyall:
        alloutput.append('%s\t%s\n'%(ReadID,'\t'.join([str(x) for x in allpenaltyall[ReadID]])))
    f1 = open('%s/%s_samresult.txt'%(args.sam,genomename),
              'w')
    f1.write(''.join(alloutput))
    f1.close()
################################################### END ########################################################
