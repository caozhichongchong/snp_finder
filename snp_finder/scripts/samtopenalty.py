# step 1 modelling SNPs and test 2 methods
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import re
from multiprocessing import Process
from multiprocessing import Queue
from itertools import islice
from datetime import datetime
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-sam",
                      help="sam file",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/SNP_model_penalty/',
                      metavar='mapper.sam')
required.add_argument("-ref",
                      help="reference genome",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/SNP_curate/test_data/penalty_test/',
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
                      type=int, default=10,
                      metavar='10')
optional.add_argument("-genome",
                      help="which genome to check",
                      type=str, default='None',
                      metavar='e.g. genome1')
optional.add_argument("-tool",
                      help="which tool to check",
                      type=str, default='None',
                      metavar='e.g. mapper')
args = parser.parse_args()
filestep = 10000
max_penalty = 1000 # no alignment
os.system('ulimit -n 2048')
thread = int(args.t)
match_penalty = args.match
mismatch_penalty = args.mismatch
gapopen_penalty = args.gapopen
gapextend_penalty = args.gapextend
n_penalty = args.n
penalty_set = {
    'bowtie':[6,5,3,0,1],
    'minimap':[6,4,2,0,1],
    'bwa':[5,6,1,0,1],
    'mapper':[1,2,0.5,0,0.1]
}

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
    #print('total pennlty=',penalty+mismatch*mismatch_penalty + nDNA*n_penalty,'mismatch = ',mismatch,'nDNA=',nDNA)
    final_penalty = (penalty + mismatch*mismatch_penalty + nDNA*n_penalty + (totallen-mismatch-nDNA)*match_penalty)
    #print('final penalty=',final_penalty,totallen,(totallen-mismatch-nDNA)*match_penalty)
    return final_penalty

def ref_read_chopbyCIGAR(CIGARtool):
    read,ref, CIGAR, soft_clip = CIGARtool
    if CIGAR == '*':
        return max_penalty
    CIGARlenset = re.split('S|D|M|I|N|H|P', CIGAR)[:-1]
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
            readset += read[kread:(kread + int(CIGARlen))]
            refset += ref[kref:(kref + int(CIGARlen))]
            kread += int(CIGARlen)
            kref += int(CIGARlen)
            #print('M',penalty,kread,kref,''.join(readset),''.join(refset))
        elif CIGARtype == 'S':
            # clipping in read treat as ambiguous matches
            if soft_clip:
                kread += int(CIGARlen)
                penalty += int(CIGARlen)*n_penalty
                #print('S', penalty,kread, kref)
            else:
                # insertion in read
                kread += int(CIGARlen)
                penalty += gapopen_penalty + (int(CIGARlen) - 1) * gapextend_penalty
                #print('S->I', penalty, kread, kref)
        elif CIGARtype == 'I':
            # insertion in read
            kread += int(CIGARlen)
            penalty += gapopen_penalty + (int(CIGARlen) - 1) * gapextend_penalty
            #print('I', penalty,kread, kref)
        elif CIGARtype == 'D':
            # insertion in ref
            kref += int(CIGARlen)
            penalty += gapopen_penalty + (int(CIGARlen) - 1) * gapextend_penalty
            #print('D', penalty,kread, kref)
        elif CIGARtype == 'N':
            # skipped region from the reference
            kref += int(CIGARlen)
            penalty += int(CIGARlen)*n_penalty
            #print('N', penalty,kread, kref)
        elif CIGARtype == 'H' or CIGARtype == 'P':
            # hard clipping or padding
            pass
        CIGARlocus += 1
    return mismatch_check(readset,refset,penalty)

def load_ref(database):
    Ref = dict()
    for record in SeqIO.parse(database, 'fasta'):
        Ref.setdefault(str(record.id),str(record.seq))
    return Ref

def process_lines(lines_gen,samlist,Ref):
    for lines in lines_gen:
        if not lines.startswith('@'):
            lines_set = lines.split('\t')
            CHR, POS = lines_set[2:4]
            read = lines_set[9]
            readID = lines_set[0] +'\t'+ min(str(Seq(read).reverse_complement()),read)
            CHR = CHR.split(' ')[0]
            if CHR!='*':
                # read matches somewhere
                CIGAR = lines_set[5]
                reflen = len(Ref[CHR])
                readlen = len(read)
                POS = int(POS)
                POSend = POS + 2 * readlen
                soft_clip = False
                if POS <= readlen or POSend >= reflen:
                    # possible soft clip when it's the end of the read
                    soft_clip = True
                ref = Ref[CHR][(POS - 1):min(POSend, reflen)]
                samlist.setdefault(readID,[read,ref, CIGAR,soft_clip])
            else:
                # no alignment
                samlist.setdefault(readID, [read,'*', '*', False])
    return samlist

def read_file(samfile,start,end,queue,Ref):
    samlist = {}
    #print(datetime.now(),'reading %s-%s lines of %s'%(start,end,samfile))
    with open(samfile, 'r+') as f:
        lines_gen = islice(f, start,end)
        samlist = process_lines(lines_gen,samlist,Ref)
    queue.put(samlist)
    #print(datetime.now(), 'finish reading %s-%s lines of %s' % (start, end, samfile))

def compare_readlistsub(allreadsub,alllistfile,queue):
    allreaddiffsub = dict()
    for readID in allreadsub:
        CIGAR0 = alllistfile[0].get(readID, ['', '*', '*', False])  # read,ref, CIGAR,soft_clip
        CIGAR1 = alllistfile[1].get(readID, ['', '*', '*', False])
        CIGAR2 = alllistfile[2].get(readID, ['', '*', '*', False])
        CIGAR3 = alllistfile[3].get(readID, ['', '*', '*', False])
        # print(readID,CIGAR0[2],CIGAR1[2],CIGAR2[2],CIGAR3[2])
        if CIGAR3[2] != CIGAR0[2] or CIGAR3[2] != CIGAR1[2] or CIGAR3[2] != CIGAR2[2]:
            # mapper CIGAR results different
            if CIGAR3[1] != CIGAR0[1] or CIGAR3[1] != CIGAR1[1] or CIGAR3[1] != CIGAR2[1]:
                # mapping ref different
                # print(readID,'stored in allreaddiff')
                allreaddiffsub.setdefault(readID, [CIGAR0, CIGAR1, CIGAR2, CIGAR3])
    queue.put(allreaddiffsub)

def compare_readlist(alllistfile,allreaddiff):
    allreads = list(alllistfile[0].keys()) + list(alllistfile[1].keys()) + list(alllistfile[2].keys()) + list(alllistfile[3].keys())
    allreads = list(set(allreads))
    allreadslen = len(allreads)
    subsetlen = int(allreadslen/thread)
    procs = dict()
    for k in range(0,thread):
        allreadssub = allreads[int(k * subsetlen):int((k + 1) * subsetlen)]
        queue = Queue()
        proc = Process(target=compare_readlistsub,
                       args=(allreadssub, alllistfile, queue))
        procs.setdefault(proc, queue)
        proc.start()
    # for last set
    allreadssub = allreads[int(thread * subsetlen):]
    queue = Queue()
    proc = Process(target=compare_readlistsub,
                   args=(allreadssub, alllistfile, queue))
    procs.setdefault(proc, queue)
    proc.start()
    for proc in procs:
        queue = procs[proc]
        allreaddiffsub = queue.get()
        allreaddiff.update(allreaddiffsub)
        proc.join()
    procs = dict()

def compare_CIGAR_sub(allreadsub,allreaddiff,allreaddiffkeep,queue):
    allreaddiffkeepsub = dict()
    for readID in allreadsub:
        #print('process', readID)
        CIGARset = allreaddiff[readID] # ref, CIGAR,soft_clip
        oldscore = allreaddiffkeep.get(readID, [])
        if oldscore == []:
            allreaddiffkeepsub.setdefault(readID,[ref_read_chopbyCIGAR(x) for x in CIGARset])
        else:
            # only keeping the min penalty for each tool, same read ID same read seq
            newscore = [min(ref_read_chopbyCIGAR(x),y)for x,y in zip(CIGARset,oldscore)]
            allreaddiffkeepsub[readID] = newscore
    queue.put(allreaddiffkeepsub)


def compare_CIGAR(allreaddiff,allreaddiffkeep):
    allreaddifflen = len(allreaddiff)
    allreads = list(allreaddiff.keys())
    subsetlen = int(allreaddifflen / thread)
    procs = dict()
    for k in range(0,thread):
        allreadssub = allreads[int(k * subsetlen):int((k + 1) * subsetlen)]
        queue = Queue()
        proc = Process(target=compare_CIGAR_sub,
                       args=(allreadssub,allreaddiff,allreaddiffkeep,queue))
        procs.setdefault(proc, queue)
        proc.start()
        # for last set
    allreadssub = allreads[int(thread * subsetlen):]
    queue = Queue()
    proc = Process(target=compare_CIGAR_sub,
                   args=(allreadssub, allreaddiff, allreaddiffkeep, queue))
    procs.setdefault(proc, queue)
    proc.start()
    for proc in procs:
        queue = procs[proc]
        allreaddiffkeepsub = queue.get()
        allreaddiffkeep.update(allreaddiffkeepsub)
        proc.join()
    procs = dict()
    return allreaddiffkeep

def parse_samfile(bowtiesam,bwasam,minimapsam,mappersam,Ref):
    allreaddiff = dict()
    allreaddiffkeep = dict()
    allreaddiffkeepnotoutput = dict()
    allpenalty = []
    i = 0
    newlinelen = 0
    # set output
    allpenalty.append('readID\tread\treadlen\tbowtie\tbwa\tminimap\tmapper\n')
    allsamfile = [bowtiesam,bwasam,minimapsam,mappersam]
    while True:
        alllistfile = []
        # read files
        procs = dict()
        for j in range(0,4):
            queue = Queue()
            proc = Process(target=read_file,
                           args=(allsamfile[j], i, i+filestep,queue,Ref))
            procs.setdefault(proc,queue)
            proc.start()
        i += filestep
        for proc in procs: # get items before joining processes
            alllistfile.append(procs[proc].get())
            proc.join()
        procs = dict()
        newlinelen += sum([len(x) for x in alllistfile])
        if newlinelen > 0:
            # compare reads, keep those with CIGAR difference
            compare_readlist(alllistfile, allreaddiff)
            print(datetime.now(),'processed %s lines' % (i))
            print('found %s new reads with diff CIGAR' % (len(allreaddiff)))
            if len(allreaddiff) > 0:
                # process CIGAR for different tools
                allreaddiffkeep = compare_CIGAR(allreaddiff,allreaddiffkeep)
            print('in total %s reads with diff CIGAR' % (len(allreaddiffkeep)))
            # output to files
            for readID in allreaddiffkeep:
                newscore = allreaddiffkeep[readID]
                if len(set(newscore)) > 1:
                    # different scores by different tools
                    if max_penalty not in newscore:
                        # read aligned by all tools compared
                        newscore = ['%.2f' % (x) for x in newscore]
                        allpenalty.append('%s\t%s\t%s\n' % (readID, len(readID.split('\t')[-1]), '\t'.join(newscore)))
                    else:
                        # read by some tools not imported yet
                        allreaddiffkeepnotoutput.setdefault(readID, newscore)
            f1.write(''.join(allpenalty))
            # empty lists
            allreaddiff = dict()
            newlinelen = 0
            allpenalty = []
            allreaddiffkeep = allreaddiffkeepnotoutput
            allreaddiffkeepnotoutput = dict()
        else:
            if len(allreaddiff) > 0:
                # process CIGAR for different tools, output all leftover reads
                allreaddiffkeep = compare_CIGAR(allreaddiff, allreaddiffkeep)
            print(datetime.now(),'end of the sam files: %s lines'%(i-filestep))
            break
    for readID in allreaddiffkeep:
        newscore = allreaddiffkeep[readID]
        if len(set(newscore)) > 1:
            # different scores by different tools
            newscore = ['%.2f'%(x) for x in newscore]
            allpenalty.append('%s\t%s\t%s\n'%(readID,len(readID.split('\t')[-1]),'\t'.join(newscore)))
    f1.write(''.join(allpenalty))

def sort_file(samfile):
    try:
        f2 = open(samfile + '.sort')
    except FileNotFoundError:
        print(datetime.now(),'sort file %s'%(samfile))
        os.system('sort %s > %s.sort'%(samfile,samfile))
    return samfile + '.sort'

################################################### MAIN ########################################################
Toolset = ['bowtie','bwa','minimap','mapper']
if args.tool != 'None':
    Toolset = [args.tool]
# load database
if args.genome != 'None':
    allref = glob.glob('%s/%s*.fasta' % (args.ref,args.genome))
else:
    allref = glob.glob('%s/*.fasta'%(args.ref))
if args.genome != 'None':
    allfastq = glob.glob('%s/%s*_1.fastq' % (args.ref, args.genome))
else:
    allfastq = glob.glob('%s/*_1.fastq'%(args.ref))
for fastq in allfastq:
    samplename = os.path.split(fastq)[-1]
    for refgenome in allref:
        Ref = load_ref(refgenome)
        # set input
        genomename = os.path.split(refgenome)[-1]
        print(datetime.now(), 'process ref %s' % (genomename))
        for tool in Toolset:
            try:
                f2 = open('%s/%s/%s_%s.mapper1.sam' % (args.sam, tool,samplename,genomename),'r')
                f2.close()
                print(datetime.now(), 'process sams for %s %s' % (samplename, genomename))
                f1 = open('%s/%s/%s_%s.%s.samcompare.txt' % (args.sam, tool, samplename,genomename, tool), 'w')
                mismatch_penalty, gapopen_penalty, gapextend_penalty, match_penalty, n_penalty = penalty_set[tool]
                # load all sam files
                bowtiesam = sort_file('%s/%s/%s_%s.bowtie.sam' % (args.sam, tool, samplename,genomename))
                bwasam = sort_file('%s/%s/%s_%s.bwa.sam' % (args.sam, tool, samplename,genomename))
                minimapsam = sort_file('%s/%s/%s_%s.minimap.sam' % (args.sam, tool, samplename,genomename))
                mappersam = sort_file('%s/%s/%s_%s.mapper1.sam' % (args.sam, tool, samplename,genomename))
                parse_samfile(bowtiesam, bwasam, minimapsam, mappersam, Ref)
                f1.close()
            except IOError:
                print('no output for %s %s'%(samplename,genomename))


################################################### END ########################################################
