import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

#input_genome = '/scratch/users/anniz44/genomes/donor_species/SNP_curate/am_BA_g0003.fasta'
#input_genome = '/scratch/users/anniz44/genomes/donor_species/SNP_curate/human/SRR2842672.fasta'
input_genome_dir = '/scratch/users/anniz44/genomes/donor_species/SNP_curate/test_data/new'
allintput_genome = glob.glob('%s/*.fasta'%(input_genome_dir))
kmer_size = 100
interval = kmer_size * 2 # generates kmers for interval bp then stops for interval bp
intervalmerge = interval * 2 # merge neighbouring kmers matching the similar positions within intervalmerge bp
total_match_cutoff = 5 # output POS with total_match_cutoff k-mers matching

def kmer_seq(CHR, sequence,kmer_set):
    print('start spliting %s into kmers'%(CHR))
    for i in range(0,len(sequence) - kmer_size):
        if int(i/interval) % 2 == 1:
            # generate kmers
            kmer = sequence[i:(i+kmer_size)]
            kmer_set.setdefault(kmer,{})
            kmer_set[kmer].setdefault(CHR,[])
            kmer_set[kmer][CHR].append(i)
    return kmer_set

def kmer_genome(input_genome):
    kmer_set = dict()
    for record in SeqIO.parse(input_genome, 'fasta'):
        kmer_set = kmer_seq(str(record.id), str(record.seq),kmer_set)
    return kmer_set

def comparePOS(POSstart,POSend,POS):
    if POS >= POSstart - intervalmerge and POS <= POSend + intervalmerge:
        return [True,min(POSstart,POS), max(POS,POSend)]
    else:
        return [False]

def compare_kmerssimple(kmer_set):
    alloutput = dict()
    i = 0
    k = 0
    allkmerlen = len(kmer_set)
    for kmer in kmer_set:
        contigmatch = kmer_set[kmer]
        totalmatch = sum([len(contigmatch[x]) for x in contigmatch])
        if totalmatch >= total_match_cutoff:
            # at least X kmer matches
            for CHR in contigmatch:
                for POS in contigmatch[CHR]:
                    alloutput.setdefault(i,[CHR,POS,totalmatch])
                    i += 1
        k += 1
        if (k%5000) == 0:
            print('processed %s of %s kmers'%(k,allkmerlen))
    print('finished processing kmer matching')
    # output replicate regions
    print('Outputing kmer matching')
    # sort CHR POS
    alloutput = pd.DataFrame.from_dict(alloutput,orient='index',columns=['CHR','POS','No.matches'])
    alloutput = alloutput.sort_values(['CHR','POS'])
    alloutput.to_csv(input_genome.replace('.fasta','.duplicate.txt'), sep='\t', index=False)

for input_genome in allintput_genome:
    kmer_set = kmer_genome(input_genome)
    print('finished splitting genome %s into kmers'%(input_genome))
    compare_kmerssimple(kmer_set)