import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

fasta_file='all.selected.gene.faa'
true_PE = 'all.species.lineage.dmrca.txt'

def load_PE(annnofile):
    CL_set = set()
    for lines in open(annnofile,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        if lines_set[-1] != 'False' and lines_set[-3] != 'False':
            # FDR True + significance
            CL = lines_set[0].replace('clustercluster','CL')
            CL_set.add(CL)
    return CL_set

def extract_PE(fasta,CL_set):
    PE_out = []
    for record in SeqIO.parse(fasta, 'fasta'):
        record_id = str(record.id)
        CL = record_id.split('__')[0]
        if CL in CL_set:
            PE_out.append('>%s\n%s\n'%(record_id,str(record.seq)))
    foutput = open(fasta + '.FDR.faa', 'w')
    foutput.write(''.join(PE_out))
    foutput.close()

CL_set=load_PE(true_PE)
print(CL_set)
extract_PE(fasta_file,CL_set)
