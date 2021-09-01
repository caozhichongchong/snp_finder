# prokka of all ref genomes
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
allreference_list = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/all.reference.list'
#all_genomes = '/scratch/users/anniz44/genomes/donor_species//vcf_round2/co-assembly/*/*.noHM.fasta'
#allprokka = '/scratch/users/anniz44/genomes/donor_species//vcf_round2/co-assembly/*/*.noHM.fasta.faa.prokka.txt'
denovo_mut = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.denovo.gene.faa'
HS_mut = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.selected.gene.faa'

def load_gene(gene_file):
    gene_set = dict()
    change_name = set()
    for record in SeqIO.parse(gene_file, 'fasta'):
        record_id = str(record.id)
        record_des = str(record.description)
        contig_id = '%s__C_%s' % (donor_species_new, record_id.split('_')[1])
        record_id = '%s__C_%s_G_%s' % (donor_species_new,record_id.split('_')[1], record_id.split('_')[-1])
        start = record_des.split('#')[1].replace(' ','')
        stop = record_des.split('#')[2].replace(' ','')
        gene_set.setdefault('%s_%s_%s'%(contig_id,start,stop),record_id)
        change_name.add('%s\t%s\n'%(record_id,len(str(record.seq))))
    f1 = open(gene_file + '.length.txt', 'w')
    f1.write(''.join(list(change_name)))
    f1.close()
    return gene_set

def contig_loci(contig):
    return [int(contig.split('_')[-3]),int(contig.split('_')[-2]),int(contig.split('_')[-1])]

def compare_contig(contig1,gene_set):
    contig1 = contig_loci(contig1)
    for contig in gene_set:
        contig2 = contig_loci(contig)
        if contig2[0] == contig1[0] and ((contig1[1]<=contig2[2] and contig1[1]>=contig2[1]) or \
                (contig1[2]>=contig2[1] and contig1[2]<=contig2[2])):
            return [True,contig]
    return [False,'']

def process_prokka(gff_file,gene_set,genome):
    # need gene COG product
    allanno = []
    for lines in open(gff_file,'r'):
        if lines.startswith(donor_species_new):
            lines_set = lines.split('\t')
            contig, notused, notsued2, start, stop, notused3, notused4, notused5, anno = lines_set
            anno = anno.split('\n')[0]
            contig = '%s_%s_%s' % (contig, start, stop)
            compare_contig_result = compare_contig(contig, gene_set)
            if compare_contig_result[0]:
                genename = gene_set[compare_contig_result[1]]
                try:
                    gene = anno.split('gene=')[1].split(';')[0]
                except IndexError:
                    gene = ''
                try:
                    product = anno.split('product=')[1].split(';')[0]
                except IndexError:
                    product = ''
                try:
                    COG = anno.split('db_xref=')[1].split(';')[0]
                except IndexError:
                    COG = ''
                allanno.append('%s\t%s\t%s\t%s\n' % (genename, gene, COG, product))
            elif 'RNA' not in lines and 'repeat_region' not in lines and 'hypothetical protein' not in lines:
                print('no match',contig,lines)
    f1 = open(genome + '.faa.prokka.txt', 'w')
    f1.write('genename\tgene\tCOG\tproduct\n')
    f1.write(''.join(allanno))
    f1.close()

def load_mut(mut_fasta):
    gene_set = dict()
    for record in SeqIO.parse(mut_fasta, 'fasta'):
        record_id = str(record.id)
        newrecord_id = '%s__%s'%(record_id.split('.donor.')[0],record_id.split('__')[1])
        gene_set.setdefault(newrecord_id,[record_id,3*len(str(record.seq))])# DNA level
    return gene_set

def load_prokka(prokka):
    for lines in open(prokka,'r'):
        lines_set = lines.split('\t')
        genename = lines_set[0]
        if genename in allmut:
            alloutput.add(genename)
            allmut_out.append('%s\t%s\t%s'%(allmut[genename][0],allmut[genename][1],lines))
        if genename in allHS:
            allHS_out.append('%s\t%s\t%s'%(allHS[genename][0],allHS[genename][1],lines))

# load genomes and prokka
all_genomes = []
allprokka = []
for lines in open(allreference_list, 'r'):
    if lines.startswith('##reference=file:'):
        # set database
        database_file = lines.split('##reference=file:')[1].split('\n')[0]
        all_genomes.append(database_file)
        allprokka.append(database_file + '.faa.prokka.txt')

# sum prokka from reports
for genome in all_genomes:
    if 'prokka_' not in genome:
        genomefolder, filename = os.path.split(genome)
        output_dir = '%s/prokka_%s' % (genomefolder, filename)
        prokka_result = glob.glob('%s/PROKKA_*.gff'%(output_dir))
        if prokka_result != []:
            try:
                f1 = open(genome + '.faa.prokka.txt','r')
            except IOError:
                donor_species = os.path.split(genome)[-1].split('.all.spades')[0]
                donor_species_new = donor_species.replace('_clustercluster', '_CL').replace('_PB_', '_PaDi_')
                gene_set = load_gene(genome + '.faa')
                process_prokka(prokka_result[0],gene_set,genome)
        else:
            print('no prokka for',genome)

# sum all prokka and denovo mutations
allmut = load_mut(denovo_mut)
allHS = load_mut(HS_mut)
allmut_out = []
allHS_out = []
alloutput = set()

# load prokka
for prokka in allprokka:
    load_prokka(prokka)

for genename in allmut:
    if genename not in alloutput:
        allmut_out.append('%s\t%s\t\t\t\t\n' % (allmut[genename][0], allmut[genename][1]))

f1 = open(denovo_mut + '.prokka.txt', 'w')
f1.write('genename\tgenelength\tshort_genename\tgene\tCOG\tproduct\n')
f1.write(''.join(allmut_out))
f1.close()

f1 = open(HS_mut + '.prokka.txt', 'w')
f1.write('genename\tgenelength\tshort_genename\tgene\tCOG\tproduct\n')
f1.write(''.join(allHS_out))
f1.close()