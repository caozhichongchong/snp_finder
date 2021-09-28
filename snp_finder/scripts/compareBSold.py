import os,glob
from Bio import SeqIO
import statistics
import numpy as np
from Bio.Seq import Seq

input_bs_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS/binding_results_ccpA.txt'
ref_BS = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS/ccpA_BS_RegPrecise_difflength.fa'
vcf_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
output_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS/'

No_BS_pick = 10# top 10 BS
mut_cutoff = 0.1 # 10% -> 5bp
mut_cutoff2 = 5

def find_strains(vcf_file,genomewithSNP):
    mut_strains = []
    for linesvcf in open(vcf_file, 'r'):
        if linesvcf.startswith('CHR'):
            linesvcf_set = linesvcf.split('\n')[0].split('\t')
            allgenome = linesvcf_set[9:]
            i = 1
            # find mutated strains
            for genome in allgenome:
                if str(i) in genomewithSNP:
                    mut_strains.append(genome)
                i += 1
            break
    return [mut_strains,allgenome]
 # compare BS SNPs
def compare_BS(seq, seq2, mut_cutoff_set=0):
    alldiff = 0
    for i in range(0, len(seq)):
        if seq2[i] != seq[i]:
            alldiff += 1
        if mut_cutoff_set != 0 and alldiff > mut_cutoff_set:
            break
    return alldiff

def load_genes(input_faa):
    Mapping_loci_all = dict()
    for record in SeqIO.parse(input_faa, 'fasta'):
        record_id = str(record.id)
        contig = '_'.join(record_id.split('_')[0:-1])
        description = str(record.description).replace(' ', '').split('#')
        Mapping_loci_all.setdefault(contig, [])
        Mapping_loci_all[contig].append([int(description[1]) - 1,
                                         int(description[2]) - 1, record_id])
    return Mapping_loci_all

def load_BS(BS_file,Mapping_loci_all):
    allBS = []
    allBS.append('BS\tpvalue\tlocus\tcontig\tstrand\ttargetgane\tlocusgene\n')
    target_gene_list = dict()
    for lines in open(BS_file, 'r'):
        i = 0
        if not lines.startswith('#') and not lines.startswith('motif_id') and lines != '\n':
            lines_set = lines.split('\n')[0].split('\t')
            if i < No_BS_pick:
                i+=1
                pvalue = lines_set[7]
                contig, locus1, locus2, strand = lines_set[2:6]
                locus1 = int(locus1)
                locus2 = int(locus2)
                targetgene = ''
                locus_target = 0
                if contig in Mapping_loci_all:
                    for locus in Mapping_loci_all[contig]:
                        locusre1, locusref2, genename = locus
                        if locus2 <= locusref2 and targetgene == '':
                            targetgene = genename
                            locus_target = locusre1
                    seq = lines_set[9]
                    allBSset.setdefault(seq, [set(), set()])
                    if genomename in mut_strains:
                        allBSset[seq][-1].add(genomename)
                    else:
                        allBSset[seq][0].add(genomename)
                    if targetgene != '':
                        if strand == '-':
                            # the gene before
                            gene_locus = int(targetgene.split('_')[-1])
                            if gene_locus > 1:
                                targetgene = '_'.join(targetgene.split('_')[0:-1]) + '_%s' % (
                                        int(targetgene.split('_')[-1]) - 1)
                    else:
                        targetgene='%s_1'%(contig)
                    allBS.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                            seq, pvalue, locus1, contig, strand, targetgene, locus_target))
                    target_gene_list.setdefault(targetgene, set())
                    target_gene_list[targetgene].add(seq)
    f1 = open('%s/%s/%s.BS.txt' % (output_file, genomename, genomename), 'w')
    f1.write(''.join(list(set(allBS))))
    f1.close()
    aa_output = []
    genomename_short = genomename.replace('_BL_', '_S')
    for record in SeqIO.parse(input_faa, 'fasta'):
        record_id = str(record.id)
        if record_id in target_gene_list:
            for seq in target_gene_list[record_id]:
                aa_output.append('>%s_%s_C_%s_G_%s\n%s\n' % (
                    seq, genomename_short, record_id.split('_')[1], record_id.split('_')[-1], str(record.seq)))
                select_seq_faa.setdefault(seq,'>%s_%s_C_%s_G_%s\n%s\n' % (
                        seq, genomename_short, record_id.split('_')[1], record_id.split('_')[-1], str(record.seq)))
    f1 = open('%s/%s/%s.BS.faa' % (output_file, genomename, genomename), 'w')
    f1.write(''.join(aa_output))
    f1.close()

def compareBS():
    BS_diff = dict()
    alldiff_set = []
    for seq in allBSset:
        BS_diff.setdefault(seq, set())
        if mut_cutoff2 == 0:
            # set cutoff as 10% top similar -> 5bp
            for seq2 in allBSset:
                if seq2 != seq:
                    alldiff = compare_BS(seq, seq2)
                    alldiff_set.append(alldiff)
            newmut_cutoff = np.quantile(alldiff_set, [0.1])[0]
        else:
            # preset cutoff
            newmut_cutoff = mut_cutoff2
        for seq2 in allBSset:
            if seq2 != seq:
                alldiff = compare_BS(seq, seq2, newmut_cutoff)
                if alldiff <= newmut_cutoff:
                    BS_diff[seq].add(seq2)
    return [BS_diff,alldiff_set]
# whether BS in some mut, not in all wt
def select_BS(list_seq):
    selected = False
    no_mut = len(list_seq[-1])
    no_wt = len(list_seq[0])
    if no_mut > 0 and no_wt < (len(allgenome)-len(mut_strains))*0.5:
        selected = True
    return [no_mut, no_wt, selected]

def select_reversecomplement(seq):
    seq_rc = str(Seq(seq).reverse_complement())
    seq_set = [seq,seq_rc]
    seq_set.sort()
    return seq == seq_set[0]

def find_candidate_mut_BS():
    allBS_all = dict()
    allseq = list(allBSset.keys())
    allBS_select = dict()
    for seq in allBSset:
        inref = False
        if seq in Ref:
            inref = True
        no_mut, no_wt, selected = select_BS(allBSset[seq])
        withsim_wt = ''
        if selected:
            if BS_diff[seq] != set():
                for seq2 in BS_diff[seq]:
                    if len(allBSset[seq2][0]) > 0 and not any(mut in allBSset[seq][-1] for mut in allBSset[seq2][-1]):
                        # does not share mut strains, some wt has it
                        # potential mutated BS from wt BS
                        # wt BS similar to mutated BS
                        withsim_wt += '%s;' % (allseq.index(seq2))
                        allBS_select.setdefault(seq, set())
                        allBS_select[seq].add(seq2)
            if withsim_wt == '':
                # no similar wt
                allBS_select[seq] = set()
        allBS_all.setdefault(seq, ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
            allseq.index(seq),
            seq, no_wt, no_mut, withsim_wt, selected,
            inref, ';'.join(BS_diff[seq]),
            ';'.join(allBSset[seq][0]), ';'.join(allBSset[seq][-1])
        )))
    # output result
    allBS = []
    allBS.append('SNPdiff\tBS_order\tBS\tNo.wt\tNo.mut\tmut_hit\twithsim_wt\tref\tsimilarseq\twt\tmut\n')
    allseqout = []
    for seq in allBS_select:
        if select_reversecomplement(seq):
            # one orientation
            if allBS_select[seq] != set():
                allBS.append('%s\t%s' % (0, allBS_all[seq]))
                for seq2 in allBS_select[seq]:
                    alldiff = compare_BS(seq, seq2)
                    allBS.append('%s\t%s' % (alldiff, allBS_all[seq2]))
                    allseqout.append(select_seq_faa.get(seq2, ''))
                allBS.append('\n')
                allseqout.append(select_seq_faa.get(seq,''))
    f1 = open('%s.BS.txt' % (output_file_BS), 'w')
    f1.write(''.join(allBS))
    f1.close()
    if allseqout!=[] and not all(gene =='' for gene in allseqout):
        fasta_output = '%s.BS.faa' % (output_file_BS)
        f1 = open(fasta_output, 'w')
        f1.write(''.join(allseqout))
        f1.close()
        # run eggnog
        annotate(fasta_output)

def annotate(fasta_output):
    cutoff = 0.7
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                           % ('usearch', fasta_output, cutoff, fasta_output,
                              fasta_output, 40))
    os.system(cmd_cluster)
    fasta_output = fasta_output + '.cluster.aa'
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/eggnog/xaa.hmm'
    cmds = ('hmmsearch --tblout %s.eggnog.1.txt --cpu 40 -E %s %s %s\n') % (
    fasta_output, cutoff, database, fasta_output)
    database = '/scratch/users/mit_alm/database/eggnog/xab.hmm'
    cmds += ('hmmsearch --tblout %s.eggnog.2.txt --cpu 40 -E %s %s %s\n') % (
        fasta_output, cutoff, database, fasta_output)
    database = '/scratch/users/mit_alm/database/eggnog/xac.hmm'
    cmds += ('hmmsearch --tblout %s.eggnog.3.txt --cpu 40 -E %s %s %s\n') % (
        fasta_output, cutoff, database, fasta_output)
    f1 = open(output_file_BS + '.eggnog.sh', 'w')
    f1.write(
        '#!/bin/bash\nsource ~/.bashrc\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (
            cmds))
    f1.close()

# load ref
Ref = []
if ref_BS != 'None':
    for record in SeqIO.parse(ref_BS, 'fasta'):
        Ref.append(str(record.seq))

# process each SNP
for lines in open(input_bs_file,'r'):
    if not lines.startswith('AA_POS_ref'):
        lines_set = lines.split('\t')
        lineage = lines_set[4].split('__')[0]
        species = lines_set[4].split('_')[0]
        donor = lines_set[5]
        SNP = lines_set[3]
        if SNP not in ['K226*','A23V','G12R','A112V']:
            # find genome names
            vcf_file = '%s/%s%s'%(vcf_folder,lineage.replace('CL','clustercluster'),'.all.parsi.fasta.linktrunc.sum.txt')
            print(vcf_file)
            mut_strains, allgenome = find_strains(vcf_file,lines_set[-9].split(';'))
            print(mut_strains)
            # process fino results
            output_file = '%s/%s_%s'%(output_folder,species,donor)
            output_file_BS = '%s/%s_%s_%s'%(output_folder,species,donor,SNP)
            print(output_file_BS)
            allBSset = dict()
            # load BS
            select_seq_faa = dict()
            for BS_folder in glob.glob('%s/*' % (output_file)):
                genomename = os.path.split(BS_folder)[-1]
                if genomename in allgenome:
                    # load BS file and target genes
                    BS_file = glob.glob('%s/fimo.tsv' % (BS_folder))[0]
                    input_faa = '%s/%s/%s.faa' % (output_file, genomename,genomename)
                    # load all gene position
                    Mapping_loci_all = load_genes(input_faa)
                    # load BS
                    load_BS(BS_file,Mapping_loci_all)
            # compare BS differences
            BS_diff,alldiff_set = compareBS()
            # find candidate mut BS
            find_candidate_mut_BS()

f1 = open(os.path.join(output_folder, 'allanno.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(output_folder, '*eggnog.sh')):
    f1.write('jobmit %s %s small1\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/allanno.sh'%(output_folder))
