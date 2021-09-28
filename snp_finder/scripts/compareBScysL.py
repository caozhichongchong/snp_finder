import os,glob
from Bio import SeqIO
import statistics
import numpy as np
from Bio.Seq import Seq

input_bs_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS_cysL/binding_results_cysL.txt'
ref_BS = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS_cysL/BS.fasta'
vcf_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
output_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS_cysL/'

No_BS_pick = 50# notused

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

def contig_length(CHR):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        try:
            total_length = CHR.split('length_')[1].split('_cov')[0]
        except IndexError:
            total_length = 10000
    return int(total_length)

def load_BS(BS_file,Mapping_loci_all):
    allBS = []
    allBS.append('BS\tpvalue\tlocus\tcontig\tstrand\ttargetgane\tlocusgene\n')
    target_gene_list = dict()
    BS_loci = dict()
    i = 0
    for lines in open(BS_file, 'r'):
        if not lines.startswith('#') and not lines.startswith('motif_id') and lines != '\n':
            lines_set = lines.split('\n')[0].split('\t')
            pvalue = lines_set[7]
            contig, locus1, locus2, strand = lines_set[2:6]
            if contig_length(contig) >= 5000:
                i += 1
                BS_loci.setdefault(contig, [])
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
                    if targetgene != '':
                        if strand == '-':
                            # the gene before
                            gene_locus = int(targetgene.split('_')[-1])
                            if gene_locus > 1:
                                targetgene = '_'.join(targetgene.split('_')[0:-1]) + '_%s' % (
                                        int(targetgene.split('_')[-1]) - 1)
                    else:
                        targetgene = '%s_1' % (contig)
                    allBS.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        seq, pvalue, locus1, contig, strand, targetgene, locus_target))
                    target_gene_list.setdefault(targetgene, set())
                    target_gene_list[targetgene].add(seq)
                    BS_loci[contig].append([locus1, locus2, seq, targetgene,pvalue])
    f1 = open('%s.BS.txt' % (output_file), 'w')
    f1.write(''.join(allBS))
    f1.close()
    aa_output = []
    for record in SeqIO.parse(input_faa, 'fasta'):
        record_id = str(record.id)
        if record_id in target_gene_list:
            for seq in target_gene_list[record_id]:
                aa_output.append('>%s_%s_C_%s_G_%s\n%s\n' % (
                    seq, lineage, record_id.split('_')[1], record_id.split('_')[-1], str(record.seq)))
                select_seq_faa.setdefault(seq,'>%s_%s_C_%s_G_%s\n%s\n' % (
                        seq, lineage, record_id.split('_')[1], record_id.split('_')[-1], str(record.seq)))
    f1 = open('%s.BS.faa' % (output_file), 'w')
    f1.write(''.join(aa_output))
    f1.close()
    return BS_loci

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

def allele_freq_to_allele(genotype, ALT_set):
    genotype = [int(i) for i in genotype.split(':')[-1].split(',')]
    Major_ALT = '-'
    if sum(genotype) > 0:
        ALT_set_sample = dict()
        ALT_frq_set = set()
        for i in range(0,len(ALT_set)):
            ALT_frq = int(genotype[i])
            ALT_set_sample.setdefault(ALT_frq, set())
            ALT_set_sample[ALT_frq].add(ALT_set[i])
            ALT_frq_set.add(ALT_frq)
        ALT_frq_set = sorted(ALT_frq_set, reverse=True)
        for ALT_frq in ALT_frq_set:
            for alleles in ALT_set_sample[ALT_frq]:
                if Major_ALT == '-':
                    Major_ALT = alleles
    return Major_ALT

def find_mutations_on_BS(vcf_file2,mut_strains,BS_loci):
    BS_SNP_output = []
    mut_strains_set = []
    wild_type_set = []
    allseqout = []
    genomewithSNPvcf = []
    genomewithnoSNPvcf = []
    samplefile = vcf_file2.replace('.all.raw.vcf.all','.all.raw.vcf.filtered.samplename.txt')
    for lines in open(samplefile, 'r'):
        # set sample name
        sample_set = lines.split('\n')[0].split('\t')
        i = 9
        for samples in sample_set:
            genome = samples
            if genome in mut_strains:
                genomewithSNPvcf.append(i)
                mut_strains_set.append('M.%s' % (genome))
            elif genome in allgenome:
                genomewithnoSNPvcf.append(i)
                wild_type_set.append('W.%s' % (genome))
            i += 1
    BS_SNP_output.append('BS\tBS_in_ref\tCHR\tPOS\tPOS_on_BS\tTarget\tpvalue\tMutated_alleles\tWildtype_alleles\t%s\t%s\n' % (
        '\t'.join(mut_strains_set),
        '\t'.join(
            wild_type_set)
    ))
    print(mut_strains,mut_strains_set,wild_type_set)
    for linesvcf in open(vcf_file2, 'r'):
        linesvcf_set = linesvcf.split('\n')[0].split('\t')
        if not linesvcf.startswith('#'):
            contig, POS = linesvcf_set[0:2]
            if contig in BS_loci:
                POS = int(POS)
                for BSpos in BS_loci[contig]:
                    BSpos1, BSpos2, seq, targetgene,pvalue = BSpos
                    if POS >= BSpos1 and POS <= BSpos2:
                        REF = linesvcf_set[3]
                        ALT_set = [REF]
                        ALT = linesvcf_set[4]
                        for a_ALT in ALT.split(','):
                            if a_ALT != '.':
                                ALT_set.append(a_ALT)
                        # BS site
                        Mut_allele = []
                        Wild_allele = []
                        # Major_ALT in mutated strains
                        for i in genomewithSNPvcf:
                            genotype = linesvcf_set[i]
                            Major_ALT = allele_freq_to_allele(genotype, ALT_set)
                            Mut_allele.append(Major_ALT)
                        # Major_ALT in wild type strains
                        for i in genomewithnoSNPvcf:
                            genotype = linesvcf_set[i]
                            Major_ALT = allele_freq_to_allele(genotype, ALT_set)
                            Wild_allele.append(Major_ALT)
                        Mut_allele_set = sorted(set(Mut_allele))
                        Wild_allele_set = sorted(set(Wild_allele))
                        if any(i not in Wild_allele_set for i in Mut_allele_set):
                            # different alleles
                            print(POS, BSpos1, BSpos2)
                            print('M', Mut_allele_set, 'W', Wild_allele_set)
                            inref = False
                            if seq in Ref:
                                inref = True
                            BS_SNP_output.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                                seq, inref, contig, POS, POS - BSpos1 + 1, targetgene,pvalue,
                                ';'.join(list(Mut_allele_set)),
                                ';'.join(list(Wild_allele_set)),
                                '\t'.join(Mut_allele),
                                '\t'.join(Wild_allele)
                            ))
                            allseqout.append(select_seq_faa.get(seq, ''))
    f1 = open('%s.BSsum.txt' % (output_file_BS),'w')
    f1.write(''.join(BS_SNP_output))
    f1.close()
    if len(allseqout)>0:
        f1 = open('%s.BS.faa' % (output_file_BS), 'w')
        f1.write(''.join(allseqout))
        f1.close()
        annotate('%s.BS.faa' % (output_file_BS))

# load ref
Ref = []
if ref_BS != 'None':
    for record in SeqIO.parse(ref_BS, 'fasta'):
        Ref.append(str(record.seq))

# process each SNP
for lines in open(input_bs_file,'r'):
    if not lines.startswith('AA_POS_ref'):
        lines_set = lines.split('\t')
        lineage = lines_set[4].split('__')[0].replace('CL', 'clustercluster')
        species = lines_set[4].split('_')[0]
        donor = lines_set[5]
        SNP = lines_set[3]
        select_seq_faa = dict()
        # find genome names
        vcf_file = '%s/%s%s' % (
            vcf_folder, lineage, '.all.parsi.fasta.linktrunc.sum.txt')
        genomewithSNP = lines_set[-9].split(';')
        mut_strains, allgenome = find_strains(vcf_file, genomewithSNP)
        # process fino results
        lineage = lineage.split('.donor')[0]
        output_file = '%s/%s' % (output_folder, lineage)
        BS_file = '%s/%s.fimo.tsv' % (output_folder, lineage)
        input_faa = '%s/co-assembly/%s.all.spades2.fasta.faa' % (output_folder, lineage)
        # load all gene position
        Mapping_loci_all = load_genes(input_faa)
        # load BS
        BS_loci = load_BS(BS_file, Mapping_loci_all)
        # find mutations on BS in mutated strains
        vcf_file2 = '%s/moredetails/%s.all.raw.vcf.all' % (
            vcf_folder, lineage)
        print(vcf_file2)
        output_file_BS = '%s/%s_%s_%s' % (output_folder, species, donor, SNP)
        find_mutations_on_BS(vcf_file2, mut_strains, BS_loci)

f1 = open(os.path.join(output_folder, 'allanno.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(output_folder, '*eggnog.sh')):
    f1.write('jobmit %s %s small1\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/allanno.sh'%(output_folder))