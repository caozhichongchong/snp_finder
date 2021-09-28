import os,glob
from Bio import SeqIO
import statistics
import numpy as np
from Bio.Seq import Seq
import re

vcf_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
output_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS/'
target_TF = '%s/target.TF.faa'%(output_folder)
ref_BS = '%s/ccpA_BS.fa'%(output_folder)
input_bs_file = '%s/binding_results_ccpA.txt'%(output_folder)
BS_folder = ''
distance = 5000 # extract 5kb neighbourhood
distance_TF = 5000 # TF with BS nearby

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

def load_BS_old(BS_file,Mapping_loci_all):
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

def load_BS(BS_file,Mapping_loci_all,BS_loci,TF_locus):
    allBS = []
    allBS.append('BS\tpvalue\tlocus\tcontig\tstrand\ttargetgane\tlocusgene\n')
    target_gene_list = dict()
    i = 0
    genomename = os.path.split(BS_file)[-1].split('.fimo.tsv')[0]
    check_TF_BS = False
    TF_locus_genome = TF_locus.get(genomename,[])
    TF_contig = [i[0] for i in TF_locus_genome]
    if TF_locus_genome != []:
        # check TF
        for lines in open(BS_file, 'r'):
            if not check_TF_BS and not lines.startswith('#') and not lines.startswith('motif_id') and lines != '\n':
                lines_set = lines.split('\n')[0].split('\t')
                contig, locus1, locus2, strand = lines_set[2:6]
                if contig in TF_contig:
                    locus1 = int(locus1)
                    locus2 = int(locus2)
                    for contig_TF,locus1_TF,locus2_TF,gene in TF_locus_genome:
                        if contig_TF == contig and locus1_TF - distance_TF <= locus2 and locus2_TF + distance_TF >= locus1:
                            #print(locus1_TF - distance_TF,locus2_TF + distance_TF,locus1,locus2,gene)
                            check_TF_BS = True
                            break
        if not check_TF_BS:
            print('not BS near TF for genome %s'%(genomename))
        else:
            # BS near TF
            print('found BS near TF for genome %s' % (genomename))
            pass_genome.append(genomename)
            # load BS
            for lines in open(BS_file, 'r'):
                if not lines.startswith('#') and not lines.startswith('motif_id') and lines != '\n':
                    lines_set = lines.split('\n')[0].split('\t')
                    pvalue = lines_set[7]
                    contig, locus1, locus2, strand = lines_set[2:6]
                    i += 1
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
                        BS_loci.setdefault(seq, [])
                        BS_loci[seq].append([genomename, locus1, locus2, '%s_C_%s_G_%s' % (genomename,
                                                                                           targetgene.split('_')[1],
                                                                                           targetgene.split('_')[-1])])
                # output BS and gene target
            # output BS and gene target
            f1 = open('%s.BS.txt' % (output_file), 'w')
            f1.write(''.join(allBS))
            f1.close()
            aa_output = []
            for record in SeqIO.parse(input_faa, 'fasta'):
                record_id = str(record.id)
                if record_id in target_gene_list:
                    for seq in target_gene_list[record_id]:
                        record_seq = str(record.seq)
                        aa_output.append('>%s_%s_C_%s_G_%s\n%s\n' % (
                            seq, genomename, record_id.split('_')[1], record_id.split('_')[-1], record_seq))
                        select_seq_faa.setdefault(seq,['',set()])
                        if record_seq not in select_seq_faa[seq][-1]:
                            select_seq_faa[seq][-1].add(record_seq)
                            select_seq_faa[seq][0]+=('>%s_%s_C_%s_G_%s\n%s\n' % (
                                    seq, genomename, record_id.split('_')[1], record_id.split('_')[-1], record_seq))
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

def find_mutations_on_BS_old(vcf_file2,mut_strains,BS_loci):
    BS_SNP_output = []
    mut_strains_set = []
    wild_type_set = []
    allseqout = set()
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
    #print(mut_strains,mut_strains_set,wild_type_set)
    for linesvcf in open(vcf_file2, 'r'):
        linesvcf_set = linesvcf.split('\n')[0].split('\t')
        if not linesvcf.startswith('#'):
            contig, POS = linesvcf_set[0:2]
            if contig in BS_loci:
                POS = int(POS)
                for BSpos in BS_loci[contig]:
                    BSpos1, BSpos2, seq, targetgene,pvalue = BSpos
                    if POS >= BSpos1 - distance and POS <= BSpos2 + distance :
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
                            allseqout.add(select_seq_faa.get(seq, ''))
    f1 = open('%s.BSsum.txt' % (output_file_BS),'w')
    f1.write(''.join(BS_SNP_output))
    f1.close()
    if len(allseqout)>0:
        f1 = open('%s.BS.faa' % (output_file_BS), 'w')
        f1.write(''.join(list(allseqout)))
        f1.close()
        annotate('%s.BS.faa' % (output_file_BS))

def find_mutations_on_BS(vcf_file2,mut_strains,BS_loci_co_assembly):
    BS_SNP_output = []
    mut_strains_set = []
    wild_type_set = []
    allseqout = set()
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
    BS_SNP_output.append('BS\tBS_in_ref\tCHR\tPOS\tPOS_on_BS\tMutated_alleles\tWildtype_alleles\t%s\t%s\n' % (
        '\t'.join(mut_strains_set),
        '\t'.join(
            wild_type_set)
    ))
    #print(mut_strains,mut_strains_set,wild_type_set)
    finished_set = []
    for linesvcf in open(vcf_file2, 'r'):
        linesvcf_set = linesvcf.split('\n')[0].split('\t')
        if not linesvcf.startswith('#'):
            contig, POS = linesvcf_set[0:2]
            if contig in BS_loci_co_assembly:
                POS = int(POS)
                for seq_loci_sub,seq in BS_loci_co_assembly[contig]:
                    if [seq_loci_sub,seq] not in finished_set:
                        finished_set.append([seq_loci_sub,seq])
                        BSpos1, BSpos2 = seq_loci_sub
                        if POS >= BSpos1 - distance and POS <= BSpos2 + distance:
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
                                inref = False
                                if seq in Ref:
                                    inref = True
                                #targetgene = [i[-1] for i in BS_loci[seq]]
                                BS_SNP_output.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                                    seq, inref, contig, POS, POS - BSpos1 + 1,
                                    ';'.join(list(Mut_allele_set)),
                                    ';'.join(list(Wild_allele_set)),
                                    '\t'.join(Mut_allele),
                                    '\t'.join(Wild_allele)
                                ))
                                allseqout.add(select_seq_faa.get(seq, '')[0])
    f1 = open('%s.BSsum.txt' % (output_file_BS), 'w')
    f1.write(''.join(BS_SNP_output))
    f1.close()
    if len(allseqout)>0:
        f1 = open('%s.BS.faa' % (output_file_BS), 'w')
        f1.write(''.join(list(allseqout)))
        f1.close()
        annotate('%s.BS.faa' % (output_file_BS))

def find_TF(input_faa,blast_out,alloutput):
    withTF = False
    genomename = os.path.split(input_faa)[-1].split('.faa')[0]
    #donor = genomename.split('_')[0]
    TF_set = dict()
    for lines in open(blast_out):
        lines_set = lines.split('\t')
        query,reference = lines_set[0:2]
        TF_set.setdefault(query,[])
        TF_set[query].append(reference)
    for record in SeqIO.parse(input_faa, 'fasta'):
        record_id = str(record.id)
        if record_id in TF_set:
            for reference in TF_set[record_id]:
                #if donor in reference:
                    description = str(record.description).replace(' ', '').split('#')
                    alloutput.append('%s\t%s\t%s\t%s\t%s\n'%(genomename,record_id,description[1],description[2],reference))
                    withTF = True
    if not withTF:
        print('TF not found for %s'%(input_faa))
    return alloutput

def load_TF(TF_locus_file):
    TF_locus = dict()
    for lines in open(TF_locus_file,'r'):
        lines_set = lines.split('\t')
        genome, gene,locus1,locus2 = lines_set[0:4]
        genome = genome.split('.all')[0]
        contig = '_'.join(gene.split('_')[0:-1])
        TF_locus.setdefault(genome,[])
        TF_locus[genome].append([contig,int(locus1),int(locus2),gene])
    return TF_locus

def find_target_BS(mut_strains, BS_loci,allgenome):
    BS_select = set()
    wild_strains = [i for i in allgenome if i not in mut_strains]
    if any(genome in mut_strains for genome in pass_genome) and any(genome in wild_strains for genome in pass_genome):
        # mutated genome and wild strains all passed
        print('found enough genomes passing the criterion of BS near TF')
        for seq in BS_loci:
            allgenomewithBS = [i[0] for i in BS_loci[seq]]
            # only in mutated
            if (any(genome in mut_strains for genome in allgenomewithBS) and not any(
                    genome in wild_strains for genome in allgenomewithBS)):
                BS_select.add(seq)
    else:
        print('not enough genomes passing the criterion of BS near TF')
        BS_select = set()
    return BS_select

def find_BS_coassembly(co_assembly,BS_select):
    BS_loci_co_assembly = dict()
    for record in SeqIO.parse(co_assembly, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        for seq in BS_select:
            seq_loci = [[m.start()+1, m.end()+1] for m in re.finditer(seq, record_seq)]
            if seq_loci!=[]:
                BS_loci_co_assembly.setdefault(record_id,[])
                for seq_loci_sub in seq_loci:
                    BS_loci_co_assembly[record_id].append([seq_loci_sub,seq])
            seq_rev = str(Seq(seq).reverse_complement())
            seq_loci = [[m.start()+1, m.end()+1] for m in re.finditer(seq_rev, record_seq)]
            if seq_loci!=[]:
                if seq_rev in select_seq_faa:
                    seq = seq_rev
                BS_loci_co_assembly.setdefault(record_id, [])
                for seq_loci_sub in seq_loci:
                    BS_loci_co_assembly[record_id].append([seq_loci_sub, seq])
    return BS_loci_co_assembly

def TF_near_BS(TF_locus,BS_file_all):
    BS_TF = []
    for BS_file in BS_file_all:
        check_TF_BS = False
        genome = os.path.split(BS_file)[-1].split('.BS.txt')[0]
        TF_locus_genome = TF_locus.get(genome, [])
        TF_contig = [i[0] for i in TF_locus_genome]
        for lines in open(BS_file,'r'):
            if not lines.startswith('BS'):
                lines_set = lines.split('\n')[0].split('\t')
                seq,contig,locus = lines_set[0:3]
                if contig in TF_contig:
                    locus = int(locus)
                    for contig_TF, locus1_TF, locus2_TF,gene in TF_locus_genome:
                        if contig_TF == contig and locus1_TF - distance_TF <= locus and locus2_TF + distance_TF >= locus:
                            check_TF_BS = True
                            BS_TF.append('%s\t%s\t%s\t%s\t%s\n'%(genome,contig,seq,locus,gene))
        if not check_TF_BS:
            print('BS not found near TF %s'%(genome))
    f1 = open('%s/allBSnearTF.txt' % (output_folder), 'w')
    f1.write(''.join(BS_TF))
    f1.close()

# find all TF
try:
    f1 = open('%s/allTFloci.txt' % (output_folder), 'r')
except FileNotFoundError:
    cmds = 'export LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n'
    os.system(cmds)
    alloutput = []
    allinput_faa = glob.glob('%s/*/*.faa'%(output_folder))
    for input_faa in allinput_faa:
        cmds = (
            "diamond blastp --query %s --db %s.dmnd --out %s.TF.txt --id 80 --outfmt 6 --max-target-seqs 10 --evalue 1e-1 --threads 40\n"%
            (
                input_faa,target_TF,input_faa,
            ))
        os.system(cmds)
        alloutput = find_TF(input_faa, input_faa + '.TF.txt', alloutput)
    f1 = open('%s/allTFloci.txt' % (output_folder), 'w')
    f1.write(''.join(alloutput))
    f1.close()

# load TF loci
TF_locus = load_TF('%s/allTFloci.txt' % (output_folder))

# find BS predicted in genomes near TF
try:
    f1 = open('%s/allBSnearTF.txt' % (output_folder), 'r')
except FileNotFoundError:
    BS_file_all = glob.glob('%s/*/*.BS.txt'%(output_folder))
    #BS_file_all = glob.glob('%s/*.BS.txt' % (output_folder))
    TF_near_BS(TF_locus, BS_file_all)
    print('finished processing TF near BS')

# load ref
Ref = []
if ref_BS != 'None':
    for record in SeqIO.parse(ref_BS, 'fasta'):
        Ref.append(str(record.seq))

# process each SNP
for lines in open(input_bs_file,'r'):
    if not lines.startswith('AA_POS_ref'):
        try:
            lines_set = lines.split('\t')
            lineage = lines_set[4].split('__')[0].replace('CL', 'clustercluster')
            species = lines_set[4].split('_')[0]
            donor = lines_set[5]
            SNP = lines_set[3]
            output_file_BS = '%s/%s/%s_%s_%s' % (output_folder, BS_folder,species, donor, SNP)
            print(output_file_BS)
            try:
                f1 = open('%s.BSsum.txt' % (output_file_BS), 'r')
            except FileNotFoundError:
                select_seq_faa = dict()
                # find genome names
                vcf_file = '%s/%s%s' % (
                    vcf_folder, lineage, '.all.parsi.fasta.linktrunc.sum.txt')
                genomewithSNP = lines_set[-9].split(';')
                mut_strains, allgenome = find_strains(vcf_file, genomewithSNP)
                print('process %s mutated strains %s'%(lineage,mut_strains))
                lineage = lineage.split('.donor')[0]
                # process fino results
                donor_species = '%s_%s'%(species,donor)
                BS_file_all = glob.glob('%s/%s/%s/*.fimo.tsv'%(output_folder,BS_folder,donor_species))
                BS_loci = dict()
                pass_genome = []
                for BS_file in BS_file_all:
                    output_file = BS_file.split('.fimo.tsv')[0]
                    input_faa = '%s/%s/%s.faa'%(output_folder,donor_species,os.path.split(output_file)[-1])
                    # load all gene position
                    Mapping_loci_all = load_genes(input_faa)
                    # load BS
                    BS_loci = load_BS(BS_file, Mapping_loci_all,BS_loci,TF_locus)
                # find targeted BS only in mutated or only in WT
                BS_select = find_target_BS(mut_strains, BS_loci,allgenome)
                if BS_select!= set():
                    # find genome locus for BS_select
                    co_assembly = '%s/co-assembly/%s.all.spades2.fasta'%(output_folder,lineage)
                    BS_loci_co_assembly = find_BS_coassembly(co_assembly,BS_select)
                    print(BS_loci_co_assembly)
                    # find mutations on BS in mutated strains
                    vcf_file2 = '%s/moredetails/%s.all.raw.vcf.all' % (
                        vcf_folder, lineage)
                    print(vcf_file2)
                    find_mutations_on_BS(vcf_file2, mut_strains, BS_loci_co_assembly)
            # filter BS
            BS_file = '%s.BSsum.txt' % (output_file_BS)
            # filter BS
            BS_set = set()
            for lines in open(BS_file, 'r'):
                if not lines.startswith('BS'):
                    lines_set = lines.split('\t')
                    BS = lines_set[0]
                    POS_on_BS = int(lines_set[4])
                    Mut = lines_set[5]
                    if POS_on_BS >= -20 and POS_on_BS <= 20:
                        # BS with changes on BS itself and not just loss
                        BS_set.add(BS)
            BS_file_out = []
            for lines in open(BS_file, 'r'):
                if lines.startswith('BS'):
                    BS_file_out.append(lines)
                elif lines != '':
                    seq = lines.split('\t')[0]
                    if seq in BS_set:
                        # BS with changes on BS itself and not just loss
                        BS_file_out.append(lines)
            f1 = open('%s.BSsum.filtered.txt' % (output_file_BS), 'w')
            f1.write(''.join(BS_file_out))
            f1.close()
        except FileNotFoundError:
            pass


alleggnog = glob.glob(os.path.join(output_folder, '%s/*eggnog.sh'%(BS_folder)))
if alleggnog!= []:
    f1 = open(os.path.join(output_folder, '%s/allanno.sh'%(BS_folder)), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in alleggnog:
        f1.write('jobmit %s %s small1\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()
    print('please run %s/%s/allanno.sh'%(output_folder,BS_folder))
