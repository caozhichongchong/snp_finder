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
                    for contig_TF,locus1_TF,locus2_TF,gene,TF in TF_locus_genome:
                        if contig_TF == contig and locus1_TF - distance_TF <= locus2 and locus2_TF + distance_TF >= locus1:
                            print(int(locus1 - locus1_TF), int(locus2_TF - locus2),contig_TF,lines_set[-1])
                            check_TF_BS = True
                            break
        if not check_TF_BS:
            print('no BS near TF for genome %s'%(genomename))
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
                        BS_loci[seq].append([genomename, locus1, locus2, contig,'%s__C_%s_G_%s' % (genomename,
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

def find_mutations_on_BS(allvcf_file):
    BS_select_info_out = dict()
    # output BS select info
    for lines in open(output_file_BS + '.target.BS.info.txt', 'r'):
        if not lines.startswith('seq'):
            lines_set = lines.split('\n')[0].split('\t')
            seq, locus1, locus2, genomename, contig,targetgene,genelocus1,genelocus2 = lines_set
            BS_select_info_out.setdefault(genomename,dict())
            BS_select_info_out[genomename].setdefault(contig,[])
            BS_select_info_out[genomename][contig].append([int(locus1), int(locus2),seq,
                                                           min(int(genelocus1),int(locus1)),
                                                           max(int(genelocus2),int(locus2))])
    BS_SNP_output = []
    allseqout = set()
    BS_SNP_output.append('BS\tBS_in_ref\tCHR\tPOS\tPOS_on_BS\tMutated_alleles\tWildtype_alleles\tMut_strain\tWild_type_strains\n')
    for vcf_file in allvcf_file:
        mut_strain_name = os.path.split(vcf_file)[-1].split('.target.BS.raw.vcf')[0]
        if mut_strain_name in BS_select_info_out:
            allBSloci = BS_select_info_out[mut_strain_name]
            for linesvcf in open(vcf_file, 'r'):
                linesvcf_set = linesvcf.split('\n')[0].split('\t')
                if not linesvcf.startswith('#'):
                    contig, POS = linesvcf_set[0:2]
                    if contig in allBSloci:
                        POS = int(POS)
                        for BSpos1,BSpos2,seq,BSgenelocus1,BSgenelocus2 in allBSloci[contig]:
                            if POS >= BSgenelocus1 - distance and POS <= BSgenelocus2 + distance:
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
                                genotype = linesvcf_set[9]
                                Major_ALT = allele_freq_to_allele(genotype, ALT_set)
                                Mut_allele.append(Major_ALT)
                                # Major_ALT in wild type strains
                                for genotype in linesvcf_set[10:]:
                                    Major_ALT = allele_freq_to_allele(genotype, ALT_set)
                                    Wild_allele.append(Major_ALT)
                                Mut_allele_set = sorted(set(Mut_allele))
                                Wild_allele_set = sorted(set(Wild_allele))
                                if any(i not in Wild_allele_set for i in Mut_allele_set):
                                    # different alleles
                                    inref = False
                                    if seq in Ref:
                                        inref = True
                                    BS_SNP_output.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                                        seq, inref, contig, POS, POS - BSpos1,
                                        ';'.join(list(Mut_allele_set)),
                                        ';'.join(list(Wild_allele_set)),
                                        '\t'.join(Mut_allele),
                                        '\t'.join(Wild_allele)
                                    ))
    f1 = open('%s.BSsum.txt' % (output_file_BS), 'w')
    f1.write(''.join(BS_SNP_output))
    f1.close()
    BS_set = BS_sum_filter('%s.BSsum.txt' % (output_file_BS))
    for seq in BS_set:
        allseqout.add(select_seq_faa.get(seq, [''])[0])
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
        genome, gene,locus1,locus2,TF = lines_set
        donor1 = genome.split('_')[0]
        donor2 = TF.split('donor.')[1].split('__')[0]
        if donor1 == donor2:
            genome = genome.split('.all')[0]
            contig = '_'.join(gene.split('_')[0:-1])
            TF_locus.setdefault(genome,[])
            TF_locus[genome].append([contig,int(locus1),int(locus2),gene,TF])
    return TF_locus

def find_target_BS(mut_strains, BS_loci,allgenome):
    BS_select = set()
    BS_select_wt = set()
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
            if (any(
                genome in wild_strains for genome in allgenomewithBS) and not any(
                    genome in mut_strains for genome in allgenomewithBS)):
                BS_select_wt.add(seq)
    else:
        print('not enough genomes passing the criterion of BS near TF')
        BS_select = set()
        BS_select_wt = set()
    return [BS_select,BS_select_wt]

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

def genome_mapping(BS_loci,BS_select):
    # find BS_neighbour
    try:
        f1 = open(output_file_BS + '.target.BS.info.txt', 'r')
    except FileNotFoundError:
        BS_select_info = dict()
        BS_select_info_out = dict()
        for seq in BS_select:
            for allinfo in BS_loci[seq]:
                genomename, locus1, locus2, contig, targetgene = allinfo
                BS_select_info.setdefault(genomename, [])
                BS_select_info[genomename].append(contig)
                BS_select_info_out.setdefault(targetgene,'%s\t%s\t%s\t%s\t%s\t%s'%(seq,locus1,locus2,genomename,contig,targetgene))
        # find gene locus
        for genomename in BS_select_info:
            faa = '%s/%s/%s/%s.faa'%(output_folder,BS_folder,donor_species,genomename)
            for record in SeqIO.parse(faa, 'fasta'):
                record_id = str(record.id)
                targetgene = '%s__C_%s_G_%s' % (genomename,
                                                record_id.split('_')[1],
                                                record_id.split('_')[-1])
                if targetgene in BS_select_info_out:
                    description = str(record.description).replace(' ', '').split('#')
                    BS_select_info_out[targetgene] += '\t%s\t%s\n'%(description[1],description[2])
        BS_select_info_outall = []
        BS_select_info_outall.append('seq\tlocus1\tlocus2\tgenomename\tcontig\ttargetgene\tlocusgene1\tlocusgene2\n')
        for targetgene in BS_select_info_out:
            templine = BS_select_info_out[targetgene]
            if '\n' not in templine:
                templine+='\t0\t0\n'
            BS_select_info_outall.append(templine)
        f1 = open(output_file_BS + '.target.BS.info.txt', 'w')
        f1.write(''.join(BS_select_info_outall))
        f1.close()
        # extract target contigs
        for genomename in BS_select_info:
            target_contig = []
            fasta = '%s/%s/%s.origin.fasta' % (output_folder, donor_species, genomename)
            for record in SeqIO.parse(fasta, 'fasta'):
                contig = str(record.id)
                record_seq = str(record.seq)
                allBSloci = BS_select_info[genomename]
                if contig in allBSloci:
                    target_contig.append('>%s\n%s\n'%(contig,record_seq))
            f1 = open('%s/%s/%s.origin.select.fasta' % (output_folder, donor_species, genomename), 'w')
            f1.write(''.join(target_contig))
            f1.close()
    allgenome = glob.glob('%s/%s/*.origin.fasta'%(output_folder,donor_species))
    os.system('mkdir %s/%s/bwa/' % (output_folder, donor_species))
    for mut_strain_name in mut_strains:
        database = '%s/%s/%s.origin.select.fasta' % (output_folder, donor_species, mut_strain_name)
        cmds = ('#!/bin/bash\nsource ~/.bashrc\nminimap2 -d %s.mmi %s\n'%(database,database))
        vcfoutput = '%s/%s/%s.target.BS' % (output_folder, donor_species, mut_strain_name)
        try:
            f1 = open('%s.raw.vcf' % (vcfoutput), 'r')
        except FileNotFoundError:
            allsam = []
            genome_file = '%s/%s/%s.origin.fasta' % (output_folder, donor_species, mut_strain_name)
            tempbamoutput = '%s/%s/bwa/%s.to.%s' % (output_folder, donor_species, mut_strain_name,mut_strain_name)
            try:
                f1 = open('%s.sorted.bam' % (tempbamoutput), 'r')
            except FileNotFoundError:
                cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                    40, database, genome_file, 'samtools', 40,
                    tempbamoutput, 'samtools', 40, tempbamoutput, tempbamoutput, 'samtools', 40,
                    tempbamoutput)
                cmds += 'rm -r %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
            allsam.append('%s.sorted.bam' % (tempbamoutput))
            for genome_file in allgenome:
                genomename = os.path.split(genome_file)[-1].split('.origin.fasta')[0]
                if genomename not in mut_strains and genomename in pass_genome:
                    tempbamoutput = '%s/%s/bwa/%s.to.%s' % (output_folder, donor_species, genomename, mut_strain_name)
                    try:
                        f1 = open('%s.sorted.bam' % (tempbamoutput), 'r')
                    except FileNotFoundError:
                        cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                            40, database, genome_file, 'samtools', 40,
                            tempbamoutput, 'samtools', 40, tempbamoutput, tempbamoutput, 'samtools', 40,
                            tempbamoutput)
                        cmds += 'rm -r %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
                    allsam.append('%s.sorted.bam' % (tempbamoutput))
            cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
                        'bcftools', 40, database,
                        ' '.join(allsam), 'bcftools', 40, vcfoutput)
            try:
                f1 = open('%s.flt.snp.vcf' % (vcfoutput))
            except FileNotFoundError:
                cmds += '%s view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
                        'bcftools', vcfoutput, vcfoutput)
            f1 = open(vcfoutput + '.sh','w')
            f1.write(cmds)
            f1.close()

def genome_mapping_reverse(BS_loci,BS_select):
    # find BS_neighbour
    try:
        f1 = open(output_file_BS + '.target.BS.loss.info.txt', 'r')
    except FileNotFoundError:
        BS_select_info = dict()
        BS_select_info_out = dict()
        for seq in BS_select:
            sefref = False
            for allinfo in BS_loci[seq]:
                genomename, locus1, locus2, contig, targetgene = allinfo
                if not sefref:
                    BS_select_info.setdefault(genomename, [])
                    BS_select_info[genomename].append(contig)
                    sefref = True
                BS_select_info_out.setdefault(targetgene,'%s\t%s\t%s\t%s\t%s\t%s'%(seq,locus1,locus2,genomename,contig,targetgene))
        print(BS_select_info)
        # find gene locus
        for genomename in BS_select_info:
            faa = '%s/%s/%s/%s.faa'%(output_folder,BS_folder,donor_species,genomename)
            for record in SeqIO.parse(faa, 'fasta'):
                record_id = str(record.id)
                targetgene = '%s__C_%s_G_%s' % (genomename,
                                                record_id.split('_')[1],
                                                record_id.split('_')[-1])
                if targetgene in BS_select_info_out:
                    description = str(record.description).replace(' ', '').split('#')
                    BS_select_info_out[targetgene] += '\t%s\t%s\n'%(description[1],description[2])
        BS_select_info_outall = []
        BS_select_info_outall.append('seq\tlocus1\tlocus2\tgenomename\tcontig\ttargetgene\tlocusgene1\tlocusgene2\n')
        for targetgene in BS_select_info_out:
            templine = BS_select_info_out[targetgene]
            if '\n' not in templine:
                templine+='\t0\t0\n'
            BS_select_info_outall.append(templine)
        f1 = open(output_file_BS + '.target.BS.loss.info.txt', 'w')
        f1.write(''.join(BS_select_info_outall))
        f1.close()
        # extract target contigs
        for genomename in BS_select_info:
            target_contig = []
            fasta = '%s/%s/%s.origin.fasta' % (output_folder, donor_species, genomename)
            for record in SeqIO.parse(fasta, 'fasta'):
                contig = str(record.id)
                record_seq = str(record.seq)
                allBSloci = BS_select_info[genomename]
                if contig in allBSloci:
                    target_contig.append('>%s\n%s\n'%(contig,record_seq))
            f1 = open('%s/%s/%s.origin.loss.select.fasta' % (output_folder, donor_species, genomename), 'w')
            f1.write(''.join(target_contig))
            f1.close()
    if False:
        allgenome = glob.glob('%s/%s/*.origin.fasta'%(output_folder,donor_species))
        os.system('mkdir %s/%s/bwa/' % (output_folder, donor_species))
        allreference = glob.glob('%s/%s/*.origin.loss.select.fasta' % (output_folder, donor_species))
        for database in allreference:
            cmds = ('#!/bin/bash\nsource ~/.bashrc\nminimap2 -d %s.mmi %s\n' % (database, database))
            ref_name = os.path.split(database)[-1].split('.origin')[0]
            vcfoutput = '%s/%s/%s.loss.target.BS' % (output_folder, donor_species, ref_name)
            try:
                f1 = open('%s.raw.vcf' % (vcfoutput), 'r')
            except FileNotFoundError:
                    allsam = []
                    for mut_strain_name in mut_strains:
                        genome_file = '%s/%s/%s.origin.fasta' % (output_folder, donor_species, mut_strain_name)
                        tempbamoutput = '%s/%s/bwa/%s.to.%s' % (
                        output_folder, donor_species, mut_strain_name, ref_name)
                        try:
                            f1 = open('%s.sorted.bam' % (tempbamoutput), 'r')
                        except FileNotFoundError:
                            cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                                40, database, genome_file, 'samtools', 40,
                                tempbamoutput, 'samtools', 40, tempbamoutput, tempbamoutput, 'samtools', 40,
                                tempbamoutput)
                            cmds += 'rm -r %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
                        allsam.append('%s.sorted.bam' % (tempbamoutput))
                    for genome_file in allgenome:
                        genomename = os.path.split(genome_file)[-1].split('.origin.fasta')[0]
                        if genomename not in mut_strains and genomename in pass_genome:
                            tempbamoutput = '%s/%s/bwa/%s.to.%s' % (output_folder, donor_species, genomename, ref_name)
                            try:
                                f1 = open('%s.sorted.bam' % (tempbamoutput), 'r')
                            except FileNotFoundError:
                                cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                                    40, database, genome_file, 'samtools', 40,
                                    tempbamoutput, 'samtools', 40, tempbamoutput, tempbamoutput, 'samtools', 40,
                                    tempbamoutput)
                                cmds += 'rm -r %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
                            allsam.append('%s.sorted.bam' % (tempbamoutput))
                    cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
                                'bcftools', 40, database,
                                ' '.join(allsam), 'bcftools', 40, vcfoutput)
                    try:
                        f1 = open('%s.flt.snp.vcf' % (vcfoutput))
                    except FileNotFoundError:
                        cmds += '%s view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
                                'bcftools', vcfoutput, vcfoutput)
                    f1 = open(vcfoutput + '.sh','w')
                    f1.write(cmds)
                    f1.close()

def BS_sum_filter(BS_file):
    # filter BS
    BS_set = set()
    for lines in open(BS_file, 'r'):
        if not lines.startswith('BS'):
            lines_set = lines.split('\t')
            BS = lines_set[0]
            POS_on_BS = int(lines_set[4])
            if POS_on_BS >= -20 and POS_on_BS <= 20:
                # BS with changes on BS itself
                BS_set.add(BS)
    BS_file_out = []
    for lines in open(BS_file, 'r'):
        if lines.startswith('BS'):
            BS_file_out.append(lines)
        elif lines != '':
            seq = lines.split('\t')[0]
            if seq in BS_set:
                # BS with changes on BS itself
                BS_file_out.append(lines)
    f1 = open('%s.BSsum.filtered.txt' % (output_file_BS), 'w')
    f1.write(''.join(BS_file_out))
    f1.close()
    return BS_set

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
            if SNP in ['A23V']:#'A23V','G12R',
                print(SNP)
                output_file_BS = '%s/%s/%s_%s_%s' % (output_folder, BS_folder, species, donor, SNP)
                print(output_file_BS)
                select_seq_faa = dict()
                # find genome names
                vcf_file = '%s/%s%s' % (
                    vcf_folder, lineage, '.all.parsi.fasta.linktrunc.sum.txt')
                genomewithSNP = lines_set[-9].split(';')
                mut_strains, allgenome = find_strains(vcf_file, genomewithSNP)
                print('process %s mutated strains %s' % (lineage, mut_strains))
                lineage = lineage.split('.donor')[0]
                # process fimo results
                donor_species = '%s_%s' % (species, donor)
                BS_file_all = glob.glob('%s/%s/%s/*.fimo.tsv' % (output_folder, BS_folder, donor_species))
                BS_loci = dict()
                pass_genome = []
                for BS_file in BS_file_all:
                    output_file = BS_file.split('.fimo.tsv')[0]
                    input_faa = '%s/%s/%s.faa' % (output_folder, donor_species, os.path.split(output_file)[-1])
                    # load all gene position
                    Mapping_loci_all = load_genes(input_faa)
                    # load BS
                    BS_loci = load_BS(BS_file, Mapping_loci_all, BS_loci, TF_locus)
                try:
                    f1 = open('%s.BSsum.txt' % (output_file_BS), 'r')
                except FileNotFoundError:
                    allvcf_file = glob.glob('%s/%s/*.target.BS.raw.vcf' % (output_folder, donor_species))
                    allvcf_file2 = glob.glob('%s/%s/*.loss.target.BS.raw.vcf' % (output_folder, donor_species))
                    if allvcf_file == [] or allvcf_file2 == []:
                        # find targeted BS only in mutated or only in WT
                        BS_select, BS_select_wt = find_target_BS(mut_strains, BS_loci, allgenome)
                        if BS_select != set():
                            # mapping wild type to mutated strains
                            genome_mapping(BS_loci, BS_select)
                        print(BS_select_wt)
                        # if BS_select_wt!= set():
                        # mutated strains mapping to wild type
                        # genome_mapping_reverse(BS_loci, BS_select_wt)
                    else:
                        # process VCF
                        find_mutations_on_BS(allvcf_file)
        except FileNotFoundError:
            pass

alleggnog = glob.glob(os.path.join(output_folder, '*eggnog.sh'))
if alleggnog!= []:
    f1 = open(os.path.join(output_folder, 'allanno.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in alleggnog:
        f1.write('jobmit %s %s small1\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()
    print('please run %s/allanno.sh'%(output_folder))

allminimap = glob.glob(os.path.join(output_folder, '*/*target.BS.sh'))
if allminimap!= []:
    f1 = open(os.path.join(output_folder, 'alltargetBS.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in allminimap:
        f1.write('jobmit %s %s small1\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()
    print('please run %s/alltargetBS.sh'%(output_folder))
