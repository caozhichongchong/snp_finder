import os,glob
from Bio import SeqIO
import statistics
import numpy as np
from Bio.Seq import Seq
import re

output_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS/'
distance = 0

def runmapping(donor_species,mut_strains,database):
    allgenome = glob.glob('%s/%s/*.origin.fasta' % (output_folder, donor_species))
    cmds = ('#!/bin/bash\nsource ~/.bashrc\nminimap2 -d %s.mmi %s\n' % (database, database))
    cmds += 'mkdir %s/%s/bwa/\n' % (output_folder, donor_species)
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
                cmds += 'rm -r  %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
            allsam.append('%s.sorted.bam' % (tempbamoutput))
        for genome_file in allgenome:
            genomename = os.path.split(genome_file)[-1].split('.origin.fasta')[0]
            if genomename not in mut_strains:
                tempbamoutput = '%s/%s/bwa/%s.to.%s' % (output_folder, donor_species, genomename, ref_name)
                try:
                    f1 = open('%s.sorted.bam' % (tempbamoutput), 'r')
                except FileNotFoundError:
                    cmds += 'minimap2' + ' -ax asm5 -t %s %s.mmi %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
                        40, database, genome_file, 'samtools', 40,
                        tempbamoutput, 'samtools', 40, tempbamoutput, tempbamoutput, 'samtools', 40,
                        tempbamoutput)
                    cmds += 'rm -r  %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
                allsam.append('%s.sorted.bam' % (tempbamoutput))
        cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
            'bcftools', 40, database,
            ' '.join(allsam), 'bcftools', 40, vcfoutput)
        #try:
        #    f1 = open('%s.flt.snp.vcf' % (vcfoutput))
        #except FileNotFoundError:
            #cmds += '%s view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
            #    'bcftools', vcfoutput, vcfoutput)
        cmds += '#rm -r %s/%s/bwa/\n'%(output_folder, donor_species)
        f1 = open(vcfoutput + '.sh', 'w')
        f1.write(cmds)
        f1.close()
        print(vcfoutput + '.sh')

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

def find_mutations_on_BS(vcf_file,mut_strains,database):
    BS_select_info_out = dict()
    # output BS select info
    for lines in open(database.split('.origin')[0] + '.fimo.tsv', 'r'):
        if not lines.startswith('#') and not lines.startswith('motif_id') and lines != '\n':
            lines_set = lines.split('\n')[0].split('\t')
            contig, locus1, locus2 = lines_set[2:5]
            seq = lines_set[-1]
            BS_select_info_out.setdefault(contig,dict())
            BS_select_info_out[contig].setdefault(seq,[int(locus1), int(locus2),'',0])
    for lines in open(database.split('.origin')[0] + '.BS.txt', 'r'):
        if not lines.startswith('BS'):
            lines_set = lines.split('\n')[0].split('\t')
            seq, pvalue, locus, contig, strand, targetgene, locusgene = lines_set
            BS_select_info_out[contig][seq][-2]=targetgene
            BS_select_info_out[contig][seq][-1]=int(locusgene)
    BS_SNP_output = []
    BS_SNP_output.append('BS\tCHR\tPOS\tPOS_on_BS\tTargetgene\tM_num\tWT_num\tMutated_alleles\tWildtype_alleles\n')
    num_MT = len(mut_strains)
    for linesvcf in open(vcf_file, 'r'):
        linesvcf_set = linesvcf.split('\n')[0].split('\t')
        if not linesvcf.startswith('#'):
            contig, POS = linesvcf_set[0:2]
            if contig in BS_select_info_out:
                POS = int(POS)
                for seq in BS_select_info_out[contig]:
                    BSpos1, BSpos2, targetgene,locusgene = BS_select_info_out[contig][seq]
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
                        for genotype in linesvcf_set[9: 9 + num_MT]:
                            Major_ALT = allele_freq_to_allele(genotype, ALT_set)
                            Mut_allele.append(Major_ALT)
                        # Major_ALT in wild type strains
                        for genotype in linesvcf_set[9 + num_MT:]:
                            Major_ALT = allele_freq_to_allele(genotype, ALT_set)
                            Wild_allele.append(Major_ALT)
                        Mut_allele_set = sorted(set(Mut_allele))
                        Wild_allele_set = sorted(set(Wild_allele))
                        if (Wild_allele_set == ['-'] and Mut_allele_set != ['-']) or (Mut_allele_set == ['-'] and Wild_allele_set != ['-']) or (
                                any(i not in Wild_allele_set and i != '-' for i in Mut_allele_set)
                        ):
                            print(Mut_allele_set, Wild_allele_set)
                            # different alleles
                            BS_SNP_output.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tM\t%s\tWT\t%s\n' % (
                                seq, contig, POS, POS - BSpos1,targetgene,
                                num_MT,len(linesvcf_set)-num_MT-9,
                                ';'.join(list(Mut_allele_set)),
                                ';'.join(list(Wild_allele_set)),
                                '\t'.join(Mut_allele),
                                '\t'.join(Wild_allele)
                            ))
    f1 = open('%s.BSsum.txt' % (vcf_file), 'w')
    f1.write(''.join(BS_SNP_output))
    f1.close()

donor_species = 'BL_D109'
mut_strains = ['D109_BL_03','D109_BL_11','D109_BL_10','D109_BL_14','D109_BL_13','D109_BL_02']
# gain
database = '%s/%s/D109_BL_03.origin.fasta' % (output_folder, donor_species)
#runmapping(donor_species,mut_strains,database)
vcf_file = '%s/%s/D109_BL_03.loss.target.BS.raw.vcf' % (output_folder, donor_species)
print(vcf_file)
find_mutations_on_BS(vcf_file,mut_strains,database)
# gain
database = '%s/%s/D109_BL_11.origin.fasta' % (output_folder, donor_species)
#runmapping(donor_species,mut_strains,database)
vcf_file = '%s/%s/D109_BL_11.loss.target.BS.raw.vcf' % (output_folder, donor_species)
print(vcf_file)
find_mutations_on_BS(vcf_file,mut_strains,database)
# gain
database = '%s/%s/D109_BL_10.origin.fasta' % (output_folder, donor_species)
#runmapping(donor_species,mut_strains,database)
vcf_file = '%s/%s/D109_BL_10.loss.target.BS.raw.vcf' % (output_folder, donor_species)
print(vcf_file)
find_mutations_on_BS(vcf_file,mut_strains,database)
# gain
database = '%s/%s/D109_BL_14.origin.fasta' % (output_folder, donor_species)
#runmapping(donor_species,mut_strains,database)
vcf_file = '%s/%s/D109_BL_14.loss.target.BS.raw.vcf' % (output_folder, donor_species)
print(vcf_file)
find_mutations_on_BS(vcf_file,mut_strains,database)
# gain
database = '%s/%s/D109_BL_02.origin.fasta' % (output_folder, donor_species)
#runmapping(donor_species,mut_strains,database)
vcf_file = '%s/%s/D109_BL_02.loss.target.BS.raw.vcf' % (output_folder, donor_species)
print(vcf_file)
find_mutations_on_BS(vcf_file,mut_strains,database)

if True:
    # finished
    donor_species = 'BL_D109'
    mut_strains = ['D109_BL_03', 'D109_BL_11', 'D109_BL_10', 'D109_BL_14', 'D109_BL_13', 'D109_BL_02']
    # loss
    database = '%s/%s/D109_BL_01.origin.fasta' % (output_folder, donor_species)
    # runmapping(donor_species,mut_strains,database)
    vcf_file = '%s/%s/D109_BL_01.loss.target.BS.raw.vcf' % (output_folder, donor_species)
    print(vcf_file)
    find_mutations_on_BS(vcf_file,mut_strains,database)
    # gain
    database = '%s/%s/D109_BL_13.origin.fasta' % (output_folder, donor_species)
    # runmapping(donor_species,mut_strains,database)
    vcf_file = '%s/%s/D109_BL_13.loss.target.BS.raw.vcf' % (output_folder, donor_species)
    print(vcf_file)
    find_mutations_on_BS(vcf_file,mut_strains,database)
    donor_species = 'BL_D134'
    mut_strains = ['D134_BL_22', 'D134_BL_05']
    # gain
    database = '%s/%s/D134_BL_05.origin.fasta' % (output_folder, donor_species)
    # runmapping(donor_species,mut_strains,database)
    vcf_file = '%s/%s/D134_BL_05.loss.target.BS.raw.vcf' % (output_folder, donor_species)
    print(vcf_file)
    find_mutations_on_BS(vcf_file,mut_strains,database)
    # loss
    database = '%s/%s/D134_BL_01.origin.fasta' % (output_folder, donor_species)
    # runmapping(donor_species,mut_strains,database)
    vcf_file = '%s/%s/D134_BL_01.loss.target.BS.raw.vcf' % (output_folder, donor_species)
    print(vcf_file)
    find_mutations_on_BS(vcf_file,mut_strains,database)
    # gain
    database = '%s/%s/D134_BL_22.origin.fasta' % (output_folder, donor_species)
    # runmapping(donor_species,mut_strains,database)
    vcf_file = '%s/%s/D134_BL_22.loss.target.BS.raw.vcf' % (output_folder, donor_species)
    print(vcf_file)
    find_mutations_on_BS(vcf_file,mut_strains,database)
    donor_species = 'BiPs_am'
    mut_strains = ['am_BiPs_g0029']
    # loss
    database = '%s/%s/am_BiPs_g0112.origin.fasta' % (output_folder, donor_species)
    # runmapping(donor_species,mut_strains,database)
    vcf_file = '%s/%s/am_BiPs_g0112.loss.target.BS.raw.vcf' % (output_folder, donor_species)
    print(vcf_file)
    find_mutations_on_BS(vcf_file,mut_strains,database)
    # gain
    database = '%s/%s/am_BiPs_g0029.origin.fasta' % (output_folder, donor_species)
    # runmapping(donor_species,mut_strains,database)
    vcf_file = '%s/%s/am_BiPs_g0029.loss.target.BS.raw.vcf' % (output_folder, donor_species)
    print(vcf_file)
    find_mutations_on_BS(vcf_file,mut_strains,database)
if False:
    # not used
    donor_species = 'BL_bj'
    mut_strains = ['bj_BL_g0035']
    # loss
    database = '%s/%s/bj_BL_g0003.origin.fasta' % (output_folder, donor_species)
    #runmapping(donor_species,mut_strains,database)
    vcf_file = '%s/%s/bj_BL_g0003.loss.target.BS.raw.vcf' % (output_folder, donor_species)
    print(vcf_file)
    #find_mutations_on_BS(vcf_file,mut_strains,database)
    # gain
    database = '%s/%s/bj_BL_g0035.origin.fasta' % (output_folder, donor_species)
    #runmapping(donor_species,mut_strains,database)
    vcf_file = '%s/%s/bj_BL_g0035.loss.target.BS.raw.vcf' % (output_folder, donor_species)
    print(vcf_file)
    #find_mutations_on_BS(vcf_file,mut_strains,database)
