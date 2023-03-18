import glob,os
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/removerec_SNP/*.norecom.gooddepth.txt')
fastqs = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round1/*/fastq/*_1.fastq') +\
glob.glob('/scratch/users/mit_alm/gutevo/2016_09_20_Bfragilis_TS1/*/*_1.fastq')
fastas = glob.glob('/scratch/users/anniz44/genomes/donor_species/selected_species/round1/allfasta/*/fasta/*_final.scaffolds.fasta') +\
glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/moredetails/BaFr_jay/*.fasta')
fastqs2 = glob.glob('/scratch/users/mit_alm/IBD_Evo/BA/*/*/sickle2050/filter_reads_1.fastq')+\
glob.glob('/scratch/users/mit_alm/IBD_Evo/BL/*/*/sickle2050/filter_reads_1.fastq')
fastas2 = glob.glob('/scratch/users/mit_alm/IBD_Evo/BA/Assembly_for_gene_flow/*/scaffolds.fasta')+\
glob.glob('/scratch/users/mit_alm/IBD_Evo/BL/Assembly_for_gene_flow/*/scaffolds.fasta')
fastqs3 = glob.glob('/scratch/users/mit_alm/IBD_Evo/PB/*/*/sickle2050/filter_reads_1.fastq')
fastas3 = glob.glob('/scratch/users/mit_alm/IBD_Evo/PB/Assembly_for_gene_flow/*/scaffolds.fasta')

fastqs_need_assembly = glob.glob('/scratch/users/mit_alm/gutevo/2016_11_22_Bfrag_D26/*/*_1.fastq')+\
glob.glob('/scratch/users/mit_alm/gutevo/2016_11_22_Bfrag_D52/*/D*_1.fastq')+\
glob.glob('/scratch/users/mit_alm/gutevo/2016_11_22_Bfrag_D66/*/D*_1.fastq')+\
glob.glob('/scratch/users/mit_alm/gutevo/2016_11_22_Bfrag_D440/*/D*_1.fastq')+\
glob.glob('/scratch/users/mit_alm/gutevo/2016_11_22_Bfrag_D57/*/D*_1.fastq')+\
glob.glob('/scratch/users/mit_alm/gutevo/2016_11_22_Bfrag_D77/*/D*_1.fastq')+\
glob.glob('/scratch/users/mit_alm/gutevo/2016_11_24_Bfrag_D128/*/D*_1.fastq')+\
glob.glob('/scratch/users/mit_alm/gutevo/2016_11_24_Bfrag_D131/*/D*_1.fastq')+\
glob.glob('/scratch/users/mit_alm/gutevo/2016_11_24_Bfrag_D14/*/D*_1.fastq')+\
glob.glob('/scratch/users/mit_alm/gutevo/2016_11_24_Bfrag_D55/*/D*_1.fastq')+\
glob.glob('/scratch/users/mit_alm/gutevo/2016_11_24_Bfrag_D97/*/D*_1.fastq')

assembly_output = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/moredetails/BaFr_jay/'
assembly_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/BaFr_assembly/'

def genome_fastq1(fastq):
    return os.path.basename(fastq).split('_1.fastq')[0]

def genome_fasta1(fasta):
    return os.path.basename(fasta).split('_final')[0].split('.fasta')[0]

def genome_fastq2(fastq):
    return fastq.split('/')[7]

def genome_fasta2(fasta):
    return fasta.split('/')[7]

def genome_fastq3(fastq):
    fastq_temp = fastq.split('/')[7].split('_')
    return '%s_PB_%s'%(fastq_temp[0],fastq_temp[1])

def genome_fasta3(fasta):
    fasta_temp = fasta.split('/')[7].split('_')
    return '%s_PB_%s' % (fasta_temp[0], fasta_temp[1])

# check all fastq path and fasta path
allfastqs = dict()
allfastas = dict()

print(len(fastqs))
for fastq in fastqs:
    allfastqs.setdefault(genome_fastq1(fastq),fastq)
print(len(fastas))
for fasta in fastas:
    allfastas.setdefault(genome_fasta1(fasta),fasta)
print(len(fastqs2))
for fastq in fastqs2:
    allfastqs.setdefault(genome_fastq2(fastq),fastq)
print(len(fastas2))
for fasta in fastas2:
    allfastas.setdefault(genome_fasta2(fasta),fasta)
print(len(fastqs3))
for fastq in fastqs3:
    allfastqs.setdefault(genome_fastq3(fastq),fastq)
print(len(fastas3))
for fasta in fastas3:
    allfastas.setdefault(genome_fasta3(fasta),fasta)

# load all lineage information
print('processing vcfs',len(allvcffiles))
genome_lineage = dict()
for vcf in allvcffiles:
    lineage_name = os.path.basename(vcf).split('.norecom.gooddepth.txt')[0]
    for lines in open(vcf):
        if lines.startswith('CHR'):
            for genome in lines.split('\n')[0].split('\t'):
                genome_lineage.setdefault(genome,lineage_name)

# match path
allfastq_output = ['genome\tlineage\tfastq\tfasta\n']
for genome in genome_lineage:
    allfastq_output.append('%s\t%s\t%s\t%s\n'%(genome,genome_lineage[genome],allfastqs.get(genome,''),allfastas.get(genome,'')))
    if genome not in allfastqs:
        print('missing fastq infor for %s'%(genome))
    if genome not in allfastas:
        print('missing fasta infor for %s'%(genome))

def runspades_single(file1,file2,output_name):
    temp_output = output_name.split('.fasta')[0]
    cmds = '%s --careful -1 %s -2 %s -o %s --threads %s --memory 100 --cov-cutoff 7\n' % \
            ('spades.py',file1, file2, temp_output,40)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_output, output_name)
    cmds += 'rm -r %s\n' % (temp_output)
    return cmds

# fastq needs assembly
i = 7
donor_old = ''
for fastq in fastqs_need_assembly:
    genome = genome_fastq1(fastq)
    donor = genome.split('_')[0].split('N')[0]
    if donor != donor_old:
        print(donor, donor_old,fastq,i)
        i += 1
    lineage_name = 'BaFr_clustercluster%s.donor.%s'%(i,donor)
    fasta = os.path.join(assembly_output,genome + '.fasta')
    f1 = open('%s/%s.sh'%(assembly_script,genome),'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n%s\n'%(
        runspades_single(fastq,fastq.replace('_1.fastq','_2.fastq'),fasta)))
    f1.close()
    allfastq_output.append(
        '%s\t%s\t%s\t%s\n' % (genome, lineage_name, fastq, fasta))
    donor_old = donor

f1 = open('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/allgenome_path.txt','w')
f1.write(''.join(allfastq_output))
f1.close()
