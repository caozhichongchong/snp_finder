#plans
#1. same mutations across donors (simulation?)
#2. parallel change setting?
#3. simulation PE significance

################################################### SET PATH ########################################################
# pipeline #newcluster.sh
# step 1 clonal population clustering
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round1/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/'
ref_dir = '/scratch/users/anniz44/genomes/donor_species/WGS_old/WGS/vcf_round1/co-assembly'
print('python mapping_WGS_cluster.py -i %s -s %s -o %s -fa %s -ref %s '%(genome_root, input_script,
                                                                               output_dir,'.fasta.corrected.fasta',ref_dir))
print('sh %s/allWGS_cluster.sh'%(input_script))
print('python vcf_process.py  -addqual False -i %s -s %s -o %s -rd 1'%(genome_root, input_script, output_dir))
# python SNPfilter_WGS_cluster_all.py
print('python SNPfilter_WGS_cluster.py -i %s/vcf_round1/merge/ -vcf .filtered.vcf -s %'%(output_dir,input_script))
print('python cluster_CP_mapping.py -i %s -s %s -clustering 1'%(output_dir,input_script))
# check files
print('python cluster_CP_mapping.py -i %s -s %s -clustering 2'%(output_dir,input_script))
#print('python cluster_CP.py -i %s -s %s -o %s'%(pangenome_dir,input_script,output_dir))
#print('python cluster_CP.py -i %s -s %s -o %s -clustering 2'%(pangenome_dir,input_script,output_dir))

# step 2 co-assembly and mapping
# Tag inside mapping_WGS.py sh mapping_WGS.sh
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round*/'
pangenome_dir = '/scratch/users/anniz44/genomes/donor_species/vcf_round1/'
output_dir = '/scratch/users/anniz44/genomes/donor_species/'
ref_genome = '/scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/reference.genome.txt'

print('python mapping_WGS.py -i %s -s %s -o %s -fa %s -ref %s -cl %s'%(genome_root, input_script,
                                                                       output_dir,'.fasta.corrected.fasta','None', pangenome_dir))
print('sh %s/allWGS.sh'%(input_script))
# use clusters! -> vcfprocess.py
print('python vcf_process.py -i %s -s %s -o %s'%(genome_root, input_script, output_dir + '/WGS/'))

# step 3 remove recombination and bad-quality SNPs
vcf_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge'
vcf_format = '.filtered.vcf'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
# use clusters! SNPfilter_WGS_all.py
print('python SNPfilter_WGS.py -i %s -vcf %s -s %s'%(vcf_folder, vcf_format,input_script))
# check error
#'ls -lh *.err | grep -v 'K''
#'ls -lh *.all.raw.vcf | grep -v 'G' | grep -v 'M'  '
#'ls -lh *.all.flt.snp.vcf | grep -v 'K' | grep -v 'M' '
# -> move finished clonal population to finished -> run mapping_WGS again for unfinished ones
# annotate_allgenome.sh extract genome fasta and cluster ratio
# step 3.5 check hypermutators check_sub_clusters.py + rm_hypermutators.py
# step 4 calculate dnds, simulat PE and PE pathway
'grep \'reference=file\' *.raw.vcf >> all.reference.list'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
coassembly_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/all.reference.list'
cutoff_file = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/total_SNP_cutoff.txt'
core_flexible_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/allgenome/all.genome.gene.faa.uc.species.sum'
#contig = 10000
contig = 5000
print('python trunc_link.py -o %s'%(output_dir_merge))
print('python dnds.py -s %s -o %s -cutoff %s -contig %s -linktrunc True'%(input_script, output_dir_merge,cutoff_file,contig))
# simulation for PE genes within lineage
print('python PE_sim.py -snp %s/details/summary -co %s -cutoff %s'%(output_dir_merge,coassembly_file,cutoff_file))
# gene based pvalue for PE within lineage
print('python PE_gene_pvalue.py -snp %s/details/summary -co %s -core %s'%(output_dir_merge,coassembly_file,core_flexible_file))
# gene based pvalue for PE across lineages
print('python PE_gene_pvalue_across.py -snp %s/details/summary -co %s -core %s'%(output_dir_merge,coassembly_file,core_flexible_file))
# simulation for PE pathways
print('python PE_sim_pathway.py -faa %s/details/summary/all.denovo.gene.faa -pe %s/details/summary/all.selected.gene.faa.cluster.aa.all.eggnog.sum.species.sum '+\
      '-all %s/details/summary/all.denovo.gene.faa.cluster.aa.all.eggnog.sum.species.sum'%(output_dir_merge,output_dir_merge,output_dir_merge))
#extract PE genes with FDR + significance pass ->  all.selected.gene.faa.FDR.faa
print('python PE_extract.py')
# step 5 parallel evolution + PE_trunc.py + sameSNP.py
co_assembly_dir = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
cutoff_file = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/total_SNP_cutoff_lineage.txt'

print('python parallel_evolution.py -i %s -s %s -o %s -cutoff %s'%(co_assembly_dir,input_script, output_dir_merge,cutoff_file))
print('please run %s/%s'%(input_script,'allannotate.sh'))
print('please run %s/%s'%(input_script,'allannotate_all.sh'))

# step 6 dmrca
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
print('python dnds.py -s %s -o %s'%(input_script, output_dir_merge))

# sum annotation + runprokka.py + sumprokka.py + snpsum.py
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
print('python annotate_sum.py -o %s'%(output_dir_merge))

# step 7 core or flexible
# core_flexible2.py core_flexible_eggnog.py or gene core flexible in R cross-genus
# all.genome.gene.faa.uc usearch 0.9 cluster from /scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/*/*.noHM.fasta.faa
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
cutoff_file = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/total_SNP_cutoff.txt'
core_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.species.High_select2.all.gene.speciesnum.core_flexible.txt'
print('python dnds.py -s %s -o %s -cutoff %s -contig 9440 -core %s'%(input_script, output_dir_merge,cutoff_file,core_file))

# stepn co-assembly mapping to metagenomes
assembly_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/withPE/'
fastq1 = '/scratch/users/anniz44/Metagenomes/public_metagenomes/fecal_human'
output_dir = '/scratch/users/anniz44/genomes/donor_species/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/human_MG_ESSG/'

## ESSG
print('python mapping_meta.31marker.py -i %s -m %s -mfq .fasta -o %s -s %s'%(assembly_folder, fastq1, output_dir, input_script))
print('sh %s/allMGvcf.sh'%(input_script))

print('python SNPfilter_meta.31marker.py -i %s -mfq _1.fasta -o %s'%(assembly_folder,output_dir))

assembly_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/withPE/'
fastq1 = '/scratch/users/anniz44/Metagenomes/public_metagenomes/fecal_human'
output_dir = '/scratch/users/anniz44/genomes/donor_species/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/human_MG'
snp_dir = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge'
## all genomes map to MG
print('python mapping_meta.py -i %s -m %s -mfq .fasta -o %s -s %s -snp %s'%(assembly_folder, fastq1, output_dir, input_script,snp_dir))
print('sh %s/allMGvcf.sh'%(input_script))

print('python SNPsum_meta.py -o %s'%(output_dir))
################################################### END ########################################################
################################################### SET PATH ########################################################
# PE general how significant
from scipy.stats import poisson
pvalue = poisson.pmf(105,0.05*136)
print(pvalue) #2.6683837909748785e-84
################################################### END ########################################################
################################################### SET PATH ########################################################
# merge annotation and snpfiles -> mergeanno.py
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

allsnpfiles = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/*.all.parsi.fasta.sum.txt')
annotation_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.denovo.gene.faa.cluster.aa.all.eggnog.sum.species.sum'
changename = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.denovo.gene.faa.changename.txt'

genename = dict()
for lines in open(changename,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    donorspecies,nouse,oldgene,newgene = lines_set[0:4]
    genename.setdefault('%s\t%s'%(donorspecies,oldgene),newgene)

annotation_set = dict()
for lines in open(annotation_file,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    newgene,eggnog,anno,cog,cog1,cog2 = lines_set[4:10]
    annotation_set.setdefault(newgene,'%s\t%s\t%s\t%s'%(eggnog,anno[:40],cog1,cog2))

for snpfile in allsnpfiles:
    alloutput = []
    donorspecies = os.path.split(snpfile)[-1].split('.all.parsi.fasta.sum.txt')[0]
    donorspecies = donorspecies.replace('_PB_','_PaDi_')
    POS_all_0 = 0
    oldCHR = ''
    for lines in open(snpfile, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        oldgene = lines_set[5]
        if oldgene == 'Gene':
            lines = 'Newgene\tEggnog\tAnno\tCOG1\tCOG2\tPOSnew\t' + lines
        else:
            CHR = lines_set[0]
            POS = int(lines_set[1])
            if oldCHR == '':
                oldCHR = CHR
            if oldCHR!=CHR:
                try:
                    total_length = oldCHR.split('size')[1]
                except IndexError:
                    try:
                        total_length = oldCHR.split('length_')[1].split('_cov')[0]
                    except IndexError:
                        total_length = 0
                oldCHR = CHR
                POS_all_0 += int(total_length)
            newgene = genename.get('%s\t%s'%(donorspecies,oldgene),'None')
            annoall = annotation_set.get(newgene,'None\tNone\tNone\tNone')
            lines = '%s\t%s\t%s\t'%(newgene,annoall,POS+POS_all_0) + lines
        alloutput.append(lines)
    f1 = open(snpfile.replace('.all.parsi.fasta.sum.txt','.snpsum.txt'), 'w')
    f1.write(''.join(alloutput))
    f1.close()

os.system('mv /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/*snpsum.txt /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/snpsummary/')

################################################### END ########################################################
################################################### SET PATH ########################################################
# snpsum.py
# merge annotation and snpfiles IBD
################################################### END ########################################################
################################################### SET PATH ########################################################
# re cluster and remove rec
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

vcf_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/allvcfdetails/'
vcf_folderold = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/allvcfdetailsold/'
vcf_format = '.filtered.vcf'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
donor_species = '1_PB_IBD_0_clustercluster1'

#os.system('mv %s/%s*.donor.* %s/'%(vcf_folder, donor_species,vcf_folderold))

print('python SNPfilter_WGS.py -i %s -vcf %s -s %s -cluster %s'%(vcf_folder, vcf_format,input_script,donor_species))

input_script_sub = input_script + '/tree'
def run_parsi(a_parsi_file):
    if 'RAxML_parsimonyTree' not in a_parsi_file:
        donor_species = os.path.split(a_parsi_file)[-1]
        input_script_sub_temp = input_script_sub + '/' + donor_species
        try:
            os.mkdir(input_script_sub_temp)
        except IOError:
            pass
        os.system('rm -rf %s %s' % (a_parsi_file + '.out.txt',
                                    a_parsi_file + '.out.tree'))
        SNP_tree_cmd3 = ('%s\n5\nV\n1\ny\n' % (a_parsi_file))
        f1 = open(os.path.join(input_script_sub_temp, '%s.parsi.optionfile.txt'%(donor_species)), 'w')
        f1.write(SNP_tree_cmd3)
        f1.close()
        os.system('rm -rf outfile outtree')
        os.system('dnapars < %s/%s.parsi.optionfile.txt > %s/%s.parsi.output\n' % (
            input_script_sub_temp, donor_species,input_script_sub_temp, donor_species))
        os.system('mv outfile %s' % (a_parsi_file + '.out.txt'))
        os.system('mv outtree %s' % (a_parsi_file + '.out.tree'))

def outputtree_parsi(fastafile,output_file):
    SNP_alignment = dict()
    for record in SeqIO.parse(fastafile, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        SNP_alignment.setdefault(record_id,record_seq)
    SNP_alignment_output_parsi = []
    seq_num = 0
    seq_len_max = 0
    for genomename in SNP_alignment:
        seq_len = len(SNP_alignment[genomename])
        newgenomename = genomename
        if len(genomename) > 8:
            newgenomename = genomename[0:4] + '_' + genomename[-4:]
        if seq_len > 0:
            SNP_alignment_output_parsi.append('S%s    %s\n' % (newgenomename, SNP_alignment[genomename]))
            seq_num += 1
            seq_len_max = max(seq_len_max,seq_len)
    temp_line = ('   %s   %s\n' % (seq_num, seq_len_max))
    vcf_file_filtered = open(output_file, 'w')
    vcf_file_filtered.write(temp_line + ''.join(SNP_alignment_output_parsi))
    vcf_file_filtered.close()

fastafile = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/allvcfdetails/1_PB_IBD_0_clustercluster1.all_subcluster0.donor.D14.raw.vcf.filtered.vcf.final.fasta'
output_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/allvcfdetails/tree/1_PB_IBD_0_clustercluster1.all_subcluster0.donor.D14.raw.vcf.filtered.vcf.final.parsi.fasta'
outputtree_parsi(fastafile,output_file)
run_parsi(output_file)
################################################### END ########################################################
################################################### SET PATH ########################################################
# check all trees -> checkalltree.py
import glob,os
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
input_script_sub = input_script + '/tree'
os.system('rm -f %s/*'%(input_script_sub))
parsi_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
tree_folder = '%s/tree/'%(parsi_folder)

def run_parsi(a_parsi_file):
    if 'RAxML_parsimonyTree' not in a_parsi_file:
        donor_species = os.path.split(a_parsi_file)[-1]
        input_script_sub_temp = input_script_sub + '/' + donor_species
        os.system('rm -r %s' % input_script_sub_temp)
        try:
            os.mkdir(input_script_sub_temp)
        except IOError:
            pass
        SNP_tree_cmd3 = ('%s\n5\nV\n1\ny\n' % (a_parsi_file))
        f1 = open(os.path.join(input_script_sub_temp, '%s.parsi.optionfile.txt'%(donor_species)), 'w')
        f1.write(SNP_tree_cmd3)
        f1.close()
        os.system('cd %s'%(input_script_sub_temp))
        os.system('rm -r outfile outtree')
        os.system('dnapars < %s/%s.parsi.optionfile.txt > %s/%s.parsi.output\n' % (
            input_script_sub_temp, donor_species,input_script_sub_temp, donor_species))
        os.system('mv outfile %s' % (a_parsi_file + '.out.txt'))
        os.system('mv outtree %s' % (a_parsi_file + '.out.tree'))

for fastafile in glob.glob('%s/*.all.parsi.fasta'%(parsi_folder)):
    output_file = '%s/%s.out.tree'%(tree_folder,os.path.split(fastafile)[-1])
    #print(fastafile)
    filesize = 0
    try:
        filesize = int(os.path.getsize(output_file))
    except FileNotFoundError:
        pass
    if filesize == 0:
        print('run tree %s'%(fastafile))
        run_parsi(fastafile)

os.system('mv %s/*.parsi.fasta.out* %s/'%(parsi_folder,tree_folder))
os.system('python treetopdf.py')
os.system('sh treetopdf.sh')
################################################### END ########################################################
################################################### SET PATH ########################################################
# merge MG results
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/withinHS'
alloutput = []
os.system('rm %s'%('%s/all.withinHS.snp.sum'%(outputdir)))
for files in glob.glob('%s/*withinHS.snp.sum'%(outputdir)):
    donor = os.path.split(files)[-1].split('.withinHS.snp.sum')[0]
    for lines in open(files, 'r'):
        alloutput.append('%s\t' % (donor) + lines)

f1 = open('%s/all.withinHS.snp.sum'%(outputdir),'w')
f1.write(''.join(alloutput))
f1.close()

outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/other'
alloutput = []
os.system('rm %s'%('%s/all.other.snp.sum'%(outputdir)))
for files in glob.glob('%s/*other.snp.sum'%(outputdir)):
    donor = os.path.split(files)[-1].split('.other.snp.sum')[0]
    for lines in open(files, 'r'):
        alloutput.append('%s\t' % (donor) + lines)

f1 = open('%s/all.other.snp.sum'%(outputdir),'w')
f1.write(''.join(alloutput))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# vcfprocess.py
# vcf processing per cluster
import glob
import os
input_script_vcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcfpro_round2'
output_dir = '/scratch/users/anniz44/genomes/donor_species/'
genome_root = '/scratch/users/anniz44/genomes/donor_species/selected_species/round*/'

os.system('rm -rf %s'%(input_script_vcf))
vcf_name = '.all*raw.vcf'
try:
    os.mkdir(input_script_vcf)
except IOError:
    pass

all_vcf_file=glob.glob(os.path.join(output_dir + '/vcf_round2/merge/','*%s'%(vcf_name)))
print(all_vcf_file)
for vcf_file in all_vcf_file:
    filesize = 0
    try:
        filesize = int(os.path.getsize(vcf_file + '.filtered.snpfreq.txt'))
        if filesize == 0:
            donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
            f1 = open(os.path.join(input_script_vcf, '%s.sh' % (donor_species)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n')
            f1.write('python %s/../vcf_process.py -i \"%s\" -s %s -o %s -cluster %s\n' % (
                input_script_vcf, genome_root, input_script_vcf, output_dir, donor_species))
            f1.close()
    except FileNotFoundError:
        donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
        f1 = open(os.path.join(input_script_vcf, '%s.sh'%(donor_species)), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n')
        f1.write('python %s/../vcf_process.py -i \"%s\" -s %s -o %s -cluster %s\n'%(input_script_vcf,genome_root, input_script_vcf, output_dir,donor_species))
        f1.close()

f1 = open(os.path.join(input_script_vcf, '../allvcfprocessing.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run: sh %s/../allvcfprocessing.sh'%(input_script_vcf))

################################################### END ########################################################
################################################### SET PATH ########################################################
# remove rec per cluster
import glob
import os
vcf_format = '.filtered.vcf'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
vcf_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge'

input_script_vcf = os.path.join(input_script,'vcfremoverec')

rec_set = ['EsCo','PB','BL','BA','TuSa','BiPs','PaDi']
os.system('rm -rf %s'%(input_script_vcf))
try:
    os.mkdir(input_script_vcf)
except IOError:
    pass

all_vcf_file=glob.glob(os.path.join(vcf_folder,'*%s'%(vcf_format)))
Donor_species = set()
for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
    if donor_species not in Donor_species:
        Donor_species.add(donor_species)
        parsi_output = os.path.join(os.path.split(vcf_file)[0], '%s.all.parsi.fasta' % (donor_species))
        filesize = 0
        try:
            filesize = int(os.path.getsize(parsi_output))
        except FileNotFoundError:
            pass
        if filesize == 0:
            if any(item in donor_species for item in rec_set):
                rec_cutoff = 50000
                contig_cutoff = 10000
            else:
                rec_cutoff = 5000
                contig_cutoff = 5000
            if donor_species == '1_PB_IBD_0_clustercluster1':
                # lots of big MGEs, small contigs are main contigs
                contig_cutoff = 5000
            f1 = open(os.path.join(input_script_vcf, '%s.sh'%(donor_species)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n')
            f1.write('python %s/SNPfilter_WGS_withdonor.py -i %s -vcf %s -cluster %s -s %s -rec %s -contig %s\n'%(input_script,vcf_folder, vcf_format, donor_species,input_script,rec_cutoff,contig_cutoff))
            f1.close()

f1 = open(os.path.join(input_script, 'allvcfremoverec.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    #f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.write('sh %s\n' % (sub_scripts))

f1.write('jobmit dnds.sh dnds.sh\n')
f1.close()
print('please run: sh %s/allvcfremoverec.sh'%(input_script))

################################################### END ########################################################
################################################### SET PATH ########################################################
# remove rec per cluster -> SNPfilter_WGS_all.py
import glob
import os
vcf_format = '.filtered.vcf'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
vcf_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge'
input_script_vcf = os.path.join(input_script,'vcfremoverec')
rec_cutoff = 5000
contig_cutoff = 10000

#rec_set = ['BL','BA','PaDi']
rec_set = []
os.system('rm -rf %s'%(input_script_vcf))
try:
    os.mkdir(input_script_vcf)
except IOError:
    pass

all_vcf_file=glob.glob(os.path.join(vcf_folder,'*%s'%(vcf_format)))
Donor_species = set()

for vcf_file in all_vcf_file:
    donor_species = os.path.split(vcf_file)[-1].split('.all')[0]
    if donor_species not in Donor_species:
        Donor_species.add(donor_species)
        parsi_output = os.path.join(os.path.split(vcf_file)[0], '%s.all.parsi.fasta' % (donor_species))
        filesize = 0
        try:
            filesize = int(os.path.getsize(parsi_output))
        except FileNotFoundError:
            pass
        if filesize == 0:
            f1 = open(os.path.join(input_script_vcf, '%s.sh'%(donor_species)), 'w')
            f1.write('#!/bin/bash\nsource ~/.bashrc\npy37\n')
            f1.write('python %s/SNPfilter_WGS.py -i %s -vcf %s -cluster %s -s %s -rec %s -contig %s\n'%(input_script,vcf_folder, vcf_format, donor_species,input_script,rec_cutoff,contig_cutoff))
            #f1.write('python %s/SNPfilter_WGS_cluster.py -i %s -vcf %s -cluster %s -s %s\n'%(input_script,vcf_folder, vcf_format, donor_species,input_script))
            f1.close()

f1 = open(os.path.join(input_script, 'allvcfremoverec.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    #f1.write('sh %s\n' % (sub_scripts))

f1.write('#jobmit dnds.sh dnds.sh\n')
f1.close()
print('please run: sh %s/allvcfremoverec.sh'%(input_script))

################################################### END ########################################################
################################################### SET PATH ########################################################

# check clonal population
from Bio import SeqIO
from Bio.Seq import Seq
import os
fasta = 'core_gene_alignment.aln'
os.system('prodigal -q -i %s -d %s.fna' % (fasta, fasta))
output = []
for record in SeqIO.parse(fasta, 'fasta'):
    record_id = str(record.id)
    record_seq = str(record.seq)
    output.append('>%s\n%s\n' % (record_id, record_seq))
    break

f1 = open(fasta + '.subset', 'w')
f1.write(''.join(output))
f1.close()
# check SNP position
#os.system('snp-sites -v %s.subset > %s.subset.vcf'%(output,output))
# check homologous in core genome
#os.system('blastn -subject %s.subset  -query %s.subset -out %s.subset.homologous.txt -outfmt 6 -max_target_seqs 100'%(output,output,output))
# check gene level SNP distribution
os.system('prodigal -q -i %s.subset -d %s.subset.fna' % (fasta, fasta))
#os.system('blastn -subject %s  -query %s -out %s.gene.divergence.txt -outfmt 6 -window_size 500 -perc_identity 80 -max_target_seqs 100'%(fasta,fasta,fasta))
os.system('%s -makeudb_usearch %s.subset.fna -output %s.subset.fna.udb' %
          ('usearch', fasta, fasta))
os.system('%s -ublast %s.fna -db %s.subset.fna.udb  -evalue 1e-2 -accel 0.5 -blast6out %s -threads 5 -strand plus -maxaccepts 500' %
          ('usearch', fasta, fasta, fasta + '.gene.divergence.txt'))

# check gene level SNP distribution of a suspicious subset
fasta = 'core_gene_alignment.aln'
subset = ['ao_BiPs_g0003','ao_BiPs_g0006']
output = []
for record in SeqIO.parse(fasta, 'fasta'):
    record_id = str(record.id)
    if record_id in subset:
        record_seq = str(record.seq)
        output.append('>%s\n%s\n' % (record_id, record_seq))

f1 = open(fasta + '.check.subset', 'w')
f1.write(''.join(output))
f1.close()
os.system('prodigal -q -i %s.check.subset -d %s.check.subset.fna' % (fasta, fasta))
os.system('%s -ublast %s.check.subset.fna -db %s.subset.fna.udb  -evalue 1e-2 -accel 0.5 -blast6out %s -threads 5 -strand plus -maxaccepts 500' %
          ('usearch', fasta, fasta, fasta + '.check.gene.divergence.txt'))
################################################### END ########################################################
################################################### SET PATH ########################################################
# set cutoff for recombination windows to calculate NS ratio
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio.Phylo import BaseTree
import statistics

input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly'
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge'
vcf_name = '.raw.vcf.filtered.snp.txt'
windowsize_set = [810000, 729000, 656100, 590490, 531441, 478296, 430466, 387419, 348677, 313809,
                  282428, 254185, 228766, 205889, 185300, 166770, 150093, 135083, 121574, 109416,
                  98474, 88626, 79763, 71786, 64607, 58146, 52331, 47097, 42387, 38148, 34333, 30899,
                  27809, 25028, 22525, 20272, 18244, 16419, 14777, 13299, 11969, 10772, 9694, 8724,
                  7851, 7065, 6358, 5722, 5149, 4634, 4170, 3753, 3377, 3039, 2735, 2461, 2214, 1992,
                  1792, 1612, 1450, 1305, 1174, 1056, 950, 855, 769, 692, 622, 559, 503, 452, 406, 365,
                  328, 295, 265, 238, 214, 192, 172, 154, 138, 124, 111, 99, 89, 80, 72, 64, 57, 51, 45,
                  40, 36, 32, 28, 25, 22, 19, 17, 15, 13, 11, 9, 8, 7, 6, 5, 4, 3, 2, 1]

try:
    os.mkdir(output_dir_merge + '/summary')
except IOError:
    pass

def windowing_concat(Seq_N):
    for windowsize in windowsize_set:
        N_S_set = [0, 0, 0]
        total_NS = ['']
        POS_old = 0
        for CHR in Seq_N:
            # calculate interval
            try:
                total_length = CHR.split('size')[1]
            except IndexError:
                total_length = CHR.split('length_')[1].split('_cov')[0]
            total_length = int(total_length)
            if POS_old + total_length >= windowsize:
                total_interval = int(total_length / windowsize) + 1
                total_NS += [''] * (total_interval)
            # windowing SNPs
            POS_set = Seq_N[CHR][0]
            NS_set = Seq_N[CHR][1]
            for i in range(0, len(POS_set)):
                POS = POS_set[i] + POS_old
                loci_POS = int(POS / windowsize)
                if total_NS[loci_POS] == '':
                    total_NS[loci_POS] = NS_set[i]
            POS_old += total_length
        N_S_set[0] += total_NS.count('N')
        N_S_set[1] += total_NS.count('S')
        try:
            N_S_set[2] = N_S_set[0] / N_S_set[1]
        except ZeroDivisionError:
            N_S_set[2] = 'N_only'
        Output.append('%s\t%s\t%s\t%s\t%s\t\n'%(donor_species,windowsize,
                                                N_S_set[0],N_S_set[1],N_S_set[2]))

def windowing(Seq_N):
    for windowsize in windowsize_set:
        N_S_set = [0, 0, 0]
        total_NS = ['']
        for CHR in Seq_N:
            # calculate interval
            try:
                total_length = CHR.split('size')[1]
            except IndexError:
                total_length = CHR.split('length_')[1].split('_cov')[0]
            total_length = int(total_length)
            if total_length >= windowsize:
                total_interval = int(total_length / windowsize) + 1
                total_NS += [''] * (total_interval)
            # windowing SNPs
            POS_set = Seq_N[CHR][0]
            NS_set = Seq_N[CHR][1]
            for i in range(0, len(POS_set)):
                POS = POS_set[i]
                loci_POS = int(POS / windowsize)
                if total_NS[loci_POS] == '':
                    total_NS[loci_POS] = NS_set[i]
        N_S_set[0] += total_NS.count('N')
        N_S_set[1] += total_NS.count('S')
        try:
            N_S_set[2] = N_S_set[0] / N_S_set[1]
        except ZeroDivisionError:
            N_S_set[2] = 'N_only'
        Output.append('%s\t%s\t%s\t%s\t%s\t\n'%(donor_species,windowsize,
                                                N_S_set[0],N_S_set[1],N_S_set[2]))

def readSNPfile(vcf_file):
    Seq_N = dict()
    for lines in open(vcf_file,'r'):
        lines_set = lines.replace('\n', '').replace('\r', '').split('\t')
        CHR = lines_set[0]
        POS = int(lines_set[1])
        N_S = lines_set[-2]
        if 'None' not in N_S and 'S' not in N_S:
            N_S = 'N'
        Seq_N.setdefault(CHR,[[],[]])
        Seq_N[CHR][0].append(POS)
        Seq_N[CHR][1].append(N_S)
    return Seq_N

# before remove rec
Output = []
all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
for vcf_file in all_vcf_file:
    if 'BaFr_clustercluster2' not in vcf_file:
        print(vcf_file)
        donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0]
        Seq_N = readSNPfile(vcf_file)
        windowing(Seq_N)

foutput = open(output_dir_merge + '/summary/all.donor.species.NSratio.txt', 'w')
foutput.write('donor_species\twindowsize\tN\tS\tNS_ratio\t\n')
foutput.write(''.join(Output))
foutput.close()

# after remove rec
vcf_name = '.raw.vcf.filtered.vcf.final.snp.txt'
Output = []
all_vcf_file = glob.glob(os.path.join(output_dir_merge, '*%s' % (vcf_name)))
for vcf_file in all_vcf_file:
    if 'BaFr_clustercluster2' not in vcf_file:
        print(vcf_file)
        donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0]
        Seq_N = readSNPfile(vcf_file)
        windowing(Seq_N)

foutput = open(output_dir_merge + '/summary/all.donor.species.NSratio.removerec.txt', 'w')
foutput.write('donor_species\twindowsize\tN\tS\tNS_ratio\t\n')
foutput.write(''.join(Output))
foutput.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# cluster ratio -> cluster_ratio.py
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

species_list = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/species_order.new.txt'
denovo_fasta = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.genome.gene.faa'
output_list = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/cluster_ratio.txt'

def split_line(lines):
    return lines.split('\n')[0].split('\t')

def find_species(record_id):
    species = record_id.split('_')[0]
    if species == '1':
        species = record_id.split('_')[1]
    return species

# load all species
All_species = dict()
for lines in open(species_list,'r'):
    if not lines.startswith('species_short'):
        lines_set = split_line(lines)
        All_species.setdefault(lines_set[2],[0,0])

# load denovo fasta
for record in SeqIO.parse(denovo_fasta, 'fasta'):
    record_id = str(record.id)
    All_species[find_species(record_id)][0]+=1

# load clustered denovo fasta
for record in SeqIO.parse(denovo_fasta + '.cluster.aa', 'fasta'):
    record_id = str(record.id)
    All_species[find_species(record_id)][1] += 1

# output cluster ratio
Output = []
Output.append('species\tseq_before_cluster\tseq_after_cluster\tcluster_ratio\t\n')
for species in All_species:
    before_cluster,aftercluster = All_species[species]
    cluster_ratio = 0
    if aftercluster > 0:
        cluster_ratio = before_cluster/aftercluster
    Output.append('%s\t%s\t%s\t%s\t\n'%(species,before_cluster,aftercluster,cluster_ratio))

foutput = open(output_list , 'w')
foutput.write(''.join(Output))
foutput.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# annotate all genomes -> annotate_allgenome.py, annotate_allgenome.sh
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

allreference_list = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/all.reference.list'
#allfasta = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/*/*.noHM.fasta.faa')
output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'
input_script_sub = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/annotate_genome'

#os.system('head %s/*.all*.raw.vcf | grep \'##reference=\' > %s/all.reference.list'%(output_dir_merge,output_dir_merge))
try:
    os.mkdir(input_script_sub)
except IOError:
    pass

def annotation(all_filter_gene_fasta_file,pre_cluster = ''):
    all_filter_gene_fasta_file = all_filter_gene_fasta_file + '.cluster.aa'
    if pre_cluster!= '':
        os.system('#%s -makeudb_usearch %s -output %s.udb' %
                  ('usearch', pre_cluster, pre_cluster))
        os.system('%s -ublast %s -db %s.udb  -evalue 1e-2 -accel 0.5 -blast6out %s -threads 2'%
                  ('usearch', all_filter_gene_fasta_file,pre_cluster, all_filter_gene_fasta_file + '.ref.out.txt'))
    # run metacyc
    cutoff = 50
    cutoff2 = 80
    database = '/scratch/users/mit_alm/database/metacyc/protseq.fsa'
    cmds = ("%s blastp --query %s --db %s.dmnd --out %s.metacyc.txt --id %s --query-cover %s --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 40\n"
            %('diamond',all_filter_gene_fasta_file,database,all_filter_gene_fasta_file,cutoff,cutoff2))
    f1 = open(os.path.join(input_script_sub, 'metacyc.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    # run eggnog
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/eggnog/xaa.hmm'
    cmds = ('%s --tblout %s.eggnog.1.txt --cpu 40 -E %s %s %s\n') %('hmmsearch', all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.1.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xab.hmm'
    cmds = ('%s --tblout %s.eggnog.2.txt --cpu 40 -E %s %s %s\n') % (
        'hmmsearch',
        all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.2.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmds))
    f1.close()
    database = '/scratch/users/mit_alm/database/eggnog/xac.hmm'
    cmds = ('%s --tblout %s.eggnog.3.txt --cpu 40 -E %s %s %s\n') % (
        'hmmsearch',
        all_filter_gene_fasta_file, cutoff, database, all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'eggnog.3.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s' % (cmds))
    f1.close()
    # run kegg
    cutoff = 0.01
    database = '/scratch/users/mit_alm/database/kegg/kofam/profiles/prokaryote/prokaryote.hmm'
    cmds = ('%s --tblout %s.kegg.txt --cpu 40 -E %s %s %s\n') %('hmmsearch', all_filter_gene_fasta_file,cutoff,database,all_filter_gene_fasta_file)
    f1 = open(os.path.join(input_script_sub, 'kegg.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n%s'%(cmds))
    f1.close()
    # run prokka
    cmdsprokka = 'py37\nexport LD_LIBRARY_PATH=/scratch/users/anniz44/bin/pro/lib/gsl-2.6:/scratch/users/anniz44/bin/pro/lib/glibc-2.14-build:/scratch/users/anniz44/bin/pro/lib/:/scratch/users/anniz44/bin/miniconda3/lib:$LD_LIBRARY_PATH\n' + \
                 'prokka --kingdom Bacteria --force --outdir %s/prokka_%s  --protein %s --locustag Bacter %s/%s\n' % \
                 (output_dir_merge + '/summary', os.path.split(all_filter_gene_fasta_file)[-1],
                  all_filter_gene_fasta_file,
                  output_dir_merge + '/summary',
                  os.path.split(all_filter_gene_fasta_file)[-1].replace('.faa', '.fna'))
    f1 = open(os.path.join(input_script_sub, 'prokka.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (cmdsprokka))
    f1.close()
    # all scripts
    f1 = open(os.path.join(input_script, 'allannotategenome.sh'), 'w')
    f1.write('#!/bin/bash\nsource ~/.bashrc\n')
    for sub_scripts in glob.glob(os.path.join(input_script_sub, '*.sh')):
        f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))
    f1.close()
    print('please run %s/%s'%(input_script,'allannotategenome.sh'))

allfasta = []
for lines in open(allreference_list, 'r'):
    if lines.startswith('##reference=file:'):
        # set database
        database_file = lines.split('##reference=file:')[1].split('\n')[0]
        allfasta.append(database_file)
        input_fasta = database_file + '.faa'
        input_fasta_dna = database_file + '.fna'
        try:
            f1 = open(input_fasta, 'r')
        except FileNotFoundError:
            os.system('prodigal -q -i %s -a %s' % (database_file, input_fasta))
        try:
            f1 = open(input_fasta_dna, 'r')
        except FileNotFoundError:
            os.system('prodigal -q -i %s -d %s' % (database_file, input_fasta_dna))

# merge all genes
output_gene = output_dir_merge + '/summary/all.genome.gene.faa'
output_fasta = []
try:
    f1 = open(output_gene,'r')
except IOError:
    change_name = set()
    for database in allfasta:
        database += '.faa'
        donor_species = os.path.split(database)[-1].split('.all')[0]
        donor_species_new = donor_species.replace('_clustercluster', '_CL')
        print(donor_species,donor_species_new)
        for record in SeqIO.parse(database, 'fasta'):
            record_id = str(record.id)
            temp_line = '%s\t%s\t%s\t' % (donor_species, donor_species_new, record_id)
            record_id = 'C_%s_G_%s' % (record_id.split('_')[1], record_id.split('_')[-1])
            output_fasta.append('>%s__%s\n%s\n' % (donor_species_new, record_id, str(record.seq)))
            temp_line += '%s__%s\t\n' % (donor_species_new, record_id)
            change_name.add(temp_line)
    f1 = open(output_gene, 'w')
    f1.write(''.join(output_fasta))
    f1.close()
    f1 = open(output_gene + '.changename.txt', 'w')
    f1.write(''.join(list(change_name)))
    f1.close()
    # run cluster
    cutoff = 0.9
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                   % ('usearch', output_fasta, cutoff, output_fasta,
                      output_fasta, 40))
    os.system(cmd_cluster)

annotation(output_gene)

output_gene = output_dir_merge + '/summary/all.genome.gene.fna'
output_fasta = []
try:
    f1 = open(output_gene,'r')
except IOError:
    change_name = set()
    for database in allfasta:
        database += '.fna'
        donor_species = os.path.split(database)[-1].split('.all')[0]
        donor_species_new = donor_species.replace('_clustercluster', '_CL')
        print(donor_species,donor_species_new)
        for record in SeqIO.parse(database, 'fasta'):
            record_id = str(record.id)
            temp_line = '%s\t%s\t%s\t' % (donor_species, donor_species_new, record_id)
            record_id = 'C_%s_G_%s' % (record_id.split('_')[1], record_id.split('_')[-1])
            output_fasta.append('>%s__%s\n%s\n' % (donor_species_new, record_id, str(record.seq)))
            temp_line += '%s__%s\t\n' % (donor_species_new, record_id)
            change_name.add(temp_line)
    f1 = open(output_gene, 'w')
    f1.write(''.join(output_fasta))
    f1.close()
    cutoff = 0.9
    cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.cluster.aa -uc %s.uc -threads %s\n'
                   % ('usearch', output_fasta, cutoff, output_fasta,
                      output_fasta, 40))
    os.system(cmd_cluster)

################################################### END ########################################################
################################################### SET PATH ########################################################
# core flexible of de novo genes, split HS annotation within lineage and across lineages -> core_flexible2.py
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'

allgenome = output_dir_merge + '/summary/all.genome.gene.faa.uc'
allgenome_denovo = output_dir_merge + '/summary/all.denovo.gene.faa.uc'
allgenome_HS = output_dir_merge + '/summary/all.selected.gene.faa.uc'

def cluster_uc(cluster_input):
    Clusters = dict()
    for lines in open(cluster_input, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        cluster = line_set[1]
        record_name = line_set[8].split(' ')[0]
        Clusters.setdefault(cluster, set())
        Clusters[cluster].add(record_name)
    return Clusters

# clusters based on all genes
Clusters = cluster_uc(allgenome)

def genus_list(species_list):
    return set([i.replace('BL','BiLo').replace('BA','BiAd')[0:2] for i in species_list])

def species_cluster(Clusters):
    Clusters_species = dict()
    Gene_cluster = dict()
    Output = []
    Output.append('record\tspecies_num\tallspecies\tgenus_num\t\n')
    for cluster in Clusters:
        allrecord = Clusters[cluster]
        Clusters_species.setdefault(cluster,set())
        for record in allrecord:
            record_new = record.split('_')[0]
            if record_new == '1':
                record_new = record.split('_')[1]
            record_new = record_new.replace('PB','PaDi').replace('Bfragilis','BaFr')
            Clusters_species[cluster].add(record_new)
        for record in allrecord:
            Gene_cluster.setdefault(record,Clusters_species[cluster])
    for record in Gene_cluster:
        allspecies = Gene_cluster[record]
        genus = genus_list(allspecies)
        Output.append('%s\t%s\t%s\t%s\t\n'%(record,len(allspecies),','.join(list(allspecies)),len(genus)
                                        ))
    f1 = open(allgenome + '.species.sum','w')
    f1.write(''.join(Output))
    f1.close()
    return [Clusters_species,Gene_cluster]

Clusters_species,Gene_cluster = species_cluster(Clusters)

def cluster_species_sum(fastaname,Gene_cluster):
    Clusters = cluster_uc(fastaname)
    Output = []
    Output.append('record\tspecies_num\tallspecies\tgenus_num\t\n')
    for cluster in Clusters:
        allrecord = Clusters[cluster]
        for record in allrecord:
            record_new = record.replace('.donor.'+record.split('.')[-1].split('__')[0],'')
            if '1_PaDi_IBD' in record_new:
                record_new = record_new.replace('1_PaDi_IBD','1_PB_IBD')
            allspecies = Gene_cluster[record_new]
            genus = genus_list(allspecies)
            Output.append('%s\t%s\t%s\t%s\t\n' % (record, len(allspecies), ','.join(list(allspecies)),len(genus)
                                                  ))
    f1 = open(fastaname + '.species.sum', 'w')
    f1.write(''.join(Output))
    f1.close()

cluster_species_sum(allgenome_denovo,Gene_cluster)
cluster_species_sum(allgenome_HS,Gene_cluster)
################################################### END ########################################################
################################################### SET PATH ########################################################
# core flexible of de novo genes eggnog ->core_flexible_eggnog.py
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq

output_dir_merge = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details'
genus_cutoff = 2 # more than 2 genera as core
allgenome = output_dir_merge + '/summary/all.genome.gene.faa'
allgenome_denovo = output_dir_merge + '/summary/all.denovo.gene.faa'
allgenome_HS = output_dir_merge + '/summary/all.selected.gene.faa'

def genus_list(species_list):
    return set([i.replace('BL','BiLo').replace('BA','BiAd')[0:2] for i in species_list])

def function_species(fastaname):
    species_file = fastaname + '.uc.species.sum'
    eggnog_file = fastaname + '.cluster.aa.all.eggnog.sum'
    Species_fun = dict()
    for lines in open(species_file, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        record_name = line_set[0]
        species_num = line_set[1]
        species_name = line_set[2].split(',')
        tag_genus = 'False'
        if len(species_name) > 1:
            genus = genus_list(species_name)
            if len(genus) > genus_cutoff:
                tag_genus = 'True'
        Species_fun.setdefault(record_name, [species_num,tag_genus])
    Output = []
    for lines in open(eggnog_file, 'r'):
        line_set = lines.split('\n')[0].split('\t')
        if lines.startswith('cluster'):
            Output.append('\t'.join(line_set) + '\tspecies_num\tcross_genus\n')
        else:
            if line_set[6] == '':
                line_set.pop(6)
            record_name = line_set[4]
            species_num, tag_genus = Species_fun[record_name]
            Output.append('\t'.join(line_set) + '\t%s\t%s\n'%(species_num, tag_genus))
    f1 = open(fastaname + '.cluster.aa.all.eggnog.sum.species.sum', 'w')
    f1.write(''.join(Output))
    f1.close()

function_species(allgenome_HS)
function_species(allgenome_denovo)
function_species(allgenome)

# split HS annotation within lineage and across lineages
def HS_lineage(filename,HS_lineagefasta):
    HS_lineage_set = set()
    for record in SeqIO.parse(HS_lineagefasta, 'fasta'):
        record_id = str(record.id)
        HS_lineage_set.add(record_id)
    Output_within = set()
    Output_across = set()
    for lines in open(filename,'r'):
        if lines.split('\t')[4] in HS_lineage_set:
            Output_within.add(lines)
        else:
            Output_across.add(lines)
    f1 = open(filename + '.within.sum', 'w')
    f1.write(''.join(list(Output_within)))
    f1.close()
    f1 = open(filename + '.across.sum', 'w')
    f1.write(''.join(list(Output_across)))
    f1.close()

HS_annotation_sum = output_dir_merge + '/summary/all.selected.gene.faa.cluster.aa.all.eggnog.sum.species.sum'
allgenome_HS_lineage = output_dir_merge + '/summary/all.selected.gene.faa'

HS_lineage(HS_annotation_sum,allgenome_HS_lineage)

################################################### END ########################################################
################################################### SET PATH ########################################################
# tree phylogenetic diversity
from Bio import Phylo
from Bio.Phylo import BaseTree
import numpy

def to_distance(tree):
    # phylogenetic diversity of a tree, sum of branch_length
    dis_all = 0
    for parent in tree.find_clades(terminal=False, order="level"):
        for child in parent.clades:
            if child.branch_length:
                dis_all += child.branch_length
    return dis_all

filename_tree = '1_BL_IBD_0_clustercluster1.donor.D77.all.parsi.fasta.out.tree'
tree = Phylo.read(filename_tree, "newick")
to_distance(tree)

################################################### END ########################################################
################################################### SET PATH ########################################################
# move MG results
import os,glob
result_dir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
snpsumdir = '/scratch/users/anniz44/genomes/donor_species/MG/snpsum/'
covsumdir = '/scratch/users/anniz44/genomes/donor_species/MG/covsum/'

try:
    os.mkdir(snpsumdir)
except IOError:
    pass

try:
    os.mkdir(covsumdir)
except IOError:
    pass

for snpfile in glob.glob('%s/*.snp.sum'%(result_dir)):
    try:
        outputdir = os.path.split(snpfile)[-1].split('.IN.')[1].split('.snp.sum')[0]
        outputdir = os.path.join('%s/%s'%(snpsumdir,outputdir))
        os.mkdir(outputdir)
    except IOError:
        pass
    os.system('mv %s %s/'%(snpfile,outputdir))


for snpfile in glob.glob('%s/*.coverage.sum'%(result_dir)) + glob.glob('%s/*.depth.MV.sum'%(result_dir)) + glob.glob('%s/*.depth.sum'%(result_dir)):
    try:
        outputdir = os.path.split(snpfile)[-1].split('.fasta.')[1].split('.raw.vcf')[0]
        outputdir = os.path.join('%s/%s'%(covsumdir,outputdir))
        os.mkdir(outputdir)
    except IOError:
        pass
    os.system('mv %s %s/'%(snpfile,outputdir))



for snpfile in glob.glob('%s/finished/*.zip'%(result_dir)):
    try:
        outputdir = os.path.split(snpfile)[-1].split('.fasta.')[1].split('.raw.vcf')[0]
        outputdir = os.path.join('%s/finished/%s'%(result_dir,outputdir))
        os.mkdir(outputdir)
    except IOError:
        pass
    os.system('mv %s %s/'%(snpfile,outputdir))
################################################### END ########################################################
################################################### SET PATH ########################################################
# PE truncation -> PE_trunc.py
# all SNPs truncation
import os,glob
snp_dir = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'

def loadHS(sum_file):
    HS_gene = dict()
    for lines in open(sum_file, 'r'):
        if not lines.startswith("#donor_species"):
            lines_set = lines.split('\n')[0].split('\t')
            if lines_set[-1] == 'True':
                lineage, gene = lines_set[0:2]
                HS_gene.setdefault(lineage,set())
                HS_gene[lineage].add(gene)
    return HS_gene

def load_snp(snpfile,alloutput):
    donor_species = os.path.split(snpfile)[-1].split('.raw.vcf.filtered.vcf.final.snp.txt')[0].replace('.all','')
    #if donor_species in HS_gene_within:
    for lines in open(snpfile,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            genename, POS,N_or_S,AAchange = lines_set[-4:]
            genenamenew, POS = lines_set[:2]
            if '*' in AAchange:
                if '*' == AAchange[0]:
                    Trunc_SNP = lines_set[2]
                else:
                    Trunc_SNP = lines_set[3]
                if donor_species in HS_gene_within and genename in HS_gene_within[donor_species]:
                    alloutput.add('%s\t%s\t%s\t%s\t%s\twithinHS\t%s\t%s\n'%(donor_species,genenamenew,POS,AAchange[0],AAchange[1],Trunc_SNP,genename))
                #elif genename in HS_gene_all[donor_species]:
                    #alloutput.add('%s\t%s\t%s\t%s\t%s\tacrossHS\t%s\t%s\n'%(donor_species,genenamenew,POS,AAchange[0],AAchange[1],Trunc_SNP,genename))
                else:
                    alloutput.add('%s\t%s\t%s\t%s\t%s\tothers\t%s\t%s\n' % (
                    donor_species, genenamenew, POS, AAchange[0], AAchange[1],Trunc_SNP,genename))
    else:
        print('no HS in %s'%(donor_species))
    return alloutput

# load HS genes
HS_gene_within = loadHS('%s/summary/all.species.txt'%(snp_dir))
#HS_gene_all = loadHS('%s/summary/all.species.txt.High_select2.txt'%(snp_dir))

alloutput=set()
snp_folder = glob.glob('%s/*.donor.*.raw.vcf.filtered.vcf.final.snp.txt' % (snp_dir))
print(snp_folder)
for snpfile in snp_folder:
    alloutput = load_snp(snpfile,alloutput)

f1 = open('%s/summary/all.species.txt.Truncated.sum'%(snp_dir), 'w')
f1.write(''.join(list(alloutput)))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# genomic diversity among strains carried different PE alleles of the same lineage
import os,glob
allmultisnp = '/scratch/users/anniz44/genomes/donor_species/MG/summary/all.withinHS.snp.genename.depthcheck.count.multigenotype.sum'
changename = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.selected.gene.faa.changename.txt'
genome_dir1 = '/scratch/users/anniz44/genomes/donor_species/selected_species/round1/'
genome_dir2 = '/scratch/users/mit_alm/IBD_Evo/BA/Assembly_for_gene_flow/'
genome_dir3 = '/scratch/users/mit_alm/IBD_Evo/BL/Assembly_for_gene_flow/'
genome_dir4 = '/scratch/users/mit_alm/IBD_Evo/PB/Assembly_for_gene_flow/'
snpfolder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/'
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/fastani'
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/fastani'
try:
    os.mkdir(outputdir)
except IOError:
    pass

try:
    os.mkdir(outputscript)
except IOError:
    pass

def loadsnpfile(snpfile, positionset):
    genomeset = dict()
    genomeorder = dict()
    for lines in open(snpfile,'r'):
        lines_set = lines.split('\n')[0].split('\t')
        if lines.startswith('CHR'):
            i = 9
            for genome in lines_set[9:]:
                genomeorder.setdefault(i,genome)
                i += 1
        else:
            CHRPOS = '%s\t%s'%(lines_set[0],lines_set[1])
            if CHRPOS in positionset:
                major, minor = lines_set[2:4]
                genomeset.setdefault(CHRPOS,[[major,set()],[minor,set()]])
                i = 9
                for alleles in lines_set[9:]:
                    if alleles == major:
                        genomeset[CHRPOS][0][-1].add(genomeorder[i])
                    elif alleles == minor:
                        genomeset[CHRPOS][1][-1].add(genomeorder[i])
                    i += 1
    return genomeset

def findfasta(lineage,genomename):
    if lineage.startswith('1_'):
        fasta = glob.glob('%s/%s/scaffolds.fasta' % (genome_dir2, genomename))+\
        glob.glob('%s/%s/scaffolds.fasta' % (genome_dir3, genomename))+\
        glob.glob('%s/%s/scaffolds.fasta' % (genome_dir4, genomename))
    else:
        fasta = glob.glob('%s/%s/fasta/%s_final.scaffolds.fasta' % (genome_dir1, genomename.split('_g')[0], genomename))
    return fasta[0]

def run_fastani(genomeset,lineage):
    for CHRPOS in genomeset:
        genomeset1,genomeset2=genomeset[CHRPOS]
        genomeset1out = set()
        genomeset2out = set()
        CHRPOSset = CHRPOS.split('\t')
        for genome in genomeset1[-1]:
            genomeset1out.add(findfasta(lineage,genome))
        for genome in genomeset2[-1]:
            genomeset2out.add(findfasta(lineage,genome))
        outputname1 = '%s__%s__%s__%s'%(lineage,CHRPOSset[0],CHRPOSset[1],genomeset1[0])
        outputname2 = '%s__%s__%s__%s' % (lineage, CHRPOSset[0], CHRPOSset[1], genomeset2[0])
        f1 = open('%s/%s.list'%(outputdir,outputname1),'w')
        f1.write('\n'.join(list(genomeset1out))+'\n')
        f1.close()
        f1 = open('%s/%s.list' % (outputdir, outputname2), 'w')
        f1.write('\n'.join(list(genomeset2out))+'\n')
        f1.close()
        command = ('fastANI --rl %s --ql %s -o %s \n' %
                    ('%s/%s.list'%(outputdir,outputname1), '%s/%s.list'%(outputdir,outputname1),
                     '%s/%s.fastaniout' % (outputdir, outputname1)))
        command += ('fastANI --rl %s --ql %s -o %s \n' %
                   ('%s/%s.list' % (outputdir, outputname2), '%s/%s.list' % (outputdir, outputname2),
                    '%s/%s.fastaniout' % (outputdir, outputname2)))
        f1 = open('%s/%s.sh' % (outputscript, outputname1), 'w')
        f1.write('#!/bin/bash\nsource ~/.bashrc\n%s' % (command))
        f1.close()

newgenename = dict()
for lines in open(changename,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    newgenename.setdefault(lines_set[3],lines_set[2])

donorspecies = dict()
for lines in open(allmultisnp,'r'):
    if not lines.startswith('subgroupnew'):
        lines_set = lines.split('\t')
        lineage = lines_set[0].split('__')[0].replace('_CL','_clustercluster')
        donorspecies.setdefault(lineage,
                                set())
        donorspecies[lineage].add('%s\t%s'%('_'.join(newgenename[lines_set[0].split(' ')[0]].split('_')[:-1]),
                                            lines_set[0].split(' ')[1]))
for lineage in donorspecies:
    snpfile = '%s/%s.all.parsi.fasta.sum.txt'%(snpfolder,lineage)
    positionset = donorspecies[lineage]
    genomeset = loadsnpfile(snpfile, positionset)
    run_fastani(genomeset,lineage)

f1 = open(os.path.join(outputscript, '../allfastani.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(outputscript, '*.sh')):
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/../allfastani.sh'%(outputscript))

################################################### END ########################################################
################################################### SET PATH ########################################################
# sum up fastani for genomic diversity among strains carried different PE alleles of the same lineage
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/fastani'

alloutput = []
for files in glob.glob('%s/*.fastaniout'%(outputdir)):
    donor_species_allele = os.path.split(files)[-1].split('.fastaniout')[0]
    for lines in open(files,'r'):
        alloutput.append(
            '%s\t%s'%(donor_species_allele,lines)
        )

f1 = open('%s/allfastaniout.sum'%(outputdir),'w')
f1.write('allele\tgenome1\tgenome2\tani\tgenemap\ttotalgenes\n')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# sum up ref genome gene length and non-ORF length -> cluster_length.py
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
import statistics

assemblyfolder1 = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/'
assemblyfolder2 = '/scratch/users/anniz44/genomes/donor_species/vcf_round1/co-assembly/'
assemblyfolder3 = '/scratch/users/anniz44/genomes/donor_species/WGS_old/WGS/vcf_round1/co-assembly/'
fasta = '.fasta.noHM.fasta'
outputfile = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/clonal_genelength.txt'

alloutput = []
for files in glob.glob('%s/*/*%s'%(assemblyfolder1,fasta))+\
             glob.glob('%s/*/*%s'%(assemblyfolder2,fasta))+\
             glob.glob('%s/*/AlOn*%s'%(assemblyfolder3,fasta))+\
             glob.glob('%s/*/BaFr*%s'%(assemblyfolder3,fasta)):
    donor_species = os.path.split(files)[-1].split(fasta)[0]
    species = donor_species.split('_')[0]
    print(files)
    if species!= 'prokka':
        total_length = 0
        total_genelength = []
        gene_num = 0
        for record in SeqIO.parse(files, 'fasta'):
            total_length += len(str(record.seq))
        filesize = 0
        try:
            filesize = int(os.path.getsize(files + '.fna'))
        except FileNotFoundError:
            pass
        if filesize == 0:
            os.system('prodigal -q -i %s -d %s -a %s' % (files, files + '.fna',files + '.faa'))
        for record in SeqIO.parse(files + '.fna', 'fasta'):
            total_genelength.append(len(str(record.seq)))
            gene_num += 1
        alloutput.append('%s\t%s\t%s\t%.1f\t%s\n'%(species,donor_species.replace('_PB_','_PaDi_'),
                                           total_length - sum(total_genelength),
                                           statistics.mean(total_genelength),gene_num))


f1 = open('%s'%(outputfile),'w')
f1.write('species\tdonor_species\tnonORFlength\tavgORFlength\ttotal_gene_num\n')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# sum up same SNPs across lineages -> sameSNP.py
import os,glob
SNPfiles = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/*.all.parsi.fasta.sum.txt')
HSfileout = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.species.txt.sameSNP.multipledonor.txt'
all_filter_gene_fasta_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/summary/all.denovo.gene.faa'

cutoff = 0.7
cmd_cluster = ('%s -sort length -cluster_fast %s -id %s -centroids %s.%s.cluster.aa -uc %s.%s.uc -threads %s\n'
                   % ('usearch', all_filter_gene_fasta_file, cutoff, all_filter_gene_fasta_file,cutoff,
                      all_filter_gene_fasta_file, cutoff,40))
#os.system(cmd_cluster)
cluster_file = all_filter_gene_fasta_file + '.%s.uc'%(cutoff)
#cluster_file = all_filter_gene_fasta_file + '.uc'
Clusters = dict()
for lines in open(cluster_file, 'r'):
    line_set = lines.split('\n')[0].split('\t')
    cluster = line_set[1]
    record_name = line_set[8]
    Clusters.setdefault(record_name, cluster)

HSsum = dict()
genetocluster = dict()
for aSNPfile in SNPfiles:
    donor_species = os.path.split(aSNPfile)[-1].split('.all.parsi.fasta.sum.txt')[0]
    CL = donor_species.split('.donor.')[0]
    donor = donor_species.split('.donor.')[1]
    donor_species_new = donor_species.replace('_clustercluster', '_CL').replace('PB_', 'PaDi_')
    for lines in open(aSNPfile,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\n')[0].split('\t')
            genename = lines_set[5]
            POS = lines_set[6]
            N_S = lines_set[4]
            if genename == 'None':
                genename = lines_set[0]
                POS = lines_set[1]
                newgenename = '%s\t%s\t%s\tNone\t%s' % (CL, genename, POS,N_S)
            else:
                genename2 = '%s__C_%s_G_%s' % (donor_species_new, genename.split('_')[1], genename.split('_')[-1])
                cluster = Clusters.get(genename2)
                newgenename = '%s\t%s' % (cluster,POS)
                genetocluster.setdefault(newgenename,set())
                genetocluster[newgenename].add('%s\t%s\t%s\t%s\t%s' % (CL, genename, POS,genename2,N_S))
            HSsum.setdefault(newgenename, set())
            HSsum[newgenename].add(donor_species)

alloutput = []
for newgenename in HSsum:
    donorcutoff = 3
    donornum = len(HSsum[newgenename])
    if donornum >= donorcutoff:
        if newgenename in genetocluster:
            for oldgenename in genetocluster[newgenename]:
                alloutput.append('%s\t%s\t%s\t%s\t\n' % (oldgenename,newgenename, donornum,
                                                         ';'.join(list(HSsum[newgenename]))
                                                         ))
        else:
            alloutput.append('%s\tNone\tNone\t%s\t%s\t\n'%(newgenename,donornum,
                                                         ';'.join(list(HSsum[newgenename]))
                                                           ))

f1 = open('%s'%(HSfileout),'w')
f1.write('Refgenome\tgenename\tPOS\tnewgenename\tN_or_S\tcluster\tPOS2\tnum_donor\tdonor_list\t\n')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# test MG filter
import os,glob
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf')
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/testMG.sh'
fastaformat = '.fasta'
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/\n'%(vcffile,fastaformat)

cmds += 'python SNPsum_meta.py  -o /scratch/users/anniz44/genomes/donor_species/\n'

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmds))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# re run MG filtering
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/refilterMG/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/refilterMG.sh'
jub_num = 200
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    i += 1
    if (i%jub_num) == 0:
        f1 = open('%s/refilterMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
        f1.write(''.join(cmds))
        f1.close()
        cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
        cmdsall += 'jobmit %s %s\n'%('%s/refilterMG.%s.sh'%(outputscript,int(i/jub_num)),
                                     'refilterMG.%s.sh' % (int(i / jub_num)))
    tempoutputfolder = '%s/refilterMG_%s'%(outputscript,int(i/jub_num))
    try:
        os.mkdir(tempoutputfolder)
    except IOError:
        pass
    cmds += 'unzip %s -d %s\n'%(vcffile,tempoutputfolder)
    cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n'%(tempoutputfolder,outputdir)
    cmds += 'rm -r %s/scratch/\n'%(tempoutputfolder)
    newvcffile = os.path.join(outputdir,os.path.split(vcffile)[-1].split('.zip')[0])
    if '.fasta' in newvcffile:
        fastaformat = '.fasta'
    else:
        fastaformat = '.1.fq'
    cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/\n'%(newvcffile,fastaformat)
    cmds += 'rm %s\n'%(newvcffile)

f1 = open('%s/refilterMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
f1.write(''.join(cmds))
f1.close()
cmdsall += 'jobmit %s %s\n'%('%s/refilterMG.%s.sh'%(outputscript,int(i/jub_num)),
                                     'refilterMG.%s.sh' % (int(i / jub_num)))

cmdsall += '#python SNPsum_meta.py  -o /scratch/users/anniz44/genomes/donor_species/\n'
f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# MG trunc filtering
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncMG/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncMG.sh'
jub_num = 200
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    if '.fasta' in vcffile:
        fastaformat = '.fasta'
    else:
        fastaformat = '.1.fq'
    samplename = os.path.split(vcffile)[-1].split(fastaformat)[0]
    donorspecies = os.path.split(vcffile)[-1].split(fastaformat + '.')[1].split('.raw.vcf')[0]
    try:
        f1 = open('%s/../truncsum/%s/%s.IN.%s.snp.sum'%(outputdir,donorspecies,samplename,donorspecies),'r')
    except IOError:
        try:
            f1 = open('%s/../truncsumHS/%s/%s.IN.%s.snp.sum' % (outputdir, donorspecies, samplename, donorspecies), 'r')
        except IOError:
            i += 1
            if (i%jub_num) == 0:
                f1 = open('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
                f1.write(''.join(cmds))
                f1.close()
                cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
                cmdsall += 'jobmit %s %s\n'%('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),
                                             'truncMG.%s.sh' % (int(i / jub_num)))
            tempoutputfolder = '%s/truncMG_%s'%(outputscript,int(i/jub_num))
            try:
                os.mkdir(tempoutputfolder)
            except IOError:
                pass
            cmds += 'unzip %s -d %s\n'%(vcffile,tempoutputfolder)
            cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n'%(tempoutputfolder,outputdir)
            cmds += 'rm -r %s/scratch/\n'%(tempoutputfolder)
            newvcffile = os.path.join(outputdir,os.path.split(vcffile)[-1].split('.zip')[0])
            cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta_trunc.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/\n'%(newvcffile,fastaformat)
            cmds += 'rm %s\n'%(newvcffile)

f1 = open('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
f1.write(''.join(cmds))
f1.close()
cmdsall += 'jobmit %s %s\n'%('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),
                                     'truncMG.%s.sh' % (int(i / jub_num)))

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# MG trunc sum
import os,glob
inputdir = glob.glob('/scratch/users/anniz44/genomes/donor_species/MG/truncsumHS/*/')+\
           glob.glob('/scratch/users/anniz44/genomes/donor_species/MG/truncsum/*/')
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/'

alloutput = []
for a_folder in inputdir:
    donorspecies = os.path.split(os.path.split(a_folder)[-2])[-1]
    for a_file in glob.glob('%s/*.snp.sum'%(a_folder)):
        samplename = os.path.split(a_file)[-1].split('.IN.')[0]
        for lines in open(a_file):
            if not lines.startswith('CHR'):
                alloutput.append('%s\t%s\t%s'%(donorspecies,samplename,lines))

f1 = open('%s/all.truncsum.txt'%(outputdir),'w')
f1.write('Lineage\tSample\tCHR\tPOS\tGene\tGenePOS\twithinHS\tMajor_ALT\tMinor_ALT\tMajor_ALT_freq\tMinor_ALT_freq\tRef_aa\tSNP_aa\n')
f1.write(''.join(alloutput))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# MG PE trunc all samples at trunc snp
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncMGHS/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncMGHS.sh'
snp_file = '/scratch/users/anniz44/genomes/donor_species/MG/summary/all.truncsum.filterfreq.trunclist.withHS.txt'
jub_num = 50
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

def load_snp(snpfile):
    SNP_list = dict()
    SNP_list2 = dict()
    if snpfile != 'None':
        for lines in open(snpfile,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            lineage,CHRPOS,Major,Minor,Gene,Genepos = lines_set[:6]
            SNP_list.setdefault(lineage,set())
            SNP_list[lineage].add(CHRPOS)
            SNP_list2.setdefault(CHRPOS, [Minor,Gene,Genepos])
    return[SNP_list,SNP_list2]

SNP_list,SNP_list2 = load_snp(snp_file)

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    if '.fasta' in vcffile:
        fastaformat = '.fasta'
    else:
        fastaformat = '.1.fq'
    samplename = os.path.split(vcffile)[-1].split(fastaformat)[0]
    donorspecies = os.path.split(vcffile)[-1].split(fastaformat + '.')[1].split('.raw.vcf')[0]
    if donorspecies in SNP_list:
        try:
            f1 = open('%s/../truncsumHSnew/%s/%s.IN.%s.snp.sum'%(outputdir,donorspecies,samplename,donorspecies),'r')
        except IOError:
            i += 1
            if (i%jub_num) == 0:
                f1 = open('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
                f1.write(''.join(cmds))
                f1.close()
                cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
                cmdsall += 'jobmit %s %s\n'%('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),
                                             'truncMG.%s.sh' % (int(i / jub_num)))
            tempoutputfolder = '%s/truncMG_%s'%(outputscript,int(i/jub_num))
            try:
                os.mkdir(tempoutputfolder)
            except IOError:
                pass
            cmds += 'unzip %s -d %s\n'%(vcffile,tempoutputfolder)
            cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n'%(tempoutputfolder,outputdir)
            cmds += 'rm -r %s/scratch/\n'%(tempoutputfolder)
            newvcffile = os.path.join(outputdir,os.path.split(vcffile)[-1].split('.zip')[0])
            cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta_trunc.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/ -snplist %s\n'%(newvcffile,fastaformat,snp_file)
            cmds += 'rm %s\n'%(newvcffile)

f1 = open('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
f1.write(''.join(cmds))
f1.close()

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# MG trunc all SNPs before trunc snp
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncallsnp/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncallsnp.sh'
snp_file = '/scratch/users/anniz44/genomes/donor_species/MG/summary/all.truncsum.txt'
jub_num = 100
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

def load_snp(snpfile):
    SNP_list = dict()
    SNP_list2 = dict()
    if snpfile != 'None':
        for lines in open(snpfile,'r'):
            lines_set = lines.split('\n')[0].split('\t')
            lineage, sample, CHR, POS, Gene, Genepos = lines_set[:6]
            SNP_list.setdefault(lineage, set())
            SNP_list[lineage].add(sample)
            lineagesample = '%s\t%s' % (lineage, sample)
            SNP_list2.setdefault(lineagesample, [])
            SNP_list2[lineagesample].append(' '.join([CHR, POS, Gene, Genepos]))
    return[SNP_list,SNP_list2]

SNP_list,SNP_list2 = load_snp(snp_file)

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    if '.fasta' in vcffile:
        fastaformat = '.fasta'
    else:
        fastaformat = '.1.fq'
    samplename = os.path.split(vcffile)[-1].split(fastaformat)[0]
    donorspecies = os.path.split(vcffile)[-1].split(fastaformat + '.')[1].split('.raw.vcf')[0]
    if samplename in SNP_list.get(donorspecies):
        lineagesample = '%s\t%s' % (donorspecies, samplename)
        snplistall = ';'.join(SNP_list2[lineagesample])
        try:
            f1 = open('%s/../truncallsnps/%s/%s.IN.%s.snp.sum'%(outputdir,donorspecies,samplename,donorspecies),'r')
        except IOError:
            i += 1
            if (i%jub_num) == 0:
                f1 = open('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
                f1.write(''.join(cmds))
                f1.close()
                cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
                cmdsall += 'jobmit %s %s\n'%('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),
                                             'truncMG.%s.sh' % (int(i / jub_num)))
            tempoutputfolder = '%s/truncMG_%s'%(outputscript,int(i/jub_num))
            try:
                os.mkdir(tempoutputfolder)
            except IOError:
                pass
            cmds += 'unzip %s -d %s\n'%(vcffile,tempoutputfolder)
            cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n'%(tempoutputfolder,outputdir)
            cmds += 'rm -r %s/scratch/\n'%(tempoutputfolder)
            newvcffile = os.path.join(outputdir,os.path.split(vcffile)[-1].split('.zip')[0])
            cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta_trunc.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/ -trunclist \'%s\'\n'%(
                newvcffile,fastaformat,snplistall)
            cmds += 'rm %s\n'%(newvcffile)

f1 = open('%s/truncMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
f1.write(''.join(cmds))
f1.close()

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# MG trunc allSNPs sum (before the truncation sites)
import os,glob
inputdir = glob.glob('/scratch/users/anniz44/genomes/donor_species/MG/truncallsnps/*/')
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/'
snpfile = '/scratch/users/anniz44/genomes/donor_species/MG/summary/all.truncsum.withdetails.nosingle.txt'

lineagesample = set()
lineagesample2=set()
for lines in open(snpfile,'r'):
    lines_set = lines.split('\n')[0].split('\t')
    lineage_sample = lines_set[1]
    gene = lines_set[7]
    lineagesample2.add('%s %s'%(lineage_sample,gene))
    lineagesample.add(lineage_sample)

alloutput = []
for a_folder in inputdir:
    donorspecies = os.path.split(os.path.split(a_folder)[-2])[-1]
    for a_file in glob.glob('%s/*.snp.sum'%(a_folder)):
        samplename = os.path.split(a_file)[-1].split('.IN.')[0]
        lineage_sample = '%s %s'%(donorspecies,samplename)
        if lineage_sample in lineagesample:
            for lines in open(a_file):
                if not lines.startswith('CHR'):
                    lines_set = lines.split('\n')[0].split('\t')
                    gene = lines_set[2]
                    if '%s %s'%(lineage_sample,gene) in lineagesample2:
                        alloutput.append('%s\t%s\t%s'%(donorspecies,samplename,lines))
        else:
            pass

f1 = open('%s/all.truncsum.allsnps.txt'%(outputdir),'w')
f1.write('Lineage\tSample\tCHR\tPOS\tGene\tGenePOS\twithinHS\tMajor_ALT\tMinor_ALT\tMajor_ALT_freq\tMinor_ALT_freq\tRef_aa\tSNP_aa\n')
f1.write(''.join(alloutput))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# MG trunc HS sum
import os,glob
inputdir = glob.glob('/scratch/users/anniz44/genomes/donor_species/MG/truncsumHSnew/*/')
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/'

alloutput = []
for a_folder in inputdir:
    donorspecies = os.path.split(os.path.split(a_folder)[-2])[-1]
    for a_file in glob.glob('%s/*.snp.sum'%(a_folder)):
        samplename = os.path.split(a_file)[-1].split('.IN.')[0]
        for lines in open(a_file):
            if not lines.startswith('CHR'):
                alloutput.append('%s\t%s\t%s'%(donorspecies,samplename,lines))

f1 = open('%s/all.truncsum.HS.txt'%(outputdir),'w')
f1.write('Lineage\tSample\tCHR\tPOS\tGene\tGenePOS\twithinHS\tMajor_ALT\tMinor_ALT\tMajor_ALT_freq\tMinor_ALT_freq\tRef_aa\tSNP_aa\n')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# MG indel filtering
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/indelMG/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/indelMG.sh'
jub_num = 200
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    if '.fasta' in vcffile:
        fastaformat = '.fasta'
    else:
        fastaformat = '.1.fq'
    samplename = os.path.split(vcffile)[-1].split(fastaformat)[0]
    donorspecies = os.path.split(vcffile)[-1].split(fastaformat + '.')[1].split('.raw.vcf')[0]
    try:
        f1 = open('%s/../indelsum/%s/%s.IN.%s.snp.sum'%(outputdir,donorspecies,samplename,donorspecies),'r')
    except IOError:
        i += 1
        if (i%jub_num) == 0:
            f1 = open('%s/indelMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
            f1.write(''.join(cmds))
            f1.close()
            cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
            cmdsall += 'jobmit %s %s\n'%('%s/indelMG.%s.sh'%(outputscript,int(i/jub_num)),
                                         'indelMG.%s.sh' % (int(i / jub_num)))
        tempoutputfolder = '%s/indelMG_%s'%(outputscript,int(i/jub_num))
        try:
            os.mkdir(tempoutputfolder)
        except IOError:
            pass
        cmds += 'unzip %s -d %s\n'%(vcffile,tempoutputfolder)
        cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n'%(tempoutputfolder,outputdir)
        cmds += 'rm -r %s/scratch/\n'%(tempoutputfolder)
        newvcffile = os.path.join(outputdir,os.path.split(vcffile)[-1].split('.zip')[0])
        cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta_indel.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/\n'%(newvcffile,fastaformat)
        cmds += 'rm %s\n'%(newvcffile)

f1 = open('%s/indelMG.%s.sh'%(outputscript,int(i/jub_num)),'a')
f1.write(''.join(cmds))
f1.close()
cmdsall += 'jobmit %s %s\n'%('%s/indelMG.%s.sh'%(outputscript,int(i/jub_num)),
                                     'indelMG.%s.sh' % (int(i / jub_num)))

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# MG indel HS sum
import os,glob
inputdir = glob.glob('/scratch/users/anniz44/genomes/donor_species/MG/indelsum/*/')
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/'
HS_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.species.txt'

HS_set = dict()

for lines in open(HS_file):
    if 'True' in lines:
        lines_set=lines.split('\t')
        donorspecies,gene = lines_set[0:2]
        donorspecies=donorspecies.split('.donor')[0]
        HS_set.setdefault(donorspecies,set())
        HS_set[donorspecies].add(gene)

alloutput = []
for a_folder in inputdir:
    donorspecies = os.path.split(os.path.split(a_folder)[-2])[-1]
    donorspecies = donorspecies.replace('_PB_','_PaDi_')
    if donorspecies in HS_set:
        HS_set_gene = HS_set[donorspecies]
        for a_file in glob.glob('%s/*.snp.sum'%(a_folder)):
            samplename = os.path.split(a_file)[-1].split('.IN.')[0]
            for lines in open(a_file):
                if not lines.startswith('CHR'):
                    lines_set = lines.split('\t')
                    gene = lines_set[2]
                    if gene in HS_set_gene:
                        alloutput.append('%s\t%s\t%s'%(donorspecies,samplename,lines))
    else:
        print(donorspecies)

f1 = open('%s/all.indelsum.HS.txt'%(outputdir),'w')
f1.write('Lineage\tSample\tCHR\tPOS\tGene\tGenePOS\tMajor_ALT\tMinor_ALT\tMajor_ALT_freq\tMinor_ALT_freq\n')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# MG PE mapping length + all gene mapping results
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncPEsum/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/truncPEsum.sh'
outputsumdir = '/scratch/users/anniz44/genomes/donor_species/MG/summary/'
jub_num = 200
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

trunc_output  = outputsumdir + '/all.lineage.PESNP.sum'
f1 = open(trunc_output,'w')
f1.write('Lineage\tSample\tNo.PESNPs\n')
f1.close()

trunc_output_gene  = outputsumdir + '/all.lineage.genemapping.sum'
f1 = open(trunc_output_gene,'w')
f1.write('Lineage\tSample\tcluster_count_PE\tcluster_count_noPE\tgene_count_PE\tgene_count_noPE\n')
f1.close()

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    if '.fasta' in vcffile:
        fastaformat = '.fasta'
    else:
        fastaformat = '.1.fq'
    samplename = os.path.split(vcffile)[-1].split(fastaformat)[0]
    donorspecies = os.path.split(vcffile)[-1].split(fastaformat + '.')[1].split('.raw.vcf')[0]
    i += 1
    if (i % jub_num) == 0:
        f1 = open('%s/truncPEsum.%s.sh' % (outputscript, int(i / jub_num)), 'a')
        f1.write(''.join(cmds))
        f1.close()
        cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
        cmdsall += 'jobmit %s %s\n' % ('%s/truncPEsum.%s.sh' % (outputscript, int(i / jub_num)),
                                       'truncPEsum.%s.sh' % (int(i / jub_num)))
    tempoutputfolder = '%s/truncPEsum_%s' % (outputscript, int(i / jub_num))
    try:
        os.mkdir(tempoutputfolder)
    except IOError:
        pass
    cmds += 'unzip %s -d %s\n' % (vcffile, tempoutputfolder)
    cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n' % (tempoutputfolder, outputdir)
    cmds += 'rm -r %s/scratch/\n' % (tempoutputfolder)
    newvcffile = os.path.join(outputdir, os.path.split(vcffile)[-1].split('.zip')[0])
    cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta_genesum.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/\n' % (
    newvcffile, fastaformat)
    cmds += 'rm %s\n' % (newvcffile)

f1 = open('%s/truncPEsum.%s.sh'%(outputscript,int(i/jub_num)),'a')
f1.write(''.join(cmds))
f1.close()

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# MG snp median depth
import os,glob
outputfile = '/scratch/users/anniz44/genomes/donor_species/MG/summary/all.lineage.depth.sum'
alldepthfiles = glob.glob('/scratch/users/anniz44/genomes/donor_species/MG/covsum/*/*.raw.vcf.depth.sum')

alloutput = []
for files in alldepthfiles:
    sample = os.path.split(files)[-1].split('.fasta')[0].split('.1.fq')[0]
    lineage = os.path.split(os.path.split(files)[0])[-1]
    totalloci = 0
    mediandepth = 0
    depth_to_loci = dict()
    alldepth = list()
    for lines in open(files,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\n')[0].split('\t')
            depth = int(lines_set[1])
            if depth > 0:
                num_loci = int(lines_set[2])
                alldepth.append(depth)
                depth_to_loci.setdefault(depth,num_loci)
                totalloci += num_loci
    newtotalloci = 0
    for depth in sorted(alldepth):
        newtotalloci += depth_to_loci[depth]
        if newtotalloci >= totalloci/2:
            mediandepth = depth
            break
    alloutput.append('%s\t%s\t%.1f\t%s\n'%(sample,lineage,mediandepth,totalloci))

f1 = open('%s'%(outputfile),'w')
f1.write('sample\tlineage\tmediandepth\ttotalloci\n')
f1.write(''.join(alloutput))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# check all depth
import os,glob
outputdir = '/scratch/users/anniz44/genomes/donor_species/MG/bwa/'
allvcffiles = glob.glob('/scratch/users/anniz44/genomes/donor_species//MG/bwa/finished/*/*.vcf.zip')
outputscript = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/MGdepth/'
outputvcf = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/MGdepth.sh'
subset = ['SRR2938117.fasta.BA_clustercluster6','SRR2912781.fasta.1_BL_IBD_1_clustercluster4',
          'SRR2938117.fasta.1_BL_IBD_1_clustercluster4','SRR2938079.fasta.1_BL_IBD_1_clustercluster4']
os.system('rm -r %s'%(outputscript))
try:
    os.mkdir(outputscript)
except IOError:
    pass

cmdsall = '#!/bin/bash\nsource ~/.bashrc\n'
i = 0
cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
for vcffile in allvcffiles:
    if any(vcfname in vcffile for vcfname in subset):
        i += 1
        if (i%2000) == 0:
            f1 = open('%s/MGdepth.%s.sh'%(outputscript,int(i/2000)),'a')
            f1.write(''.join(cmds))
            f1.close()
            cmds = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
            cmdsall += 'jobmit %s %s\n'%('%s/MGdepth.%s.sh'%(outputscript,int(i/2000)),
                                         'MGdepth.%s.sh' % (int(i / 2000)))
        tempoutputfolder = '%s/MGdepth_%s'%(outputscript,int(i/2000))
        try:
            os.mkdir(tempoutputfolder)
        except IOError:
            pass
        cmds += 'unzip %s -d %s\n'%(vcffile,tempoutputfolder)
        cmds += 'mv %s/scratch/users/anniz44/genomes/donor_species//MG/bwa/*.vcf %s/\n'%(tempoutputfolder,outputdir)
        cmds += 'rm -r %s/scratch/\n'%(tempoutputfolder)
        newvcffile = os.path.join(outputdir,os.path.split(vcffile)[-1].split('.zip')[0])
        if '.fasta' in newvcffile:
            fastaformat = '.fasta'
        else:
            fastaformat = '.1.fq'
        cmds += 'python /scratch/users/anniz44/scripts/1MG/donor_species/assembly/SNPfilter_meta.allloci.py -vcf %s -snp /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/ -mfq %s -o /scratch/users/anniz44/genomes/donor_species/\n'%(newvcffile,fastaformat)
        cmds += 'rm %s\n'%(newvcffile)

f1 = open('%s/MGdepth.%s.sh'%(outputscript,int(i/2000)),'a')
f1.write(''.join(cmds))
f1.close()
cmdsall += 'jobmit %s %s\n'%('%s/MGdepth.%s.sh'%(outputscript,int(i/2000)),
                                     'MGdepth.%s.sh' % (int(i / 2000)))

f1 = open('%s'%(outputvcf),'w')
f1.write(''.join(cmdsall))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# simulation to loci info
import os
alloutput = []
for lines in open('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.genome.gene.fnatruncation.simulation.withmutfreq.txt','r'):
    GC = lines.split('\t')[0]
    truncation_loci = lines.split('\t')[-1].split('\n')[0]
    truncation_loci = truncation_loci.replace('[','').replace(']','')
    for loci in truncation_loci.split(','):
        alloutput.append('%s\t%s\n'%(GC,loci.replace(' ','')))

f1 = open('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.genome.gene.fnatruncation.simulation.withmutfreq.loci.txt','w')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# start codon near truncations not used
from Bio import SeqIO
from Bio.Seq import Seq
import glob
import os

fasta = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.genome.gene.fna'
truncfile = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/summary/all.truncsum.withdetails.addtruncpercentage.uniqueGC.txt'
codon_check = 1 # search 10 codons
gene_trunc = dict()
for lines in open(truncfile,'r'):
    if not lines.startswith('gene_name'):
        lines_set = lines.split('\n')[0].split('\t')
        gene_name = lines_set[0]
        GC,gene_length,trunc_loci = lines_set[-3:]
        trunc_loci = float(trunc_loci)
        gene_trunc.setdefault(gene_name,[GC,int(gene_length),set()])
        gene_trunc[gene_name][-1].add(trunc_loci)

outputfasta = []
outputstartcodon = []
outputstartcodonsum = dict()
for record in SeqIO.parse(fasta, 'fasta'):
    record_id = str(record.id)
    record_seq = str(record.seq)
    if record_id in gene_trunc:
        outputfasta.append('>%s\n%s\n'%(record_id,record_seq))
        GC, gene_length, trunc_lociset = gene_trunc[record_id]
        for trunc_loci in trunc_lociset:
            trunc_loci2 = int(trunc_loci*gene_length/100)
            codon_loci = int(trunc_loci2/3)
            start_loci = int(max(1,codon_loci-codon_check)) # start with first codon
            stop_loci = min(int(gene_length/3),codon_loci + codon_check)
            record_seq_sub = [record_seq[i*3:(i*3+3)] for i in range(start_loci,stop_loci+1)]
            trunc_loci_interset = int(trunc_loci/5)
            outputstartcodonsum.setdefault(trunc_loci_interset,[0,0])
            if 'ATG' in record_seq_sub or 'CAT' in record_seq_sub:
                outputstartcodon.append('%s\t%s\t%s\tTrue\n' % (record_id, GC, trunc_loci))
                outputstartcodonsum[trunc_loci_interset][0] += 1
            else:
                outputstartcodon.append('%s\t%s\t%s\tFalse\n' % (record_id, GC, trunc_loci))
                outputstartcodonsum[trunc_loci_interset][1] += 1

outputsum = []
for trunc_loci_interset in outputstartcodonsum:
    outputsum.append('%s\t%s\t%s\n'%(trunc_loci_interset,
                                 outputstartcodonsum[trunc_loci_interset][0],
                                 outputstartcodonsum[trunc_loci_interset][1]))

f1 = open('%s.trunc.start.codon.sum.txt'%(fasta),'w')
f1.write('truncloci_per\twith_near_start_codon\twithout_near_start_codon\n')
f1.write(''.join(outputsum))
f1.close()

f1 = open('%s.trunc.fna'%(fasta),'w')
f1.write(''.join(outputfasta))
f1.close()

f1 = open('%s.trunc.start.codon.txt'%(fasta),'w')
f1.write('gene_name\tGC\ttruncloci_per\tnear_start_codon\n')
f1.write(''.join(outputstartcodon))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# check wrong vcf_round2
import os,glob
shfileall = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/*.sh')
filesize = 29

cmds = '#!/bin/bash\nsource ~/.bashrc\n'
for shfile in shfileall:
    try:
        filesize = int(os.path.getsize(shfile))
    except FileNotFoundError:
        pass
    if filesize == 29:
        donor_species = os.path.split(shfile)[-1].split('.all.flt.snp.vcf.sh')[0]
        print(donor_species)
        cmds += ('rm /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/%s.*all.parsi.fasta\n'%(
            donor_species
        ))
        ref_file = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/%s.*.ref.raw.vcf'%(donor_species))
        if ref_file!= []:
            cmds += ('bcftools filter -T /scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/%s.all.flt.snp.vcf.CHRPOS %s > %s.filter\n'%(donor_species,ref_file[0],ref_file[0]))
        else:
            print(donor_species,'no ref')
        cmds += ('sh /scratch/users/anniz44/scripts/1MG/donor_species/assembly/vcfremoverec/%s.sh\n'%(donor_species))

f1 = open('checkalltree2.sh','w')
f1.write(cmds)
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# check finished vcf_round2
# vcfround2finish.py
import os,glob
parsiall = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/*.all.parsi.fasta')

for parsi in parsiall:
    donor_species = os.path.split(parsi)[-1].split('_cluster')[0]
    clonal_file = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/%s.*'%(donor_species))
    if clonal_file!=[]:
        os.system('mv %s /scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/finished/'%(' '.join(clonal_file)))
################################################### END ########################################################
################################################### SET PATH ########################################################
#convert tree to tree.pdf -> treetopdf.py
import os,glob
alltree = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/tree/*.tree')
cmds = '#!/bin/bash\nsource ~/.bashrc\n'
cmds += 'cd /scratch/users/anniz44/bin/pro/FigTree_v1.4.4/\n'
for tree in alltree:
    try:
        f1 = open('%s.pdf'%(tree),'r')
    except FileNotFoundError:
        cmds += ('./figtree -graphic PDF %s %s.pdf\n'%(tree,tree))

f1 = open('/scratch/users/anniz44/scripts/1MG/donor_species/assembly/treetopdf.sh','w')
f1.write(cmds)
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
#compare cluster difference -> clustercheck.py
import os,glob
newcluster = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/*.genome.cluster.txt')
oldclusterfolder = '/scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/finished/'

def donor_num(clusterfile):
    donor_cluster = dict()
    for lines in open(clusterfile, 'r'):
        lines_set = lines.split('\t')
        donor = lines_set[1].split('_')[0]
        donor_cluster.setdefault(donor, set())
        donor_cluster[donor].add(lines_set[2])
    return donor_cluster

newredo = []
for cluster in newcluster:
    species = os.path.split(cluster)[-1].split('.')[0]
    oldcluster = glob.glob('%s/%s'%(oldclusterfolder,os.path.split(cluster)[-1]))[0]
    donor_clusterold = donor_num(oldcluster)
    donor_clusternew = donor_num(cluster)
    for donor in donor_clusternew:
        if len(donor_clusternew[donor])<len(donor_clusterold[donor]):
            print('something wrong with %s %s'%(species,donor))
        elif len(donor_clusternew[donor])>len(donor_clusterold[donor]):
            newredo.append('%s\t%s\t\n'%(species,donor))

f1 = open('/scratch/users/anniz44/genomes/donor_species/vcf_round1/clonal_population/needredo.txt','w')
f1.write(''.join(newredo))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
#check whether sim cutoff works
import os,glob
donor_species_cutoff = []
for lines in open('total_SNP_cutoff.txt','r'):
    donor_species_cutoff.append(lines.split('\t')[0])

output = []
for lines in open('all.species.txt','r'):
    if lines.split('\t')[0] in donor_species_cutoff and int(lines.split('\t')[4])<=2:
        output.append(lines.replace('False', 'FalseFalse').replace('True', 'FalseFalse'))
    else:
        output.append(lines)

f1 = open('all.species.cutoff.txt','w')
f1.write(''.join(output))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
#check sub_clusters -> check_sub_clusters.py
import os,glob
allsumfile = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/*.sum.txt')

allout = []
for sumfile in allsumfile:
    total_snp = 0
    snp_type = dict()
    donor_species = os.path.split(sumfile)[-1].split('.all.parsi.fasta.sum.txt ')[0]
    for lines in open(sumfile,'r'):
        if not lines.startswith('CHR'):
            snptype = lines.split('\t')[8]
            snp_type.setdefault(snptype,[0,0,0])
            total_snp += 1
            if lines.split('\t')[4] == 'N':
                snp_type[snptype][0] += 1
            elif lines.split('\t')[4] == 'S':
                snp_type[snptype][1] += 1
            else:
                snp_type[snptype][2] += 1
    allout.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (donor_species, snptype, snp_type[snptype][0],
                                                snp_type[snptype][1], snp_type[snptype][2], total_snp))

f1 = open('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/all.genotype.sum.txt','w')
f1.write('donor_species\tsnptype\tN\tS\tothers\ttotal_SNP\n')
f1.write(''.join(allout))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
#remove hypermutators -> rm_hypermutators.py
# quantile 50% total 24 SNPs of all sum.txt + more than 20% SNPs (tails, 95%)
import os,glob
input_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
sum_file = '%s/all.genotype.sum.txt'%(input_folder)
allsumfile = glob.glob('%s/*.sum.txt'%(input_folder))

total_SNP_cutoff = 24
ratio_cutoff = 0.2
try:
    os.mkdir(input_folder+'/withhypermutator/')
except IOError:
    pass

snp_type = dict()
for lines in open(sum_file,'r'):
    donor_species,snptype,N,S,others,total_SNP =lines.split('\n')[0].split('\t')
    total_genotype = int(N)+int(S)+int(others)
    if int(total_SNP) > total_SNP_cutoff and total_genotype/total_SNP > ratio_cutoff:
        snp_type.setdefault(donor_species,set())
        snp_type[donor_species].add(snptype)

print(snp_type)
for sumfile in allsumfile:
    donor_species = os.path.split(sumfile)[-1].split('.all.parsi.fasta.sum.txt ')[0]
    if donor_species in snp_type:
        remove_snptype = snp_type[donor_species]
        os.system('mv %s %s/withhypermutator/' % (sumfile, input_folder))
        newoutput = []
        for lines in open(sumfile,'r'):
            if not lines.startswith('CHR'):
                snptype = lines.split('\t')[8]
                if snptype not in remove_snptype:
                    newoutput.append(lines)
            else:
                newoutput.append(lines)
        f1 = open(sumfile, 'w')
        f1.write(''.join(newoutput))
        f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
#check false HS -> check_HS_false.py
import os,glob
input_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
sum_file = '%s/summary/all.species.txt'%(input_folder)
Min_SNP_highselect_cutoff = 1/3000

HS_false = dict()
for lines in open(sum_file,'r'):
    if 'NODE_' in lines and 'False' in lines:
        lines_set = lines.split('\t')
        if int(lines_set[4])>=2 and int(lines_set[9]) == 0 and int(lines_set[5])/int(lines_set[3])>=Min_SNP_highselect_cutoff:
            HS_false.setdefault(lines_set[0],set())
            HS_false[lines_set[0]].add(lines_set[1])

alloutput = []
for donor_species in HS_false:
    print(donor_species)
    vcf_file = '%s/moredetails/%s.raw.vcf.filtered.vcf.noremoverec.snpfreq.txt'%(
        input_folder,donor_species.replace('.donor','.all.donor'))
    allgene = HS_false[donor_species]
    oldline = ''
    oldline2 = ''
    nextline = False
    for lines in open(vcf_file,'r'):
        lines_set = lines.split('\t')
        if lines_set[7] in allgene:
            alloutput.append('%s\t%s\t%s' % (donor_species, lines_set[7], oldline2))
            alloutput.append('%s\t%s\t%s' % (donor_species, lines_set[7], oldline))
            alloutput.append('%s\t%s\t%s' % (donor_species, lines_set[7], lines))
            nextline = True
            targetgene = lines_set[7]
        elif nextline == True and targetgene !='':
            alloutput.append('%s\t%s\t%s' % (donor_species, targetgene, lines))
            nextline = False
            targetgene = ''
        oldline2 = oldline
        oldline = lines

f1 = open('%s/summary/all.HS.false.txt'%(input_folder), 'w')
f1.write(''.join(alloutput))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
#check core flexible -> check_core_flexible.py
import os,glob
input_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/'
sum_file = '%s/summary/all.species.txt'%(input_folder)
core_flexible = '%s/summary/all.denovo.gene.faa.uc.species.sum'%(input_folder)

core_flexible_set = dict()
for lines in open(core_flexible,'r'):
    lines_set = lines.split('\t')
    core_flexible_set.setdefault(lines_set[0],lines_set[3].split('\n')[0])

newoutput = []
for lines in open(sum_file,'r'):
    lines = lines.split('\n')[0]
    if lines.startswith('#donor_species'):
        newoutput.append('%s\tgenus_num\n' % (lines))
    else:
        lines_set = lines.split('\t')
        donor_species = lines_set[0]
        donor_species_new = donor_species.replace('_clustercluster', '_CL').replace('_PB_', '_PaDi_')
        record_id = lines_set[1]
        if 'other' not in record_id and record_id!='0' and 'allspecies' not in record_id:
            #print(record_id)
            record_id = 'C_%s_G_%s' % (record_id.split('_')[1], record_id.split('_')[-1])
            newrecord_id = '%s__%s' % (donor_species_new, record_id)
            if newrecord_id in core_flexible_set:
                newoutput.append('%s\t%s\n'%(lines,core_flexible_set[newrecord_id]))
            else:
                print('%s not in cluster'%(newrecord_id))

f1 = open('%s/summary/all.species.genusnum.txt'%(input_folder), 'w')
f1.write(''.join(newoutput))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
#extract sequences of interests -> extract_fasta.py
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

fasta_file='all.denovo.gene.faa'
SNP_file = glob.glob('SNP_profile.*')

genes = dict()
for aSNP_file in SNP_file:
    genename = aSNP_file.split('.txt')[0].split('.')[-1]
    print(genename)
    for lines in open(aSNP_file,'r'):
        lines_set = lines.split('\t')
        genes.setdefault(lines_set[0],genename)

outputgene = set()
alloutput=[]
seqall = set()
alloutput2 = []
for record in SeqIO.parse(fasta_file, 'fasta'):
    record_id = str(record.id)
    if record_id in genes:
        record_seq = str(record.seq)
        if record_seq not in seqall:
            alloutput2.append('>%s_%s\n%s\n' % (genes[record_id], record_id, record_seq.replace('*','')))
            seqall.add(record_seq)
        alloutput.append('>%s_%s\n%s\n'%(genes[record_id],record_id,record_seq))
        outputgene.add(record_id)

alloutput.sort()
f1 = open(fasta_file + '.select.faa', 'w')
f1.write(''.join(alloutput))
f1.close()
alloutput2.sort()
f1 = open(fasta_file + '.unique.select.faa', 'w')
f1.write(''.join(alloutput2))
f1.close()
print([genename for genename in genes if genename not in outputgene])
################################################### END ########################################################
################################################### SET PATH ########################################################
# find SNPs loci on the structure -> locus_aligned.py
import os,glob
from Bio import SeqIO

#SNP_file = 'SNP_profile.cysL.txt'
#fasta_file = 'SNP_profile.cysL.alignment.fa'
SNP_file = 'SNP_profile.degA.to.ccpA.txt'
fasta_file = 'SNP_profile.degA.to.ccpA.alignment.fa'

def find_SNP_locus(seq,locus):
    for i in range(0,len(seq)):
        if seq[i] != '-':
            locus -= 1
        if locus == 0:
            return i+1

def find_ref_locus(Ref,locus):
    if Ref[locus-1]=='-':
        return ['None','None']
    locus_ref = 0
    for i in range(0,locus):
        if Ref[i] != '-':
            locus_ref += 1
    return [locus_ref,Ref[locus-1]]

f1 = open(SNP_file + '.aligned.txt', 'w')

SNP_loci = dict()
for lines in open(SNP_file, 'r'):
    if not lines.startswith('Newgene'):
        lines_set = lines.split('\n')[0].split('\t')
        newSNP = lines_set[2].replace(lines_set[1],lines_set[-1])
        lines = lines.replace(lines_set[2],newSNP)
        newname = '%s_%s'%(lines_set[-2],lines_set[0])
        SNP_loci.setdefault(newname,[])
        SNP_loci[newname].append([int(lines_set[-1]),lines])
    else:
        f1.write('%s\t%s\t%s\t%s'%('AA_POS_ref','AA_ref','AA_POS_aligned',lines))

newoutput = []
for record in SeqIO.parse(fasta_file, 'fasta'):
    record_id = str(record.id)
    if not record_id.startswith('ccpA'):
        Ref = str(record.seq)

for record in SeqIO.parse(fasta_file, 'fasta'):
    record_id = str(record.id)
    if record_id in SNP_loci:
        record_seq = str(record.seq)
        for temp_locus in SNP_loci[record_id]:
            locus,lines = temp_locus
            locus_new = find_SNP_locus(record_seq,locus)
            locus_Ref,Ref_AA = find_ref_locus(Ref,locus_new)
            newoutput.append('%s\t%s\t%s\t%s'%(locus_Ref,Ref_AA,locus_new,lines))

f1.write(''.join(newoutput))
f1.close()

################################################### END ########################################################
################################################### SET PATH ########################################################
# find ccpA BS -> BSfinding.py
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

input_fasta='/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/BL_clustercluster33/BL_clustercluster33.all.spades2.fasta.noHM.fasta'
input_fna='/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/BL_clustercluster33/BL_clustercluster33.all.spades2.fasta.noHM.fasta.fna'

length_halfBS = 7
allBS = []
allBS.append('BS\tcontig\tlocus\tregion\n')
Mapping_loci = dict()
for record in SeqIO.parse(input_fna, 'fasta'):
    record_id = str(record.id)
    description = str(record.description).replace(' ', '').split('#')
    contig = '_'.join(record_id.split('_')[0:-1])
    Mapping_loci.setdefault(contig, [])
    Mapping_loci[contig].append([int(description[1])-1,
                                 int(description[2])-1])

for record in SeqIO.parse(input_fasta, 'fasta'):
    record_id = str(record.id)
    record_seq = str(record.seq)
    for i in range(0,(len(record_seq)-length_halfBS*2)):
        seq1 = record_seq[i : (i+length_halfBS)]
        seq2 = record_seq[(i + length_halfBS):(i+2*length_halfBS)]
        if seq1 == str(Seq(seq2).reverse_complement()) and seq1!='NNNNNNN':
            if not any(i >=locus[0] and i <= locus[1] for locus in Mapping_loci[record_id]):
                allBS.append('%s%s%s%s\t%s\t%s\tnon-ORF\n'%(record_seq[(i-5):(i)].lower(),seq1,seq2,
                                                   record_seq[(i+2*length_halfBS):(i+2*length_halfBS+5)].lower(),
                                                   record_id,i))
            else:
                allBS.append('%s%s%s%s\t%s\t%s\tORF\n' % (record_seq[(i - 5):(i)].lower(), seq1, seq2,
                                                              record_seq[(i + 2 * length_halfBS):(
                                                                          i + 2 * length_halfBS + 5)].lower(),
                                                              record_id, i))

f1 = open(input_fasta + '.BS.txt', 'w')
f1.write(''.join(allBS))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# find ccpA BS by motif -> BSfinding2.py
import os,glob
from Bio import SeqIO

pvalue_cutoff=1E-5
input_fasta='/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/BL_clustercluster33/BL_clustercluster33.all.spades2.fasta.noHM.fasta'
input_fna='/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/BL_clustercluster33/BL_clustercluster33.all.spades2.fasta.fna'
input_faa='/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/BL_clustercluster33/BL_clustercluster33.all.spades2.fasta.faa'
input_motif = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/BL_clustercluster33/BS_candidate.tsv'
Mapping_loci = dict()
for record in SeqIO.parse(input_fna, 'fasta'):
    record_id = str(record.id)
    description = str(record.description).replace(' ', '').split('#')
    contig = '_'.join(record_id.split('_')[0:-1])
    Mapping_loci.setdefault(contig, [])
    Mapping_loci[contig].append([int(description[1])-1,
                                 int(description[2])-1,record_id])

allBS = []
allBS.append('BS\tpvalue\tlocus\tcontig\tstrand\ttargetgane\tlocusgene\n')
target_gene_list=[]
for lines in open(input_motif,'r'):
    if not lines.startswith('#') and not lines.startswith('motif_id') and lines!='\n':
        lines_set = lines.split('\n')[0].split('\t')
        pvalue = lines_set[7]
        if float(pvalue) <= pvalue_cutoff:
            contig,locus1,locus2,strand = lines_set[2:6]
            locus1=int(locus1)
            locus2=int(locus2)
            targetgene = ''
            notongene = True
            locus_target = 0
            if contig in Mapping_loci:
                for locus in Mapping_loci[contig]:
                    locusre1,locusref2,genename = locus
                    if locus2 < locusre1 and targetgene=='':
                        targetgene = genename
                        locus_target = locusre1
                    if locusre1<locus1 and locus1<locusref2:
                        notongene = False
                if notongene and targetgene!='':
                    seq = lines_set[9]
                    allBS.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                    seq, pvalue, locus1, contig, strand, targetgene, locus_target))
                    target_gene_list.append(targetgene)
                    if strand == '-':
                        # the gene before
                        targetgene = '_'.join(targetgene.split('_')[0:-1]) + '_%s'% (int(targetgene.split('_')[-1])-1)
                        allBS.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(seq,pvalue,locus1,contig,strand,targetgene,locus_target))
                        target_gene_list.append(targetgene)


f1 = open(input_fasta + '.BStarget.txt', 'w')
f1.write(''.join(allBS))
f1.close()

aa_output = []
for record in SeqIO.parse(input_faa, 'fasta'):
    record_id = str(record.id)
    if record_id in target_gene_list:
        aa_output.append('>%s\n%s\n'%(record_id,str(record.seq)))

f1 = open(input_fasta + '.BStarget.faa', 'w')
f1.write(''.join(aa_output))
f1.close()
################################################### END ########################################################
################################################### SET PATH ########################################################
# run fimo and predict aa -> runfimo.py
################################################### END ########################################################
################################################### SET PATH ########################################################
# compare BS -> compareBS.py
################################################### END ########################################################
################################################### SET PATH ########################################################
# BS gene annotation -> annotate_eggnog.py
################################################### END ########################################################
# compare BS loss -> compareBS_loss.py
# 2-3 BSs unique in each genome, no matter whether it's a mutant or wild type


