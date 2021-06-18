################################################### END ########################################################
################################################### SET PATH ########################################################
import glob
import os
input_script = '/scratch/users/anniz44/scripts/1MG/vg'
vg_WGS = '/scratch/users/amyxiao/projects/lacto_clock/data/demux/cut_trim/'
vg_WGS2 = '/scratch/users/amyxiao/projects/lacto_clock/data/demux/'
vg_genome = '/scratch/users/amyxiao/projects/lacto_clock/data/demux/spades_results/'
vg_checkm = glob.glob('/scratch/users/amyxiao/projects/lacto_clock/data/demux/spades_results/*/results_checkM_remove1000.txt')
outout_dir = '/scratch/users/anniz44/genomes/vg/WGS'

try:
    os.mkdir('%s/gtdb'%(outout_dir))
except IOError:
    pass

try:
    os.mkdir('%s/fasta'%(outout_dir))
except IOError:
    pass

try:
    os.mkdir('%s/fastq'%(outout_dir))
except IOError:
    pass

def readcheckm(checkm):
    for lines in open(checkm):
        if lines.startswith('genome_'):
            try:
                return lines.split('\t')[1].split('__')[1][0:4]
            except IndexError:
                return lines.split('\t')[1][0:4]

def rungrdb():
    cmds = '#!/bin/bash\nsource ~/.bashrc\nexport GTDBTK_DATA_PATH=/scratch/users/needham/databases/release89/\npy37\n'
    cmds += 'python -m pip install gtdbtk\n'
    cmds += 'gtdbtk classify_wf --genome_dir %s/fasta -x fasta --prefix vg_taxon  --cpus 40 --out_dir %s/gtdb'%(outout_dir,outout_dir)
    f1 = open('%s/gtdb.sh'%(input_script),'w')
    f1.write(cmds)
    f1.close()

Taxon = dict()
cmds = '#!/bin/bash\nsource ~/.bashrc\n'
for checkm in vg_checkm:
    tagnum = os.path.split(os.path.split(checkm)[0])[-1]
    taxon = readcheckm(checkm)
    Taxon.setdefault(taxon,0)
    newname = '%s__g%.4d' % (taxon, Taxon[taxon])
    Taxon[taxon] += 1
    print(taxon,newname)
    genome = glob.glob('%s/%s/genome_final.scaffolds.fa'%(vg_genome,tagnum))
    if genome == []:
        genome = glob.glob('%s/%s/scaffolds.fasta' % (vg_genome, tagnum))
    fastq1 = glob.glob('%s/%s_1.trim.fastq'%(vg_WGS,tagnum))
    fastq2 = glob.glob('%s/%s_2.trim.fastq' % (vg_WGS, tagnum))
    if fastq1 == [] or fastq2 == []:
        fastq1 = glob.glob('%s/out_%s/1_1.fastq' % (vg_WGS2, tagnum))
        fastq2 = glob.glob('%s/out_%s/1_2.fastq' % (vg_WGS2, tagnum))
    print(genome,fastq1,fastq2)
    #os.system('ln -s %s %s/fasta/%s.fasta'%(genome[0],outout_dir,newname))
    #os.system('ln -s %s %s/fastq/%s_1.fastq' % (fastq1[0], outout_dir, newname))
    #os.system('ln -s %s %s/fastq/%s_2.fastq' % (fastq2[0], outout_dir, newname))
    cmds += ('cp %s %s/fasta/%s.fasta\n'%(genome[0],outout_dir,newname))
    cmds += ('cp %s %s/fastq/%s_1.fastq\n' % (fastq1[0], outout_dir, newname))
    cmds += ('cp %s %s/fastq/%s_2.fastq\n' % (fastq2[0], outout_dir, newname))

f1 = open('%s/coppyfile.sh'%(input_script),'w')
f1.write(cmds)
f1.close()
rungrdb()