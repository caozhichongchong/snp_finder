import os,glob
from Bio import SeqIO
from Bio.Seq import Seq

output_folder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/BS/'
input_script = '/scratch/users/anniz44/scripts/1MG/donor_species/assembly/'

input_coassembly = '%s/co-assembly'%(output_folder)
input_bs_file = '%s/binding_results_ccpA.txt'%(output_folder)
motif_file = '%s/T17A.meme'%(output_folder)
script_file = '%s/fimo/'%(input_script)
script_file2 = '%s/predictaa/'%(input_script)
length_halfBS = 7

os.system('rm -r %s'%(script_file))
try:
    os.mkdir(script_file)
except IOError:
    pass
os.system('rm -r %s'%(script_file2))
try:
    os.mkdir(script_file2)
except IOError:
    pass

CL = dict()
for lines in open(input_bs_file,'r'):
    if not lines.startswith('AA_POS_ref'):
        lines_set = lines.split('\t')
        species = lines_set[4].split('_')[0]
        donor = lines_set[5]
        SNP = lines_set[3]
        if SNP == 'T17A':
            CL.setdefault(species,set())
            CL[species].add(donor)

def find_BS(input_fasta):
    allBS = []
    allBS.append('BS\tcontig\tlocus\n')
    for record in SeqIO.parse(input_fasta, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        for i in range(0, (len(record_seq) - length_halfBS * 2)):
            seq1 = record_seq[i: (i + length_halfBS)]
            seq2 = record_seq[(i + length_halfBS):(i + 2 * length_halfBS)]
            if seq1 == str(Seq(seq2).reverse_complement()) and seq1 != 'NNNNNNN':
                allBS.append('%s%s\t%s\t%s\n' % (seq1, seq2,
                                                              record_id, i))
    f1 = open('%s/%s.BS.txt'%(output_file, genomename), 'w')
    f1.write(''.join(allBS))
    f1.close()

def exclude_aa_region(fasta,output_file,genomename):
    faa = '%s/%s.faa' % (output_file, genomename)
    newfasta = '%s/%s.fasta' % (output_file, genomename)
    return newfasta
    os.system('cp %s %s/%s.origin.fasta'%(fasta,output_file, genomename))
    try:
        f1 = open(faa,'r')
    except IOError:
        cmd = 'prodigal -q -i %s -a %s/%s.faa' % (
            fasta, output_file, genomename)
        os.system(cmd)
    AA_loci = dict()
    for record in SeqIO.parse(faa, 'fasta'):
        record_id = str(record.id)
        contig = '_'.join(record_id.split('_')[0:-1])
        description = str(record.description).replace(' ', '').split('#')
        AA_loci.setdefault(contig, [])
        AA_loci[contig].append([int(description[1]) - 1,
                                         int(description[2]) - 1])
    fasta_no_aa = []
    for record in SeqIO.parse(fasta, 'fasta'):
        contig = str(record.id)
        record_seq = str(record.seq)
        record_seq = list(record_seq)
        if contig in AA_loci:
            for loci in AA_loci[contig]:
                # excluding gene regions
                record_seq[loci[0]:loci[1]] = 'N'*(loci[1]-loci[0])
        fasta_no_aa.append('>%s\n%s\n'%(contig,''.join(record_seq)))
    f1 = open(newfasta, 'w')
    f1.write(''.join(fasta_no_aa))
    f1.close()
    return newfasta

for species in CL:
    for donor in CL[species]:
        output_file = '%s/%s_%s/' % (output_folder, species, donor)
        try:
            os.mkdir(output_file)
        except IOError:
            pass
        if donor.startswith('D') or donor.startswith('H'):
            input_fasta=glob.glob('/scratch/users/mit_alm/IBD_Evo/%s/Assembly_for_gene_flow/%s_*/scaffolds.fasta'%(species,donor))
            for fasta in input_fasta:
                cmds = '#!/bin/bash\nsource ~/.bashrc\n'
                genomename = os.path.split(os.path.split(fasta)[0])[-1]
                fasta = exclude_aa_region(fasta, output_file, genomename)
                cmds += 'fimo -o %s/%s %s %s\n'%(output_file,genomename,motif_file,fasta)
                cmds += 'mv %s/%s/fimo.tsv %s/%s.fimo.tsv\n' % (output_file, genomename, output_file, genomename)
                cmds += 'rm -r %s/%s\n' % (output_file, genomename)
                f1 = open('%s/%s.sh'%(script_file,genomename), 'w')
                f1.write(cmds)
                f1.close()
                try:
                    f1 = open('%s/%s.BS.txt'%(output_file, genomename),'r')
                except FileNotFoundError:
                    find_BS(fasta)
        else:
            input_fasta = glob.glob(
                '/scratch/users/anniz44/genomes/donor_species/selected_species/round1/%s_%s/fasta/*.origin.fasta' % (donor,species))
            print(input_fasta)
            for fasta in input_fasta:
                cmds = '#!/bin/bash\nsource ~/.bashrc\n'
                genomename = os.path.split(fasta)[-1].split('.origin.fasta')[0]
                fasta = exclude_aa_region(fasta, output_file, genomename)
                cmds += 'fimo -o %s/%s %s %s\n' % (output_file, genomename, motif_file, fasta)
                cmds += 'mv %s/%s/fimo.tsv %s/%s.fimo.tsv\n' % (output_file, genomename, output_file, genomename)
                cmds += 'rm -r %s/%s\n' % (output_file, genomename)
                f1 = open('%s/%s.sh' % (script_file, genomename), 'w')
                f1.write(cmds)
                f1.close()
                try:
                    f1 = open('%s/%s.BS.txt'%(output_file, genomename),'r')
                except FileNotFoundError:
                    find_BS(fasta)

f1 = open(os.path.join(script_file, '../allfimo.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(script_file, '*.sh')):
    f1.write('jobmit %s %s small1\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run %s/../allfimo.sh'%(script_file))

