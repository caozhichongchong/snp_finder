# sum up aligned region ORF length and non-ORF length -> cluster_length_new.py
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
assemblyfolder = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/co-assembly/'
fasta = '.all.spades2.fasta.noHM.fasta.faa'
allcovfiles = glob.glob(
'/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/MV_cov/*.cov.MW.txt')
outputfile = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/removerec_SNP/clonal_genelength_new.txt'
max_fraction_ambigious_samples = .4 #If more than % of the samples have ambiguous NT, discard the candidate location
min_fraction_of_good_coverage = 1-max_fraction_ambigious_samples
Cov_dis_overall = 1000 # calculate coverage per 1000 bp
min_depth = 6 #Remove candidate locations have lower than this depth
old_genelength_file = '/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/removerec_SNP/clonal_genelength.txt'

def load_oldgenelengthfile(old_genelength_file):
    lost_assembly_data = dict()
    for lines in open(old_genelength_file, 'r'):
        if '.all.spades1' in lines:
            lines_set = lines.split('\n')[0].split('\t')
            species = lines_set[0].split('.all.spades1')[0]
            nonORFlength = int(float(lines_set[3])*int(lines_set[4]))
            num_genes = int(lines_set[4])
            lost_assembly_data.setdefault(species,[nonORFlength,num_genes])
    return lost_assembly_data

def load_covfile(covfile,fnafile):
    allgenes = dict()
    for record in SeqIO.parse(fnafile, 'fasta'):
        record_id = str(record.id)
        chr = '_'.join(record_id.split('_')[:-1])
        allgenes.setdefault(chr,{})
        record_descriptionset = str(record.description).split(' # ')
        allgenes[chr].setdefault(int(record_descriptionset[1]),int(record_descriptionset[2])-int(record_descriptionset[1])+1) # start, gene length
    covresult = [0,0,0] # genomic POS covered, ORF-region POS, number of genes
    for lines in open(covfile,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\n')[0].split('\t')
            if len([int(x)>=min_depth for x in lines_set[3:]]) >= min_fraction_of_good_coverage*(len(lines_set)-3):
                # >= min_fraction_of_good_coverage isolates with >= min_depth depth
                # good coverage
                CHR, POS = lines_set[:2]
                POS = int(POS)
                if POS > 1:
                    covresult[0] += Cov_dis_overall # genomic POS covered
                    # check genes
                    allgeneschr = allgenes.get(CHR,{})
                    checked_genes = []
                    for genestart in allgeneschr:
                        if genestart<=POS and genestart >= POS - Cov_dis_overall:
                            # genes within this moving window region
                            covresult[1] += allgeneschr[genestart] # ORF-region POS covered
                            covresult[2] += 1 # number of genes covered
                            checked_genes.append(genestart)
                        elif genestart > POS:
                            break
                    for genestart in checked_genes:
                        # remove genes that already checked
                        allgeneschr.pop(genestart)
                    allgenes[chr] = allgeneschr
    return [str(x) for x in covresult]

def load_covfile_nofna(covfile,nonORFlength, num_genes):
    covresult = [0,0,0] # genomic POS covered, ORF-region POS, number of genes
    for lines in open(covfile,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\n')[0].split('\t')
            if len([int(x)>=min_depth for x in lines_set[3:]]) >= min_fraction_of_good_coverage*(len(lines_set)-3):
                # >= min_fraction_of_good_coverage isolates with >= min_depth depth
                # good coverage
                CHR, POS = lines_set[:2]
                POS = int(POS)
                if POS > 1:
                    covresult[0] += Cov_dis_overall # genomic POS covered
    covresult[1] = covresult[0] # assuming no nonORF region
    covresult[2] = int(covresult[0]/nonORFlength*num_genes) # assuming gene distribution is even
    return [str(x) for x in covresult]

alloutput = ['cluster\tgenome_coverage\tORF_region_coverage\tnum_genes_covered\n']
lost_assembly_data = load_oldgenelengthfile(old_genelength_file)
for covfile in allcovfiles:
    covfilename = os.path.split(covfile)[-1]
    assembly_name = covfilename.split('.all.donor')[0]
    fnafile = glob.glob('%s/%s/%s%s'%(assemblyfolder,
                                      assembly_name,assembly_name,
                                   fasta))
    lineagename = covfilename.split('.raw.vcf.filtered.cov.MW.txt')[0].replace('.all.donor', '.donor')
    if len(fnafile) > 0:
        print('process %s'%(covfile))
        covresult = load_covfile(covfile,fnafile[0])
        alloutput.append('%s\t%s\n'%(lineagename,
                         '\t'.join(covresult)))
        print(covresult)
    else:

        nonORFlength, num_genes = lost_assembly_data.get(covfilename.split('_')[0],[0,0])
        print('process %s with no fasta but nonORFlength %s num_genes %s' % (covfile,nonORFlength, num_genes ))
        covresult = load_covfile_nofna(covfile,nonORFlength, num_genes)
        alloutput.append(
            '%s\t%s\n' % (lineagename,
                          '\t'.join(covresult)))
        print(covresult)
    if lineagename in ['AlOn_clustercluster2.donor.am','BaFr_clustercluster3.donor.am']:
        alloutput.append(
            '%s\t%s\n' % (lineagename.replace('AlOn_clustercluster2.donor.am','AlOn_clustercluster1.donor.am').replace(
                'BaFr_clustercluster3.donor.am','BaFr_clustercluster7.donor.am'),
                          '\t'.join(covresult)))

f1 = open(outputfile,'w')
f1.write(''.join(alloutput))
f1.close()

