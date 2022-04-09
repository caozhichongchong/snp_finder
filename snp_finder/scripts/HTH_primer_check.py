import os,glob

genome_dict = {
    'K1':'P1',
'K2':'P1',
'R1':'P2',
'R2':'P2',
'L1':'P3',
'L2':'P3'
}
allgenomes = glob.glob('/scratch/users/anniz44/genomes/donor_species/vcf_round2/HTH/*.fa')
for genomes in allgenomes:
    genomesname = os.path.split(genomes)[-1].replace('.fa','')
    if genomesname in genome_dict:
        os.system('%s -makeudb_usearch %s -output %s.udb' %
                  ('usearch', genomes, genomes))
        mapping_genome = genomes.replace(genomesname,genome_dict.get(genomesname,genomesname)).replace('.fa','.genome.fasta')
        print(mapping_genome)
        os.system('%s -ublast %s -db %s.udb  -evalue 1 -accel 0.5 -blast6out %s -threads 2 -strand both -maxaccepts 0 -maxrejects 0' %
                  ('usearch', mapping_genome, genomes, genomes + '.genome.txt'))
