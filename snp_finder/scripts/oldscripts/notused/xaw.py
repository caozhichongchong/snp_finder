# start
# clean up
import os
os.system('rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/MG/bwa')
os.system('rm -rf /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/MG/merge/*.raw.vcf')
os.system('mv /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/MG/merge/*.all.flt.snp.vcf /scratch/users/anniz44/genomes/donor_species/selected_species/vcf_round4/MG/merge/rawdata/')

################################################### END ########################################################
