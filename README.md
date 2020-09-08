# snp_finder
## Introduction
* snp_finder identifies and summarizes SNPs in clonal populations
* snp_finder can easily handel hundreds of clonal populations
* input: WGS (whole genome sequences)
* requirement: python >= 3.0, samtools, bcftools, spades, minimap2, usearch
## Install
`pip install snp_finder`
### latest version (unstable though)
`git clone https://github.com/caozhichongchong/snp_finder.git `\
`cd snp_finder`\
`python setup.py build`\
`python setup.py install`
## A little preparation
### organize files
organize folders of WGS (whole genome sequences) into each species isolated from each subject\
rename folders as 'subject_species'
## Start finding SNPs!
### prepare and correct assembled genomes for each WGS
`snp_finder prepare -i input_folder -fq _1.fastq`\
`sh scripts/snp_finder.sh`
* Output SNVs curated in assemblies: scripts/SNP_currate1.assembly.sum.curated
### cluster clonal populations, identify and summarize SNPs
`snp_finder snp -i input_folder -fq _1.fastq`\
`sh scripts/snp_finder.sh`
* Output the clustering and SNPs of each WGS: scripts/SNP_round*.sum (Round 1-4)
* Output clustering cutoff: scripts/SNP_round*.sum.cutoff
* Output pair-wise SNPs between WGSs: scripts/SNP_round*.allpair.sum
* Output pair-wise SNPs between WGSs: scripts/SNP_round*.allpair.sum
* Output trees across clonal populations: snp_output/vcf_round*/merge/tree/*.tree
* Output final trees within clonal populations: snp_output/vcf_round4/merge_genome/tree/*.tree
* Output (Matlab input) final SNP summary of each clonal population: snp_output/vcf_round4/merge_genome/clonal_pop.all.flt.snp.vcf.filtered.vcf.final.vcf.removerec.*
### calculate dN/dS ratio, select genes under parallel evolution and extract gene sequences
`snp_finder select -i input_folder -fq _1.fastq`\
`sh scripts/snp_finder.sh`
* Output dnds of all clonal populations and genes: snp_output/vcf_round4/merge_genome/summary/all.donor.species.dnds.txt.High_select2.txt
* Extract all genes with de novo mutations: snp_output/vcf_round4/merge_genome/summary/all.denovo.gene.faa
* Extract all genes with parallel evolution: snp_output/vcf_round4/merge_genome/summary/all.selected.gene.faa.High_select2.faa
### track strains in metagenomes
`snp_finder meta -i input_folder -fq _1.fastq -m input_meta -mfq _1.fasta`\
`sh scripts/snp_finder.sh`
* Output coverage/prevalence and of each strain across metagenomes: snp_output/vcf_round4/MG/merge/allcov.sum.txt
### c3ddb users only: annotate genes under parallel evolution and genes with de novo mutaitons
`snp_finder anno -i input_folder -fq _1.fastq`\
`sh scripts/snp_finder.sh`
* Annotate all genes with parallel evolution: snp_output/vcf_round4/merge_genome/summary/all.selected.gene.faa.High_select2.faa.cluster.sum
* Annotate all genes with de novo mutations: snp_output/vcf_round4/merge_genome/summary/all.denovo.gene.faa.cluster.sum
* Annotate all genes with parallel evolution only by kegg: snp_output/vcf_round4/merge_genome/summary/all.selected.gene.faa.High_select2.faa.cluster.kegg.sum.txt
* Annotate all genes with de novo mutations only by kegg: snp_output/vcf_round4/merge_genome/summary/all.denovo.gene.faa.cluster.kegg.sum.txt
