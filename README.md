# snp_finder
## Introduction
* The snpfinder is a stringent and convenient bioinformatics pipeline for population genetics analysis
* The snpfinder can cluster lineages, call SNPs, filter SNPs, and screen for PE genes (genes under parallel evolution). 
* The snpfinder can easily handle tens of species, hundreds of lineages and process them in parallel. 
* simple input: WGS samples (whole genome sequences)
* environment: python >= 3.0, 
* required tools: [samtools](https://github.com/samtools/samtools/releases/tag/1.11), [bcftools](https://github.com/samtools/bcftools/releases/tag/1.11), [spades](http://cab.spbu.ru/files/release3.13.0/manual.html#sec2.1), [minimap2](https://github.com/lh3/minimap2), [usearch](https://www.drive5.com/usearch/download.html), [bowtie2](https://anaconda.org/bioconda/bowtie2), [prodigal](https://github.com/hyattpd/Prodigal), [roary](https://sanger-pathogens.github.io/Roary/), [blastn](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), [snp-sites](https://github.com/sanger-pathogens/snp-sites), [prokka](https://anaconda.org/bioconda/prokka)\
samtools: `curl -O https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2`\
bcftools `curl -O https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2`\
spades: `wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz`, `tar -xzf SPAdes-3.13.0-Linux.tar.gz` and `cd SPAdes-3.13.0-Linux/bin/`\
minimap2: `git clone https://github.com/lh3/minimap2.git` and `cd minimap2 && make`\
usearch: `curl -O https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz`\
bowtie2: `conda install -c bioconda bowtie2`\
prodigal: `git clone https://github.com/hyattpd/Prodigal.git` and `cd Prodigal && make install`\
roary: `sudo apt-get install roary` or `curl -O https://github.com/sanger-pathogens/Roary/tarball/master`\
blastn: `curl -O https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz`
snp-sites: `apt-get install snp-sites` or `conda install snp-sites`\
prokka: `conda install -c bioconda prokka`
## Install snp_finder
`pip install snp_finder`
## Preparation
### organize files
* Please move WGS samples (whole genome sequences) isolated from one a species in one subject into one folder
* Please rename the folder as 'subject_species' (avoid using '_' in donor names and species names)
* Please rename WGS samples as 'donor_species_gXXXX_1.fq' and 'donor_species_gXXXX_2.fq' (XXXX can be any number and/or common letters)
* Please put all folders ('subject_species') into one big input folder (`-i`)
## Start calling SNPs!
### super lazy auto version
1. Simply run `snp_finder auto -i input_folder`\
or `snp_finder auto -i input_folder -fq _1.fastq`, if the file extension of WGS #1 files is not _1.fq\
To identify **truncated genes** caused by SNPs, run `snp_finder auto -i input_folder -trunc True`\
2. Then please run `sh snpfinder_scripts/snp_finder.sh`
### fast manual version
If you would like to run snp_finder in the most efficient and fast mode, try `snp_finder manual`\
1. Simply run `snp_finder manual -i input_folder`\
or `snp_finder manual -i input_folder -fq _1.fastq`, if the file extension of WGS #1 files is not _1.fq\
To identify **truncated genes** caused by SNPs, run `snp_finder manual -i input_folder -trunc True`\
2. Then please run the scripts in `snpfinder_scripts/snp_finder.sh` step by step, and please make sure all tasks of one step are finished before running the next step!
### tune cutoff for PE identification
As a default, snp_finder identifies PE genes within a lineage by the criteria of `>=2 SNPs, >= 1 SNP per 2 kb, and > 2 unique genotypes` in a lineage.\
snp_finder identifies PE genes across lineages by the criteria of `>=2 SNPs and > 2 unique genotypes` in all lineages of a species.\
To set different cutoff for how many SNPs on a gene to call PE, you could provide files containing cutoff for specific lineages (`-cutoff`) and/or species (`-cutoffsp`).\
Please refer to the format of *example/SNP_cutoff_species.txt* and *example/SNP_cutoff_lineage.txt*.
### reference genomes
As a default, snp_finder uses co-assembly of a species as the reference genome.\
You could also use your preferred reference genomes by providing a file containing the path of genomes to specific folders (`-ref`).\
Please refer to the format of *example/reference.genome.txt*
### pre-assembly
As a default, snp_finder automatically assembles genomes from WGS samples.\
To use pre-assembled genomes, you could put them in the same folder of WGS samples and rename them into 'donor_species_gXXXX.fasta'.
## Output!
### PE genes and sequences
* lineages are named as *specis_cluster_number.donor.donor_name*
* all PE within a lineage and across lineages: *snpfinder_output/vcf_round1/merge/summary/all.species.txt.High_select2.txt*
* PE within a lineage: *snpfinder_output/vcf_round1/merge/summary/all.species.txt*
* Sequences of PE genes: *snpfinder_output/vcf_round1/merge/summary/all.selected.gene.faa*
* Sequences of genes with SNPs: *snpfinder_output/vcf_round1/merge/summary/all.denovo.gene.faa*
* Sequences of truncated genes caused by SNPs (`-trunc True`): *snpfinder_output/vcf_round1/merge/summary/all.trunc.gene.faa*
* Clustered sequences: *snpfinder_output/vcf_round1/merge/summary/\*.faa.cluster.aa*
* dmrca of each lineage: *snpfinder_output/vcf_round1/merge/summary/alldmrca.txt*
### SNPs and trees
* '\*.donor.\*' is the name of a lineage
* all high-quality SNPs in a lineage: *snpfinder_output/vcf_round1/merge/\*.donor.\*.all.parsi.fasta.sum*
* SNP tree of a lineage: *snpfinder_output/vcf_round1/merge/tree/\*.donor.\*.all.parsi.fasta.out.tree*
* SNP tree report of a lineage: *snpfinder_output/vcf_round1/merge/tree/\*.donor.\*.all.parsi.fasta.out.txt*
* vcf of high-quality SNPs: *snpfinder_output/vcf_round1/merge/\*.donor.\*.raw.vcf.filtered.vcf.final.removerec*
### cluster of lineages
* clustering results: *snpfinder_output/pangenome/clonal_population/\*.genome.cluster.txt*
* core-genome alignment: *snpfinder_output/pangenome/clonal_population/\*.tree*
* all pan-genomes: *snpfinder_output/pangenome/pangenome/roary_species*
### co-assembly
* all co-assemblies after removing redundant regions: *snpfinder_output/vcf_round1/co-assembly/\*/*.all.spades1.fasta.noHM.fasta*
## Citations
## Contact
Dr. An-Ni Zhang, anniz44@mit.edu, Alm Lab, Biological Engineering, MIT