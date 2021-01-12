library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(stringr)
library(pheatmap)
###---- clonal population finished----
# test core gene cutoff
SNPfiles = list.files(path = 'final_SNP/clonal_pop/BN10/',pattern="*.SNP.pair.sum.txt")
length(SNPfiles)
i = 11
library(stringr)
library(ggplot2)
library(RColorBrewer)
width_set = 0.1
{
  SNP = read.delim(paste('final_SNP/clonal_pop/BN10/',SNPfiles[i],sep=''),header=T)
  SNP$SNP_ratio = SNP$SNPs_curated_bytree/SNP$Core_gene_length*1000
  SNP$donor1 = str_split_fixed(SNP$Genome1, "_", 2)[,1]
  SNP$donor2 = str_split_fixed(SNP$Genome2, "_", 2)[,1]
  SNP$genome_pair = 'intra-subject'
  SNP$genome_pair[which(as.matrix(SNP$donor1)!=as.matrix(SNP$donor2))]='inter-subject'
  color_set = c(1)
  if (length(unique(SNP$genome_pair))>1)
    color_set = c(5,1)
  species = str_split_fixed(SNPfiles[i], "\\.", 2)[,1]
  ggplot(SNP,aes(
    x=SNP_ratio,
    group = genome_pair
  ))+
    geom_histogram(aes(
      # y=..density..,
      color=genome_pair,
      fill=genome_pair
    ),
    binwidth = width_set,
    position="dodge",
    alpha = 1,
    size=0
    )+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1)
      # axis.text.x=element_blank(),
      # axis.ticks.x=element_blank()#,
      # axis.ticks.y=element_blank(),
      # axis.text.y=element_blank()
    ) +
    scale_fill_manual(values =brewer.pal(5, "PRGn")[color_set])+
    scale_color_manual(values =brewer.pal(5, "PRGn")[color_set])+
    scale_x_continuous(breaks = seq(0,max(SNP$SNP_ratio)+width_set,width_set))+
    xlab(paste('SNPs per 1 kbp (core genes) in ',species,sep=''))+
    ylab('total pairs')
}
ggsave(paste(species,'.pdf',sep = ''),path = 'final_SNP/clonal_pop/BN10/',
       width = 15, height = 12, units = "in")
SNPfiles[i]
i = i+1

# test core gene cutoff
# tree distance
SNPfiles = list.files(path = 'final_SNP/clonal_pop/',pattern="*.SNP.pair.sum.txt")
length(SNPfiles)
i = 1
width_set = 0.01
library(stringr)
library(ggplot2)
library(RColorBrewer)
{
  SNP = read.delim(paste('final_SNP/clonal_pop/',SNPfiles[i],sep=''),header=T)
  SNP$SNP_ratio = SNP$tree_distance
  SNP$donor1 = str_split_fixed(SNP$Genome1, "_", 2)[,1]
  SNP$donor2 = str_split_fixed(SNP$Genome2, "_", 2)[,1]
  SNP$genome_pair = 'intra-subject'
  SNP$genome_pair[which(as.matrix(SNP$donor1)!=as.matrix(SNP$donor2))]='inter-subject'
  color_set = c(1)
  if (length(unique(SNP$genome_pair))>1)
    color_set = c(5,1)
  species = str_split_fixed(SNPfiles[i], "\\.", 2)[,1]
  ggplot(SNP,aes(
    x=SNP_ratio,
    group = genome_pair
  ))+
    geom_histogram(aes(
      # y=..density..,
      color=genome_pair,
      fill=genome_pair
    ),
    binwidth = width_set,
    position="dodge",
    alpha = 1,
    size=0
    )+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1)
      # axis.text.x=element_blank(),
      # axis.ticks.x=element_blank()#,
      # axis.ticks.y=element_blank(),
      # axis.text.y=element_blank()
    ) +
    scale_fill_manual(values =brewer.pal(5, "PRGn")[color_set])+
    scale_color_manual(values =brewer.pal(5, "PRGn")[color_set])+
    scale_x_continuous(breaks = seq(0,max(SNP$SNP_ratio)+width_set,width_set))+
    xlab(paste('Tree distance (core genes) ',species,sep=''))+
    ylab('total pairs')
}
ggsave(paste(species,'.pdf',sep = ''),path = 'final_SNP/clonal_pop/',
       width = 15, height = 12, units = "in")
SNPfiles[i]
i = i+1

#----recombination length cutoff finished ----
NS_ratio = read.delim('final_SNP/vcf/all.donor.species.NSratio.txt',header=T)
NS_ratio=NS_ratio[,-ncol(NS_ratio)]
NS_ratio=NS_ratio[,-ncol(NS_ratio)]
NS_ratio$S_cu =NS_ratio$S
NS_ratio$S_cu[which(NS_ratio$S_cu==0)]=0.5
NS_ratio$NS_ratio = NS_ratio$N/NS_ratio$S_cu
NS_ratio$remove = 'before_remove'
NS_ratio2 = NS_ratio
NS_ratio = read.delim('final_SNP/vcf/all.donor.species.NSratio.removerec.txt',header=T)
NS_ratio=NS_ratio[,-ncol(NS_ratio)]
NS_ratio=NS_ratio[,-ncol(NS_ratio)]
NS_ratio$S_cu =NS_ratio$S
NS_ratio$S_cu[which(NS_ratio$S_cu==0)]=0.5
NS_ratio$NS_ratio = NS_ratio$N/NS_ratio$S_cu
NS_ratio$remove = 'after_remove'
NS_ratio=rbind(NS_ratio,NS_ratio2)
NS_ratio_all = NS_ratio
library(stringr)
NS_ratio_all$species = str_split_fixed(NS_ratio_all$donor_species, "_", 4)[,1]

NS_ratio_all$species[which(NS_ratio_all$species=='1')]=str_split_fixed(NS_ratio_all$donor_species[which(NS_ratio_all$species=='1')], "_", 4)[,2]
allspecies = unique(NS_ratio_all$species)
NS_ratio_all2 = NS_ratio_all
subset_species = allspecies[sample(1:length(allspecies),as.integer(length(allspecies)/2))]

# all species
library(ggplot2)
library(RColorBrewer)
# group 1
NS_ratio_all = NS_ratio_all2[which(NS_ratio_all2$species
                                   %in% subset_species),]
NS_ratio = NS_ratio_all[!duplicated(paste(NS_ratio_all$windowsize,
                                          NS_ratio_all$remove)),c(2,3,4,7)]

for (i in 1: nrow(NS_ratio)){
  temp_window = NS_ratio_all[which(NS_ratio_all$windowsize==
                                     NS_ratio$windowsize[i]&
                                     NS_ratio_all$remove==
                                     NS_ratio$remove[i] ),]
  NS_ratio$N[i] = sum(temp_window$N)
  NS_ratio$S[i] = sum(temp_window$S)
}
NS_ratio$S_cu =NS_ratio$S
NS_ratio$S_cu[which(NS_ratio$S_cu==0)]=0.5
NS_ratio$NS_ratio = NS_ratio$N/NS_ratio$S_cu
Cutoff = NS_ratio$windowsize[which(NS_ratio$NS_ratio==
                                     max(NS_ratio$NS_ratio[which(NS_ratio$windowsize<=1e+5&
                                                                   NS_ratio$windowsize>=1e+3&
                                                                   NS_ratio$remove=='after_remove')]))]
ggplot(NS_ratio,aes(
  x= windowsize,
  y = NS_ratio,
  group = remove))+
  geom_point(aes(
    color=remove,
    fill=remove
  ),
  alpha = 1,
  size=2
  )+
  geom_line(aes(
    color=remove
  ),
  alpha = 0.5,
  size=1
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_x_continuous(breaks = c(1,10,100,1000,10000,100000,1000000),trans = 'log10')+
    geom_vline(xintercept = Cutoff,alpha=0.3)

  # group 2
NS_ratio_all = NS_ratio_all2[which(!(NS_ratio_all2$species
                                   %in% subset_species)),]
NS_ratio = NS_ratio_all[!duplicated(paste(NS_ratio_all$windowsize,
                                          NS_ratio_all$remove)),c(2,3,4,7)]
for (i in 1:nrow(NS_ratio)){
  temp_window = NS_ratio_all[which(NS_ratio_all$windowsize==
                                     NS_ratio$windowsize[i]&
                                     NS_ratio_all$remove==
                                     NS_ratio$remove[i] ),]
  NS_ratio$N[i] = sum(temp_window$N)
  NS_ratio$S[i] = sum(temp_window$S)
}
NS_ratio$S_cu =NS_ratio$S
NS_ratio$S_cu[which(NS_ratio$S_cu==0)]=0.5
NS_ratio$NS_ratio = NS_ratio$N/NS_ratio$S_cu
Cutoff2 = NS_ratio$windowsize[which(NS_ratio$NS_ratio==
                                     max(NS_ratio$NS_ratio[which(NS_ratio$windowsize<=1e+5&
                                                                   NS_ratio$windowsize>=1e+3&
                                                                   NS_ratio$remove=='after_remove')]))]

ggplot(NS_ratio,aes(
  x= windowsize,
  y = NS_ratio,
  group = remove))+
  geom_point(aes(
    color=remove,
    fill=remove
  ),
  alpha = 1,
  size=2
  )+
  geom_line(aes(
    color=remove
  ),
  alpha = 0.5,
  size=1
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_x_continuous(breaks = c(1,10,100,1000,10000,100000,1000000),trans = 'log10')+
  geom_vline(xintercept = Cutoff2,alpha=0.3)

# all species
NS_ratio_all = NS_ratio_all2
NS_ratio = NS_ratio_all[!duplicated(paste(NS_ratio_all$windowsize,
                                          NS_ratio_all$remove)),c(2,3,4,7)]
for (i in 1:nrow(NS_ratio)){
  temp_window = NS_ratio_all[which(NS_ratio_all$windowsize==
                                     NS_ratio$windowsize[i]&
                                     NS_ratio_all$remove==
                                     NS_ratio$remove[i] ),]
  NS_ratio$N[i] = sum(temp_window$N)
  NS_ratio$S[i] = sum(temp_window$S)
}
NS_ratio$S_cu =NS_ratio$S
NS_ratio$S_cu[which(NS_ratio$S_cu==0)]=0.5
NS_ratio$NS_ratio = NS_ratio$N/NS_ratio$S_cu
Cutoffnew = (Cutoff + Cutoff2)/2

ggplot(NS_ratio,aes(
  x= windowsize,
  y = NS_ratio,
  group = remove))+
  geom_point(aes(
    color=remove,
    fill=remove
  ),
  alpha = 1,
  size=2
  )+
  geom_line(aes(
    color=remove
  ),
  alpha = 0.5,
  size=1
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_x_continuous(breaks = c(1,10,100,1000,10000,100000,1000000),trans = 'log10')+
  geom_vline(xintercept = c(Cutoff,Cutoffnew,Cutoff2),alpha=0.3)
Cutoffnew#9440
# for each species
NS_ratio_all = NS_ratio_all2[which(!(NS_ratio_all2$species
                                     %in% subset_species)),]

NS_ratio = NS_ratio_all[!duplicated(paste(NS_ratio_all$species,paste(NS_ratio_all$windowsize,
                                          NS_ratio_all$remove))),]
for (i in 1:nrow(NS_ratio)){
  temp_window = NS_ratio_all[which(NS_ratio_all$windowsize==
                                     NS_ratio$windowsize[i]&
                                     NS_ratio_all$remove==
                                     NS_ratio$remove[i]&
                                     NS_ratio_all$species==
                                     NS_ratio$species[i] ),]
  NS_ratio$N[i] = sum(temp_window$N)
  NS_ratio$S[i] = sum(temp_window$S)
}
NS_ratio$S_cu =NS_ratio$S
NS_ratio$S_cu[which(NS_ratio$S_cu==0)]=0.5
NS_ratio$NS_ratio = NS_ratio$N/NS_ratio$S_cu


allspecies = unique(NS_ratio$species)
i = 1

{species = allspecies[i]
  NS_ratio2=NS_ratio[which(NS_ratio$species==species),]
ggplot(NS_ratio2,aes(
  x= windowsize,
  y = NS_ratio,
  group = remove))+
  geom_point(aes(
    color=remove,
    fill=remove
  ),
  alpha = 1,
  size=2
  )+
  geom_line(aes(
    color=remove
  ),
  alpha = 0.5,
  size=1
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_x_continuous(breaks = c(1,10,100,1000,10000,100000,1000000),trans = 'log10')
}
species
ggsave(paste(species,'.pdf',sep = ''),path = 'final_SNP/vcf/plot/',
       width = 15, height = 12, units = "in")
i = i + 1

###---- sum up HS table finished----
# sum up metadata
allspecies=read.delim('final_SNP/HSgene/all.species.txt',header=F)
colnames(allspecies)=as.matrix(allspecies[1,])
colnames(allspecies)[1]='X.donor_species'
allspecies=allspecies[-1,]
allspecies_temp = allspecies[which(allspecies$No.genome==''),]
allspecies_temp=as.matrix(allspecies_temp[1,])
allspecies_temp[1,c(1,2)]=c('BaCa_clustercluster1','0')
allspecies = rbind(allspecies,allspecies_temp)
allspecies_temp[1,c(1,2)]=c('FaPr_clustercluster1','0')
allspecies = rbind(allspecies,allspecies_temp)
allspecies_temp[1,c(1,2)]=c('EuRe_clustercluster1','0')
allspecies = rbind(allspecies,allspecies_temp)
library(stringr)
allspecies$species = str_split_fixed(allspecies$X.donor_species, "_", 4)[,1]
allspecies$cluster = str_split_fixed(allspecies$X.donor_species, "_cluster", 2)[,2]

allspecies$species[which(allspecies$species=='1')]=str_split_fixed(allspecies$X.donor_species[which(allspecies$species=='1')], "_", 4)[,2]
unique(allspecies$species)

library(plyr)
allspecies_HS = count(allspecies,c('X.donor_species','High_selected'))
allspecies_HS=allspecies_HS[which(allspecies_HS$High_selected=='True'),]
colnames(allspecies_HS)[3]='HS_NO'
allspecies=merge(allspecies,allspecies_HS[,c(1,3)],
                 by='X.donor_species',all.x=T)
allspecies_HS=allspecies[which(allspecies$High_selected=='True'),c(1,6)]
colnames(allspecies_HS)[2]='freq'
allspecies_HS$freq=as.numeric(as.matrix(allspecies_HS$freq))
allspecies_HS = count(allspecies_HS,c('X.donor_species'))
colnames(allspecies_HS)[2]='HS_SNP_NO'
allspecies=merge(allspecies,allspecies_HS,
                 by='X.donor_species',all.x=T)

neworder=read.delim('species_order.new.txt',header=T)
allspecies2=allspecies
allspecies2=merge(allspecies2,neworder,
                  by='species',all=T)

allspecies2$gene_name = str_split_fixed(allspecies2$gene, "_", 4)[,1]
allspecies1=allspecies2[which(allspecies2$gene_name!='NODE'),]
write.table(allspecies2,'final_SNP/HSgene/all.species.all.txt',
            quote=F,sep='\t',row.names=F)

write.table(allspecies1,'final_SNP/HSgene/all.species.all.short.withHS.txt',
            quote=F,sep='\t',row.names=F)

allspecies1=allspecies2[which(allspecies2$gene_name!='NODE' &
                                !is.na(allspecies2$neworder) &
                                (as.matrix(allspecies2$X.donor_species) !=
                                as.matrix(allspecies2$species) |
                                (allspecies2$gene == '0'))),]
write.table(allspecies1,'final_SNP/HSgene/all.species.lineage.txt',
            quote=F,sep='\t',row.names=F)

allspecies1=allspecies2[which(allspecies2$gene_name!='NODE' &
                                !is.na(allspecies2$neworder) &
                                as.matrix(allspecies2$X.donor_species) ==
                                   as.matrix(allspecies2$species)),]

write.table(allspecies1,'final_SNP/HSgene/all.species.all.species.txt',
            quote=F,sep='\t',row.names=F)

allspecies1=allspecies2[which(allspecies2$gene_name=='NODE'),]
write.table(allspecies1,'final_SNP/HSgene/all.species.all.gene.txt',
            quote=F,sep='\t',row.names=F)

#---- simulate adaptive genes within a species finished----
Species_set = c()
allspecies=read.delim('final_SNP/HSgene/all.species.all.short.withHS.txt',header=T)
allspecies2=read.delim('final_SNP/HSgene/all.species.all.gene.txt',header=T)

allspecies=allspecies[which(allspecies$gene!='0' &
                              allspecies$No.genome!='0' &
                              allspecies$X.donor_species!=
                              as.matrix(allspecies$species)),]
allspecies = allspecies[which(!grepl('highselect',allspecies$gene)),]
library(binom)
library(ggplot2)
library(RColorBrewer)
clonal = read.delim('clonal_pop_gene_num.sum.txt',header=T)
cluster_ratio = read.delim('cluster_ratio.txt',header=T)
clonal = merge(clonal,cluster_ratio,by = 'species',all.x=T)
library(stringr)
clonal=clonal[!duplicated(clonal$cluster),]
allspecies$donor_species = str_split_fixed(allspecies$X.donor_species, ".donor", 2)[,1]
allspecies=merge(allspecies,clonal[,c(2,3,6)],by='donor_species',
                 by.y='cluster',
                 all.x=T)
allspecies_set = unique(allspecies$species)

allspecies_geneonly = allspecies2[which(allspecies2$gene_length !=1000),]
Min_SNP_highselect_cutoff = 1/3000
SNP_cutoff_original = 2
#Max_SNP_highselect_cutoff = 0.02
simround = 100
library(plyr)
simulate_HS <- function(total_gene,total_N,total_S,avg_gene_length,SNP_cutoff){
  mut_N = sample(1:total_gene, total_N, replace=T)# repeats allowed
  mut_S = sample(1:total_gene, total_S, replace=T)
  mut_N = count(mut_N)
  mut_S = count(mut_S)
  mut_N = merge(mut_N,mut_S,by='x',all.x=T)
  mut_N$freq.y[is.na(mut_N$freq.y)]=0
  mut_N$totalSNP_ratio = (mut_N$freq.x + mut_N$freq.y)/avg_gene_length
  return(length(which(mut_N$freq.x >= SNP_cutoff & mut_N$totalSNP_ratio >= Min_SNP_highselect_cutoff )))
}
temp_simulate = matrix(0,nrow=length(allspecies_set)*simround,ncol=10)
k = 1
for(species in allspecies_set)
{
  if (species %in% Species_set)
    SNP_cutoff = 3
  else
    SNP_cutoff = SNP_cutoff_original
  species = toString(species)
  allspecies_sub = allspecies[which(allspecies$species == species),]
  allspecies_sub2 = allspecies_geneonly[which(allspecies_geneonly$X.donor_species %in%
                                                as.matrix(allspecies_sub$X.donor_species)),]
  total_N = sum(allspecies_sub$N)
  total_S = sum(allspecies_sub$S)
  if (nrow(allspecies_sub)>=2)
    {total_gene = sum(as.numeric(as.matrix(allspecies_sub$tota_gene_num))/
                       as.numeric(as.matrix(allspecies_sub$cluster_ratio)))
    avg_gene =sum(as.numeric(as.matrix(allspecies_sub$tota_gene_num)))/nrow(allspecies_sub)
    }else
    {total_gene = allspecies_sub$tota_gene_num[1]
    avg_gene = allspecies_sub$tota_gene_num[1]
    
    }  
  if (nrow(allspecies_sub2)>=5)
    avg_gene_length = sum(allspecies_sub2$gene_length)/nrow(allspecies_sub2)
  else
    avg_gene_length = sum(allspecies_geneonly$gene_length)/nrow(allspecies_geneonly)
  HS_NO = sum(allspecies_sub$HS_NO[!is.na(allspecies_sub$HS_NO)])
  speciesname = toString(species)
  for(i in 1:simround)
  {temp_simulate[k,]=c(speciesname,i,
                       simulate_HS(as.integer(total_gene),as.integer(total_N),
                                   as.integer(total_S),as.integer(avg_gene_length),SNP_cutoff),
                       as.integer(total_gene),as.integer(total_N),
                       as.integer(total_S),as.integer(avg_gene_length),
                       as.integer(HS_NO),nrow(allspecies_sub),avg_gene)
  k = k + 1}
}

colnames(temp_simulate)=c('species','simulate_time','adaptive_genes',
                          'total_gene','total_N','total_S','avg_gene_length','HS_NO','total_lineage_num','avg_gene')
temp_simulate=data.frame(temp_simulate)
temp_simulate=merge(temp_simulate,allspecies[!duplicated(allspecies$species),c(2,32)],
                    by='species',all.x=T)
write.table(temp_simulate,'final_SNP/HSgene/HS_simulation.total.txt',
            quote=F,sep='\t',row.names=F)

#----set adaptive genes cutoff within a species finished----
# Hypergeometric distribution
temp_simulate=read.delim('final_SNP/HSgene/HS_simulation.total.txt',header=T)
temp_simulate=data.frame(temp_simulate)
temp_simulate_sum = temp_simulate[!duplicated(temp_simulate$species),]
# all SNPs in a species
temp_simulate_sum2=matrix(0,nrow=0,ncol=4)
colnames(temp_simulate_sum2)=c(
  'species','total_gene','No_SNP_gene','p'
)
temp_simulate_sum2=data.frame(temp_simulate_sum2)
SNP_gene_cutoff = 2
for (i in 1:nrow(temp_simulate_sum))
{
  m = as.integer(temp_simulate_sum$total_N[i]+temp_simulate_sum$total_S[i]) #SNP bp
  n = as.integer(temp_simulate_sum$total_gene[i]*temp_simulate_sum$avg_gene_length[i]-m)# no-SNP bp
  k = as.integer(temp_simulate_sum$avg_gene_length[i])# gene length
  x <- 0:5
  temp_result = dhyper(x, m, n, k)
  temp_result2 = data.frame(species=temp_simulate_sum$species[i],
                            total_gene=temp_simulate_sum$total_gene[i],
                            No_SNP_gene = x,
                            p=temp_result)
  temp_simulate_sum2=rbind(temp_simulate_sum2,
                           temp_result2)
}
temp_simulate_sum2$gene_num=temp_simulate_sum2$total_gene*temp_simulate_sum2$p
temp_simulate_sum2=merge(temp_simulate_sum2,
                         temp_simulate[!duplicated(temp_simulate$species),c(1,8,9)],
                         by='species',all.x=T)

write.table(temp_simulate_sum2,'final_SNP/HSgene/HS_hypergeometric_distribution.total.txt',
            quote=F,sep='\t',row.names=F)
temp_simulate_sum2=read.delim('final_SNP/HSgene/HS_hypergeometric_distribution.total.txt',header=T)
temp_simulate_sum3=temp_simulate_sum2[which(temp_simulate_sum2$No_SNP_gene == 2),]
temp_simulate_sum3$species[which(temp_simulate_sum3$gene_num>=SNP_gene_cutoff )]
temp_simulate_sum3=temp_simulate_sum2[which(temp_simulate_sum2$No_SNP_gene == 3),]
temp_simulate_sum3$species[which(temp_simulate_sum3$gene_num>=SNP_gene_cutoff )]
temp_simulate_sum3=temp_simulate_sum2[which(temp_simulate_sum2$No_SNP_gene == 4),]
temp_simulate_sum3$species[which(temp_simulate_sum3$gene_num>=SNP_gene_cutoff )]

library(ggplot2)
library(RColorBrewer)
ggplot(temp_simulate_sum2,
       aes(y=gene_num,
           x=No_SNP_gene,
           group = species))+
  geom_point(aes(
    color='1',
    fill = '1'
  ),
  size=2,
  alpha = 0.5
  )+
  geom_line(aes(
    color='1'
  ),
  size=1,
  alpha = 0.5
  )+ 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_y_continuous(breaks = c(1e-10,1e-5,1e-1,1,10,100,1000,10000,100000),trans = 'log10')+
  scale_x_continuous(breaks = x)+
  geom_hline(aes(yintercept = SNP_gene_cutoff,color = '2')) 

###---- sum up final HS table finished----
# sum up metadata
allspecies=read.delim('final_SNP/HSgene/all.species.txt.High_select2.txt',header=F)
colnames(allspecies)=as.matrix(allspecies[4,])
colnames(allspecies)[1]='X.donor_species'
allspecies=allspecies[-4,]
allspecies_temp = allspecies[which(allspecies$No.genome==''),]
allspecies_temp=as.matrix(allspecies_temp[1,])
allspecies_temp[1,c(1,2)]=c('BaCa_clustercluster1','0')
allspecies = rbind(allspecies,allspecies_temp)
allspecies_temp[1,c(1,2)]=c('FaPr_clustercluster1','0')
allspecies = rbind(allspecies,allspecies_temp)
allspecies_temp[1,c(1,2)]=c('EuRe_clustercluster1','0')
allspecies = rbind(allspecies,allspecies_temp)
library(stringr)
allspecies$species = str_split_fixed(allspecies$X.donor_species, "_", 4)[,1]
allspecies$cluster = str_split_fixed(allspecies$X.donor_species, "_cluster", 2)[,2]

allspecies$species[which(allspecies$species=='1')]=str_split_fixed(allspecies$X.donor_species[which(allspecies$species=='1')], "_", 4)[,2]
unique(allspecies$species)

library(plyr)
allspecies_HS = count(allspecies,c('X.donor_species','High_selected'))
allspecies_HS = allspecies_HS[which(allspecies_HS$High_selected == 'True'),]

colnames(allspecies_HS)[3]='HS_NO'
allspecies=merge(allspecies,allspecies_HS[,c(1,3)],
                 by='X.donor_species',all.x=T)

allspecies_HS = allspecies[which(allspecies$High_selected == 'True'),]
allspecies_HS$contiglength = str_split_fixed(allspecies_HS$gene, "length_", 2)[,2]
allspecies_HS$contiglength = str_split_fixed(allspecies_HS$contiglength, "_cov", 2)[,1]
allspecies_HS$contiglength[which(allspecies_HS$contiglength=='')]='250468'
allspecies_HS2=allspecies_HS[which(as.numeric(as.matrix(allspecies_HS$contiglength))>=5000),]

allspecies_HS=allspecies_HS2[,c(1,6)]
colnames(allspecies_HS)[2]='freq'
allspecies_HS$freq=as.numeric(as.matrix(allspecies_HS$freq))
allspecies_HS = count(allspecies_HS,c('X.donor_species'))
colnames(allspecies_HS)[2]='HS_SNP_NO'
allspecies=merge(allspecies,allspecies_HS,
                 by='X.donor_species',all.x=T)

neworder=read.delim('species_order.new.txt',header=T)
allspecies2=allspecies
allspecies2=merge(allspecies2,neworder,
                  by='species',all=T)

allspecies2$gene_name = str_split_fixed(allspecies2$gene, "_", 4)[,1]
allspecies1=allspecies2[which(allspecies2$gene_name!='NODE'),]
write.table(allspecies2,'final_SNP/HSgene/all.species.High_select2.all.txt',
            quote=F,sep='\t',row.names=F)

write.table(allspecies1,'final_SNP/HSgene/all.species.High_select2.all.short.withHS.txt',
            quote=F,sep='\t',row.names=F)

allspecies1=allspecies2[which(allspecies2$gene_name!='NODE' &
                                !is.na(allspecies2$neworder) &
                                (as.matrix(allspecies2$X.donor_species) !=
                                   as.matrix(allspecies2$species) |
                                   (allspecies2$gene == '0'))),]
write.table(allspecies1,'final_SNP/HSgene/all.species.High_select2.lineage.txt',
            quote=F,sep='\t',row.names=F)

allspecies1=allspecies2[which(allspecies2$gene_name!='NODE' &
                                !is.na(allspecies2$neworder) &
                                as.matrix(allspecies2$X.donor_species) ==
                                as.matrix(allspecies2$species)),]

write.table(allspecies1,'final_SNP/HSgene/all.species.High_select2.all.species.txt',
            quote=F,sep='\t',row.names=F)

allspecies1=allspecies2[which(allspecies2$gene_name!='NODE' &
                                is.na(allspecies2$neworder) &
                                as.matrix(allspecies2$X.donor_species) ==
                                as.matrix(allspecies2$species)),]

write.table(allspecies1,'final_SNP/HSgene/all.species.High_select2.all.genus.txt',
            quote=F,sep='\t',row.names=F)

allspecies1=allspecies2[which(allspecies2$gene_name=='NODE'),]
write.table(allspecies1,'final_SNP/HSgene/all.species.High_select2.all.gene.txt',
            quote=F,sep='\t',row.names=F)

allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.all.short.withHS.txt',header=T)
allspecies=allspecies[which(allspecies$X.donor_species!=
                              as.matrix(allspecies$species)),]
colnames(allspecies)
library(stringr)
library(plyr)
library(ggplot2)
library(RColorBrewer)
# dmrca
dmrca=read.delim("final_SNP/dmrca/alldmrca.txt",header=T)
dmrca$donor_species = paste(str_split_fixed(str_split_fixed(dmrca$donor_species, "_subcluster", 2)[,1]
                                            , ".donor", 2)[,1],
                            str_split_fixed(dmrca$donor_species, "donor", 2)[,2],sep='.donor')
allspecies$dmrca = 0
for (i in 1:nrow(allspecies))
{
  tempdmrca = dmrca[which(dmrca$donor_species==toString(allspecies$X.donor_species[i])),]
  allspecies$dmrca[i]=mean(tempdmrca$dmrca)
}

write.table(allspecies,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.txt',
            quote=F,sep='\t',row.names=F)
allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.all.short.withHS.txt',header=T)
allspecies=allspecies[which(allspecies$X.donor_species==
                              as.matrix(allspecies$species) &
                                          !is.na(allspecies$species_short)),]
allspecies$No.SNP = as.numeric(as.matrix(allspecies$No.SNP))

allspecies$No.SNP[which(is.na(allspecies$No.SNP))]=0
allspecies=allspecies[rev(order(allspecies$No.SNP)),]
allspecies=allspecies[!duplicated(paste(allspecies$gene)),]
colnames(allspecies)
library(stringr)
library(plyr)
library(ggplot2)
library(RColorBrewer)
# dmrca
dmrca=read.delim("final_SNP/dmrca/alldmrca.txt",header=T)
allspecies$dmrca = 0
dmrca$donor_species = paste(str_split_fixed(str_split_fixed(dmrca$donor_species, "_subcluster", 2)[,1]
                                            , ".donor", 2)[,1],
                            str_split_fixed(dmrca$donor_species, "donor", 2)[,2],sep='.donor')

for (i in 1:nrow(allspecies))
{
  tempdmrca = dmrca[which(grepl(toString(allspecies$X.donor_species[i]),dmrca$donor_species)),]
  allspecies$dmrca[i]=mean(tempdmrca$dmrca)
}

write.table(allspecies,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.species.txt',
            quote=F,sep='\t',row.names=F)

allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.txt',header=T)
allspecies2=read.delim('final_SNP/HSgene/all.species.High_select2.all.gene.txt',header=T)
allspecies$highselect = grepl('highselect',allspecies$gene)
allspecies = allspecies[, c(1,2,4:14,29:35)]
allspecies$No.SNP[is.na(allspecies$No.SNP)]=0
allspecies$HS_NO[is.na(allspecies$HS_NO)]=0
allspecies$HS_SNP_NO[is.na(allspecies$HS_SNP_NO)]=0
allspecies$dmrca[is.na(allspecies$dmrca)]=0
allspecies$No.genomelog10 = log10(allspecies$No.genome)
allspecies$donor =str_split_fixed(allspecies$X.donor_species, "donor.", 2)[,2]
allspecies$genus =str_split_fixed(allspecies$species_short, "_", 2)[,1]
allspecies2_count=count(allspecies2,'X.donor_species')
colnames(allspecies2_count)[2]='No.gene_SNPs'
allspecies=merge(allspecies,allspecies2_count,
                 by='X.donor_species',all.x=T)
# normalization
# Hypergeometric distribution for lineage
temp_simulate=read.delim('final_SNP/HSgene/HS_simulation.total.txt',header=T)
temp_simulate=data.frame(temp_simulate)
temp_simulate_sum = temp_simulate[!duplicated(temp_simulate$species),]
allspecies_noHS = allspecies[which(allspecies$highselect==FALSE),]
allspecies_noHS=merge(allspecies_noHS,temp_simulate_sum,by='species',all.y=T)
# all SNPs in a lineage
temp_simulate_sum2=matrix(0,nrow=0,ncol=5)
SNP_gene_cutoff = 2
colnames(temp_simulate_sum2)=c(
  'X.donor_species','total_gene','No_SNP_gene','p','species'
)
temp_simulate_sum2=data.frame(temp_simulate_sum2)
for (i in 1:nrow(allspecies_noHS))
{
  m = as.integer(allspecies_noHS$N[i]+allspecies_noHS$S[i]) #SNP bp
  n = as.integer((allspecies_noHS$avg_gene[i]*allspecies_noHS$avg_gene_length[i]-m))# no-SNP bp
  k = as.integer(allspecies_noHS$avg_gene_length[i])# gene length
  x <- 0:5
  temp_result = dhyper(x, m, n, k)
  temp_result2 = data.frame(X.donor_species=allspecies_noHS$X.donor_species[i],
                            total_gene=allspecies_noHS$avg_gene[i],
                            No_SNP_gene = x,
                            p=temp_result,
                            species = allspecies_noHS$species[i])
  temp_simulate_sum2=rbind(temp_simulate_sum2,
                           temp_result2)
}
temp_simulate_sum2$gene_num=temp_simulate_sum2$total_gene*temp_simulate_sum2$p
write.table(temp_simulate_sum2,'final_SNP/HSgene/HS_hypergeometric_distribution.total.lineage.txt',
            quote=F,sep='\t',row.names=F)

accumulate_genenum <- function(hyper_lineage2,hyper_lineage)
{
  for(i in 1:nrow(hyper_lineage2))
  {
    sub_hyper_lineage = hyper_lineage[which(hyper_lineage$species == hyper_lineage2$species[i] &
                                              hyper_lineage$No_SNP_gene >=  hyper_lineage2$No_SNP_gene[i]),]
    hyper_lineage2$gene_num[i] = sum(sub_hyper_lineage$gene_num)
    hyper_lineage2$SNP_num[i] = sum(sub_hyper_lineage$gene_num*sub_hyper_lineage$No_SNP_gene)
  }
  return(hyper_lineage2)
}
accumulate_genenum_lineage <- function(hyper_lineage2,hyper_lineage)
{
  for(i in 1:nrow(hyper_lineage2))
  {
    sub_hyper_lineage = hyper_lineage[which(hyper_lineage$X.donor_species == 
                                              hyper_lineage2$X.donor_species[i] &
                                              hyper_lineage$No_SNP_gene >=  
                                              hyper_lineage2$No_SNP_gene[i]),]
    hyper_lineage2$gene_num[i] = sum(sub_hyper_lineage$gene_num)
    hyper_lineage2$SNP_num[i] = sum(sub_hyper_lineage$gene_num*sub_hyper_lineage$No_SNP_gene)
  }
  return(hyper_lineage2)
}

hyper_lineage=read.delim('final_SNP/HSgene/HS_hypergeometric_distribution.total.lineage.txt',header=T)
hyper_species=read.delim('final_SNP/HSgene/HS_hypergeometric_distribution.total.txt',header=T)
species_cutoff = read.delim('total_SNP_cutoff_species.txt',header=F)
hyper_lineage2 = hyper_lineage[which(hyper_lineage$No_SNP_gene==2),]
hyper_lineage2$SNP_num = 0
hyper_lineage2 = accumulate_genenum_lineage(hyper_lineage2,hyper_lineage)
hyper_species2 = hyper_species[which(hyper_species$No_SNP_gene==3 & 
                                       hyper_species$species %in% species_cutoff$V1),]
hyper_species2=rbind(hyper_species2,
                     hyper_species[which(hyper_species$No_SNP_gene==2 & 
                                           !(hyper_species$species %in% species_cutoff$V1)),])
hyper_species2$SNP_num = 0
hyper_species2 = accumulate_genenum(hyper_species2,hyper_species)

hyper_lineage = hyper_lineage2
hyper_lineage = merge(hyper_lineage,hyper_species2[,c(1,5,8)],
                      by='species')
hyper_lineage$sim_HS = hyper_lineage$gene_num.x+hyper_lineage$gene_num.y
hyper_lineage$sim_SNP = hyper_lineage$SNP_num.x+hyper_lineage$SNP_num.y

allspecies=merge(allspecies,
                 hyper_lineage[,c(2,10,11)],by='X.donor_species')
allspecies$HS_norm_gene = allspecies$HS_NO/(
  allspecies$HS_NO+
    allspecies$sim_HS)
allspecies$HS_norm_SNP = allspecies$HS_SNP_NO/(
  allspecies$HS_SNP_NO+
    allspecies$sim_SNP)
write.table(allspecies,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.txt',
            quote=F,sep='\t',row.names=F)

allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.txt',header=T)
allspecies_HS = allspecies[which(allspecies$highselect==TRUE),]
allspecies = allspecies[which(allspecies$highselect==FALSE),]
SNP_set = read.delim('final_SNP/HSgene/SNP_set.txt',header=F)
colnames(SNP_set)=c('X.donor_species','SNP_set')
HS_pop = read.delim('final_SNP/HSgene/HS_withinpop.txt',header=T)


# check data
allspecies$genomelog10=log10(as.numeric(as.matrix(allspecies$No.genome)))
allspecies$dmrcalog10=log10(as.numeric(as.matrix(allspecies$dmrca)))
allspecies$No.SNPlog10=log10(as.numeric(as.matrix(allspecies$No.SNP))+1)
allspecies$No.HS_SNP_NOlog10=log10(as.numeric(as.matrix(allspecies$HS_SNP_NO))+1)
allspecies$S_cut = as.numeric(as.matrix(allspecies$S))
allspecies$S_cut[which(allspecies$S_cut==0)]=0.5
allspecies$dNdS_cut = as.numeric(as.matrix(allspecies$dNdS))
allspecies$dNdS_cut[which(allspecies$S_cut==0.5)]=allspecies$N[which(allspecies$S_cut==0.5)]/
  allspecies$S_cut[which(allspecies$S_cut==0.5)]/
  as.numeric(as.matrix(allspecies$expected_ratio[which(allspecies$S_cut==0.5)]))
allspecies$HS_ratio = allspecies$HS_NO/allspecies$No.gene_SNPs
allspecies$HS_ratio2 = allspecies$HS_SNP_NO/allspecies$No.SNP
allspecies_HS$S_cut = as.numeric(as.matrix(allspecies_HS$S))
allspecies_HS$S_cut[which(allspecies_HS$S_cut==0)]=0.5
allspecies_HS$dNdS_cut = as.numeric(as.matrix(allspecies_HS$dNdS))
allspecies_HS$dNdS_cut[which(allspecies_HS$S_cut==0.5)]=allspecies_HS$N[which(allspecies_HS$S_cut==0.5)]/
  allspecies_HS$S_cut[which(allspecies_HS$S_cut==0.5)]/
  as.numeric(as.matrix(allspecies_HS$expected_ratio[which(allspecies_HS$S_cut==0.5)]))
allspecies_HS=allspecies_HS[,which(colnames(allspecies_HS) %in% c(
  'X.donor_species','dNdS_cut'
)),]
colnames(allspecies_HS)[2]='dNdS_cut_highselect'
allspecies=merge(allspecies,allspecies_HS,by='X.donor_species',all.x=T)
allspecies$short_residence = FALSE
allspecies=merge(allspecies,SNP_set,by='X.donor_species',all.x=T)
allspecies$SNP_set = allspecies$SNP_set
allspecies$HS_SNP_set=allspecies$HS_SNP_NO/(allspecies$SNP_set)^0.5
allspecies$No.SNP_set=allspecies$No.SNP/(allspecies$SNP_set)^0.5
allspecies$HS_ratio3 = allspecies$HS_SNP_set/allspecies$No.SNP_set
allspecies$SNP_set2 = (allspecies$SNP_set)^0.5
allspecies=merge(allspecies,HS_pop,by='X.donor_species',all.x=T)
write.table(allspecies,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.txt',
            quote=F,sep='\t',row.names=F)
#
#----statistic test dNdS confidence finished figure 1a----
allspecies1 = read.delim('final_SNP/HSgene/all.species.High_select2.all.txt',header=T)
allspecies1=allspecies1[which(allspecies1$X.donor_species=='allspecies'),]
# dN dS
library(binom)
Sum_table = matrix(0,nrow=0,ncol=4)
colnames(Sum_table)=c('selected','dNdS','lower','higher')
Sum_table2 = matrix(0,nrow=0,ncol=5)
colnames(Sum_table2)=c('selected','p-value','conf.int1','conf.int2','alternative')

unique(allspecies1$gene)
# all genes de novo mutation
N_total = as.numeric(as.matrix(allspecies1$N[which(allspecies1$gene=='allspecies')]))
S_total =as.numeric(as.matrix(allspecies1$S[which(allspecies1$gene=='allspecies')]))
Total_SNP = N_total + S_total
expected_ratio = as.numeric(as.matrix(allspecies1$expected_ratio[which(allspecies1$gene=='allspecies')]))
observed_ratio=N_total/S_total	 
Pall = N_total/Total_SNP
Pref = 1/(1+1/expected_ratio)

result=binom.confint(N_total, (N_total+S_total), 
                     conf.level = 0.95, 
                     methods = "all")
dNdS_low = (Total_SNP*mean(result$lower)/(Total_SNP*(1-mean(result$lower))))/expected_ratio
dNdS_high = (Total_SNP*mean(result$upper)/(Total_SNP*(1-mean(result$upper))))/expected_ratio
dNdS = (Total_SNP*mean(result$mean)/(Total_SNP*(1-mean(result$mean))))/expected_ratio
c(dNdS_low,dNdS_high,dNdS) #0.7871940 0.8506457 0.8182684
Sum_table=rbind(Sum_table,c('de_novo_mutation',dNdS,dNdS_low,dNdS_high))
result=binom.test(N_total,Total_SNP,Pref, alternative="two.sided")
Sum_table2=rbind(Sum_table2,c('All',result$p.value,result$conf.int[1],result$conf.int[2],result$alternative))


# highly selected genes
N_total = as.numeric(as.matrix(allspecies1$N[which(allspecies1$gene=='allspecies_highselect')]))#908
S_total = as.numeric(as.matrix(allspecies1$S[which(allspecies1$gene=='allspecies_highselect')])) #109
Total_SNP = N_total + S_total
expected_ratio = as.numeric(as.matrix(allspecies1$expected_ratio[which(allspecies1$gene=='allspecies_highselect')]))
observed_ratio=N_total/S_total 
PHS = N_total/Total_SNP
result=binom.test(N_total,Total_SNP,Pref, alternative="greater")
Sum_table2=rbind(Sum_table2,c('Highly_selected',result$p.value,result$conf.int[1],result$conf.int[2],result$alternative))

result=binom.confint(N_total, (N_total+S_total), 
                     conf.level = 0.95, 
                     methods = "all")
dNdS_low = (Total_SNP*mean(result$lower)/(Total_SNP*(1-mean(result$lower))))/expected_ratio
dNdS_high = (Total_SNP*mean(result$upper)/(Total_SNP*(1-mean(result$upper))))/expected_ratio
dNdS = (Total_SNP*mean(result$mean)/(Total_SNP*(1-mean(result$mean))))/expected_ratio
c(dNdS_low,dNdS_high,dNdS) #5.594611 8.044727 6.683171
Sum_table=rbind(Sum_table,c('Highly_selected',dNdS,dNdS_low,dNdS_high))

# flexible genes
N_total = as.numeric(as.matrix(allspecies1$N[which(allspecies1$gene=='allspecies_flexible')]))
S_total = as.numeric(as.matrix(allspecies1$S[which(allspecies1$gene=='allspecies_flexible')]))
Total_SNP = N_total + S_total
expected_ratio = as.numeric(as.matrix(allspecies1$expected_ratio[which(allspecies1$gene=='allspecies_flexible')]))
observed_ratio=N_total/S_total	 
Pfle = N_total/Total_SNP
result=binom.test(N_total,Total_SNP,Pref,alternative = 'two.sided')
Sum_table2=rbind(Sum_table2,c('flexible',result$p.value,
                              result$conf.int[1],result$conf.int[2],result$alternative))
result=binom.confint(N_total, (N_total+S_total), 
                     conf.level = 0.95, 
                     methods = "all")
dNdS_low = (Total_SNP*mean(result$lower)/(Total_SNP*(1-mean(result$lower))))/expected_ratio
dNdS_high = (Total_SNP*mean(result$upper)/(Total_SNP*(1-mean(result$upper))))/expected_ratio
dNdS = (Total_SNP*mean(result$mean)/(Total_SNP*(1-mean(result$mean))))/expected_ratio
c(dNdS_low,dNdS_high,dNdS) #0.7871940 0.8506457 0.8182684
Sum_table=rbind(Sum_table,c('flexible',dNdS,dNdS_low,dNdS_high))

# core genes
N_total = as.numeric(as.matrix(allspecies1$N[which(allspecies1$gene=='allspecies_core')]))
S_total = as.numeric(as.matrix(allspecies1$S[which(allspecies1$gene=='allspecies_core')]))
Total_SNP = N_total + S_total
expected_ratio = as.numeric(as.matrix(allspecies1$expected_ratio[which(allspecies1$gene=='allspecies_core')]))
observed_ratio=N_total/S_total	 
Pcore = N_total/Total_SNP
result=binom.test(N_total,Total_SNP,Pref,alternative = 'less')
Sum_table2=rbind(Sum_table2,c('core',result$p.value,result$conf.int[1],result$conf.int[2],result$alternative))
result=binom.confint(N_total, (N_total+S_total), 
                     conf.level = 0.95, 
                     methods = "all")
dNdS_low = (Total_SNP*mean(result$lower)/(Total_SNP*(1-mean(result$lower))))/expected_ratio
dNdS_high = (Total_SNP*mean(result$upper)/(Total_SNP*(1-mean(result$upper))))/expected_ratio
dNdS = (Total_SNP*mean(result$mean)/(Total_SNP*(1-mean(result$mean))))/expected_ratio
c(dNdS_low,dNdS_high,dNdS) #0.7871940 0.8506457 0.8182684
Sum_table=rbind(Sum_table,c('core',dNdS,dNdS_low,dNdS_high))

# plot
Sum_table=data.frame(Sum_table)
Sum_table2=data.frame(Sum_table2)
Sum_table2$q.value = p.adjust(as.matrix(Sum_table2$p.value), method = 'fdr')
write.table(Sum_table2,'final_SNP/HSgene/dNdS.binomtest.txt',
            quote=F,sep='\t',row.names=F)
write.table(Sum_table,'final_SNP/HSgene/dNdS.confi.txt',
            quote=F,sep='\t',row.names=F)
Sum_table=read.delim('final_SNP/HSgene/dNdS.confi.txt',header=T)
library(ggplot2)
library(RColorBrewer)
ggplot(Sum_table,aes(x=selected,
                     y=as.numeric(as.matrix(dNdS))
))+
  #y=..scaled..,
  geom_bar(aes(
    color=selected,
    fill=selected
  ),
  stat = 'identity',
  size=0,
  alpha = 0.5,
  width=0.9
  )+ 
  geom_errorbar(aes(ymin=as.numeric(as.matrix(lower)), 
                    ymax=as.numeric(as.matrix(higher)),
                    color=selected), 
                width=0.6,
                position=position_dodge(.9))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  geom_hline(yintercept = 1,alpha=0.3)+
  scale_y_continuous(labels = c(0,1,2,3,4,5,6,7,8),
                     breaks = c(0,1,2,3,4,5,6,7,8))+
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(4,5,2,1)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(4,5,2,1)])
#
#---- simulate adaptive SNPs within a species after tuning cutoff finished----
#total_SNP_cutoff = read.delim('total_SNP_cutoff.txt',header=F)
total_SNP_cutoff = c()
total_SNP_cutoff_species = read.delim('total_SNP_cutoff_species.txt',header=F)
allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.all.short.withHS.txt',header=T)
allspecies2=read.delim('final_SNP/HSgene/all.species.High_select2.all.gene.txt',header=T)

allspecies=allspecies[which(allspecies$gene!='0' &
                              allspecies$No.genome!='0' &
                              allspecies$X.donor_species!=
                              as.matrix(allspecies$species)),]
allspecies = allspecies[which(!grepl('highselect',allspecies$gene)),]
library(binom)
library(ggplot2)
library(RColorBrewer)
clonal = read.delim('clonal_pop_gene_num.sum.txt',header=T)
cluster_ratio = read.delim('cluster_ratio.txt',header=T)
clonal = merge(clonal,cluster_ratio,by = 'species',all.x=T)
library(stringr)
clonal=clonal[!duplicated(clonal$cluster),]
allspecies$donor_species = str_split_fixed(allspecies$X.donor_species, ".donor", 2)[,1]
allspecies=merge(allspecies,clonal[,c(2,3,6)],by='donor_species',
                 by.y='cluster',
                 all.x=T)
allspecies_set = unique(allspecies$species)

allspecies_geneonly = allspecies2[which(allspecies2$gene_length !=1000),]
Min_SNP_highselect_cutoff = 1/3000
SNP_cutoff_original = 2
cluster_ratio = 6604/4086
#Max_SNP_highselect_cutoff = 0.02
simround = 1000
library(plyr)
simulate_HS <- function(total_gene,total_N,total_S,avg_gene_length,SNP_cutoff){
  mut_N = sample(1:total_gene, total_N, replace=T)# repeats allowed
  mut_S = sample(1:total_gene, total_S, replace=T)
  mut_N = count(mut_N)
  mut_S = count(mut_S)
  mut_N = merge(mut_N,mut_S,by='x',all.x=T)
  mut_N$freq.y[is.na(mut_N$freq.y)]=0
  mut_N$totalSNP_ratio = (mut_N$freq.x + mut_N$freq.y)/avg_gene_length
  mut_N_HS = mut_N[which(mut_N$freq.x >= SNP_cutoff & 
                mut_N$totalSNP_ratio >= Min_SNP_highselect_cutoff ),]
  return(c(sum(mut_N_HS$freq.x+mut_N_HS$freq.y),
           nrow(mut_N_HS),sum(mut_N_HS$freq.x),sum(mut_N_HS$freq.y)))
}
temp_simulate = matrix(0,nrow=length(allspecies_set)*simround,ncol=17)
k = 1
for(species in allspecies_set)
{
  if (species %in% total_SNP_cutoff_species$V1)
    SNP_cutoff_species = total_SNP_cutoff_species$V2[which(total_SNP_cutoff_species$V1 == species)]
  else
    SNP_cutoff_species = SNP_cutoff_original
  if (species %in% total_SNP_cutoff$V1)
    SNP_cutoff = total_SNP_cutoff$V2[which(total_SNP_cutoff$V1 == species)]
  else
    SNP_cutoff = SNP_cutoff_original
  species = toString(species)
  allspecies_sub = allspecies[which(allspecies$species == species),]
  allspecies_sub2 = allspecies_geneonly[which(allspecies_geneonly$X.donor_species %in%
                                                as.matrix(allspecies_sub$X.donor_species)),]
  total_N = sum(allspecies_sub$N)
  total_S = sum(allspecies_sub$S)
  total_gene = sum(as.numeric(as.matrix(allspecies_sub$tota_gene_num))/
                     as.numeric(as.matrix(allspecies_sub$cluster_ratio)))
  avg_gene =sum(as.numeric(as.matrix(allspecies_sub$tota_gene_num)))/nrow(allspecies_sub)
  
  if (nrow(allspecies_sub2)>=5)
    avg_gene_length = sum(allspecies_sub2$gene_length)/nrow(allspecies_sub2)
  else
    avg_gene_length = sum(allspecies_geneonly$gene_length)/nrow(allspecies_geneonly)
  HS_NO = sum(allspecies_sub$HS_NO[!is.na(allspecies_sub$HS_NO)])
  speciesname = toString(species)
  lineage_num = nrow(allspecies_sub)
  N_average = as.integer(total_N)/lineage_num
  S_average = as.integer(total_S)/lineage_num
  print(species)
  print(SNP_cutoff_species)
  print(SNP_cutoff)
  for(i in 1:simround)
  {temp_simulate[k,]=c(speciesname,i,
                       simulate_HS(as.integer(total_gene),as.integer(total_N),
                                   as.integer(total_S),as.integer(avg_gene_length),SNP_cutoff_species),
                       simulate_HS(as.integer(avg_gene),as.integer(N_average),
                                   as.integer(S_average),as.integer(avg_gene_length),SNP_cutoff),
                       as.integer(total_gene),as.integer(total_N),
                       as.integer(total_S),as.integer(avg_gene_length),
                       as.integer(HS_NO),lineage_num,avg_gene)
  k = k + 1}
}

colnames(temp_simulate)=c('species','simulate_time','adaptive_genes_species_SNPs','adaptive_genes_species',
                          'adaptive_genes_species_N','adaptive_genes_species_S',
                          'adaptive_genes_lineage_SNPs','adaptive_genes_lineage',
                          'adaptive_genes_lineage_N','adaptive_genes_lineage_S',
                          'total_gene','total_N','total_S','avg_gene_length','HS_NO','total_lineage_num','avg_gene')
temp_simulate=data.frame(temp_simulate)
temp_simulate=merge(temp_simulate,allspecies[!duplicated(allspecies$species),c(2,33)],
                    by='species',all.x=T)
write.table(temp_simulate,'final_SNP/HSgene/HS_simulation.total.afterturningcutoff.clusterratio.withSNPs.txt',
            quote=F,sep='\t',row.names=F)

temp_simulate = read.delim('final_SNP/HSgene/HS_simulation.total.afterturningcutoff.clusterratio.withSNPs.txt',
                           header=T)
temp_simulate$adaptive_genes = temp_simulate$adaptive_genes_species+temp_simulate$adaptive_genes_lineage
temp_simulate2 = data.frame(
  neworder = unique(temp_simulate$neworder),
  lower = 0,
  higher = 0,
  Mean = 0
)
for(i in 1:nrow(temp_simulate2)){
  neworder=temp_simulate2$neworder[i]
  temp_simulate_sub = temp_simulate[which(temp_simulate$neworder == neworder),]
  temp_simulate2$lower[i]=quantile(temp_simulate_sub$adaptive_genes,c(0.05))
  temp_simulate2$higher[i]=quantile(temp_simulate_sub$adaptive_genes,c(0.95))
  temp_simulate2$Mean[i]=mean(temp_simulate_sub$adaptive_genes)
}
temp_simulate2$lower[which(temp_simulate2$lower<0)]=0
ggplot(temp_simulate2,aes(x=neworder,
                          y=as.numeric(as.matrix(Mean))
))+
  #y=..scaled..,
  geom_point(aes(
    color='1',
    fill='1'
  ),
  size=2,
  alpha = 1
  )+ 
  geom_errorbar(aes(ymin=as.numeric(as.matrix(lower)), 
                    ymax=as.numeric(as.matrix(higher)),
                    color='1'), 
                width=0.6,
                position=position_dodge(.9))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(values =c(brewer.pal(5, "PuOr")[c(1)]))+
  scale_color_manual(values =c(brewer.pal(5, "PuOr")[c(1)]))+
  scale_x_continuous(limits = c(0,37),breaks = c(1,36),labels = c(1,36))

#---- simulate adaptive SNPs within a lineage after tuning cutoff finished----
total_SNP_cutoff = c()
total_SNP_cutoff_species = read.delim('total_SNP_cutoff_species.txt',header=F)
allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.all.short.withHS.txt',header=T)
allspecies2=read.delim('final_SNP/HSgene/all.species.High_select2.all.gene.txt',header=T)
allspecies=allspecies[which(allspecies$gene!='0' &
                              allspecies$No.genome!='0' &
                              allspecies$X.donor_species!=
                              as.matrix(allspecies$species)),]
allspecies = allspecies[which(!grepl('highselect',allspecies$gene)),]
library(binom)
library(ggplot2)
library(RColorBrewer)
clonal = read.delim('clonal_pop_gene_num.sum.txt',header=T)
cluster_ratio = read.delim('cluster_ratio.txt',header=T)
clonal = merge(clonal,cluster_ratio,by = 'species',all.x=T)
library(stringr)
clonal=clonal[!duplicated(clonal$cluster),]
allspecies$donor_species = str_split_fixed(allspecies$X.donor_species, ".donor", 2)[,1]
allspecies=merge(allspecies,clonal[,c(2,3,6)],by='donor_species',
                 by.y='cluster',
                 all.x=T)
allspecies_set = unique(allspecies$species)
allspecies_geneonly = allspecies2[which(allspecies2$gene_length !=1000),]
Min_SNP_highselect_cutoff = 1/3000
SNP_cutoff_original = 2
cluster_ratio = 6604/4086
#Max_SNP_highselect_cutoff = 0.02
simround = 100
library(plyr)
simulate_HS <- function(total_gene,total_N,total_S,avg_gene_length,SNP_cutoff){
  mut_N = sample(1:total_gene, total_N, replace=T)# repeats allowed
  mut_S = sample(1:total_gene, total_S, replace=T)
  mut_N = count(mut_N)
  mut_S = count(mut_S)
  mut_N = merge(mut_N,mut_S,by='x',all.x=T)
  mut_N$freq.y[is.na(mut_N$freq.y)]=0
  mut_N$totalSNP_ratio = (mut_N$freq.x + mut_N$freq.y)/avg_gene_length
  mut_N_HS = mut_N[which((mut_N$freq.x+mut_N$freq.y) >= SNP_cutoff & 
                           mut_N$totalSNP_ratio >= Min_SNP_highselect_cutoff ),]
  # No. HS SNPs, No. HS genes, No. HS N, No. HS S 
  return(c(sum(mut_N_HS$freq.x+mut_N_HS$freq.y),
           nrow(mut_N_HS),sum(mut_N_HS$freq.x),sum(mut_N_HS$freq.y)))
}
sample_species <- function(Num,Name){
  temp_name = c()
  for(i in 1:length(Num))
  {
   temp_name=c(temp_name, rep(toString(Name[i]),Num[i]))
  }
  return(temp_name)
}
species_HS <- function(temp_species_HS,SNP_cutoff_species,allspecies_sub)
{
  # No. HS SNPs, No. HS genes, No. HS N, No. HS S 
  temp_lineage = sample(sample_species(allspecies_sub$No.SNP,
                                       allspecies_sub$X.donor_species),
         temp_species_HS[1],replace=FALSE)
  temp_lineage_N = sample(temp_lineage,
                          temp_species_HS[3],
                          replace=FALSE)
  temp_lineage_count = count(temp_lineage)
  temp_lineage_count_N = count(temp_lineage_N)
  temp_lineage_count_gene = temp_lineage_count
  temp_lineage_count=merge(temp_lineage_count,temp_lineage_count_gene,
                           by='x',all=T)
  temp_lineage_count=merge(temp_lineage_count,temp_lineage_count_N,
                           by='x',all=T)
  colnames(temp_lineage_count)=c('X.donor_species','SNP','gene','N')
  temp_lineage_count$N[is.na(temp_lineage_count$N)]=0
  temp_lineage_count$S = temp_lineage_count$SNP -  temp_lineage_count$N
  return(temp_lineage_count)
}
temp_simulate = matrix(0,nrow=length(unique(allspecies$X.donor_species))*simround,ncol=18)
k = 1
for(species in allspecies_set)
{
  if (species %in% total_SNP_cutoff_species$V1)
    SNP_cutoff_species = total_SNP_cutoff_species$V2[which(total_SNP_cutoff_species$V1 == species)]
  else
    SNP_cutoff_species = SNP_cutoff_original
  if (species %in% total_SNP_cutoff$V1)
    SNP_cutoff = total_SNP_cutoff$V2[which(total_SNP_cutoff$V1 == species)]
  else
    SNP_cutoff = SNP_cutoff_original
  species = toString(species)
  allspecies_sub = allspecies[which(allspecies$species == species),]
  allspecies_sub2 = allspecies_geneonly[which(allspecies_geneonly$X.donor_species %in%
                                                as.matrix(allspecies_sub$X.donor_species)),]
  total_N = sum(allspecies_sub$N)
  total_S = sum(allspecies_sub$S)
  total_gene = sum(as.numeric(as.matrix(allspecies_sub$tota_gene_num))/
                     as.numeric(as.matrix(allspecies_sub$cluster_ratio)))
  if (nrow(allspecies_sub2)>=5)
    avg_gene_length = sum(allspecies_sub2$gene_length)/nrow(allspecies_sub2)
  else
    avg_gene_length = sum(allspecies_geneonly$gene_length)/nrow(allspecies_geneonly)
  for(i in 1:simround)
  {
    # species simulation
    temp_species_HS = simulate_HS(as.integer(total_gene),as.integer(total_N),
                                  as.integer(total_S),as.integer(avg_gene_length),SNP_cutoff_species)
    temp_lineage_count = species_HS(temp_species_HS,SNP_cutoff_species,allspecies_sub)
    # lineage simulation
    for(rownum in 1:nrow(allspecies_sub))
    {
      speciesname = toString(allspecies_sub$X.donor_species[rownum])
      N_average = allspecies_sub$N[rownum]
      S_average = allspecies_sub$S[rownum]
      avg_gene =allspecies_sub$tota_gene_num[rownum]
      HS_NO = allspecies_sub$HS_NO[rownum]
      print(speciesname)
      print(SNP_cutoff_species)
      print(SNP_cutoff)
      temp_lineage = c(0,0,0,0)
      if(speciesname %in% temp_lineage_count$X.donor_species){
        temp_lineage=as.matrix(temp_lineage_count[
          which(temp_lineage_count$X.donor_species==speciesname),c(2:5)])
      }
      {temp_simulate[k,]=c(speciesname,i,
                           temp_lineage,
                           simulate_HS(as.integer(avg_gene),as.integer(N_average),
                                       as.integer(S_average),as.integer(avg_gene_length),SNP_cutoff),
                           as.integer(total_gene),as.integer(total_N),
                           as.integer(total_S),as.integer(avg_gene_length),
                           as.integer(HS_NO),avg_gene,N_average,S_average)
        k = k + 1}
    }
  }
}

colnames(temp_simulate)=c('X.donor_species','simulate_time','adaptive_genes_species_SNPs','adaptive_genes_species',
                          'adaptive_genes_species_N','adaptive_genes_species_S',
                          'adaptive_genes_lineage_SNPs','adaptive_genes_lineage',
                          'adaptive_genes_lineage_N','adaptive_genes_lineage_S',
                          'total_gene_species','N_species','S_species','avg_gene_length','HS_NO','avg_gene','N','S')
temp_simulate=data.frame(temp_simulate)
write.table(temp_simulate,'final_SNP/HSgene/HS_simulation.total.afterturningcutoff.clusterratio.withSNPs.lineage.txt',
            quote=F,sep='\t',row.names=F)
temp_simulate = read.delim('final_SNP/HSgene/HS_simulation.total.afterturningcutoff.clusterratio.withSNPs.lineage.txt',header=T)
temp_simulate2 = data.frame(
  X.donor_species = unique(temp_simulate$X.donor_species),
  Mean_SNP_lineage = 0,
  Mean_gene_lineage = 0,
  SNP_lineage = 0,
  HS_NO = 0
)
for(i in 1:nrow(temp_simulate2)){
  donor_species=temp_simulate2$X.donor_species[i]
  temp_simulate_sub = temp_simulate[which(temp_simulate$X.donor_species == donor_species),]
  temp_simulate2$Mean_SNP_lineage[i] = mean(temp_simulate_sub$adaptive_genes_lineage_SNPs)
  temp_simulate2$Mean_gene_lineage[i] = mean(temp_simulate_sub$adaptive_genes_lineage)
  temp_simulate2$SNP_lineage[i] = mean(temp_simulate_sub$N + temp_simulate_sub$S)
  temp_simulate2$HS_NO[i] = mean(temp_simulate_sub$HS_NO)
}
temp_simulate2$lower[which(temp_simulate2$lower<0)]=0
temp_simulate2$HS_NO[is.na(temp_simulate2$HS_NO)]=0
temp_simulate3 =temp_simulate2[which(temp_simulate2$HS_NO>0
                                     & temp_simulate2$Mean_SNP_lineage == 0),]
#
#---- check data----
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(stringr)

# plot
allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.txt',header=T)
allspecies2=read.delim('final_SNP/HSgene/all.species.High_select2.all.gene.txt',header=T)

MG_abu = FALSE
if (MG_abu){
  # add abundance
  Total_MG = 103
  allotu = read.delim('final_SNP/MG/MG.abu.sum.txt',header=T)
  allotu=allotu[,-ncol(allotu)]
  allotu$avg_depth=as.numeric(as.matrix(allotu$avg_depth))
  library(stringr)
  allotu$donor = stringi::stri_sub(allotu$sample, 1,2)
  allotu$lineage_donor = paste(allotu$donor_species,allotu$donor,sep='.donor.')
  MG_size = read.delim('final_SNP/MG/MG.size.txt',header=F)
  allotu=merge(allotu,MG_size,by.x='sample',by.y='V1')
  allotu$avg_depth[which(allotu$avg_depth==0)]=0.5
  allotu$avg_depth=allotu$avg_depth/allotu$V2
  allotu2=allotu[which(allotu$lineage_donor %in% as.matrix(allspecies$X.donor_species)),]
  allspecies=merge(allspecies,allotu2[,c(3,7)],by.x='X.donor_species',
                   by.y='lineage_donor',all.x=T)
  allspecies$avg_depth=log10(allspecies$avg_depth)
  allspecies2=allspecies[which(!is.na(allspecies$avg_depth)),]
  allspecies2=allspecies2[which(!is.na(allspecies2$dNdS_cut_highselect)),]
  
  }

# de novo mutations
ggscatter(allspecies, x = "genomelog10", y = "No.SNPlog10",
          color = 'short_residence',
          # add = "reg.line", conf.int = TRUE,
          palette = brewer.pal(5, "PRGn")[c(5,1)],
          size = 3,
          #label = 'species_short',
          title = 'de novo mutations',
          xlab = 'No. genomes per population (log10)',
          ylab = 'No. SNPs per population (log10)',
          add.params = list(color = brewer.pal(5, "PRGn")[c(5)]
                            , fill =brewer.pal(5, "PRGn")[c(4)]),
          cor.coef = TRUE, cor.method = "spearman")
ggscatter(allspecies, x = "genomelog10", y = "No.HS_SNP_NOlog10",
          color = 'short_residence',
          # add = "reg.line", conf.int = TRUE,
          palette = brewer.pal(5, "PRGn")[c(5,1)],
          size = 3,
          #label = 'species_short',
          title = 'de novo mutations',
          xlab = 'No. genomes per population (log10)',
          ylab = 'No. SNPs under parallel evolution (log10)',
          add.params = list(color = brewer.pal(5, "PRGn")[c(5)]
                            , fill =brewer.pal(5, "PRGn")[c(4)]),
          cor.coef = TRUE, cor.method = "spearman")


# species sum
allspecies_pop = count(allspecies,c('species'))
colnames(allspecies_pop)[2]='num_pop'
allspecies_pop$No.SNP = 0
allspecies_pop$HS_NO = 0
allspecies_pop$HS_SNP_NO = 0
allspecies_pop$No.genome = 0
for(i in 1:nrow(allspecies))
{
  j = which(allspecies_pop$species ==allspecies$species[i])
  allspecies_pop$No.SNP[j]=allspecies_pop$No.SNP[j] + allspecies$No.SNP[i]
  allspecies_pop$HS_NO[j]=allspecies_pop$HS_NO[j] + allspecies$HS_NO[i]
  allspecies_pop$HS_SNP_NO[j]=allspecies_pop$HS_SNP_NO[j] + allspecies$HS_SNP_NO[i]
  allspecies_pop$No.genome[j]=allspecies_pop$No.genome[j] + allspecies$No.genome[i]
  
}
allspecies_pop$No.genome[1]=44
allspecies_pop$num_poplog10 = log10(allspecies_pop$num_pop)
allspecies_pop=merge(allspecies_pop,allspecies[!duplicated(allspecies$species),c(2,17)],
                     by='species',all.x=T)
allspecies_pop$No.SNPlog10 = log10(allspecies_pop$No.SNP)
allspecies_pop$HS_NOlog10 = log10(allspecies_pop$HS_NO)
allspecies_pop$HS_SNP_NOlog10 = log10(allspecies_pop$HS_SNP_NO)
allspecies2_count=count(allspecies2,'species')
colnames(allspecies2_count)[2]='No.gene_SNPs'
allspecies_pop=merge(allspecies_pop,allspecies2_count,
                     by='species')
write.table(allspecies_pop,'final_SNP/HSgene/all.species.all.species.lineagenum.txt',
            quote=F,sep='\t',row.names=F)

# No. SNPs versus NO. HS, HS ratio per species line
allspecies2 = allspecies[which(allspecies$species %in% 
                                 allspecies_pop$species[which(allspecies_pop$num_pop > 2)]),]
allspecies2 = allspecies2[which(allspecies2$species %in% 
                                 unique(allspecies_pop$species[which(allspecies_pop$HS_NO > 1)])),]
allspecies2$logSNP = log10(as.numeric(as.matrix(allspecies2$No.SNP))+1)
allspecies2$logHS_SNP = log10(as.numeric(as.matrix(allspecies2$HS_SNP_NO))+1)
allspecies2$logHS = log10(as.numeric(as.matrix(allspecies2$HS_NO))+1)
ggplot(allspecies2,aes(y=logHS,
                      x=logSNP,
                      group = species))+
  geom_point(aes(
    color=species,
    fill=species
  ),
  alpha = 0.6,
  size=2
  )+ 
  geom_smooth(
    aes(
      group = species,
      color=species,
      fill=species
    ),
    se = FALSE
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =c(brewer.pal(9, "Set1")[c(1:9)],brewer.pal(12, "Set3")[c(1:12)]))+
  scale_color_manual(values =c(brewer.pal(9, "Set1")[c(1:9)],brewer.pal(12, "Set3")[c(1:12)]))

# No. genome versus No.SNPs, species color
allspecies2 = allspecies[which(allspecies$species %in% 
                                 allspecies_pop$species[which(allspecies_pop$num_pop > 2)]),]
allspecies2$logSNP = log10(as.numeric(as.matrix(allspecies2$No.SNP))+1)
allspecies2$logHS_SNP = log10(as.numeric(as.matrix(allspecies2$HS_SNP_NO))+1)
allspecies2$logHS = log10(as.numeric(as.matrix(allspecies2$HS_NO))+1)
allspecies2$No.genomelog10 = log10(as.numeric(as.matrix(allspecies2$No.genome))+1)

ggplot(allspecies2,aes(y=logSNP,
                       x=No.genomelog10,
                       group = species))+
  geom_point(aes(
    color=species,
    fill=species
  ),
  alpha = 0.6,
  size=2
  )+ 
  geom_smooth(
    aes(
      group = species,
      color=species,
      fill=species
    ),
    se = FALSE
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =c(brewer.pal(9, "Set1")[c(1:9)],brewer.pal(12, "Set3")[c(1:12)]))+
  scale_color_manual(values =c(brewer.pal(9, "Set1")[c(1:9)],brewer.pal(12, "Set3")[c(1:12)]))

# donor sum
allspecies_pop = count(allspecies,c('donor'))
colnames(allspecies_pop)[2]='num_pop'
allspecies_pop$No.SNP = 0
allspecies_pop$HS_NO = 0
allspecies_pop$HS_SNP_NO = 0
allspecies_pop$No.genome = 0
for(i in 1:nrow(allspecies))
{
  j = which(allspecies_pop$donor ==allspecies$donor[i])
  allspecies_pop$No.SNP[j]=allspecies_pop$No.SNP[j] + allspecies$No.SNP[i]
  allspecies_pop$HS_NO[j]=allspecies_pop$HS_NO[j] + allspecies$HS_NO[i]
  allspecies_pop$HS_SNP_NO[j]=allspecies_pop$HS_SNP_NO[j] + allspecies$HS_SNP_NO[i]
  allspecies_pop$No.genome[j]=allspecies_pop$No.genome[j] + allspecies$No.genome[i]
  
}

allspecies_pop$num_poplog10 = log10(allspecies_pop$num_pop)
allspecies_pop$No.SNPlog10 = log10(allspecies_pop$No.SNP)

allspecies_pop=allspecies_pop[which(allspecies_pop$donor!=''),]
write.table(allspecies_pop,'final_SNP/HSgene/all.species.all.donor.lineagenum.txt',
            quote=F,sep='\t',row.names=F)
# No. SNPs versus NO. HS, HS ratio per donor line
allspecies2 = allspecies[which(allspecies$donor %in% 
                                 allspecies_pop$donor[which(allspecies_pop$num_pop > 2)]),]
allspecies2 = allspecies2[which(allspecies2$donor %in% 
                                  unique(allspecies_pop$donor[which(allspecies_pop$HS_NO > 1)])),]
allspecies2$logSNP = log10(as.numeric(as.matrix(allspecies2$No.SNP))+1)
allspecies2$logHS_SNP = log10(as.numeric(as.matrix(allspecies2$HS_SNP_NO))+1)
allspecies2$logHS = log10(as.numeric(as.matrix(allspecies2$HS_NO))+1)
ggplot(allspecies2,aes(y=logHS_SNP,
                       x=logSNP,
                       group = donor))+
  geom_point(aes(
    color=donor,
    fill=donor
  ),
  alpha = 0.6,
  size=2
  )+ 
  geom_smooth(
    aes(
      group = donor,
      color=donor,
      fill=donor
    ),
    se = FALSE
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =c(brewer.pal(9, "Set1")[c(1:9)],brewer.pal(12, "Set3")[c(1:12)]))+
  scale_color_manual(values =c(brewer.pal(9, "Set1")[c(1:9)],brewer.pal(12, "Set3")[c(1:12)]))

#---- general figure 2----
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(stringr)

# plot
allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.txt',header=T)
allspecies2=read.delim('final_SNP/HSgene/all.species.High_select2.all.gene.txt',header=T)

# number of highly selected lineage
{temp_simulate_sum3=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted2.txt',header=T)
temp_simulate_sum3 = temp_simulate_sum3[!duplicated(temp_simulate_sum3$lineage),]
neworder=read.delim('species_order.new.txt',header=T)
temp_simulate_sum3=merge(temp_simulate_sum3,neworder,
                         by='species',all=T)
temp_simulate_sum3_sub = temp_simulate_sum3[which(!is.na(temp_simulate_sum3$HS_norm_gene) &
                                                    temp_simulate_sum3$HS_norm_gene >0),]

temp_simulate_sum3=rbind(
  temp_simulate_sum3_sub,
  temp_simulate_sum3[which(!(temp_simulate_sum3$species %in% temp_simulate_sum3_sub$species)),]
)
qualify_species = unique(temp_simulate_sum3$species[
  which(!is.na(temp_simulate_sum3$HS_norm_gene))
  ])
temp_simulate = read.delim('final_SNP/HSgene/HS_simulation.total.afterturningcutoff.clusterratio.withSNPs.txt',
                           header=T)
multiple_lineage = unique(temp_simulate_sum3$species[which(temp_simulate_sum3$total_lineage_num>=6)])
temp_simulate = temp_simulate[which(temp_simulate$species
                                    %in% qualify_species),]
temp_simulate$adaptive_genes_species[which(!(temp_simulate$species
                                           %in% multiple_lineage))]=0
temp_simulate$adaptive_genes = temp_simulate$adaptive_genes_lineage
temp_simulate2 = data.frame(
  species = unique(temp_simulate$species),
  lower = 0,
  higher = 0,
  Mean = 0,
  Mean_SNP = 0
)
for(i in 1:nrow(temp_simulate2)){
  species=temp_simulate2$species[i]
  temp_simulate_sub = temp_simulate[which(temp_simulate$species == species),]
  temp_simulate2$lower[i]=quantile(temp_simulate_sub$adaptive_genes,c(0.05))
  temp_simulate2$higher[i]=quantile(temp_simulate_sub$adaptive_genes,c(0.95))
  temp_simulate2$Mean[i]=quantile(temp_simulate_sub$adaptive_genes,c(0.5))
  temp_simulate2$Mean_SNP[i]=quantile(temp_simulate_sub$adaptive_genes_lineage_SNPs,
                                      c(0.5))
 }

temp_simulate2$lower[which(temp_simulate2$lower<0)]=0

temp_simulate_sum3 = merge(temp_simulate_sum3,temp_simulate2,
                       by='species',all.x=T)

temp_simulate_sum3$HS_gene[is.na(temp_simulate_sum3$HS_gene)]=0

#Plot
ggplot()+
  geom_boxplot(data = temp_simulate_sum3,
               aes(
                 y=as.numeric(as.matrix(HS_lineage)),
                 x=neworder,
                 group = species_short,
                 color='1',
                 fill='1'
               ),
               outlier.shape = NA,
               alpha = 0.5,
               size=0.8
  )+ stat_summary(data = temp_simulate_sum3,fun.y = mean, geom = "point", 
                  aes( y=as.numeric(as.matrix(HS_lineage)),
                       x=neworder,
                       group = species_short,
                       color='1',fill = '1'),alpha = 1,
                  size=3)+
  geom_point(data = temp_simulate_sum3,
             aes(
               x=neworder,
               y=as.numeric(as.matrix(Mean)),
               color='2',
               fill='2'
             ),
             size=2,
             alpha = 1
  )+ 
  geom_errorbar(data = temp_simulate_sum3,
                aes(x=neworder,
                    ymin=as.numeric(as.matrix(lower)), 
                    ymax=as.numeric(as.matrix(higher)),
                    color='2'), 
                width=0.6,
                position=position_dodge(.9))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(values =c(brewer.pal(5, "PRGn")[c(1)],brewer.pal(5, "PuOr")[c(1)]))+
  scale_color_manual(values =c(brewer.pal(5, "PRGn")[c(1)],brewer.pal(5, "PuOr")[c(1)]))+
  scale_x_continuous(limits = c(0,37),breaks = c(1,36),labels = c(1,36))+
  scale_y_continuous(breaks = c(0,1,2,3,4,5,10,15,20))
}
# HS norm lineage
ggplot(temp_simulate_sum3,aes(y=HS_norm_gene_lineage,
                      x=neworder,
                      group = species_short))+
  geom_boxplot(aes(
    color='1',
    fill='1'
  ),
  outlier.shape = NA,
  alpha = 0.5,
  size=0.8
  )+ stat_summary(fun.y = mean, geom = "point", 
                  aes(color='1',fill = '1'),alpha = 1,
                  size=3)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1)])+
  scale_x_continuous(limits = c(0,37),breaks = c(1,36),labels = c(1,36))+
  geom_hline(yintercept = 0.6)
# HS norm ratio lineage
ggplot(temp_simulate_sum3,aes(y=HSratio_gene_diff_lineage,
                              x=neworder,
                              group = species_short))+
  geom_boxplot(aes(
    color='1',
    fill='1'
  ),
  outlier.shape = NA,
  alpha = 0.5,
  size=0.8
  )+ stat_summary(fun.y = mean, geom = "point", 
                  aes(color='1',fill = '1'),alpha = 1,
                  size=3)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1)])+
  scale_x_continuous(limits = c(0,37),breaks = c(1,36),labels = c(1,36))

ggplot(temp_simulate_sum3,aes(y=HSratio_gene_diff_lineage,
                              x=neworder,
                              group = species_short))+
  geom_boxplot(aes(
    color='1',
    fill='1'
  ),
  outlier.shape = NA,
  alpha = 0.5,
  size=0.8
  )+ stat_summary(fun.y = mean, geom = "point", 
                  aes(color='1',fill = '1'),alpha = 1,
                  size=3)+
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               aes(color='1',
                   width = 0.8,
                   linetype = "solid"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1)])+
  scale_x_continuous(limits = c(0,37),breaks = c(1,36),labels = c(1,36))+
  scale_y_continuous(limits = c(-0.05,0.2))+
geom_hline(yintercept = 0.01)

# dmrca
ggplot(allspecies,aes(y=dmrca,
                      x=neworder,
                      group = species_short))+
  geom_boxplot(aes(
    color='1',
    fill='1'
  ),
  outlier.shape = NA,
  alpha = 0.5,
  size=0.8
  )+ stat_summary(fun.y = mean, geom = "point", 
                  aes(color='1',fill = '1'),alpha = 1,
                  size=3)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(5)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(5)])+
  scale_x_continuous(limits = c(0,37),breaks = c(1,36),labels = c(1,36))

ggplot(allspecies,aes(y=dmrca,
                      x=neworder,
                      group = species_short))+
  geom_boxplot(aes(
    color='1',
    fill='1'
  ),
  outlier.shape = NA,
  alpha = 0.5,
  size=0.8
  )+ stat_summary(fun.y = mean, geom = "point", 
                  aes(color='1',fill = '1'),alpha = 1,
                  size=3)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(5)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(5)])+
  scale_x_continuous(limits = c(0,37),breaks = c(1,36),labels = c(1,36))+
  scale_y_continuous(limits = c(0,20))
# number of highly selected species
{temp_simulate_sum3=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted2.txt',header=T)
  temp_simulate_sum3 = temp_simulate_sum3[!duplicated(temp_simulate_sum3$lineage),]
  neworder=read.delim('species_order.new.txt',header=T)
  temp_simulate_sum3=merge(temp_simulate_sum3,neworder,
                           by='species',all=T)
  temp_simulate_sum3_sub = temp_simulate_sum3[which(!is.na(temp_simulate_sum3$HS_norm_gene) &
                                                      temp_simulate_sum3$HS_norm_gene >0),]
  
  temp_simulate_sum3=rbind(
    temp_simulate_sum3_sub,
    temp_simulate_sum3[which(!(temp_simulate_sum3$species %in% temp_simulate_sum3_sub$species)),]
  )
  qualify_species = unique(temp_simulate_sum3$species[
    which(!is.na(temp_simulate_sum3$HS_norm_gene))
    ])
  temp_simulate = read.delim('final_SNP/HSgene/HS_simulation.total.afterturningcutoff.clusterratio.withSNPs.txt',
                             header=T)
  multiple_lineage = unique(temp_simulate_sum3$species[which(temp_simulate_sum3$total_lineage_num>5)])
  temp_simulate = temp_simulate[which(temp_simulate$species
                                      %in% qualify_species),]
  temp_simulate$adaptive_genes_species[which(!(temp_simulate$species
                                               %in% multiple_lineage))]=0
  temp_simulate$adaptive_genes = temp_simulate$adaptive_genes_species
  temp_simulate2 = data.frame(
    species = unique(temp_simulate$species),
    lower = 0,
    higher = 0,
    Mean = 0,
    Mean_SNP = 0
  )
  for(i in 1:nrow(temp_simulate2)){
    species=temp_simulate2$species[i]
    temp_simulate_sub = temp_simulate[which(temp_simulate$species == species),]
    temp_simulate2$lower[i]=quantile(temp_simulate_sub$adaptive_genes,c(0.05))
    temp_simulate2$higher[i]=quantile(temp_simulate_sub$adaptive_genes,c(0.95))
    temp_simulate2$Mean[i]=quantile(temp_simulate_sub$adaptive_genes,c(0.5))
    temp_simulate2$Mean_SNP[i]=quantile(temp_simulate_sub$adaptive_genes_lineage_SNPs,
                                        c(0.5))
  }
  
  temp_simulate2$lower[which(temp_simulate2$lower<0)]=0
  
  temp_simulate_sum3 = merge(temp_simulate_sum3,temp_simulate2,
                             by='species',all.x=T)
  
  temp_simulate_sum3$HS_gene[is.na(temp_simulate_sum3$HS_gene)]=0
  temp_simulate_sum3 = temp_simulate_sum3[which(temp_simulate_sum3$species
                                                %in% multiple_lineage),]
  #Plot
  ggplot()+
    geom_boxplot(data = temp_simulate_sum3,
                 aes(
                   y=as.numeric(as.matrix(HS_species)),
                   x=species,
                   group = species_short,
                   color='1',
                   fill='1'
                 ),
                 outlier.shape = NA,
                 alpha = 0.5,
                 size=0.8
    )+ stat_summary(data = temp_simulate_sum3,fun.y = mean, geom = "point", 
                    aes( y=as.numeric(as.matrix(HS_species)),
                         x=species,
                         group = species_short,
                         color='1',fill = '1'),alpha = 1,
                    size=3)+
    geom_point(data = temp_simulate_sum3,
               aes(
                 x=species,
                 y=as.numeric(as.matrix(Mean)),
                 color='2',
                 fill='2'
               ),
               size=2,
               alpha = 1
    )+ 
    geom_errorbar(data = temp_simulate_sum3,
                  aes(x=species,
                      ymin=as.numeric(as.matrix(lower)), 
                      ymax=as.numeric(as.matrix(higher)),
                      color='2'), 
                  width=0.6,
                  position=position_dodge(.9))+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(values =c(brewer.pal(5, "PRGn")[c(1)],brewer.pal(5, "PuOr")[c(1)]))+
    scale_color_manual(values =c(brewer.pal(5, "PRGn")[c(1)],brewer.pal(5, "PuOr")[c(1)]))+
    scale_y_continuous(breaks = c(0,1,2,3,4,5,10,15,20))
  }
# HS norm species
ggplot(temp_simulate_sum3,aes(y=HS_norm_gene_species,
                              x=species,
                              group = species_short))+
  geom_boxplot(aes(
    color='1',
    fill='1'
  ),
  outlier.shape = NA,
  alpha = 0.5,
  size=0.8
  )+ stat_summary(fun.y = mean, geom = "point", 
                  aes(color='1',fill = '1'),alpha = 1,
                  size=3)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1)])
# HS norm ratio species
ggplot(temp_simulate_sum3,aes(y=HSratio_gene_diff_species,
                              x=species,
                              group = species_short))+
  geom_boxplot(aes(
    color='1',
    fill='1'
  ),
  outlier.shape = NA,
  alpha = 0.5,
  size=0.8
  )+ stat_summary(fun.y = mean, geom = "point", 
                  aes(color='1',fill = '1'),alpha = 1,
                  size=3)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1)])

#
#----adaptive SNPs ratio vesus abundance species----
# read abundance
#allotu = read.delim('final_SNP/MG/all.otu_table.metagenome.txt',header=T)
#unique(allotu$STArea)
#allotu=allotu[which(allotu$STArea == 'Gut'),]
#write.table(allotu,'final_SNP/MG/all.otu_table.gut.metagenome.txt',
#            quote=F,sep='\t',row.names=F)

# adaptive SNPs ratio!!! adaptive genes ratio!!!
# exclude species that are not significant!!!
# abundance of 16S + ubiquity!!!
allotu = read.delim('final_SNP/MG/all.otu_table.gut.metagenome.txt',header=T)
allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.species.txt',header=T)
allspecies=allspecies[which(allspecies$No.SNP>0),-c(15:26)]
library(stringr)
library(psych)
library(ggpubr)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
allspecies_set = unique(allspecies$species_short)
temp_abu = allotu[,c(10:ncol(allotu))]
temp_abu=as.matrix(temp_abu)
temp_abu = temp_abu[which(temp_abu>0)]
temp_abu=log10(temp_abu)
hist(temp_abu)
quantile(temp_abu,c(0.001,0.01,0.05,0.1)) #-6.154902 -5.744727 -5.130768 -4.772113 
cutoff_prevalence = 10^quantile(temp_abu,c(0.001))
allotunew = data.frame(species_short = allspecies_set,
                       abu_mean = 0,
                       abu_median = 0,
                       abu_harmonic_mean = 0,
                       abu_CV = 0,
                       prevalence = 0)
for(i in 1:nrow(allotunew))
{
  species = allotunew$species_short[i]
  genus = str_split_fixed(species, "_", 2)[,1]
  species_sub = str_split_fixed(species, "_", 2)[,2]
  col_species = which(grepl(genus, colnames(allotu)) & grepl(species_sub, colnames(allotu)))
  colnames(allotu)[col_species]
  Set = FALSE
  if (length(col_species)>0)
  {abu_temp = allotu[,col_species[1]]
  Set = TRUE}else
  {
    col_species = which(grepl(genus, colnames(allotu)))
    if (length(col_species)>0)
      {abu_temp = allotu[,col_species[1]]
      Set = TRUE}
  }
  if(Set)
  {allotunew$abu_mean[i]=mean(abu_temp[which(abu_temp>0)])
  allotunew$abu_median[i]=median(abu_temp[which(abu_temp>0)])
  allotunew$abu_harmonic_mean[i]=harmonic.mean(abu_temp[which(abu_temp>0)],zero = F)
  allotunew$abu_CV[i]=sd(abu_temp[which(abu_temp>0)])/mean(abu_temp[which(abu_temp>0)])
  allotunew$prevalence[i]=length(abu_temp[which(abu_temp>cutoff_prevalence)])/length(abu_temp)
  }else
    print(species)
}
allotunew$logabu = log10(allotunew$abu_mean)
allotunew=allotunew[which(allotunew$species_short!='Ruthen_lactat'),]
library(ggpubr)
ggscatter(allotunew, x = "logabu", y = "prevalence",
          add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(1)],
          size = 3,
          label = 'species_short',
          title = 'HMP',
          xlab = 'abundance (log10)',
          ylab = 'prevalence (%)',
          add.params = list(color = brewer.pal(5, "PRGn")[c(1)]
                            , fill =brewer.pal(5, "PRGn")[c(2)]),
          cor.coef = TRUE, cor.method = "pearson")

allspecies=merge(allspecies,allotunew,by='species_short',all=T)
allspecies=allspecies[which(allspecies$species_short!='Ruthen_lactat'),]
temp_species=allspecies[!duplicated(allspecies$species),c(1,2)]
allotunew=merge(allotunew,
                temp_species,
                by='species_short')
allspecies2=read.delim('final_SNP/HSgene/all.species.all.species.lineagenum.txt',header=T)
allspecies2=merge(allspecies2,allotunew,
                 by='species',all=T)
# read abundance of bn10 metagenomes
mgabu = FALSE
if(mgabu)
{# read abundance of metagenomes
Total_MG = 103
allotu = read.delim('final_SNP/MG/MG.abu.sum.txt',header=T)
allotu=allotu[,-ncol(allotu)]
allotu$avg_depth=as.numeric(as.matrix(allotu$avg_depth))
library(stringr)
allotu$donor = stringi::stri_sub(allotu$sample, 1,2)
allotu$donor_species = paste(allotu$donor_species,allotu$donor,sep = '.donor.')
MG_size = read.delim('final_SNP/MG/MG.size.txt',header=F)
allotu=merge(allotu,MG_size,by.x='sample',by.y='V1')
allotu$avg_depth=allotu$avg_depth/allotu$V2
allotu2=allotu
allotu2$lineage = str_split_fixed(allotu2$donor_species, ".donor", 2)[,1]
allspecies_set = unique(allspecies$species)
allotunew = data.frame(species = allspecies_set,
                       abu_mean_MG = 0,
                       abu_median_MG = 0,
                       abu_CV_MG = 0,
                       prevalence_MG = 0)
for(i in 1:nrow(allotunew))
{
  species = allotunew$species[i]
  abu_temp=allotu2[which(grepl(species, allotu2$donor_species)),]
  temp_line = matrix(0,ncol = ncol(abu_temp),
                     nrow = Total_MG*length(unique(abu_temp$lineage)) - nrow(abu_temp))
  colnames(temp_line)=colnames(abu_temp)
  abu_temp = rbind(as.matrix(abu_temp), temp_line)[,3]
  abu_temp=as.numeric(abu_temp)
  allotunew$abu_mean_MG[i]=mean(abu_temp[which(abu_temp>0)])
  allotunew$abu_median_MG[i]=median(abu_temp[which(abu_temp>0)])
  allotunew$abu_CV_MG[i]=sd(abu_temp[which(abu_temp>0)])/mean(abu_temp[which(abu_temp>0)])
  abu_temp=allotu2[which(grepl(species, allotu2$donor_species)),]
  temp_line = matrix(0,ncol = ncol(abu_temp),
                     nrow = Total_MG*length(unique(abu_temp$lineage)) - nrow(abu_temp))
  colnames(temp_line)=colnames(abu_temp)
  abu_temp = rbind(as.matrix(abu_temp), temp_line)[,3]
  abu_temp=as.numeric(abu_temp)
  allotunew$prevalence_MG[i]=length(abu_temp[which(abu_temp>0)])/length(abu_temp)
  
}
allotunew=allotunew[which(!is.nan(allotunew$abu_mean)),]
allspecies2=merge(allspecies2,allotunew,
                  by='species',all=T)}
# read abundance of bn10 16S
bn10abu = T
if(bn10abu)
{# read abundance of bn10
  allotu2 = read.delim('final_SNP/MG/BN10_16S_bacterial_taxa.tab',header=T)
  allotu2_colsum = colSums(allotu2)
  y = quantile(allotu2_colsum,c(0.001,0.001,0.01,0.05,0.1,0.2,0.5,1))
  x = c(0.001,0.001,0.01,0.05,0.1,0.2,0.5,1)
  plot(x,log10(y))
  sample_cutoff = y[3]
  allotu2=allotu2[,which(allotu2_colsum>=sample_cutoff)]# 1 inside
  allotu2_colsum=allotu2_colsum[which(allotu2_colsum>=sample_cutoff)]
  for(i in 2:ncol(allotu2))
  {
    allotu2[,i]=allotu2[,i]/allotu2_colsum[i]
  }
  allotu2_taxa = read.delim('final_SNP/MG/BN10_16S_taxid_to_taxonomy.tab',header=T)
  allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.species.txt',header=T)
  allspecies=allspecies[which(allspecies$No.SNP>0),-c(15:26)]
  
  library(stringr)
  library(psych)
  library(ggpubr)
  allspecies_set = unique(allspecies$species_short)
  temp_abu = allotu2[,c(2:ncol(allotu2))]
  temp_abu=as.matrix(temp_abu)
  temp_abu = temp_abu[which(temp_abu>0)]
  temp_abu=log10(temp_abu)
  hist(temp_abu)
  quantile(temp_abu,c(0.001,0.01,0.05,0.1)) #-6.810950 -6.607194 -6.392091 -6.215239
  cutoff_prevalence = 10^quantile(temp_abu,c(0.001))
  allotunew = data.frame(species_short = allspecies_set,
                         abu_mean = 0,
                         abu_median = 0,
                         abu_harmonic_mean = 0,
                         abu_CV = 0,
                         prevalence = 0)
  for(i in 1:nrow(allotunew))
  {
    species = allotunew$species_short[i]
    genus = str_split_fixed(species, "_", 2)[,1]
    species_sub = str_split_fixed(species, "_", 2)[,2]
    col_species = which(grepl(genus, allotu2_taxa$taxonomy) & 
                          grepl(species_sub, allotu2_taxa$taxonomy))
    Set = FALSE
    if (length(col_species)>0)
    {abu_temp = allotu2[which(allotu2$OTU.ID %in% allotu2_taxa$X.OTU.ID[col_species]),c(2:ncol(allotu2))]
    Set = TRUE}else
    {
      col_species = which(grepl(genus, allotu2_taxa$taxonomy))
      if (length(col_species)>0)
      {abu_temp = allotu2[which(allotu2$OTU.ID %in% allotu2_taxa$X.OTU.ID[col_species]),c(2:ncol(allotu2))]
      Set = TRUE}
    }
    allotu2_taxa$taxonomy[col_species]
    if(Set)
    {
      abu_temp=as.matrix(abu_temp)
      allotunew$abu_mean[i]=mean(abu_temp[which(abu_temp>0)])
    allotunew$abu_median[i]=median(abu_temp[which(abu_temp>0)])
    allotunew$abu_harmonic_mean[i]=harmonic.mean(abu_temp[which(abu_temp>0)],zero = F)
    allotunew$abu_CV[i]=sd(abu_temp[which(abu_temp>0)])/mean(abu_temp[which(abu_temp>0)])
    allotunew$prevalence[i]=length(abu_temp[which(abu_temp>cutoff_prevalence)])/length(abu_temp)
    }else
      print(species)
  }
  allotunew$logabu = log10(allotunew$abu_mean)
  allspecies=merge(allspecies,allotunew,by='species_short',all=T)
  temp_species=allspecies[!duplicated(allspecies$species),c(1,2)]
  allotunew=merge(allotunew,
                  temp_species,
                  by='species_short')
  library(ggpubr)
  ggscatter(allotunew, x = "logabu", y = "prevalence",
            add = "reg.line", conf.int = TRUE,
            color = brewer.pal(5, "PRGn")[c(1)],
            size = 3,
            label = 'species_short',
            title = 'HMP',
            xlab = 'abundance (log10)',
            ylab = 'prevalence (%)',
            add.params = list(color = brewer.pal(5, "PRGn")[c(1)]
                              , fill =brewer.pal(5, "PRGn")[c(2)]),
            cor.coef = TRUE, cor.method = "pearson")
    allspecies2=merge(allspecies2,allotunew,
                    by='species',all=T)
    }
allspecies2=allspecies2[rev(order(allspecies2$abu_mean.y)),]
allspecies2=allspecies2[!duplicated(allspecies2$species),]
write.table(allspecies2,'final_SNP/HSgene/all.species.all.species.lineagenum.abu.txt',
            quote=F,sep='\t',row.names=F)
# add dnds and dmrca
allspecies=read.delim('final_SNP/HSgene/all.species.all.species.lineagenum.abu.txt',header=T)
allspecies$HSraio = allspecies$HS_SNP_NO/allspecies$No.SNP
allspecies$HSraio_gene = allspecies$HS_NO/allspecies$No.gene_SNPs
allspecies2 =read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.species.txt',header=T)
allspecies2$S_cut = as.numeric(as.matrix(allspecies2$S))
allspecies2$S_cut[which(allspecies2$S_cut==0)]=0.5
allspecies2$dNdS_cut = 0
allspecies2$dNdS_cut=allspecies2$N/allspecies2$S_cut/as.numeric(as.matrix(allspecies2$expected_ratio))
allspecies=allspecies[order(allspecies$species),]
allspecies2=allspecies2[order(allspecies2$species),]
allspecies$dNdS_cut_all = c(allspecies2$dNdS_cut[which(!(grepl('select', allspecies2$gene)))])
allspecies$dNdS_cut_HS = c(allspecies2$dNdS_cut[which((grepl('select', allspecies2$gene)))])
allspecies$dmrca = c(allspecies2$dmrca[which(!(grepl('select', allspecies2$gene)))])

write.table(allspecies,'final_SNP/HSgene/all.species.all.species.lineagenum.abu.dnds.txt',
            quote=F,sep='\t',row.names=F)
# add simulation
allspecies=read.delim('final_SNP/HSgene/all.species.all.species.lineagenum.abu.dnds.txt',header=T)
temp_simulate = read.delim('final_SNP/HSgene/HS_simulation.total.afterturningcutoff.clusterratio.withSNPs.txt',
                           header=T)
temp_simulate$adaptive_genes = temp_simulate$adaptive_genes_species+temp_simulate$adaptive_genes_lineage
temp_simulate$adaptive_genes_SNPs = temp_simulate$adaptive_genes_species_SNPs+temp_simulate$adaptive_genes_lineage_SNPs
temp_simulate2 = data.frame(
  species = unique(temp_simulate$species),
  Mean_sim_species = 0,
  Mean_sim_species_SNPs = 0,
  Mean_sim_lineage = 0,
  Mean_sim_lineage_SNPs = 0,
  Mean_sim = 0,
  Mean_sim_SNPs = 0,
  total_gene = 0,
  avg_gene = 0,
  avg_gene_length = 0
)
for(i in 1:nrow(temp_simulate2)){
  species=temp_simulate2$species[i]
  temp_simulate_sub = temp_simulate[which(temp_simulate$species == species),]
  temp_simulate2$Mean_sim[i]=mean(temp_simulate_sub$adaptive_genes)
  temp_simulate2$Mean_sim_SNPs[i]=mean(temp_simulate_sub$adaptive_genes_SNPs)
  temp_simulate2$Mean_sim_species[i]=mean(temp_simulate_sub$adaptive_genes_species)
  temp_simulate2$Mean_sim_species_SNPs[i]=mean(temp_simulate_sub$adaptive_genes_species_SNPs)
  temp_simulate2$Mean_sim_lineage[i]=mean(temp_simulate_sub$adaptive_genes_lineage)
  temp_simulate2$Mean_sim_lineage_SNPs[i]=mean(temp_simulate_sub$adaptive_genes_lineage_SNPs)
  temp_simulate2$total_gene[i]=mean(temp_simulate_sub$total_gene)
  temp_simulate2$avg_gene[i]=mean(temp_simulate_sub$avg_gene)
  temp_simulate2$avg_gene_length[i]=mean(temp_simulate_sub$avg_gene_length)
}

allspecies=merge(allspecies,temp_simulate2,by='species',all.x=T)
write.table(allspecies,'final_SNP/HSgene/all.species.all.species.lineagenum.abu.dnds.sim.txt',
            quote=F,sep='\t',row.names=F)

# linear regression
allspecies=read.delim('final_SNP/HSgene/all.species.all.species.lineagenum.abu.dnds.sim.txt',header=T)
allspecies$dNdS_cut_HS[is.na(allspecies$dNdS_cut_HS)]=0
allspecies2=allspecies
Not_qualify = c()#c('BaVu','BaTh')#'ClBe','BaVu'
Out_lier = c('ClBe')
# plot add error bar!!!
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(plyr)
library(stringr)
library(rpart)
library(randomForest)
allspecies = allspecies2
allspecies$logabuy = log10(allspecies$abu_mean.y)
allspecies$species[which(allspecies$logabuy< -3)]
allspecies$logabu = log10(allspecies$abu_mean.x)
allspecies$species[which(allspecies$logabu< -3)]
allspecies=allspecies[which(!(allspecies$species %in% Not_qualify |
                                allspecies$species %in% Out_lier )),]
allspecies$HSraio_sim = allspecies$Mean_sim_SNPs/allspecies$No.SNP
allspecies$HSraio_sim2 = allspecies$Mean_sim/allspecies$No.gene_SNPs
allspecies$NO.SNP_lineage = allspecies$No.SNP/allspecies$num_pop
allspecies$simSNP_lineage = allspecies$Mean_sim_lineage_SNPs/allspecies$num_pop
allspecies$simSNP_species = allspecies$Mean_sim_species_SNPs/allspecies$num_pop
allspecies$simSNP_all = allspecies$Mean_sim_SNPs/allspecies$num_pop
allspecies$HS_SNP_NO_all = allspecies$HS_SNP_NO/allspecies$num_pop
allspecies$simgene_lineage = allspecies$Mean_sim_lineage/allspecies$num_pop
allspecies$simgene_species = allspecies$Mean_sim_species/allspecies$num_pop
allspecies$simgene_all = allspecies$Mean_sim/allspecies$num_pop
allspecies$HS_NO_all = allspecies$HS_NO/allspecies$num_pop
allspecies=allspecies[!is.na(allspecies$logabu),]

ggscatter(allspecies, x = "NO.SNP", y = "HS_SNP_NO_all",
          # add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(1)],
          size = 3,
          label = 'species_short.y',
          add.params = list(color = brewer.pal(5, "PRGn")[c(1)]
                            , fill =brewer.pal(5, "PRGn")[c(2)]),
          cor.coef = TRUE, cor.method = "spearman")

ggscatter(allspecies, x = "HSraio", y = "HSraio_gene",
          # add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(1)],
          size = 3,
          label = 'species_short.y',
          add.params = list(color = brewer.pal(5, "PRGn")[c(1)]
                            , fill =brewer.pal(5, "PRGn")[c(2)]),
          cor.coef = TRUE, cor.method = "spearman")

fit <- rpart(HSraio ~ logabu + NO.SNP_lineage + avg_gene + 
        avg_gene_length + total_gene + num_pop, data=allspecies,method = 'anova') 

printcp(fit) # display the results
plotcp(fit) # visualize cross-validation results
summary(fit) # detailed summary of splits
par(mfrow=c(1,2)) # two plots on one page
rsq.rpart(fit) # visualize cross-validation results  

# plot tree
plot(fit, uniform=TRUE)
text(fit, use.n=TRUE, all=TRUE, cex=.8)

# prune the tree
#pfit<- prune(fit, cp=fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]) # from cptable   
#plot(pfit, uniform=TRUE)
#text(pfit, use.n=TRUE, all=TRUE, cex=.8)

# random forests
library(randomForest)
fit <- randomForest(HS_SNP_NO_all ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
                      avg_gene_length + total_gene + num_pop, data=allspecies)
print(fit) # view results
importance(fit, scale = FALSE) # importance of each predictor
allspecies$HS

# linear regression
fit_sum = matrix(0,nrow=0,ncol=9)

fit <- lm(simSNP_lineage ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
            avg_gene_length + total_gene, data=allspecies)
summary_fit = summary(fit) # show results
fit_sum=rbind(fit_sum,
              c('simSNP_lineage',coefficients(fit),summary_fit$adj.r.squared))

fit <- lm(simSNP_species ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
            avg_gene_length + total_gene, data=allspecies)
summary_fit = summary(fit) # show results
fit_sum=rbind(fit_sum,
              c('simSNP_species',coefficients(fit),summary_fit$adj.r.squared))

fit <- lm(simSNP_all ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
            avg_gene_length + total_gene, data=allspecies)
summary_fit = summary(fit) # show results
fit_sum=rbind(fit_sum,
              c('simSNP_all',coefficients(fit),summary_fit$adj.r.squared))

fit <- lm(HS_SNP_NO_all ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
            avg_gene_length + total_gene, data=allspecies)
summary_fit = summary(fit) # show results
fit_sum=rbind(fit_sum,
              c('HS_SNP_NO_all',coefficients(fit),summary_fit$adj.r.squared))

fit <- lm(simgene_lineage ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
            avg_gene_length + total_gene, data=allspecies)
summary_fit = summary(fit) # show results
fit_sum=rbind(fit_sum,
              c('simgene_lineage',coefficients(fit),summary_fit$adj.r.squared))

fit <- lm(simgene_species ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
            avg_gene_length + total_gene, data=allspecies)
summary_fit = summary(fit) # show results
fit_sum=rbind(fit_sum,
              c('simgene_species',coefficients(fit),summary_fit$adj.r.squared))

fit <- lm(simgene_all ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
            avg_gene_length + total_gene, data=allspecies)
summary_fit = summary(fit) # show results
fit_sum=rbind(fit_sum,
              c('simgene_all',coefficients(fit),summary_fit$adj.r.squared))

fit <- lm(HS_NO_all ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
            avg_gene_length + total_gene, data=allspecies)
summary_fit = summary(fit) # show results
fit_sum=rbind(fit_sum,
              c('HS_NO_all',coefficients(fit),summary_fit$adj.r.squared))
fit_HSratio = TRUE
if(fit_HSratio)
{fit <- lm(HSraio_sim ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
            avg_gene_length + total_gene, data=allspecies)
summary_fit = summary(fit) # show results
fit_sum=rbind(fit_sum,
              c('HSraio_sim',coefficients(fit),summary_fit$adj.r.squared))

fit <- lm(HSraio_sim2 ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
            avg_gene_length + total_gene, data=allspecies)
summary_fit = summary(fit) # show results
fit_sum=rbind(fit_sum,
              c('HSraio_sim2',coefficients(fit),summary_fit$adj.r.squared))
fit <- lm(HSraio ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
            avg_gene_length + total_gene, data=allspecies)
summary_fit = summary(fit) # show results
fit_sum=rbind(fit_sum,
              c('HSraio',coefficients(fit),summary_fit$adj.r.squared))

fit <- lm(HSraio_gene ~ logabu + NO.SNP_lineage + No.SNP + avg_gene + 
            avg_gene_length + total_gene, data=allspecies)
summary_fit = summary(fit) # show results
fit_sum=rbind(fit_sum,
              c('HSraio_gene',coefficients(fit),summary_fit$adj.r.squared))

}
colnames(fit_sum)[1]='parameter'
colnames(fit_sum)[9]='adj.r.squared'
write.table(fit_sum,'final_SNP/HSgene/all.species.all.species.lineagenum.abu.dnds.sim.fitmodel.txt',
            quote=F,sep='\t',row.names=F)

#plot
allspecies=read.delim('final_SNP/HSgene/all.species.all.species.lineagenum.abu.dnds.sim.txt',header=T)
allspecies$dNdS_cut_HS[is.na(allspecies$dNdS_cut_HS)]=0
allspecies2=allspecies
Not_qualify = c()#c('BaVu','BaTh')#'ClBe','BaVu'
Out_lier = c('ClBe')
# plot add error bar!!!
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(plyr)
library(stringr)

sample1 = 553
sample2 = 564
allspecies = allspecies2[which(allspecies2$HS_NO > 0),]
allspecies=allspecies[which(!(allspecies$species %in% Not_qualify |
                                allspecies$species %in% Out_lier )),]

ggscatter(allspecies, x = "HSraio", y = "HSraio_gene",
          # add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(1)],
          size = 3,
          label = 'species_short.y',
          title = 'genes-target parallel evolution',
          xlab = 'No. SNPs under parallel evolution (%)',
          ylab = 'No. genes under parallel evolution (%)',
          add.params = list(color = brewer.pal(5, "PRGn")[c(1)]
                            , fill =brewer.pal(5, "PRGn")[c(2)]),
          cor.coef = TRUE, cor.method = "spearman")

allspecies = allspecies2
allspecies$logabuy = log10(allspecies$abu_mean.y)
allspecies$species[which(allspecies$logabuy< -3)]
allspecies$logabu = log10(allspecies$abu_mean.x)
allspecies$species[which(allspecies$logabu< -3)]
allspecies=allspecies[which(!(allspecies$species %in% Not_qualify |
                                allspecies$species %in% Out_lier )),]
allspecies$HSraio_sim = allspecies$Mean_sim_SNPs/allspecies$No.SNP
allspecies$HSraio_sim2 = allspecies$Mean_sim/allspecies$No.gene_SNPs
allspecies$NO.SNP_lineage = allspecies$No.SNP/allspecies$num_pop
allspecies$simSNP_lineage = allspecies$Mean_sim_lineage_SNPs/allspecies$num_pop
allspecies$simSNP_species = allspecies$Mean_sim_species_SNPs/allspecies$num_pop
allspecies$simSNP_all = allspecies$Mean_sim_SNPs/allspecies$num_pop
allspecies$HS_SNP_NO_all = allspecies$HS_SNP_NO/allspecies$num_pop

allspecies$simgene_lineage = allspecies$Mean_sim_lineage/allspecies$num_pop
allspecies$simgene_species = allspecies$Mean_sim_species/allspecies$num_pop
allspecies$simgene_all = allspecies$Mean_sim/allspecies$num_pop
allspecies$HS_NO_all = allspecies$HS_NO/allspecies$num_pop

fit_sum=read.delim('final_SNP/HSgene/all.species.all.species.lineagenum.abu.dnds.sim.fitmodel.txt',header=T)
fit_sum=fit_sum[which(fit_sum$adj.r.squared > 0.5),]
diff = c(
  quantile(allspecies$logabu[!is.na(allspecies$logabu)],c(0.95))-quantile(allspecies$logabu[!is.na(allspecies$logabu)],c(0.05)),
  #quantile(allspecies$NO.SNP_lineage,c(0.95))-quantile(allspecies$NO.SNP_lineage,c(0.05)),
  #quantile(allspecies$No.SNP,c(0.95))-quantile(allspecies$No.SNP,c(0.05)),
  quantile(allspecies$avg_gene,c(0.95))-quantile(allspecies$avg_gene,c(0.05)),
  quantile(allspecies$avg_gene_length,c(0.95))-quantile(allspecies$avg_gene_length,c(0.05)),
  quantile(allspecies$total_gene,c(0.95))-quantile(allspecies$total_gene,c(0.05))
)
library(pheatmap)
mat = fit_sum[1:3,c(3,6:8)]
for(i in 1:ncol(mat))
{
  mat[,i]=mat[,i]*diff[i]
}
row.names(mat)=fit_sum$parameter[1:3]
mat=as.matrix(mat)
#mat = log10(mat+2)
color_set = colorRampPalette(brewer.pal(n = 9, name = "Purples"))(9)[c(1:9)]
pheatmap(mat= mat,
                color = color_set,
                show_rownames= TRUE,
                show_colnames = TRUE,
                cluster_cols = FALSE,
               #clustering_distance_rows = 'correlation',
                cluster_rows = TRUE)

mat = fit_sum[4:6,c(3,6:8)]
for(i in 1:ncol(mat))
{
  mat[,i]=mat[,i]*diff[i]
}
row.names(mat)=fit_sum$parameter[4:6]
mat=as.matrix(mat)
#mat = log10(mat+2)
color_set = colorRampPalette(brewer.pal(n = 9, name = "Purples"))(9)[c(1:9)]
pheatmap(mat= mat,
         color = color_set,
         show_rownames= TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         #clustering_distance_rows = 'correlation',
         cluster_rows = TRUE)
# plot
ggscatter(allspecies, x = "Mean_sim_lineage_SNPs", y = "No.SNP",
          # add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(1)],
          size = 3,
          label = 'species_short.y',
          title = 'genes-target parallel evolution',
          xlab = 'median abundance (log10)',
          ylab = 'median abundance BN10 (log10)',
          add.params = list(color = brewer.pal(5, "PRGn")[c(1)]
                            , fill =brewer.pal(5, "PRGn")[c(2)]),
          cor.coef = TRUE, cor.method = "spearman")

ggscatter(allspecies, x = "logabu", y = "logabuy",
          # add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(1)],
          size = 3,
          label = 'species_short.y',
          title = 'genes-target parallel evolution',
          xlab = 'median abundance (log10)',
          ylab = 'median abundance BN10 (log10)',
          add.params = list(color = brewer.pal(5, "PRGn")[c(1)]
                            , fill =brewer.pal(5, "PRGn")[c(2)]),
          cor.coef = TRUE, cor.method = "spearman")

ggscatter(allspecies, x = "logabu", y = "HSraio",
         # add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(1)],
          size = 3,
          label = 'species_short.y',
          title = 'genes-target parallel evolution',
          xlab = 'median abundance (log10)',
          ylab = 'No. SNPs under parallel evolution (%)',
          add.params = list(color = brewer.pal(5, "PRGn")[c(1)]
                            , fill =brewer.pal(5, "PRGn")[c(2)]),
          cor.coef = TRUE, cor.method = "spearman")

ggscatter(allspecies, x = "logabu", y = "HSraio_sim",
          # add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(5)],
          size = 3,
          label = 'species_short.y',
          title = 'genes-target parallel evolution',
          xlab = 'median abundance (log10)',
          ylab = 'No. SNPs under parallel evolution simulation (%)',
          add.params = list(color = brewer.pal(5, "PRGn")[c(1)]
                            , fill =brewer.pal(5, "PRGn")[c(2)]),
          cor.coef = TRUE, cor.method = "spearman")

ggscatter(allspecies, x = "logabu", y = "HSraio_gene",
          # add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(1)],
          size = 3,
          label = 'species_short.y',
          title = 'genes-target parallel evolution',
          xlab = 'median abundance (log10)',
          ylab = 'No. genes under parallel evolution (%)',
          add.params = list(color = brewer.pal(5, "PRGn")[c(1)]
                            , fill =brewer.pal(5, "PRGn")[c(2)]),
          cor.coef = TRUE, cor.method = "spearman")

ggscatter(allspecies, x = "logabu", y = "HSraio_sim2",
          # add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(5)],
          size = 3,
          label = 'species_short.y',
          title = 'genes-target parallel evolution',
          xlab = 'median abundance (log10)',
          ylab = 'No. genes under parallel evolution simulation (%)',
          add.params = list(color = brewer.pal(5, "PRGn")[c(1)]
                            , fill =brewer.pal(5, "PRGn")[c(2)]),
          cor.coef = TRUE, cor.method = "spearman")

ggscatter(allspecies, x = "logabu", y = "dNdS_cut_HS",
          #add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(1)],
          size = 3,
          label = 'species_short.y',
          title = 'genes-target parallel evolution',
          xlab = 'mean abundance (log10)',
          ylab = 'dN/dS of genes under parallel evolution',
          add.params = list(color = brewer.pal(5, "PRGn")[c(1)]
                            , fill =brewer.pal(5, "PRGn")[c(2)]),
          cor.coef = TRUE, cor.method = "spearman")

ggscatter(allspecies, x = "logabu", y = "dNdS_cut_all",
          #add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(5)],
          size = 3,
          label = 'species_short.y',
          title = 'de novo mutations',
          xlab = 'mean abundance (log10)',
          ylab = 'dN/dS of all genes with de novo mutations',
          add.params = list(color = brewer.pal(5, "PRGn")[c(5)]
                            , fill =brewer.pal(5, "PRGn")[c(4)]),
          cor.coef = TRUE, cor.method = "spearman")

ggscatter(allspecies, x = "logabu", y = "dmrca",
          #add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(5)],
          size = 3,
          label = 'species_short.y',
          title = 'de novo mutations',
          xlab = 'mean abundance (log10)',
          ylab = 'dmrca',
          add.params = list(color = brewer.pal(5, "PRGn")[c(5)]
                            , fill =brewer.pal(5, "PRGn")[c(4)]),
          cor.coef = TRUE, cor.method = "spearman")

ggscatter(allspecies[which(allspecies$dNdS_cut_HS>0),], 
          x = "dmrca", y = "dNdS_cut_HS",
          # add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(1)],
          size = 3,
          label = 'species_short.y',
          title = 'genes-target parallel evolution',
          xlab = 'dmrca',
          ylab = 'dNdS',
          add.params = list(color = brewer.pal(5, "PRGn")[c(5)]
                            , fill =brewer.pal(5, "PRGn")[c(4)]),
          cor.coef = TRUE, cor.method = "spearman")

ggscatter(allspecies, x = "dmrca", y = "dNdS_cut_all",
          # add = "reg.line", conf.int = TRUE,
          color = brewer.pal(5, "PRGn")[c(5)],
          size = 3,
          label = 'species_short.y',
          title = 'de novo mutations',
          xlab = 'dmrca',
          ylab = 'dNdS',
          add.params = list(color = brewer.pal(5, "PRGn")[c(5)]
                            , fill =brewer.pal(5, "PRGn")[c(4)]),
          cor.coef = TRUE, cor.method = "spearman")

allspecies = allspecies2#[which(allspecies2$HS_NO > 0),]
allspecies=allspecies[which(!(allspecies$species %in% Not_qualify)),]
allspecies$genus =str_split_fixed(allspecies$species_short.y, "_", 2)[,1]
allspecies3=allspecies[which(allspecies$genus %in% c('Bifido','Bacter','Paraba')),]
allspecies3$logabu = log10(allspecies3$abu_mean.x)

ggplot(allspecies3,aes(
  y=dNdS_cut_HS,
  x = logabu))+
  geom_point(aes(
    color=genus,
    fill=genus
  ),
  alpha = 1,
  size=2
  )+
  geom_text(aes(
              color=genus,
              fill=genus,
              label=species_short.y
            ),
            hjust = 0,nudge_x = 0.02, nudge_y = 0)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  scale_fill_manual(values =brewer.pal(5, "Set1")[c(1,2,3)])+
  scale_color_manual(values =brewer.pal(5, "Set1")[c(1,2,3)])+
  xlab('mean abundance (log10)') + ylab('dN/dS of genes under parallel evolution')

ggplot(allspecies3,aes(
  y=HSraio,
  x = logabu))+
  geom_point(aes(
    color=genus,
    fill=genus
  ),
  alpha = 1,
  size=2
  )+
  geom_text(aes(
    color=genus,
    fill=genus,
    label=species_short.y
  ),
  hjust = 0,nudge_x = 0.02, nudge_y = 0)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  scale_fill_manual(values =brewer.pal(5, "Set1")[c(1,2,3)])+
  scale_color_manual(values =brewer.pal(5, "Set1")[c(1,2,3)])+
  xlab('mean abundance (log10)') + ylab('No. SNPs under parallel evolution (%)')

ggplot(allspecies3,aes(
  y=HSraio_gene,
  x = logabu))+
  geom_point(aes(
    color=genus,
    fill=genus
  ),
  alpha = 1,
  size=2
  )+
  geom_text(aes(
    color=genus,
    fill=genus,
    label=species_short.y
  ),
  hjust = 0,nudge_x = 0.02, nudge_y = 0)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  scale_fill_manual(values =brewer.pal(5, "Set1")[c(1,2,3)])+
  scale_color_manual(values =brewer.pal(5, "Set1")[c(1,2,3)])+
  xlab('mean abundance (log10)') + ylab('No. genes under parallel evolution (%)')


#----gene core flexible prep----
allgene = read.delim('final_SNP/annotation/all.genome.gene.faa.uc.species.sum',header=T)
genecount = count(allgene,c('genus_num'))
genecount$tag = 'all genes'
genecount$percentage = genecount$freq/nrow(allgene)

denovo = read.delim('final_SNP/annotation/all.denovo.gene.faa.uc.species.sum',header=T)
denovocount = count(denovo,c('genus_num'))
denovocount$tag = 'de novo mutations'
denovocount$percentage = denovocount$freq/nrow(denovo)

HS = read.delim('final_SNP/annotation/all.selected.gene.faa.High_select2.faa.uc.species.sum',header=T)
HScount = count(HS,c('genus_num'))
HScount$tag = 'parallel evolution'
HScount$percentage = HScount$freq/nrow(HS)

genecount=rbind(genecount,denovocount,HScount)

ggplot(genecount,aes(y=percentage,
                     x=genus_num,
                     group = tag))+
  geom_point(aes(
    color=tag,
    fill=tag
  ),
  alpha = 1,
  size=2
  )+ 
  geom_line(aes(
    color=tag,
    fill=tag
  ),
  alpha = 0.5,
  size=1
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =c(brewer.pal(9, "Set1")[c(9,3,4)]))+
  scale_color_manual(values =c(brewer.pal(9, "Set1")[c(9,3,4)]))

# N and S
allspecies = read.delim('final_SNP/HSgene/all.species.High_select2.all.gene.txt',header=T)
allchangename = read.delim('final_SNP/annotation/all.denovo.gene.faa.changename.txt',header=F)
colnames(allchangename)=c('X.donor_species','short_name_donor','gene','gene_newname','V5')
length(unique(allchangename$gene))==nrow(allchangename)
allchangename$gene_donor = paste(allchangename$X.donor_species,allchangename$gene)
allspecies$gene_donor = paste(allspecies$X.donor_species,allspecies$gene)
allspecies=merge(allspecies,allchangename[,c(4,6)],
                 by = 'gene_donor')
allspecies=merge(allspecies,denovo[,c(-5)],
                 by.x='gene_newname',
                 by.y='record')
write.table(allspecies,'final_SNP/HSgene/all.species.High_select2.all.gene.speciesnum.txt',
            quote=F,sep='\t',row.names=F)

allspecies = read.delim('final_SNP/HSgene/all.species.High_select2.all.gene.speciesnum.txt',header=T)
allspecies$S_cut = as.numeric(as.matrix(allspecies$S))
allspecies$S_cut[which(allspecies$S_cut==0)]=0.5
allspecies$cross_genus = 'False'
allspecies$cross_genus[which(allspecies$genus_num>1)]='True'
# N
allspecies_count= count(allspecies,c('N','cross_genus'))
allspecies_count2 = count(allspecies,c('cross_genus'))
allspecies_count=merge(allspecies_count,allspecies_count2,
                       by='cross_genus')
allspecies_count$percentage = allspecies_count$freq.x/allspecies_count$freq.y
ggplot(allspecies_count,aes(x=N,
                     y=percentage,
                     group = cross_genus
))+
  geom_bar(aes(
    color=cross_genus,
    fill=cross_genus
  ),
  stat = 'identity',
  position = 'dodge')+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1,5)])

# S
allspecies_count= count(allspecies,c('S','cross_genus'))
allspecies_count2 = count(allspecies,c('cross_genus'))
allspecies_count=merge(allspecies_count,allspecies_count2,
                       by='cross_genus')
allspecies_count$percentage = allspecies_count$freq.x/allspecies_count$freq.y
ggplot(allspecies_count,aes(x=S,
                            y=percentage,
                            group = cross_genus
))+
  geom_bar(aes(
    color=cross_genus,
    fill=cross_genus
  ),
stat = 'identity',
position = 'dodge')+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1,5)])
# overall
ggplot(allspecies,aes(x=cross_genus,
                      y=N/S_cut
))+
  geom_boxplot(aes(
    color=cross_genus,
    fill=cross_genus
  ),
  outlier.shape = NA,
  alpha = 0,
  size=0
  )+ stat_summary(fun = mean, geom = "bar", 
                  aes(color=cross_genus,fill = cross_genus),alpha = 0.5,
                  size=0)+
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               aes(color=cross_genus,
                   width = 0.8,
                   linetype = "solid"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_y_continuous(breaks = seq(0,1.6,0.2))

allspecies$cross_genus = 'flexible'
allspecies$cross_genus[which(allspecies$genus_num>1)]='core'
write.table(allspecies,'final_SNP/HSgene/all.species.High_select2.all.gene.speciesnum.core_flexible.txt',
            quote=F,sep='\t',row.names=F)

# eggnog all
eggnogall = read.delim('final_SNP/annotation/all.genome.gene.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
eggnogall_count = count(eggnogall,c('COG','COG1','COG2','cross_genus'))
eggnogall_count2 = count(eggnogall,c('COG','COG1','COG2'))
eggnogall_count$COG1COG2 = paste(eggnogall_count$COG1,eggnogall_count$COG2)
eggnogall_count2$COG1COG2 = paste(eggnogall_count2$COG1,eggnogall_count2$COG2)
eggnogall_count = merge(eggnogall_count,eggnogall_count2[,c(4,5)],by='COG1COG2')
eggnogall_count$percentage = eggnogall_count$freq.x/eggnogall_count$freq.y
eggnogall_count$tag = 'all genes'

eggnogall_countall = count(eggnogall,c('COG1','cross_genus'))
eggnogall_countall2=count(eggnogall,'COG1')
eggnogall_countall = merge(eggnogall_countall,eggnogall_countall2,by='COG1')
eggnogall_countall$percentage = eggnogall_countall$freq.x/eggnogall_countall$freq.y

# eggnog de novo
eggnog = read.delim('final_SNP/annotation/all.denovo.gene.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
eggnog2 = eggnog
eggnog = eggnog2
eggnog_count = count(eggnog,c('COG','COG1','COG2','cross_genus'))
eggnog_count2 = count(eggnog,c('COG','COG1','COG2'))
eggnog_count$COG1COG2 = paste(eggnog_count$COG1,eggnog_count$COG2)
eggnog_count2$COG1COG2 = paste(eggnog_count2$COG1,eggnog_count2$COG2)
eggnog_count = merge(eggnog_count,eggnog_count2[,c(4,5)],by='COG1COG2')
eggnog_count$percentage = eggnog_count$freq.x/eggnog_count$freq.y
eggnog_count$tag = 'de novo mutations'
eggnogall_count=rbind(eggnog_count,eggnogall_count)
# eggnog HS
eggnog = read.delim('final_SNP/annotation/all.selected.gene.faa.High_select2.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
eggnog2 = eggnog
eggnog = eggnog2
eggnog_count = count(eggnog,c('COG','COG1','COG2','cross_genus'))
eggnog_count2 = count(eggnog,c('COG','COG1','COG2'))
eggnog_count$COG1COG2 = paste(eggnog_count$COG1,eggnog_count$COG2)
eggnog_count2$COG1COG2 = paste(eggnog_count2$COG1,eggnog_count2$COG2)
eggnog_count = merge(eggnog_count,eggnog_count2[,c(4,5)],by='COG1COG2')
eggnog_count$percentage = eggnog_count$freq.x/eggnog_count$freq.y
eggnog_count$tag = 'parallel evolution'
eggnogall_count=rbind(eggnog_count,eggnogall_count)
write.table(eggnogall_count,'final_SNP/annotation/all.core.flexible.sum',
            quote=F,sep='\t',row.names=F)
#
#----gene core flexible plot----
eggnogall_count=read.delim('final_SNP/annotation/all.core.flexible.sum',header=T)
eggnogall_count2=eggnogall_count[which(eggnogall_count$cross_genus=='True'),]
temp_row = eggnogall_count[which(!(eggnogall_count$COG1COG2 %in% eggnogall_count2$COG1COG2)),]
temp_row$percentage = 0
temp_row$cross_genus = 'True'
eggnogall_count2 = rbind(eggnogall_count2,temp_row)
eggnogall_count2 = eggnogall_count2[order(eggnogall_count2$COG1COG2),]
eggnogall_count2 = eggnogall_count2[which(eggnogall_count2$COG != 'S'),]
temp_row = eggnogall_count2[which(eggnogall_count2$tag == 'all genes'),]
temp_row$percentage = 0
temp_row$tag = 'parallel evolution'
eggnogall_count2 = rbind(eggnogall_count2,temp_row)
eggnogall_count2 = eggnogall_count2[!duplicated(paste(eggnogall_count2$COG1COG2,
                                                      eggnogall_count2$tag)),]
temp_row = eggnogall_count2[which(eggnogall_count2$tag == 'all genes'),]
temp_row$percentage = 0
temp_row$tag = 'de novo mutations'
eggnogall_count2 = rbind(eggnogall_count2,temp_row)
eggnogall_count2 = eggnogall_count2[!duplicated(paste(eggnogall_count2$COG1COG2,
                                                      eggnogall_count2$tag)),]

eggnogall_count2 = eggnogall_count2[order(eggnogall_count2$COG1COG2),]
eggnogall_count2 = eggnogall_count2[which(eggnogall_count2$COG != 'S'),]

ggplot(eggnogall_count2,aes(y=percentage,
                     x=COG1COG2,
                     group = tag))+
  geom_point(aes(
    color=tag,
    fill=tag
  ),
  alpha = 1,
  size=2
  )+ 
  geom_line(aes(
    color=tag,
    fill=tag
  ),
  alpha = 0.5,
  size=1
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =c(brewer.pal(9, "Set1")[c(9,3,4)]))+
  scale_color_manual(values =c(brewer.pal(9, "Set1")[c(9,3,4)]))

alltag = unique(eggnogall_count2$tag)
alltag2 = unique(eggnogall_count2$COG1)
eggnogall_count3 = eggnogall_count[which(eggnogall_count$cross_genus == 'False'),]
eggnogall_count3 = eggnogall_count3[which(eggnogall_count3$COG !='S'),]
tag1 = alltag[2]
{tag1
sum(eggnogall_count3$freq.y[which(eggnogall_count3$tag==tag1)]
    ) - sum(eggnogall_count3$freq.x[which(eggnogall_count3$tag==tag1)])
sum(eggnogall_count3$freq.y[which(eggnogall_count3$tag==tag1)])
100 - sum(eggnogall_count3$freq.x[which(eggnogall_count3$tag==tag1)])/sum(
  eggnogall_count3$freq.y[which(eggnogall_count3$tag==tag1)])*100
}
for(tag2 in alltag2)
{
  
  eggnogall_count3_temp = eggnogall_count3[which(eggnogall_count3$tag==tag1
                                                 &
                                                   eggnogall_count3$COG1==tag2 ),]
 temp_result =  c(toString(tag1),tag2,sum(eggnogall_count3_temp$freq.y
  ) - sum(eggnogall_count3_temp$freq.x),
  sum(eggnogall_count3_temp$freq.y),
  100 - sum(eggnogall_count3_temp$freq.x)/sum(
    eggnogall_count3_temp$freq.y)*100)
  
  print(temp_result)
  
}
#
#----eggnog enrichment test prep----
tomatrix <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=mat[rownum,colnum] + 1/dm[i,colnom]
  }
  return(mat)
}
tomatrix_norm <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=mat[rownum,colnum] + dm[i,colnom]
  }
  return(mat)
}
tomatrix_norm_test <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  dm[which(dm[,colnom]=='TRUE'),colnom]=1
  dm[which(dm[,colnom]=='FALSE'),colnom]=0
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=as.numeric(dm[i,colnom])
  }
  return(mat)
}
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(stringr)
#against all genes!!!
# load all genes
keggall = read.delim('final_SNP/annotation/all.genome.gene.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
keggall$Tag = 'All'
# HS genes
kegg = read.delim('final_SNP/annotation/all.selected.gene.faa.High_select2.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
kegg$Tag = 'High_select'
kegg=rbind(kegg,keggall)
kegg$speciesnew = kegg$species
kegg$speciesnew = str_split_fixed(kegg$speciesnew, "_", 2)[,1]
kegg$speciesnew[which(kegg$speciesnew=='95')]='EsCo'
kegg$speciesnew[which(kegg$speciesnew=='Bfragilis')]='BaFr'
kegg$species = kegg$speciesnew
kegg=kegg[which(kegg$COG!=''),]
kegg=kegg[which(kegg$COG1!='POORLY CHARACTERIZED'),]
#normalize against all
newTag='High_select'
kegg2=kegg[which(kegg$Tag==newTag),]
kegg_HS = kegg2
kegg2_count = count(kegg2,c('species','COG1','COG','COG2','cross_genus'))
colnames(kegg2_count)[c(ncol(kegg2_count))]=c('num_HS')
kegg2_count$speciesfun = paste(kegg2_count$species,
                               kegg2_count$cross_genus,
                               kegg2_count$COG)
kegg2new=kegg[which(kegg$Tag!=newTag),]
kegg2_countall = count(kegg2new,c('species','COG1','COG','COG2','cross_genus'))
kegg4=count(kegg2new,c('species','cross_genus'))
kegg2_countall$speciesgenus = paste(kegg2_countall$species,kegg2_countall$cross_genus)
kegg4$speciesgenus = paste(kegg4$species,kegg4$cross_genus)
kegg2_countall=merge(kegg2_countall,kegg4[,c(3,4)],by='speciesgenus')
colnames(kegg2_countall)[c(ncol(kegg2_countall)-1,
                           ncol(kegg2_countall))]=c(
                             'num_gene','num_allgene')
kegg4=count(kegg2,c('species','cross_genus'))
kegg4$speciesgenus = paste(kegg4$species,kegg4$cross_genus)
kegg2_countall=merge(kegg2_countall,kegg4[,c(3,4)],by='speciesgenus')

colnames(kegg2_countall)[ncol(kegg2_countall)]=c('num_allHS')
kegg2_countall$expectedHS = kegg2_countall$num_gene/kegg2_countall$num_allgene*kegg2_countall$num_allHS
kegg2_countall$speciesfun = paste(kegg2_countall$species,
                                  kegg2_countall$cross_genus,
                                  kegg2_countall$COG)
m=ncol(kegg2_count)
kegg2_count=merge(kegg2_count[,c(m-1,m)],kegg2_countall,
                  by='speciesfun',all=T)
kegg2_count$num_HS[is.na(kegg2_count$num_HS)]=0
kegg2_count$diff = (kegg2_count$num_HS/kegg2_count$expectedHS)
kegg2_count=kegg2_count[which(kegg2_count$species
                              %in% kegg2$species),]
for(i in 1:nrow(kegg2_count))
{
  a = kegg2_count$num_HS[i]
  b = kegg2_count$expectedHS[i]
  kegg2_count$diff[i] = max(a,b)/min(a,b)*(a-b)/abs((a-b))
}
write.table(kegg2_count,'final_SNP/annotation/all.HS.eggnog.norm.species.txt',
            quote=F,sep='\t',row.names=F)

# de novo genes
kegg = read.delim('final_SNP/annotation/all.denovo.gene.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
kegg$Tag = 'denovo'
kegg=rbind(kegg,keggall)
kegg$speciesnew = kegg$species
kegg$speciesnew = str_split_fixed(kegg$speciesnew, "_", 2)[,1]
kegg$speciesnew[which(kegg$speciesnew=='95')]='EsCo'
kegg$speciesnew[which(kegg$speciesnew=='Bfragilis')]='BaFr'
kegg$species = kegg$speciesnew
kegg=kegg[which(kegg$COG!=''),]
kegg=kegg[which(kegg$COG1!='POORLY CHARACTERIZED'),]
#normalize against all
newTag='denovo'
kegg2=kegg[which(kegg$Tag==newTag),]
kegg_HS = kegg2
kegg2_count = count(kegg2,c('species','COG1','COG','COG2','cross_genus'))
colnames(kegg2_count)[c(ncol(kegg2_count))]=c('num_HS')
kegg2_count$speciesfun = paste(kegg2_count$species,
                               kegg2_count$cross_genus,
                               kegg2_count$COG)
kegg2new=kegg[which(kegg$Tag!=newTag),]
kegg2_countall = count(kegg2new,c('species','COG1','COG','COG2','cross_genus'))
kegg4=count(kegg2new,c('species','cross_genus'))
kegg2_countall$speciesgenus = paste(kegg2_countall$species,kegg2_countall$cross_genus)
kegg4$speciesgenus = paste(kegg4$species,kegg4$cross_genus)
kegg2_countall=merge(kegg2_countall,kegg4[,c(3,4)],by='speciesgenus')
colnames(kegg2_countall)[c(ncol(kegg2_countall)-1,ncol(kegg2_countall))]=c('num_gene','num_allgene')
kegg4=count(kegg2,c('species','cross_genus'))
kegg4$speciesgenus = paste(kegg4$species,kegg4$cross_genus)
kegg2_countall=merge(kegg2_countall,kegg4[,c(3,4)],by='speciesgenus')
colnames(kegg2_countall)[ncol(kegg2_countall)]=c('num_allHS')
kegg2_countall$expectedHS = kegg2_countall$num_gene/kegg2_countall$num_allgene*kegg2_countall$num_allHS
kegg2_countall$speciesfun = paste(kegg2_countall$species,
                                  kegg2_countall$cross_genus,
                                  kegg2_countall$COG)
m=ncol(kegg2_count)
kegg2_count=merge(kegg2_count[,c(m-1,m)],kegg2_countall,
                  by='speciesfun',all=T)
kegg2_count$num_HS[is.na(kegg2_count$num_HS)]=0
kegg2_count$diff = (kegg2_count$num_HS/kegg2_count$expectedHS)
kegg2_count=kegg2_count[which(kegg2_count$species
                              %in% kegg2$species),]
for(i in 1:nrow(kegg2_count))
{
  a = kegg2_count$num_HS[i]
  b = kegg2_count$expectedHS[i]
  kegg2_count$diff[i] = max(a,b)/min(a,b)*(a-b)/abs((a-b))
}
write.table(kegg2_count,'final_SNP/annotation/all.denovo.eggnog.norm.species.txt',
            quote=F,sep='\t',row.names=F)
#plot
library(pheatmap)
kegg2_count=read.delim('final_SNP/annotation/all.denovo.eggnog.norm.species.txt',header=T)
hist(kegg2_count$diff)
kegg2_count$diff[is.infinite(kegg2_count$diff)]=-2000
kegg2_count2 = kegg2_count
kegg2_count2=kegg2_count2[order(kegg2_count2$COG2),]
kegg2_count2=kegg2_count2[order(kegg2_count2$COG1),]
kegg2_count2$COG3_new = paste(kegg2_count2$COG1,kegg2_count2$COG2)
kegg2_count2=kegg2_count2[order(kegg2_count2$COG3_new),]
kegg2_mat = tomatrix_norm(kegg2_count2,14,4,13)
hist(kegg2_mat)
kegg2_mat[which(kegg2_mat==0)]=-2000
max_mat = as.integer(max(kegg2_mat))+1
min_mat =  as.integer(min(kegg2_mat))
hist(kegg2_count2$diff[which(kegg2_count2$diff>=-20)]
     ,breaks=seq(-20,max(kegg2_mat)+1,1))
mat_breaks = c(-10,-5,-3,-2,-1,1,2,3,5,10)
color_set = colorRampPalette(brewer.pal(n = 9, name = "PRGn"))(9)[c(1:9)]
pplot<-pheatmap(mat= kegg2_mat,
                color = color_set,
                breaks = mat_breaks,
                legend_breaks = mat_breaks,
                legend_labels = mat_breaks,
                show_rownames= TRUE,
                show_colnames = TRUE,
                cluster_cols = FALSE,
               # clustering_distance_cols = 'correlation',
                cluster_rows = FALSE)
# fisher test, binomial test, one way anova test
kegg2_count=read.delim('final_SNP/annotation/all.HS.eggnog.norm.species.txt',header=T)
#kegg2_count=read.delim('final_SNP/annotation/all.denovo.eggnog.norm.species.txt',header=T)
result_test=matrix(0,nrow=nrow(kegg2_count),ncol=6)
colnames(result_test)=c('species',
                        'fun',
                        'cross_genus',
                        'fisher_test',
                        'chisq_test','binomial_test')
k=0
for(i in 1:nrow(kegg2_count))
{
  k=k+1
  temp=matrix(0,nrow=2,ncol=2)
  temp[1,1]=kegg2_count$num_HS[i]
  temp[1,2]=kegg2_count$num_gene[i]-temp[1,1]
  temp[2,1]=kegg2_count$num_allHS[i]-temp[1,1]
  temp[2,2]=kegg2_count$num_allgene[i]-kegg2_count$num_gene[i]
  result2 = fisher.test(temp,alternative='greater',
                        simulate.p.value =TRUE,B = 1000,
                        conf.int = TRUE, conf.level = 0.9)
  result3 = chisq.test(temp,
                       simulate.p.value =TRUE,B = 1000)
  result4 = binom.test(temp[1,1], kegg2_count$num_gene[i], 
                       kegg2_count$num_allHS[i]/kegg2_count$num_allgene[i], 
                       alternative="greater")
  result_test[k,]=c(
    as.matrix(kegg2_count$species[i]),
    as.matrix(kegg2_count$COG[i]),
    as.matrix(kegg2_count$cross_genus[i]),
    result2$p.value,
    result3$p.value,
    result4$p.value
  )
}
result_test=data.frame(result_test)
# one way anova test: one COG compared to other COG (average value)
species_fun_sig = c()
#temp_kegg_anova_all = data.frame(
#  fun = c(rep('all',time =sum(kegg2_count$num_gene))),
#  gene_num = c(rep(1,time = sum(kegg2_count$num_HS)),
#               rep(0,time = sum(kegg2_count$num_gene) - sum(kegg2_count$num_HS))
#  ))
for(species in unique(kegg2_count$species))
{
  for(cross_genus in c('False','True'))
  {
    temp_kegg = kegg2_count[which(kegg2_count$species==species &
                                  kegg2_count$cross_genus==cross_genus ),]
  if(nrow(temp_kegg[which(temp_kegg$num_HS>0),])>10)
  {for(i in 1:nrow(temp_kegg))
  {
    if(temp_kegg$num_allgene[i]-temp_kegg$num_HS[i] > 0)
   { temp_kegg_anova = data.frame(
      fun = c(rep(toString(temp_kegg$COG[i]),time = temp_kegg$num_gene[i]),
              rep('other',time = temp_kegg$num_allgene[i]-temp_kegg$num_gene[i])),
      gene_num = c(rep(1,time = temp_kegg$num_HS[i]),
                    rep(0,time = (temp_kegg$num_allgene[i]-temp_kegg$num_HS[i]))
                   ))
   temp_kegg_anova = data.frame(
     fun = c(rep(toString(temp_kegg$COG[i]),time = temp_kegg$num_gene[i]),
             rep('other',time = temp_kegg$num_allgene[i]-temp_kegg$num_gene[i])),
     gene_num = c(rep(1,time = temp_kegg$num_HS[i]),
                  rep(0,time = (temp_kegg$num_gene[i]-temp_kegg$num_HS[i])),
                  rep(1,time = temp_kegg$num_allHS[i]-temp_kegg$num_HS[i]),
                  rep(0,time = (temp_kegg$num_allgene[i]-temp_kegg$num_gene[i]
                                -temp_kegg$num_allHS[i]+temp_kegg$num_HS[i]))
     ))
    res.aov <- aov(gene_num ~ fun, data = temp_kegg_anova)
    summary(res.aov)
    result4 = summary(res.aov)
    if(!is.null(result4[[1]][["Pr(>F)"]][1]))
      if(!is.nan(result4[[1]][["Pr(>F)"]][1]))
      if(result4[[1]][["Pr(>F)"]][1]<=0.05)
      species_fun_sig=c(species_fun_sig,paste(
        species,
        cross_genus,
        toString(temp_kegg$COG[i])
      ))}
  }}
  else if(nrow(temp_kegg)>1)
    for(i in 1:nrow(temp_kegg))
    {
      if(abs(temp_kegg$diff[i])>=1.5)
      species_fun_sig=c(species_fun_sig,paste(
        species,
        cross_genus,
        toString(temp_kegg$COG[i])
      ))
    }}
}
kegg2_count$anova_test = FALSE
kegg2_count$anova_test[which(kegg2_count$speciesfun %in% species_fun_sig)]=TRUE
result_test$speciesfun = paste(result_test$species,result_test$cross_genus,result_test$fun)
kegg2_count=merge(kegg2_count,
                  result_test[,c(4:7)],
                  by = 'speciesfun')
kegg2_count$anova_test[which(as.numeric(as.matrix(kegg2_count$fisher_test))<=0.05)]=TRUE
kegg2_count$anova_test[which(as.numeric(as.matrix(kegg2_count$diff))<1)]=FALSE
write.table(kegg2_count,'final_SNP/annotation/all.HS.eggnog.norm.species.test.txt',
            quote=F,sep='\t',row.names=F)
#write.table(kegg2_count,'final_SNP/annotation/all.denovo.eggnog.norm.species.test.txt',
#            quote=F,sep='\t',row.names=F) # has no true anova

#----eggnog enrichment test heatmap ----
tomatrix <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=mat[rownum,colnum] + 1/dm[i,colnom]
  }
  return(mat)
}
tomatrix_norm <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=mat[rownum,colnum] + dm[i,colnom]
  }
  return(mat)
}
tomatrix_norm_test <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  dm[which(dm[,colnom]=='TRUE'),colnom]=1
  dm[which(dm[,colnom]=='FALSE'),colnom]=0
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=as.numeric(dm[i,colnom])
  }
  return(mat)
}
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(stringr)
#heatmap
kegg2_count=read.delim('final_SNP/annotation/all.HS.eggnog.norm.species.test.txt',header=T)
#kegg2_count=read.delim('final_SNP/annotation/all.denovo.eggnog.norm.species.test.txt',header=T)
kegg2_count$anova_test[which(as.numeric(as.matrix(kegg2_count$diff))<1)]='FALSE'
hist(kegg2_count$diff)
kegg2_count$diff[is.infinite(kegg2_count$diff)]=-10
#Not_qualify = c('BaVu','BaTh')#'ClBe'
Not_qualify = c()
#kegg2_count2 = kegg2_count[which(!(kegg2_count$species %in% Not_qualify)),]
# Not core
kegg2_count2 = kegg2_count[which(kegg2_count$cross_genus=='False'),]
kegg2_count2=kegg2_count2[order(kegg2_count2$COG2),]
kegg2_count2=kegg2_count2[order(kegg2_count2$COG1),]
kegg2_count2$COG3_new = paste(kegg2_count2$COG1,kegg2_count2$COG2)
kegg2_count2=kegg2_count2[order(kegg2_count2$COG3_new),]
library(pheatmap)
kegg2_mat = tomatrix_norm(kegg2_count2,18,4,13)
max_mat = as.integer(max(kegg2_mat))+1
min_mat =  as.integer(min(kegg2_mat))
hist(kegg2_count2$diff,breaks=seq(min(kegg2_mat)-1,max(kegg2_mat)+1,1))
mat_breaks = c(min_mat,-5,-3,-2,-1,1,2,3,5,10)
color_set = colorRampPalette(brewer.pal(n = 9, name = "PRGn"))(9)[c(1:3,5,5,5,7:9)]
pplot<-pheatmap(mat= kegg2_mat,
         color = color_set,
         breaks = mat_breaks,
         legend_breaks = mat_breaks,
         legend_labels = mat_breaks,
         show_rownames= TRUE,
         show_colnames = TRUE,
         cluster_cols = TRUE,
         #clustering_distance_cols = 'correlation',
         cluster_rows = FALSE)

# core
new_order = pplot$tree_col[["order"]]
kegg2_count2 = kegg2_count[which(kegg2_count$cross_genus=='True'),]
kegg2_count2=kegg2_count2[order(kegg2_count2$COG2),]
kegg2_count2=kegg2_count2[order(kegg2_count2$COG1),]
kegg2_count2$COG3_new = paste(kegg2_count2$COG1,kegg2_count2$COG2)
kegg2_count2=kegg2_count2[order(kegg2_count2$COG3_new),]
kegg2_mat = tomatrix_norm(kegg2_count2,18,4,13)
max_mat = as.integer(max(kegg2_mat))+1
min_mat =  as.integer(min(kegg2_mat))
#hist(kegg2_count2$diff,breaks=seq(min(kegg2_mat)-1,max(kegg2_mat)+1,1))
mat_breaks = c(min_mat,-5,-3,-2,-1,1,2,3,5,10)
color_set = colorRampPalette(brewer.pal(n = 9, name = "PRGn"))(9)[c(1:3,5,5,5,7:9)]
kegg2_mat = kegg2_mat[,new_order]
pheatmap(mat= kegg2_mat,
                color = color_set,
                breaks = mat_breaks,
                legend_breaks = mat_breaks,
                legend_labels = mat_breaks,
                show_rownames= TRUE,
                show_colnames = TRUE,
                cluster_cols = FALSE,
                cluster_rows = FALSE)
#significant
kegg2_mat = tomatrix_norm_test(kegg2_count2,18,4,14)
kegg2_mat = kegg2_mat[,new_order]

pheatmap(mat= kegg2_mat,
         show_rownames= TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE)

#---- species enrichment gene figure 4 prep----
# HS genes
# anova test + fisher exact: one COG compared to other COG (average value) in the same species
# diff > 1 compared to expected 
kegg2_count=read.delim('final_SNP/annotation/all.HS.eggnog.norm.species.test.txt',header=T)
kegg2_count = kegg2_count[which(kegg2_count$cross_genus!='True'),]
kegg2_count$anova_test[which(as.numeric(as.matrix(kegg2_count$diff))<1)]='FALSE'
kegg2_count=kegg2_count[which(kegg2_count$anova_test=='TRUE'),]
kegg = read.delim('final_SNP/annotation/all.selected.gene.faa.High_select2.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
kegg$Tag = 'High_select'
kegg$speciesnew = kegg$species
kegg$speciesnew = str_split_fixed(kegg$speciesnew, "_", 2)[,1]
kegg$speciesnew[which(kegg$speciesnew=='95')]='EsCo'
kegg$speciesnew[which(kegg$speciesnew=='Bfragilis')]='BaFr'
kegg$species = kegg$speciesnew
kegg=kegg[which(kegg$COG!=''),]
kegg=kegg[which(kegg$COG1!='POORLY CHARACTERIZED'),]
kegg_count = count(
  kegg,c('species','COG','COG1','COG2') 
)
kegg_count$speciesfun = paste(kegg_count$species,kegg_count$COG)
kegg2_count$speciesfun = paste(kegg2_count$species,kegg2_count$COG)
kegg_count=kegg_count[which(kegg_count$freq>1),]
unique(kegg_count$species)
kegg$speciesfun = paste(kegg$species,kegg$COG)
kegg_count_fun = kegg[which(kegg$speciesfun %in% 
                              as.matrix(kegg_count$speciesfun)),]
kegg_count_fun$enrichment_anova = FALSE
kegg_count_fun$enrichment_anova[which(kegg_count_fun$speciesfun %in%
                                    as.matrix(kegg2_count$speciesfun))]=TRUE
allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.all.gene.txt',header=T)
allspecies=allspecies[,c(2,3,7,9,10)]
change_name = read.delim('final_SNP/annotation/all.denovo.gene.faa.changename.txt',header=F)
change_name=change_name[which(change_name$V4 %in%
                                kegg_count_fun$gene_name),]
allspecies$V5 = paste(allspecies$X.donor_species,allspecies$gene)
change_name$V5=paste(change_name$V1,change_name$V3)
change_name=merge(change_name,allspecies,by='V5',all.x=T)
change_name = change_name[,c(5,8:10)]
kegg_count_fun=merge(kegg_count_fun,change_name,by.x='gene_name',
                     by.y='V4',all=T)
kegg_count_fun=kegg_count_fun[,c(1,2,4,6:10,16:19)]
kegg_count_fun=kegg_count_fun[order(kegg_count_fun$enrichment_anova),]
kegg_count_fun=kegg_count_fun[order(kegg_count_fun$COG),]
kegg_count_fun=kegg_count_fun[order(kegg_count_fun$COG1),]
kegg_count_fun=kegg_count_fun[order(kegg_count_fun$species),]
kegg_count_fun$gene_name=str_split_fixed(kegg_count_fun$gene_name, "donor.", 2)[,2]
kegg_count_fun$gene_name=str_split_fixed(kegg_count_fun$gene_name, "__", 2)[,1]
colnames(kegg_count_fun)[c(1,2)]=
  c('donor','gene_cluster_num')
colnames(kegg_count_fun)[c(11,12)]=
  c('No.N_SNPs','No.S_SNPs')
colnames(kegg_count_fun)[c(9)]=
  c('enrichment_anova_test')
write.table(kegg_count_fun[,c(3,1,7,8,5,4,2,9:12)],'final_SNP/annotation/all.parallel.evolution.multiple.function.withdonor.txt',
            quote=F,sep='\t',row.names=F)
kegg_count_fun_count = count(kegg_count_fun[which(kegg_count_fun$enrichment_anova_test==TRUE),],c('species','COG2'))
write.table(kegg_count_fun_count,'final_SNP/annotation/all.parallel.evolution.multiple.function.withdonor.sum.txt',
            quote=F,sep='\t',row.names=F)
# split into within donor, across donor
kegg = read.delim('final_SNP/annotation/all.parallel.evolution.multiple.function.withdonor2.txt',header=T)
kegg_donor =count(kegg,c('species','COG2'))
kegg_donor=kegg_donor[which(kegg_donor$freq>1),]
kegg$speciesfun = paste(kegg$species,kegg$COG2)
kegg_donor$speciesfun = paste(kegg_donor$species,kegg_donor$COG2)
kegg=kegg[which(kegg$speciesfun %in% kegg_donor$speciesfun),]
kegg$species_gene = paste(kegg$species,
                                kegg$gene_cluster_num)
kegg$donor_species_gene = paste(paste(kegg$donor,kegg$species),
                          kegg$gene_cluster_num)

kegg2=kegg[!duplicated(paste(kegg$donor_species_gene)),]
kegg_donor =count(kegg2,c('species','gene_cluster_num'))
kegg_donor$species_gene = paste(kegg_donor$species,
                                      kegg_donor$gene_cluster_num)

kegg=merge(kegg,kegg_donor,by='species_gene')
write.table(kegg,'final_SNP/annotation/all.parallel.evolution.multiple.function.withdonor2.txt',
            quote=F,sep='\t',row.names=F)

#KEGG
kegg = read.delim('final_SNP/annotation/all.selected.gene.faa.High_select2.faa.cluster.aa.all.kegg.sum',header=T)
kegg$Tag = 'High_select'
kegg$speciesnew = kegg$species
kegg$speciesnew = str_split_fixed(kegg$speciesnew, "_", 2)[,1]
kegg$species = kegg$speciesnew
kegg=kegg[which(kegg$BRITE_KO3!=''),]
kegg=kegg[which(kegg$BRITE_KO3!='None'),]
kegg_count = count(
  kegg,c('species','BRITE_KO1','BRITE_KO2','BRITE_KO3') 
)
kegg_count=kegg_count[which(kegg_count$freq>1),]
kegg_count2= kegg_count[which(!(kegg_count$species %in% c('BA','BL'))),]
#
#----species fun eggnog figure 4----
Not_qualify = c()#c('BaVu','BaTh')#'ClBe'
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(plyr)
# start plotting
tomatrix_norm <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=mat[rownum,colnum] + dm[i,colnom]
  }
  return(mat)
}
tomatrix_text <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix('',nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=toString(dm[i,colnom])
  }
  return(mat)
}
tomatrix_text_single <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix('',nrow=length(allrow),ncol=1)
  row.names(mat)=allrow
  colnames(mat)='alltext'
  for (i in 1:nrow(dm))
  {
    if(dm[i,colnom]!='')
    {rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = 1
    mat[rownum,colnum]=toString(dm[i,colnom])}
  }
  return(mat)
}
add_extra_colomn <-function(mat)
{
  newmat = matrix(0,nrow = nrow(mat),ncol = 2*ncol(mat))
  row.names(newmat)=row.names(mat)
  colname_temp=c()
  for(i in 1:ncol(mat))
  {
    newmat[,i*2-1]=mat[,i]
    colname_temp = c(colname_temp,colnames(mat)[i],'')
  }
  colnames(newmat)=colname_temp
  return(newmat)
}

# within donor
kegg2=read.delim('final_SNP/annotation/all.parallel.evolution.multiple.function.withdonor.withindonor.txt',header=T)
kegg=kegg2
kegg$newgene = paste(kegg$COG1,kegg$COG2)
kegg$newgene2 = ''
kegg_count = count(kegg,c('newgene','species.x'))
for(i in 1:nrow(kegg_count))
{
  m = which(kegg$newgene == kegg_count$newgene[i] &
              kegg$species.x == kegg_count$species.x[i] )
  temp_name = c(1:kegg_count$freq[i])
  kegg$newgene2[m]=paste(kegg_count$species.x[i],temp_name)
}
kegg=kegg[order(kegg$newgene2),]
kegg=kegg[order(kegg$newgene),]
kegg2_mat = tomatrix_norm(kegg,16,15,12)
kegg2_mat=kegg2_mat[,order(colnames(kegg2_mat))]
kegg2_mat=kegg2_mat[order(row.names(kegg2_mat)),]
mat_breaks =c(seq(-0.5,4.5),max(kegg2_mat))
Metadata2=data.frame(as.matrix(kegg$COG1)[!duplicated(kegg$newgene),])
row.names(Metadata2)=colnames(kegg2_mat)
colnames(Metadata2) ='Annotation'
mat_colors <- list(Annotation = brewer.pal(5, "RdBu")[c(5,4,1)])
names(mat_colors$Annotation) <- unique(Metadata2$Annotation)
kegg2_matnew = add_extra_colomn(kegg2_mat)
pplot <- pheatmap(mat= kegg2_matnew,
                  color = c('#fcfbfd',
                            colorRampPalette(brewer.pal(n = 9, name = "Purples"))(
                              9)[c(-1,-2,-3)]),
                  breaks = mat_breaks,
                  annotation_col    = Metadata2,
                  annotation_colors = mat_colors,
                  show_rownames= TRUE,
                  show_colnames = TRUE,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE)

kegg2_mat = tomatrix_norm(kegg,16,15,13)
kegg2_mat=kegg2_mat[,order(colnames(kegg2_mat))]
kegg2_mat=kegg2_mat[order(row.names(kegg2_mat)),]

kegg2_matnew = add_extra_colomn(kegg2_mat)
kegg2_matnew = add_extra_colomn(kegg2_matnew)

mat_breaks =c(seq(-0.5,5.5),max(kegg2_mat))
pplot <- pheatmap(mat= kegg2_matnew,
                  color = c('#fcfbfd',
                            colorRampPalette(brewer.pal(n = 9, name = "Oranges"))(length(mat_breaks))[c(-1,-2)]),
                  breaks = mat_breaks,
                  #annotation_col    = Metadata2,
                  #annotation_colors = mat_colors,
                  show_rownames= TRUE,
                  show_colnames = TRUE,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE)

kegg2_mat = tomatrix_text(kegg,16,15,7)
kegg2_mat=kegg2_mat[,order(colnames(kegg2_mat))]
kegg2_mat=kegg2_mat[order(row.names(kegg2_mat)),]

ncol(kegg2_mat)
grid.table(kegg2_mat)

# across donor
kegg2=read.delim('final_SNP/annotation/all.parallel.evolution.multiple.function.withdonor.acrossdonor.txt',header=T)
kegg=kegg2
kegg$newgene = paste(kegg$COG1,kegg$COG2)
kegg$newgene = paste(kegg$newgene,kegg$EGGNOG)
kegg$newgene2 = ''
kegg_count = count(kegg,c('newgene','species.x'))
for(i in 1:nrow(kegg_count))
{
  m = which(kegg$newgene == kegg_count$newgene[i] &
              kegg$species.x == kegg_count$species.x[i] )
  temp_name = c(1:kegg_count$freq[i])
  temp_name2 = kegg_count$species.x[i]
  kegg$newgene2[m]=paste(temp_name2,temp_name)
}
kegg=kegg[order(kegg$newgene2),]
kegg=kegg[order(kegg$newgene),]
kegg2_mat = tomatrix_norm(kegg,16,15,12)
kegg2_mat=kegg2_mat[,order(colnames(kegg2_mat))]
kegg2_mat=kegg2_mat[order(row.names(kegg2_mat)),]
mat_breaks =c(seq(-0.5,4.5),max(kegg2_mat))
Metadata2=data.frame(as.matrix(kegg$COG1)[!duplicated(kegg$newgene),])
row.names(Metadata2)=colnames(kegg2_mat)
colnames(Metadata2) ='Annotation'
mat_colors <- list(Annotation = brewer.pal(5, "RdBu")[c(5,4,1)])
names(mat_colors$Annotation) <- unique(Metadata2$Annotation)
kegg2_matnew = kegg2_mat
pplot <- pheatmap(mat= kegg2_matnew,
                  color = c('#fcfbfd',
                            colorRampPalette(brewer.pal(n = 9, name = "Purples"))(
                              9)[c(-1,-2,-3)]),
                  breaks = mat_breaks,
                  annotation_col    = Metadata2,
                  annotation_colors = mat_colors,
                  show_rownames= TRUE,
                  show_colnames = TRUE,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE)

kegg$newgene = paste(kegg$newgene,kegg$annotation_short)
kegg2_mat = tomatrix_norm(kegg,16,15,13)
kegg2_mat=kegg2_mat[,order(colnames(kegg2_mat))]
kegg2_mat=kegg2_mat[order(row.names(kegg2_mat)),]

kegg2_matnew = add_extra_colomn(kegg2_mat)
#kegg2_matnew = add_extra_colomn(kegg2_matnew)

mat_breaks =c(seq(-0.5,5.5),max(kegg2_mat))
pplot <- pheatmap(mat= kegg2_matnew,
                  color = c('#fcfbfd',
                            colorRampPalette(brewer.pal(n = 9, name = "Oranges"))(length(mat_breaks))[c(-1,-2)]),
                  breaks = mat_breaks,
                  #annotation_col    = Metadata2,
                  #annotation_colors = mat_colors,
                  show_rownames= TRUE,
                  show_colnames = TRUE,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE)

# finished


#---- species fun fisher test figure 4 prep----
tomatrix <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=mat[rownum,colnum] + 1/dm[i,colnom]
  }
  return(mat)
}
tomatrix_norm <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=mat[rownum,colnum] + dm[i,colnom]
  }
  return(mat)
}
tomatrix_norm_test <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  dm[which(dm[,colnom]=='TRUE'),colnom]=1
  dm[which(dm[,colnom]=='FALSE'),colnom]=0
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=as.numeric(dm[i,colnom])
  }
  return(mat)
}
# prep
kegg2_count=read.delim('final_SNP/annotation/all.HS.eggnog.norm.species.test.txt',header=T)
kegg2_count = kegg2_count[which(kegg2_count$cross_genus!='True'),]
kegg2_count$anova_test[which(as.numeric(as.matrix(kegg2_count$diff))<1)]='FALSE'
kegg2_count=kegg2_count#[which(kegg2_count$anova_test=='TRUE'),]
kegg = read.delim('final_SNP/annotation/all.selected.gene.faa.High_select2.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
kegg$Tag = 'High_select'
kegg$speciesnew = kegg$species
kegg$speciesnew = str_split_fixed(kegg$speciesnew, "_", 2)[,1]
kegg$speciesnew[which(kegg$speciesnew=='95')]='EsCo'
kegg$speciesnew[which(kegg$speciesnew=='Bfragilis')]='BaFr'
kegg$species = kegg$speciesnew
kegg=kegg[which(kegg$COG!=''),]
kegg=kegg[which(kegg$COG1!='POORLY CHARACTERIZED'),]
kegg_count = count(
  kegg,c('species','COG','COG1','COG2') 
)
kegg_count$speciesfun = paste(kegg_count$species,kegg_count$COG)
kegg2_count$speciesfun = paste(kegg2_count$species,kegg2_count$COG)
kegg_count=kegg_count#[which(kegg_count$freq>1),]
unique(kegg_count$species)
kegg$speciesfun = paste(kegg$species,kegg$COG)
kegg_count_fun = kegg[which(kegg$speciesfun %in% 
                              as.matrix(kegg_count$speciesfun)),]
kegg_count_fun$enrichment_anova = FALSE
kegg2_count=kegg2_count#[which(kegg2_count$anova_test=='TRUE'),]
kegg_count_fun$enrichment_anova[which(kegg_count_fun$speciesfun %in%
                                        as.matrix(kegg2_count$speciesfun))]=TRUE
allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.all.gene.txt',header=T)
allspecies=allspecies[,c(2,3,7,9,10)]
change_name = read.delim('final_SNP/annotation/all.denovo.gene.faa.changename.txt',header=F)
change_name=change_name[which(change_name$V4 %in%
                                kegg_count_fun$gene_name),]
allspecies$V5 = paste(allspecies$X.donor_species,allspecies$gene)
change_name$V5=paste(change_name$V1,change_name$V3)
change_name=merge(change_name,allspecies,by='V5',all.x=T)
change_name = change_name[,c(5,8:10)]
kegg_count_fun=merge(kegg_count_fun,change_name,by.x='gene_name',
                     by.y='V4',all=T)
kegg_count_fun=kegg_count_fun[,c(1,2,4,6:10,16:19)]
kegg_count_fun=kegg_count_fun[order(kegg_count_fun$enrichment_anova),]
kegg_count_fun=kegg_count_fun[order(kegg_count_fun$COG),]
kegg_count_fun=kegg_count_fun[order(kegg_count_fun$COG1),]
kegg_count_fun=kegg_count_fun[order(kegg_count_fun$species),]
kegg_count_fun$gene_name=str_split_fixed(kegg_count_fun$gene_name, "donor.", 2)[,2]
kegg_count_fun$gene_name=str_split_fixed(kegg_count_fun$gene_name, "__", 2)[,1]
colnames(kegg_count_fun)[c(1,2)]=
  c('donor','gene_cluster_num')
colnames(kegg_count_fun)[c(11,12)]=
  c('No.N_SNPs','No.S_SNPs')
colnames(kegg_count_fun)[c(9)]=
  c('enrichment_anova_test')
write.table(kegg_count_fun[,c(3,1,7,8,5,4,2,9:12)],'final_SNP/annotation/all.parallel.evolution.all.function.withdonor.txt',
            quote=F,sep='\t',row.names=F)

# split into within donor, across donor
kegg = read.delim('final_SNP/annotation/all.parallel.evolution.all.function.withdonor.txt',header=T)
kegg_donor =count(kegg,c('species','COG2'))
kegg$speciesfun = paste(kegg$species,kegg$COG2)
kegg_donor$speciesfun = paste(kegg_donor$species,kegg_donor$COG2)
kegg$species_gene = paste(kegg$species,
                          kegg$gene_cluster_num)
kegg$donor_species_gene = paste(paste(kegg$donor,kegg$species),
                                kegg$gene_cluster_num)

kegg2=kegg[!duplicated(paste(kegg$donor_species_gene)),]
kegg_donor =count(kegg2,c('species','gene_cluster_num'))
kegg_donor$species_gene = paste(kegg_donor$species,
                                kegg_donor$gene_cluster_num)

kegg=merge(kegg,kegg_donor,by='species_gene')
write.table(kegg,'final_SNP/annotation/all.parallel.evolution.all.function.withdonor2.txt',
            quote=F,sep='\t',row.names=F)

#----simulation fisher test figure 4 prep----
#Not_qualify = c('BaVu','BaTh')#'ClBe'
Not_qualify = c()#'ClBe'
library(pheatmap)
tomatrix <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=mat[rownum,colnum] + 1/dm[i,colnom]
  }
  return(mat)
}
tomatrix_mat <- function(dm,colnum1,colnum2)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=mat[rownum,colnum] + 1
  }
  return(mat)
}
#fisher test
fisher_test_mat <- function(mat)
{
  result_test=matrix(0,nrow=nrow(mat),ncol=ncol(mat))
  colnames(result_test)=colnames(mat)
  row.names(result_test) = row.names(mat)
  #test
  k=0
  Total = sum(mat)
  for (i in 1:nrow(mat))
  {
    for (j in 1:ncol(mat))
    {
      k=k+1
      temp=matrix(0,nrow=2,ncol=2)
      temp[1,1]=as.integer(mat[i,j])
      temp[1,2]=as.integer(sum(mat[i,])-temp[1,1])
      temp[2,1]=as.integer(sum(mat[,j])-temp[1,1])
      temp[2,2]=as.integer(Total -sum(mat[i,])-sum(mat[,j])+ temp[1,1])
      result2 = fisher.test(temp,alternative='greater',
                            simulate.p.value =TRUE,B = 1000,
                            conf.int = TRUE, conf.level = 0.95)
      result_test[i,j]=result2$p.value
    }
  }
  return(result_test)
}
#pheatmap
fisher_heatmap <- function(mat,p = 0.01)
{
  mat=mat[,order(colnames(mat))]
  mat=mat[order(row.names(mat)),]
  mat[which(mat <= p)]=0
  mat[which(mat > p)] = 1
  pheatmap(mat= mat,
           color = c(colorRampPalette(brewer.pal(n = 5, name = "PRGn"))(
             5)[c(1,3)]),
           legend_breaks = c(0,1), 
           legend_labels = c('sig','non_sig'),
           show_rownames= TRUE,
           show_colnames = TRUE,
           cluster_cols = FALSE,
           cluster_rows = FALSE)
}
heatmap_mat <- function(mat)
{
  mat=mat[,order(colnames(mat))]
  mat=mat[order(row.names(mat)),]
  
  mat[which(mat >= 5)] = 5
  pheatmap(mat= mat,
           color = c(colorRampPalette(brewer.pal(n = 9, name = "Purples"))(
             9)[c(1,2,5:9)]),
           show_rownames= TRUE,
           show_colnames = TRUE,
           cluster_cols = FALSE,
           cluster_rows = FALSE)
}
heatmap_mat_percentile <- function(mat)
{
  mat=mat[,order(colnames(mat))]
  mat=mat[order(row.names(mat)),]
  pheatmap(mat= mat,
           color = c(colorRampPalette(brewer.pal(n = 9, name = "PRGn"))(
             9)[c(1:9)]),
           breaks = seq(0.1,1,0.1),
           show_rownames= TRUE,
           show_colnames = TRUE,
           cluster_cols = FALSE,
           cluster_rows = FALSE)
}
tomatrix <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=mat[rownum,colnum] + 1/dm[i,colnom]
  }
  return(mat)
}
tomatrix_norm <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=mat[rownum,colnum] + dm[i,colnom]
  }
  return(mat)
}
tomatrix_norm_test <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  dm[which(dm[,colnom]=='TRUE'),colnom]=1
  dm[which(dm[,colnom]=='FALSE'),colnom]=0
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=as.numeric(dm[i,colnom])
  }
  return(mat)
}

# HS load
keggHS = read.delim('final_SNP/annotation/all.parallel.evolution.all.function.withdonor2.txt',header=T)
temp_simulate_sum3=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted2.txt',header=T)
temp_simulate_sum3 = temp_simulate_sum3[!duplicated(temp_simulate_sum3$lineage),]
temp_simulate_sum3=temp_simulate_sum3[which((temp_simulate_sum3$species %in% as.matrix(keggHS$species.x))),]
species_HS = data.frame(species=unique(temp_simulate_sum3$species),
                        HS_lineage = 0
                        )
for(i in 1:nrow(species_HS))
{
  temp_all = temp_simulate_sum3[which(temp_simulate_sum3$species==
                                        species_HS$species[i]),]
  species_HS$HS_lineage[i]=sum(temp_all$HS_lineage)
}
# all genes in genomes
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(stringr)
keggall = read.delim('final_SNP/annotation/all.genome.gene.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
keggall$Tag = 'All'
keggall$speciesnew = keggall$species
keggall$speciesnew = str_split_fixed(keggall$speciesnew, "_", 2)[,1]
keggall$speciesnew[which(keggall$speciesnew=='95')]='EsCo'
keggall$speciesnew[which(keggall$speciesnew=='Bfragilis')]='BaFr'
keggall$species = keggall$speciesnew
keggall=keggall[which(keggall$COG!=''),]
keggall=keggall[which(keggall$COG1!='POORLY CHARACTERIZED'),]
keggall$species_genus = paste(keggall$species,keggall$cross_genus)

keggall$fun = paste(keggall$COG1,keggall$COG2)
keggall_sum = tomatrix_mat(keggall,3,16)
keggall_sum2 = tomatrix_mat(keggall,15,16)
test1 = fisher_test_mat(keggall_sum)
test2 = fisher_test_mat(keggall_sum2)
test1=test1[row.names(test1) %in% species_HS$species,]

fisher_heatmap(test1)
write.table(keggall_sum,'final_SNP/annotation/all.genome.eggnog.mat.txt',
            quote=F,sep='\t',row.names=T)
write.table(keggall_sum2,'final_SNP/annotation/all.genome.eggnog.mat.withcore.txt',
            quote=F,sep='\t',row.names=T)
write.table(test1,'final_SNP/annotation/all.genome.eggnog.mat.fisher.txt',
            quote=F,sep='\t',row.names=T)
write.table(test2,'final_SNP/annotation/all.genome.eggnog.mat.withcore.fisher.txt',
            quote=F,sep='\t',row.names=T)

# Denovo
kegg = read.delim('final_SNP/annotation/all.denovo.gene.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
kegg$Tag = 'Denovo'
kegg$speciesnew = kegg$species
kegg$speciesnew = str_split_fixed(kegg$speciesnew, "_", 2)[,1]
kegg$speciesnew[which(kegg$speciesnew=='95')]='EsCo'
kegg$speciesnew[which(kegg$speciesnew=='Bfragilis')]='BaFr'
kegg$species = kegg$speciesnew
kegg=kegg[which(kegg$COG!=''),]
kegg=kegg[which(kegg$COG1!='POORLY CHARACTERIZED'),]
kegg$species_genus = paste(kegg$species,kegg$cross_genus)
kegg = kegg[which(kegg$species %in% species_HS$species),]
kegg$fun = paste(kegg$COG1,kegg$COG2)
kegg_sum = tomatrix_mat(kegg,3,16)
kegg_sum2 = tomatrix_mat(kegg,15,16)
test1_denovo = fisher_test_mat(kegg_sum)
test2_denovo = fisher_test_mat(kegg_sum2)
fisher_heatmap(test1_denovo, 0.05)
write.table(kegg_sum,'final_SNP/annotation/all.denovo.eggnog.mat.txt',
            quote=F,sep='\t',row.names=T)
write.table(kegg_sum2,'final_SNP/annotation/all.denovo.eggnog.mat.withcore.txt',
            quote=F,sep='\t',row.names=T)
write.table(test1_denovo,'final_SNP/annotation/all.denovo.eggnog.mat.fisher.txt',
            quote=F,sep='\t',row.names=T)
write.table(test2_denovo,'final_SNP/annotation/all.denovo.eggnog.mat.withcore.fisher.txt',
            quote=F,sep='\t',row.names=T)

test1_denovo2 = matrix(
  p.adjust(test1_denovo, method = 'fdr'),
  ncol = ncol(test1_denovo))
row.names(test1_denovo2)=row.names(test1_denovo)
colnames(test1_denovo2)=colnames(test1_denovo)
fisher_heatmap(test1_denovo2, 0.05)
write.table(test1_denovo2,'final_SNP/annotation/all.denovo.eggnog.mat.fisher.fdr.txt',
            quote=F,sep='\t',row.names=T)

# HS simulation
kegg2new=kegg
keggsim_sum = kegg_sum
keggsim_sum2 =kegg_sum2
temp_list = list()
keggsim_sum[] = 0
keggsim_sum2[] = 0
merge_mat <- function(keggsim_sum,mat)
{
  for(i in 1:nrow(mat))
  {
    for(j in 1:ncol(mat))
    {
      keggsim_sum[which(row.names(keggsim_sum)==row.names(mat)[i]),
                  which(colnames(keggsim_sum)==colnames(mat)[j])
                  ] = keggsim_sum[which(row.names(keggsim_sum)==row.names(mat)[i]),
                                  which(colnames(keggsim_sum)==colnames(mat)[j])
                                  ] + mat[i,j]
    }
  }
  return(keggsim_sum)
}
format_mat <- function(keggsim_sum,mat)
{
  keggsim_sum[]=0
  for(i in 1:nrow(mat))
  {
    for(j in 1:ncol(mat))
    {
      keggsim_sum[which(row.names(keggsim_sum)==row.names(mat)[i]),
                  which(colnames(keggsim_sum)==colnames(mat)[j])
                  ] = keggsim_sum[which(row.names(keggsim_sum)==row.names(mat)[i]),
                                  which(colnames(keggsim_sum)==colnames(mat)[j])
                                  ] + mat[i,j]
    }
  }
  return(keggsim_sum)
}
# simulate against all genes
for(k in 1:100)
{
  kegg2sim=data.frame()
  for(m in 1:nrow(species_HS))
  {
    species=toString(species_HS$species[m])
    HS_gene = species_HS$HS_lineage[m]
    temp_kegg = kegg2new[which(kegg2new$species==species),]
    temp_kegg2 = temp_kegg[!duplicated(temp_kegg$gene_name),]
    subset_gene = temp_kegg2$gene_name[sample(1:nrow(temp_kegg2), HS_gene,replace = TRUE)]
    kegg2_control=as.matrix(temp_kegg[which(temp_kegg$gene_name %in% subset_gene),])
    kegg2sim=rbind(kegg2sim,kegg2_control)
  }
  kegg2sim=data.frame(kegg2sim)
  keggsim_sum = merge_mat(keggsim_sum,tomatrix_mat(kegg2sim,3,16))
  temp_list[[k]]=format_mat(keggsim_sum,tomatrix_mat(kegg2sim,3,16))
  keggsim_sum2 = merge_mat(keggsim_sum2,tomatrix_mat(kegg2sim,15,16))
}
keggsim_sum = keggsim_sum/100
keggsim_sum2 = keggsim_sum2/100
test1_sim = fisher_test_mat(keggsim_sum)
test2_sim = fisher_test_mat(keggsim_sum2)
heatmap_mat(keggsim_sum)
write.table(keggsim_sum,'final_SNP/annotation/all.HS.sim.eggnog.mat.txt',
            quote=F,sep='\t',row.names=T)
write.table(keggsim_sum2,'final_SNP/annotation/all.HS.sim.eggnog.mat.withcore.txt',
            quote=F,sep='\t',row.names=T)
write.table(test1_sim,'final_SNP/annotation/all.HS.sim.eggnog.mat.fisher.txt',
            quote=F,sep='\t',row.names=T)
write.table(test2_sim,'final_SNP/annotation/all.HS.sim.eggnog.mat.withcore.fisher.txt',
            quote=F,sep='\t',row.names=T)

# real HS load
keggHS = read.delim('final_SNP/annotation/all.parallel.evolution.all.function.withdonor2.txt',header=T)
CI_eva <- function(dis, value)
{
  dis = as.numeric(dis)
  return(length(which(value>dis))/length(dis)+
           length(which(value==dis))/length(dis)/2)
}
# single donor
keggHS = keggHS[which(keggHS$freq == 1),]
kegg2_count2=keggHS
kegg2_count2=kegg2_count2[order(kegg2_count2$COG2),]
kegg2_count2=kegg2_count2[order(kegg2_count2$COG1),]
kegg2_count2$COG3_new = paste(kegg2_count2$COG1,kegg2_count2$COG2)
kegg2_count2=kegg2_count2[order(kegg2_count2$COG3_new),]
kegg2_count2=kegg2_count2[!duplicated(paste(kegg2_count2$donor_species_gene,
                                            kegg2_count2$COG3_new)),]
keggHS_mat = tomatrix_norm(kegg2_count2,2,18,17)
heatmap_mat(keggHS_mat)
keggHS_mat_percentile = keggHS_mat
keggHS_mat_percentile[]=0
#calculate percentile
for(i in 1:nrow(keggHS_mat))
  for(j in 1:ncol(keggHS_mat))
{
    observe_num = keggHS_mat[i,j]
    expect_num_set = c()
    for(k in 1:100)
{
      if(k == 1)
      {
        colnum = which(colnames(temp_mat)==colnames(keggHS_mat)[j])
        rownum = which(row.names(temp_mat)==row.names(keggHS_mat)[i])
      }
  temp_mat = temp_list[[k]]
  expect_num = temp_mat[rownum,
                        colnum
                        ]
  expect_num_set=c(expect_num_set,expect_num)
    }
    keggHS_mat_percentile[i,j]=CI_eva(expect_num_set,observe_num)
  }
heatmap_mat_percentile(keggHS_mat_percentile)
write.table(keggHS_mat_percentile,'final_SNP/annotation/all.HS.sim.eggnog.mat.quantile.txt',
            quote=F,sep='\t',row.names=T)

#
#----simulation eggnog figure 3 prep----
#Not_qualify = c('BaVu','BaTh')#'ClBe'
Not_qualify = c()#'ClBe'

tomatrix <- function(dm,colnum1,colnum2,colnom)
{
  allrow = unique(dm[,colnum1])
  allcol = unique(dm[,colnum2])
  mat=matrix(0,nrow=length(allrow),ncol=length(allcol))
  row.names(mat)=allrow
  colnames(mat)=allcol
  for (i in 1:nrow(dm))
  {
    rownum = which(row.names(mat)==dm[i,colnum1])
    colnum = which(colnames(mat)==dm[i,colnum2])
    mat[rownum,colnum]=mat[rownum,colnum] + 1/dm[i,colnom]
  }
  return(mat)
}
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(stringr)
keggall = read.delim('final_SNP/annotation/all.genome.gene.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
keggall$Tag = 'All'
keggall$speciesnew = keggall$species
keggall$speciesnew = str_split_fixed(keggall$speciesnew, "_", 2)[,1]
keggall$speciesnew[which(keggall$speciesnew=='95')]='EsCo'
keggall$speciesnew[which(keggall$speciesnew=='Bfragilis')]='BaFr'
keggall$species = keggall$speciesnew
keggall=keggall[which(keggall$COG!=''),]
keggall=keggall[which(keggall$COG1!='POORLY CHARACTERIZED'),]
keggall$species_genus = paste(keggall$species,keggall$cross_genus)

# HS genes
kegg = read.delim('final_SNP/annotation/all.selected.gene.faa.High_select2.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
kegg$Tag = 'High_select'
kegg$speciesnew = kegg$species
kegg$speciesnew = str_split_fixed(kegg$speciesnew, "_", 2)[,1]
kegg$speciesnew[which(kegg$speciesnew=='95')]='EsCo'
kegg$speciesnew[which(kegg$speciesnew=='Bfragilis')]='BaFr'
kegg$species = kegg$speciesnew
kegg=kegg[which(kegg$COG!=''),]
kegg=kegg[which(kegg$COG1!='POORLY CHARACTERIZED'),]
kegg$species_genus = paste(kegg$species,kegg$cross_genus)
kegg_nom=read.delim('final_SNP/annotation/all.HS.eggnog.norm.species.test.txt',header=T)
allpop = read.delim('final_SNP/HSgene/all.species.all.species.lineagenum.txt')
allpop=allpop[which(!(allpop$species %in% Not_qualify)),]
allpop=allpop[which((allpop$species %in% c(as.matrix(kegg$species),
                                           'BaSa'))),]
# COG1 level
# HS genes
kegg2=kegg
kegg4=count(kegg2,c('species'))
kegg2=merge(kegg2,kegg4,by='species')
kegg_mat = tomatrix(kegg2,15,9,16)
row.names(kegg_mat)=paste(row.names(kegg_mat),'High_select')
unique(kegg2$species)
# control simulation
kegg2new=keggall
kegg4=count(kegg2new,c('species'))
kegg2new=merge(kegg2new,kegg4,by='species')
kegg2new=kegg2new[which(kegg2new$species %in% unique(kegg2$species)),]
species_HS = allpop[,c(1,4)]
# simulate against all genes
for(k in 1:100)
{
  kegg2sim=data.frame()
  for(m in 1:nrow(species_HS))
  {
    species=toString(species_HS$species[m])
    HS_gene = species_HS$HS_NO[m]
    temp_kegg = kegg2new[which(kegg2new$species==species),]
    temp_kegg2 = temp_kegg[!duplicated(temp_kegg$gene_name),]
    subset_gene = temp_kegg2$gene_name[sample(1:nrow(temp_kegg2), HS_gene,replace = TRUE)]
    kegg2_control=as.matrix(temp_kegg[which(temp_kegg$gene_name %in% subset_gene),])
    kegg2_control[,1]=toString(species)
    kegg2_control[,16] = nrow(kegg2_control)
    kegg2sim=rbind(kegg2sim,kegg2_control)
  }
  kegg2sim$freq=as.numeric(as.matrix(kegg2sim$freq))
  kegg_mat2 = tomatrix(kegg2sim,15,9,16)
  row.names(kegg_mat2)=paste(row.names(kegg_mat2),'simulation HS')
  kegg_mat=rbind(kegg_mat,kegg_mat2)
}
write.table(kegg_mat,'final_SNP/annotation/all.HS.eggnog.all.simulation.txt',
            quote=F,sep='\t',row.names=T)

# denovo genes
kegg = read.delim('final_SNP/annotation/all.denovo.gene.faa.cluster.aa.all.eggnog.sum.species.sum',header=T)
kegg$Tag = 'Denovo'
kegg$speciesnew = kegg$species
kegg$speciesnew = str_split_fixed(kegg$speciesnew, "_", 2)[,1]
kegg$speciesnew[which(kegg$speciesnew=='95')]='EsCo'
kegg$speciesnew[which(kegg$speciesnew=='Bfragilis')]='BaFr'
kegg$species = kegg$speciesnew
kegg=kegg[which(kegg$COG!=''),]
kegg=kegg[which(kegg$COG1!='POORLY CHARACTERIZED'),]
kegg$species_genus = paste(kegg$species,kegg$cross_genus)
kegg2=kegg
kegg4=count(kegg2,c('species'))
kegg2=merge(kegg2,kegg4,by='species')
kegg_mat2 = tomatrix(kegg2,15,9,16)
row.names(kegg_mat2)=paste(row.names(kegg_mat2),'Denovo')
kegg_mat=kegg_mat2#rbind(kegg_mat,kegg_mat2)
unique(kegg2$species)
# control simulation
allpop = read.delim('final_SNP/HSgene/all.species.all.species.lineagenum.txt')
allpop=allpop[which(!(allpop$species %in% Not_qualify)),]
allpop=allpop[which((allpop$species %in% as.matrix(kegg$species))),]
species_HS = allpop[which(allpop$HS_NO>0),c(1,12)]

kegg2new=keggall
kegg4=count(kegg2new,c('species'))
kegg2new=merge(kegg2new,kegg4,by='species')
kegg2new=kegg2new[which(kegg2new$species %in% unique(kegg2$species)),]

# simulate against all genes
for(k in 1:100)
{
  kegg2sim=data.frame()
  for(m in 1:nrow(species_HS))
  {
    species=toString(species_HS$species[m])
    HS_gene = species_HS$No.gene_SNPs[m]
    temp_kegg = kegg2new[which(kegg2new$species==species),]
    temp_kegg2 = temp_kegg[!duplicated(temp_kegg$gene_name),]
    subset_gene = temp_kegg2$gene_name[sample(1:nrow(temp_kegg2), HS_gene,replace = TRUE)]
    kegg2_control=as.matrix(temp_kegg[which(temp_kegg$gene_name %in% subset_gene),])
    kegg2_control[,1]=toString(species)
    kegg2_control[,16] = nrow(kegg2_control)
    kegg2sim=rbind(kegg2sim,kegg2_control)
  }
  kegg2sim$freq=as.numeric(as.matrix(kegg2sim$freq))
  kegg_mat2 = tomatrix(kegg2sim,15,9,16)
  row.names(kegg_mat2)=paste(row.names(kegg_mat2),'simulation Denovo')
  kegg_mat=rbind(kegg_mat,kegg_mat2)
}
write.table(kegg_mat,'final_SNP/annotation/all.HS.eggnog.all.simulation.2.txt',
            quote=F,sep='\t',row.names=T)

# form table
library(stringr)
CI_eva <- function(dis, value)
{
  dis = as.numeric(dis)
  return(length(which(value>dis))/length(dis)+
           length(which(value==dis))/length(dis)/2)
}
kegg_mat = read.delim('final_SNP/annotation/all.HS.eggnog.all.simulation.txt',header=F)
colnames(kegg_mat)[2:4]=as.matrix(kegg_mat[1,1:3])
colnames(kegg_mat)[1]='species_tag'
kegg_mat=kegg_mat[-1,]
kegg_mat2 = kegg_mat
kegg_mat = read.delim('final_SNP/annotation/all.HS.eggnog.all.simulation.2.txt',header=F)
colnames(kegg_mat)[2:4]=as.matrix(kegg_mat[1,1:3])
colnames(kegg_mat)[1]='species_tag'
kegg_mat=kegg_mat[-1,]
kegg_mat2=kegg_mat2[,c(1,which(colnames(kegg_mat2)==colnames(kegg_mat)[2]),
                       which(colnames(kegg_mat2)==colnames(kegg_mat)[3]),
                       which(colnames(kegg_mat2)==colnames(kegg_mat)[4])
                       )]
kegg_mat=rbind(as.matrix(kegg_mat),as.matrix(kegg_mat2))
write.table(kegg_mat,'final_SNP/annotation/all.HS.eggnog.all.simulation.all.txt',
            quote=F,sep='\t',row.names=F)
kegg_mat = read.delim('final_SNP/annotation/all.HS.eggnog.all.simulation.all.txt',header=T)

kegg_mat2=data.frame(species_tag = kegg_mat$species_tag,
                    value = kegg_mat$METABOLISM,
                    tag2 = 'M')
temp_mat = data.frame(species_tag = kegg_mat$species_tag,
                      value = kegg_mat$INFORMATION.STORAGE.AND.PROCESSING,
                      tag2 = 'I')
kegg_mat2$value=as.numeric(as.matrix(kegg_mat2$value))
kegg_mat2=rbind(kegg_mat2,temp_mat)
temp_mat = data.frame(species_tag = kegg_mat$species_tag,
                      value = kegg_mat$CELLULAR.PROCESSES.AND.SIGNALING,
                      tag2 = 'C')
kegg_mat2$value=as.numeric(as.matrix(kegg_mat2$value))
kegg_mat2=rbind(kegg_mat2,temp_mat)
kegg_mat=kegg_mat2
kegg_mat$species=str_split_fixed(kegg_mat$species_tag, " ", 2)[,1]
kegg_mat$genus=str_split_fixed(kegg_mat$species_tag, " ", 3)[,2]
kegg_mat$tag=str_split_fixed(kegg_mat$species_tag, " ", 3)[,3]
kegg_mat$label = kegg_mat$species
kegg_mat$label[which(kegg_mat$tag=='simulation HS')]=''
kegg_mat$label[which(kegg_mat$tag=='simulation Denovo')]=''
kegg_mat_sim_denovo = kegg_mat[which(kegg_mat$tag=='simulation Denovo'
                                     
                                     & kegg_mat$genus == 'False'),]
kegg_mat_sim_HS = kegg_mat[which(kegg_mat$tag=='simulation HS' 
                                 & kegg_mat$genus == 'False'),]
alltag = unique(kegg_mat$tag2)
kegg_mat_sim_mean = matrix(0,nrow = 6,ncol = 4)
colnames(kegg_mat_sim_mean)=c('tag2','tag','mean','sd')
k = 0
for(a_tag in alltag)
{
  k = k + 1
  kegg_mat_sim_mean[k,]=c(a_tag,'Genes of parallel evolution',
                          mean(kegg_mat_sim_HS$value[which(kegg_mat_sim_HS$tag2 == a_tag)]),
                          sd(kegg_mat_sim_HS$value[which(kegg_mat_sim_HS$tag2 == a_tag)])
                          )
  k = k + 1
  kegg_mat_sim_mean[k,]=c(a_tag,'All mutations',
                          mean(kegg_mat_sim_denovo$value[which(kegg_mat_sim_denovo$tag2 == a_tag)]),
                          sd(kegg_mat_sim_denovo$value[which(kegg_mat_sim_denovo$tag2 == a_tag)])
                          )
}
write.table(kegg_mat_sim_mean,'final_SNP/annotation/all.HS.eggnog.all.simulation.form.ref.txt',
            quote=F,sep='\t',row.names=F)

kegg_mat$tagall = paste(paste(kegg_mat$species,kegg_mat$tag2),kegg_mat$genus)
kegg_mat=kegg_mat[which(!(kegg_mat$species %in% Not_qualify)),]
# genes
kegg_pca3 = kegg_mat[which(kegg_mat$tag %in% c(
  'Denovo'
)),]
kegg_pca2 = kegg_mat[which(kegg_mat$tag %in% c(
  'High_select'
)),]
kegg_pca3 = merge(kegg_pca3,kegg_pca2[,c(2,8)],by='tagall',all=T)

# de novo simulation
kegg_pca1 = kegg_mat[which(kegg_mat$tag %in% c(
  'simulation Denovo'
)),]
kegg_pca2=kegg_pca3
kegg_pca2$quantile.x = 0
num_col = ncol(kegg_pca2)
for(i in 1:nrow(kegg_pca2))
{
  temp_kegg = kegg_pca1[which(kegg_pca1$tagall == kegg_pca2$tagall[i]),]
  kegg_pca2[i,num_col]=CI_eva(temp_kegg$value,kegg_pca2$value.x[i])
}
kegg_pca3 = kegg_pca2
# HS simulation
kegg_pca1 = kegg_mat[which(kegg_mat$tag %in% c(
  'simulation HS'
)),]
kegg_pca2=kegg_pca3
kegg_pca2$quantile.y = 0
num_col = ncol(kegg_pca2)
for(i in 1:nrow(kegg_pca2))
{
  temp_kegg = kegg_pca1[which(kegg_pca1$tagall == kegg_pca2$tagall[i]),]
  kegg_pca2[i,num_col]=CI_eva(temp_kegg$value,kegg_pca2$value.y[i])
}
allpop = read.delim('final_SNP/HSgene/all.species.all.species.lineagenum.txt')
allpop = allpop[which(allpop$HS_NO>0),]

kegg_pca3 = kegg_pca2[which(kegg_pca2$species %in% allpop$species),]
# output
write.table(kegg_pca3,'final_SNP/annotation/all.HS.eggnog.all.simulation.form.txt',
            quote=F,sep='\t',row.names=F)

#----plot simulation eggnog figure 3----
# plot each catogory!!! no PCA
# each species a p-value
#HS species cog1
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(stringr)

kegg_pca = read.delim('final_SNP/annotation/all.HS.eggnog.all.simulation.form.txt',header=T)
kegg_mat_sim_mean = read.delim('final_SNP/annotation/all.HS.eggnog.all.simulation.form.ref.txt',header=T)

temp_simulate_sum3=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted2.txt',header=T)
temp_simulate_sum3_sub = temp_simulate_sum3[which(!is.na(temp_simulate_sum3$HS_norm_gene) &
                                                    temp_simulate_sum3$HS_norm_gene >0  ),]
species_order = read.delim('species_order.new.fullname.txt',header=T)
#qualified species
qualified_species = unique(temp_simulate_sum3_sub$species)
kegg_pca=kegg_pca[which((kegg_pca$species %in% qualified_species
                           )),]
kegg_pca=merge(kegg_pca,species_order,by='species',all.x=T)
# start plot
kegg_pca=kegg_pca[order(kegg_pca$neworder),]
kegg_pca$value.y[is.na(kegg_pca$value.y)]=0
kegg_pca$quantile.y[is.na(kegg_pca$quantile.y)]=0
alltag = unique(kegg_pca$tag2)
i = 2 # 1-3
{Tag2 = alltag[i] 
kegg_pca1=kegg_pca[which(kegg_pca$tag2==Tag2
                         &kegg_pca$genus!='True'),]
kegg_pca1=kegg_pca1[rev(order(kegg_pca1$neworder)),]
kegg_pca1$neworder = 1:nrow(kegg_pca1)
kegg_pca2 = kegg_pca1
kegg_pca1 = data.frame(
  neworder = kegg_pca2$neworder,
  quantile.x =  kegg_pca2$quantile.x,
  tag = 'All mutations',
  tag2 = 'All mutations observe',
  label = kegg_pca2$label
)
kegg_pca1 = rbind(kegg_pca1,
                  data.frame(
                    neworder = kegg_pca2$neworder,
                    quantile.x =  kegg_pca2$quantile.y,
                    tag = 'Genes of parallel evolution',
                    tag2 = 'Genes of parallel evolution observe',
                    label = kegg_pca2$label
                  )
                  )
kegg_pca1 = rbind(kegg_pca1,
                  data.frame(
                    neworder = kegg_pca2$neworder,
                    quantile.x =  0.5,
                    tag = 'Genes of parallel evolution',
                    tag2 = 'Genes of parallel evolution',
                    label = kegg_pca2$label
                  )
)
kegg_pca1 = rbind(kegg_pca1,
                  data.frame(
                    neworder = kegg_pca2$neworder,
                    quantile.x =  0.5,
                    tag = 'All mutations',
                    tag2 = 'All mutations',
                    label = kegg_pca2$label
                  )
)
ggplot()+
  geom_point(kegg_pca1,mapping = aes(
    y=neworder,
    x=(quantile.x-0.5)*100,
    color=tag2,
    fill=tag2
  ),
  alpha = 1,
  size=4
  )+
  geom_line(kegg_pca1,mapping = aes(
    y=neworder,
    x=(quantile.x-0.5)*100,
    color=tag,
    fill=tag,
    group = paste(neworder,tag)
  ),
  alpha = 0.3,
  size=1
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  scale_fill_manual(values =c(brewer.pal(5, "PRGn")[c(5,1,2,4)]))+
  scale_color_manual(values =c(brewer.pal(5, "PRGn")[c(5,1,2,4)]))+
  scale_x_continuous(limit = c(-50,50))+
  scale_y_continuous(breaks = 1:nrow(kegg_pca1),
                     labels = kegg_pca1$label)+
  geom_vline(xintercept=0, color='grey', size=1)+
  xlab("percentile compared to the neutral model (%)") + 
  labs(title=Tag2)+
  scale_x_continuous(limits = c(0,50))
}
{Tag2 = alltag[i] 
  kegg_pca1=kegg_pca[which(kegg_pca$tag2==Tag2
                           &kegg_pca$genus!='True'),]
  kegg_pca1=kegg_pca1[rev(order(kegg_pca1$neworder)),]
  kegg_pca1$neworder = 1:nrow(kegg_pca1)
  kegg_pca2 = kegg_pca1
  ref_all = kegg_mat_sim_mean$mean[which(kegg_mat_sim_mean$tag2==Tag2
                                         &
                                           kegg_mat_sim_mean$tag=='All mutations')]
  ref_parallel = kegg_mat_sim_mean$mean[which(kegg_mat_sim_mean$tag2==Tag2
                                              &
                                                kegg_mat_sim_mean$tag=='Genes of parallel evolution')]
  
  kegg_pca1 = data.frame(
    neworder = kegg_pca2$neworder,
    value.x =  kegg_pca2$value.x,
    tag = 'All mutations',
    tag2 = 'All mutations observed',
    label = kegg_pca2$label
  )
  kegg_pca1 = rbind(kegg_pca1,
                    data.frame(
                      neworder = kegg_pca2$neworder,
                      value.x =  kegg_pca2$value.y,
                      tag = 'Genes of parallel evolution',
                      tag2 = 'Genes of parallel evolution observed',
                      label = kegg_pca2$label
                    )
  )
  kegg_pca1 = rbind(kegg_pca1,
                    data.frame(
                      neworder = kegg_pca2$neworder,
                      value.x =  ref_parallel,
                      tag = 'Genes of parallel evolution',
                      tag2 = 'Genes of parallel evolution',
                      label = kegg_pca2$label
                    )
  )
  kegg_pca1 = rbind(kegg_pca1,
                    data.frame(
                      neworder = kegg_pca2$neworder,
                      value.x =  ref_all,
                      tag = 'All mutations',
                      tag2 = 'All mutations',
                      label = kegg_pca2$label
                    )
  )
  ggplot()+
    geom_point(kegg_pca1,mapping = aes(
      y=neworder,
      x=value.x*100,
      color=tag2,
      fill=tag2
    ),
    alpha = 1,
    size=4
    )+
    geom_line(kegg_pca1,mapping = aes(
      y=neworder,
      x=value.x*100,
      color=tag,
      fill=tag,
      group = paste(neworder,tag)
    ),
    alpha = 0.3,
    size=1
    )+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1)
    ) +
    scale_fill_manual(values =c(brewer.pal(5, "PRGn")[c(5,1,2,4)]))+
    scale_color_manual(values =c(brewer.pal(5, "PRGn")[c(5,1,2,4)]))+
    scale_x_continuous(limit = c(-0,100))+
    scale_y_continuous(breaks = 1:nrow(kegg_pca1),
                       labels = kegg_pca1$label)+
    geom_vline(xintercept=ref_parallel*100, size=1,color = c(brewer.pal(5, "PRGn")[c(1)]))+
    geom_vline(xintercept=ref_all*100, size=1, color = c(brewer.pal(5, "PRGn")[c(5)]))+
    xlab("percentage compared to all functions (%)") + 
    labs(title=Tag2)
}
{Tag2 = alltag[i]
  kegg_pca1=kegg_pca[which(kegg_pca$tag2==Tag2
                           &kegg_pca$genus!='True'),]
  
  temp_result = c(toString(Tag2),
                  mean(kegg_pca1$value.x[which(kegg_pca1$tag2 == Tag2)]),
                  kegg_mat_sim_mean$mean[which(kegg_mat_sim_mean$tag2==Tag2
                                               &
                                                 kegg_mat_sim_mean$tag=='All mutations')],
                  mean(kegg_pca1$value.y[which(kegg_pca1$tag2 == Tag2)]),
                  kegg_mat_sim_mean$mean[which(kegg_mat_sim_mean$tag2==Tag2
                                               &
                                                 kegg_mat_sim_mean$tag=='Genes of parallel evolution')]
                  
  )
  print(c('Fun','mean_denovo','mean_simulation',
          'mean_HS','mean_simulation_HS'))
  print(temp_result)
}
#
#---- mutation freq----
Mut = read.delim('final_SNP/HSgene/all.species.High_select2.txt',header=T)
# calculate select and nonselect mut freq
Mut_select = Mut[which(as.matrix(Mut$X.donor_species)==as.matrix(Mut$gene)),
                 c(2,24,22,20,18,16,14)]
Mut_select[,c(2:7)]=0
Mut_nonselect=Mut_select
for(species in Mut_select$gene)
{
  temp_mut=Mut[which(as.matrix(Mut$X.donor_species)==species),
               c(2,24,22,20,18,16,14,26)]
  for(i in 1:nrow(temp_mut))
  {
    if(temp_mut$gene[i]!=species)
    {
      if (temp_mut$High_selected[i] == 'True'){
        Mut_select[which(as.matrix(Mut_select$gene)==species),c(2:7)]=
          Mut_select[which(as.matrix(Mut_select$gene)==species),c(2:7)]+
          temp_mut[i,c(2:7)]
      }else{
        Mut_nonselect[which(as.matrix(Mut_nonselect$gene)==species),c(2:7)]=
          Mut_nonselect[which(as.matrix(Mut_nonselect$gene)==species),c(2:7)]+
          temp_mut[i,c(2:7)]
      }
    }
  }
}

Mut_nonselect[nrow(Mut_nonselect),c(2:7)]=Mut[which(as.matrix(Mut$gene)=='allspecies'),
                                        c(24,22,20,18,16,14)]

Mut_select[nrow(Mut_select),c(2:7)]=Mut[which(as.matrix(Mut$gene)=='allspecies_highselect'),
                                        c(24,22,20,18,16,14)]
Mut_nonselect[nrow(Mut_nonselect),c(2:7)] = Mut_nonselect[nrow(Mut_nonselect),c(2:7)] -
  Mut_select[nrow(Mut_select),c(2:7)]
write.table(Mut_select,'mutfreq/species_mut_freq.select.txt',
            quote=F,sep='\t',row.names=F)
write.table(Mut_nonselect,'mutfreq/species_mut_freq.nonselect.txt',
            quote=F,sep='\t',row.names=F)

# sum up for species
select = 'non' #'all' or 'select' or 'non'
Mut2 = read.delim('final_SNP/HSgene/all.species.High_select2.all.short.txt',header=T)
if (select == 'all'){
 Mut=Mut[,c(2,24,22,20,18,16,14)] #all
}else if (select == 'select'){
 Mut = Mut_select #select
}else{
  Mut = Mut_nonselect 
}
Mut2=merge(Mut,Mut2,by.y='gene',by.x='gene')
Mut2$No.SNP=rowSums(as.matrix(Mut2[,c(2:7)]))

Mut2=Mut2[which(Mut2$No.SNP > 0),]

Mut3 = data.frame(
  population = Mut2$gene,
  freq = Mut2$G.A_freq.x,
  Tag='G.A_freq',
  species = Mut2$species_short,
  No.SNP = Mut2$No.SNP,
  genus = Mut2$genus
)
Mut3 =rbind(Mut3,
            data.frame(
              population = Mut2$gene,
              freq = Mut2$A.G_freq.x,
              Tag='A.G_freq',
              species = Mut2$species_short,
              No.SNP = Mut2$No.SNP,
              genus = Mut2$genus
            ))
Mut3 =rbind(Mut3,
            data.frame(
              population = Mut2$gene,
              freq = Mut2$G.T_freq.x,
              Tag='G.T_freq',
              species = Mut2$species_short,
              No.SNP = Mut2$No.SNP,
              genus = Mut2$genus
            ))
Mut3 =rbind(Mut3,
            data.frame(
              population = Mut2$gene,
              freq = Mut2$G.C_freq.x,
              Tag='G.C_freq',
              species = Mut2$species_short,
              No.SNP = Mut2$No.SNP,
              genus = Mut2$genus
            ))
Mut3 =rbind(Mut3,
            data.frame(
              population = Mut2$gene,
              freq = Mut2$A.C_freq.x,
              Tag='A.C_freq',
              species = Mut2$species_short,
              No.SNP = Mut2$No.SNP,
              genus = Mut2$genus
            ))
Mut3 =rbind(Mut3,
            data.frame(
              population = Mut2$gene,
              freq = Mut2$A.T_freq.x,
              Tag='A.T_freq',
              species = Mut2$species_short,
              No.SNP = Mut2$No.SNP,
              genus = Mut2$genus
            ))
if (select == 'all'){
  write.table(Mut3,'mutfreq/species_mut_freq.txt',
                          quote=F,sep='\t',row.names=F) #all
}else if (select == 'select'){
  write.table(Mut3,'mutfreq/species_mut_freq.select.txt',
              quote=F,sep='\t',row.names=F) #select
}else{
  write.table(Mut3,'mutfreq/species_mut_freq.nonselect.txt',
              quote=F,sep='\t',row.names=F)
}

library(ggplot2)
library(RColorBrewer)

# statistic test
Mut2$G.A_freq=as.numeric(as.matrix(Mut2$G.A_freq.x))/
  as.numeric(as.matrix(Mut2$No.SNP))
Mut2$A.G_freq=as.numeric(as.matrix(Mut2$A.G_freq.x))/
  as.numeric(as.matrix(Mut2$No.SNP))
Mut2$G.T_freq=as.numeric(as.matrix(Mut2$G.T_freq.x))/
  as.numeric(as.matrix(Mut2$No.SNP))
Mut2$G.C_freq=as.numeric(as.matrix(Mut2$G.C_freq.x))/
  as.numeric(as.matrix(Mut2$No.SNP))
Mut2$A.C_freq=as.numeric(as.matrix(Mut2$A.C_freq.x))/
  as.numeric(as.matrix(Mut2$No.SNP))
Mut2$A.T_freq=as.numeric(as.matrix(Mut2$A.T_freq.x))/
  as.numeric(as.matrix(Mut2$No.SNP))
colnames(Mut2)
if (select == 'all'){
  write.table(Mut2,'mutfreq/species_mut_freq.all.txt',
                          quote=F,sep='\t',row.names=F) #all
}else if (select == 'select'){
  write.table(Mut2,'mutfreq/species_mut_freq.all.select.txt',
              quote=F,sep='\t',row.names=F) #select
}else{
  write.table(Mut2,'mutfreq/species_mut_freq.all.nonselect.txt',
              quote=F,sep='\t',row.names=F)
}

Mut2=read.delim('final_SNP/mutfreq/species_mut_freq.all.select.txt',header=T)
Mut2$HS='HS'
Mut2_2=read.delim('final_SNP/mutfreq/species_mut_freq.all.nonselect.txt',header=T)
Mut2_2$HS='non_HS'
Mut2=rbind(Mut2,Mut2_2)
Mut2_2=Mut2
library(plyr)
genus_set = count(Mut2_2,c('genus'))
allgenus_set = genus_set$genus[which(genus_set$freq>5)]
testall=matrix(0,ncol=5,nrow=0)

for(genus in allgenus_set)
{
  Mut2=Mut2_2[which(Mut2_2$genus==genus),]
  if(nrow(Mut2)>1)
{Tag = 'G.A_freq'
library(effsize)
#Nonparametric Tests of Group Difference
#wilcox test
# independent 2-group Mann-Whitney U Test 
temp=wilcox.test(G.A_freq ~ HS, Mut2) 
testall=rbind(testall,c(temp$statistic,temp$p.value,'wilcos',Tag,genus))

#kruskal test
# Kruskal Wallis Test One Way Anova by Ranks 
temp=kruskal.test(Mut2$G.A_freq, as.factor(Mut2$HS) ) 
testall=rbind(testall,c(temp$statistic,temp$p.value,'kruskal',Tag,genus))

#calculate effect size
temp=cliff.delta(Mut2$G.A_freq, as.matrix(Mut2$HS) )
testall=rbind(testall,c(temp$estimate,levels(temp$magnitude)[temp$magnitude],
                        'effect size',Tag,genus))

Tag = 'A.G_freq'
#wilcox test
# independent 2-group Mann-Whitney U Test 
temp=wilcox.test(A.G_freq ~ HS, Mut2) 
testall=rbind(testall,c(temp$statistic,temp$p.value,'wilcos',Tag,genus))

#kruskal test
# Kruskal Wallis Test One Way Anova by Ranks 
temp=kruskal.test(Mut2$A.G_freq, as.factor(Mut2$HS) ) 
testall=rbind(testall,c(temp$statistic,temp$p.value,'kruskal',Tag,genus))

#calculate effect size
temp=cliff.delta(Mut2$A.G_freq,as.matrix(Mut2$HS))
testall=rbind(testall,c(temp$estimate,levels(temp$magnitude)[temp$magnitude],
                        'effect size',Tag,genus))

Tag = 'G.C_freq'
#wilcox test
# independent 2-group Mann-Whitney U Test 
temp=wilcox.test(G.C_freq ~ HS, Mut2) 
testall=rbind(testall,c(temp$statistic,temp$p.value,'wilcos',Tag,genus))

#kruskal test
# Kruskal Wallis Test One Way Anova by Ranks 
temp=kruskal.test(Mut2$G.C_freq, as.factor(Mut2$HS) ) 
testall=rbind(testall,c(temp$statistic,temp$p.value,'kruskal',Tag,genus))

#calculate effect size
temp=cliff.delta(Mut2$G.C_freq,as.matrix(Mut2$HS))
testall=rbind(testall,c(temp$estimate,levels(temp$magnitude)[temp$magnitude],
                        'effect size',Tag,genus))


Tag = 'G.T_freq'
#wilcox test
# independent 2-group Mann-Whitney U Test 
temp=wilcox.test(G.T_freq ~ HS, Mut2) 
testall=rbind(testall,c(temp$statistic,temp$p.value,'wilcos',Tag,genus))

#kruskal test
# Kruskal Wallis Test One Way Anova by Ranks 
temp=kruskal.test(Mut2$G.T_freq, as.factor(Mut2$HS) ) 
testall=rbind(testall,c(temp$statistic,temp$p.value,'kruskal',Tag,genus))

#calculate effect size
temp=cliff.delta(Mut2$G.T_freq,as.matrix(Mut2$HS))
testall=rbind(testall,c(temp$estimate,levels(temp$magnitude)[temp$magnitude],
                        'effect size',Tag,genus))

Tag = 'A.C_freq'
#wilcox test
# independent 2-group Mann-Whitney U Test 
temp=wilcox.test(A.C_freq ~ HS, Mut2) 
testall=rbind(testall,c(temp$statistic,temp$p.value,'wilcos',Tag,genus))

#kruskal test
# Kruskal Wallis Test One Way Anova by Ranks 
temp=kruskal.test(Mut2$A.C_freq, as.factor(Mut2$HS) ) 
testall=rbind(testall,c(temp$statistic,temp$p.value,'kruskal',Tag,genus))

#calculate effect size
temp=cliff.delta(Mut2$A.C_freq,as.matrix(Mut2$HS))
testall=rbind(testall,c(temp$estimate,levels(temp$magnitude)[temp$magnitude],
                        'effect size',Tag,genus))

Tag = 'A.T_freq'
#wilcox test
# independent 2-group Mann-Whitney U Test 
temp=wilcox.test(A.T_freq ~ HS, Mut2) 
testall=rbind(testall,c(temp$statistic,temp$p.value,'wilcos',Tag,genus))

#kruskal test
# Kruskal Wallis Test One Way Anova by Ranks 
temp=kruskal.test(Mut2$A.T_freq, as.factor(Mut2$HS) ) 
testall=rbind(testall,c(temp$statistic,temp$p.value,'kruskal',Tag,genus))

#calculate effect size
temp=cliff.delta(Mut2$A.T_freq,as.matrix(Mut2$HS))
testall=rbind(testall,c(temp$estimate,levels(temp$magnitude)[temp$magnitude],
                        'effect size',Tag,genus))
  }
  }

colnames(testall)=c('estimate','p-value',
                    'test_method','test_project','genus')
write.table(testall,'mutfreq/species_mut_freq.select.nonselect.t-test.txt',
            quote=F,sep='\t',row.names=F)

# re-plot
Mut3=read.delim('final_SNP/mutfreq/species_mut_freq.select.txt',header=T)
Mut3$HS='HS'
Mut3_2=read.delim('final_SNP/mutfreq/species_mut_freq.nonselect.txt',header=T)
Mut3_2$HS='non_HS'
Mut3=rbind(Mut3_2,Mut3)
Mut3=Mut3[which(Mut3$genus %in% allgenus_set[c(1,2,3,4)]),]

ggplot(Mut3,aes(y=freq/No.SNP,
                x=paste(genus,Tag,HS),
                group = genus))+
  geom_boxplot(aes(
    color=HS,
    fill=HS
  ),
  outlier.shape = NA,
  alpha = 0,
  size=0
  )+ stat_summary(fun = mean, geom = "bar", 
                  aes(color=HS,fill = HS),alpha = 0.5,
                  size=0)+
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               aes(color=HS,
                   width = 0.8,
                   linetype = "solid"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1,5)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1,5)])

library(ggtern)

Mut2=Mut2[which(Mut2$genus %in% genus_set),]
ggtern(data = Mut2, aes(x=as.numeric(as.matrix(G.A_freq)),
                        y=(as.numeric(as.matrix(A.T_freq))),
                        z=1-as.numeric(as.matrix(G.A_freq))-(as.numeric(as.matrix(A.T_freq)))
)) +
  
  #the layers
  geom_point(aes(fill = genus,
                 color = genus),
             size = 3, 
             alpha = 0.5,
             shape = 21)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "Set1")[c(5,4,3,2)])+
  scale_color_manual(values =brewer.pal(5, "Set1")[c(5,4,3,2)])

#----species enrichment gene----
# HS genes
anno = read.delim('final_SNP/annotation/all.selected.gene.faa.High_select2.faa.cluster.sum.new.txt',header=F)
anno=anno[,-9]
colnames(anno)=as.matrix(anno[1,])
colnames(anno)[9:16]=c('fungroup',
                       'description',
                       'anno1',
                       'anno2',
                       'anno3',
                       'note',
                       'gene',
                       'description2')
anno=anno[-1,]
anno=anno[,c(1,7:16)]
write.table(anno[,c(1,2,4,6,8)],'HSgene/all.cluster.anno.txt',
            quote=F,sep='\t',row.names=F)

kegg=read.delim('final_SNP/annotation/all.selected.gene.faa.High_select2.faa.cluster.aa.kegg.sum.txt',header=T)

changename = read.delim('final_SNP/pangenome/all.denovo.gene.faa.changename.txt',header=F)
changename=changename[,c(1:4)]
colnames(changename)=c('donor_species','short_species','gene1','gene2')
NS_ratio = read.delim('final_SNP/HSgene/all.species.High_select2.txt',header=T)
#NS_ratio2 = NS_ratio[,c(1,2,3,5,6,7,8,9,26)]
NS_ratio2 = NS_ratio[,c(1,2,3,4,5,6,7,8,9,26)]

NS_ratio2$allname = paste(NS_ratio2$X.donor_species,NS_ratio2$gene)
changename$allname =paste(changename$donor_species,changename$gene1)
NS_ratio2=merge(NS_ratio2,changename,by='allname')
library(stringr)

NS_ratio2$allname = paste(str_split_fixed(NS_ratio2$short_species, "_CL", 2)[,1],NS_ratio2$gene2)
NS_ratio2 = NS_ratio2[which(NS_ratio2$High_selected=='False' &
                            NS_ratio2$total_SNP_genomeset > 1),-6]
#NS_ratio2 = NS_ratio2[which(NS_ratio2$High_selected=='True'),]

kegg$allname = paste(kegg$donor_species,kegg$gene_name)
kegg=merge(kegg,anno[,c(1,2)],by='cluster',all.x=T)
kegg=kegg[!duplicated(paste(kegg$allname,kegg$BRITE_KO3)),]
kegg2=merge(kegg[,c(10,1,6:11)],
           NS_ratio2[,c(1,5,6,8:14)],by='allname',all.y=T)
kegg2$species = str_split_fixed(kegg2$short_species, "_CL", 2)[,1]
kegg2$species = paste(str_split_fixed(kegg2$species, "_", 4)[,2],
                      str_split_fixed(kegg2$species, "_", 4)[,3],
                      sep='_')
unique(kegg2$species)
kegg2$species[which(kegg2$species=='BA_BA')]='Bifido_adoles'
kegg2$species[which(kegg2$species=='BL_BL')]='Bifido_longum'
kegg2$species[which(kegg2$species=='BA_cluste')]='Bifido_adoles'
kegg2$species[which(kegg2$species=='BL_cluste')]='Bifido_longum'
kegg2$species[which(kegg2$species=='PB_cluste')]='PB_PB'
unique(kegg2$species)
kegg2$BRITE_KO3[which(is.na(kegg2$BRITE_KO3))]='None'
kegg2$BRITE_KO3[which(kegg2$BRITE_KO3=='')]='None'
kegg2$BRITE_KO2[which(is.na(kegg2$BRITE_KO2))]='None'
kegg2$BRITE_KO2[which(kegg2$BRITE_KO2=='')]='None'
kegg2$BRITE_KO1[which(is.na(kegg2$BRITE_KO1))]='None'
kegg2$BRITE_KO1[which(kegg2$BRITE_KO1=='')]='None'
kegg2=kegg2[order(kegg2$BRITE_KO3),]
kegg2=kegg2[!duplicated(paste(kegg2$allname,kegg2$BRITE_KO3)),]
kegg2=kegg2[which(kegg2$N <= kegg2$S &
                    kegg2$No.SNP/kegg2$gene_length >= 1/2000),]
#write.table(kegg2[,c(17,2,8,18,15,3:6)],'HSgene/all.gene.kegg.HS.anno.final.txt',
#            quote=F,sep='\t',row.names=F)
write.table(kegg2[,c(17,2,8,18,15,3:6)],'HSgene/all.gene.kegg.HS.anno.final.FP.txt',
            quote=F,sep='\t',row.names=F)
# enriched gene
enrichment=read.delim('final_SNP/enrichment/all.gene.kegg.enrichment.norm.species.txt',header=T)
enrichment=enrichment[which(enrichment$diff>=1),]
kegg=read.delim('final_SNP/enrichment/all.gene.kegg.enrichment.txt',header=T)
kegg$speciesfun = paste(kegg$species,kegg$BRITE_KO2)
kegg=merge(kegg,enrichment[,c(1,2,9,10)],by='speciesfun')
anno = read.delim('final_SNP/annotation/all.selected.gene.faa.High_select2.faa.cluster.sum.new.txt',header=F)
anno=anno[,-9]
colnames(anno)=as.matrix(anno[1,])
colnames(anno)[9:16]=c('fungroup',
                        'description',
                        'anno1',
                        'anno2',
                        'anno3',
                        'note',
                        'gene',
                        'description2')
anno=anno[-1,]
anno=anno[,c(1,7:16)]
kegg=merge(kegg,anno,by='cluster',all.x=T)
kegg=kegg[which(kegg$Tag == 'High_select'),]
write.table(kegg,'enrichment/all.gene.kegg.enrichment.species.norm.enriched.anno.final.txt',
            quote=F,sep='\t',row.names=F)
kegg=read.delim('final_SNP/enrichment/all.gene.kegg.enrichment.species.norm.enriched.anno.final.txt',header=T)
kegg_count = count(kegg,'species')

#----heatmap selected genes ----
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
i = 1
species_set = c('Bacter.Bacter_sterco',
'Paraba.Paraba_distas',
                'Paraba.Paraba_merdae',
                'Clostr.Clostr_innocu',
                'Clostr.Clostr_beijer',
                'Ruthen.Ruthen_lactat',
                'Lactob.Lactob_rumini',
                'Bacter.Bacter_fragil',
                'Blauti.Blauti_wexler',
                'Collin.Collin_aerofa',
                'Parasu.Parasu_excrem',
                'Pseudo.Pseudo_sp.',
                'Turici.Turici_sangui',
                'Bacter.Bacter_salyer',
                'Bacter.Bacter_ovatus',
                'Bacter.Bacter_vulgat',
                'Escher.Escher_coli',
                'Bacter.Bacter_thetai',
                'PB.PB_PB',
                'Bifido.Bifido_pseudo',
                'Bifido.Bifido_longum',
                'Bifido.Bifido_adoles')

species = species_set[i]
species
{
  Metadata1 = read.delim(
  paste(paste('speciesgenes/all.vcf.frq.snp.',species,sep=''),'.fun.txt',sep=''),header=F)
Mat_snp=read.delim(
  paste(paste('speciesgenes/all.vcf.frq.snp.',species,sep=''),'.matrix.txt',sep=''),header=T)
row.names(Mat_snp)=Mat_snp[,1]
Mat_snp=Mat_snp[,-c(1,ncol(Mat_snp))]
Mat_snp=Mat_snp[order(row.names(Mat_snp)),]
Mat_snp2=data.frame(t(Mat_snp))
if (nrow(Metadata1)>90)
{
  Mat_snp2=Mat_snp2[which(Metadata1$V2!='None'),]
  Metadata1=Metadata1[which(Metadata1$V2!='None'),]
}
Color_set = c('#e5e5e5','#ffffff','#ca0020',
              '#f4a582','#92c5de','#0571b0')

# heatmap pre
Mat_snp = Mat_snp2
mat=as.matrix(Mat_snp)
mat[which(mat=='')]=-1
mat[which(mat=='REF')]=0
mat[which(mat=='A')]=11
mat[which(mat=='T')]=12
mat[which(mat=='G')]=13
mat[which(mat=='C')]=14
mat=matrix(as.numeric(mat),nrow = nrow(mat))
colnames(mat)=colnames(Mat_snp)
row.names(mat)=row.names(Mat_snp)

Metadata2=data.frame(group = c(as.matrix(Metadata1$V2)))
row.names(Metadata2)=c(row.names(mat))
Metadata2$newname = paste(Metadata2$group,row.names(mat),sep='.')
row.names(mat)=c(Metadata2$newname)
#heatmap2
Metadata2=data.frame(SNP_type = c(as.matrix(Metadata1$V3),'N','S'))
row.names(Metadata2)=c(row.names(mat),'Ref1','Ref2')
mat=mat[order(row.names(mat)),]
Metadata2=data.frame(as.matrix(Metadata2)[order(Metadata2$SNP_type),])
colnames(Metadata2)='SNP_type'
mat_colors <- list(SNP_type = brewer.pal(5, "PuOr")[c(5,1)])
names(mat_colors$SNP_type) <- unique(Metadata2$SNP_type)
mat=mat[order(row.names(mat)),]
mat_breaks <- c(-2,-0.5,10.5,11.5,12.5,13.5,14.5)
pheatmap(mat= mat,
         color = Color_set,
         show_rownames= TRUE,
         show_colnames = TRUE,
         breaks= mat_breaks,
         annotation_row    = Metadata2,
         annotation_colors = mat_colors,
         clustering_distance_rows = "correlation",
        # cluster_cols = TRUE,
         cluster_rows = FALSE,
         border_color = Color_set[2],
         legend_breaks =   c(-1,0,11,12,13,14),
         legend_labels = c('Not present','Major allele','A','T','G','C'),
         main = species)
}
i=i+1
#----heatmap not selected genes ----
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
i = 1
species_set = c('Bacter.Bacter_caccae',
                'Bacter.Bacter_fragil',
                'Bacter.Bacter_ovatus',
                'Bacter.Bacter_thetai',
                'Bacter.Bacter_vulgat',
                'Bacter.Bacter_xylani',
                'Bifido.Bifido_adoles',
                'Bifido.Bifido_longum',
                'Bifido.Bifido_pseudo',
                'Blauti.Blauti_wexler',
                'Clostr.Clostr_innocu',
                'Collin.Collin_aerofa',
                'Entero.Entero_hirae',
                'Entero.Entero_mundti',
                'Escher.Escher_coli',
                'Eubact.Eubact_rectal',
                'Faecal.Faecal_prausn',
                'Parasu.Parasu_excrem',
                'Turici.Turici_sangui')

species = species_set[i]
species = 'Bifido.Bifido_adoles'
{
  Metadata1 = read.delim(
    paste(paste('nospeciesgenes/all.vcf.frq.snp.',species,sep=''),'.fun.txt',sep=''),header=F)
  Mat_snp=read.delim(
    paste(paste('nospeciesgenes/all.vcf.frq.snp.',species,sep=''),'.matrix.txt',sep=''),header=T)
  row.names(Mat_snp)=Mat_snp[,1]
  Mat_snp=Mat_snp[,-c(1,ncol(Mat_snp))]
  Mat_snp=Mat_snp[order(row.names(Mat_snp)),]
  Mat_snp2=data.frame(t(Mat_snp))
  Color_set = c('#e5e5e5','#ffffff','#ca0020',
                '#f4a582','#92c5de','#0571b0')
  
  # heatmap pre
  Mat_snp = Mat_snp2
  mat=as.matrix(Mat_snp)
  mat[which(mat=='')]=-1
  mat[which(mat=='REF')]=0
  mat[which(mat=='A')]=11
  mat[which(mat=='T')]=12
  mat[which(mat=='G')]=13
  mat[which(mat=='C')]=14
  mat=matrix(as.numeric(mat),nrow = nrow(mat))
  colnames(mat)=colnames(Mat_snp)
  row.names(mat)=row.names(Mat_snp)
  
  Metadata2=data.frame(group = c(as.matrix(Metadata1$V2)))
  row.names(Metadata2)=c(row.names(mat))
  Metadata2$newname = paste(Metadata2$group,row.names(mat),sep='.')
  row.names(mat)=c(Metadata2$newname)
  #heatmap2
  Metadata2=data.frame(SNP_type = c(as.matrix(Metadata1$V3),'N','S'))
  row.names(Metadata2)=c(row.names(mat),'Ref1','Ref2')
  mat=mat[order(row.names(mat)),]
  Metadata2=data.frame(as.matrix(Metadata2)[order(Metadata2$SNP_type),])
  colnames(Metadata2)='SNP_type'
  mat_colors <- list(SNP_type = brewer.pal(5, "PuOr")[c(5,1)])
  names(mat_colors$SNP_type) <- unique(Metadata2$SNP_type)
  mat=mat[order(row.names(mat)),]
  mat_breaks <- c(-2,-0.5,10.5,11.5,12.5,13.5,14.5)
  pheatmap(mat= mat,
           color = Color_set,
           show_rownames= TRUE,
           show_colnames = TRUE,
           breaks= mat_breaks,
           annotation_row    = Metadata2,
           annotation_colors = mat_colors,
           clustering_distance_rows = "correlation",
           # cluster_cols = TRUE,
           cluster_rows = FALSE,
           border_color = Color_set[2],
           legend_breaks =   c(-1,0,11,12,13,14),
           legend_labels = c('Not present','Major allele','A','T','G','C'),
           main = species)
}
i=i+1
#----metagenome depth versus SNP----
# load cov
MG=read.delim('final_SNP/MG/allcov.sum.txt',header=T)
MG=MG[which(MG$donor_species!='donor_species'),-7]
MG2=MG[which(MG$avg_depth>0),]
total_sample = 205
library(plyr)
MG_pre=count(MG2,c('donor_species'))
MG_pre=MG_pre[which(MG_pre$freq>1),]
MG2=merge(MG2,MG_pre,by='donor_species')
MG2$prevalence=MG2$freq/total_sample
MG2=MG2[!duplicated(MG2$donor_species),]
MG2$avg_depth = 0
MG2$CV_depth = 0
MG2$avg_depth2 = 0
MG2$CV_depth2 = 0
for(i in 1:nrow(MG2))
{
  MG_sub = MG[which(MG$donor_species == as.matrix(MG2$donor_species)[i]),]
  MG2$avg_depth[i] = mean(MG_sub$avg_depth)
  MG2$CV_depth[i] = sd(MG_sub$avg_depth)/mean(MG_sub$avg_depth)    
  MG_sub = MG_sub[which(MG_sub$avg_depth>0),]
  MG2$avg_depth2[i] = mean(MG_sub$avg_depth)
  # corrected by qualified sample number
  MG2$CV_depth2[i] = sd(MG_sub$avg_depth)/mean(MG_sub$avg_depth)*(1+1/4/nrow(MG_sub))  
}
# prevalence depth SNP
library(ggplot2)
library(RColorBrewer)
allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.all.short.txt',header=T)
#allspecies=allspecies[which(allspecies$donor=='am'),]
allspecies=merge(allspecies,MG2[,c(1,3,4,8:11)],
                 by.x='X.donor_species',by.y='donor_species')

ggplot(allspecies,aes(x=log10(avg_depth),
                      y=prevalence
                      ))+
  geom_point(aes(
    color=genus,
    fill=genus
  ),
  alpha = 0.8,
  size=2
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(7, "Set1"))+
  scale_color_manual(values =brewer.pal(7, "Set1"))

ggplot(allspecies,aes(x=log10(avg_depth2),
                      y=prevalence
))+
  geom_point(aes(
    color=genus,
    fill=genus
  ),
  alpha = 0.8,
  size=2
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(7, "Set1"))+
  scale_color_manual(values =brewer.pal(7, "Set1"))

ggplot(allspecies,aes(y=log10(CV_depth2),
                      x=log10(avg_depth2)
))+
  geom_point(aes(
    color=genus,
    fill=genus
  ),
  alpha = 0.8,
  size=2
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(7, "Set1"))+
  scale_color_manual(values =brewer.pal(7, "Set1"))

allspecies$HS_NO[is.na(allspecies$HS_NO)]=0
allspecies$HS_SNP_NO[is.na(allspecies$HS_SNP_NO)]=0
ggplot(allspecies,aes(x=(CV_depth2),
                      y=(No.SNP)
))+
  geom_point(aes(
    color=genus,
    fill=genus
  ),
  alpha = 1,
  size=2
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(7, "Set1"))+
  scale_color_manual(values =brewer.pal(7, "Set1"))

#---- metagenome----
# run simple species
species_set = c('am_Bacteroides_fragilis_cluster1',
                'am_Bacteroides_fragilis_cluster2',
                'am_Bacteroides_salyersiae_cluster1', 
                'am_Bacteroides_vulgatus_cluster1', 
                'am_Bacteroides_vulgatus_cluster2', 
                'am_Bifidobacterium_adolescentis_cluster1', 
                'am_Bifidobacterium_adolescentis_cluster2', 
                'am_Bifidobacterium_longum_cluster1', 
                'am_Bifidobacterium_longum_cluster2', 
                'am_Bifidobacterium_longum_cluster3', 
                'am_Bifidobacterium_longum_cluster4', 
                'am_Bifidobacterium_longum_cluster5', 
                'am_Bifidobacterium_longum_cluster6', 
                'am_Bifidobacterium_pseudocatenulatum_cluster1', 
                'am_Collinsella_aerofaciens_cluster1', 
                'am_Escherichia_coli_cluster1', 
                'am_Escherichia_coli_cluster2', 
                'am_Escherichia_coli_cluster3', 
                'am_Escherichia_coli_cluster4', 
                'am_Escherichia_coli_cluster5', 
                'am_Turicibacter_sanguinis_cluster1', 
                'am_Turicibacter_sanguinis_cluster4')
library(ggplot2)
library(RColorBrewer)

i = 3
species = species_set[i]
metaspecies = read.delim(paste('MG/',
                               paste(species,'.all.flt.snp.vcf.freq',sep=''),sep=''
                               ),
                         header = F)
colnames(metaspecies)=as.matrix(metaspecies[1,])
colnames(metaspecies)[1]='seq_type'
metaspecies=metaspecies[-1,-ncol(metaspecies)]
metaspecies=data.frame(seq_type = metaspecies$seq_type,
                       timepoint=as.numeric(as.matrix(metaspecies$timepoint)),
                       freq = as.numeric(as.matrix(metaspecies$freq)),
                       closest_relative=metaspecies$closest_relative,
                       SNP = as.numeric(as.matrix(metaspecies$SNP)))
metaspecies_snptype = metaspecies[!duplicated(metaspecies$seq_type),]

# correct MLF
timepoint = unique(metaspecies$timepoint)
for(time_one in timepoint)
{
  metaspecies_multitime = metaspecies[which(metaspecies$timepoint==time_one),]
  if(nrow(metaspecies_multitime)>1)
  {metaspecies$freq[which(metaspecies$timepoint==time_one &
                      metaspecies$seq_type == metaspecies_multitime$seq_type[1])]=
    1-as.numeric(as.matrix(metaspecies_multitime$freq[2]))
  }else
  {
    metaspecies$freq[which(metaspecies$timepoint==time_one &
                             metaspecies$seq_type == metaspecies_multitime$seq_type[1])]=1
  }
  for(seq_type in metaspecies_snptype$seq_type){
    if(!(seq_type %in% metaspecies_multitime$seq_type))
    {
      metaspecies=rbind(metaspecies,
                        c(seq_type,time_one,0,
                          as.matrix(metaspecies_snptype$closest_relative[which(metaspecies_snptype$seq_type == seq_type)]),
                          metaspecies_snptype$SNP[which(metaspecies_snptype$seq_type == seq_type)])
                        )
    }
  }
}
metaspecies=data.frame(seq_type = metaspecies$seq_type,
                       timepoint=as.numeric(as.matrix(metaspecies$timepoint)),
                       freq = as.numeric(as.matrix(metaspecies$freq)),
                       closest_relative=metaspecies$closest_relative,
                       SNP = as.numeric(as.matrix(metaspecies$SNP)))

# re plot
# set color accordding to origin
ggplot(metaspecies,aes(x=timepoint,
                       y=freq
))+
  geom_area(
    aes(fill = seq_type,
        color = seq_type),
    alpha = 1,
    size=0
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  )+
  scale_fill_manual(values =brewer.pal(9, "PRGn")[c(1,2,3,4,6,7,8,9)])+
  scale_color_manual(values =brewer.pal(9, "PRGn")[c(1,2,3,4,6,7,8,9)])+
  scale_x_continuous(breaks = seq(0,max(metaspecies$timepoint),
                                  by = 20))
#----metagenome complicated species----
species_set = c('am_Bacteroides_fragilis','am_Bacteroides_ovatus',
                'am_Bacteroides_stercoris','am_Bacteroides_vulgatus',
                'am_Bacteroides_vulgatus_cluster2','am_Bifidobacterium_longum')
library(ggplot2)
library(RColorBrewer)

i = 3 # 2-6
species = species_set[i]
metaspecies = read.delim(paste('meta/',
                               paste(species,'.all.flt.snp.vcf.freq',sep=''),sep=''
),
header = F)
colnames(metaspecies)=as.matrix(metaspecies[1,])
colnames(metaspecies)[1]='seq_type'
metaspecies=metaspecies[-1,-ncol(metaspecies)]
metaspecies=metaspecies[order(metaspecies$seq_type),]
metaspecies=data.frame(seq_type = metaspecies$seq_type,
                       timepoint=as.numeric(as.matrix(metaspecies$timepoint)),
                       freq = as.numeric(as.matrix(metaspecies$freq)),
                       closest_relative=metaspecies$closest_relative,
                       SNP = as.numeric(as.matrix(metaspecies$SNP)))
metaspecies2=metaspecies
max(metaspecies$SNP)
metaspecies_snptype = metaspecies[!duplicated(metaspecies$seq_type),]
if (nrow(metaspecies_snptype)>10)
{
  metaspecies$SNP[which(metaspecies$SNP>6&metaspecies$SNP<=8)]=8
  metaspecies$SNP[which(metaspecies$SNP>4&metaspecies$SNP<=6)]=6
  metaspecies$SNP[which(metaspecies$SNP>2&metaspecies$SNP<=4)]=4
  metaspecies$SNP[which(metaspecies$SNP>0&metaspecies$SNP<=2)]=2
  metaspecies=data.frame(seq_type = paste(metaspecies$closest_relative,
                                          metaspecies$SNP,sep='_'),
                         timepoint=as.numeric(as.matrix(metaspecies$timepoint)),
                         freq = as.numeric(as.matrix(metaspecies$freq)),
                         closest_relative=metaspecies$closest_relative,
                         SNP = as.numeric(as.matrix(metaspecies$SNP)))
  metaspecies_snptype = metaspecies[!duplicated(metaspecies$seq_type),]
}

metaspecies=metaspecies[rev(order(metaspecies$freq)),]
metaspecies=metaspecies[!duplicated(paste(metaspecies$seq_type,metaspecies$timepoint)),]
metaspecies$percentage = 0
# correct MLF
timepoint = unique(metaspecies$timepoint)
for(time_one in timepoint)
{
  metaspecies_multitime = metaspecies[which(metaspecies$timepoint==time_one),]
  for(seq_type in metaspecies_snptype$seq_type){
    if(!(seq_type %in% metaspecies_multitime$seq_type))
    {
      metaspecies=rbind(metaspecies,
                        c(seq_type,time_one,0,
                          as.matrix(metaspecies_snptype$closest_relative[which(metaspecies_snptype$seq_type == seq_type)]),
                          metaspecies_snptype$SNP[which(metaspecies_snptype$seq_type == seq_type)])
      )
    }
  }
  metaspecies_multitime = metaspecies[which(metaspecies$timepoint==time_one),]
  
  metaspecies$percentage[which(metaspecies$timepoint==time_one)
                         ]=as.numeric(as.matrix(metaspecies_multitime$freq))/sum(as.numeric(as.matrix(metaspecies_multitime$freq)))
}
metaspecies=metaspecies[order(metaspecies$seq_type),]
metaspecies=data.frame(seq_type = metaspecies$seq_type,
                       timepoint=as.numeric(as.matrix(metaspecies$timepoint)),
                       freq = as.numeric(as.matrix(metaspecies$percentage)),
                       closest_relative=metaspecies$closest_relative,
                       SNP = as.numeric(as.matrix(metaspecies$SNP)))


# trial plot
ggplot(metaspecies,aes(x=timepoint,
                       y=freq
))+
  geom_area(
    aes(fill = seq_type,
        color = seq_type),
    alpha = 0.6,
    size=0.5
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_brewer(palette='Set1')+
  scale_color_brewer(palette='Set1')+
  scale_x_continuous(breaks = seq(0,max(metaspecies$timepoint),
                                  by = 20))

# re plot
metaspecies_snptype = metaspecies[!duplicated(metaspecies$seq_type),]

metaspecies_snptype$color = c(
  brewer.pal(9, "PRGn")[1],
  brewer.pal(9, "PRGn")[2],
  brewer.pal(9, "PRGn")[3],
  brewer.pal(9, "PRGn")[4],
  brewer.pal(9, "PRGn")[5]
)
# set color accordding to origin
ggplot(metaspecies,aes(x=timepoint,
                       y=freq
))+
  geom_area(
    aes(fill = seq_type,
        color = seq_type),
    alpha = 1,
    size=0
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  )+
  scale_fill_manual(values =c(metaspecies_snptype$color))+
  scale_color_manual(values =c(metaspecies_snptype$color))+
  scale_x_continuous(breaks = seq(0,max(metaspecies$timepoint),
                                  by = 20))

c(
  brewer.pal(9, "PRGn")[9],
  brewer.pal(9, "PRGn")[8],
  brewer.pal(9, "PRGn")[7],
  brewer.pal(9, "PRGn")[6],
  brewer.pal(9, "PRGn")[1],
  brewer.pal(9, "PRGn")[2],
  brewer.pal(9, "PRGn")[3],
  brewer.pal(9, "PRGn")[4]
)
#---- BN10 genome time tag finished----
genomename = read.delim('final_SNP/metadata/BN10_WGS_newname.txt',header=F)
colnames(genomename)=c('oldname','newnameold','newname')
BN10time = read.delim('final_SNP/metadata/table_BN10_WGS_metadata_07052018.txt',header=T)
BN10time=BN10time[,c(1:5)]
genomename=merge(genomename[,c(1:3)],BN10time,
                 by.x='oldname',by.y='Genome_name',all.x=T)

write.table(genomename,'metadata/BN10_WGS_newname_meta_time.txt',
            quote=F,sep='\t',row.names=F)
genomename = read.delim('final_SNP/BN10_WGS_newname_meta_time.txt',header=T)
genomename_multitime =genomename[!duplicated(paste(paste(genomename$Species,
                                                         genomename$Individual),
                                                   genomename$Time_point)),]
library(plyr)
genomename_multitime = count(genomename_multitime,c('Species','Individual'))
genomename_multitime = genomename_multitime[which(genomename_multitime$freq>1),]

# generate tree color file
genomename2=genomename[which(paste(genomename$Species,
                                   genomename$Individual) %in% 
                               paste(genomename_multitime$Species,
                                     genomename_multitime$Individual)),]

write.table(genomename2,'BN10_WGS_newname_meta_multitime.txt',
            quote=F,sep='\t',row.names=F)

alltime = unique(genomename2$Time_point)
alltime=data.frame(Time_point=alltime,
                   color = c('#54278f','#756bb1','#9e9ac8','#cbc9e2','#f2f0f7',
                             '#ffffb2','#fecc5c','#fd8d3c','#f03b20','#bd0026'))


genomename2=merge(genomename2,alltime,by='Time_point')

genomename2$shortname = ''
for(i in 1:nrow(genomename2))
{
  genomename2$shortname[i]=paste('S',str_sub(genomename2$newname[i],-8,-1),
                                 sep = '')
}
genomename2$donor_species = paste(genomename2$Individual,
                                  genomename2$Species,sep='_')
alldonor_species = unique(genomename2$donor_species)

write.table(genomename2,'metadata/BN10_WGS_newname_meta_multitime_color.txt',
            quote=F,sep='\t',row.names=F)
temptitle = 'DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL	Timepoint
COLOR	#ff0000
LEGEND_TITLE	Dataset_legend
LEGEND_SHAPES	1	1	1	1	1	1	1	1	1	1
LEGEND_COLORS	#54278f	#756bb1	#9e9ac8	#cbc9e2	#f2f0f7	#ffffb2	#fecc5c	#fd8d3c	#f03b20	#bd0026
LEGEND_LABELS	1	3	6	14	25	52	70	111	171	224
DATA'
for(donor_species in alldonor_species)
{
  tempgenomename2 = genomename2[which(genomename2$donor_species==donor_species),]
  tempfilename = paste(paste('metadata/',donor_species,sep = ''),
                       '.tree.COLORSTRIP.txt',sep='')
  write.table(temptitle,tempfilename,
              quote=F,sep='\t',row.names=F,append = FALSE,col.names = F)
  write.table(tempgenomename2[,c(9,8)],tempfilename,
              quote=F,sep='\t',row.names=F,append = TRUE,col.names = F)
}

#----summarize all SNPs finished----
SNP1 = read.delim('final_SNP/SNP_round1.allpair.all.sum',header=T)
SNP1=SNP1[,c(1:3)]
SNP1=SNP1[which(SNP1$Genome1!='reference' & 
                  SNP1$Genome2!='reference'),]
SNP1=SNP1[which(SNP1$Genome1!='Genome1' & 
                  SNP1$Genome1!='G1'),]
SNP1$SNP_total=as.numeric(
  as.matrix(SNP1$SNP_total)
) 
SNP2=SNP1
library(ggplot2)
library(RColorBrewer)
# same donor all SNPs
SNP1 =SNP2
ggplot(SNP1,aes(x=SNP_total))+
  geom_histogram(aes(#y=(..count..)/sum(..count..),
    color='1',
    fill='1'
  ),
  position='dodge',
  binwidth = 1000,
  alpha = 0.2,
  size=0.5
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1)])+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8)*10000)+
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14)*1E+4)

# less than 1000 SNPs
SNP1=SNP2[which(SNP2$SNP_total<=200),]
nrow(SNP1)/nrow(SNP2)
SNP1=SNP2[which(SNP2$SNP_total<=2000),]
ggplot(SNP1,aes(x=SNP_total))+
  geom_histogram(aes(#y=(..count..)/sum(..count..),
    color='1',
    fill='1'
  ),
  position='dodge',
  binwidth = 100,
  alpha = 0.2,
  size=0.5
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(1)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(1)])+
  scale_x_continuous(
    breaks = c(0,1,2,3,4,5,6,7,8,9,10)*200)

# same donor diff clusters
SNP1 = read.delim('final_SNP/SNP_round1.allpair.all.diffcluster.sumnew.txt',header=T)
SNP1$SNP_count=as.numeric(
  as.matrix(SNP1$SNP_count)
) 
SNP3=SNP1[which(SNP1$SNP_count>=10000),]
nrow(SNP3)/nrow(SNP1)
ggplot(SNP1,aes(x=SNP_count))+
  geom_histogram(aes(#y=(..count..)/sum(..count..),
    color='1',
    fill='1'
  ),
  position='dodge',
  binwidth = 1000,
  alpha = 0.2,
  size=0.5
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()#,
    # axis.ticks.y=element_blank(),
    # axis.text.y=element_blank()
  ) +
  scale_fill_manual(values =brewer.pal(5, "PRGn")[c(5)])+
  scale_color_manual(values =brewer.pal(5, "PRGn")[c(5)])+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8)*10000)

#---- predict HS ecology figure 1b prep----
library(binom)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(plyr)
library(rpart)
library(randomForest)
# C20,n*(10 genes, 90 untargeted genes, 20 SNPs)*(2 SNPs, n SNPs,)
#p = 1-phyper(SNP_cutoff - 1, m, n, k)
#binom.test(10 genes,100 genes,p, alternative="greater")

# sum table
allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.txt',header=T)
allspecies2=read.delim('final_SNP/HSgene/all.species.High_select2.all.gene.txt',header=T)
clonal = read.delim('clonal_pop_gene_num.sum.txt',header=T)
cluster_ratio = read.delim('cluster_ratio.txt',header=T)
clonal = merge(clonal,cluster_ratio,by = 'species',all.x=T)
clonal=clonal[!duplicated(clonal$cluster),]
allspecies$donor_species = str_split_fixed(allspecies$X.donor_species, ".donor", 2)[,1]
allspecies=merge(allspecies,clonal[,c(2,3,6)],by='donor_species',
                 by.y='cluster',
                 all.x=T)
allspecies_geneonly = allspecies2[which(!grepl('_other',allspecies2$gene)),]
allspecies$avg_gene_length = 0
for(i in 1:nrow(allspecies))
{
  allspecies_sub2 = allspecies_geneonly[which(allspecies_geneonly$X.donor_species %in%
                                                as.matrix(allspecies$X.donor_species[i])),]
  if (nrow(allspecies_sub2)>=4)
    avg_gene_length = sum(allspecies_sub2$gene_length)/nrow(allspecies_sub2)
  else
    avg_gene_length = sum(allspecies_geneonly$gene_length)/nrow(allspecies_geneonly)
  
  allspecies$avg_gene_length[i]=avg_gene_length
}
write.table(allspecies,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.txt',
            quote=F,sep='\t',row.names=F)

# significant calculating for species and lineage
allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.txt',header=T)
allspecies=allspecies[which(allspecies$No.SNP>0),]
allspecies$HS_pop[is.na(allspecies$HS_pop)]=0

quantile_cutoff = 0.01
# return the probability of higher than HS_NO
# The quantile is defined as the smallest value x SNPs per gene 
# such that F(x) > p, where F is the distribution function (accumulative).
quantile_hyper <- function(HS_NO,m, n, k, start_time = 0, startpoint=0){
  quantile_test =1- 2^-(startpoint:(startpoint+1000))
  temp_test = qhyper(quantile_test, m, n, k)
  small_list = which(temp_test <= HS_NO)
  endpoint = small_list[length(small_list)]
  if(endpoint==1001 & start_time <= 5)
  {
  startpoint = startpoint+1000
  start_time = start_time + 1
  quantile_hyper(HS_NO,m, n, k, start_time, startpoint)
  }
  else
  return(1-quantile_test[endpoint])
}
temp_simulate_sum2=matrix(0,nrow=0,ncol=14)
for (i in 1:nrow(allspecies))
{
  lineage = allspecies$X.donor_species[i]
  HS_NO = allspecies$HS_NO[i]
  HS_pop = allspecies$HS_pop[i]
  HS_NO_species = HS_NO - HS_pop 
  species = allspecies$species[i]
  SNP_set = allspecies$SNP_set[i]
  allspecies_sub = allspecies[which(allspecies$species==species),]
  Small_lineage = 'F'
  Small_species = 'NA'
  # lineage
  m = as.integer(allspecies$No.SNP[i]) #SNP bp
  total_gene = allspecies$tota_gene_num[i]
  n = as.integer(allspecies$tota_gene_num[i]*allspecies$avg_gene_length[i]-m)# no-SNP bp
  k = as.integer(allspecies$avg_gene_length[i])# gene length
  if (HS_pop == 0)
  {if(SNP_set < 2 | quantile_hyper(0,m, n, k, 1, 0)<=quantile_cutoff)
    Small_lineage = 'T'
  }
  quantile_out = quantile_hyper(HS_pop,m, n, k, 1, 0)
  # species
  if(nrow(allspecies_sub)>1)
  {
    Small_species = 'F'
    m2 = sum(as.integer(allspecies_sub$No.SNP))
    n2 = as.integer(sum((allspecies_sub$tota_gene_num)/(allspecies_sub$cluster_ratio)*
                        allspecies_sub$avg_gene_length-allspecies_sub$No.SNP))
  k2 = as.integer(mean(allspecies_sub$avg_gene_length))
  if (HS_NO_species == 0)
  {if(quantile_hyper(0,m2, n2, k2, 1, 0)<=quantile_cutoff)
    Small_species = 'T'
  }
    quantile_out_species = quantile_hyper(HS_NO_species,m2, n2, k2, 1, 0)
}  
  temp_result2 = data.frame(lineage = lineage,
    species=species,
    quantile_lineage = quantile_out,
    quantile_species = quantile_out_species,
    HS_lineage = HS_pop,
    HS_species = HS_NO_species,
    Small_lineage=Small_lineage,
    Small_species=Small_species,
    SNP_set = SNP_set,
    No.SNP = allspecies$No.SNP[i],
    No.gene_SNP = allspecies$No.gene_SNPs[i],
    HS_SNP_NO = allspecies$HS_SNP_NO[i],
    dNdS_cut = allspecies$dNdS_cut_highselect[i],
    expected_ratio = allspecies$expected_ratio[i])
  temp_simulate_sum2=rbind(temp_simulate_sum2,
                           temp_result2)
}

allspecies2=read.delim('final_SNP/HSgene/all.species.all.species.lineagenum.abu.dnds.sim.txt',header=T)
allspecies2=allspecies2[,c(1,14:19,21:26)]
temp_simulate_sum2=data.frame(temp_simulate_sum2)
temp_simulate_sum2=merge(temp_simulate_sum2,allspecies2,by='species')

# add simulation results
temp_simulate = read.delim('final_SNP/HSgene/HS_simulation.total.afterturningcutoff.clusterratio.withSNPs.lineage.txt',
                           header=T)
temp_simulate$adaptive_genes = temp_simulate$adaptive_genes_species+temp_simulate$adaptive_genes_lineage
temp_simulate$sim_N = temp_simulate$adaptive_genes_species_N + temp_simulate$adaptive_genes_lineage_N
temp_simulate$sim_S = temp_simulate$adaptive_genes_species_S + temp_simulate$adaptive_genes_lineage_S
temp_simulate$sim_S_cut = temp_simulate$sim_S
temp_simulate$sim_S_cut[which(temp_simulate$sim_S_cut==0)]=0.5
temp_simulate2 = matrix(0,nrow = 0,ncol= 11)
for(donor_species in unique(temp_simulate$X.donor_species)){
  temp_simulate_sub = temp_simulate[which(temp_simulate$X.donor_species == donor_species),]
  temp_simulate2 = rbind(temp_simulate2,
                         c(donor_species,
                           quantile(temp_simulate_sub$adaptive_genes,c(0.5)),
                           quantile(temp_simulate_sub$adaptive_genes_species_SNPs+
                                      temp_simulate_sub$adaptive_genes_lineage_SNPs,c(0.5,0.95,0.99)),
                           quantile(temp_simulate_sub$sim_N/temp_simulate_sub$sim_S_cut,c(0.5)),
                           mean(temp_simulate_sub$N_species+temp_simulate_sub$S_species),
                           mean(temp_simulate_sub$N+temp_simulate_sub$S),
                           mean(temp_simulate_sub$total_gene_species),
                           mean(temp_simulate_sub$avg_gene),
                           mean(temp_simulate_sub$avg_gene_length)
                         ) )
}

colnames(temp_simulate2)=c('X.donor_species','sim_gene_median',
                           'sim_SNP_median','sim_SNP_higher','sim_SNP_higher0.99','sim_NS_median',
                           'No.NS_species',
                           'No.NS','total_gene_species','total_gene','avg_gene_length')
temp_simulate_sum2=merge(temp_simulate_sum2,temp_simulate2,
                         by.x = 'lineage',by.y='X.donor_species',all.x=T)
write.table(temp_simulate_sum2,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.txt',
            quote=F,sep='\t',row.names=F)

# quality control
temp_simulate_sum2=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.txt',header=T)
temp_simulate_sum2$Small_species[is.na(temp_simulate_sum2$Small_species)]=FALSE
temp_simulate_sum2$HSratio = temp_simulate_sum2$HS_SNP_NO/temp_simulate_sum2$No.SNP
temp_simulate_sum2$HSratio_gene = (temp_simulate_sum2$HS_lineage+
                                     temp_simulate_sum2$HS_species
)/temp_simulate_sum2$No.gene_SNP
temp_simulate_sum2$dNdS_cut[is.na(temp_simulate_sum2$dNdS_cut)]=0
temp_simulate_sum2$HSratio_diff = temp_simulate_sum2$HSratio - temp_simulate_sum2$sim_SNP_median/temp_simulate_sum2$No.SNP
temp_simulate_sum2$HSratio_gene_diff = temp_simulate_sum2$HSratio_gene - temp_simulate_sum2$sim_gene_median/temp_simulate_sum2$No.gene_SNP
temp_simulate_sum2$dNdS_diff = temp_simulate_sum2$dNdS_cut - temp_simulate_sum2$sim_NS_median/temp_simulate_sum2$expected_ratio

colnames(temp_simulate_sum2)
temp_simulate_sum2$Small = FALSE
temp_simulate_sum2$Small[which(temp_simulate_sum2$HSratio == 0 & 
                                 (temp_simulate_sum2$sim_SNP_higher==0 | 
                                    temp_simulate_sum2$SNP_set == 1) )]=TRUE
count(temp_simulate_sum2,'Small')
write.table(temp_simulate_sum2,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.quality.txt',
            quote=F,sep='\t',row.names=F)

temp_simulate_sum2=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.quality.txt',header=T)
temp_simulate_sum3 = temp_simulate_sum2[which(temp_simulate_sum2$Small == FALSE &
                                                temp_simulate_sum2$HSratio == 0 &
                                                temp_simulate_sum2$sim_SNP_higher >= 1 &
                                                temp_simulate_sum2$SNP_set >= 3),]
temp_simulate_sum4 = temp_simulate_sum2[which(temp_simulate_sum2$species %in% temp_simulate_sum3$species),]
species_withHS = temp_simulate_sum4$species[which(temp_simulate_sum4$HSratio>0)]
temp_simulate_sum3 = temp_simulate_sum3[which(!(temp_simulate_sum3$species %in% species_withHS)),]
temp_simulate_sum4 = temp_simulate_sum2[which(temp_simulate_sum2$species %in% temp_simulate_sum3$species),]
write.table(temp_simulate_sum3,'final_SNP/HSgene/lineage_confident_noHS.txt',
            quote=F,sep='\t',row.names=F)
temp_simulate_sum4 = temp_simulate_sum2[which(!(temp_simulate_sum2$species 
                                              %in% temp_simulate_sum3$species)),]
temp_simulate_sum3 = temp_simulate_sum4[which(temp_simulate_sum4$HSratio > 0),]
length(unique(temp_simulate_sum3$species))
temp_simulate_sum4 = temp_simulate_sum4[which(!(temp_simulate_sum4$species %in% temp_simulate_sum3$species)),]
length(unique(temp_simulate_sum4$species))
length(unique(temp_simulate_sum2$species))

# sum up in species level
temp_simulate_sum2=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.quality.txt',header=T)
temp_simulate_sum2=temp_simulate_sum2[which(!is.na(temp_simulate_sum2$logabu.x)),]
temp_simulate_sum2$Small_species[is.na(temp_simulate_sum2$Small_species)]=FALSE
temp_simulate_sum2$HSratio = temp_simulate_sum2$HS_SNP_NO/temp_simulate_sum2$No.SNP
temp_simulate_sum2$HSratio_gene = (temp_simulate_sum2$HS_lineage+
                                     temp_simulate_sum2$HS_species
)/temp_simulate_sum2$No.gene_SNP
temp_simulate_sum2$dNdS_cut[is.na(temp_simulate_sum2$dNdS_cut)]=0
temp_simulate_sum2$HSratio_diff = temp_simulate_sum2$HSratio - temp_simulate_sum2$sim_SNP_median/temp_simulate_sum2$No.SNP
temp_simulate_sum2$HSratio_gene_diff = temp_simulate_sum2$HSratio_gene - temp_simulate_sum2$sim_gene_median/temp_simulate_sum2$No.gene_SNP
temp_simulate_sum2$dNdS_diff = temp_simulate_sum2$dNdS_cut - temp_simulate_sum2$sim_NS_median/temp_simulate_sum2$expected_ratio
temp_simulate_sum_sum = matrix(0,nrow=0,ncol=17)
colnames(temp_simulate_sum_sum)=c('species',
                                  'tag',
                                  'quantile_lineage_median',
                                  'quantile_lineage_lower',
                                  'quantile_lineage_higher',
                                  'quantile_species_median',
                                  'quantile_species_lower',
                                  'quantile_species_higher',
                                  'HS_ratio_diff_median',
                                  'HS_ratio_diff_lower',
                                  'HS_ratio_diff_higher',
                                  'HS_ratio_gene_diff_median',
                                  'HS_ratio_gene_diff_lower',
                                  'HS_ratio_gene_diff_higher',
                                  'dNdS_diff_median',
                                  'dNdS_diff_lower',
                                  'dNdS_diff_higher')
allspecies_set = unique(temp_simulate_sum2$species)
for(species in allspecies_set)
{
  temp_simulate_sum2_sub = temp_simulate_sum2[which(temp_simulate_sum2$species==species),]
  for(tag in unique(temp_simulate_sum2_sub$Small))
  {
    temp_simulate_sum_sum = rbind(temp_simulate_sum_sum,
                                  c(species,
                                    tag,
                                    quantile(temp_simulate_sum2_sub$quantile_lineage[
                                      which(temp_simulate_sum2_sub$Small == tag)],
                                      c(0.5,0.25,0.75)),
                                    quantile(temp_simulate_sum2_sub$quantile_species[
                                      which(temp_simulate_sum2_sub$Small == tag)],
                                      c(0.5,0.25,0.75)),
                                    quantile(temp_simulate_sum2_sub$HSratio_diff[
                                      which(temp_simulate_sum2_sub$Small == tag)],
                                      c(0.5,0.25,0.75)),
                                    quantile(temp_simulate_sum2_sub$HSratio_gene_diff[
                                      which(temp_simulate_sum2_sub$Small == tag)],
                                      c(0.5,0.25,0.75)),
                                    quantile(temp_simulate_sum2_sub$dNdS_diff[
                                      which(temp_simulate_sum2_sub$Small == tag)],
                                      c(0.5,0.25,0.75))
                                  )
    )
  }
}
temp_simulate_sum_sum = data.frame(temp_simulate_sum_sum)
colnames(temp_simulate_sum2)
temp_simulate_sum_sum=merge(temp_simulate_sum_sum,
                            temp_simulate_sum2[!duplicated(temp_simulate_sum2$species),
                                               c(2,15:26)],
                            by='species',all.x=T)

write.table(temp_simulate_sum_sum,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.species.txt',
            quote=F,sep='\t',row.names=F)

#----random forest preparation, sample size to max(freq)*2 figure 1b prep----
# prep
temp_simulate_sum2=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.quality.txt',header=T)
temp_simulate_sum2$Small = FALSE
temp_simulate_sum2$Small[which(temp_simulate_sum2$HSratio == 0 & 
                                 (temp_simulate_sum2$sim_SNP_higher==0 | 
                                    temp_simulate_sum2$SNP_set == 1) )]=TRUE
count(temp_simulate_sum2,'Small')

temp_simulate_sum2=temp_simulate_sum2[which(temp_simulate_sum2$Small!='TRUE'),]
temp_simulate_sum2_count = count(temp_simulate_sum2,c('species'))
temp_simulate_sum2_count=temp_simulate_sum2_count[order(temp_simulate_sum2_count$freq),]
maxfreq = max(temp_simulate_sum2_count$freq)*2
temp_simulate_sum2_count$freqnew = as.integer(maxfreq/(temp_simulate_sum2_count$freq))
temp_simulate_sum3_sub = matrix(0,
                                ncol = ncol(temp_simulate_sum2),
                                nrow = sum(temp_simulate_sum2_count$freqnew*temp_simulate_sum2_count$freq))
colnames(temp_simulate_sum3_sub)=colnames(temp_simulate_sum2)
k = 1
for(i in 1:nrow(temp_simulate_sum2_count))
{
  copy_time = temp_simulate_sum2_count$freqnew[i]
  for(j in 1:copy_time)
  {
    newmatrix = as.matrix(temp_simulate_sum2[
      which(temp_simulate_sum2$species==toString(temp_simulate_sum2_count$species[i])
      ),])
    m = nrow(newmatrix)
    temp_simulate_sum3_sub[c(k:(k+m-1)),]=newmatrix
  k = k + m
  }
}
temp_simulate_sum3_sub=data.frame(temp_simulate_sum3_sub)
write.table(temp_simulate_sum3_sub,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted.txt',
            quote=F,sep='\t',row.names=F)

# add sim_HS
accumulate_genenum <- function(hyper_lineage2,hyper_lineage)
{
  for(i in 1:nrow(hyper_lineage2))
  {
    sub_hyper_lineage = hyper_lineage[which(hyper_lineage$species == hyper_lineage2$species[i] &
                                              hyper_lineage$No_SNP_gene >=  hyper_lineage2$No_SNP_gene[i]),]
    hyper_lineage2$gene_num[i] = sum(sub_hyper_lineage$gene_num)
    hyper_lineage2$SNP_num[i] = sum(sub_hyper_lineage$gene_num*sub_hyper_lineage$No_SNP_gene)
  }
  return(hyper_lineage2)
}
accumulate_genenum_lineage <- function(hyper_lineage2,hyper_lineage)
{
  for(i in 1:nrow(hyper_lineage2))
  {
    sub_hyper_lineage = hyper_lineage[which(hyper_lineage$X.donor_species == 
                                              hyper_lineage2$X.donor_species[i] &
                                              hyper_lineage$No_SNP_gene >=  
                                              hyper_lineage2$No_SNP_gene[i]),]
    hyper_lineage2$gene_num[i] = sum(sub_hyper_lineage$gene_num)
    hyper_lineage2$SNP_num[i] = sum(sub_hyper_lineage$gene_num*sub_hyper_lineage$No_SNP_gene)
  }
  return(hyper_lineage2)
}

hyper_lineage=read.delim('final_SNP/HSgene/HS_hypergeometric_distribution.total.lineage.txt',header=T)
hyper_species=read.delim('final_SNP/HSgene/HS_hypergeometric_distribution.total.txt',header=T)
species_cutoff = read.delim('total_SNP_cutoff_species.txt',header=F)
hyper_lineage2 = hyper_lineage[which(hyper_lineage$No_SNP_gene==2),]
hyper_lineage2$SNP_num = 0
hyper_lineage2 = accumulate_genenum_lineage(hyper_lineage2,hyper_lineage)
hyper_species2 = hyper_species[which(hyper_species$No_SNP_gene==3 & 
                                       hyper_species$species %in% species_cutoff$V1),]
hyper_species2=rbind(hyper_species2,
                     hyper_species[which(hyper_species$No_SNP_gene==2 & 
                                           !(hyper_species$species %in% species_cutoff$V1)),])
hyper_species2$SNP_num = 0
hyper_species2 = accumulate_genenum(hyper_species2,hyper_species)

hyper_lineage = hyper_lineage2
hyper_lineage = merge(hyper_lineage,hyper_species2[,c(1,5,8)],
                      by='species')
hyper_lineage$sim_HS_lineage = hyper_lineage$gene_num.x
hyper_lineage$sim_SNP_lineage = hyper_lineage$SNP_num.x
hyper_lineage$sim_HS_species = hyper_lineage$gene_num.y
hyper_lineage$sim_SNP_species = hyper_lineage$SNP_num.y

temp_simulate_sum3_sub=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted.txt',header=T)
temp_simulate_sum3_sub=merge(temp_simulate_sum3_sub,
                             hyper_lineage[,c(2,10:13)],by.x='lineage',by.y='X.donor_species')
temp_simulate=read.delim('final_SNP/HSgene/HS_simulation.total.txt',header=T)
temp_simulate=data.frame(temp_simulate)
temp_simulate_sum = temp_simulate[!duplicated(temp_simulate$species),c(1,9)]
temp_simulate_sum3_sub=merge(temp_simulate_sum3_sub,temp_simulate_sum,by='species')

temp_simulate_sum3_sub$sim_HS_species[which(temp_simulate_sum3_sub$total_lineage_num==1)]=0
temp_simulate_sum3_sub$sim_SNP_species[which(temp_simulate_sum3_sub$total_lineage_num==1)]=0
temp_simulate_sum3_sub$HS_gene = temp_simulate_sum3_sub$HS_lineage + temp_simulate_sum3_sub$HS_species
temp_simulate_sum3_sub$HS_norm_gene = temp_simulate_sum3_sub$HS_gene/(
  temp_simulate_sum3_sub$HS_gene+
    temp_simulate_sum3_sub$sim_HS_lineage+temp_simulate_sum3_sub$sim_HS_species)
temp_simulate_sum3_sub$HS_norm_gene_lineage = temp_simulate_sum3_sub$HS_lineage/(
  temp_simulate_sum3_sub$HS_lineage+
    temp_simulate_sum3_sub$sim_HS_lineage)
temp_simulate_sum3_sub$HS_norm_gene_species = temp_simulate_sum3_sub$HS_species/(
  temp_simulate_sum3_sub$HS_species+temp_simulate_sum3_sub$sim_HS_species)
temp_simulate_sum3_sub$HSratio_gene_diff_lineage = (
  temp_simulate_sum3_sub$HS_lineage- temp_simulate_sum3_sub$sim_HS_lineage
)/temp_simulate_sum3_sub$No.gene_SNP

temp_simulate_sum3_sub$HSratio_gene_diff_species = (
  temp_simulate_sum3_sub$HS_species- temp_simulate_sum3_sub$sim_HS_species
)/temp_simulate_sum3_sub$No.gene_SNP

write.table(temp_simulate_sum3_sub,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted2.txt',
            quote=F,sep='\t',row.names=F)

allspecies=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.txt',header=T)
temp_simulate_sum3_sub=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted2.txt',header=T)
temp_simulate_sum3_sub=temp_simulate_sum3_sub[!duplicated(temp_simulate_sum3_sub$lineage),]
temp_simulate_sum3_sub=merge(temp_simulate_sum3_sub,allspecies[,c(1,3)],
                             by.x='lineage',by.y='X.donor_species',all.x=T)
write.table(temp_simulate_sum3_sub,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted2.addgenome.txt',
            quote=F,sep='\t',row.names=F)
# re sampling for lineage and species separately
temp_simulate_sum3=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted2.txt',header=T)
temp_simulate_sum3_sub = temp_simulate_sum3[which(!is.na(temp_simulate_sum3$HS_norm_gene) &
                                                    temp_simulate_sum3$HS_norm_gene >0  ),]
temp_simulate_sum3_sub_lineage =  temp_simulate_sum3_sub[which((temp_simulate_sum3_sub$HS_norm_gene_lineage>0) | 
                                                                 temp_simulate_sum3_sub$sim_HS_lineage>=0.05 ),]
length(unique(temp_simulate_sum3_sub_lineage$lineage))

temp_simulate_sum2=temp_simulate_sum3_sub_lineage[!duplicated(temp_simulate_sum3_sub_lineage$lineage),]
temp_simulate_sum2_count = count(temp_simulate_sum2,c('species'))
temp_simulate_sum2_count=temp_simulate_sum2_count[order(temp_simulate_sum2_count$freq),]
maxfreq = max(temp_simulate_sum2_count$freq)*2
temp_simulate_sum2_count$freqnew = as.integer(maxfreq/(temp_simulate_sum2_count$freq))
temp_simulate_sum3_sub = matrix(0,
                                ncol = ncol(temp_simulate_sum2),
                                nrow = sum(temp_simulate_sum2_count$freqnew*temp_simulate_sum2_count$freq))
colnames(temp_simulate_sum3_sub)=colnames(temp_simulate_sum2)
k = 1
for(i in 1:nrow(temp_simulate_sum2_count))
{
  copy_time = temp_simulate_sum2_count$freqnew[i]
  for(j in 1:copy_time)
  {
    newmatrix = as.matrix(temp_simulate_sum2[
      which(temp_simulate_sum2$species==toString(temp_simulate_sum2_count$species[i])
      ),])
    m = nrow(newmatrix)
    temp_simulate_sum3_sub[c(k:(k+m-1)),]=newmatrix
    k = k + m
  }
}
temp_simulate_sum3_sub=data.frame(temp_simulate_sum3_sub)
write.table(temp_simulate_sum3_sub,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted.forlineage.txt',
            quote=F,sep='\t',row.names=F)

temp_simulate_sum3_sub = temp_simulate_sum3[which(!is.na(temp_simulate_sum3$HS_norm_gene) &
                                                    temp_simulate_sum3$HS_norm_gene >0  ),]

temp_simulate_sum3_sub_species=temp_simulate_sum3_sub[which((temp_simulate_sum3_sub$HS_norm_gene_species>0) | 
                                                              temp_simulate_sum3_sub$sim_HS_species>=0.05 ),]
length(unique(temp_simulate_sum3_sub_species$lineage))

temp_simulate_sum2=temp_simulate_sum3_sub_species[!duplicated(temp_simulate_sum3_sub_species$lineage),]
temp_simulate_sum2_count = count(temp_simulate_sum2,c('species'))
temp_simulate_sum2_count=temp_simulate_sum2_count[order(temp_simulate_sum2_count$freq),]
maxfreq = max(temp_simulate_sum2_count$freq)*2
temp_simulate_sum2_count$freqnew = as.integer(maxfreq/(temp_simulate_sum2_count$freq))
temp_simulate_sum3_sub = matrix(0,
                                ncol = ncol(temp_simulate_sum2),
                                nrow = sum(temp_simulate_sum2_count$freqnew*temp_simulate_sum2_count$freq))
colnames(temp_simulate_sum3_sub)=colnames(temp_simulate_sum2)
k = 1
for(i in 1:nrow(temp_simulate_sum2_count))
{
  copy_time = temp_simulate_sum2_count$freqnew[i]
  for(j in 1:copy_time)
  {
    newmatrix = as.matrix(temp_simulate_sum2[
      which(temp_simulate_sum2$species==toString(temp_simulate_sum2_count$species[i])
      ),])
    m = nrow(newmatrix)
    temp_simulate_sum3_sub[c(k:(k+m-1)),]=newmatrix
    k = k + m
  }
}
temp_simulate_sum3_sub=data.frame(temp_simulate_sum3_sub)
write.table(temp_simulate_sum3_sub,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted.forspecies.txt',
            quote=F,sep='\t',row.names=F)
#
#----run random forest figure 1b prep----
# calculate importance
library(randomForest)
temp_simulate_sum3=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted2.txt',header=T)
#temp_simulate_sum3_sub_lineage = read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted.forlineage.txt',header=T)
#temp_simulate_sum3_sub_species = read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted.forspecies.txt',header=T)
{
  fit <- randomForest(HS_norm_gene_lineage ~ logabu.x +prevalence.x + No.SNP + No.gene_SNP +SNP_set+total_lineage_num+  
                       avg_gene_length + total_gene, data=
                        temp_simulate_sum3, importance=TRUE)
  print(fit) # view results, Mean of squared residuals: 1.925882e-05, % Var explained: 99.96
  fit_im = importance(fit, scale = FALSE) # importance of each predictor
  fit_im[,1] = fit_im[,1]/max(fit_im[,1])
  fit_im[,2] = fit_im[,2]/max(fit_im[,2])
  colnames(fit_im)=paste('HS_norm_gene_lineage',colnames(fit_im))
  test_randomForest = fit_im
  
  fit <- randomForest(HS_norm_gene_species ~ logabu.x +prevalence.x + No.SNP + No.gene_SNP +SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, data=
                        temp_simulate_sum3[which(!is.na(temp_simulate_sum3$HS_norm_gene_species)),], importance=TRUE)
  print(fit) # view results, Mean of squared residuals: 3.064234e-05, % Var explained: 99.94
  fit_im = importance(fit, scale = FALSE) # importance of each predictor
  fit_im[,1] = fit_im[,1]/max(fit_im[,1])
  fit_im[,2] = fit_im[,2]/max(fit_im[,2])
  colnames(fit_im)=paste('HS_norm_gene_species',colnames(fit_im))
  test_randomForest = cbind(test_randomForest,fit_im)
  
  
  fit <- randomForest(HSratio_gene_diff_lineage ~ logabu.x +prevalence.x+
                        No.SNP + No.gene_SNP +
                        SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, 
                      data=temp_simulate_sum3, importance=TRUE)
  print(fit) # view results, Mean of squared residuals: 0.0001382924, % Var explained: 99.76
  # the prediction error on the out-of-bag portion of the data, %IncMSE
  # IncNodePurity, total decrease in node impurities from splitting on the variable, averaged over all trees
  fit_im = importance(fit, scale = FALSE) # importance of each predictor
  fit_im[,1] = fit_im[,1]/max(fit_im[,1])
  fit_im[,2] = fit_im[,2]/max(fit_im[,2])
  colnames(fit_im)=paste('HSratio_gene_diff_lineage',colnames(fit_im))
  test_randomForest = cbind(test_randomForest,fit_im)
  
  
  fit <- randomForest(HSratio_gene_diff_species ~ logabu.x +prevalence.x + No.SNP + No.gene_SNP +SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, data=temp_simulate_sum3[which(!is.na(temp_simulate_sum3$HS_norm_gene_species)),], importance=TRUE)
  print(fit) # view results, Mean of squared residuals: 0.0001478777, % Var explained: 99.75
  fit_im = importance(fit, scale = FALSE) # importance of each predictor
  fit_im[,1] = fit_im[,1]/max(fit_im[,1])
  fit_im[,2] = fit_im[,2]/max(fit_im[,2])
  colnames(fit_im)=paste('HSratio_gene_diff_species',colnames(fit_im))
  test_randomForest = cbind(test_randomForest,fit_im)
  
  
  fit <- randomForest(dNdS_cut ~ logabu.x +prevalence.x + No.SNP + No.gene_SNP +SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, data=temp_simulate_sum3, importance=TRUE)
  print(fit) # view results, Mean of squared residuals: 0.5757298, % Var explained: 89.9
  fit_im = importance(fit, scale = FALSE) # importance of each predictor
  fit_im[,1] = fit_im[,1]/max(fit_im[,1])
  fit_im[,2] = fit_im[,2]/max(fit_im[,2])
  colnames(fit_im)=paste('dNdS_cut',colnames(fit_im))
  test_randomForest = cbind(test_randomForest,fit_im)
  
  # double check
  #sim_HS_lineage, sim_HS_species, sim_NS_median
  
  fit <- randomForest(sim_HS_species ~ logabu.x +prevalence.x + No.SNP + No.gene_SNP +SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, data=temp_simulate_sum3[which(!is.na(temp_simulate_sum3$HS_norm_gene_species)),], importance=TRUE)
  print(fit) # view results, Mean of squared residuals: 0.125985, % Var explained: 97.34
  fit_im = importance(fit, scale = FALSE) # importance of each predictor
  fit_im[,1] = fit_im[,1]/max(fit_im[,1])
  fit_im[,2] = fit_im[,2]/max(fit_im[,2])
  colnames(fit_im)=paste('sim_HS_species',colnames(fit_im))
  test_randomForest = cbind(test_randomForest,fit_im)
  
  fit <- randomForest(sim_HS_lineage ~ logabu.x +prevalence.x + No.SNP + No.gene_SNP +SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, data=temp_simulate_sum3, importance=TRUE)
  print(fit) # view results, Mean of squared residuals: 0.125985, % Var explained: 97.34
  fit_im = importance(fit, scale = FALSE) # importance of each predictor
  fit_im[,1] = fit_im[,1]/max(fit_im[,1])
  fit_im[,2] = fit_im[,2]/max(fit_im[,2])
  colnames(fit_im)=paste('sim_HS_lineage',colnames(fit_im))
  test_randomForest = cbind(test_randomForest,fit_im)
  
  fit <- randomForest(sim_NS_median ~ logabu.x +prevalence.x + No.SNP + No.gene_SNP +SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, data=temp_simulate_sum3, importance=TRUE)
  print(fit) # view results, Mean of squared residuals: 0.007751538, % Var explained: 99.49
  fit_im = importance(fit, scale = FALSE) # importance of each predictor
  fit_im[,1] = fit_im[,1]/max(fit_im[,1])
  fit_im[,2] = fit_im[,2]/max(fit_im[,2])
  colnames(fit_im)=paste('sim_NS_median',colnames(fit_im))
  test_randomForest = cbind(test_randomForest,fit_im)
}
write.table(test_randomForest,
            'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted.randomforest.importance.txt',
            quote=F,sep='\t',row.names=F)
# testing accuracy
allspecies_set = unique(temp_simulate_sum3$lineage)
total_species = length(allspecies_set)
accuracy_model = matrix(0,nrow = total_species*8,
                        ncol = 4)
k = 1
for(m in c(1:(total_species)))
{
  test_time = allspecies_set[m]
  test_time_all = toString(test_time)
  # subset
  test_set = which(temp_simulate_sum3$lineage == test_time)
  temp_simulate_sum3_test = temp_simulate_sum3[
    test_set,]
  temp_simulate_sum3_training = temp_simulate_sum3[
    setdiff(c(1:nrow(temp_simulate_sum3)),test_set) ,]
  # test model
  #HS_norm_gene_lineage
  fit <- randomForest(HS_norm_gene_lineage ~ logabu.x +prevalence.x+
                        No.SNP + No.gene_SNP +
                        SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, 
                      data=temp_simulate_sum3_training, importance=TRUE)
  temp_simulate_sum3_test$predict = predict(fit, temp_simulate_sum3_test)
  temp_simulate_sum3_test$accuracy = (temp_simulate_sum3_test$HS_norm_gene_lineage-
                                        temp_simulate_sum3_test$predict)^2
  accuracy_model[k,]=c(test_time_all,
                       mean(temp_simulate_sum3_test$accuracy)^0.5,
                       mean(fit$mse)^0.5,
                       'HS_norm_gene_lineage')
  k = k + 1
  #HS_norm_gene_species
  if(!is.na(temp_simulate_sum3_test$HS_norm_gene_species[1]))
  { fit <- randomForest(HS_norm_gene_species ~ logabu.x +prevalence.x+
                          No.SNP + No.gene_SNP +
                          SNP_set+total_lineage_num+  
                          avg_gene_length + total_gene, 
                        data=temp_simulate_sum3_training[which(!is.na(temp_simulate_sum3_training$HS_norm_gene_species)),], importance=TRUE)
  temp_simulate_sum3_test$predict = predict(fit, temp_simulate_sum3_test)
  temp_simulate_sum3_test$accuracy = (temp_simulate_sum3_test$HS_norm_gene_species-
                                        temp_simulate_sum3_test$predict)^2
  accuracy_model[k,]=c(test_time_all,
                       mean(temp_simulate_sum3_test$accuracy)^0.5,
                       mean(fit$mse)^0.5,
                       'HS_norm_gene_species')
  k = k + 1}
  #HSratio_gene_diff_lineage
  fit <- randomForest(HSratio_gene_diff_lineage ~ logabu.x +prevalence.x+
                        No.SNP + No.gene_SNP +
                        SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, 
                      data=temp_simulate_sum3_training, importance=TRUE)
  temp_simulate_sum3_test$predict = predict(fit, temp_simulate_sum3_test)
  temp_simulate_sum3_test$accuracy = (temp_simulate_sum3_test$HSratio_gene_diff_lineage-
                                        temp_simulate_sum3_test$predict)^2
  accuracy_model[k,]=c(test_time_all,
                       mean(temp_simulate_sum3_test$accuracy)^0.5,
                       mean(fit$mse)^0.5,
                       'HSratio_gene_diff_lineage')
  k = k + 1
  #HSratio_gene_diff_species
  fit <- randomForest(HSratio_gene_diff_species ~ logabu.x +prevalence.x+
                        No.SNP + No.gene_SNP +
                        SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, 
                      data=temp_simulate_sum3_training, importance=TRUE)
  temp_simulate_sum3_test$predict = predict(fit, temp_simulate_sum3_test)
  temp_simulate_sum3_test$accuracy = (temp_simulate_sum3_test$HSratio_gene_diff_species-
                                        temp_simulate_sum3_test$predict)^2
  accuracy_model[k,]=c(test_time_all,
                       mean(temp_simulate_sum3_test$accuracy)^0.5,
                       mean(fit$mse)^0.5,
                       'HSratio_gene_diff_species')
  k = k + 1
  #dNdS_cut
  fit <- randomForest(dNdS_cut ~ logabu.x +prevalence.x+
                        No.SNP + No.gene_SNP +
                        SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, 
                      data=temp_simulate_sum3_training, importance=TRUE)
  temp_simulate_sum3_test$predict = predict(fit, temp_simulate_sum3_test)
  temp_simulate_sum3_test$accuracy = (temp_simulate_sum3_test$dNdS_cut-
                                        temp_simulate_sum3_test$predict)^2
  accuracy_model[k,]=c(test_time_all,
                       mean(temp_simulate_sum3_test$accuracy)^0.5,mean(fit$mse)^0.5,
                       'dNdS_cut')
  k = k + 1
  #sim_HS_lineage
  fit <- randomForest(sim_HS_lineage ~ logabu.x +prevalence.x+
                        No.SNP + No.gene_SNP +
                        SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, 
                      data=temp_simulate_sum3_training, importance=TRUE)
  temp_simulate_sum3_test$predict = predict(fit, temp_simulate_sum3_test)
  temp_simulate_sum3_test$accuracy = (temp_simulate_sum3_test$sim_HS_lineage-
                                        temp_simulate_sum3_test$predict)^2
  accuracy_model[k,]=c(test_time_all,
                       mean(temp_simulate_sum3_test$accuracy)^0.5,mean(fit$mse)^0.5,
                       'sim_HS_lineage')
  k = k + 1
  #sim_HS_species
  fit <- randomForest(sim_HS_species ~ logabu.x +prevalence.x+
                        No.SNP + No.gene_SNP +
                        SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, 
                      data=temp_simulate_sum3_training, importance=TRUE)
  temp_simulate_sum3_test$predict = predict(fit, temp_simulate_sum3_test)
  temp_simulate_sum3_test$accuracy = (temp_simulate_sum3_test$sim_HS_species-
                                        temp_simulate_sum3_test$predict)^2
  accuracy_model[k,]=c(test_time_all,
                       mean(temp_simulate_sum3_test$accuracy)^0.5,mean(fit$mse)^0.5,
                       'sim_HS_species')
  k = k + 1
  #sim_NS_median
  fit <- randomForest(sim_NS_median ~ logabu.x +prevalence.x+
                        No.SNP + No.gene_SNP +
                        SNP_set+total_lineage_num+  
                        avg_gene_length + total_gene, 
                      data=temp_simulate_sum3_training, importance=TRUE)
  temp_simulate_sum3_test$predict = predict(fit, temp_simulate_sum3_test)
  temp_simulate_sum3_test$accuracy = (temp_simulate_sum3_test$sim_NS_median-
                                        temp_simulate_sum3_test$predict)^2
  accuracy_model[k,]=c(test_time_all,
                       mean(temp_simulate_sum3_test$accuracy)^0.5,mean(fit$mse)^0.5,
                       'sim_NS_median')
  k = k + 1
}
colnames(accuracy_model)=c(
  'species','mse_squared_test',
  'mse_squared_training','tag'
)
write.table(accuracy_model,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted.randomforest.accuracy.txt',
            quote=F,sep='\t',row.names=F)

#----plot random forest figure 1b----
library(pheatmap)
temp_simulate_sum3_sub=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted2.txt',header=T)
temp_simulate_sum3_sub = temp_simulate_sum3[which(!is.na(temp_simulate_sum3$HS_norm_gene) &
                                                    temp_simulate_sum3$HS_norm_gene >0  ),]
accuracy_model=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted.randomforest.accuracy.txt',
                          header=T)
allvariable = unique(accuracy_model$tag)[c(1:5,7,6,8)]
mat = matrix(0,nrow = 1, ncol = length(allvariable))
colnames(mat)=allvariable
k = 1
for(tag in allvariable)
{
  range_tag = temp_simulate_sum3_sub[,which(colnames(temp_simulate_sum3_sub)==tag)]
  range_tag = range_tag[which(!is.na(range_tag))]
  print(tag)
  # error rate
  error_rate = mean(accuracy_model$mse_squared_test[
    which(accuracy_model$tag==tag)
    ])/(max(range_tag)-min(range_tag))
  print(error_rate
    )
  mat[1,k]=error_rate
  k = k+1
}
mat = 1-mat
mat = rbind(mat,
            c(97.83,99.37,99.94,92.81,90.77,99.97,97,99.55)/100)
row.names(mat)=c('Accuracy','Var explained')
mat_accuracy = mat
write.table(mat_accuracy,'final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted.randomforest.accuracy.sum.txt',
            quote=F,sep='\t',row.names=T)

test_randomForest = read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.adjusted.randomforest.importance.txt',
                               header=T)
row.names(test_randomForest)=
  c('abundance (log10)','prevalence (%)',
    'No. SNPs per pop','No. gene with SNPs per pop',
    'No. unique genotypes',"No. lineages per species",
   'Avg gene length','No. genes per pop')

subsettag = c(1,3,5,7,9,11,13,15)
mat = test_randomForest[,subsettag]
row.names(mat)=row.names(test_randomForest)
mat=as.matrix(mat)
colnames(mat)=c('Norm PE within a lineage',
                'Norm PE across lineages',
                'Norm PE ratio within a lineage',
                'Norm PE ratio across lineages',
                'Observed dN/dS of PE',
                'Expected PE across lineages',
                'Expected PE within a lineage',
                'Expected dN/dS of PE'
)
Metadata2=data.frame(
  Parameters = c('ecology',
             'ecology',
             'genome collection',
             'genome collection',
             'genome collection',
             'genome collection',
             'population genetics',
             'population genetics'
             )
  )
row.names(Metadata2)=row.names(test_randomForest)
mat_colors <- list(Parameters = brewer.pal(9, "Set1")[c(1,5,6)])
names(mat_colors$Parameters) <- unique(Metadata2$Parameters)

color_set = colorRampPalette(brewer.pal(n = 9, name = "Purples"))(9)[c(1:9)]
neworder <- pheatmap(mat= mat,
         color = color_set,
         show_rownames= TRUE,
         show_colnames = TRUE,
         cluster_cols = TRUE,
         annotation_row    = Metadata2,
         annotation_colors = mat_colors,
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         cluster_rows = TRUE)

breaks_mat = seq(0.5,1,0.05)
color_set = colorRampPalette(brewer.pal(n = 9, name = "OrRd"))(9)[c(1:9)]
mat_accuracy = mat_accuracy[,neworder$tree_col$order]
pheatmap(mat= mat_accuracy,
         color = color_set,
         breaks = breaks_mat,
         show_rownames= TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE)

mat = test_randomForest[,subsettag+1]
row.names(mat)=row.names(test_randomForest)
mat=as.matrix(mat)
colnames(mat)=c('Norm PE within a lineage',
                'Norm PE across lineages',
                'Norm PE ratio within a lineage',
                'Norm PE ratio across lineages',
                'Observed dN/dS of PE',
                'Expected PE across lineages',
                'Expected PE within a lineage',
                'Expected dN/dS of PE'
)
color_set = colorRampPalette(brewer.pal(n = 9, name = "Purples"))(9)[c(1:9)]
pheatmap(mat= mat,
         color = color_set,
         show_rownames= TRUE,
         show_colnames = TRUE,
         cluster_cols = TRUE,
         annotation_row    = Metadata2,
         annotation_colors = mat_colors,
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         cluster_rows = TRUE)
#
#----examples ecology----
library(randomForest)
temp_simulate_sum2=read.delim('final_SNP/HSgene/all.species.High_select2.dmrca.all.short.lineage.genenum.allsum.withgenenum.HS_hypergeometric_distribution.quality.txt',header=T)
temp_simulate_sum2$No.SNP2 = as.integer(log10(temp_simulate_sum2$No.SNP))
temp_simulate_sum2$SNP_set2 = as.integer(log2(temp_simulate_sum2$SNP_set))
plot(log10(temp_simulate_sum2$No.SNP),log2(temp_simulate_sum2$SNP_set),
     col=temp_simulate_sum2$species)
# add simulation results
temp_simulate = read.delim('final_SNP/HSgene/HS_simulation.total.afterturningcutoff.clusterratio.withSNPs.lineage.txt',
                           header=T)
temp_simulate2 = matrix(0,nrow = 0,ncol= 3)
temp_simulate$adaptive_genes_lineage_SNPs
for(donor_species in unique(temp_simulate$X.donor_species)){
  temp_simulate_sub = temp_simulate[which(temp_simulate$X.donor_species == donor_species),]
  temp_simulate2 = rbind(temp_simulate2,
                         c(donor_species,
                           quantile(temp_simulate_sub$adaptive_genes_lineage,c(0.5)),
                           quantile(temp_simulate_sub$adaptive_genes_lineage_SNPs,c(0.5))
                         ) )
}

colnames(temp_simulate2)=c('X.donor_species','sim_gene_lineage_median',
                           'sim_SNP_lineage_median')
temp_simulate_sum2=merge(temp_simulate_sum2,temp_simulate2,
                         by.x = 'lineage',by.y='X.donor_species',all.x=T)

count_species = count(temp_simulate_sum2,c('No.SNP2','SNP_set2'))

SNP_range = c(50:120)
#SNP_set_range = c(3:2^5)
temp_simulate_sum2_example = temp_simulate_sum2[
  which(temp_simulate_sum2$No.SNP %in% SNP_range & 
          temp_simulate_sum2$SNP_set >=40 &
          temp_simulate_sum2$Small_lineage!='TRUE'),
  ]
# example 1
SNP_range1 = c(10:40)
SNP_set_range1 = c(8:12)
temp_simulate_sum2_example1 = temp_simulate_sum2[
  which(temp_simulate_sum2$No.SNP %in% SNP_range1 & 
          temp_simulate_sum2$SNP_set %in% SNP_set_range1 &
          temp_simulate_sum2$Small_lineage!='TRUE'&
          temp_simulate_sum2$lineage!='1_BL_IBD_0_clustercluster1.donor.D77'),
  ]

# example 2
SNP_range1 = c(115:250)
SNP_set_range1 = c(25:45)
temp_simulate_sum2_example = temp_simulate_sum2[
  which(temp_simulate_sum2$No.SNP %in% SNP_range1 & 
          temp_simulate_sum2$SNP_set %in% SNP_set_range1 
        & temp_simulate_sum2$Small_lineage!='TRUE'
        & !(temp_simulate_sum2$lineage %in%
        c('1_BA_IBD_2_clustercluster5.donor.D469','1_BA_IBD_1_clustercluster2.donor.D544'))
  ),
  ]

plot(log10(temp_simulate_sum2_example$No.SNP),temp_simulate_sum2_example$HS_lineage,
     col=temp_simulate_sum2_example$species)
plot(log2(temp_simulate_sum2_example$SNP_set),temp_simulate_sum2_example$HS_lineage,
     col=temp_simulate_sum2_example$species)
plot(temp_simulate_sum2_example$logabu.x,log10(temp_simulate_sum2_example$HS_lineage),
     col=temp_simulate_sum2_example$species)


####
