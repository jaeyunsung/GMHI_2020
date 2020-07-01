
# This script shows how to reproduce the results in Figure 6 of Gupta et. al (version June 2020)
# Author: Vinod K. Gupta, PhD

install.packages('vegan')
install.packages('ggplot2')
install.packages('ggpubr')
install.packages('rcompanion')

library(vegan)
library(ggplot2)
library(ggpubr)
library(rcompanion)

fig6<-read.csv("./validation_abundance.csv", sep = ",", header = TRUE,row.names = 1,check.names = F) # change path to the directory wherein this file is located

# Metadata for validation datasets
fig6_metadata<-read.csv("./validation_metadata.csv", sep = ",", header = TRUE,row.names = 1,check.names = F) # change path to the directory wherein this file is located
fig6_1 <- fig6[-grep('unclassified', row.names(fig6)),] # unclassified species removed
fig6_2 <- fig6_1[-grep('virus', row.names(fig6_1)),] # viruses removed
fig6_3<-sweep(fig6_2,2,colSums(fig6_2),`/`) # re-normalization to get relative abundances after removing unclassified and viruses
fig6_4<-(fig6_3)*100
fig6_4[fig6_4 < 0.001] <- 0

# GMHI calculation for validation datasets
H_species<-c("s__Alistipes_senegalensis", "s__Bacteroidales_bacterium_ph8",
             "s__Bifidobacterium_adolescentis", "s__Bifidobacterium_angulatum",
             "s__Bifidobacterium_catenulatum", "s__Lachnospiraceae_bacterium_8_1_57FAA",
             "s__Sutterella_wadsworthensis") # from Table 2
NH_species<-c("s__Anaerotruncus_colihominis", "s__Atopobium_parvulum", "s__Bifidobacterium_dentium",
              "s__Blautia_producta", "s__candidate_division_TM7_single_cell_isolate_TM7c",
              "s__Clostridiales_bacterium_1_7_47FAA", "s__Clostridium_asparagiforme",
              "s__Clostridium_bolteae", "s__Clostridium_citroniae", "s__Clostridium_clostridioforme",
              "s__Clostridium_hathewayi", "s__Clostridium_nexile", "s__Clostridium_ramosum",
              "s__Clostridium_symbiosum", "s__Eggerthella_lenta", "s__Erysipelotrichaceae_bacterium_2_2_44A",
              "s__Flavonifractor_plautii", "s__Fusobacterium_nucleatum", "s__Gemella_morbillorum",
              "s__Gemella_sanguinis", "s__Granulicatella_adiacens", "s__Holdemania_filiformis",
              "s__Klebsiella_pneumoniae", "s__Lachnospiraceae_bacterium_1_4_56FAA",
              "s__Lachnospiraceae_bacterium_2_1_58FAA", "s__Lachnospiraceae_bacterium_3_1_57FAA_CT1",
              "s__Lachnospiraceae_bacterium_5_1_57FAA", "s__Lachnospiraceae_bacterium_9_1_43BFAA",
              "s__Lactobacillus_salivarius", "s__Peptostreptococcus_stomatis",
              "s__Ruminococcaceae_bacterium_D16", "s__Ruminococcus_gnavus",
              "s__Solobacterium_moorei", "s__Streptococcus_anginosus", "s__Streptococcus_australis",
              "s__Streptococcus_gordonii", "s__Streptococcus_infantis", "s__Streptococcus_mitis_oralis_pneumoniae",
              "s__Streptococcus_sanguinis", "s__Streptococcus_vestibularis",
              "s__Subdoligranulum_sp_4_3_54A2FAA", "s__Subdoligranulum_variabile",
              "s__Veillonella_atypica") # from Table 2
H_constant<-7 # Healthy normalization constant calculated based on training dataset
NH_constant<-31 # Nonhealthy normalization constant calculated based on training dataset
validation_H <- fig6_4[row.names(fig6_4) %in% H_species, ] # extracting H-plus species
validation_NH <- fig6_4[row.names(fig6_4) %in% NH_species, ] # extracting H-minus species
alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
H_shannon<-apply((validation_H/100), 2, alpha)
NH_shannon<-apply((validation_NH/100), 2, alpha)
H_sig_count<-apply(validation_H, 2, function(i) (sum(i > 0)))
NH_sig_count<-apply(validation_NH, 2, function(i) (sum(i > 0)))
H_GMHI<-((H_sig_count/H_constant)*H_shannon)
NH_GMHI<-((NH_sig_count/NH_constant)*NH_shannon)
GMHI<-data.frame(log10((H_GMHI+0.00001)/(NH_GMHI+0.00001)))
colnames(GMHI) <- c("GMHI")

# Shannon diversity calculation for validation datasets
alpha_shannon_fig6<-data.frame(diversity(fig6_4, 2, index="shannon"))
fig6_alpha_GMHI_data<-merge(as.data.frame(GMHI),as.data.frame(alpha_shannon_fig6), by='row.names', all=F)
colnames(fig6_alpha_GMHI_data)<-c("Sample_id","GMHI","alpha")
row.names(fig6_alpha_GMHI_data)<-fig6_alpha_GMHI_data$Sample_id
fig6_alpha_GMHI_data1<-merge(as.data.frame(fig6_alpha_GMHI_data),as.data.frame(fig6_metadata), by='row.names', all=F)
fig6_alpha_GMHI_data_final<-fig6_alpha_GMHI_data1[-grep(':OB|:OW|:UW', fig6_alpha_GMHI_data1$Phenotype_all),] # remove obesity, overweight and underweight samples
H_Un_count<-data.frame(table(fig6_alpha_GMHI_data_final$Phenotype)) # healthy_nonhealthy_validation_sample_counts
All_phenotype_counts<-data.frame(table(fig6_alpha_GMHI_data_final$Phenotype_all)) # phenotype_wise_validation_sample_counts
fig6_alpha_GMHI_data_final$Phenotype <- gsub(x = fig6_alpha_GMHI_data_final$Phenotype, pattern = "Unhealthy", replacement = "Nonhealthy")

# Figure 6a (GMHI for combined phenotype as Healthy vs Nonhealthy)
p_gmhi_violin<-ggplot(fig6_alpha_GMHI_data_final, aes(x=Phenotype, y=GMHI, fill=Phenotype)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
p_gmhi_violin+rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_color_manual(values=c('steelblue','orange2'))+stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+scale_fill_manual(values=c('steelblue','orange2'))+labs(x = "",y="Gut Microbiome Health Index (GMHI)",title = "Fig6a")+theme_bw()+theme(panel.grid = element_blank())+rremove("legend")
cliffDelta(data =fig6_alpha_GMHI_data_final,GMHI~Phenotype)

# Figure 6b (GMHI for all study-wise phenotypes)
fig6_GMHI_boxplot <- ggboxplot(fig6_alpha_GMHI_data_final, x = "Phenotype_all", y = "GMHI",
                               color = "Phenotype_all",
                               add = "jitter",
                               short.panel.labs = FALSE,outlier.size=0,xlab = "", ylab = "Gut Microbiome Health Index (GMHI)",title = "Fig6b")
fig6_GMHI_boxplot +rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_color_manual(values=c("#FB8072","steelblue","orange2","steelblue","orange2","orange2","#BEBADA","#BEBADA","#FFFFB3","#BEBADA","steelblue","#D9D9D9"))
Significant_GMHI_H1 <- compare_means(GMHI ~ Phenotype_all, data = fig6_alpha_GMHI_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "Healthy_P")
Significant_GMHI_H2 <- compare_means(GMHI ~ Phenotype_all, data = fig6_alpha_GMHI_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "N-1:Healthy")
Significant_GMHI_H3 <- compare_means(GMHI ~ Phenotype_all, data = fig6_alpha_GMHI_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "N-8:Healthy")

# Figure 6c (Shannon diversity for combined phenotype as Healthy vs Nonhealthy)
p_alpha_violin<-ggplot(fig6_alpha_GMHI_data_final, aes(x=Phenotype, y=alpha, fill=Phenotype)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
p_alpha_violin+rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_color_manual(values=c('steelblue','orange2'))+stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+scale_fill_manual(values=c('steelblue','orange2'))+labs(x = "",y="Shannon Diversity",title = "Fig6c")+theme_bw()+theme(panel.grid = element_blank())+rremove("legend")
cliffDelta(data =fig6_alpha_GMHI_data_final,alpha~Phenotype)

# Fig6d (Shannon diversity for all study-wise phenotypes)
fig6d <- ggboxplot(fig6_alpha_GMHI_data_final, x = "Phenotype_all", y = "alpha",
                                color = "Phenotype_all",
                                add = "jitter",
                                short.panel.labs = FALSE,outlier.size=0,xlab = "", ylab = "Shannon Diversity",title = "Fig6d")
fig6d +rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_color_manual(values=c("#FB8072","steelblue","orange2","steelblue","orange2","orange2","#BEBADA","#BEBADA","#FFFFB3","#BEBADA","steelblue","#D9D9D9"))
Significant_alpha_H1 <- compare_means(alpha ~ Phenotype_all, data = fig6_alpha_GMHI_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "Healthy_P")
Significant_alpha_H2 <- compare_means(alpha ~ Phenotype_all, data = fig6_alpha_GMHI_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "N-1:Healthy")
Significant_alpha_H3 <- compare_means(alpha ~ Phenotype_all, data = fig6_alpha_GMHI_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "N-8:Healthy")

# Figures 6e and 6f (80% abundance coverage for combined phenotype as Healthy vs Nonhealthy)
richness_fig6<-data.frame((apply(fig6_4>0, 2, sum)))
cum<-data.frame(apply(apply(fig6_4, 2, sort,decreasing=T), 2, cumsum))
cum_cutoff2 <- function(x){max(which(x < 80))+1}
prevalence2<-data.frame(apply(cum, 2, cum_cutoff2))
prevalence2[prevalence2 =="-Inf"] <- 1
fig6_richness_abcov_data<-merge(as.data.frame(richness_fig6),as.data.frame(prevalence2), by='row.names', all=F)
colnames(fig6_richness_abcov_data)<-c("Sample_id","richness","abun_cov")
row.names(fig6_richness_abcov_data)<-fig6_richness_abcov_data$Sample_id
fig6_richness_abcov_data1<-merge(as.data.frame(fig6_richness_abcov_data),as.data.frame(fig6_metadata), by='row.names', all=F)
fig6_richness_abcov_data_final<-fig6_richness_abcov_data1[-grep(':OB|:OW|:UW', fig6_richness_abcov_data1$Phenotype_all),] # remove obesity, overweight and underweight samples
p_gmhi_violin<-ggplot(fig6_richness_abcov_data_final, aes(x=Phenotype, y=abun_cov, fill=Phenotype)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
p_gmhi_violin+rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_color_manual(values=c('steelblue','orange2'))+stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+scale_fill_manual(values=c('steelblue','orange2'))+labs(x = "",y="80% Abundance Coverage (# of Species)")+theme_bw()+theme(panel.grid = element_blank())+rremove("legend")
cliffDelta(data =fig6_richness_abcov_data_final,abun_cov~Phenotype)
fig6f<- ggboxplot(fig6_richness_abcov_data_final, x = "Phenotype_all", y = "abun_cov",
                               color = "Phenotype_all",
                               add = "jitter",
                               short.panel.labs = FALSE,outlier.size=0,xlab = "", ylab = "80% Abundance Coverage (# of Species)")
fig6f +rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_color_manual(values=c("#FB8072","steelblue","orange2","steelblue","orange2","orange2","#BEBADA","#BEBADA","#FFFFB3","#BEBADA","steelblue","#D9D9D9"))
Significant_abundance_coverage_H1 <- compare_means(abun_cov ~ Phenotype_all, data = fig6_richness_abcov_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "Healthy_P")
Significant_abundance_coverage_H2 <- compare_means(abun_cov ~ Phenotype_all, data = fig6_richness_abcov_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "N-1:Healthy")
Significant_abundance_coverage_H3 <- compare_means(abun_cov ~ Phenotype_all, data = fig6_richness_abcov_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "N-8:Healthy")

# Figures 6g and 6h (species richness for combined phenotype as Healthy vs Nonhealthy)
p_alpha_violin<-ggplot(fig6_richness_abcov_data_final, aes(x=Phenotype, y=richness, fill=Phenotype)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
p_alpha_violin+rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_color_manual(values=c('steelblue','orange2'))+stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+scale_fill_manual(values=c('steelblue','orange2'))+labs(x = "",y="Species richness")+theme_bw()+theme(panel.grid = element_blank())+rremove("legend")
cliffDelta(data =fig6_richness_abcov_data_final,richness~Phenotype)
fig6h <- ggboxplot(fig6_richness_abcov_data_final, x = "Phenotype_all", y = "richness",
                                color = "Phenotype_all",
                                add = "jitter",
                                short.panel.labs = FALSE,outlier.size=0,xlab = "", ylab = "Species richness")
fig6h +rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_color_manual(values=c("#FB8072","steelblue","orange2","steelblue","orange2","orange2","#BEBADA","#BEBADA","#FFFFB3","#BEBADA","steelblue","#D9D9D9"))
Significant_richness_H1 <- compare_means(richness ~ Phenotype_all, data = fig6_richness_abcov_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "Healthy_P")
Significant_richness_H2 <- compare_means(richness ~ Phenotype_all, data = fig6_richness_abcov_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "N-1:Healthy")
Significant_richness_H3 <- compare_means(richness ~ Phenotype_all, data = fig6_richness_abcov_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "N-8:Healthy")

# End
