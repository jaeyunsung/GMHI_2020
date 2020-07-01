
# This script shows how to reproduce the results in Figure 5 of Gupta et. al (version June 2020)
# Author: Vinod K. Gupta, PhD

install.packages('dplyr')
install.packages('vegan')
install.packages('ggpubr')
install.packages('ggplot2')
install.packages('reshape2')

library(dplyr)
library(vegan)
library(ggpubr)
library(ggplot2)
library(reshape2)

Final_microbiome_data_4347<-read.csv("./Final_metadata_4347.csv", sep = ",", header = T,row.names = 1,check.names = F)
fig5_1<-data.frame(t(Final_microbiome_data_4347),check.rows = F,check.names = F)

# Studies having samples from healthy as well as from Nonhealthy individuals
study_wise_data<-c("S-7_Healthy",
                   "S-7_Crohns-disease",
                   "V-12_Healthy",
                   "V-12_Obesity",
                   "V-16_Crohns-disease",
                   "V-16_Overweight",
                   "V-16_Healthy",
                   "V-16_Ulcerative-colitis",
                   "V-16_Obesity",
                   "V-19_Obesity",
                   "V-19_Overweight",
                   "V-19_Healthy",
                   "V-19_advanced-adenoma",
                   "V-19_CRC",
                   "V-2_ACVD",
                   "V-2_Overweight",
                   "V-2_Healthy",
                   "V-29_Overweight",
                   "V-29_Healthy",
                   "V-29_Rheumatoid-Arthritis",
                   "V-43_Overweight",
                   "V-43_Healthy",
                   "V-43_CRC",
                   "V-46_Healthy",
                   "V-46_Obesity",
                   "V-48_Healthy",
                   "V-48_Ulcerative-colitis",
                   "V-48_Crohns-disease",
                   "V-6_Overweight",
                   "V-6_Healthy",
                   "V-6_CRC",
                   "V-6_advanced-adenoma",
                   "V-8_Underweight",
                   "V-8_Overweight",
                   "V-8_Healthy",
                   "V-8_T2D",
                   "V-9_Overweight",
                   "V-9_Healthy",
                   "V-9_IGT",
                   "V-9_T2D")
final_dataset<-fig5_1[fig5_1$study %in% study_wise_data,]
final_dataset1<-final_dataset[,c(1:2,4,15,33:ncol(final_dataset))]
colnames(final_dataset1)[2]<-"study1"
metadata<-data.frame(final_dataset1[,c(1:4)])
colnames(metadata)[3]<-"author"
metadata$author_phenotype = paste(metadata$author,"",metadata$Phenotype)
metadata$author_phenotype<-as.factor(metadata$author_phenotype)
final_dataset1[,-c(1:4)] <- lapply(final_dataset1[,-c(1:4)], function(x) as.numeric(as.character(x)))
final_dataset2<-data.frame(t(final_dataset1),check.rows = F,check.names = F)
final_dataset3<-data.frame(final_dataset2[-c(1:4),])
final_dataset3[] <- lapply(final_dataset3, function(x) as.numeric(as.character(x)))

# Figure 5a
final_dataset3[final_dataset3 < 0.001] <- 0
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
H_constant<-7 # from Figure 2
NH_constant<-31 # from Figure 2
sp_H <- final_dataset3[row.names(final_dataset3) %in% H_species, ]
sp_NH <- final_dataset3[row.names(final_dataset3) %in% NH_species, ]
alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
H_shannon<-apply((sp_H/100), 2, alpha)
NH_shannon<-apply((sp_NH/100), 2, alpha)
H_sig_count<-apply(sp_H, 2, function(i) (sum(i > 0)))
NH_sig_count<-apply(sp_NH, 2, function(i) (sum(i > 0)))
H_GMHI<-((H_sig_count/H_constant)*H_shannon)
NH_GMHI<-((NH_sig_count/NH_constant)*NH_shannon)
GMHI<-data.frame(log10((H_GMHI+0.00001)/(NH_GMHI+0.00001)))
colnames(GMHI) <- c("GMHI")
Result<-data.frame(cbind(GMHI,H_sig_count,NH_sig_count,H_shannon,NH_shannon,H_GMHI,NH_GMHI))
GMHI_final<-data.frame(metadata,Result)
colnames(GMHI_final)[6]<-"GMHI_final"
Significant_GMHI <- compare_means(GMHI_final ~ Phenotype, group.by = "author", data = GMHI_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.001,0.05,1), symbols = c("***", "*","ns")),ref.group = "Healthy")
GMHI_plot<-ggboxplot(GMHI_final, x = "author_phenotype", y = "GMHI_final",
                     color = 'black',fill = "Phenotype",
                     title = "Fig5a",
                     short.panel.labs = FALSE,outlier.shape=NA,xlab = "", ylab = "Gut Microbiome Health Index (GMHI)",ylim = c(-6, 6))
gp<-GMHI_plot+theme_linedraw()+scale_fill_brewer(palette="Set3")+rremove("legend")
gp + coord_flip()

# Figure 5b
alpha_shannon<-data.frame(metadata,(diversity(final_dataset1[,-c(1:4)], 1, index="shannon")))
colnames(alpha_shannon)[6]<-"alpha"
Significant_alpha <- compare_means(alpha ~ Phenotype, group.by = "author", data = alpha_shannon, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.001,0.05,1), symbols = c("***", "*","ns")),ref.group = "Healthy")
alpha_plot<-ggboxplot(alpha_shannon, x = "author_phenotype", y = "alpha",
                      color = 'black',title = "Fig5b",fill = "Phenotype",
                      short.panel.labs = FALSE,outlier.shape=NA,xlab = "", ylab = "Shannon Diversity",ylim = c(0, 5))
ap<-alpha_plot+theme_linedraw()+scale_fill_brewer(palette="Set3")+rremove("legend")
ap+coord_flip()

# Figure 5c
cum<-data.frame(apply(apply(final_dataset3, 2, sort,decreasing=T), 2, cumsum))
cum_cutoff2 <- function(x){max(which(x < 80))+1}
prevalence2<-data.frame(apply(cum, 2, cum_cutoff2))
prevalence2[prevalence2 =="-Inf"] <- 1
colnames(prevalence2)<-"prevalence"
prevalence<-data.frame(metadata,prevalence2$prevalence)
colnames(prevalence)[6]<-"prevalence"
Significant_prevalence <- compare_means(prevalence ~ Phenotype, group.by = "author", data = prevalence, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.001,0.05,1), symbols = c("***", "*","ns")),ref.group = "Healthy")
prevalence_plot<-ggboxplot(prevalence, x = "author_phenotype", y = "prevalence",
                           color = 'black',fill = "Phenotype",
                           title = "Fig5c",
                           short.panel.labs = FALSE,outlier.shape=NA,xlab = "", ylab = "80% abundance coverage (# of Species)",ylim = c(0, 40))
pp<-prevalence_plot+theme_linedraw()+scale_fill_brewer(palette="Set3")+rremove("legend")
pp+coord_flip()

# Figure 5d
richness<-data.frame(metadata,(apply(final_dataset1[,-c(1:4)]>0, 1, sum)))
colnames(richness)[6]<-"richness"
Significant_richness <- compare_means(richness ~ Phenotype, group.by = "author", data = richness, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.001,0.05,1), symbols = c("***", "*","ns")),ref.group = "Healthy")
richness_plot<-ggboxplot(richness, x = "author_phenotype", y = "richness",
                         color = 'black',fill = "Phenotype",
                         title = "Fig5d",
                         short.panel.labs = FALSE,outlier.shape=NA,xlab = "", ylab = "Species richness",ylim = c(0, 180))
rp<-richness_plot+theme_linedraw()+scale_fill_brewer(palette="Set3")+rremove("legend")
rp+coord_flip()

# End
