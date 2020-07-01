
# This script shows how to reproduce the results in Figure 3 of Gupta et. al (version June 2020)
# Author: Vinod K. Gupta, PhD

install.packages('vegan')
install.packages('ape')
install.packages('ade4')
install.packages('ggplot2')
install.packages('ggpubr')
install.packages('easyGgplot2')
install.packages('dplyr')
install.packages('ade4')
install.packages('rcompanion')

library(vegan)
library(ape)
library(ade4)
library(ggplot2)
library(ggpubr)
library(easyGgplot2)
library(dplyr)
library(ade4)
library(rcompanion)

Final_microbiome_data_4347<-read.csv("./Final_metadata_4347.csv", sep = ",", header = TRUE,row.names = 1,check.names = F) # change path to the directory wherein this file is located
names(Final_microbiome_data_4347) <- as.matrix(Final_microbiome_data_4347[15, ])
baseline <- Final_microbiome_data_4347[grep('s__', row.names(Final_microbiome_data_4347)),] # species data
baseline[] <- lapply(baseline, function(x) type.convert(as.character(x)))
baseline[baseline < 0.001] <- 0
Healthy <- baseline[,grep('Healthy', names(baseline))]  # sample with columns Healthy
Nonhealthy <- baseline[,-grep('Healthy', names(baseline))]  # sample with columns Nonhealthy
PH<-apply(Healthy, 1, function(i) (sum(i > 0))*100/2636)
PNH<-apply(Nonhealthy, 1, function(i) (sum(i > 0))*100/1711)
PH_diff<-(PH-PNH)
PH_fold<-(PH/PNH)
PNH_fold<-(PNH/PH)
all_matrix<-data.frame(cbind(baseline,PH_diff,PH_fold,PNH_fold))
H_signature<-data.frame(subset(all_matrix, all_matrix$PH_fold >= 1.4 & all_matrix$PH_diff >=10))
NH_signature<-data.frame(subset(all_matrix, all_matrix$PNH_fold >= 1.4 & all_matrix$PH_diff <= -10))
alpha_gmhi <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
H_shannon<-apply((H_signature[,-c(4348:4350)]/100), 2, alpha_gmhi)
NH_shannon<-apply((NH_signature[,-c(4348:4350)]/100), 2, alpha_gmhi)
H_sig_count<-apply(H_signature[,-c(4348:4350)], 2, function(i) (sum(i > 0)))
NH_sig_count<-apply(NH_signature[,-c(4348:4350)], 2, function(i) (sum(i > 0)))
constant<-data.frame(cbind(H_sig_count,NH_sig_count))
HC1<-constant[with(constant, order(-H_sig_count, NH_sig_count)), ]
H_constant<-median(HC1$H_sig_count[1:26])
NHC1<-constant[with(constant, order(H_sig_count, -NH_sig_count)), ]
NH_constant<-median(NHC1$NH_sig_count[1:17])
H_GMHI<-((H_sig_count/H_constant)*H_shannon)
NH_GMHI<-((NH_sig_count/NH_constant)*NH_shannon)
GMHI<-data.frame(log10((H_GMHI+0.00001)/(NH_GMHI+0.00001)))
Healthy_GMHI <- data.frame(GMHI[grep('Healthy', row.names(GMHI)),])
Healthy_GMHI$Phenotype<-"Healthy"
Nonhealthy_GMHI <- data.frame(GMHI[-grep('Healthy', row.names(GMHI)),])
Nonhealthy_GMHI$Phenotype<-"Nonhealthy"
colnames(Healthy_GMHI)[1] <- "GMHI"
colnames(Nonhealthy_GMHI)[1] <- "GMHI"
GMHI_20<-data.frame(rbind(Healthy_GMHI,Nonhealthy_GMHI))
Healthy_accuracy<-sum(Healthy_GMHI$GMHI>0)*100/2636
Nonhealthy_accuracy<-sum(Nonhealthy_GMHI$GMHI<0)*100/1711
total_accuracy<-(Healthy_accuracy+Nonhealthy_accuracy)
report<-cbind(nrow(H_signature),nrow(NH_signature),Healthy_accuracy,Nonhealthy_accuracy,total_accuracy)
report

# Figures 3a and 3e (GMHI)
fig3a<-ggplot(GMHI_20, aes(x=Phenotype, y=GMHI, fill=Phenotype)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")+theme_classic()
fig3a +rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_colour_manual(values=c("Healthy"="steelblue","Nonhealthy"="orange2"))+stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+scale_fill_manual(values=c('steelblue','orange2'))+labs(x = "",y="Gut Microbiome Health Index (GMHI)",title = "Fig3a")
cliffDelta(data = GMHI_20,GMHI~Phenotype) # effect size
GMHI_phenotype_data<-GMHI
GMHI_phenotype_data$Phenotype<-row.names(GMHI_phenotype_data)
colnames(GMHI_phenotype_data)<-c("GMHI","Phenotype_all")
GMHI_phenotype_data$Phenotype_all<-gsub(x = GMHI_phenotype_data$Phenotype_all, pattern = "[.]", replacement = " ")
GMHI_phenotype_data$Phenotype_all<-gsub(x = GMHI_phenotype_data$Phenotype_all, pattern = " \\d+", replacement = "")
fig3e<-ggplot(GMHI_phenotype_data, aes(x=Phenotype_all, y=GMHI, fill=Phenotype_all)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig3e+scale_fill_brewer(palette="Set3")+theme(axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"))+labs(x = "",y="Gut Microbiome Health Index (GMHI)",title = "Fig3e")+theme_bw()+theme(panel.grid = element_blank())+rremove("legend")
Significance_GMHI <- compare_means(GMHI ~ Phenotype_all, data = GMHI_phenotype_data, p.adjust.method = "fdr",ref.group = "Healthy",symnum.args = list(cutpoints = c(0,0.001,1), symbols = c("***", "ns")))

# Figures 3b and 3f (shannon diversity)
fig3_dataset1<-data.frame(t(Final_microbiome_data_4347),check.rows = F,check.names = F)
fig3_dataset<-fig3_dataset1[,c(1,7,15,15,33:ncol(fig3_dataset1))]
fig3_dataset[,-c(1:4)] <- lapply(fig3_dataset[,-c(1:4)], function(x) as.numeric(as.character(x)))
colnames(fig3_dataset)[4]<-"Phenotype_all"
fig3_dataset$Phenotype <- gsub(x = fig3_dataset$Phenotype, pattern = "[^Healthy].+|advanced adenoma", replacement = "Nonhealthy")
alpha_shannon<-data.frame(fig3_dataset[,c(1:4)],(diversity(fig3_dataset[,-c(1:4)], 1, index="shannon")))
colnames(alpha_shannon)[c(2,5)]<-c("Sample Accession","alpha")
fig3b<-ggplot(alpha_shannon, aes(x=Phenotype, y=alpha, fill=Phenotype)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")+theme_classic()
fig3b +rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_colour_manual(values=c("Healthy"="steelblue","Nonhealthy"="orange2"))+stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+scale_fill_manual(values=c('steelblue','orange2'))+labs(x = "",y="Shannon diversity",title = "fig3b")
cliffDelta(data = alpha_shannon,alpha~Phenotype)
fig3f<-ggplot(alpha_shannon, aes(x=Phenotype_all, y=alpha, fill=Phenotype_all)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig3f+scale_fill_brewer(palette="Set3")+theme(axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"))+labs(x = "",y="Shannon diversity",title = "fig3f")+rremove("legend")
Significance_alpha_diversity <- compare_means(alpha ~ Phenotype_all, data = alpha_shannon, p.adjust.method = "fdr",ref.group = "Healthy",symnum.args = list(cutpoints = c(0,0.001,1), symbols = c("***", "ns")))

# Figures 3c and 3g (80% abundance coverage for species)
dominant_species1<-data.frame(t(fig3_dataset),check.rows = F,check.names = F)
dominant_species2<-data.frame(dominant_species1[-c(1:4),])
dominant_species2[] <- lapply(dominant_species2, function(x) as.numeric(as.character(x)))
cum<-data.frame(apply(apply(dominant_species2, 2, sort,decreasing=T), 2, cumsum))
cum_cutoff2 <- function(x){max(which(x < 80))+1}
prevalence2<-data.frame(apply(cum, 2, cum_cutoff2))
prevalence2[prevalence2 =="-Inf"] <- 1
colnames(prevalence2)<-"80_prev"
prevalence<-data.frame(fig3_dataset[,c(1:4)],prevalence2$`80_prev`)
colnames(prevalence)[5]<-"prevalence"
fig3c<-ggplot(prevalence, aes(x=Phenotype, y=prevalence, fill=Phenotype)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")+theme_classic()
fig3c +rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_colour_manual(values=c("Healthy"="steelblue","Nonhealthy"="orange2"))+stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+scale_fill_manual(values=c('steelblue','orange2'))+labs(x = "",y="80% abundance coverage",title = "fig3c")
cliffDelta(data = prevalence,prevalence~Phenotype)
fig3h<-ggplot(prevalence, aes(x=Phenotype_all, y=prevalence, fill=Phenotype_all)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig3h+scale_fill_brewer(palette="Set3")+theme(axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"))+labs(x = "",y="80% abundance coverage",title = "fig3g")+rremove("legend")
Significance_prevalence <- compare_means(prevalence ~ Phenotype_all, data = prevalence, p.adjust.method = "fdr",ref.group = "Healthy",symnum.args = list(cutpoints = c(0,0.001,1), symbols = c("***", "ns")))

# Figures 3d and 3h (species richness)
richness<-data.frame(fig3_dataset[,c(1:4)],(apply(fig3_dataset[,-c(1:4)]>0, 1, sum)))
colnames(richness)[5]<-"richness"
fig3d<-ggplot(richness, aes(x=Phenotype, y=richness, fill=Phenotype)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")+theme_classic()
fig3d +rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_colour_manual(values=c("Healthy"="steelblue","Nonhealthy"="orange2"))+stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+scale_fill_manual(values=c('steelblue','orange2'))+labs(x = "",y="Species richness",title = "fig3d")
cliffDelta(data = richness,richness~Phenotype)
fig3i<-ggplot(richness, aes(x=Phenotype_all, y=richness, fill=Phenotype_all)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")
fig3i+scale_fill_brewer(palette="Set3")+theme(axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"))+labs(x = "",y="Species richness",title = "fig3h")+rremove("legend")
Significance_richness <- compare_means(richness ~ Phenotype_all, data = richness, p.adjust.method = "fdr",ref.group = "Healthy",symnum.args = list(cutpoints = c(0,0.001,1), symbols = c("***", "ns")))

#End
