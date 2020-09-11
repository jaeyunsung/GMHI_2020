
# This script shows how to reproduce the results in Figure 5 of Gupta et. al (Nature Communications, 2020)
# Author: Vinod K. Gupta, PhD

unavailable_pkg <- setdiff(c("dplyr","vegan","ggpubr","ggplot2","reshape2"),
                           rownames(installed.packages()))
install.packages(unavailable_pkg, repos = "http://cran.us.r-project.org")

library(dplyr)
library(vegan)
library(ggpubr)
library(ggplot2)
library(reshape2)


# Functions---------------------------------------------
parse_file_list <- function(data_file){
    data_list <- scan(data_file, what="", sep="\n")
    return (data_list)
}


# I/O
Final_microbiome_data_4347 <- read.csv("./Final_metadata_4347.csv", sep = ",", 
                                     header = T,row.names = 1,check.names = F)
study_wise_file <- "study_wise_data.txt"
MH_species_file <- "./MH_species.txt"
MN_species_file <- "./MN_species.txt"

fig5a_out <- "Fig5a.pdf"
fig5b_out <- "Fig5b.pdf"
fig5c_out <- "Fig5c.pdf"
fig5d_out <- "Fig5d.pdf"


# Studies having samples from healthy as well as from Nonhealthy individuals
study_wise_data <- parse_file_list(study_wise_file)
fig5_1 <- data.frame(t(Final_microbiome_data_4347),check.rows = F,check.names = F)

final_dataset <- fig5_1[fig5_1$study %in% study_wise_data,]
final_dataset1 <- final_dataset[,c(1:2,4,15,33:ncol(final_dataset))]

colnames(final_dataset1)[2] <- "study1"
metadata <- data.frame(final_dataset1[,c(1:4)])
colnames(metadata)[3] <- "author"

metadata$author_phenotype = paste(metadata$author,"",metadata$Phenotype)
metadata$author_phenotype <- as.factor(metadata$author_phenotype)

final_dataset1[,-c(1:4)] <- lapply(final_dataset1[,-c(1:4)], function(x) as.numeric(as.character(x)))
final_dataset2 <- data.frame(t(final_dataset1),check.rows = F,check.names = F)
final_dataset3 <- data.frame(final_dataset2[-c(1:4),])
final_dataset3[] <- lapply(final_dataset3, function(x) as.numeric(as.character(x)))
final_dataset3[final_dataset3 < 0.001] <- 0


# Figure 5a
H_species <- parse_file_list(MH_species_file)
NH_species <- parse_file_list(MN_species_file) # from table 2

H_constant <- 7 # from Figure 2
NH_constant <- 31 # from Figure 2

sp_H <- final_dataset3[row.names(final_dataset3) %in% H_species, ]
sp_NH <- final_dataset3[row.names(final_dataset3) %in% NH_species, ]

alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
H_shannon <- apply((sp_H/100), 2, alpha)
NH_shannon <- apply((sp_NH/100), 2, alpha)

H_sig_count <- apply(sp_H, 2, function(i) (sum(i > 0)))
NH_sig_count <- apply(sp_NH, 2, function(i) (sum(i > 0)))
                      
H_GMHI <- ((H_sig_count/H_constant)*H_shannon)
NH_GMHI <- ((NH_sig_count/NH_constant)*NH_shannon)
                      
GMHI <- data.frame(log10((H_GMHI+0.00001)/(NH_GMHI+0.00001)))
colnames(GMHI) <- c("GMHI")
Result <- data.frame(cbind(GMHI,H_sig_count,NH_sig_count,H_shannon,NH_shannon,H_GMHI,NH_GMHI))
GMHI_final <- data.frame(metadata,Result)
colnames(GMHI_final)[6] <- "GMHI_final"
                      
Significant_GMHI <- compare_means(GMHI_final ~ Phenotype, group.by = "author", 
                                  data = GMHI_final, p.adjust.method = "fdr",
                                  symnum.args = list(cutpoints = c(0,0.001,0.05,1), 
                                  symbols = c("***", "*","ns")),ref.group = "Healthy")
                      
pdf(fig5a_out)
GMHI_plot <- ggboxplot(GMHI_final, x = "author_phenotype", y = "GMHI_final",
                     color = 'black',fill = "Phenotype",
                     title = "Fig5a",
                     short.panel.labs = FALSE,outlier.shape=NA,
                     xlab = "", ylab = "Gut Microbiome Health Index (GMHI)",ylim = c(-6, 6))
gp <- GMHI_plot+theme_linedraw()+scale_fill_brewer(palette="Set3")+rremove("legend")
gp + coord_flip()
dev.off()


# Figure 5b
pdf(fig5b_out)
alpha_shannon<-data.frame(metadata,(diversity(final_dataset1[,-c(1:4)], 1, index="shannon")))
colnames(alpha_shannon)[6]<-"alpha"
Significant_alpha <- compare_means(alpha ~ Phenotype, group.by = "author", 
                                   data = alpha_shannon, p.adjust.method = "fdr",
                                   symnum.args = list(cutpoints = c(0,0.001,0.05,1), symbols = c("***", "*","ns")),
                                   ref.group = "Healthy")
alpha_plot<-ggboxplot(alpha_shannon, x = "author_phenotype", y = "alpha",
                      color = 'black',title = "Fig5b",fill = "Phenotype",
                      short.panel.labs = FALSE,outlier.shape=NA,xlab = "", 
                      ylab = "Shannon Diversity",ylim = c(0, 5))
ap<-alpha_plot+theme_linedraw()+scale_fill_brewer(palette="Set3")+rremove("legend")
ap+coord_flip()
dev.off()


# Figure 5c
pdf(fig5c_out)
cum<-data.frame(apply(apply(final_dataset3, 2, sort,decreasing=T), 2, cumsum))
cum_cutoff2 <- function(x){max(which(x < 80))+1}
prevalence2<-data.frame(apply(cum, 2, cum_cutoff2))
prevalence2[prevalence2 =="-Inf"] <- 1
colnames(prevalence2)<-"prevalence"
prevalence<-data.frame(metadata,prevalence2$prevalence)
colnames(prevalence)[6]<-"prevalence"

Significant_prevalence <- compare_means(prevalence ~ Phenotype, group.by = "author", 
                                        data = prevalence, p.adjust.method = "fdr",
                                        symnum.args = list(cutpoints = c(0,0.001,0.05,1), symbols = c("***", "*","ns")),
                                        ref.group = "Healthy")

prevalence_plot<-ggboxplot(prevalence, x = "author_phenotype", y = "prevalence",
                           color = 'black',fill = "Phenotype",title = "Fig5c",
                           short.panel.labs = FALSE,outlier.shape=NA,
                           xlab = "", ylab = "80% abundance coverage (# of Species)",ylim = c(0, 40))
pp<-prevalence_plot+theme_linedraw()+scale_fill_brewer(palette="Set3")+rremove("legend")
pp+coord_flip()
dev.off()


# Figure 5d
pdf(fig5d_out)
richness<-data.frame(metadata,(apply(final_dataset1[,-c(1:4)]>0, 1, sum)))
colnames(richness)[6]<-"richness"
Significant_richness <- compare_means(richness ~ Phenotype, group.by = "author", data = richness, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.001,0.05,1), symbols = c("***", "*","ns")),ref.group = "Healthy")
richness_plot<-ggboxplot(richness, x = "author_phenotype", y = "richness",
                         color = 'black',fill = "Phenotype",
                         title = "Fig5d",
                         short.panel.labs = FALSE,outlier.shape=NA,xlab = "", ylab = "Species richness",ylim = c(0, 180))
rp<-richness_plot+theme_linedraw()+scale_fill_brewer(palette="Set3")+rremove("legend")
rp+coord_flip()
dev.off()

# End
