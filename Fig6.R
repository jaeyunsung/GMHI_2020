
# This script shows how to reproduce the results in Figure 6 of Gupta et. al (Nature Communications, 2020)
# Author: Vinod K. Gupta, PhD

unavailable_pkg <- setdiff(c("vegan","ggplot2","ggpubr","rcompanion"),
                           rownames(installed.packages()))
install.packages(unavailable_pkg, repos = "http://cran.us.r-project.org")

library(vegan)
library(ggplot2)
library(ggpubr)
library(rcompanion)


# Functions---------------------------------------------
parse_file_list <- function(data_file){
    data_list <- scan(data_file, what="", sep="\n")
    return (data_list)
}


# I/O
Final_microbiome_data_4347_file <- "Final_metadata_4347.csv"
validation_abundance_file <- "validation_abundance.csv"
validation_meta_file <- "validation_metadata.csv"
study_wise_file <- "study_wise_data.txt"
MH_species_file <- "./MH_species.txt"
MN_species_file <- "./MN_species.txt"

fig6a_out <- "Fig6a.pdf"
fig6b_out <- "Fig6b.pdf"

Final_microbiome_data_4347 <- read.csv(Final_microbiome_data_4347_file, sep = ",", 
                                       header = T,row.names = 1,check.names = F)

fig6 <- read.csv(validation_abundance_file, sep = ",", header = TRUE,row.names = 1,check.names = F) 
fig6_metadata <- read.csv(validation_meta_file, sep = ",", header = TRUE,row.names = 1,check.names = F) 


# GMHI calculation for validation datasets
H_species <- parse_file_list(MH_species_file) # from Table 2
NH_species <- parse_file_list(MN_species_file) # from Table 2


# Processing of meta-data from validation cohort
fig6_1 <- fig6[-grep('unclassified', row.names(fig6)),] # unclassified species removed
fig6_2 <- fig6_1[-grep('virus', row.names(fig6_1)),] # viruses removed
fig6_3<-sweep(fig6_2,2,colSums(fig6_2),`/`) # re-normalization to get relative abundances after removing unclassified and viruses
fig6_4<-(fig6_3)*100
fig6_4[fig6_4 < 0.001] <- 0


# Healthy normalization constant calculated based on training dataset
# Nonhealthy normalization constant calculated based on training dataset
H_constant <- 7 
NH_constant <- 31 


# Extracting H-plus species
# Extracting H-minus species
validation_H <- fig6_4[row.names(fig6_4) %in% H_species, ] 
validation_NH <- fig6_4[row.names(fig6_4) %in% NH_species, ] 
alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}

H_shannon <- apply((validation_H/100), 2, alpha)
NH_shannon <- apply((validation_NH/100), 2, alpha)

H_sig_count <- apply(validation_H, 2, function(i) (sum(i > 0)))
NH_sig_count <- apply(validation_NH, 2, function(i) (sum(i > 0)))
                      
H_GMHI <- ((H_sig_count/H_constant)*H_shannon)
NH_GMHI <- ((NH_sig_count/NH_constant)*NH_shannon)
GMHI <- data.frame(log10((H_GMHI+0.00001)/(NH_GMHI+0.00001)))
colnames(GMHI) <- c("GMHI")


# Shannon diversity calculation for validation cohort
alpha_shannon_fig6 <- data.frame(diversity(fig6_4, 2, index="shannon"))
fig6_alpha_GMHI_data <- merge(as.data.frame(GMHI),as.data.frame(alpha_shannon_fig6), by='row.names', all=F)
colnames(fig6_alpha_GMHI_data) <- c("Sample_id","GMHI","alpha")
row.names(fig6_alpha_GMHI_data) <- fig6_alpha_GMHI_data$Sample_id
fig6_alpha_GMHI_data1 <- merge(as.data.frame(fig6_alpha_GMHI_data),
                               as.data.frame(fig6_metadata), by='row.names', all=F)
    
                  
# Remove obesity, overweight and underweight samples
fig6_alpha_GMHI_data_final <- fig6_alpha_GMHI_data1[-grep(':OB|:OW|:UW', fig6_alpha_GMHI_data1$Phenotype_all),] 


# healthy_nonhealthy_validation_sample_counts                      
H_Un_count <- data.frame(table(fig6_alpha_GMHI_data_final$Phenotype))


# phenotype_wise_validation_sample_counts
All_phenotype_counts <- data.frame(table(fig6_alpha_GMHI_data_final$Phenotype_all)) 
fig6_alpha_GMHI_data_final$Phenotype <- gsub(x = fig6_alpha_GMHI_data_final$Phenotype, 
                                             pattern = "Unhealthy", replacement = "Nonhealthy")


# Figure 6a (GMHI for combined phenotype as Healthy vs Nonhealthy)
pdf(fig6a_out)
p_gmhi_violin <- ggplot(fig6_alpha_GMHI_data_final, aes(x=Phenotype, y=GMHI, fill=Phenotype)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1,fill="white")

p_gmhi_violin+rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),
                                      axis.title=element_text(size=14,face="bold"))+
                scale_color_manual(values=c('steelblue','orange2'))+
                stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+
                scale_fill_manual(values=c('steelblue','orange2'))+
                labs(x = "",y="Gut Microbiome Health Index (GMHI)",title = "Fig6a")+
                theme_bw()+theme(panel.grid = element_blank())+rremove("legend")
dev.off()
cliffDelta(data =fig6_alpha_GMHI_data_final,GMHI~Phenotype)


# For re-ordering purposes
fig6_alpha_GMHI_data_final$Phenotype_all <- factor(fig6_alpha_GMHI_data_final$Phenotype_all, 
        levels=c('Healthy_P','N-1:Healthy','N-8:Healthy','N-5:AS','N-7:Colorectal adenoma',
                 'N-1:OB','N-8:CRC','N-1:OW','N-7: Italy:CRC','N-7: Japan:CRC',
                 'N-1:UW','CD','Liver cirrhosis','N-4:NAFLD','N-7:OB','N-7:OW','N-7:UW',
                 'N-8:OB','N-8:OW','N-8:UW','Rheumatoid arthritis'))
    

# Figure 6b (GMHI for all study-wise phenotypes)
pdf(fig6b_out)
fig6_GMHI_boxplot <- ggboxplot(fig6_alpha_GMHI_data_final, x = "Phenotype_all", y = "GMHI",
                               color = "Phenotype_all",
                               add = "jitter",
                               short.panel.labs = FALSE,outlier.size=0,xlab = "", ylab = "Gut Microbiome Health Index (GMHI)",title = "Fig6b")
# fig6_GMHI_boxplot +rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))+scale_color_manual(values=c("#FB8072","steelblue","orange2","steelblue","orange2","orange2","#BEBADA","#BEBADA","#FFFFB3","#BEBADA","steelblue","#D9D9D9"))
# Significant_GMHI_H1 <- compare_means(GMHI ~ Phenotype_all, data = fig6_alpha_GMHI_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "Healthy_P")
# Significant_GMHI_H2 <- compare_means(GMHI ~ Phenotype_all, data = fig6_alpha_GMHI_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "N-1:Healthy")
# Significant_GMHI_H3 <- compare_means(GMHI ~ Phenotype_all, data = fig6_alpha_GMHI_data_final, p.adjust.method = "fdr",symnum.args = list(cutpoints = c(0,0.01,1), symbols = c("*","ns")),ref.group = "N-8:Healthy")
fig6_GMHI_boxplot
dev.off()

# End
