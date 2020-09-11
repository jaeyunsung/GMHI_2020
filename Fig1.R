
# This script shows how to reproduce the results in Figure 1 of Gupta et. al (Nature Communications, 2020)
# Author: Vinod K. Gupta, PhD

unavailable_pkg <- setdiff(c("vegan","ape","ade4","ggplot2","ggpubr","easyGgplot2","dplyr"),rownames(installed.packages()))
install.packages(unavailable_pkg, repos = "http://cran.us.r-project.org")

library(vegan)
library(ape)
library(ade4)
library(ggplot2)
library(ggpubr)
library(easyGgplot2)
library(dplyr)
library(ade4)


# I/O
relative_abundance_file <- "./4347_final_relative_abundances.txt"
meta_data_file <- "./Final_metadata_4347.csv"

fig1b_out <- "Fig1b.pdf"
fig1c_out <- "Fig1c.pdf"
fig1d_out <- "Fig1d.pdf"


# Figure 1b
Fig1b_input <- read.table(relative_abundance_file, sep = "\t", header = TRUE,row.names = 1)
Fig1c_d_input <- read.csv(meta_data_file, sep = ",", header = T,row.names = 1,check.names = F)
Fig1b_classified <- Fig1b_input[-grep('unclassified', row.names(Fig1b_input)),] # unclassified species removed
Fig1b_classified[Fig1b_classified < 0.001] <- 0 # threshold of being present for species
Fig1b_classified[Fig1b_classified > 0] <- 1
species_coverage<-apply(Fig1b_classified,1,sum)*100/4347

pdf(file = fig1b_out)
Fig1b <- hist(species_coverage,ylim = c(0,500),
            col = 'orange',xlab = "Species prevalence across 4,347 samples (%)",
            ylab = "Number of species",breaks=100,main = "Fig1b")
dev.off()


# Figures 1c and d
# count for species to be excluded from the analyses due to low coverage
# count for species to be included for the further analysis 
low_coverage_species_count<-Fig1b$counts[1] 
final_species_count<-sum(Fig1b$counts[2:99]) 

Fig1c_d_dataset1 <- data.frame(t(Fig1c_d_input),check.rows = F,check.names = F)
Fig1c_d_dataset <- Fig1c_d_dataset1[,c(1,7,15,15,33:ncol(Fig1c_d_dataset1))]
Fig1c_d_dataset[,-c(1:4)] <- lapply(Fig1c_d_dataset[,-c(1:4)], function(x) as.numeric(as.character(x)))
                                    
colnames(Fig1c_d_dataset)[4] <- "Phenotype_all"
                                    
Fig1c_d_dataset$Phenotype <- gsub(x = Fig1c_d_dataset$Phenotype, 
                                  pattern = "[^Healthy].+|advanced adenoma", replacement = "Nonhealthy")

# Figures 1c and 1d 
# Note: PERMANOVA process is a bit time consuming
# PCoA based on Bray-Curtis distances using asin(sqrt) transformed relative abundances of microbial species

# Multiplied by 0.01 because relative abundances in input file are in percentange (%)
trs <- function(x) asin(sqrt(x*0.01)) 


# Preprocess for Figures 1c and d: START
permanova1 <- Fig1c_d_dataset %>% mutate_each(funs(trs), colnames(Fig1c_d_dataset[,-c(1:4)]))
permanova1$study <- gsub(x = permanova1$study, pattern = "_.+", replacement = "")

BC.dist=vegdist(permanova1[,c(5:317)], method = 'bray')
adonis2(BC.dist ~ Phenotype, data = permanova1, permutations = 999,strata = permanova1$study)

pcoa_all <- dudi.pco(BC.dist, scannf=FALSE, nf=3)
evals <- eigenvals(pcoa_all)
Variance <- evals / sum(evals)
Variance1 <- 100 * signif(Variance[1], 2)
Variance2 <- 100 * signif(Variance[2], 2)
Variance3 <- 100 * signif(Variance[3], 2)

pc_plot_data<-data.frame(pcoa_all$li$A1, pcoa_all$li$A2, permanova1$Phenotype,
                         Fig1c_d_dataset$Phenotype_all, Fig1c_d_dataset$`Sample Accession or Sample ID`)
colnames(pc_plot_data)<-c("x","y","Phenotype","Phenotype_all","Sample_ID")
pc_plot_data1 <- merge(pc_plot_data,aggregate(cbind(mean.x=x,mean.y=y)~Phenotype, 
                                              pc_plot_data,mean), by="Phenotype")
# Preprocess for Figures 1c and d: END


# Draw Figure 1c
pdf(file = fig1c_out)
Fig1c <- ggplot(pc_plot_data1, aes(x,y,color=factor(Phenotype)))+geom_point(size=2)+
        geom_point(aes(x=mean.x,y=mean.y),size=0)+stat_ellipse(level = 0.95)+
        labs(x = paste("PCoA1 (", Variance1, "%)", sep=""), y = paste("PCoA2 (", Variance2, "%)", sep=""),
        title = 'Fig1c')+theme_bw()+theme(legend.title=element_blank())+
        scale_colour_manual(values=c("Healthy"="steelblue","Nonhealthy"="orange2"))
Fig1c+annotate("text", label = paste0(c("R-squared = ",0.01705,'\n',"p-value = ",0.001),collapse = ''), 
               x = -0.35, y = 0.35, color = "black")
dev.off()


# Draw Figure 1d
pdf(file = fig1d_out)
Fig1d <- ggplot(pc_plot_data1, aes(x,y,color=factor(Phenotype_all)))+geom_point(size=2)+
        geom_point(aes(x=mean.x,y=mean.y),size=0)+ 
        labs(x = paste("PCoA1 (", Variance1, "%)", sep=""),y = paste("PCoA2 (", Variance2, "%)", sep=""),title = 'Fig1d')+
        theme_bw()+theme(legend.title=element_blank())+
        scale_colour_manual(values=c("Healthy"="#80B1D3","ACVD"="#FFFFB3","advanced adenoma"="#8DD3C7","CRC"="#FB8072","Crohns disease"="#BEBADA","Obesity"="#B3DE69","Overweight"="#FCCDE5","Rheumatoid Arthritis"="#D9D9D9","Symptomatic atherosclerosis"="#BC80BD","IGT"="#ED1E79","T2D"="#CCEBC5","Ulcerative colitis"="#FFED6F","Underweight"="white"))
Fig1d+annotate("text", label = paste0(c("R-squared = ",0.01705,'\n',"p-value = ",0.001),collapse = ''), 
               x = -0.35, y = 0.35, color = "black")
dev.off()

# End
