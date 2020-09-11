
# This script shows how to reproduce the results in Figure 4 of Gupta et. al (Nature Communications, 2020)
# Author: Vinod K. Gupta, PhD

unavailable_pkg <- setdiff(c("ggplot2","ggpubr","vegan","ggplot2","rcompanion"),
                           rownames(installed.packages()))
install.packages(unavailable_pkg, repos = "http://cran.us.r-project.org")

library(ggplot2)
library(ggpubr)
library(vegan)
library(rcompanion)


# I/O
microbiome_data_file <- "Final_metadata_4347.csv"
fig4a_out <- "Fig4a.pdf"
fig4b_out <- "Fig4b.pdf"


# Preprocess: Start
Final_microbiome_data_4347 <- read.csv(microbiome_data_file, sep = ",", 
                                       header = TRUE,row.names = 1,check.names = F) 
names(Final_microbiome_data_4347) <- as.matrix(Final_microbiome_data_4347[15, ])

baseline <- Final_microbiome_data_4347[grep('s__', row.names(Final_microbiome_data_4347)),] # sample with columns Healthy
baseline[] <- lapply(baseline, function(x) type.convert(as.character(x)))
baseline[baseline < 0.001] <- 0


# sample with columns Healthy
# sample with columns Nonhealthy
Healthy <- baseline[,grep('Healthy', names(baseline))]  
Nonhealthy <- baseline[,-grep('Healthy', names(baseline))]  

PH <- apply(Healthy, 1, function(i) (sum(i > 0))*100/2636)
PNH <- apply(Nonhealthy, 1, function(i) (sum(i > 0))*100/1711)
PH_diff <- (PH-PNH)
PH_fold <- (PH/PNH)
PNH_fold <- (PNH/PH)
           
all_matrix <- data.frame(cbind(baseline,PH_diff,PH_fold,PNH_fold))
H_signature <- data.frame(subset(all_matrix, all_matrix$PH_fold >= 1.4 & all_matrix$PH_diff >=10))
NH_signature <- data.frame(subset(all_matrix, all_matrix$PNH_fold >= 1.4 & all_matrix$PH_diff <= -10))
           
alpha_gmhi <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
H_shannon <- apply((H_signature[,-c(4348:4350)]/100), 2, alpha_gmhi)
NH_shannon <- apply((NH_signature[,-c(4348:4350)]/100), 2, alpha_gmhi)
             
H_sig_count <- apply(H_signature[,-c(4348:4350)], 2, function(i) (sum(i > 0)))
NH_sig_count <- apply(NH_signature[,-c(4348:4350)], 2, function(i) (sum(i > 0)))

constant <- data.frame(cbind(H_sig_count,NH_sig_count))
HC1 <- constant[with(constant, order(-H_sig_count, NH_sig_count)), ]
H_constant <- median(HC1$H_sig_count[1:26])
                      
NHC1 <- constant[with(constant, order(H_sig_count, -NH_sig_count)), ]
NH_constant <- median(NHC1$NH_sig_count[1:17])
                      
H_GMHI <- ((H_sig_count/H_constant)*H_shannon)
NH_GMHI <- ((NH_sig_count/NH_constant)*NH_shannon)
                      
GMHI <- data.frame(log10((H_GMHI+0.00001)/(NH_GMHI+0.00001)))
Healthy_GMHI <- data.frame(GMHI[grep('Healthy', row.names(GMHI)),])
Healthy_GMHI$Phenotype<-"Healthy"
                      
Nonhealthy_GMHI <- data.frame(GMHI[-grep('Healthy', row.names(GMHI)),])
Nonhealthy_GMHI$Phenotype<-"Nonhealthy"
                      
colnames(Healthy_GMHI)[1] <- "GMHI"
colnames(Nonhealthy_GMHI)[1] <- "GMHI"
                      
GMHI_20 <- data.frame(rbind(Healthy_GMHI,Nonhealthy_GMHI))
Healthy_accuracy <- sum(Healthy_GMHI$GMHI>0)*100/2636
Nonhealthy_accuracy <- sum(Nonhealthy_GMHI$GMHI<0)*100/1711
                      
total_accuracy <- (Healthy_accuracy+Nonhealthy_accuracy)
report <- cbind(nrow(H_signature),nrow(NH_signature),Healthy_accuracy,Nonhealthy_accuracy,total_accuracy)

Healthy_associated_species <- row.names(H_signature)
Nonhealthy_associated_species <- row.names(NH_signature)
# Preprocess: END


# Figure 4a
GMHI_hist_h <- hist(Healthy_GMHI$GMHI,breaks = 10)
GMHI_hist_h$breaks
GMHI_hist_Un <- hist(Nonhealthy_GMHI$GMHI,breaks = 10)

intervals1 <- GMHI_hist_Un$breaks[1:10]
intervals2 <- GMHI_hist_h$breaks[2:11]
intervals <- paste(intervals1,intervals2,sep = " - ")

Healthy_count <- data.frame(GMHI_hist_h$counts[1:10])
Nonhealthy_count <- data.frame(GMHI_hist_Un$counts[2:11])

line <- data.frame(cbind(intervals,intervals1,intervals2,Healthy_count,Nonhealthy_count))

count_bin <- (Healthy_count) + (Nonhealthy_count)
count_bin_H_percentage <- (Healthy_count)*100/count_bin
count_bin_Un_percentage <- (Nonhealthy_count*100)/count_bin
count_diff_percentage <- count_bin_H_percentage-count_bin_Un_percentage
count_table <- data.frame(intervals,count_bin_H_percentage,count_bin_Un_percentage,count_diff_percentage)
colnames(count_table) <- c("Intervals","Healthy","Nonhealthy","Healthy - Nonhealthy")

pdf(fig4a_out)

par(mar = c(5,5,2,5))
plot(count_table$Healthy, type="o", col="steelblue", pch="o", lty=1,xaxt='n',
     yaxt='n',xlab = "GMHI Bins",ylab = "Proportion of Samples (%)",ylim=c(0,100))
par(new=T)
plot(count_table$Nonhealthy, type="o", col="orange2", pch="o", lty=1,xaxt='n',
     yaxt='n',xlab = '',ylab = '',ylim=c(0,100),main = "Fig4a")
axis(2, seq(100, 0, by = -10), seq(100, 0, by = -10), las = 2)
axis(1,1:20,seq(-4.5,5,by= 0.5))
par(new=T)
hist(GMHI_20$GMHI,breaks = 10,xaxt='n',yaxt='n',main = "",axes=F, xlab=NA, ylab=NA)
axis(side = 4)
mtext(side = 4, line = 3, 'Number of Samples')
dev.off()


# Figure 4b
alpha_shannon <- diversity(baseline, 2, index="shannon")
GMHI_alpha_plot_H_U <- data.frame(cbind(row.names(GMHI),GMHI,alpha_shannon))
colnames(GMHI_alpha_plot_H_U) <- c("Phenotype","GMHI","alpha")

GMHI_alpha_plot_H_U$Phenotype <- gsub(x = GMHI_alpha_plot_H_U$Phenotype, pattern = "[.]", replacement = " ")
GMHI_alpha_plot_H_U$Phenotype <- gsub(x = GMHI_alpha_plot_H_U$Phenotype, pattern = " \\d+", replacement = "")
GMHI_alpha_plot_H_U$Phenotype_all <- GMHI_alpha_plot_H_U$Phenotype
GMHI_alpha_plot_H_U$Phenotype <- gsub(x = GMHI_alpha_plot_H_U$Phenotype, 
                                      pattern = "[^Healthy].+|advanced adenoma", replacement = "Nonhealthy")

pdf(fig4b_out)
fig4b <- ggscatterhist(
  GMHI_alpha_plot_H_U, x = "alpha", y = "GMHI",
  color = "Phenotype",
  palette = c("steelblue", "orange2"),
  margin.plot = "histogram",
  ggtheme = theme_bw(),
  show.legend = FALSE,xlab = "Shannon diversity",ylab = "Gut Microbiome Health Index (GMHI)",title = "Fig4b",
  margin.params = list(fill = "Phenotype", color = "black", size = 0.1)
)
dev.off()
cor.test(GMHI_alpha_plot_H_U$alpha,GMHI_alpha_plot_H_U$GMHI,method = 'spearman')

# End
