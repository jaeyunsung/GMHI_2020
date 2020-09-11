
# This script shows how to reproduce the results in Figure 2, Table 2, and Supplementary Table 1 of Gupta et. al (Nature Communications, 2020)
# Author: Vinod K. Gupta, PhD

unavailable_pkg <- setdiff(c("ggplot2","ggpubr","vegan","rcompanion"),rownames(installed.packages()))
install.packages(unavailable_pkg, repos = "http://cran.us.r-project.org")

library('ggplot2')
library('ggpubr')
library('vegan')
library('rcompanion')


#I/O
microbiome_data_file <- "Final_metadata_4347.csv"
fig2a_out <- "Fig2a.pdf"
fig2b_out <- "Fig2b.pdf"

Final_microbiome_data_4347 <- read.csv(microbiome_data_file, sep = ",", header = TRUE,
                                       row.names = 1,check.names = F)
names(Final_microbiome_data_4347) <- as.matrix(Final_microbiome_data_4347[15, ])
baseline <- Final_microbiome_data_4347[grep('s__', row.names(Final_microbiome_data_4347)),] # species data
baseline[] <- lapply(baseline, function(x) type.convert(as.character(x)))


# Supplementary Table 1. 
# Identification of best accuracy thresholds for fold change in prevalence and difference in prevalence between Healthy and Nonhealthy
baseline[baseline < 0.001] <- 0
Healthy <- baseline[,grep('Healthy', names(baseline))]
Nonhealthy <- baseline[,-grep('Healthy', names(baseline))]


# PH: species prevalence among healthy
# PNH: species prevalence among non-healthy
PH <- apply(Healthy, 1, function(i) (sum(i > 0))*100/2636) 
PNH <- apply(Nonhealthy, 1, function(i) (sum(i > 0))*100/1711) 
           
PH_diff <- (PH-PNH)
PH_fold <- (PH/PNH)
PNH_fold <- (PNH/PH)
             
all_matrix<-data.frame(cbind(baseline,PH_diff,PH_fold,PNH_fold))
             
report <- list()
GMHI_20_list<- list()
             
cutoffs1 <- seq(1.2, 2, by=0.1)
cutoffs2 <- c(5,10,15,20)

alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}

for(i in 1:length(cutoffs1)) {
    
  H_signature_sublist <- list()
  NH_signature_sublist <- list()
  report_sublist <- list()
  GMHI_list <- list()
    
  for(j in 1:length(cutoffs2)) {
      
    H_signature_sublist[[j]] <- data.frame(subset(as.matrix(all_matrix), 
                                                  all_matrix$PH_fold >= cutoffs1[i] & all_matrix$PH_diff >= cutoffs2[j]))
    NH_signature_sublist[[j]] <- data.frame(subset(as.matrix(all_matrix), 
                                                   all_matrix$PNH_fold >= cutoffs1[i] & all_matrix$PH_diff <= -cutoffs2[j]))
       
    H_shannon <- apply((H_signature_sublist[[j]][,-c(4348:4350)]/100), 2, alpha)
    NH_shannon <- apply((NH_signature_sublist[[j]][,-c(4348:4350)]/100), 2, alpha)
      
    H_sig_count <- apply(H_signature_sublist[[j]][,-c(4348:4350)], 2, function(i) (sum(i > 0)))
    NH_sig_count <- apply(NH_signature_sublist[[j]][,-c(4348:4350)], 2, function(i) (sum(i > 0)))
                          
    constant <- data.frame(cbind(H_sig_count,NH_sig_count))
                          
    HC1 <- constant[with(constant, order(-H_sig_count, NH_sig_count)), ]
    H_constant <- median(HC1$H_sig_count[1:26])
                          
    NHC1 <- constant[with(constant, order(H_sig_count, -NH_sig_count)), ]
    NH_constant <- median(NHC1$NH_sig_count[1:17])
                          
    H_GMHI <- ((H_sig_count/H_constant)*H_shannon)
    NH_GMHI <- ((NH_sig_count/NH_constant)*NH_shannon)
                          
    GMHI_list[[j]] <- data.frame(log10((H_GMHI+0.00001)/(NH_GMHI+0.00001)))
    GMHI_20_list[[i]] <- GMHI_list
    GMHI<-data.frame(log10((H_GMHI+0.00001)/(NH_GMHI+0.00001)))
                          
    Healthy_GMHI <- data.frame(GMHI[grep('Healthy', row.names(GMHI)),])
    Healthy_GMHI$Phenotype <- "Healthy"
                          
    Nonhealthy_GMHI <- data.frame(GMHI[-grep('Healthy', row.names(GMHI)),])
    Nonhealthy_GMHI$Phenotype <- "Nonhealthy"
                          
    colnames(Healthy_GMHI)[1] <- "GMHI_20"
    colnames(Nonhealthy_GMHI)[1] <- "GMHI_20"
                          
    GMHI_20<-data.frame(rbind(Healthy_GMHI,Nonhealthy_GMHI))
                          
    Healthy_accuracy<-sum(Healthy_GMHI$GMHI_20>0)*100/2636
    Nonhealthy_accuracy<-sum(Nonhealthy_GMHI$GMHI_20<0)*100/1711
                          
    total_accuracy<-(Healthy_accuracy+Nonhealthy_accuracy)
    total_average_accuracy<-(Healthy_accuracy+Nonhealthy_accuracy)/2
                          
    report_sublist[[j]]<-cbind(cutoffs1[i],cutoffs2[j],nrow(H_signature_sublist[[j]]),
                               nrow(NH_signature_sublist[[j]]),Healthy_accuracy,Nonhealthy_accuracy,
                               total_accuracy,total_average_accuracy)
    report[[i]] <- report_sublist
  }
}

Accuracy_table <- matrix(unlist(report), ncol = 8, byrow = TRUE)
colnames(Accuracy_table) <- c("Fold change", "Difference in coverage b/w healthy and non-healthy",
                              "H+ count","H- count","Healthy accuracy","Non-healthy accuracy",
                              "Total accuracy","Total average accuracy")

Supplementary_Table3<-data.frame(Accuracy_table)
Supplementary_Table3<-Supplementary_Table3[,c(1:4,8)]

colnames(Supplementary_Table3)<-c("Fold change", "Difference in prevalence b/w healthy and Nonhealthy",
                                  "Health-prevalent species count","Health-scarce species count",
                                  "Classification accuracy")


# GMHI calculation for training dataset based on best classification accuracy using 1.4 fold change of prevalence and 10% of difference in prevalences between healthy and non-healthy
final_dataset <- data.frame(t(Final_microbiome_data_4347),check.rows = F,check.names = F)
final_dataset1 <- final_dataset[,c(1,2,11:13,15:20,33:ncol(final_dataset))]

colnames(final_dataset1)[c(2,3,4,6:11)] <- c("study1","age","BMI","Phenotype","FBG",
                                             "TRIG","LDLC","CHOL","HDLC")

final_dataset1[,-c(1:11)] <- lapply(final_dataset1[,-c(1:11)], function(x) as.numeric(as.character(x)))
final_dataset2<-data.frame(t(final_dataset1),check.rows = F,check.names = F)
final_dataset3<-data.frame(final_dataset2[-c(1:11),])
                                    
final_dataset3[] <- lapply(final_dataset3, function(x) as.numeric(as.character(x)))
final_dataset3[final_dataset3 < 0.001] <- 0


# Identifying signature species
H_signature <- data.frame(subset(all_matrix, all_matrix$PH_fold >= 1.4 & all_matrix$PH_diff >=10))
NH_signature <- data.frame(subset(all_matrix, all_matrix$PNH_fold >= 1.4 & all_matrix$PH_diff <= -10))

H_species <- row.names(H_signature)
NH_species <- row.names(NH_signature)

H_constant <- 7
NH_constant <- 31

sp_H <- final_dataset3[row.names(final_dataset3) %in% H_species, ]
sp_NH <- final_dataset3[row.names(final_dataset3) %in% NH_species, ]

alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}

H_shannon <- apply((sp_H/100), 2, alpha)
NH_shannon <- apply((sp_NH/100), 2, alpha)

H_sig_count <- apply(sp_H, 2, function(i) (sum(i > 0)))
NH_sig_count <- apply(sp_NH, 2, function(i) (sum(i > 0)))
                    
H_GMHI <- ((H_sig_count/H_constant)*H_shannon)
NH_GMHI <- ((NH_sig_count/NH_constant)*NH_shannon)
                    
GMHI<-data.frame(log10((H_GMHI+0.00001)/(NH_GMHI+0.00001)))
colnames(GMHI) <- c("GMHI")
                    
Result<-data.frame(cbind(GMHI,H_sig_count,NH_sig_count,H_shannon,NH_shannon,H_GMHI,NH_GMHI))
metadataset<-data.frame(final_dataset1[,c(1:11)])
GMHI_final<-merge(as.data.frame(metadataset),as.data.frame(Result), by='row.names', all=F)
row.names(GMHI_final)<-GMHI_final$Row.names


# Table 2. Microbial species of the Health-prevalent and Health-scarce groups
table2_input <- data.frame(cbind(PH,PNH,PH_diff,
                                 PH_fold,PNH_fold))
H_plus_species <- table2_input[row.names(table2_input) %in% H_species,]
H_plus_species$PNH_fold <- NULL

colnames(H_plus_species)[4] <- "Fold-change (PH/PNH or PNH/PH)"

H_minus_species <- table2_input[row.names(table2_input) %in% NH_species,]
H_minus_species$PH_fold <- NULL

colnames(H_minus_species)[4] <- "Fold-change (PH/PNH or PNH/PH)"

table2<-data.frame(rbind(H_plus_species,H_minus_species))
colnames(table2)<-c("Prevalence in Healthy samples (%)","Prevalence in Nohealthy samples (%)",
                    "Difference (PH-PNH)","Fold-change (PH/PNH or PNH/PH)")


# Figure 2a
pdf(fig2a_out)

HDLC_plot_data<-GMHI_final[grep("[[:digit:]]", GMHI_final$HDLC), ]
HDLC_plot_data$HDLC <- as.numeric(as.character(HDLC_plot_data$HDLC))
HDLC_plot_data<-data.frame(subset(HDLC_plot_data, HDLC_plot_data$HDLC > 0))
HDLC_plot_data<-data.frame(subset(HDLC_plot_data, HDLC_plot_data$HDLC < 5)) # remove outliers
HDLC_plot_data$Phenotype<-gsub(x = HDLC_plot_data$Phenotype, pattern = "[^Healthy].+|advanced adenoma", replacement = "Non-healthy")
table(HDLC_plot_data$Phenotype)
HDLC_plot_data$HDLC<-(HDLC_plot_data$HDLC)*38.67 # unit conversion from mmol/L to mg/dL
HDLC_plot<-ggscatter(HDLC_plot_data, x = "HDLC", y = "GMHI",
                     color = "black", shape = 21, size = 3,
                     add = "reg.line",
                     add.params = list(color = "blue", fill = "lightgray"),
                     conf.int = TRUE,ylim = c(-5, 5),
                     cor.coef = TRUE, ylab = "Gut Microbiome Health Index (GMHI)",xlab = "High-Density Lipoproteins Cholesterol (mg/dL)",
                     cor.method = "spearman",cor.coeff.args = list(method = "spearman", label.sep = "\n"))+geom_point(aes(color=Phenotype))
HDLC_plot+scale_colour_manual(values=c("Healthy"="steelblue","Non-healthy"="orange2"))
dev.off()


# Figure 2b
pdf(fig2b_out)

HDLC_plot_data<-data.frame(subset(HDLC_plot_data, HDLC_plot_data$GMHI > 0|HDLC_plot_data$GMHI < 0))
HDLC_plot_data$GMHI_group = cut(HDLC_plot_data$GMHI,c(0,-6,6))
levels(HDLC_plot_data$GMHI_group) = c("GMHI_neg","GMHI_pos")
level_order <- c('GMHI_pos', 'GMHI_neg')
HDLC_plot<-ggplot(HDLC_plot_data, aes(x = factor(GMHI_group,level = level_order), y=HDLC, fill=GMHI_group)) +
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1, fill="white")+theme_classic()
HDLC_plot +rremove("legend")+theme(axis.text=element_text(size=14,face="bold"),
     axis.title=element_text(size=14,face="bold"))+
	 scale_colour_manual(values=c("GMHI_pos"="steelblue","GMHI_neg"="orange2"))+
	 stat_compare_means(label = "p.format",method = "wilcox.test",label.x.npc = "middle")+
	 scale_fill_manual(values=c("GMHI_pos"="steelblue","GMHI_neg"="orange2"))+
	 labs(x = "",y="High-Density Lipoproteins Cholesterol (mg/dL)")
table(HDLC_plot_data$GMHI_group)
cliffDelta(HDLC~GMHI_group,data = HDLC_plot_data) # effect size
dev.off()

#End
