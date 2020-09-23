
# This script shows how to calculate the Gut Microbiome Health Index (GMHI) of a MetaPhlAn2 species-level abundance profile, as demonstrated in Gupta et. al (Nature Communications, 2020)
# Author: Vinod K. Gupta, PhD

#Functions---------------------------------------------
parse_species_list <- function(species_file){
	species_list <- scan(species_file, what="", sep="\n")
	return (species_list)
}


#I/O---------------------------------------------------
species_relative_abundance_file <- "./species_relative_abundances.csv"
MH_species_file <- "./MH_species.txt"
MN_species_file <- "./MN_species.txt"
output_file = 'GMHI_output.csv'
#------------------------------------------------------


# Overall Steps
# Step-1: Run MetaPhlAn2 on a stool metagenome using '--tax_lev s' argument.
# Step-2: Merge outputs using 'merge_metaphlan_tables.py' script provided in the MetaPhlAn2 pipeline (see MetaPhlAn2's online tutorial).
# Step-3: Make sure of the following: The merged species' relative abundance profile should be arranged as shown in 'species_relative_abundances.csv'. Accordingly, the first column should contain names of the species-level clades (i.e., taxonomic names with 's__' flag). Subsequent columns should contain the species' relative abundances corresponding to each metagenome sample.
# Step-4: Save input data from Step-3 as a '.csv' file, and run the following script to calculate GMHI for each stool metagenome. GMHI values for each sample in 'species_relative_abundances.csv' are shown in 'GMHI_output.csv'.

# Pre-processing data matrix of species' relative abundances:
# species_profile_1: Result after removing unclassified and virus species
# species_profile_2: Result after transposing species_profile_1
# species_profile_3: Re-normalization of species' relative abundances after removing unclassified and virus species
species_profile <- read.csv(species_relative_abundance_file, sep = ",",
							header = TRUE,row.names = 1,check.names = F) 
library(tidyverse)
tmp1 <- data.frame(t(species_profile),check.rows = F,check.names = F) # transposing data matrix
species_profile_1 <- tmp1 %>% select(-contains(c("unclassified","virus")))
species_profile_2 <- data.frame(t(species_profile_1),check.rows = F,check.names = F)
species_profile_3 <- sweep(species_profile_2,2,colSums(species_profile_2),`/`)
species_profile_3[species_profile_3 < 0.00001] <- 0

MH_species <- parse_species_list(MH_species_file) # Health-prevalent species (7 in total)
MN_species <- parse_species_list(MN_species_file) # Health-scarce species (43 in total)


# Extracting Health-prevalent species present in metagenome
# Extracting Health-scarce species present in metagenome
MH_species_metagenome <- species_profile_3[row.names(species_profile_3) %in% MH_species, ]
MN_species_metagenome <- species_profile_3[row.names(species_profile_3) %in% MN_species, ]


# Diversity among Health-prevalent species
# Diversity among Health-scarce species
alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
MH_shannon <- apply((MH_species_metagenome), 2, alpha) 
MN_shannon <- apply((MN_species_metagenome), 2, alpha) 


# Richness of Health-prevalent species
# Richness of Health-scarce species
R_MH <- apply(MH_species_metagenome, 2, function(i) (sum(i > 0))) 
R_MN <- apply(MN_species_metagenome, 2, function(i) (sum(i > 0)))


# Median RMH from 1% of the top-ranked samples (see Methods)
# Median RMN from 1% of the bottom-ranked samples (see Methods)
MH_prime <- 7
MN_prime <- 31


# Collective abundance of Health-prevalent species
# Collective abundance of Health-scarce species
psi_MH <- ((R_MH/MH_prime)*MH_shannon) 
psi_MN <- ((R_MN/MN_prime)*MN_shannon)

GMHI <- data.frame(log10((psi_MH+0.00001)/(psi_MN+0.00001))) # 0.00001 added to avoid having the denominator as 0
colnames(GMHI) <- c("GMHI")

if (file.exists(output_file)){
	file.remove(output_file)
}
write.csv(GMHI, file=output_file) # Saving GMHI results as 'GMHI_output.csv'. User should change the path to the appropriate directory.

# End
