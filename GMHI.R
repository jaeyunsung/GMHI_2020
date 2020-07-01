
# This script shows how to calculate the Gut Microbiome Health Index (GMHI) of a MetaPhlAn2 species-level abundance profile, as demonstrated in Gupta et. al (version June 2020)
# Author: Vinod K. Gupta, PhD

# Step-1: Run MetaPhlAn2 on a stool metagenome using '--tax_lev s' argument.
# Step-2: Merge outputs using 'merge_metaphlan_tables.py' script provided in the MetaPhlAn2 pipeline (see MetaPhlAn2's online tutorial).
# Step-3: Make sure of the following: The merged species' relative abundance profile should be arranged as shown in 'species_relative_abundances.csv'. Accordingly, the first column should contain names of the species-level clades (i.e., taxonomic names with 's__' flag). Subsequent columns should contain the species' relative abundances corresponding to each metagenome sample.
# Step-4: Save input data from Step-3 as a '.csv' file, and run the following script to calculate GMHI for each stool metagenome. GMHI values for each sample in 'species_relative_abundances.csv' are shown in 'GMHI_output.csv'.

species_profile <- read.csv("./species_relative_abundances.csv", sep = ",", header = TRUE,row.names = 1,check.names = F) # User should change the path to the appropriate directory
species_profile_1 <- species_profile[-grep('unclassified', row.names(species_profile)),]  # Removing unclassified species
species_profile_2 <- species_profile_1[-grep('virus', row.names(species_profile_1)),]  # Removing virus species
species_profile_3 <- sweep(species_profile_2,2,colSums(species_profile_2),`/`) # Re-normalizing to get species relative abundances after removing unclassified and virus species
species_profile_3[species_profile_3 < 0.00001] <- 0 # Discarding species that are of very low relative abundance from downstream analyses

MH_species <- c("s__Alistipes_senegalensis", "s__Bacteroidales_bacterium_ph8",
             "s__Bifidobacterium_adolescentis", "s__Bifidobacterium_angulatum",
             "s__Bifidobacterium_catenulatum", "s__Lachnospiraceae_bacterium_8_1_57FAA",
             "s__Sutterella_wadsworthensis") # Health-prevalent species

MN_species <- c("s__Anaerotruncus_colihominis", "s__Atopobium_parvulum", "s__Bifidobacterium_dentium",
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
              "s__Veillonella_atypica") # Health-scarce species

MH_species_metagenome <- species_profile_3[row.names(species_profile_3) %in% MH_species, ] # Extracting Health-prevalent species present in metagenome
MN_species_metagenome <- species_profile_3[row.names(species_profile_3) %in% MN_species, ] # Extracting Health-scarce species present in metagenome
alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
MH_shannon <- apply((MH_species_metagenome), 2, alpha) # Diversity among Health-prevalent species
MN_shannon <- apply((MN_species_metagenome), 2, alpha) # Diversity among Health-scarce species
R_MH <- apply(MH_species_metagenome, 2, function(i) (sum(i > 0))) # Richness of Health-prevalent species
R_MN <- apply(MN_species_metagenome, 2, function(i) (sum(i > 0))) # Richness of Health-scarce species
MH_prime <- 7 # Median RMH from 1% of the top-ranked samples (see Supplementary Methods for further details)
MN_prime <- 31 # Median RMN from 1% of the bottom-ranked samples (see Supplementary Methods for further details)
psi_MH <- ((R_MH/MH_prime)*MH_shannon) # Collective abundance of Health-prevalent species
psi_MN <- ((R_MN/MN_prime)*MN_shannon) # Collective abundance of Health-scarce species
GMHI <- data.frame(log10((psi_MH+0.00001)/(psi_MN+0.00001))) # 0.00001 added to avoid having the denominator as 0
colnames(GMHI) <- c("GMHI")
write.csv(GMHI,file = "/GMHI_output.csv") # Saving GMHI results as 'GMHI_output.csv'. User should change the path to the appropriate directory.
