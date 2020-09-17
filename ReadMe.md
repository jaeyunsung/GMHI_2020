A Predictive Index for Health Status using Species-level Gut Microbiome Profiling
=================================================================================
DOI: 10.1038/s41467-020-18476-8

Vinod K. Gupta, Minsuk Kim, Utpal Baksh, Kevin Y. Cunningham, John M. Davis III, Konstantinos N. Lazaridis, Heidi Nelson, Nicholas Chia, & Jaeyun Sung


# Description (GMHI)

GMHI is a robust index for evaluating health status based on the species-level taxonomic profile of a stool shotgun metagenome (gut microbiome) sample. GMHI is designed to evaluate the balance between two sets of microbial species associated with good and adverse health conditions, which are identified from a meta-analysis on 4,347 publicly available, human stool metagenomes integrated across multiple studies encompassing various phenotypes and geographies. GMHI denotes the degree to which a subject’s stool metagenome sample portrays microbial taxonomic properties associated with healthy (GMHI > 0) or non-healthy (GMHI < 0). A positive or negative GMHI allows the sample to be classified as healthy or non-healthy, respectively; a GMHI of 0 indicates an equal balance of Health-prevalent and Health-scarce species, and thereby classified as neither. Higher (more positive) and lower (more negative) values of GMHI reflects the dominant influence of Health-prevalent species over Health-scarce species in the healthy group, and vice versa in the non-healthy group, respectively.

# Basic Usage of GMHI

Step 1: Run MetaPhlAn2 on stool metagenome(s) using the '--tax_lev s' argument. It will produce a separate output file (of clade-specific relative abundances) for each metagenome.

Step 2: If there are multiple stool metagenome samples to be processed, first merge the MetaPhlAn2 outputs using 'merge_metaphlan_tables.py' provided in the MetaPhlAn2 pipeline (follow online MetaPhlAn2 tutorial), and then save the merged file in .csv format. This file should contain species' names in the first column, and corresponding relative abundances in subsequent columns for all metagenome samples (see 'species_relative_abundances.csv').

Step 3: Open the GMHI.R script and load the merged MetaPhlAn2 output file (generated in Step 2) by providing the appropriate path. Then run GMHI.R in its entirety.

# Description (Figures)

R scripts for reproducing the figures illustrated in the corresponding paper. Tested on R (3.6.1)

# Scripts

#### GMHI (Gut Microbiome Health Index)

R script to calculate the Gut Microbiome Health Index (GMHI) for a species-level gut microbiome profile.

>R script GMHI.R

#### Figure 1

Includes Figure 1 (b,c,d) generation with PERMANOVA (iterations = 999). 

>R script Fig1.R

#### Supplementary Table 1, Table 2, Figure 2

Related with Supplementary Table 1, Table 2, and Figure 2 (a,b) generation. 

>R script TableS1_Table2_Fig2.R

#### Figure 3

Related with Figure 3 (a ~ h) generation.

>R script Fig3.R

#### Figure 4

Related with Figure 4 (a ~ b) generation.

>R script Fig4.R

#### Figure 5

Related with Figure 5 (a ~ d) generation.

>R script Fig5.R

#### Figure 6

Related with Figure 6 (a ~ b).

>R script Fig6.R

#### Data

Merged version of MetaPhlAn2 output files showing only the species' relative abundances from 25 stool metagenome samples.
````
species_relative_abundances.csv
````
GMHI values (generated from GMHI.R) for all 25 stool metagenomes in 'species_relative_abundances.csv'
```
GMHI_output.csv
```
Other Meta data

```
validation_metadata.csv
study_wise_data.txt
MH_species.txt
MN_species.txt
```
