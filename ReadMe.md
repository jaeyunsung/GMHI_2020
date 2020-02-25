
# GMHI: Gut Microbiome Health Index

'GMHI.R' is an R script to calculate the Gut Microbiome Health Index (GMHI) for a species-level gut microbiome profile.

# Installation

Download the 'GMHI.R' script from https://github.com/jaeyunsung/GMHI_2020.

# Description

GMHI is a robust index for evaluating health status based on the species-level taxonomic profile of a stool shotgun metagenome (gut microbiome) sample. GMHI is designed to evaluate the balance between two sets of microbial species associated with good and adverse health conditions, which are identified from a meta-analysis on 4,347 publicly available, human stool metagenomes integrated across multiple studies encompassing various phenotypes and geographies. GMHI denotes the degree to which a subject’s stool metagenome sample portrays microbial taxonomic properties associated with healthy (GMHI > 0) or non-healthy (GMHI < 0). A positive or negative GMHI allows the sample to be classified as healthy or non-healthy, respectively; a GMHI of 0 indicates an equal balance of Health-prevalent and Health-scarce species, and thereby classified as neither. Higher (more positive) and lower (more negative) values of GMHI reflects the dominant influence of Health-prevalent species over Health-scarce species in the healthy group, and vice versa in the non-healthy group, respectively.

If you use GMHI.R, please cite or acknowledge our manuscript.

# Pre-requisites

We recommend using R (3.5.x or later) or the latest version of RStudio.

# Basic Usage

Step 1:
Run MetaPhlAn2 on stool metagenome(s) using the '--tax_lev s' argument. It will produce a separate output file (of clade-specific relative abundances) for each metagenome.

Step 2:
If there are multiple stool metagenome samples to be processed, first merge the MetaPhlAn2 outputs using 'merge_metaphlan_tables.py' provided in the MetaPhlAn2 pipeline (follow online MetaPhlAn2 tutorial), and then save the merged file in .csv format. This file should contain species' names in the first column, and corresponding relative abundances in subsequent columns for all metagenome samples (see 'species_relative_abundances.csv').

Step 3:
Open the GMHI.R script and load the merged MetaPhlAn2 output file (generated in Step 2) by providing the appropriate path. Then run GMHI.R in its entirety.
 

# Example

A simple demo of GMHI.R can be conducted using the following files:

Files:
1. 'species_relative_abundances.csv': Merged version of MetaPhlAn2 output files showing only the species' relative abundances from 25 stool metagenome samples.
2. 'GMHI_output.csv': GMHI values (generated from GMHI.R) for all 25 stool metagenomes in 'species_relative_abundances.csv'.
