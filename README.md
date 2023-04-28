# Williamson et al. 2023, Hummingbird blood traits track oxygen availability across space and time, Ecology Letters

In this study, we sought to test whether phylogenetic and population patterns for blood-trait-environment variation match. Specifically, we asked: *Does species adaptation change trait-environment relationships?* 

Data are archived on Dryad: https://doi.org/10.5061/dryad.mkkwh714z.

# DATASET 
We spent 2006–2020 collecting specimen-vouchered hematological data from 1,217 wild hummingbirds of 77 species, representing all 9 clades, and spanning ~4,600 meters in elevation in the Andes (Chile and Peru). We collected data on six blood traits known to affect O2 carrying functions (hemoglobin concentration ([Hb]), hematocrit (Hct), total red blood cell count (TRBC), mean cell volume (MCV), mean cell hemoglobin (MCH), and mean cell hemoglobin concentration (MCHC)). Using these data, we then constructed hierarchical Bayesian models to estimate the responses of six blood traits to elevation, while accounting for other sources of variation. 

Analysis code is available on GitHub: https://github.com/jlwilliamson/ComparativeHummingbirdBlood. All data are linked to vouchered specimens housed at the **Museum of Southwestern Biology** (MSB) at the University of New Mexico, the **Centro de Ornitología y Biodiversidad** in Peru, and the **Pontificia Universidad Católica de Chile** in Chile. Specimen records are accessible in the Artcos database (https://www.arctosdb.org). 

 
# R SCRIPTS

`HumBlood_DataWrangling.Rmd`: This script includes code for processing raw data, evaluating and eliminating outliers, creating new variables, combining data w/ Stotz, making individual elevational range adjustments, processing spatial data and generating sampling map, BioClim processing, reading out of final data file for modeling, plotting, and downstream analyses. Start here. 

`HumBlood_Phylogeny.Rmd`: Script to process McGuire et al. 2014 hummingbird phylogenies for our specific subsets. 

`HumBlood_PhyloElevRangeMapFig.Rmd`: Code and notes for making phylo elev range component of Figure 1 (note that final figure was constructed in Adobe Illustrator; script notes when this happens).

`HumBlood_Modeling.Rmd`: Individual-level (within species analysis) modeling script. Read in processed data from `HumBlood_DataWrangling.Rmd`, standardize predictors, make data subsets that match phylogeny subsets (so there are no NAs in any dataset for modeling), and then run sets of Bayesian hierarchical models w/ brms(). Script describes all individual-level (within species) models. Script also contains the individual contmap components for each blood trait that were used to create **Figure S1** (i.e., trait-specific contmaps were created in R and combined in Illustrator for easy adjustment of graphics).

`HumBlood_Modeling_SpeciesMeans.Rmd`: Species mean (among species analysis) modeling script. Very similar to `HumBlood_Modeling.Rmd` in terms of scope and flow, but this script outlines the procedure for processing data and running and evaluating species mean (among species) models, with necessary modifications. 

`HumBlood_CeNS.Rmd`: Cell Number Versus Size analysis script. Script contains all code and annotations required to reproduce our three-pronged cell number versus cell size analyses: 1) variation in cell number versus size within species, among species, and within single-species, respectively; 2) creation of the Cell Number-Size Index (*CeNS*) using species-specific model coefficients; 3) analysis of the factors that contribute to variation in *CeNS*. This script contains code for **Figure S1** and **Figure S6**. 

`HumBlood_Plots&Figs.Rmd`: This script contains most polished plots and figures for our paper. In most cases, I confined figures/graphics to this standalone script to not bog down analysis scripts with beefy, multi-panel plots, and for ease of future revisions. Script contains code for **Figures 2-4**, **Figure S2**, **Figure S3**, **Figure S4**, **Figure S5** (and entire breakpoint analysis), **Figure S7**, **Figure S8**, and **Figure S9**.

Questions? Contact me at jlw432 [at] cornell.edu. 