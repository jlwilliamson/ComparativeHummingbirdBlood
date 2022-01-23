# Code for Williamson et al., Oxygen availability drives blood traits and the cell number-size tradeoff across Andean hummingbirds
What drives variation in blood-O2 carrying capacity? We zoom in on one clade along a broad range of elevational gradients to get at causes of variation in detail. 

**Our preprint is now live on bioRxiv:** https://www.biorxiv.org/content/10.1101/2022.01.18.476833v1


# QUESTIONS 
In this paper we ask: 
**1)** To what extent do blood traits shift predictably with elevational changes in O2 availability, while controlling for other sources of variation?
**2)** Do blood traits vary with O2 availability differently over developmental (within species) and evolutionary (among species) time scales? 
**3)** How are blood traits affected by phylogenetic history and evolved species differences, such as body mass, metabolic intensity, and known genetic variants that affect Hb-O2 binding affinity?
**4)** To what extent are individual and species variation in O2-carrying capacity per unit volume of blood attributable to changes in cell number versus cell size? 
**5)** Does manifestation of the cell number-size tradeoff vary predictably among species?


# DATASET 
We collected specimen-vouchered hematological data from hummingbirds from 2006–2020 across a 4,578-meter elevational gradient in the Andes (Chile and Peru). Data for six blood traits (hemoglobin concentration ([Hb]), hematocrit (Hct), total red blood cell count (TRBC), mean cell volume (MCV), mean cell hemoglobin (MCH), and mean cell hemoglobin concentration (MCHC)) from 1,217 individuals of 77 species were used to model within- and among-species variation, respectively. All data are linked to vouchered specimens housed at the **Museum of Southwestern Biology** (MSB) at the University of New Mexico, the **Centro de Ornitología y Biodiversidad** in Peru, and the **Pontificia Universidad Católica de Chile** in Chile. Specimen records are accessible in the Artcos database (https://www.arctosdb.org). 

 
# R SCRIPTS

`HumBlood_DataWrangling.Rmd`: This script includes code for processing raw data, evaluating and eliminating outliers, creating new variables, combining data w/ Stotz, making individual elevational range adjustments, processing spatial data and generating sampling map, BioClim processing, reading out of final data file for modeling, plotting, and downstream analyses. Start here. 

`HumBlood_Phylogeny.Rmd`: Script to process McGuire et al. 2014 hummingbird phylogenies for our specific subsets. 

`HumBlood_PhyloElevRangeMapFig.Rmd`: Code and notes for making phylo elev range component of Figure 1 (note that final figure was constructed in Adobe Illustrator; script details when this happens).

`HumBlood_Modeling.Rmd`: Individual-level (within species analysis) modeling script. Read in processed data from "HumBlood_DataWrangling.Rmd", standardize predictors, make data subsets that match phylogeny subsets (so there are no NAs in any dataset for modeling), and then run checking, plotting sets of Bayesian phylogenetic mixed models w/ brms(). Script details all inndividual-level (within species) models. This script also contains the individual contmap components for each blood trait that were used to create **Figure S8** (i.e., trait-specific contmaps were created in R and combined in Illustrator for easy adjustment of graphics).

`HumBlood_Modeling_SpeciesMeans.Rmd`: Species mean (among species analysis) modeling script. Very similar to `HumBlood_Modeling.Rmd`, but this script outlines the procedure for processing data and running and evaluating species mean (among species) models. 

`HumBlood_CeNS.Rmd`: Cell Number Versus Size analysis script. Script contains all code and annotations required to reproduce our three-pronged cell number versus cell size analyses: 1) variation in cell number versus size within species, among species, and within single-species, respectively; 2) creation of the Cell Number versus Size Index (*CeNS*) using species-specific model coefficients; 3) analysis of the factors that contribute to variation in *CeNS*. This script also contains code for **Figure S1** and **Figure S4**. 

`HumBlood_Plots&Figs.Rmd`: This script contains most polished plots and figures for our all-hummingbird blood comparison. In most cases, I confined these graphics to this standalone script to not bog down analysis scripts with really beefy plots, and for ease of future revisions. Script contains code for **Figures 1-4** and **Figure S2**, **Figure S3**, **Figure S5**, **Figure S6**, **Figure S7**, **Figure S9**, *and* the entire breakpoint analysis associated with Figure S7. 
Questions? Feel free to email me at jwilliamson0110@gmail.com. 