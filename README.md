# Comparative Hummingbird Blood
What drives variation in blood-oxygen carrying capacity? We zoom in on one clade along a broad range of elevational gradients to get at causes of variation in detail. 

R scripts: 

`HumBlood_DataWrangling.Rmd`: This script includes: processing raw data, evaluating and eliminating outliers, creating new variables, combining data w/ Stotz, making individual elevational range adjustments, processing spatial data and generating sampling map, BioClim processing, reading out of final data file for modeling, plotting, etc.

`HumBlood_Phylogeny.Rmd`: Script to process McGuire et al. 2014 hummingbird phylogenies for our specific subsets. 

`HumBlood_Modeling.Rmd`: This script includes: Reading in pre-processed data from "HumBlood_DataWrangling.Rmd", standardizing predictors, making data subsets that match subsetted phylogenies (so there are no NAs in any dataset for modeling), and then running, checking, plotting sets of Bayesian phylogenetic mixed models w/ brms(). 

`HumBlood_Plots&Figs.Rmd`: This script contains exploratory plots, polished plots, and figures for our all-hummingbird blood comparison. These are separated into a standalone script to not bog down main script.
