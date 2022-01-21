# Comparative Hummingbird Blood
What drives variation in blood-oxygen carrying capacity? We zoom in on one clade along a broad range of elevational gradients to get at causes of variation in detail. TEST EDIT 1/21/22


# QUESTIONS 

**1) How do blood parameters vary with abiotic aspects of the environment?** 
Predictors: 
1) elevation 
2) elevation position (position relative to the width of the elevational range; scale from 0-1)
3) temperature index (PC1 of 19 BioClim variables)
4) aridity index (PC1 of 19 BioClim variables)
5) latitude

**2) How to blood parameters vary with aspects of species biology?** 
Predictors: 
5) body mass
6) hemoglobin genotype (beta13-beta83 genotype) 

**3) Are patterns in blood phylogenetically constrained?**

**4) What relationships exist among blood parameters?**
6 total: [Hb], Hct, TRBC, MCV, MCH, MCHC


# DATASET 
Dataset includes blood data for 77 species of hummingbirds from Chile and Peru. Sampling was conducted from 2006-2020 by students, affiliates, researchers, and professors at the Museum of Southwestern Biology at the University of New Mexico. All data are linked to vouchered specimens*.
*Note: Because some tissues and specimens haven't been exported from Peru, and because Chile collection isn't yet catalogued, not all are archived in Arctos. 

# APPROACH 
We will build 6 sets of phylogenetic mixed models in brms (x 10 models per set)

Response variables: 
1) Hb
2) Hct
3) TRBC
4) MCV
5) MCH
6) MCHC

Predictor set for each set of models: 
1) elevation 
2) elevation position (position relative to the width of the elevational range; scale from 0-1)
3) temperature index
4) precip index
5) latitude 
6) mass
7) Hb genotype    
8) intraspecific variation: elevation
9) intraspecific variation: elevation position
10) intraspecific variation: latitude
11) intraspecific variation: mass
12) intraspecific vatiation: precip
13) intraspecific variation: temp 
+ (1|species)
+ (1|phylogeny)

For each response, we will include: 
1) Null model (intercept-only)
2) Random effects-only model (intercept + phylogenetic random effect + species random effect)
3) Full model with no random effects (just predictors) 
4) Full model with just phylo random effect
5) Full model with just species random effect 
6) Full model with phylo and species random effects 
7) Reduced model with no random effects (just predictors whose CIs do NOT overlap zero)
8) Reduced model with just phylo random effect 
9) Reduced model with just sepcies random effect 
10) Reduced model with phylogenetic and species random effects 

We'll use LOOIC and WAIC to compare model sets.



**R scripts:** 

`HumBlood_DataWrangling.Rmd`: This script includes: processing raw data, evaluating and eliminating outliers, creating new variables, combining data w/ Stotz, making individual elevational range adjustments, processing spatial data and generating sampling map, BioClim processing, reading out of final data file for modeling, plotting, etc.

`HumBlood_Phylogeny.Rmd`: Script to process McGuire et al. 2014 hummingbird phylogenies for our specific subsets. 

`HumBlood_Modeling.Rmd`: This script includes: Reading in pre-processed data from "HumBlood_DataWrangling.Rmd", standardizing predictors, making data subsets that match subsetted phylogenies (so there are no NAs in any dataset for modeling), and then running, checking, plotting sets of Bayesian phylogenetic mixed models w/ brms(). 

`HumBlood_Plots&Figs.Rmd`: This script contains exploratory plots, polished plots, and figures for our all-hummingbird blood comparison. These are separated into a standalone script to not bog down main script.
