Ecological erosion and expanding extinction risk of sharks and rays
In this repository, you will find the R codes to support the analyses in the publication.
Cite this article as XXCITATIONTOADDXX. DOI: XXDOITOADDXX READ THE FULL ARTICLE AT XXLINKTOADDXX

Time-series and trait data are publicly available on https://www.sharkipedia.org/ XXADDSharkipediacitationXX
XXADD description of folder and analysisXX



Boosted_Regression_Tree_Analysis is the folder for the Boosted Regression Tree analysis. It contains:
- Source code contains with Libraries and User-defined functions
- Data (Red List Index ([RLI] and covariates)
- Analysis script

Functional_richness_analysis is the folder for the Simulations of the loss of functional diversity analysis. It contains:
- a script to perform the functional richness analysis
- two functions

General_data is the folder with species names, IUCN categories and covariates.

Global_Chondrichthyan_Catch_and_Effort is the folder for the catch reconstructions and fishing effort analysis. It conatins:
- a script to perform the analysis
- a stan model
- two dataset of global catch and effort

Red_List_Index_and_Sankey_plot is the folder for the Red List Index (RLI) over time analysis. It contains:
- a script to calculate the RLI (as well as bootstrapped uncertainty) and do a Sankey plot
- RLI values

Red_List_Index_and_Spatial_Units is the folder for the Red List Index (RLI) over spatial units analysis. It contains:
- .Rmd script to calculate the RLI values by spatial units (equal area hex grid, EEZ, LME, MEOW, FAO)
- .csv files containing the output RLI values 
- Note: script uses package hosted at https://github.com/Marine-Biodiversity-Conservation-Lab/RLIspatial.git for calculating RLI. 
