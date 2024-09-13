# This is the source code for the Red List Index Boosted Regression Tree Analysis
# Source code contains:
  # 1. Libraries
  # 2. User-defined functions

# 1. Load libraries -----------------------------------------------------------------------
library(tidyverse) # data wrangling
library(countrycode) # add ISO codes to each model
library(patchwork) # patch plots together
library(reshape2) # to create collinearity heatmap
library(xgboost) # run BRT analyses
library(pdp) # calculate partial dependence in BRTs

# 2. Custom functions -----------------------------------------------------------------------
# calculate the total number of NAs in a dataframe
calc_na <- 
  function(data_frame) {
    
    df <- 
      lapply(data_frame, function(x) sum(is.na(x)))
    
    return(df)
    
  }

# publication theme
publication_theme <- function(axis_title_size = 13,
                              axis_text_size = 11) {
  
  theme(text = element_text(family = 'Helvetica'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(colour = 'grey20', size = 12),
        axis.line = element_line(colour = 'grey40'),
        axis.text = element_text(colour = 'grey40', size = axis_text_size),
        axis.title = element_text(colour = 'grey20', size = axis_title_size),
        axis.ticks = element_line(colour = 'grey40'),
        strip.background = element_blank())
  
}
