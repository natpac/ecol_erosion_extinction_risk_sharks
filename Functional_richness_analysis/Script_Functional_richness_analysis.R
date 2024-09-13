rm(list=ls()) # removes all objects from the environment
cat("\014") # clears the console

#### Load packages and functions ####
if(!require("pacman")) install.packages("pacman")
p_load(mFD, ape, readxl, parallel, foreach, doFuture, progressr, scales, ggplot2, xlsx)
`%notin%` <- Negate(`%in%`)

source('./show_final_pcoa_fun.R')
source('./compute_arrows_fun.R')
####

#### Load IUCN statuts df ####
status_global <- as.data.frame(read_excel("General_data/Table_IUCNstatus_traits.xlsx"))
status_global <- status_global[!is.na(as.numeric(status_global$Range_area)),] # REMOVE SP WITH NO RANGE DISTRIB
status_global <- status_global[status_global$Repro_mode != "unknown",] # REMOVE SP WITH UNKNOWN RepMode

## Removing unknown_habitat (species that do not have IUCN habitat coded)
SP_unknown_habitat <- c("Parascyllium elongatum","Chiloscyllium caeruleopunctatum", "Akheilos suwartanai", "Scyliorhinus garmani", "Dentiraja flindersi", "Dipturus amphispinus", "Dipturus ecuadoriensis", "Dipturus wuhanlingi", "Okamejei mengae", "Torpedo suessii", "Squatina caillieti", "Urobatis marmoratus", "Spinilophus armatus") 
status_global <- status_global[status_global$species %notin% SP_unknown_habitat, ]

## Make Repro_mode ordered factor column
status_global$Repro_mode <- ordered(status_global$Repro_mode, levels = c("egg-laying_NoMatrotrophy", "Live-bearing_NoMatrotrophy", "Live-bearing_Matrotrophy"))
####


#### PREP TRAIT DF ####
# Select relevant columns
df_shark_trait <- status_global[,
                                c(3 # "species"
                                  , 14 # "size"
                                  , 15 # "Range_area"
                                  , 16 # "Repro_mode"
                                  , 17:45 # IUCN habitats
                                )
]

# Set row names of 'df_shark_trait' to be the 'species' column and then remove the 'species' column
row.names(df_shark_trait) <- df_shark_trait[,1]
df_shark_trait <- df_shark_trait[,-1]

# Convert 'size' and 'Range_area' columns to numeric
df_shark_trait[,"size"] <- as.numeric(df_shark_trait[,"size"])
df_shark_trait[,"Range_area"] <- as.numeric(df_shark_trait[,"Range_area"])

# rescale col in fuzzy (sum = 1)
df_shark_trait[, 4:32] <- df_shark_trait[, 4:32]/rowSums(df_shark_trait[, 4:32])

# # Normality test for quantitative data
# hist(df_shark_trait[,"size"], col='steelblue', breaks = 100)
# qqnorm(df_shark_trait[,"size"])
# qqline(df_shark_trait[,"size"])
# shapiro.test(df_shark_trait[,"size"]) # p-value = 2.2e-16 NOT OK
# 
# hist(log(df_shark_trait[,"size"]), col='steelblue', breaks = 100)
# qqnorm(log(df_shark_trait[,"size"]))
# qqline(log(df_shark_trait[,"size"]))
# shapiro.test(log(df_shark_trait[,"size"])) # p-value = 9.136e-12 NOT OK BUT GOOD ENOUGH
# 
# # Normality test for quantitative data
# hist(df_shark_trait[,"Range_area"], col='steelblue', breaks = 100)
# qqnorm(df_shark_trait[,"Range_area"])
# qqline(df_shark_trait[,"Range_area"])
# shapiro.test(df_shark_trait[,"Range_area"]) # p-value = 2.2e-16 NOT OK
# 
# hist(log(df_shark_trait[,"Range_area"]), col='steelblue', breaks = 100)
# qqnorm(log(df_shark_trait[,"Range_area"]))
# qqline(log(df_shark_trait[,"Range_area"]))
# shapiro.test(log(df_shark_trait[,"Range_area"])) # p-value = 8.808e-14 NOT OK BUT GOOD ENOUGH

# SIZE BECOMES LOG (LN) SIZE
df_shark_trait[,"size"] <- log(df_shark_trait[,"size"])
colnames(df_shark_trait)[colnames(df_shark_trait) == 'size'] <- 'log_size'

# RANGE_AREA BECOMES LOG (LN) RANGE_AREA
df_shark_trait[,"Range_area"] <- log(df_shark_trait[,"Range_area"])
colnames(df_shark_trait)[colnames(df_shark_trait) == 'Range_area'] <- 'log_Range_area'


# NORMALIZED NUMERIC VARIABLES
df_shark_trait[,"log_size"] <- (df_shark_trait[,"log_size"]-min(df_shark_trait[,"log_size"]))/(max(df_shark_trait[,"log_size"])-min(df_shark_trait[,"log_size"]))
df_shark_trait[,"log_Range_area"] <- (df_shark_trait[,"log_Range_area"]-min(df_shark_trait[,"log_Range_area"]))/(max(df_shark_trait[,"log_Range_area"])-min(df_shark_trait[,"log_Range_area"]))

## PREP CAT FUZZY DF
shark_trait_cat <- data.frame(names(df_shark_trait)
                              , c(rep("Q", 2), "O", rep("F", 29))
                              , c(rep(NA_character_, 3), rep("Habitat", 29)) # only fuzzy have a name
)
colnames(shark_trait_cat) <- c("trait_name", "trait_type", "fuzzy_name")
####


#### START THE ANALYSIS ####
# Compute trait-based distances using the 'funct.dist' function
dist_shark <- funct.dist(
  sp_tr         = df_shark_trait,
  tr_cat        = shark_trait_cat,
  metric        = "gower",
  scale_euclid  = "noscale", # "scale_center"
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# use quality.fpscaes function to compute quality metrics:
quality_fspaces_shark <- quality.fspaces(
  sp_dist             = dist_shark,
  fdendro             = "average",
  maxdim_pcoa         = 10,
  deviation_weighting = "squared",
  fdist_scaling       = F)


# retrieve the functional space associated with minimal quality metric: 
quality_fspaces_shark$"quality_fspaces"

# apply(quality_fspaces_shark$quality_fspaces, 2, which.min)
nber_axis <- 5

# Extract PCoA eigenvalues; ensure there are no negative eigenvalues
quality_fspaces_shark$"details_fspaces"$"pc_eigenvalues"

# Extract the coordinates of species in the PCoA space
sp_faxes_coord_shark <- quality_fspaces_shark$details_fspaces$sp_pc_coord

## INFO ON PCOA AXIS
# Compute correlation between traits and PCoA axes
axis_correl <- traits.faxes.cor(
  sp_tr = df_shark_trait,
  sp_faxes_coord = sp_faxes_coord_shark[,paste0("PC",1:nber_axis)],
  plot = F,
)

# Filter correlations with p-values < 0.05 and values > 0.1
axis_correl_red <- axis_correl[axis_correl$p.value <0.05 & axis_correl$value > 0.1,]
axis_correl_red <- axis_correl_red[order(axis_correl_red$axis, axis_correl_red$value, decreasing = c(FALSE, TRUE), method="radix"),]

# Perform Principal Coordinates Analysis (PCA) on the distance matrix
pcoa_trdist <- pcoa(dist_shark) 

# Extract the first two principal coordinates (PC1 and PC2) and their variance explained
vec.pcoa <- pcoa_trdist$vectors[ ,1:2]
var.pcoa <- round(pcoa_trdist$values$Relative_eig[1:2],2)

# Convert 'Repro_mode' to numeric for further analysis
df_shark_trait_quantitative <- df_shark_trait
df_shark_trait_quantitative[,"Repro_mode"] <- as.numeric(df_shark_trait_quantitative[,"Repro_mode"])-1

# Compute arrows for PCoA plot
basic_arrow_pcoa = compute_arrows_fun(given_pcoa=pcoa_trdist, trait_df=df_shark_trait_quantitative)
# Plot PCoA with arrows for different axis combinations
show_final_pcoa_fun(basic_arrow_pcoa, sp_plot=F, axis1 = 1, axis2 = 2, max.overlaps_label = 10) # , variable_scalar = 15
show_final_pcoa_fun(basic_arrow_pcoa, sp_plot=F, axis1 = 1, axis2 = 3, max.overlaps_label = 10) # , variable_scalar = 15
show_final_pcoa_fun(basic_arrow_pcoa, sp_plot=F, axis1 = 2, axis2 = 3, max.overlaps_label = 10) # , variable_scalar = 15


# Generate a big plot showing the functional space
big_plot <- funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord_shark[ , paste0("PC",1:min(nber_axis,4))], # , "PC4", "PC5"
  faxes           = paste0("PC",1:min(nber_axis,4)), # , "PC5"
  alpha_ch        = 0.5,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
)
# Display the plot using 'patchwork'
big_plot$patchwork

# Create and plot combinations of PCoA axes
combination_axis_pcoa <- combn(1:nber_axis, 2, simplify = F)
for(comb_axis in 1:length(combination_axis_pcoa)){
  show_final_pcoa_fun(basic_arrow_pcoa, sp_plot=F, axis1 = combination_axis_pcoa[[comb_axis]][1], axis2 = combination_axis_pcoa[[comb_axis]][2], max.overlaps_label = 10)
}
####


#### create extinction scenarios ####
nber_axis <- 5 # Number of PCoA axes to use
bootstrap_rep_nber <- 500 # Number of bootstrap replications
nber_extinctin_scenar <- nrow(df_shark_trait)-nber_axis # Number of species to simulate extinction

### -> RANDOM
# Initialize an array to hold the extinction scenarios for each bootstrap replication
asb_sp_w_extinction_scenar <- array(1, dim = c(nber_extinctin_scenar,nrow(df_shark_trait),bootstrap_rep_nber))
dimnames(asb_sp_w_extinction_scenar)[[1]] <-  paste0("Ext_",0:(nber_extinctin_scenar-1))
dimnames(asb_sp_w_extinction_scenar)[[2]] <- row.names(df_shark_trait)
dimnames(asb_sp_w_extinction_scenar)[[3]] <-  paste0("boot_",1:bootstrap_rep_nber)

# List to store the order of species to be extinct in each bootstrap replicationlist_order_extinct_scenar <- list()
for(bootstrap_nber in 1:bootstrap_rep_nber){
  # Randomly sample species to be extinct
  order_extinct_scenar <- sample(1:nrow(df_shark_trait), nber_extinctin_scenar)
  list_order_extinct_scenar[[bootstrap_nber]] <- order_extinct_scenar
  # Update the extinction scenarios with the selected extinct species
  for(extinct_scenar in 2:length(order_extinct_scenar)){
    asb_sp_w_extinction_scenar[extinct_scenar,order_extinct_scenar[1:(extinct_scenar-1)],bootstrap_nber] <- 0
  }
} 


# Parallel processing setup
registerDoFuture()  ## %dopar% parallelizes via future
plan(multisession, workers = availableCores()-1)
handlers(global = TRUE)
handlers("cli")
options(future.globals.maxSize= 20000*1024^2)

# Function to compute functional richness for each bootstrap replication
bootstrap_FRic_random <- function(xs,sp_faxes_coord_shark,list_order_extinct_scenar) {
  p <- progressor(along = xs) # Progress bar for tracking computation
  foreach(bootstrap_nber = 1:bootstrap_rep_nber, .combine = cbind) %dopar% {
    # Initialize a temporary array for each bootstrap replication
    asb_sp_w_extinction_scenar_tmp <- array(1, dim = c(nber_extinctin_scenar,nrow(df_shark_trait)))
    dimnames(asb_sp_w_extinction_scenar_tmp)[[1]] <-  paste0("Ext_",0:(nber_extinctin_scenar-1))
    dimnames(asb_sp_w_extinction_scenar_tmp)[[2]] <- row.names(df_shark_trait)
    # Get the current extinction scenario
    order_extinct_scenar <- list_order_extinct_scenar[[bootstrap_nber]]
    for(extinct_scenar in 2:length(order_extinct_scenar)){
      asb_sp_w_extinction_scenar_tmp[extinct_scenar,order_extinct_scenar[1:(extinct_scenar-1)]] <- 0
    }
    # Compute alpha functional diversity indices (FRic) for the current scenario
    alpha_fr <- mFD::alpha.fd.multidim(
      sp_faxes_coord = sp_faxes_coord_shark[,paste0("PC",1:nber_axis)],
      asb_sp_w = asb_sp_w_extinction_scenar_tmp,
      ind_vect = "fric",
      scaling = TRUE,
      check_input = TRUE,
      details_returned = TRUE) # compute alpha functional diversity indices
    p(message = sprintf("Bootstrap %g", bootstrap_nber)) # Update progress
    return(alpha_fr$functional_diversity_indices$fric) # Return functional richness (FRic)
  }
} 


# Perform the bootstrap analysis to compute functional richness
df_bootstrap_FRic_RANDOM = bootstrap_FRic_random(xs = 1:bootstrap_rep_nber,sp_faxes_coord_shark,list_order_extinct_scenar)
colnames(df_bootstrap_FRic_RANDOM) <- paste0("boot_",1:bootstrap_rep_nber)
row.names(df_bootstrap_FRic_RANDOM) <- paste0("Ext_",0:(nber_extinctin_scenar-1))

# Plotting
plot(as.numeric(df_bootstrap_FRic_RANDOM[,1])
     , xlab = 'nber of species'
     , ylab = 'FRic'
     , col = alpha('#a6cee3', 0.3)
     , xaxt = 'n'
     , yaxt = 'n'
     , type = 'l')
axis(1, at=c(0, nrow(df_bootstrap_FRic_RANDOM)-1000, nrow(df_bootstrap_FRic_RANDOM)-500, nrow(df_bootstrap_FRic_RANDOM)), labels=rev(c(0, 500, 1000, nrow(df_bootstrap_FRic_RANDOM))))
axis(2, at=seq(0,1,by=.25), labels=seq(0,100,by=25), las = 1)

# Add lines for each bootstrap replication to the plot
for(bootstrap_nber in 1:bootstrap_rep_nber){
  lines(as.numeric(df_bootstrap_FRic_RANDOM[,bootstrap_nber])
        , col = alpha('#a6cee3', 0.3)
  )
}
###


### Following IUCN extinction
IUCN_cat_all <- c('CR(PE)', 'CR', 'EN', 'VU', 'NT', 'LC', 'DD')

nber_axis <- 5
bootstrap_rep_nber <- 500
nber_extinctin_scenar <- nrow(df_shark_trait)-nber_axis
asb_sp_w_extinction_scenar_IUCN <- array(1, dim = c(nber_extinctin_scenar,nrow(df_shark_trait),bootstrap_rep_nber))
dimnames(asb_sp_w_extinction_scenar_IUCN)[[1]] <-  paste0("Ext_",0:(nber_extinctin_scenar-1))
dimnames(asb_sp_w_extinction_scenar_IUCN)[[2]] <- row.names(df_shark_trait)
dimnames(asb_sp_w_extinction_scenar_IUCN)[[3]] <-  paste0("boot_",1:bootstrap_rep_nber)

list_order_extinct_scenar_IUCN <- list()
for(bootstrap_nber in 1:bootstrap_rep_nber){
  order_extinct_scenar_IUCN <- unlist(lapply(IUCN_cat_all, function(x){sample(which(status_global$y2020 == x), if(x != 'DD'){sum(status_global$y2020 == x)}else{sum(status_global$y2020 == x)-nber_axis})}))
  list_order_extinct_scenar_IUCN[[bootstrap_nber]] <- order_extinct_scenar_IUCN
  for(extinct_scenar in 2:length(order_extinct_scenar_IUCN)){
    asb_sp_w_extinction_scenar_IUCN[extinct_scenar,order_extinct_scenar_IUCN[1:(extinct_scenar-1)],bootstrap_nber] <- 0
  }
} 


registerDoFuture()  ## %dopar% parallelizes via future
plan(multisession) # availableCores() = 8
handlers(global = TRUE)
handlers("cli")
options(future.globals.maxSize= 20000*1024^2)

bootstrap_FRic_IUCN <- function(xs,sp_faxes_coord_shark,list_order_extinct_scenar_IUCN) {
  p <- progressor(along = xs)
  foreach(bootstrap_nber = 1:bootstrap_rep_nber, .combine = cbind) %dopar% {
    asb_sp_w_extinction_scenar_tmp <- array(1, dim = c(nber_extinctin_scenar,nrow(df_shark_trait)))
    dimnames(asb_sp_w_extinction_scenar_tmp)[[1]] <-  paste0("Ext_",0:(nber_extinctin_scenar-1))
    dimnames(asb_sp_w_extinction_scenar_tmp)[[2]] <- row.names(df_shark_trait)
    order_extinct_scenar <- list_order_extinct_scenar_IUCN[[bootstrap_nber]]
    for(extinct_scenar in 2:length(order_extinct_scenar)){
      asb_sp_w_extinction_scenar_tmp[extinct_scenar,order_extinct_scenar[1:(extinct_scenar-1)]] <- 0
    }
    alpha_fr <- mFD::alpha.fd.multidim(
      sp_faxes_coord = sp_faxes_coord_shark[,paste0("PC",1:nber_axis)],
      asb_sp_w = asb_sp_w_extinction_scenar_tmp,
      ind_vect = "fric",
      scaling = TRUE,
      check_input = TRUE,
      details_returned = TRUE) #compute alpha functional diversity indices
    p(message = sprintf("Bootstrap %g", bootstrap_nber))
    return(alpha_fr$functional_diversity_indices$fric)
  }
} 


df_bootstrap_FRic_IUCN = bootstrap_FRic_IUCN(xs = 1:bootstrap_rep_nber,sp_faxes_coord_shark,list_order_extinct_scenar_IUCN)
colnames(df_bootstrap_FRic_IUCN) <- paste0("boot_",1:bootstrap_rep_nber)
row.names(df_bootstrap_FRic_IUCN) <- paste0("Ext_",0:(nber_extinctin_scenar-1))


## PLOT
plot(as.numeric(df_bootstrap_FRic_IUCN[,1])
     , xlab = 'nber of species'
     , ylab = 'FRic'
     , col = alpha('#a6cee3', 0.3)
     , xaxt = 'n'
     , yaxt = 'n'
     , type = 'l')
axis(1, at=c(0, nrow(df_bootstrap_FRic_IUCN)-1000, nrow(df_bootstrap_FRic_IUCN)-500, nrow(df_bootstrap_FRic_IUCN)), labels=rev(c(0, 500, 1000, nrow(df_bootstrap_FRic_IUCN))))
axis(2, at=seq(0,1,by=.25), labels=seq(0,100,by=25), las = 1)

for(bootstrap_nber in 1:bootstrap_rep_nber){
  lines(as.numeric(df_bootstrap_FRic_IUCN[,bootstrap_nber])
        , col = alpha('#a6cee3', 0.3)
  )
}
###



### -> Following IUCN extinction with DD = CR
sp_faxes_coord_shark <- quality_fspaces_shark$details_fspaces$sp_pc_coord
IUCN_cat_all <- c('CR(PE)', 'CR', 'EN', 'VU', 'NT', 'LC')

nber_axis <- 5
bootstrap_rep_nber <- 500
nber_extinctin_scenar <- nrow(df_shark_trait)-nber_axis
asb_sp_w_extinction_scenar_IUCN <- array(1, dim = c(nber_extinctin_scenar,nrow(df_shark_trait),bootstrap_rep_nber))
dimnames(asb_sp_w_extinction_scenar_IUCN)[[1]] <-  paste0("Ext_",0:(nber_extinctin_scenar-1))
dimnames(asb_sp_w_extinction_scenar_IUCN)[[2]] <- row.names(df_shark_trait)
dimnames(asb_sp_w_extinction_scenar_IUCN)[[3]] <-  paste0("boot_",1:bootstrap_rep_nber)

# Update the conservation status by treating 'DD' (Data Deficient) as 'CR' (Critically Endangered)
status_global$y2020DDCR <- status_global$y2020
status_global$y2020DDCR <- gsub('DD', 'CR', status_global$y2020DDCR)
list_order_extinct_scenar_IUCN <- list()
for(bootstrap_nber in 1:bootstrap_rep_nber){
  order_extinct_scenar_IUCN <- unlist(lapply(IUCN_cat_all, function(x){sample(which(status_global$y2020DDCR == x), if(x != 'LC'){sum(status_global$y2020DDCR == x)}else{sum(status_global$y2020DDCR == x)-nber_axis})}))
  list_order_extinct_scenar_IUCN[[bootstrap_nber]] <- order_extinct_scenar_IUCN
  for(extinct_scenar in 2:length(order_extinct_scenar_IUCN)){
    asb_sp_w_extinction_scenar_IUCN[extinct_scenar,order_extinct_scenar_IUCN[1:(extinct_scenar-1)],bootstrap_nber] <- 0
  }
} 


registerDoFuture()  ## %dopar% parallelizes via future
plan(multisession) # availableCores() = 8
handlers(global = TRUE)
handlers("cli")
options(future.globals.maxSize= 20000*1024^2)

bootstrap_FRic_IUCN <- function(xs,sp_faxes_coord_shark,list_order_extinct_scenar_IUCN) {
  p <- progressor(along = xs)
  foreach(bootstrap_nber = 1:bootstrap_rep_nber, .combine = cbind) %dopar% {
    asb_sp_w_extinction_scenar_tmp <- array(1, dim = c(nber_extinctin_scenar,nrow(df_shark_trait)))
    dimnames(asb_sp_w_extinction_scenar_tmp)[[1]] <-  paste0("Ext_",0:(nber_extinctin_scenar-1))
    dimnames(asb_sp_w_extinction_scenar_tmp)[[2]] <- row.names(df_shark_trait)
    order_extinct_scenar <- list_order_extinct_scenar_IUCN[[bootstrap_nber]]
    for(extinct_scenar in 2:length(order_extinct_scenar)){
      asb_sp_w_extinction_scenar_tmp[extinct_scenar,order_extinct_scenar[1:(extinct_scenar-1)]] <- 0
    }
    alpha_fr <- mFD::alpha.fd.multidim(
      sp_faxes_coord = sp_faxes_coord_shark[,paste0("PC",1:nber_axis)],
      asb_sp_w = asb_sp_w_extinction_scenar_tmp,
      ind_vect = "fric",
      scaling = TRUE,
      check_input = TRUE,
      details_returned = TRUE) #compute alpha functional diversity indices
    p(message = sprintf("Bootstrap %g", bootstrap_nber))
    return(alpha_fr$functional_diversity_indices$fric)
  }
} 

df_bootstrap_FRic_IUCN = bootstrap_FRic_IUCN(xs = 1:bootstrap_rep_nber,sp_faxes_coord_shark,list_order_extinct_scenar_IUCN)
colnames(df_bootstrap_FRic_IUCN) <- paste0("boot_",1:bootstrap_rep_nber)
row.names(df_bootstrap_FRic_IUCN) <- paste0("Ext_",0:(nber_extinctin_scenar-1))




### -> Following IUCN extinction with DD = LC
load("./Analysis.R/Functional_richness/dist-quality_fspaces.RData")
sp_faxes_coord_shark <- quality_fspaces_shark$details_fspaces$sp_pc_coord
IUCN_cat_all <- c('CR(PE)', 'CR', 'EN', 'VU', 'NT', 'LC')

nber_axis <- 5
bootstrap_rep_nber <- 500
nber_extinctin_scenar <- nrow(df_shark_trait)-nber_axis
asb_sp_w_extinction_scenar_IUCN <- array(1, dim = c(nber_extinctin_scenar,nrow(df_shark_trait),bootstrap_rep_nber))
dimnames(asb_sp_w_extinction_scenar_IUCN)[[1]] <-  paste0("Ext_",0:(nber_extinctin_scenar-1))
dimnames(asb_sp_w_extinction_scenar_IUCN)[[2]] <- row.names(df_shark_trait)
dimnames(asb_sp_w_extinction_scenar_IUCN)[[3]] <-  paste0("boot_",1:bootstrap_rep_nber)

# Update the conservation status by treating 'DD' (Data Deficient) as 'LC' (Least Concern)
status_global$y2020DDLC <- status_global$y2020
status_global$y2020DDLC <- gsub('DD', 'LC', status_global$y2020DDLC)
list_order_extinct_scenar_IUCN <- list()
for(bootstrap_nber in 1:bootstrap_rep_nber){
  order_extinct_scenar_IUCN <- unlist(lapply(IUCN_cat_all, function(x){sample(which(status_global$y2020DDLC == x), if(x != 'LC'){sum(status_global$y2020DDLC == x)}else{sum(status_global$y2020DDLC == x)-nber_axis})}))
  list_order_extinct_scenar_IUCN[[bootstrap_nber]] <- order_extinct_scenar_IUCN
  for(extinct_scenar in 2:length(order_extinct_scenar_IUCN)){
    asb_sp_w_extinction_scenar_IUCN[extinct_scenar,order_extinct_scenar_IUCN[1:(extinct_scenar-1)],bootstrap_nber] <- 0
  }
} 


registerDoFuture()  ## %dopar% parallelizes via future
plan(multisession) # availableCores() = 8
handlers(global = TRUE)
handlers("cli")
options(future.globals.maxSize= 20000*1024^2)

bootstrap_FRic_IUCN <- function(xs,sp_faxes_coord_shark,list_order_extinct_scenar_IUCN) {
  p <- progressor(along = xs)
  foreach(bootstrap_nber = 1:bootstrap_rep_nber, .combine = cbind) %dopar% {
    asb_sp_w_extinction_scenar_tmp <- array(1, dim = c(nber_extinctin_scenar,nrow(df_shark_trait)))
    dimnames(asb_sp_w_extinction_scenar_tmp)[[1]] <-  paste0("Ext_",0:(nber_extinctin_scenar-1))
    dimnames(asb_sp_w_extinction_scenar_tmp)[[2]] <- row.names(df_shark_trait)
    order_extinct_scenar <- list_order_extinct_scenar_IUCN[[bootstrap_nber]]
    for(extinct_scenar in 2:length(order_extinct_scenar)){
      asb_sp_w_extinction_scenar_tmp[extinct_scenar,order_extinct_scenar[1:(extinct_scenar-1)]] <- 0
    }
    alpha_fr <- mFD::alpha.fd.multidim(
      sp_faxes_coord = sp_faxes_coord_shark[,paste0("PC",1:nber_axis)],
      asb_sp_w = asb_sp_w_extinction_scenar_tmp,
      ind_vect = "fric",
      scaling = TRUE,
      check_input = TRUE,
      details_returned = TRUE) #compute alpha functional diversity indices
    p(message = sprintf("Bootstrap %g", bootstrap_nber))
    return(alpha_fr$functional_diversity_indices$fric)
  }
} 


df_bootstrap_FRic_IUCN = bootstrap_FRic_IUCN(xs = 1:bootstrap_rep_nber,sp_faxes_coord_shark,list_order_extinct_scenar_IUCN)
colnames(df_bootstrap_FRic_IUCN) <- paste0("boot_",1:bootstrap_rep_nber)
row.names(df_bootstrap_FRic_IUCN) <- paste0("Ext_",0:(nber_extinctin_scenar-1))
####