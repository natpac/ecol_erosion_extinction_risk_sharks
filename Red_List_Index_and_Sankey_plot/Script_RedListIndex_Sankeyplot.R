rm(list=ls()) # removes all objects from the environment
cat("\014") # clears the console

#### Load packages ####
if(!require("pacman")) install.packages("pacman")
p_load(readxl, writexl, PNWColors, plyr, tidyr, ggalluvial)
`%notin%` <- Negate(`%in%`) # define a new operator '%notin%' as the negation of '%in%'
####

#### Load IUCN statuts df ####
status_global <- Table_IUCNstatus_traits <- as.data.frame(read_excel("General_data/Table_IUCNstatus_traits.xlsx"))
status_global <- status_global_withDD[status_global_withDD$y2020 %notin% c('DD', 'NE'),] # 155 DD or NE removed
####



#### Calculate Red List Index + bootstrap ####
## Calculate Red List Index with the uncertainty due to Data Deficient species where all of them are either assumed LC or CR
create_RLI_df <- function(subclass = 'all', status_df = status_global_withDD, uncertainty = T){
  if(!uncertainty){
    status_df <- status_df[status_df$y2020 %notin% c('DD', 'NE'),]
  }
  if(subclass != 'all'){
    status_global_temp <- status_df[status_df$subclass == subclass, ]
  } else{status_global_temp <- status_df}
  
  # Transform status data to a suitable format for RLI calculation
  status_sp <- as.data.frame(t(status_global_temp[,c('y1970', 'y1980', 'y2005', 'y2020')]))
  status_sp <- cbind(as.integer(sub('.', '', row.names(status_sp))),status_sp)
  colnames(status_sp) <- c('year', status_global_temp[,1])
  rownames(status_sp) <- NULL
  status_sp[,2:ncol(status_sp)] <- apply(status_sp[,2:ncol(status_sp)], 2, function(x){as.integer(mapvalues(x, from=c('CR(PE)', 'CR', 'EN', 'VU', 'NT', 'LC', 'DD'), to=c(4,4,3,2,1,0,666), warn_missing = F))})
  
  RLI_df <- data.frame('year'=status_sp$year, 'RLI'= NA_real_) # Initialize RLI dataframe
  if(uncertainty){
    # UPPER bound all DD are assumed to be LC
    status_sp_withDD_low <- status_sp
    status_sp_withDD_low[status_sp_withDD_low == 666] <- 0
    # LOWER Opposite is all DD are CR
    status_sp_withDD_high <- status_sp
    status_sp_withDD_high[status_sp_withDD_high == 666] <- 4
    # normal RLI no DD
    status_sp <- status_sp[, !apply(status_sp, 2, function(x) any(x == 666))]
    
    # Calculate RLI for each year
    RLI_df$RLI <- apply(RLI_df,1,function(df){
      xxx <- apply(as.data.frame(table(as.numeric(status_sp[status_sp$year == df['year'],2:ncol(status_sp)])), stringsAsFactors =F),2,as.numeric)
      1 - (sum(apply(xxx,1,prod)) / (5 * sum(xxx[,2])))
    })
    # Calculate lower and upper bounds
    RLI_df$RLIlow <- apply(RLI_df,1,function(df){
      xxx <- apply(as.data.frame(table(as.numeric(status_sp_withDD_low[status_sp_withDD_low$year == df['year'],2:ncol(status_sp_withDD_low)])), stringsAsFactors =F),2,as.numeric)
      1 - (sum(apply(xxx,1,prod)) / (5 * sum(xxx[,2])))
    })
    RLI_df$RLIhigh <- apply(RLI_df,1,function(df){
      xxx <- apply(as.data.frame(table(as.numeric(status_sp_withDD_high[status_sp_withDD_high$year == df['year'],2:ncol(status_sp_withDD_high)])), stringsAsFactors =F),2,as.numeric)
      1 - (sum(apply(xxx,1,prod)) / (5 * sum(xxx[,2])))
    })
  } else{
    # Calculate RLI without considering uncertainty
    RLI_df$RLI <- apply(RLI_df,1,function(df){
      xxx <- apply(as.data.frame(table(as.numeric(status_sp[status_sp$year == df['year'],2:ncol(status_sp)])), stringsAsFactors =F),2,as.numeric)
      1 - (sum(apply(xxx,1,prod)) / (5 * sum(xxx[,2])))
    })
  }
  print(RLI_df)
  return(RLI_df)
}

## Calculate Red List Index with the uncertainty due to potential classification error by assessors
create_RLI_df_with_bootstrap <- function(subclass = 'all', status_df = status_df, boot_nsp = boot_nsp, boot_type = c('all', 'NT_THR')[1]){
  if(subclass != 'all'){
    status_global_temp <- status_df[status_df$subclass == subclass, ]
  } else{status_global_temp <- status_df}
  
  # Transform status data to a suitable format for RLI calculation
  status_sp <- as.data.frame(t(status_global_temp[,c('y1970', 'y1980', 'y2005', 'y2020')]))
  status_sp <- cbind(as.integer(sub('.', '', row.names(status_sp))),status_sp)
  colnames(status_sp) <- c('year', status_global_temp[,1])
  rownames(status_sp) <- NULL
  status_sp[,2:ncol(status_sp)] <- apply(status_sp[,2:ncol(status_sp)], 2, function(x){as.integer(mapvalues(x, from=c('CR(PE)', 'CR', 'EN', 'VU', 'NT', 'LC'), to=c(4,4,3,2,1,0), warn_missing = F))})
  
  status_sp_boot=status_sp
  
  if(boot_type == 'all'){
    # Random sampling for bootstrapping
    boot_sp <- list(sample(2:ncol(status_sp), round(boot_nsp*(ncol(status_sp)-1)/100)) #1970
                    , sample(2:ncol(status_sp), round(boot_nsp*(ncol(status_sp)-1)/100)) #1980
                    , sample(2:ncol(status_sp), round(boot_nsp*(ncol(status_sp)-1)/100)) #2005
                    , sample(2:ncol(status_sp), round(boot_nsp*(ncol(status_sp)-1)/100)) #2020
    ) 
  }
  
  if(boot_type == 'NT_THR'){
    # Sampling based on threatened (Critically Endangered, Endangered, Vulnerable) and Near Threatened species resample
    sp_NT_THR <- list((2:ncol(status_sp))[status_df$y1970 != 'LC']
                      , (2:ncol(status_sp))[status_df$y1980 != 'LC']
                      , (2:ncol(status_sp))[status_df$y2005 != 'LC']
                      , (2:ncol(status_sp))[status_df$y2020 != 'LC']
    )
    boot_sp <- list(sample(sp_NT_THR[[1]], min(c(round(boot_nsp*(ncol(status_sp)-1)/100), length(sp_NT_THR[[1]])))) #1970
                    , sample(sp_NT_THR[[2]], min(c(round(boot_nsp*(ncol(status_sp)-1)/100), length(sp_NT_THR[[2]])))) #1980
                    , sample(sp_NT_THR[[3]], min(c(round(boot_nsp*(ncol(status_sp)-1)/100), length(sp_NT_THR[[3]])))) #2005
                    , sample(sp_NT_THR[[4]], min(c(round(boot_nsp*(ncol(status_sp)-1)/100), length(sp_NT_THR[[4]])))) #2020
    ) 
  }
  
  # Apply bootstrap modifications
  for(year in 1:nrow(status_sp_boot)){
    boot_spLC_year <- status_sp_boot[year, boot_sp[[year]]] %in% 0
    boot_spCR_year <- status_sp_boot[year, boot_sp[[year]]] %in% 4
    boot_sprest_year <- status_sp_boot[year, boot_sp[[year]]] %notin% c(0,4)
    
    if(boot_type == 'all'){
      status_sp_boot[year, boot_sp[[year]]][boot_spLC_year] <- sapply(status_sp_boot[year, boot_sp[[year]]][boot_spLC_year], function(x){x+1})
    }
    status_sp_boot[year, boot_sp[[year]]][boot_sprest_year] <- sapply(status_sp_boot[year, boot_sp[[year]]][boot_sprest_year], function(x){sample(c(x+1,x-1),1, prob = c(1,1))})
    status_sp_boot[year, boot_sp[[year]]][boot_spCR_year] <- sapply(status_sp_boot[year, boot_sp[[year]]][boot_spCR_year], function(x){x-1})
  }
  
  RLI_df_boot <- data.frame('year'=status_sp_boot$year, 'RLI_boot'= NA_real_) 
  RLI_df_boot$RLI_boot <- apply(RLI_df_boot,1,function(df){
    xxx <- apply(as.data.frame(table(as.numeric(status_sp_boot[status_sp_boot$year == df['year'],2:ncol(status_sp_boot)])), stringsAsFactors =F),2,as.numeric)
    1 - (sum(apply(xxx,1,prod)) / (5 * sum(xxx[,2])))
  })
  print(RLI_df_boot)
  return(RLI_df_boot)
}





## SAVE RLI values
# Create the Red List Index (RLI) for various subclasses and with/without uncertainty
RLI_df <- create_RLI_df(subclass = c('all', 'shark', 'ray', 'chimaera')[1], status_df = status_global_withDD, uncertainty = T)
RLI_df$RLI_shark <- create_RLI_df(subclass = c('all', 'shark', 'ray', 'chimaera')[2], status_df = status_global, uncertainty = F)$RLI
RLI_df$RLI_ray <- create_RLI_df(subclass = c('all', 'shark', 'ray', 'chimaera')[3], status_df = status_global, uncertainty = F)$RLI
if(any(status_global$subclass == 'chimaera')){
  RLI_df$RLI_ghost_shark <- create_RLI_df(subclass = c('all', 'shark', 'ray', 'chimaera')[4], status_df = status_global, uncertainty = F)$RLI
}

# Define bootstrap parameters
boot_nsp_vector=c(10,20) # number of boostrapped species in percent of total number of species # 10 or 20
nber_bootstrap=1000 # number of bootstrap
boot_type_vector = c('all', 'NT_THR')
bar = c('range','CI')[2]

# Loop for RLI with bootstrapping
for(boot_nsp in boot_nsp_vector){
  for(boot_type in boot_type_vector){
    RLI_df_bootstrap <- create_RLI_df_with_bootstrap(subclass = c('all', 'shark', 'ray', 'chimaera')[1], status_df = status_global, boot_nsp = boot_nsp, boot_type = boot_type)
    for(boot in 1:nber_bootstrap){
      RLI_df_bootstrap[,boot+1] <- create_RLI_df_with_bootstrap(subclass = c('all', 'shark', 'ray', 'chimaera')[1], status_df = status_global, boot_nsp = boot_nsp, boot_type = boot_type)$RLI_boot
    }
    colnames(RLI_df_bootstrap) <- c('year', paste0('RLI_', boot_type, '_', boot_nsp, 'sp_boot_', 1:nber_bootstrap))
    # RLI_df$RLI_boot_average <- t(apply(RLI_df_bootstrap[,2:ncol(RLI_df_bootstrap)],1,function(x)quantile(x, c(0.025,0.975))))
    
    if(bar == 'range'){
      RLI_df[, c(paste0('RLI_', boot_type, '_', boot_nsp, '_rangeLow'), paste0('RLI_', boot_type, '_', boot_nsp, '_rangeHigh'))] <- t(apply(RLI_df_bootstrap[,2:ncol(RLI_df_bootstrap)],1,function(x)range(x)))
    } else{
      RLI_df[, c(paste0('RLI_', boot_type, '_', boot_nsp, '_average_2.5%'), paste0('RLI_', boot_type, '_', boot_nsp, '_average_97.5%'))] <- t(apply(RLI_df_bootstrap[,2:ncol(RLI_df_bootstrap)],1,function(x)quantile(x, c(0.025,0.975))))
    }
  }
}

# Print final RLI data frame
RLI_df
####



#### Sankey plot ####
# Prepare data for Sankey plot
status_global_restricted <- status_global[, c('species', 'y1970', 'y1980', 'y2005', 'y2020')]
alluv_sp_status <- as.data.frame(pivot_longer(status_global_restricted, y1970:y2020))
colnames(alluv_sp_status) <- c('ID', 'year', 'status')
alluv_sp_status <- alluv_sp_status[, c('status', 'year', 'ID')]
alluv_sp_status$year <- as.numeric(gsub('y','', alluv_sp_status$year))
alluv_sp_status$ID <- as.numeric(as.factor(alluv_sp_status$ID))
alluv_sp_status$freq <- 1
alluv_sp_status[,1] <- as.factor(alluv_sp_status[,1])
alluv_sp_status[,1] <- factor(alluv_sp_status[,1], levels = c("LC", "NT", "VU", "EN", "CR", "CR(PE)"))
alluv_sp_status$ID <- as.integer(alluv_sp_status$ID)
# Define colors for IUCN categories
IUCN_colors <- c("#6C0F03", "#D81E05", "#FC7F3F", "#F9E814", "#CCE226", "#60C659") # CR(PE) color obtained between EX color and CR color, CR, EN, VU, NT, LC

# Create Sankey plot
ggplot_sankey <- ggplot(alluv_sp_status,
                        aes(x = year, stratum = status
                            , alluvium = ID
                            , fill = status
                        )) +
  geom_flow(stat = "alluvium", width = 4) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
  ) +
  scale_x_continuous(breaks = c(1970, 1980, 2005, 2020), labels = c('1970', '1980', '2005', '2020'), position = "top") +
  geom_stratum(width = 4) +
  geom_text(aes(label = status),
            stat = "stratum", size = 3) +
  scale_fill_manual(values=IUCN_colors, breaks=rev(levels(alluv_sp_status$status)), labels=rev(levels(alluv_sp_status$status)))

# Print Sankey plot
ggplot_sankey
####
