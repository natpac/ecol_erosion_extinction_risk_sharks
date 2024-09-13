# This script does the BRT analysis for the RLI data

source('Source.R')

# Need to change back to including CPUE
rli_data <- 
  read_csv('../Data/RLI_covariates.csv') %>% 
  # Here are the variables that we're dropping
  dplyr::select(-catch_cv, -catch, -effort, -species_richness,
                -latitude, -HDI_2021, -finsHK_2010)

names(rli_data)
str(rli_data)

# let's check the skewness of the data
pairs_skew <- 
  lapply(colnames(rli_data)[4:12], function(var_name) {
    
    the_plot <- 
      rli_data %>% 
      # Pull out the only column we care about
      dplyr::select(all_of(var_name)) %>% 
      # Remove NA values
      drop_na() %>% 
      # Rename to make it easier to plot
      dplyr::rename('variable_column' = 1) %>% 
      # make the plot
      ggplot(aes(x = variable_column)) +
      geom_histogram(bins = 30) +
      labs(x = paste(var_name)) +
      theme_light()
    
    return(the_plot)
    
  })

# Bind all the plots together to look at their skewness
pairs_plot <- 
  pairs_skew[[1]] + pairs_skew[[2]] + pairs_skew[[3]] + pairs_skew[[4]] + 
  pairs_skew[[5]] + pairs_skew[[6]] + pairs_skew[[7]] + pairs_skew[[8]] + 
  pairs_skew[[9]] + 
  plot_layout(ncol = 4)

# From here, we'll figure out which variables we'll need to log
pairs_plot

# OK, so looks like we'll have to log cpue, coastpop_2020, 
  # area_km2, subsidy_good, subsidy_bad, GDP_2019
rli_data <- 
  rli_data %>% 
  mutate_at(vars(c('area_km2', 'GDP_2019', 'cpue')), log) %>% 
  mutate_at(vars(c('coastpop_2020', 'subsidy_good', 'subsidy_bad')), log1p)

# RUN BRT -----------------------------------------------------------------------
# Prep BRT objects
# Specify the model formula
fmod <- formula(~ coastpop_2020 + PrimaryProduction + area_km2 + cpue + 
                  SST + subsidy_good + subsidy_bad + GDP_2019 + WGI_2018)

# create model matrix
modmat <- 
  stats::model.matrix(fmod, model.frame(~ ., rli_data, na.action = na.pass))[, -1]

head(modmat)

# create a vector of response
labels <- rli_data$RLI_2020_weighted

## Tune hyperparameters ----------------------------------------------------------
# We're going to find the optimal combination of hyperparameters to minimize the
# loss function

# Create a matrix of tuning variables
tune_grid <- 
  expand.grid(eta = seq(0.1, 0.9, 0.1), # learning rate
              gamma = seq(0.1, 0.9, 0.1), # minimum loss reduction
              max_depth = c(5, 10, 15), # maximum tree depth
              subsample = seq(0.1, 0.9, 0.1), # subsample ratio of the training instance
              rmse = NA) # empty column to infill

# Run a loop through all the combinations of hyperparameters to find the
# optimal combination
for(i in 1:nrow(tune_grid)) {
  
  # create a list of parameters
  tune_params <- 
    list(eta = tune_grid$eta[i],
         gamma = tune_grid$gamma[i],
         max_depth = tune_grid$max_depth[i],
         subsample = tune_grid$subsample[i])
  
  # run the model with specific set of hyperparameters
  tune_brt <- 
    xgboost::xgboost(modmat, label = labels, nrounds = 150,
                     params = tune_params,
                     objective = 'reg:logistic',
                     verbose = 0)
  
  # infill the dataframe
  tune_grid$rmse[i] <- min(tune_brt$evaluation_log$train_rmse)
  
  # Keep track of progress
  cat(paste('Bootstrapping the model, round ', i, '/', nrow(tune_grid), sep = ''), '\n')
  
}

# find the best combinations of hyperparameters
tune_grid %>% dplyr::slice_max(rmse) 
tune_grid %>% dplyr::slice_min(rmse) 

# For CPUE:
    # Max: 0.1110837
    # Min: 0.04187448
    # eta = 0.9; gamma = 0.1; max_depth = 10; subsample = 0.9

## Full model ----------------------------------------------------------------------
# Create list of output vectors
brt_rmse <- list() # root mean square error
rel_imp <- list() # relative influence of predictors
test_preds <- list() # predictions on the test set
pdp_values <- list() # partial dependence values from each BRT

for (i in 1:1000) {
  
  # first, let's randomize the entire dataframe
  set.seed(i)
  
  # now specify row number
  brt_random <- 
    rli_data %>% 
    mutate(ID = row_number())
  
  # randomly split the data into 80-20 training-test set
  brt_train <- 
    brt_random %>% 
    dplyr::sample_frac(0.8)
  
  # test set
  brt_test <- 
    anti_join(brt_random, brt_train, by = 'ID')
  
  # create model matrix
  modmat_train <- 
    stats::model.matrix(fmod, model.frame(~ ., brt_train, na.action = na.pass))[, -1]
  
  modmat_test <- 
    stats::model.matrix(fmod, model.frame(~ ., brt_test, na.action = na.pass))[, -1]
  
  # Run the model
  params <- 
    list(eta = 0.9, gamma = 0.1, max_depth = 10, subsample = 0.9)
  
  # The model
  brt_model <- 
    xgboost::xgboost(modmat_train,
                     objective = 'reg:logistic',
                     label = brt_train$RLI_2020_weighted,
                     params = params,
                     nrounds = 150, 
                     verbose = 0)
  
  # extract RMSE
  brt_rmse[[i]] <- 
    tibble(N_round = paste(i),
           rmse = min(brt_model$evaluation_log$train_rmse))
  
  # relative importance
  rel_imp[[i]] <- 
    as_tibble(xgb.importance(colnames(modmat_train), model = brt_model)) %>% 
    mutate(N_round = paste(i))
  
  # assess performance on test set
  test_preds[[i]] <- 
    brt_test %>% 
    # create column of predictions
    mutate(pred_RLI = stats::predict(brt_model, newdata = modmat_test),
           # calculate bias
           RLI_bias = RLI_2020_weighted - pred_RLI,
           N_round = paste(i)) %>% 
    dplyr::select(N_round, location, RLI_2020_weighted, pred_RLI, RLI_bias)
  
  # calculate partial dependence plot values
  pdp_dfs <- 
    lapply(colnames(modmat_train), function(variable_name) {
      
      df <- 
        as_tibble(pdp::partial(brt_model, 
                               pred.var = paste0(variable_name),
                               train = modmat_train)) %>% 
        mutate(variable = paste(variable_name)) %>% 
        dplyr::rename('value' = paste(variable_name))
      
      return(df)
      
    })
  
  pdp_values[[i]] <- 
    do.call(rbind, pdp_dfs) %>% 
    mutate(N_round = paste(i))
  
  cat(paste('Bootstrapping the model, round', i), '\n')
  
}

# Clean things from the models
# RMSE
rmse_clean <- 
  do.call(rbind, brt_rmse)

head(rmse_clean)

#write.csv(rmse_clean, '../Model_outputs/BRT_RMSE.csv', row.names = FALSE)

# Relative influence
ri_clean <- 
  do.call(rbind, rel_imp)

head(ri_clean)

#write.csv(ri_clean, '../Model_outputs/BRT_RelImp.csv', row.names = FALSE)

# Test set RLI bias
test_clean <- 
  do.call(rbind, test_preds)

head(test_clean)

#write.csv(test_clean, '../Model_outputs/BRT_TestSet.csv', row.names = FALSE)

# Partial dependence values
pdp_clean <- 
  do.call(rbind, pdp_values)

head(pdp_clean)

#write.csv(pdp_clean, '../Model_outputs/BRT_PDP.csv', row.names = FALSE)

# SUMMARY STATISTICS ----------------------------------------------------------------------
# Let's compare the RMSE and prediction bias of the models including HDI or WGI
test_clean %>% 
  group_by(N_round) %>% 
  summarise(bias = mean(RLI_bias)) %>% 
  ungroup() %>% 
  summarise(mean(bias),
            min(bias),
            max(bias))


# Relative influence
#read_csv('../Model_outputs/BRT_RelImp_WGI.csv') %>% 
  ri_clean %>% 
  group_by(Feature) %>% 
  summarise(mean = sum(Gain)/1000) %>% 
  ungroup() %>% 
  arrange(-mean)
