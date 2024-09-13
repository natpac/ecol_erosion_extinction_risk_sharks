rm(list=ls()) # removes all objects from the environment
cat("\014") # clears the console

#### Load packages ####
if(!require("pacman")) install.packages("pacman")
p_load(tidyverse, here, rstan)
`%notin%` <- Negate(`%in%`) # define a new operator '%notin%' as the negation of '%in%'
####

### Data used obtained from the Sea Around Us Project version 50 (www.seaaroundus.org) 

#### Load data and setup ####
### Global catch
catch_obs <- read.csv(here("GlobalCatch.csv"))

### Global effort
effort_obs <- read.csv(here("GlobalEffort.csv"))
effort_obs <- effort_obs[1:61, ] # only have data 1950-2010 ; dropping NAs

### Years of observations
T1 = length(catch_obs$year)
T2 = length(effort_obs$year)

### How many years to estimate effort
T_pred = 9
####

#### Run model ####
### Combine simulated data and parameters
data_list <- list(
  T1 = length(catch_obs$year),
  T2 = length(effort_obs$year),
  T_pred = 9,
  y1 = catch_obs$total_catch,
  y2 = effort_obs$total_effort
)


### Fit the Stan model modified with large prior 
stan_fit_final <- stan_fit_final_TEST_largeprior <- stan(file = here("Stan_model.stan")  # Stan program
                                                         , data = data_list
                                                         , chains = 3
                                                         , thin = 3
                                                         , iter = 11000, warmup = 2000
)

print(stan_fit_final, probs = c(0.025,0.5, 0.975))
traceplot(stan_fit_final, pars = c("alpha1","alpha2"))
traceplot(stan_fit_final, pars = "mu1")
traceplot(stan_fit_final, pars = "mu2")
plot(stan_fit_final, show_density = TRUE, pars = "mu1", fill_color = "purple")
plot(stan_fit_final, show_density = TRUE, pars = "mu2", fill_color = "purple")
plot(stan_fit_final, show_density = TRUE, pars = "alpha1", fill_color = "purple")
plot(stan_fit_final, show_density = TRUE, pars = "alpha2", fill_color = "purple")



### FOR catch and effort
# Extract posterior samples
posterior_samples <- extract(stan_fit_final)

# Plot the results
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))


# Plot for the first observed dataset
posterior_means1 <- colMeans(posterior_samples$mu1)
posterior_intervals1 <- apply(posterior_samples$mu1, 2, function(x) quantile(x, c(0.025, 0.975)))

plot(1:T1, catch_obs$total_catch, type = 'l', col = 'blue', lwd = 2, ylab = 'Value', xlab = 'Time', main = 'Catch')
lines(1:T1, posterior_means1, col = 'red', lty = 2, lwd = 2)
shade1 <- cbind(posterior_intervals1[1, ], rev(posterior_intervals1[2, ]))
polygon(c(1:T1, rev(1:T1)), shade1, col = rgb(1, 0, 0, 0.2), border = NA)

# Plot for the second observed dataset
posterior_means2 <- colMeans(posterior_samples$mu2)
posterior_intervals2 <- apply(posterior_samples$mu2, 2, function(x) quantile(x, c(0.025, 0.975)))

plot(1:T1, c(effort_obs$total_effort, rep(NA,9)), type = 'l', col = 'blue', lwd = 2, ylab = 'Value', xlab = 'Time', main = 'Effort'
     ,ylim = c(0,max(posterior_intervals2)+.1*max(posterior_intervals2))
)
lines(1:T2, posterior_means2[1:T2], col = 'red', lty = 2, lwd = 2)
shade2 <- cbind(posterior_intervals2[1, 1:T2], rev(posterior_intervals2[2, 1:T2]))
polygon(c(1:T2, rev(1:T2)), shade2, col = rgb(1, 0, 0, 0.2), border = NA)

lines((T2+1):T1, posterior_means2[(T2+1):T1], col = 'darkgreen', lty = 2, lwd = 2)
shade2 <- cbind(posterior_intervals2[1, (T2+1):T1], rev(posterior_intervals2[2, (T2+1):T1]))
polygon(c((T2+1):T1, rev((T2+1):T1)), shade2, col = alpha('darkgreen', 0.2), border = NA)

# Plot for the CPUE
par(mfrow=c(1,1))
posterior_means_CPUE <- colMeans(posterior_samples$CPUE)
posterior_intervals_CPUE <- apply(posterior_samples$CPUE, 2, function(x) quantile(x, c(0.025, 0.975)))

plot(1:T1, (catch_obs$total_catch)/c(effort_obs$total_effort, rep(NA,9)), type = 'l', col = 'blue', lwd = 2, ylab = 'Value', xlab = 'Time', main = 'CPUE',
     ylim = c(0,0.17), las=1)
lines(1:T1, posterior_means_CPUE, col = 'red', lty = 2, lwd = 2)
shade1 <- cbind(posterior_intervals_CPUE[1, ], rev(posterior_intervals_CPUE[2, ]))
polygon(c(1:T1, rev(1:T1)), shade1, col = rgb(1, 0, 0, 0.2), border = NA)