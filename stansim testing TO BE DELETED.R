library(rstansim)

library(readr); library(rstan); library(shinystan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# main loadings, identical across three factors

initf1 <- function() {
  
  loadingValues <- matrix(0.0001, 15, 3)
  ind <- 0
  
  for(i in 1:3){
    loadingValues[1 + ind,i] <- 0.8
    
    loadingValues[2 + ind,i] <- 0.7
    
    loadingValues[3 + ind,i] <- 0.6
    
    loadingValues[4 + ind,i] <- 0.5
    
    loadingValues[5 + ind,i] <- 0.4
    
    ind <- ind + 5
  }
  
  list(ML = t(loadingValues),
       CL = array(0, dim = c(3,15)))
} 

stan_arg_list <- list(file = 'sim_study/stan_models/small_guassian.stan', 
     iter = 5000, chains = 4, thin = 5, init = initf1)

data_list <- dir("sim_study/sim_data/dense/N_200/lat_no_corr/data_output", full.names = T)[1:4]

testout <- stan_sim(stan_args = stan_arg_list, sim_data = data_list, loo = FALSE, use_cores = 4, parameters = c("ML", "CL"))

library(dplyr)
temp <- testout$data

temp[grep("ML", temp$parameter),] %>% filter(estimate =="50%") %>% 
  ggplot(aes(x = value, y = parameter)) + geom_point()

