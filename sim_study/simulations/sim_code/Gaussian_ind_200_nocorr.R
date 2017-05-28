#### modelling for Using Bayesian regularisation to systematically estimate sparse cross-loading solutions
## Cross-Loading Structure: Independent
## samples size: 200
## Latent Correlations: Not Present
## CL prior Distribution: Gaussian

##//////////////////////////////////////////////////////////////////////////////
#### load required libraries ####
##//////////////////////////////////////////////////////////////////////////////

devtools::install_github("Ewan-Keith/rstansim")

library(rstan)
# STILL IN DEVELOPMENT, I/O not yet stable so code likely to break in future
library(rstansim)

##//////////////////////////////////////////////////////////////////////////////
#### set up simulation ####
##//////////////////////////////////////////////////////////////////////////////

# initial values function
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

stan_arg_list <- list(file = 'sim_study/stan_models/small_gaussian.stan', 
                      iter = 5000, chains = 4, thin = 5, init = initf1)

data_list <- dir("sim_study/sim_data/ind/N_200/lat_no_corr/data_output", full.names = T)

##//////////////////////////////////////////////////////////////////////////////
#### run simulation ####
##//////////////////////////////////////////////////////////////////////////////

stansim_out <- stan_sim(stan_args = stan_arg_list, 
                        sim_data = data_list, 
                        calc_loo = T, 
                        use_cores = 4, 
                        parameters = c("ML", "CL", "Omega"))

saveRDS(stansim_out, "sim_study/simulations/sim_results/Gaussian_ind_200_nocorr_output.rds")


