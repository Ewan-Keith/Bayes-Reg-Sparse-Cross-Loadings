library(readr); library(rstan); library(shinystan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

######### testing 200 ind no cor data ############
# WORKS

ind_nocor_data <- read_csv(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/ind/N_200/lat_no_corr/data_output/ind_200_no-corr_file29.csv")

# prep data for stan
ind_nocor_dat <- read_rds(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/ind/N_200/lat_no_corr/data_output/ind_200_no-corr_file99.rds")

#list(N = 200,
#                  P = 15,
#                  D = 3,
#                  C = 30,
#                  X = as.matrix(ind_nocor_data)
# )

# fit model
fit <- stan(file = 'sim_study/stan_models/small_hs.stan', data = ind_nocor_dat, 
            iter = 5000, chains = 4, thin = 5, init = initf1, control = list(adapt_delta = 0.9))

###### testing ind cor data ###########
# WORKS

ind_cor_data <- read_csv(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/ind/N_200/lat_corr/data_output/ind_200_corr_file14.csv")

# prep data for stan
ind_cor_dat <- list(N = 200,
                      P = 15,
                      D = 3,
                      C = 30,
                      X = as.matrix(ind_cor_data)
)

# fit model
fit <- stan(file = 'sim_study/stan_models/small_guassian.stan', data = ind_cor_dat, 
            iter = 5000, chains = 4, thin = 5, init = initf1)

###### testing sparse no_cor data ###########
# WORKS

sp_nocor_data <- read_csv(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/sparse/N_200/lat_no_corr/data_output/sparse_200_no-corr_file14.csv")

# prep data for stan
sp_nocor_dat <- list(N = 200,
                    P = 15,
                    D = 3,
                    C = 30,
                    X = as.matrix(sp_nocor_data)
)

# fit model
fit <- stan(file = 'sim_study/stan_models/gaussian_cl_prior_sem matrix mult.stan', data = sp_nocor_dat, 
            iter = 5000, chains = 4, thin = 5)

###### testing sparse cor data ###########
# WORKS

sp_cor_data <- read_rds(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/sparse/N_200/lat_corr/data_output/sparse_200_corr_file19.rds")

# prep data for stan
sp_cor_dat <- list(N = 200,
                     P = 15,
                     D = 3,
                     C = 30,
                     X = as.matrix(sp_cor_data)
)

# fit model
fit <- stan(file = 'sim_study/stan_models/small_hs.stan', data = sp_cor_data, 
            iter = 5000, chains = 4, thin = 5, init = initf1, control = list(adapt_delta = 0.9))

###### testing dense nocor data ###########
# WORKS

dense_nocor_data <- read_csv(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/dense/N_200/lat_no_corr/data_output/dense_200_no-corr_file42.csv")

# prep data for stan
dense_nocor_dat <- list(N = 200,
                   P = 15,
                   D = 3,
                   C = 30,
                   X = as.matrix(dense_nocor_data)
)

# fit model
fit <- stan(file = 'sim_study/stan_models/gaussian_cl_prior_sem matrix mult.stan', data = dense_nocor_dat, 
            iter = 5000, chains = 4, thin = 5)

###### testing dense cor data ###########
# WORKS

dense_cor_data <- read_csv(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/dense/N_200/lat_corr/data_output/dense_200_corr_file83.csv")

# prep data for stan
dense_cor_dat <- list(N = 200,
                        P = 15,
                        D = 3,
                        C = 30,
                        X = as.matrix(dense_cor_data)
)

# fit model
fit <- stan(file = 'sim_study/stan_models/gaussian_cl_prior_sem matrix mult.stan', data = dense_cor_dat, 
            iter = 5000, chains = 4, thin = 5)

###//////////////////////////////////////////////////////////////////////
#### LARGE SAMPLE TESTING ####
###//////////////////////////////////////////////////////////////////////

######### testing 1000 ind no cor data ############
# WORKS

ind_nocor_data_large <- read_csv(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/ind/N_1000/lat_no_corr/data_output/ind_1000_no-corr_file1.csv")

# prep data for stan
ind_nocor_dat_large <- list(N = 1000,
                      P = 15,
                      D = 3,
                      C = 30,
                      X = as.matrix(ind_nocor_data_large)
)

# fit model
fit <- stan(file = 'sim_study/stan_models/gaussian_cl_prior_sem matrix mult.stan', data = ind_nocor_dat_large, 
            iter = 5000, chains = 4, thin = 5, init = initf1)

##

test_data2 <- readRDS(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/sparse/N_1000/lat_corr/data_output/sparse_1000_corr_file32.rds")

test_data_hard <- readRDS(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/dense/N_1000/lat_corr/data_output/dense_1000_corr_file99.rds")

# prep data for stan
stan_dat <- list(N = 200,
                 P = 15,
                 D = 3,
                 C = 30,
                 X = as.matrix(test_data)
)

stan_dat2 <- list(N = 1000,
                 P = 15,
                 D = 3,
                 #C = 30,
                 X = as.matrix(test_data_hard)
)

# fit model
fit <- stan(file = 'sim_study/stan_models/gaussian_cl_prior_sem matrix mult.stan', data = stan_dat, 
            iter = 5000, chains = 4, thin = 5, control = list(adapt_delta = 0.99))

fit2 <- stan(file = 'sim_study/stan_models/large_hs.stan', data = test_data2, 
            iter = 5000, chains = 4, thin = 5, init = initf1, control = list(adapt_delta = 0.9))

library(shinystan)

launch_shinystan(fit)
launch_shinystan(fit2)

###### testing inits ######




initf1 <- function() {
  # main loadings, identical across three factors
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

###### testing from scratch model ###########
ind_nocor_data <- read_csv(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/ind/N_200/lat_no_corr/data_output/ind_200_no-corr_file1.csv")

scratch_dat <- list(N = 200,
                 P = 15,
                 X = as.matrix(ind_nocor_data)
)

scratch_fit <- stan(file = 'sim_study/stan_models/from_scratch.stan', data = scratch_dat, 
            iter = 4000, chains = 4, thin = 5)
