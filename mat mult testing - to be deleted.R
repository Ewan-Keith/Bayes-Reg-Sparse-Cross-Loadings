library(readr); library(rstan); library(shinystan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

test_data <- read_csv(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/ind/N_200/lat_no_corr/data_output/ind_200_no-corr_file1.csv")

test_data2 <- read_csv(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/ind/N_200/lat_corr/data_output/ind_200_corr_file32.csv")

test_data_hard <- read_csv(
  "~/Documents/github/Bayes Reg Sparse Cross Loadings/sim_study/sim_data/dense/N_1000/lat_corr/data_output/dense_1000_corr_file24.csv")

# stan_dat_mat <- list(N = 200,
#                  P = 15,
#                  D = 3,
#                  C = 30,
#                  X = as.matrix(test_data),
#                  ML_row = c(rep(1, 5), rep(2, 5), rep(3, 5)),
#                  ML_col = c(seq(1, 15)),
#                  CL_row = c(rep(1, 10), rep(2, 10), rep(3, 10)),
#                  CL_col = c(seq(6, 15), seq(1, 5), seq(11, 15), seq(1, 10))               
# )

stan_dat_mat <- list(N = 200,
                     P = 15,
                     D = 3,
                     C = 30,
                     X = as.matrix(test_data)     
)


fitmat <- stan(file = 'sim_study/stan_models/gaussian_cl_prior_sem matrix mult.stan', data = stan_dat_mat, 
            iter = 20, chains = 4, thin = 5, control = list(adapt_delta = 0.99))

launch_shinystan(fit)