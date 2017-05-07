#### data simulation for Using Bayesian regularisation to systematically estimate sparse cross-loading solutions
## Cross-Loading Structure: Independent
## samples size: 200
## Latent Correlations: Not Present

library(simsem)

# //////////////////////////////////////////////////
#### specify sparse simsem model ####
# //////////////////////////////////////////////////

#### specify loading matrix ####

## Specify varying/fixed parameters

# main loadings
loading <- matrix(0, 15, 3)

loading[1:5, 1] <- NA 

loading[6:10, 2] <- NA 

loading[11:15, 3] <- NA 

## specify population values for variable parameters
loadingValues <- matrix(0, 15, 3)

# main loadings, identical across three factors
ind <- 0

for(i in 1:3){
  loadingValues[1 + ind,i] <- 0.8
  
  loadingValues[2 + ind,i] <- 0.7
  
  loadingValues[3 + ind,i] <- 0.6
  
  loadingValues[4 + ind,i] <- 0.5
  
  loadingValues[5 + ind,i] <- 0.4
  
  ind <- ind + 5
}

# bind into simsem object
LY <- bind(loading, loadingValues) 

#### define PS as the residual correlation matrix among factors ####

latentCor <- matrix(NA, 3, 3) 
diag(latentCor) <- 1 # population factor variance of 1
PS <- binds(latentCor,0) # population correlation of factors = 0

#### define TE  as the measurement error correlation among indicators ####

manifestCor <- matrix(NA, 15, 15) 
diag(manifestCor) <- 1 # population factor variance of 1
TE <- binds(manifestCor, 0) # residual variances of indicators = 1

#### create templates for data generation and analysis using model() ####
CFA.Model <- model(LY = LY, RPS = PS, RTE = TE, modelType = "CFA")

# //////////////////////////////////////////////////
#### Simulate Sparse Data ####
# //////////////////////////////////////////////////

n_data <- 100 # number of datasets to simulate

set.seed(456789)

for(i in 1:n_data){
  
  temp_data <- generate(CFA.Model, n = 200)
  
  out_list <- list(N = 200,
                   P = 15,
                   D = 3,
                   C = 30,
                   X = as.matrix(temp_data))
  
  saveRDS(out_list, 
          file = paste0("sim_study/sim_data/ind/N_200/lat_no_corr/data_output/ind_200_no-corr_file", i, ".rds"))
}
