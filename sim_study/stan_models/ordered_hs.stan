data{
  int N; // sample size
  int P; // number of observed variables
  int D; // number of latent variables
  int C; // number of cross loadings
  matrix[N, P] X; // Observed data matrix
}

parameters{
  vector[D] Ltv_score_vect[N]; // vectorised latent variable scores, matrix [N, D]
  corr_matrix[D] Omega; // latent variable correlation matrix
  matrix<lower=0>[D, P] ML; // Main Loading Matrix
  vector<lower=0, upper=pi()/2>[P] var_p_unif; // for efficient cauchy dist
  vector[C] z; // Gaussian used for matt tricking the horseshoe
  vector<lower=0, upper=pi()/2>[C] lambda_unif; // local hs var param
	real<lower=0, upper=pi()/2> tau_unif; // global hs var param
	vector[P] cl_zeros; // zero restraints for transformed CL
}

transformed parameters{
  matrix[N, P] mu; // NxP Matrix where each cell is ind N pred score for var P
  matrix[D,D] Ld; // scale-parameter for multi-normal latent variable scores
  vector<lower=0>[P] var_p; // sd for each variable, vector
  matrix[N, D] Ltv_score; // latent variable scores, matrix [N,D]
  vector<lower=0>[D] ordering; // positive values to enforce ordering
  vector[C] lambda; // local scale param for horseshoe priors
	real tau; // global scale param for horseshoe priors
  vector[C] hs_priors; // horseshoe priors for cross-loadings
  matrix[D, P] CL; // Cross Loading Matrix
  
  // transform Ltv_score_vect in to matrix format for use in matrix algebra
  for(i in 1:N) Ltv_score[i,] = to_row_vector(Ltv_score_vect[i]);
  
  // Omega is the latent variable correlation matrix, here decomposed
  // to its cholesky and multiplied with a D-dimensional 
  // identity matrix, fixing latent variable variance to 1. Ld is the
  // scale parameter for the multi-normal cholesky sampling
  // of latent variable scores.
  Ld = diag_matrix(rep_vector(1, D)) * cholesky_decompose(Omega);
  
  // efficient statement of var_p ~ cauchy(0,1)
  for(i in 1:P) var_p[i] = tan(var_p_unif[i]);
  
  // Set up horseshoe priors by setting up both global
  // and local Cauchy distribution parameters. These are
  // then multipled by z ~ normal(0, 1) to get an efficiently
  // sampled vector of hs priors.
  for(i in 1:C) lambda[i] = tan(lambda_unif[i]);  // lambda ~ cauchy(0,1)
	tau = tan(tau_unif);  // tau ~ cauchy(0,1)
	hs_priors = z .* lambda * tau; // cross loadings horseshoe prior (Matt trick)
	
	//////// specify cross loading matrix free estimates ///////
  // equivalent to CL[r, c] ~ horseshoe(0)
  CL[2, 1] = hs_priors[1]; CL[3, 1] = hs_priors[2];
  CL[2, 2] = hs_priors[3]; CL[3, 2] = hs_priors[4];
  CL[2, 3] = hs_priors[5]; CL[3, 3] = hs_priors[6];
  CL[2, 4] = hs_priors[7]; CL[3, 4] = hs_priors[8];
  CL[2, 5] = hs_priors[9]; CL[3, 5] = hs_priors[10];
  CL[1, 6] = hs_priors[11]; CL[3, 6] = hs_priors[12];
  CL[1, 7] = hs_priors[13]; CL[3, 7] = hs_priors[14];
  CL[1, 8] = hs_priors[15]; CL[3, 8] = hs_priors[16];
  CL[1, 9] = hs_priors[17]; CL[3, 9] = hs_priors[18];
  CL[1, 10] = hs_priors[19]; CL[3, 10] = hs_priors[20];
  CL[1, 11] = hs_priors[21]; CL[2, 11] = hs_priors[22];
  CL[1, 12] = hs_priors[23]; CL[2, 12] = hs_priors[24];
  CL[1, 13] = hs_priors[25]; CL[2, 13] = hs_priors[26];
  CL[1, 14] = hs_priors[27]; CL[2, 14] = hs_priors[28];
  CL[1, 15] = hs_priors[29]; CL[2, 15] = hs_priors[30];
  
  ///////// specify cross loading matrix restrictions to 0 ///////
  // equivalent to CL[r, c] ~ normal(0, .001)
  CL[1, 1] = cl_zeros[1] * .001;
  CL[1, 2] = cl_zeros[2] * .001;
  CL[1, 3] = cl_zeros[3] * .001;
  CL[1, 4] = cl_zeros[4] * .001;
  CL[1, 5] = cl_zeros[5] * .001;
  CL[2, 6] = cl_zeros[6] * .001;
  CL[2, 7] = cl_zeros[7] * .001;
  CL[2, 8] = cl_zeros[8] * .001;
  CL[2, 9] = cl_zeros[9] * .001;
  CL[2, 10] = cl_zeros[10] * .001;
  CL[3, 11] = cl_zeros[11] * .001;
  CL[3, 12] = cl_zeros[12] * .001;
  CL[3, 13] = cl_zeros[13] * .001;
  CL[3, 14] = cl_zeros[14] * .001;
  CL[3, 15] = cl_zeros[15] * .001;

  // Factor Model Equation. ML and CL split in to 2 matrices
  // due to different constraints (ML positive, CL no constraints).
  // Equivalent to Ltv_score * full loading matrix. mu is 
  // predicted scores for all individuals on all observed vars.
  mu = (Ltv_score * ML) + (Ltv_score * CL);
  
  // values of ordering are restricted positive. The effect of 
  // the below is to fix the first main loading of each factor
  // to be greater than the last.
  ordering[1] = ML[1, 1] - ML[1, 5];
  ordering[2] = ML[2, 6] - ML[2, 10];
  ordering[3] = ML[3, 11] - ML[3, 15];
}

model{
  ///////// specify main loading matrix free estimates ///////
    // equivalent to ML[r, c] ~ normal(1, 1)
  ML[1, 1] - 1 ~ normal(0, 1);
  ML[1, 2] - 1 ~ normal(0, 1);
  ML[1, 3] - 1 ~ normal(0, 1);
  ML[1, 4] - 1 ~ normal(0, 1);
  ML[1, 5] - 1 ~ normal(0, 1);
  ML[2, 6] - 1 ~ normal(0, 1);
  ML[2, 7] - 1 ~ normal(0, 1);
  ML[2, 8] - 1 ~ normal(0, 1);
  ML[2, 9] - 1 ~ normal(0, 1);
  ML[2, 10] - 1 ~ normal(0, 1);
  ML[3, 11] - 1 ~ normal(0, 1);
  ML[3, 12] - 1 ~ normal(0, 1);
  ML[3, 13] - 1 ~ normal(0, 1);
  ML[3, 14] - 1 ~ normal(0, 1);
  ML[3, 15] - 1 ~ normal(0, 1);
  
  ///////// specify main loading matrix restrictions to 0 ///////
    // equivalent to ML[r, c] ~ normal(0, .001)
  ML[2, 1]*1000 ~ normal(0, 1); ML[3, 1]*1000 ~ normal(0, 1);
  ML[2, 2]*1000 ~ normal(0, 1); ML[3, 2]*1000 ~ normal(0, 1);
  ML[2, 3]*1000 ~ normal(0, 1); ML[3, 3]*1000 ~ normal(0, 1);
  ML[2, 4]*1000 ~ normal(0, 1); ML[3, 4]*1000 ~ normal(0, 1);
  ML[2, 5]*1000 ~ normal(0, 1); ML[3, 5]*1000 ~ normal(0, 1);
  ML[1, 6]*1000 ~ normal(0, 1); ML[3, 6]*1000 ~ normal(0, 1);
  ML[1, 7]*1000 ~ normal(0, 1); ML[3, 7]*1000 ~ normal(0, 1);
  ML[1, 8]*1000 ~ normal(0, 1); ML[3, 8]*1000 ~ normal(0, 1);
  ML[1, 9]*1000 ~ normal(0, 1); ML[3, 9]*1000 ~ normal(0, 1);
  ML[1, 10]*1000 ~ normal(0, 1); ML[3, 10]*1000 ~ normal(0, 1);
  ML[1, 11]*1000 ~ normal(0, 1); ML[2, 11]*1000 ~ normal(0, 1);
  ML[1, 12]*1000 ~ normal(0, 1); ML[2, 12]*1000 ~ normal(0, 1);
  ML[1, 13]*1000 ~ normal(0, 1); ML[2, 13]*1000 ~ normal(0, 1);
  ML[1, 14]*1000 ~ normal(0, 1); ML[2, 14]*1000 ~ normal(0, 1);
  ML[1, 15]*1000 ~ normal(0, 1); ML[2, 15]*1000 ~ normal(0, 1);
  
  /////////////////////////////////////////////////////
  
  // prior for use in Horseshoe distribution construction
  z ~ normal(0,1);
  
  // zero values for fixing non-CL entries in CL matrix
  cl_zeros ~ normal(0,1);
    
  // for each individual, their actual observed item value
  // is normally distributed around the estimated mean with
  // cauchy+ distributed variances. Item level errors are uncorrelated
  // so vectorised normal() can be used rather than multi_normal()
  for(i in 1:N) X[i] ~ normal(mu[i], var_p);
  
  // latent variable scores have a multi_normal() distribution
  // with mean vector fixed at zero, and a covariance matrix with
  // variances fixed at one and correlations freely estimated.
  // multi_normal_cholesky() is used for sampling efficiency.
  Ltv_score_vect ~ multi_normal_cholesky(rep_vector(0, D), Ld);
}

// generate the log-likelihood variable to be used in LOO calculations
generated quantities {
  vector[P] log_lik[N];  // log_likelihood matrix of order [N,P]
  
  for(n in 1:P){
    for (i in 1:N){ 
      log_lik[i, n] = normal_lpdf(X[i, n] | mu[i, n], var_p[n]); 	
    }
  }
}
