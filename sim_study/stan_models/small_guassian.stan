data{
  int N; // sample size
  int P; // number of observed variables
  int D; // number of latent variables
  matrix[N, P] X; // Observed data matrix
}

parameters{
  vector[D] Ltv_score_vect[N]; // vectorised latent variable scores, matrix [N, D]
  corr_matrix[D] Omega; // latent variable correlation matrix
  matrix<lower=0>[D, P] ML; // Main Loading Matrix
  matrix[D, P] CL; // Cross Loading Matrix
  vector<lower=0, upper=pi()/2>[P] var_p_unif; // for efficient cauchy dist
}

transformed parameters{
  matrix[N, P] mu; // NxP Matrix where each cell is ind N pred score for var P
  matrix[D,D] Ld; // scale-parameter for multi-normal latent variable scores
  vector<lower=0>[P] var_p; // sd for each variable, vector
  matrix[N, D] Ltv_score; // latent variable scores, matrix [N,D]
  //vector<lower=0>[D] ordering; // positive values to enforce ordering
  
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
    
  // Factor Model Equation. ML and CL split in to 2 matrices
  // due to different constraints (ML positive, CL no constraints).
  // Equivalent to Ltv_score * full loading matrix. mu is 
  // predicted scores for all individuals on all observed vars.
  mu = (Ltv_score * ML) + (Ltv_score * CL);
  
  // value sof ordering are restricted positive. The effect of 
  // the below is to fix the first main loading of each factor
  // to be greater than the last.
  //ordering[1] = ML[1, 1] - ML[1, 5];
  //ordering[2] = ML[2, 6] - ML[2, 10];
  //ordering[3] = ML[3, 11] - ML[3, 15];
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
    
    ///////// specify cross loading matrix restrictions to 0 ///////
    // equivalent to CL[r, c] ~ normal(0, .001)
    CL[1, 1]*1000 ~ normal(0, 1);
    CL[1, 2]*1000 ~ normal(0, 1);
    CL[1, 3]*1000 ~ normal(0, 1);
    CL[1, 4]*1000 ~ normal(0, 1);
    CL[1, 5]*1000 ~ normal(0, 1);
    CL[2, 6]*1000 ~ normal(0, 1);
    CL[2, 7]*1000 ~ normal(0, 1);
    CL[2, 8]*1000 ~ normal(0, 1);
    CL[2, 9]*1000 ~ normal(0, 1);
    CL[2, 10]*1000 ~ normal(0, 1);
    CL[3, 11]*1000 ~ normal(0, 1);
    CL[3, 12]*1000 ~ normal(0, 1);
    CL[3, 13]*1000 ~ normal(0, 1);
    CL[3, 14]*1000 ~ normal(0, 1);
    CL[3, 15]*1000 ~ normal(0, 1);
    
    ///////// specify cross loading matrix free estimates ///////
    // equivalent to CL[r, c] ~ normal(0, .05)
    CL[2, 1]*20 ~ normal(0, 1); CL[3, 1]*20 ~ normal(0, 1);
    CL[2, 2]*20 ~ normal(0, 1); CL[3, 2]*20 ~ normal(0, 1);
    CL[2, 3]*20 ~ normal(0, 1); CL[3, 3]*20 ~ normal(0, 1);
    CL[2, 4]*20 ~ normal(0, 1); CL[3, 4]*20 ~ normal(0, 1);
    CL[2, 5]*20 ~ normal(0, 1); CL[3, 5]*20 ~ normal(0, 1);
    CL[1, 6]*20 ~ normal(0, 1); CL[3, 6]*20 ~ normal(0, 1);
    CL[1, 7]*20 ~ normal(0, 1); CL[3, 7]*20 ~ normal(0, 1);
    CL[1, 8]*20 ~ normal(0, 1); CL[3, 8]*20 ~ normal(0, 1);
    CL[1, 9]*20 ~ normal(0, 1); CL[3, 9]*20 ~ normal(0, 1);
    CL[1, 10]*20 ~ normal(0, 1); CL[3, 10]*20 ~ normal(0, 1);
    CL[1, 11]*20 ~ normal(0, 1); CL[2, 11]*20 ~ normal(0, 1);
    CL[1, 12]*20 ~ normal(0, 1); CL[2, 12]*20 ~ normal(0, 1);
    CL[1, 13]*20 ~ normal(0, 1); CL[2, 13]*20 ~ normal(0, 1);
    CL[1, 14]*20 ~ normal(0, 1); CL[2, 14]*20 ~ normal(0, 1);
    CL[1, 15]*20 ~ normal(0, 1); CL[2, 15]*20 ~ normal(0, 1);
    
    /////////////////////////////////////////////////////
    
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

//generated quantities {
//    vector[P] log_lik[N];  // data matrix of order [N,P]
//    
//    for(i in 1:N) log_lik[i] = normal_lpdf(X[i] | mu[i], var_p); 	
//}
    
//generated quantities {
//    vector[P] log_lik[N];  // data matrix of order [N,P]
//    
//    for(n in 1:P){
//      for (i in 1:N){ 
//        log_lik[i, n] = normal_lpdf(X[i, n] | mu[i, n], var_p[n]); 	
//      }
//    }
//}
