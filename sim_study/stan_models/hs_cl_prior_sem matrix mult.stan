data{
  int N; // sample size
  int P; // number of variables
  int C; // number of cross loadings
  int D; // number of factors
  matrix[N, P] X; // data matrix of order [N,P]
}

parameters{
  vector[P] b; // intercepts
  matrix[N, D] FS; // factor scores, matrix of order [N,D]
  corr_matrix[D] Rho; // correlation matrix between factors
  matrix<lower=0>[D, P] ML; // Main Loading Matrix
  matrix[D, P] CL; // Cross Loading Matrix
  vector<lower=0, upper=pi()/2>[P] var_p_unif;
}

transformed parameters{
  vector[D] M;
  vector<lower=0, upper=1>[D] Sd_d; // sd of factors
  matrix[N, P] mu; // NxP Matrix where each cell is ind N mean score for var P
  matrix[D,D] Ld;
  matrix[N, P] b_dummy; // dummy intercept matrix
  vector<lower=0>[P] var_p; // sd for each variable, vector
  
  for (m in 1:D) {
    M[m] = 0;
    Sd_d[m] = 1;
  }
  
  Ld = diag_matrix(Sd_d) * cholesky_decompose(Rho);
  
  for(int_rows in 1:N){
    b_dummy[int_rows] = (rep_vector(1, P) .* b)';
  }
    
    for(i in 1:P) var_p[i] = tan(var_p_unif[i]);  // var_p ~ cauchy(0,1)
    
    mu = b_dummy + (FS * ML) + (FS * CL);
  }
    
    
    
    model{
    
    ///////// specify main loading matrix free estimates ///////
    
    ML[1, 1] ~ normal(.7, .05);
    ML[1, 2] ~ normal(0, 1);
    ML[1, 3] ~ normal(0, 1);
    ML[1, 4] ~ normal(0, 1);
    ML[1, 5] ~ normal(0, 1);
    ML[2, 6] ~ normal(0, 1);
    ML[2, 7] ~ normal(0, 1);
    ML[2, 8] ~ normal(0, 1);
    ML[2, 9] ~ normal(0, 1);
    ML[2, 10] ~ normal(0, 1);
    ML[3, 11] ~ normal(0, 1);
    ML[3, 12] ~ normal(0, 1);
    ML[3, 13] ~ normal(0, 1);
    ML[3, 14] ~ normal(0, 1);
    ML[3, 15] ~ normal(0, 1);
    
    ///////// specify main loading matrix restrictions to 0 ///////
    
    ML[2, 1]/.001 ~ normal(0, 1); ML[3, 1]/.001 ~ normal(0, 1);
    ML[2, 2]/.001 ~ normal(0, 1); ML[3, 2]/.001 ~ normal(0, 1);
    ML[2, 3]/.001 ~ normal(0, 1); ML[3, 3]/.001 ~ normal(0, 1);
    ML[2, 4]/.001 ~ normal(0, 1); ML[3, 4]/.001 ~ normal(0, 1);
    ML[2, 5]/.001 ~ normal(0, 1); ML[3, 5]/.001 ~ normal(0, 1);
    ML[1, 6]/.001 ~ normal(0, 1); ML[3, 6]/.001 ~ normal(0, 1);
    ML[1, 7]/.001 ~ normal(0, 1); ML[3, 7]/.001 ~ normal(0, 1);
    ML[1, 8]/.001 ~ normal(0, 1); ML[3, 8]/.001 ~ normal(0, 1);
    ML[1, 9]/.001 ~ normal(0, 1); ML[3, 9]/.001 ~ normal(0, 1);
    ML[1, 10]/.001 ~ normal(0, 1); ML[3, 10]/.001 ~ normal(0, 1);
    ML[1, 11]/.001 ~ normal(0, 1); ML[2, 11]/.001 ~ normal(0, 1);
    ML[1, 12]/.001 ~ normal(0, 1); ML[2, 12]/.001 ~ normal(0, 1);
    ML[1, 13]/.001 ~ normal(0, 1); ML[2, 13]/.001 ~ normal(0, 1);
    ML[1, 14]/.001 ~ normal(0, 1); ML[2, 14]/.001 ~ normal(0, 1);
    ML[1, 15]/.001 ~ normal(0, 1); ML[2, 15]/.001 ~ normal(0, 1);
    
    ///////// specify cross loading matrix restrictions to 0 ///////
    
    CL[1, 1]/.001 ~ normal(0, 1);
    CL[1, 2]/.001 ~ normal(0, 1);
    CL[1, 3]/.001 ~ normal(0, 1);
    CL[1, 4]/.001 ~ normal(0, 1);
    CL[1, 5]/.001 ~ normal(0, 1);
    CL[2, 6]/.001 ~ normal(0, 1);
    CL[2, 7]/.001 ~ normal(0, 1);
    CL[2, 8]/.001 ~ normal(0, 1);
    CL[2, 9]/.001 ~ normal(0, 1);
    CL[2, 10]/.001 ~ normal(0, 1);
    CL[3, 11]/.001 ~ normal(0, 1);
    CL[3, 12]/.001 ~ normal(0, 1);
    CL[3, 13]/.001 ~ normal(0, 1);
    CL[3, 14]/.001 ~ normal(0, 1);
    CL[3, 15]/.001 ~ normal(0, 1);
    
    ///////// specify cross loading matrix free estimates ///////
    
    CL[2, 1]/.1 ~ normal(0, 1); CL[3, 1]/.1 ~ normal(0, 1);
    CL[2, 2]/.1 ~ normal(0, 1); CL[3, 2]/.1 ~ normal(0, 1);
    CL[2, 3]/.1 ~ normal(0, 1); CL[3, 3]/.1 ~ normal(0, 1);
    CL[2, 4]/.1 ~ normal(0, 1); CL[3, 4]/.1 ~ normal(0, 1);
    CL[2, 5]/.1 ~ normal(0, 1); CL[3, 5]/.1 ~ normal(0, 1);
    CL[1, 6]/.1 ~ normal(0, 1); CL[3, 6]/.1 ~ normal(0, 1);
    CL[1, 7]/.1 ~ normal(0, 1); CL[3, 7]/.1 ~ normal(0, 1);
    CL[1, 8]/.1 ~ normal(0, 1); CL[3, 8]/.1 ~ normal(0, 1);
    CL[1, 9]/.1 ~ normal(0, 1); CL[3, 9]/.1 ~ normal(0, 1);
    CL[1, 10]/.1 ~ normal(0, 1); CL[3, 10]/.1 ~ normal(0, 1);
    CL[1, 11]/.1 ~ normal(0, 1); CL[2, 11]/.1 ~ normal(0, 1);
    CL[1, 12]/.1 ~ normal(0, 1); CL[2, 12]/.1 ~ normal(0, 1);
    CL[1, 13]/.1 ~ normal(0, 1); CL[2, 13]/.1 ~ normal(0, 1);
    CL[1, 14]/.1 ~ normal(0, 1); CL[2, 14]/.1 ~ normal(0, 1);
    CL[1, 15]/.1 ~ normal(0, 1); CL[2, 15]/.1 ~ normal(0, 1);
    
    /////////////////////////////////////////////////////
    
    
    b ~ normal(0, 1);
    //var_p ~ cauchy(0, 2);
    
    for(i in 1:N){   
    X[i] ~ normal(mu[i], var_p);
    FS[i] ~ multi_normal_cholesky(M, Ld);
    }
    
    }
    
    generated quantities {
    vector[P] log_lik[N];  // data matrix of order [N,P]
    
    for(n in 1:P){
    for (i in 1:N){ 
    log_lik[i, n] = normal_lpdf(X[i, n] | mu[i, n], var_p[n]); 	
    }
    }
    }
    