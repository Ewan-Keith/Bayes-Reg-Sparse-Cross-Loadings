data{
	int N; // sample size
	int P; // number of observed variables
	// int C; // number of cross loadings
	// int D; // number of factors
	vector[P] X[N]; // data matrix of order [N,P]
}

parameters{
  vector<lower=0>[P] items_sd; // SD of observed variables
  vector[P] item_ints; // Intercepts for observed variables
}

transformed parameters{
  vector[P] mu[N]; // matrix of individual's item means
  matrix[P,P] L;
  matrix[P,P] chol_ident;
  
  // calculate mu
  // at the moment this is just item intercept (mean)
  for(i in 1:N){
    mu[i] = item_ints;
  }

  chol_ident = cholesky_decompose(diag_matrix(rep_vector(1, P)));

  // calcs input for observed data multi_normal_cholesky cov matrix
  // all correlations fixed at 0, only item SD is estimated
  L = diag_matrix(items_sd) * chol_ident;
  
}


model{
  for(c in 1:N){
  X[c] ~ multi_normal_cholesky(mu, L);
  }

   
  
  // priors
  item_ints ~ normal(0, 1);
  items_sd ~ normal(0, 1);
}
