// model used to fit models with several phylogenetic matrices (avgC-5, multiC-100)

functions{

  matrix kronecker(matrix R, matrix C) {
    
    int r = dims(R)[1];
    int n = dims(C)[1];
    
    matrix[n * r, n * r] V;
    
    for ( i in 1:r ){
        for (j in 1:r) {
            for (k in 1:n) {
                for (l in 1:n) {
                    V[ n * (i-1)+k , n * (j-1)+l ] = R[i, j] * C[k, l];
                }
            }
        }
    }
    return V;
  }

}

data {
  // number of languages
  int<lower=0> N;
  // traits, language-wise
  int x[2 * N];
  // number of matrices
  int<lower=0> M;
  // covariance / time matrix
  matrix[N, N] C[M];
  //Design Matrix
  matrix[2 * N,2] dMat;
  // hyperparameters for prior:
  // lkj
  real<lower=0> eta;
  // sigma
  real<lower=0> lambda;
  // z priors
  real mu_z;
  real sigma_z;
}

parameters {
  // mean / assumed starting value per trait
  vector[2] z;
  // trait correlations cholesky-factored
  cholesky_factor_corr[2] R_L;
  // correlation scales
  vector<lower=0>[2] sigma;
  // not-center for brownian motion
  vector[2*N] not_center;
}

transformed parameters{
  matrix[2,2] Sigma_L = diag_pre_multiply(sigma, R_L);
  vector[2*N] p[M];
  for (m in 1:M){
    matrix[2*N, 2*N] V = kronecker(Sigma_L, cholesky_decompose(C[m]));
    p[m] = (V * not_center) + (dMat * z);
  }
}

model{
  
  z ~ normal( mu_z, sigma_z );
  R_L ~ lkj_corr_cholesky( eta );
  sigma ~ exponential( lambda );
  not_center ~ normal(0, 1);
  
  for ( m in 1:M ){
    x ~ binomial_logit( 1, p[m] );    
  }

}

generated quantities{
  // log likelihood
  vector[2*N] log_lik;
  // correlation
  corr_matrix[2] R = multiply_lower_tri_self_transpose(R_L);
  cov_matrix[2] Sigma = multiply_lower_tri_self_transpose(Sigma_L);
  // trait correlations
  { // nested, so ll_m isn't added as part of the output
    matrix[2*N,M] ll_m;
    for( m in 1:M ){
      for (n in 1:(2*N)){
        ll_m[n,m] = binomial_logit_lpmf( x[n] | 1, p[m][n] );
      }
    }
    for( n in 1:(2*N) ){
      log_lik[n] = log(mean(exp(ll_m[n])));
    }
  }
}

