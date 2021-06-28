// independence model used for avgC-1 / using only one phylogenetic covariance matrix

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
  // covariance / time matrix
  matrix[N, N] C;
  //Design Matrix
  matrix[2 * N,2] dMat;
  // hyperparameters for prior:
  // sigma
  real<lower=0> lambda;
  // z priors
  real mu_z;
  real sigma_z;
}

transformed data{
  // independence assuming fixed correlation
  corr_matrix[2] R_L = cholesky_decompose( diag_matrix(rep_vector(1, 2)) );
}

parameters {
  // mean / assumed starting value per trait
  vector[2] z;
  // correlation scales
  vector<lower=0>[2] sigma;
  // not-center for brownian motion
  vector[2*N] not_center;
}

transformed parameters{
  matrix[2,2] Sigma_L = diag_pre_multiply(sigma, R_L);
  vector[2*N] p;
  matrix[2*N, 2*N] V = kronecker( Sigma_L, cholesky_decompose(C));
  p = (V * not_center) + (dMat * z);
}

model{
  
  z ~ normal( mu_z, sigma_z );
  sigma ~ exponential( lambda );
  not_center ~ normal(0, 1);
  
  x ~ binomial_logit( 1, p );

}

generated quantities{
  // log likelihood
  vector[2*N] log_lik = rep_vector(0, 2*N);
  // correlation
  cov_matrix[2] Sigma = multiply_lower_tri_self_transpose(Sigma_L);
  // trait correlations
  for (n in 1:(2*N)){
    log_lik[n] = binomial_logit_lpmf( x[n] | 1, p[n] );
  }
}

