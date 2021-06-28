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
  int<lower=0> nMat;
  // covariance / time matrix
  matrix[N, N] C[nMat];
  //Design Matrix
  matrix[2 * N,2] dMat;
  // hyperparameters for prior:
  // lkj
  real<lower=0> eta;
  // z priors
  real mu_z;
  real sigma_z;
}

transformed data{
    matrix[N, N] L[nMat];
    for (i in 1:nMat){
      L[i] = cholesky_decompose(C[i]);
    }
}

parameters {
  // mean / assumed starting value per trait
  vector[2] z;
  // trait correlations cholesky-factored
  cholesky_factor_corr[2] R_L;
  // multi normal 'output'
  vector[2 * N] p[nMat];
}

transformed parameters{
  vector[2 * N] mu = dMat * z;
}

model{
  
  z ~ normal( mu_z, sigma_z );
  R_L ~ lkj_corr_cholesky( eta );
  
  for ( j in 1:nMat ){
    matrix[2 * N,2 * N] V = kronecker(R_L, L[j]);
    p[j] ~ multi_normal_cholesky(mu, V );
    
    target += binomial_logit_lpmf( x | 1, p[j] );    
  }

}

generated quantities{
  // log likelihood
  real log_lik[nMat];
  // posterior prediction
  int x_pred[nMat, 2*N];
  // trait correlations
  corr_matrix[2] R = R_L * R_L';

  for ( j in 1:nMat ){
    
    log_lik[j] = binomial_logit_lpmf( x | 1, p[j] );
    x_pred[j] = binomial_rng( 1, inv_logit(p[j]) );
  }
}

