// early attempt to fit a separate correlation matrix for each phylogenetic covariance matrix

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
  // z priors
  real mu_z;
  real sigma_z;
}

transformed data{
    matrix[N, N] L[M];
    for (m in 1:M){
      L[m] = cholesky_decompose(C[m]);
    }
}

parameters {
  // mean / assumed starting value per trait
  vector[2] z;
  // trait correlations cholesky-factored
  cholesky_factor_corr[2] R_L[M];
  // not-center for brownian motion
  vector[2*N] not_center;
}

transformed parameters{
  // mean for brownian motion
  vector[2*N] mu = dMat * z;
  // logit from decentered brownian motion
  vector[2*N] p[M];
  for (m in 1:M){
    matrix[2*N, 2*N] V = kronecker(R_L[M], L[m]);
    p[m] = (V * not_center) + mu;
  }
}

model{
  
  z ~ normal( mu_z, sigma_z );
  not_center ~ normal(0, 1);
  
  for ( m in 1:M ){
    R_L[m] ~ lkj_corr_cholesky( eta );
    x ~ binomial_logit( 1, p[m] );    
  }

}

generated quantities{
  // log likelihood
  vector[2*N] log_lik;
  // trait correlations
  corr_matrix[2] R[M];

  matrix[2*N,M] ll_m;
  for( m in 1:M ){
    R[m] = R_L[m] * R_L[m]';
    for (n in 1:(2*N)){
      ll_m[n,m] = binomial_logit_lpmf( x[n] | 1, p[m][n] );
    }
  }
  for( n in 1:(2*N) ){
    log_lik[n] = log(mean(exp(ll_m[n])));
  }
}

