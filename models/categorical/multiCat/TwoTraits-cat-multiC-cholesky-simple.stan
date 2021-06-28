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
  // transition rates
  cholesky_factor_corr[2] R_L;
  // brownian motion mn output
  vector[2 * N] brown[nMat];
}

transformed parameters{
  vector[2 * N] mu = dMat * z;
}

model{
  
  vector[3] s;
  vector[3] p;
  z ~ normal( mu_z, sigma_z );
  
  R_L ~ lkj_corr_cholesky( eta );
  
  for ( j in 1:nMat ){
    matrix[2 * N,2 * N] V = kronecker(R_L, L[j]);
    brown[j] ~ multi_normal_cholesky(mu, V );
    for ( i in 1:(2*N) ){
      
      s[1] = brown[j][i];
      s[2] = brown[j][i];
      s[3] = 0; // pivot
      
      p = softmax( s );
      target += categorical_lpmf( x[i] | p );
    }
  }

}

generated quantities{
  real log_lik[nMat];
  vector[2*N] x_pred[nMat];
  
  corr_matrix[2] R = R_L * R_L';
  
  for ( j in 1:nMat ){
    real log_lik_j[2*N];
    for ( i in 1:(2*N) ){
      vector[3] s_gq;
      s_gq[1] = brown[j][i];
      s_gq[2] = brown[j][i];
      s_gq[3] = 0; // pivot
      
      log_lik_j[i] = categorical_lpmf( x[i] | softmax( s_gq ) );
      x_pred[j][i] = categorical_rng( softmax( s_gq ) );
    }
    log_lik[j] = log_sum_exp(log_lik_j);
  }
  
  
}

