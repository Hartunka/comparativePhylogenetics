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
  // a mean
  real mu_a;
  // a standart deviation 
  real sigma_a;
  // b mean
  real mu_b;
  // b standart deviation 
  real sigma_b;
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
  // 'intercept' for mn output
  vector[2] a;
  // slope for mn output
  vector[2] b;
}

transformed parameters{
  vector[2 * N] mu = dMat * z;
}

model{
  
  vector[3] s;
  vector[3] p;
  z ~ normal( 2, 1 );
  a ~ normal( mu_a, sigma_a );
  b ~ normal( mu_b, sigma_b );
  R_L ~ lkj_corr_cholesky( eta );
  
  for ( j in 1:nMat ){
    matrix[2 * N,2 * N] V = kronecker(R_L, L[j]);
    brown[j] ~ multi_normal_cholesky(mu, V );
    for ( i in 1:(2*N) ){
      
      s[1] = a[1] + b[1] * brown[j][i];
      s[2] = a[2] + b[2] * brown[j][i];
      s[3] = 0; // pivot
      
      p = softmax( s );
      target += categorical_lpmf( x[i] | p );
    }
  }

}

generated quantities{
  real log_lik[nMat];
  
  corr_matrix[2] R = R_L * R_L';
  
  vector[3] s_gq;
  vector[3] p_gq;
  for ( j in 1:nMat ){
    real log_lik_j[2*N];
    for ( i in 1:(2*N) ){
      
      s_gq[1] = a[1] + b[1] * brown[j][i];
      s_gq[2] = a[2] + b[2] * brown[j][i];
      s_gq[3] = 0; // pivot
      
      p_gq = softmax( s_gq );
      log_lik_j[i] = categorical_lpmf( x[i] | p_gq );
    }
    log_lik[j] = log_sum_exp(log_lik_j);
  }
  
  
}

