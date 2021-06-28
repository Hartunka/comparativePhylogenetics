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
}

transformed data{
    int xMat[2*N, 3] = rep_array( 0, 2*N, 3);
    matrix[N, N] L[nMat];
    for (i in 1:nMat){
      L[i] = cholesky_decompose(C[i]);
    }
    
    for (i in 1:(2*N)){
      xMat[i, x[i]] = 1;
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
  vector[3] a;
  // slope for mn output
  vector[3] b;
}

transformed parameters{
  vector[2 * N] mu = dMat * z;
}

model{
  
  //vector[3] lp;
  vector[3] p;
  z ~ normal( 2, 0.5 );
  a ~ normal( 0, 0.5 );
  b ~ normal( 1, 0.5 );
  R_L ~ lkj_corr_cholesky( eta );
  
  for ( j in 1:nMat ){
    matrix[2 * N,2 * N] V = kronecker(R_L, L[j]);
    brown[j] ~ multi_normal_cholesky(mu, V );
    for ( i in 1:(2*N) ){
      
      p[1] = exp( a[1] + b[1] * brown[j][i] );
      p[2] = exp( a[2] + b[2] * brown[j][i] );
      p[3] = exp( a[3] + b[3] * brown[j][i] );
      
      //p = exp( lp );
      xMat[i,1] ~ poisson( p[1] );
      xMat[i,2] ~ poisson( p[2] );
      xMat[i,3] ~ poisson( p[3] );
    }
  }

}

generated quantities{
  real log_lik[nMat];
  corr_matrix[2] R = R_L * R_L';
  
  vector[3] p_gq;
  for ( j in 1:nMat ){
    real log_lik_j[2*N];
    for ( i in 1:(2*N) ){
      
      p_gq[1] = exp( a[1] + b[1] * brown[j][i] );
      p_gq[2] = exp( a[2] + b[2] * brown[j][i] );
      p_gq[3] = exp( a[3] + b[3] * brown[j][i] );
      
      // p_gq = softmax( lp_gq );
      // TODO:
      log_lik_j[i] = log_sum_exp( 
                         log_sum_exp( poisson_lpmf( xMat[i,1] |  p_gq[1] ),
                                      poisson_lpmf( xMat[i,2] |  p_gq[2] ) ),
                         poisson_lpmf( xMat[i,3] |  p_gq[3] )
                                );
    }
    log_lik[j] = log_sum_exp(log_lik_j);
  }
  
  
}

