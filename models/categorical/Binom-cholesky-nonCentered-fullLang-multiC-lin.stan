// model for lineage-specific correlation assumption, fitting to all languages

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
  
  matrix designMatrix(int N){
    matrix[2*N, 2] dMat;
    for (i in 1:(N*2)) {
      for (j in 1:2) {
        if ( ((j-1*N) < i) && (i<=(j*N)) ){
          dMat[i,j] = 1;
        } else {
          dMat[i,j] = 0;
        }
      }
    }
    return dMat;
  }
}

data {
  // number of families : 34
  // total number of languages : 768
  // number of languages per family
  int<lower=0> N[34];
  int<lower=0> M;
  // traits, language-wise
  int x[2*768];
  // covariance / time matrix
  matrix[768, 768] C[M];
  // hyperparameters for prior:
  // lkj
  real<lower=0> eta;
  // sigma
  real<lower=0> lambda;
  // z priors
  real mu_z;
  real sigma_z;
}

transformed data {
  // positions for each family's 'submatrix' of C
  int ranges[34,4];
  ranges[1,1] = 1;
  ranges[1,3] = 1;
  ranges[1,2] = N[1];
  ranges[1,4] = (2*N[1]);
  for (f in 2:34){
    ranges[f,1] = ranges[f-1,2]+1;
    ranges[f,2] = ranges[f,1]+N[f]-1;
    ranges[f,3] = ranges[f-1,4]+1;
    ranges[f,4] = ranges[f,3]+(2*N[f])-1;
  }
}

parameters {
  // mean / assumed starting value per trait per family
  vector[2] z[34];
  // trait correlations cholesky-factored per family
  cholesky_factor_corr[2] R_L[34];
  // correlation scales per family
  vector<lower=0>[2] sigma[34];
  // not-center for brownian motion
  vector[2*768] not_center;
}

transformed parameters{
  vector[2*768] p[M];
  for (m in 1:M){
    for (f in 1:34){
      matrix[2,2] Sigma_L = diag_pre_multiply(sigma[f], R_L[f]);
      matrix[N[f], N[f]] C_f = C[m][ ranges[f,1]:ranges[f,2], ranges[f,1]:ranges[f,2] ];
      matrix[2*N[f], 2*N[f]] V = kronecker(Sigma_L, cholesky_decompose(C_f));
      matrix[2*N[f],2] dMat = designMatrix(N[f]);
      p[m][ranges[f,3]:ranges[f,4]] = (V * not_center[ranges[f,3]:ranges[f,4]]) + (dMat * z[f]);
    }
  }
}

model{
  
  for (f in 1:34){
    z[f] ~ normal( mu_z, sigma_z );
    R_L[f] ~ lkj_corr_cholesky( eta );
    sigma[f] ~ exponential( lambda );
  }
  not_center ~ std_normal();
  
  for ( m in 1:M ){
    x ~ binomial_logit( 1, p[m] );    
  }

}

generated quantities{
  // log likelihood
  vector[2*768] log_lik = rep_vector(0, 2*768);
  // trait correlations
  //corr_matrix[2] R[34];
  //for (f in 1:34){
  //  R[f] =  multiply_lower_tri_self_transpose(R_L[f]);
  //}
  { // nested, so ll_m isn't added as part of the output
    matrix[2*768,M] ll_m;
    for (n in 1:(2*768)){
      for( m in 1:M ){
        ll_m[n,m] = binomial_logit_lpmf( x[n] | 1, p[m][n] );
      }
      log_lik[n] = log(mean(exp(ll_m[n])));
    }
  }
}

