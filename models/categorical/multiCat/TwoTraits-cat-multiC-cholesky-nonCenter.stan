data {
  // number of languages
  int<lower=0> N;
  // traits, language-wise
  int x[N];
  // number of matrices
  int<lower=0> M;
  // covariance / time matrix
  matrix[N, N] C[M];
  // hyperparameters for prior:
  // sigma_sq
  real<lower=0> lambda;
  // z priors
  real mu_z;
  real sigma_z;
  // alpha priors
  real mu_alpha;
  real sigma_alpha;
  // beta priors
  //real mu_beta;
  //real sigma_beta;

}

transformed data{
    matrix[N, N] L[M];
    for (i in 1:M){
      L[i] = cholesky_decompose(C[i]);
    }
}

parameters {
  // mean / assumed starting value per trait
  real z;
  // transition rates
  real<lower=0> sigma_sq;
  // not-center for brownian motion
  vector[N] not_center;
  
  // 'intercept' for mn output
  vector[3] alpha;
  //vector[3] beta;
}

transformed parameters{
  // logit from decentered brownian motion
  vector[N] brown[M];
  // softmax for multinomial
  simplex[4] p[M,N];
  for (m in 1:M){
    // non-centered brownian motion
    brown[m] = ( (L[m]*sigma_sq) * not_center ) + z;
    
    // softmax construction
    for ( n in 1:N ){
      p[m][n] = softmax( 
        to_vector(
          { alpha[1] + brown[m][n],
            alpha[2] + brown[m][n],
            alpha[3] + brown[m][n],
            0 } // pivot
          )
        );
    }
  }
}

model{
  
  z ~ normal( mu_z, sigma_z );
  alpha ~ normal( mu_alpha, sigma_alpha );
  //beta ~ normal( mu_beta, sigma_beta );
  sigma_sq ~ exponential(lambda);
  not_center ~ normal(0, 1);
  
  for ( m in 1:M ){
    for ( n in 1:N ){
      x[n] ~ categorical( p[m,n] );
    }
  }

}

generated quantities{
  real log_lik[N];
    
  vector[4] soft_a = softmax( append_row(alpha, 0) );
  vector[4] soft_z = softmax( rep_vector(z, N) );
  real exp_z = exp( z );
  matrix[N,M] ll_m;
  for( m in 1:M ){
    for (n in 1:(N)){
      ll_m[n,m] = categorical_lpmf( x[n] | p[m,n] );
    }
  }
  for( n in 1:(N) ){
    log_lik[n] = log(mean(exp(ll_m[n])));
  }
  
  
}

