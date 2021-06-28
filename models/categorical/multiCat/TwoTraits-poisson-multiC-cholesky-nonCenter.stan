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
  real mu_beta;
  real sigma_beta;
}

transformed data{
    matrix[N, N] L[M];
    int xMat[N, 4] = rep_array( 0, N, 4);
    for (m in 1:M){
      L[m] = cholesky_decompose(C[m]);
    }
    for ( n in 1:N){
      xMat[n, x[n]] = 1;
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
  vector[4] alpha;
  //vector[4] beta;
}

transformed parameters{
  // logit from decentered brownian motion
  vector[N] brown[M];
  // log-probabilities for poisson
  vector[4] l[M,N];
  for (m in 1:M){
    // non-centered brownian motion
    brown[m] = ( (L[m]*sigma_sq) * not_center ) + z ;
    
    // softmax construction
    for ( n in 1:N ){
      l[m][n] = to_vector(
        {alpha[1] + brown[m][n],
         alpha[2] + brown[m][n],
         alpha[3] + brown[m][n],
         alpha[4] + brown[m][n]}
        );
    }
  }
}

model{
  
  vector[3] s;
  
  z ~ normal( mu_z, sigma_z );
  alpha ~ normal( mu_alpha, sigma_alpha );
  //beta ~ normal( mu_beta, sigma_beta );
  sigma_sq ~ exponential(lambda);
  not_center ~ normal(0, 1);
  
  for ( m in 1:M ){
    for ( n in 1:N ){
      xMat[n,] ~ poisson( exp(l[m,n]) );
    }
  }

}

generated quantities{
  real log_lik[N];
    
  vector[4] p[M,N];
  vector[4] soft_a = softmax(alpha);
  
  matrix[N,M] ll_m;
  for( m in 1:M ){
    for (n in 1:(N)){
      ll_m[n,m] = poisson_lpmf( xMat[n,] | exp(l[m,n]) );
      p[m,n] = softmax(l[m,n]);
    }
  }
  for( n in 1:(N) ){
    log_lik[n] = log(mean(exp(ll_m[n])));
  }
  
  
}

