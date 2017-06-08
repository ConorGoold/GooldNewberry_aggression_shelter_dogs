data {
  int<lower=1> N;     // length of data (long format)
  int<lower=1> Nid;   // number of individuals
  int<lower=1> K;     // number of predictors
  int y[N];           // response variable (0/1)
  int id[N];          // identification variable
  matrix[N,K] X;      // matrix of predictors
}

parameters{
  real a;              // group-level intercept
  vector[Nid] zaID;    // dog-varying intercepts' scaling vector (see the uncentered parameterisations)
  vector[K] Beta;      // vector of predictor regression coefficients
  real<lower=0> sigmaID;  // standard deviation of dog-varying intercepts
}

transformed parameters{
  vector[Nid] rID;        // transform scaled varying intercepts to the original scale
  rID = zaID * sigmaID;
}

model{
  vector[N] prob;         // probability of 1 (i.e. agonistic behaviour)

  // prior distributions
  sigmaID ~ cauchy( 0 , 2 );
  a ~ normal( 0 , 5 );
  Beta ~ normal( 0 , 1 );
  zaID ~ normal( 0 , 1 );

  // likelihood
  prob = a + X * Beta;       // the model (without varying intercepts)

  for(i in 1:N) {
    prob[i] = prob[i] + rID[id[i]];   // add varying intercepts
  }
  y ~ bernoulli_logit(prob);        // sampling statement
}

generated quantities {
  // calculate log-likelihood to compute WAIC
  vector[N] log_lik;

  for(i in 1:N) {
    log_lik[i] = bernoulli_logit_lpmf( y[i] | a + X[i] * Beta + rID[id[i]]);
  }
}
