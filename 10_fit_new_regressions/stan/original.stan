// saved as ./temp/model1.stan
data{
  int<lower=1> N;
  int MG[N];
  real Mean_c[N];
  real Lag1[N];
  real Lag2[N];
  real Phylogeny[N];
  real Space[N];
}
parameters{
  real a;
  real b_sc;
  real b_l1;
  real b_l2;
  real b_ph;
  real b_sp;
}
model{
  vector[N] p;
  b_sp ~ normal( 0 , 7 );
  b_ph ~ normal( 0 , 7 );
  b_l2 ~ normal( 0 , 7 );
  b_l1 ~ normal( 0 , 7 );
  b_sc ~ normal( 0 , 7 );
  a ~ normal( 0 , 7 );
  for ( i in 1:N ) {
    p[i] = a + b_sc * Mean_c[i] + b_l1 * Lag1[i] + b_l2 * Lag2[i] + b_ph * Phylogeny[i] +    b_sp * Space[i];
  }
  MG ~ binomial_logit( 1 , p );
}
generated quantities{
  vector[N] p;
  for ( i in 1:N ) {
    p[i] = a + b_sc * Mean_c[i] + b_l1 * Lag1[i] + b_l2 * Lag2[i] + b_ph * Phylogeny[i] +    b_sp * Space[i];
  }
}

