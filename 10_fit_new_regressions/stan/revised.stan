data{
    int<lower=1> N;
    int<lower=1> N_nga;
    int MG[N];
    real Phylogeny[N];
    real Space[N];
    real Mean_c[N];
    int nga[N];
}
parameters{
    real a;
    vector[N_nga] a_nga;
    real<lower=0> a_sigma;
    real b_sc;
    real b_ph;
    real b_sp;
}
model{
    vector[N] p;
    b_sp ~ normal(0, 4);
    b_ph ~ normal(0, 4);
    b_sc ~ normal(0, 4);
    a_sigma ~ exponential(1);
    a_nga ~ normal(0, a_sigma);
    a ~ normal(0, 1);
    for (i in 1:N) {
        p[i] = a + a_nga[nga[i]] + b_sc * Mean_c[i] +
          b_ph * Phylogeny[i] + b_sp * Space[i];
    }
    MG ~ binomial_logit(1, p);
}
generated quantities{
    vector[N] p;
    for (i in 1:N) {
        p[i] = a + a_nga[nga[i]] + b_sc * Mean_c[i] +
          b_ph * Phylogeny[i] + b_sp * Space[i];
    }
}
