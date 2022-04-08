functions {
  vector SIR1(real time, vector y, real[] params) {
    vector[4] dydt;
    real S_to_I;
    real I1_to_R;
    real S_to_C;
    S_to_I = params[1]*y[1]*y[2]/1000;
    I1_to_R = 0.25*y[2];
    S_to_C = S_to_I;
    dydt[1] = -S_to_I;
    dydt[2] = S_to_I-I1_to_R;
    dydt[3] = I1_to_R;
    dydt[4] = S_to_C;
    return dydt;
  }
}
data {
  int<lower = 1> n_obs;
  int<lower = 1> n_params;
  int<lower = 1> n_difeq;
  int y[n_obs];
  real t0;
  real ts[n_obs + 1];
}
parameters {
  real<lower = 0> beta;
  real<lower = 0> I0;
}
transformed parameters{
  vector[n_difeq] o[n_obs + 1]; // Output from the ODE solver
  array[n_obs] real x;
  real pred;
  vector[n_difeq] x0;
  array[n_params] real params;
  x0[1] = 1000 - 1 * I0;
  x0[2] = I0;
  x0[3] = 0;
  x0[4] = 1 * I0;
  params[1] = beta;
  o = ode_rk45(SIR1, x0, t0, ts, params);
  x[1] =  o[1, 4]  - x0[4];
  for (i in 1:n_obs-1) {
    x[i + 1] = o[i + 1, 4] - o[i, 4] + 1e-5;
  }
  pred = o[n_obs + 1, 4] - o[n_obs, 4];
}
model {
  beta  ~ lognormal(0, 1);
  I0    ~ lognormal(0, 1);
  y     ~ poisson(x);
}
generated quantities {
  real log_lik;
  log_lik = poisson_lpmf(y | x);
}
