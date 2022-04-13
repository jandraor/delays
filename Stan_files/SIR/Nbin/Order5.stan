functions {
  vector SIR5(real time, vector y, real[] params) {
    vector[8] dydt;
    real S_to_I;
    real I2_to_I3;
    real S_to_C;
    real I1_to_I2;
    real I3_to_I4;
    real I4_to_I5;
    real I5_to_R;
    S_to_I = params[1]*y[1]*(y[2]+y[5]+y[6]+y[7]+y[8])/1000;
    I2_to_I3 = 1.25*y[5];
    S_to_C = S_to_I;
    I1_to_I2 = 1.25*y[2];
    I3_to_I4 = 1.25*y[6];
    I4_to_I5 = 1.25*y[7];
    I5_to_R = 1.25*y[8];
    dydt[1] = -S_to_I;
    dydt[2] = S_to_I-I1_to_I2;
    dydt[3] = I5_to_R;
    dydt[4] = S_to_C;
    dydt[5] = I1_to_I2-I2_to_I3;
    dydt[6] = I2_to_I3-I3_to_I4;
    dydt[7] = I3_to_I4-I4_to_I5;
    dydt[8] = I4_to_I5-I5_to_R;
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
  real<lower = 0> phi;
}
transformed parameters{
  vector[n_difeq] o[n_obs + 1]; // Output from the ODE solver
  array[n_obs] real x;
  real pred;
  vector[n_difeq] x0;
  array[n_params] real params;
  real phi_inv;
  phi_inv = 1 / phi;
  x0[1] = 1000 - 5 * I0;
  x0[2] = I0;
  x0[3] = 0;
  x0[4] = 5 * I0;
  x0[5] = I0;
  x0[6] = I0;
  x0[7] = I0;
  x0[8] = I0;
  params[1] = beta;
  o = ode_rk45(SIR5, x0, t0, ts, params);
  x[1] =  o[1, 4]  - x0[4];
  for (i in 1:n_obs-1) {
    x[i + 1] = o[i + 1, 4] - o[i, 4] + 1e-5;
  }
  pred = o[n_obs + 1, 4] - o[n_obs, 4];
}
model {
  beta  ~ lognormal(0, 1);
  I0    ~ lognormal(0, 1);
  phi   ~ exponential(6);
  y     ~ neg_binomial_2(x, phi_inv);
}
generated quantities {
  real log_lik;
  log_lik = neg_binomial_2_lpmf(y | x, phi_inv);
}
