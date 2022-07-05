functions {
  vector SE3I3R(real time, vector y, real[] params) {
    vector[9] dydt;
    real S_to_E;
    real E1_to_E2;
    real I1_to_I2;
    real I2_to_I3;
    real I3_to_R;
    real E2_to_E3;
    real E3_to_I1;
    real C_in;
    S_to_E = params[1]*y[1]*(y[3]+y[6]+y[7])/10000;
    E1_to_E2 = 1.5*y[2];
    I1_to_I2 = 1.5*y[3];
    I2_to_I3 = 1.5*y[6];
    I3_to_R = 1.5*y[7];
    E2_to_E3 = 1.5*y[8];
    E3_to_I1 = 1.5*y[9];
    C_in = params[2]*E3_to_I1;
    dydt[1] = -S_to_E;
    dydt[2] = S_to_E-E1_to_E2;
    dydt[3] = E3_to_I1-I1_to_I2;
    dydt[4] = I3_to_R;
    dydt[5] = C_in;
    dydt[6] = I1_to_I2-I2_to_I3;
    dydt[7] = I2_to_I3-I3_to_R;
    dydt[8] = E1_to_E2-E2_to_E3;
    dydt[9] = E2_to_E3-E3_to_I1;
    return dydt;
  }
}
data {
  int<lower = 1> n_obs;
  int<lower = 1> n_params;
  int<lower = 1> n_difeq;
  array[n_obs] int y;
  real t0;
  array[n_obs + 1] real ts;
}
parameters {
  real<lower = 0> beta;
  real<lower = 0> I0;
  real<lower = 0, upper = 1> rho;
  real<lower = 0> phi;
}
transformed parameters{
  array[n_obs + 1] vector[n_difeq] o; // Output from the ODE solver
  array[n_obs] real x;
  real pred;
  vector[n_difeq] x0;
  array[n_params] real params;
  real phi_inv;
  phi_inv = 1 / phi;
  x0[1] = 10000 - 3 * I0;
  x0[2] = 0;
  x0[3] = I0;
  x0[4] = 0;
  x0[5] = 0; // This is C
  x0[6] = I0;
  x0[7] = I0;
  x0[8] = 0;
  x0[9] = 0;
  params[1] = beta;
  params[2] = rho;
  o = ode_rk45(SE3I3R, x0, t0, ts, params);
  x[1] =  o[1, 5]  - x0[5];
  for (i in 1:n_obs-1) {
    x[i + 1] = o[i + 1, 5] - o[i, 5] + 1e-5;
  }
  pred = o[n_obs + 1, 5] - o[n_obs, 5];
}
model {
  beta  ~ lognormal(0, 1);
  rho   ~ beta(2, 2);
  I0    ~ lognormal(0, 1);
  phi   ~ exponential(5);
  y     ~ neg_binomial_2(x, phi_inv);
}
generated quantities {
  real log_lik;
  log_lik = neg_binomial_2_lpmf(y | x, phi_inv);
}
