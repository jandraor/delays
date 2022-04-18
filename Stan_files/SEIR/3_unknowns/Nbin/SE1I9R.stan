functions {
  vector SE1I9R(real time, vector y, real[] params) {
    vector[13] dydt;
    real S_to_E;
    real E1_to_I1;
    real I1_to_I2;
    real C_in;
    real I2_to_I3;
    real I3_to_I4;
    real I6_to_I7;
    real I4_to_I5;
    real I5_to_I6;
    real I7_to_I8;
    real I8_to_I9;
    real I9_to_R;
    S_to_E = params[1]*y[1]*(y[3]+y[6]+y[7]+y[8]+y[9]+y[10]+y[11]+y[12]+y[13])/1000;
    E1_to_I1 = 0.5*y[2];
    I1_to_I2 = 4.5*y[3];
    C_in = params[2]*E1_to_I1;
    I2_to_I3 = 4.5*y[6];
    I3_to_I4 = 4.5*y[7];
    I6_to_I7 = 4.5*y[9];
    I4_to_I5 = 4.5*y[8];
    I5_to_I6 = 4.5*y[9];
    I7_to_I8 = 4.5*y[11];
    I8_to_I9 = 4.5*y[12];
    I9_to_R = 4.5*y[13];
    dydt[1] = -S_to_E;
    dydt[2] = S_to_E-E1_to_I1;
    dydt[3] = E1_to_I1-I1_to_I2;
    dydt[4] = I9_to_R;
    dydt[5] = C_in;
    dydt[6] = I1_to_I2-I2_to_I3;
    dydt[7] = I2_to_I3-I3_to_I4;
    dydt[8] = I3_to_I4-I4_to_I5;
    dydt[9] = I4_to_I5-I5_to_I6;
    dydt[10] = I5_to_I6-I6_to_I7;
    dydt[11] = I6_to_I7-I7_to_I8;
    dydt[12] = I7_to_I8-I8_to_I9;
    dydt[13] = I8_to_I9-I9_to_R;
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
  x0[1] = 1000 - 9 * I0;
  x0[2] = 0;
  x0[3] = I0;
  x0[4] = 0;
  x0[5] = 0; // This is C
  x0[6] = I0;
  x0[7] = I0;
  x0[8] = I0;
  x0[9] = I0;
  x0[10] = I0;
  x0[11] = I0;
  x0[12] = I0;
  x0[13] = I0;
  params[1] = beta;
  params[2] = rho;
  o = ode_rk45(SE1I9R, x0, t0, ts, params);
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
