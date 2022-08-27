functions {
  vector X_model(real time, vector y, array[] real params) {
    vector[8] dydt;
    real S_to_E;
    real E1_to_I1;
    real I1_to_I2;
    real C_in;
    real I2_to_I3;
    real I3_to_I4;
    real I4_to_R;
    S_to_E = params[1]*y[1]*(y[3]+y[6]+y[7]+y[8])/10000;
    E1_to_I1 = 0.5*y[2];
    I1_to_I2 = 4*params[3]*y[3];
    C_in = params[2]*E1_to_I1;
    I2_to_I3 = 4*params[3]*y[6];
    I3_to_I4 = 4*params[3]*y[7];
    I4_to_R = 4*params[3]*y[8];
    dydt[1] = -S_to_E;
    dydt[2] = S_to_E-E1_to_I1;
    dydt[3] = E1_to_I1-I1_to_I2;
    dydt[4] = I4_to_R;
    dydt[5] = C_in;
    dydt[6] = I1_to_I2-I2_to_I3;
    dydt[7] = I2_to_I3-I3_to_I4;
    dydt[8] = I3_to_I4-I4_to_R;
    return dydt;
  }
}
data {
  int<lower = 1> n_obs;
  int<lower = 1> n_params;
  int<lower = 1> n_difeq;
  array[n_obs] int y;
  real t0;
  array[n_obs] real ts;
}
parameters {
  real<lower = 0> par_beta;
  real<lower = 0, upper = 1> par_rho;
  real<lower = 0> I0;
  real<lower = 0, upper = 1> par_gamma;
  real<lower = 0> inv_phi;
}
transformed parameters{
  array[n_obs] vector[n_difeq] x; // Output from the ODE solver
  array[n_params] real params;
  vector[n_difeq] x0; // init values
  array[n_obs] real delta_x_1;
  real phi;
  phi = 1 / inv_phi;
  x0[1] = (10000) - I0; // S
  x0[2] = 0; // E1
  x0[3] = I0; // I1
  x0[4] = 0; // R
  x0[5] = I0; // C
  x0[6] = 0; // I2
  x0[7] = 0; // I3
  x0[8] = 0; // I4
  params[1] = par_beta;
  params[2] = par_rho;
  params[3] = par_gamma;
  x = ode_rk45(X_model, x0, t0, ts, params);
  delta_x_1[1] =  x[1, 5] - x0[5] + 1e-5;
  for (i in 1:n_obs-1) {
    delta_x_1[i + 1] = x[i + 1, 5] - x[i, 5] + 1e-5;
  }
}
model {
  par_beta ~ lognormal(0, 1);
  par_rho ~ beta(2, 2);
  I0 ~ lognormal(0, 1);
  par_gamma ~ beta(2, 2);
  inv_phi ~ exponential(5);
  y ~ neg_binomial_2(delta_x_1, phi);
}
generated quantities {
  real log_lik;
  log_lik = neg_binomial_2_lpmf(y | delta_x_1, phi);
}
