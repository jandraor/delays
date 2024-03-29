functions {
  vector X_model(real time, vector y, array[] real params) {
    vector[6] dydt;
    real E1_to_I1;
    real C_in;
    real aux_j;
    real aux_tau;
    real var_beta;
    real var_gamma;
    real I2_to_R;
    real S_to_E;
    real I1_to_I2;
    E1_to_I1 = params[5]*y[2];
    C_in = params[2]*E1_to_I1;
    aux_j = (2+1)/(2.0*2);
    aux_tau = params[4]-(1/params[5]);
    var_beta = (1/params[1])*(aux_j/aux_tau);
    var_gamma = var_beta*params[1];
    I2_to_R = 2*var_gamma*y[6];
    S_to_E = var_beta*y[1]*(y[3]+y[6])/params[3];
    I1_to_I2 = 2*var_gamma*y[3];
    dydt[1] = -S_to_E;
    dydt[2] = S_to_E-E1_to_I1;
    dydt[3] = E1_to_I1-I1_to_I2;
    dydt[4] = I2_to_R;
    dydt[5] = C_in;
    dydt[6] = I1_to_I2-I2_to_R;
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
  real N;
  real par_tau;
  real par_sigma;
  real xi;
}
parameters {
  real<lower = 0, upper = 1> par_inv_R0;
  real<lower = 0, upper = 1> par_rho;
  real<lower = 0> I0;
  real<lower = 0> inv_phi;
}
transformed parameters{
  array[n_obs] vector[n_difeq] x; // Output from the ODE solver
  array[n_params] real params;
  vector[n_difeq] x0; // init values
  array[n_obs] real delta_x_1;
  real phi;
  phi = 1 / inv_phi;
  x0[1] = N * (1 -  xi) - I0; // S
  x0[2] = 0; // E1
  x0[3] = I0; // I1
  x0[4] = xi * N; // R
  x0[5] = I0; // C
  x0[6] = 0; // I2
  params[1] = par_inv_R0;
  params[2] = par_rho;
  params[3] = N;
  params[4] = par_tau;
  params[5] = par_sigma;
  x = ode_rk45(X_model, x0, t0, ts, params);
  delta_x_1[1] =  x[1, 5] - x0[5] + 1e-5;
  for (i in 1:n_obs-1) {
    delta_x_1[i + 1] = x[i + 1, 5] - x[i, 5] + 1e-5;
  }
}
model {
  par_inv_R0 ~ beta(2, 2);
  par_rho ~ beta(2, 2);
  I0 ~ lognormal(0, 1);
  inv_phi ~ exponential(5);
  y ~ neg_binomial_2(delta_x_1, phi);
}
generated quantities {
  real log_lik;
  array[n_obs] int sim_y;
  log_lik = neg_binomial_2_lpmf(y | delta_x_1, phi);
  sim_y = neg_binomial_2_rng(delta_x_1, phi);
}

