functions {
  vector X_model(real time, vector y, array[] real params) {
    vector[5] dydt;
    real E1_to_I1;
    real C_in;
    real aux_j;
    real aux_tau;
    real var_beta;
    real var_gamma;
    real S_to_E;
    real I_to_R;
    E1_to_I1 = 0.5*y[2];
    C_in = params[2]*E1_to_I1;
    aux_j = (1+1)/(2.0*1);
    aux_tau = params[3]-(1/0.5);
    var_beta = (1/params[1])*(aux_j/aux_tau);
    var_gamma = var_beta*params[1];
    S_to_E = var_beta*y[1]*y[3]/10000;
    I_to_R = 1*var_gamma*y[3];
    dydt[1] = -S_to_E;
    dydt[2] = S_to_E-E1_to_I1;
    dydt[3] = E1_to_I1-I_to_R;
    dydt[4] = I_to_R;
    dydt[5] = C_in;
    return dydt;
  }
}
data {
  int<lower = 1> n_obs;
  int<lower = 1> n_params;
  int<lower = 1> n_difeq;
  array[n_obs] real y;
  real t0;
  array[n_obs] real ts;
  real par_tau;
}
parameters {
  real<lower = 0, upper = 1> par_inv_R0;
  real<lower = 0, upper = 1> par_rho;
  real<lower = 0> I0;
  real<lower = 0> epsilon;
}
transformed parameters{
  array[n_obs] vector[n_difeq] x; // Output from the ODE solver
  array[n_params] real params;
  vector[n_difeq] x0; // init values
  array[n_obs] real delta_x_1;
  x0[1] = (10000) - I0; // S
  x0[2] = 0; // E1
  x0[3] = I0; // I1
  x0[4] = 0; // R
  x0[5] = I0; // C
  params[1] = par_inv_R0;
  params[2] = par_rho;
  params[3] = par_tau;
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
  epsilon ~ exponential(0.2);
  y ~ normal(delta_x_1, epsilon);
}
generated quantities {
  real log_lik;
  log_lik = normal_lpdf(y | delta_x_1, epsilon);
}

