functions {
  vector X_model(real time, vector y, array[] real params) {
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
    E1_to_E2 = 3*0.5*y[2];
    I1_to_I2 = 3*0.5*y[3];
    I2_to_I3 = 3*0.5*y[6];
    I3_to_R = 3*0.5*y[7];
    E2_to_E3 = 3*0.5*y[8];
    E3_to_I1 = 3*0.5*y[9];
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
  array[n_obs] real ts;
}
parameters {
  real<lower = 0> par_beta;
  real<lower = 0, upper = 1> par_rho;
  real<lower = 0> I0;
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
  x0[6] = 0; // I2
  x0[7] = 0; // I3
  x0[8] = 0; // E2
  x0[9] = 0; // E3
  params[1] = par_beta;
  params[2] = par_rho;
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
  y ~ poisson(delta_x_1);
}
generated quantities {
  real log_lik;
  log_lik = poisson_lpmf(y | delta_x_1);
}

