create_SEIR_files <- function(m, n, inv_sigma, inv_gamma, N, meas_model) {
  
  orders <- cross2(m, n)
  
  lapply(orders, function(ord_obj) {
    
    E_ord <- ord_obj[[1]]
    I_ord <- ord_obj[[2]]
    
    output <- list(E_ord         = E_ord,
                   I_ord         = I_ord,
                   n_stocks      = NULL,
                   stan_filepath = NULL)
    
    fp              <- str_glue("./models/SEIR/SE{E_ord}I{I_ord}R.stmx")
    mdl             <- read_xmile(fp)
    consts          <- sd_constants(mdl)
    stocks          <- sd_stocks(mdl)
    C_index         <- which(stocks$name == "C")    
    output$n_stocks <- nrow(stocks)
    ODE_fn          <- str_glue("SE{E_ord}I{I_ord}R")
    
    sigma_val    <- 1 / (inv_sigma /  E_ord)
    gamma_val    <- 1 / (inv_gamma /  I_ord)
    
    stan_fun     <- stan_ode_function(fp, ODE_fn, pars = consts$name[c(1,5)],
                                      const_list = list(par_sigma = sigma_val,
                                                        par_gamma = gamma_val,
                                                        N         = N))
    
    fun_exe_line <- str_glue("  o = ode_rk45({ODE_fn}, x0, t0, ts, params);") 
    
    stan_data   <- paste(
      "data {",
      "  int<lower = 1> n_obs;",
      "  int<lower = 1> n_params;",
      "  int<lower = 1> n_difeq;",
      "  array[n_obs] int y;",
      "  real t0;",
      "  array[n_obs + 1] real ts;",
      "}", sep = "\n")
    
    stan_params <- get_stan_params(meas_model)
    
    stan_tp1a <- paste(
      "transformed parameters{",
      "  array[n_obs + 1] vector[n_difeq] o; // Output from the ODE solver",
      "  array[n_obs] real x;",
      "  real pred;",
      "  vector[n_difeq] x0;",
      "  array[n_params] real params;",
      sep = "\n")
    
    if(meas_model == "nbin") {
      
      stan_tp1a <- paste(stan_tp1a,
                         "  real phi_inv;",
                         "  phi_inv = 1 / phi;", sep = "\n")
    }
    
    stan_tp1b <- paste(
      str_glue(
        "  x0[1] = {N} - {I_ord} * I0;"),
        "  x0[2] = 0;",
        "  x0[3] = I0;",
        "  x0[4] = 0;",
        "  x0[5] = 0; // This is C", 
      sep = "\n")
    
    stan_tp1 <- paste(stan_tp1a, stan_tp1b, sep = "\n")
    
    if(I_ord > 1) {
      indexes <- 5 + 1:(I_ord - 1)
      additional_stocks <- str_glue("  x0[{indexes}] = I0;") |> 
        paste(collapse = "\n")
      
      stan_tp1 <- paste(stan_tp1, additional_stocks, sep = "\n")
    }
    
    if(E_ord > 1) {
      
      indexes <- 5 + I_ord - 1 + 1:(E_ord - 1)
      
      additional_stocks <- str_glue("  x0[{indexes}] = 0;") |> 
        paste(collapse = "\n")
      
      stan_tp1 <- paste(stan_tp1, additional_stocks, sep = "\n")
    }
    
    stan_tp2 <- paste(
      "  params[1] = beta;",
      "  params[2] = rho;",
      fun_exe_line,
      str_glue("  x[1] =  o[1, {C_index}]  - x0[{C_index}];"),
      "  for (i in 1:n_obs-1) {",
      str_glue("    x[i + 1] = o[i + 1, {C_index}] - o[i, {C_index}] + 1e-5;"),
      "  }",
      str_glue("  pred = o[n_obs + 1, {C_index}] - o[n_obs, {C_index}];"),
      "}", sep = "\n")
    
    stan_tp <- paste(stan_tp1, stan_tp2, sep = "\n")
    
    stan_model <- get_stan_model(meas_model)
    
    stan_gc <- get_stan_gc(meas_model)
    
    stan_text   <- paste(stan_fun, stan_data, stan_params,
                         stan_tp, stan_model, stan_gc, sep = "\n")
    
    dist  <- str_to_title(meas_model)    
    
    stan_fldr     <- str_glue("./Stan_files/SEIR/{dist}")
    dir.create(stan_fldr, showWarnings = FALSE, recursive = TRUE)  
    stan_filepath <- file.path(stan_fldr, str_glue("SE{E_ord}I{I_ord}R.stan"))
    
    output$stan_filepath <- stan_filepath
    
    create_stan_file(stan_text, stan_filepath)
    
    output
  })
}

get_stan_params <- function(meas_model) {
  
  
  if(meas_model == "pois") {
    
    sp <- paste(
      "parameters {",
      "  real<lower = 0>            beta;",
      "  real<lower = 0>            I0;",
      "  real<lower = 0, upper = 1> rho;",
      "}", sep = "\n")
  }
  
  if(meas_model == "nbin") {
    
    sp <- paste(
      "parameters {",
      "  real<lower = 0> beta;",
      "  real<lower = 0> I0;",
      "  real<lower = 0, upper = 1> rho;",
      "  real<lower = 0> phi;",
      "}", sep = "\n")
  }
  
  sp
}

get_stan_model <- function(meas_model) {
  
  if(meas_model == "pois") {
    
    sm <- paste(
      "model {",
      "  beta  ~ lognormal(0, 1);",
      "  rho   ~ beta(2, 2);",
      "  I0    ~ lognormal(0, 1);",
      "  y     ~ poisson(x);", 
      "}",
      sep = "\n")
  }
  
  if(meas_model == "nbin") {
    
    sm <- paste(
      "model {",
      "  beta  ~ lognormal(0, 1);",
      "  rho   ~ beta(2, 2);",
      "  I0    ~ lognormal(0, 1);",
      "  phi   ~ exponential(5);",
      "  y     ~ neg_binomial_2(x, phi_inv);",
      "}",
      sep = "\n")
  }
  
  sm
}

get_stan_gc <- function(meas_model) {
  
  if(meas_model == "pois") {
    
    sgc <- paste(
      "generated quantities {",
      "  real log_lik;",
      "  log_lik = poisson_lpmf(y | x);",
      "}",
      sep = "\n")
  }
  
  if(meas_model == "nbin") {
    
    sgc <- paste(
      "generated quantities {",
      "  real log_lik;",
      "  log_lik = neg_binomial_2_lpmf(y | x, phi_inv);",
      "}", sep = "\n")
  }
  
  sgc
}