create_SIR_files <- function(model_orders, inv_gamma_val, pop_val,
                             meas_model = "pois") {
  
  lapply(model_orders, function(mdl_ord) {
    
    output <- list(mdl_ord       = mdl_ord,
                   n_stocks      = NULL,
                   stan_filepath = NULL)
    
    fp           <- str_glue("./models/SIR/SIR{mdl_ord}.stmx") # filepath
    mdl          <- read_xmile(fp)
    consts       <- sd_constants(mdl)
    stocks       <- sd_stocks(mdl)
    output$n_stocks <- nrow(stocks)
    ODE_fn       <- str_glue("SIR{mdl_ord}")
    
    gamma_val    <- 1 / (inv_gamma_val /  mdl_ord)
    stan_fun     <- stan_ode_function(fp, ODE_fn, pars = consts$name[[1]],
                                      const_list = list(par_gamma = gamma_val,
                                                        N         = pop_val))
    
    fun_exe_line <- str_glue("  o = ode_rk45({ODE_fn}, x0, t0, ts, params);") 
    
    stan_data   <- paste(
      "data {",
      "  int<lower = 1> n_obs;",
      "  int<lower = 1> n_params;",
      "  int<lower = 1> n_difeq;",
      "  int y[n_obs];",
      "  real t0;",
      "  real ts[n_obs + 1];",
      "}", sep = "\n")
    
    stan_params <- get_stan_params(meas_model)
    
    stan_tp1a <- paste(
      "transformed parameters{",
      "  vector[n_difeq] o[n_obs + 1]; // Output from the ODE solver",
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
      str_glue("  x0[1] = {pop_val} - {mdl_ord} * I0;"),
      "  x0[2] = I0;",
      "  x0[3] = 0;",
      str_glue("  x0[4] = {mdl_ord} * I0;"),
      sep = "\n")
    
    stan_tp1 <- paste(stan_tp1a, stan_tp1b, sep = "\n")
    
    if(mdl_ord > 1) {
      indexes <- 4 + 1:(mdl_ord - 1)
      additional_stocks <- str_glue("  x0[{indexes}] = I0;") |> 
        paste(collapse = "\n")
      
      stan_tp1 <- paste(stan_tp1, additional_stocks, sep = "\n")
    }
    
    stan_tp2 <- paste(
      "  params[1] = beta;",
      fun_exe_line,
      "  x[1] =  o[1, 4]  - x0[4];",
      "  for (i in 1:n_obs-1) {",
      "    x[i + 1] = o[i + 1, 4] - o[i, 4] + 1e-5;",
      "  }",
      "  pred = o[n_obs + 1, 4] - o[n_obs, 4];",
      "}", sep = "\n")
    
    stan_tp <- paste(stan_tp1, stan_tp2, sep = "\n")
    
    stan_model <- get_stan_model(meas_model)
    
    stan_gc <- get_stan_gc(meas_model)
    
    stan_text   <- paste(stan_fun, stan_data, stan_params,
                         stan_tp, stan_model, stan_gc, sep = "\n")
    
    dist  <- str_to_title(meas_model)    
    
    stan_fldr     <- str_glue("./Stan_files/SIR/{dist}")
    dir.create(stan_fldr, showWarnings = FALSE, recursive = TRUE)  
    stan_filepath <- file.path(stan_fldr, str_glue("Order{mdl_ord}.stan"))
    
    output$stan_filepath <- stan_filepath
    
    create_stan_file(stan_text, stan_filepath)
    
    output
  })
}

get_stan_params <- function(meas_model) {
  
  
  if(meas_model == "pois") {
    
    sp <- paste(
      "parameters {",
      "  real<lower = 0> beta;",
      "  real<lower = 0> I0;",
      "}", sep = "\n")
  }
  
  if(meas_model == "nbin") {
    
    sp <- paste(
      "parameters {",
      "  real<lower = 0> beta;",
      "  real<lower = 0> I0;",
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
      "  I0    ~ lognormal(0, 1);",
      "  y     ~ poisson(x);", 
      "}",
      sep = "\n")
  }
  
  if(meas_model == "nbin") {
    
    sm <- paste(
      "model {",
      "  beta  ~ lognormal(0, 1);",
      "  I0    ~ lognormal(0, 1);",
      "  phi   ~ exponential(6);",
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