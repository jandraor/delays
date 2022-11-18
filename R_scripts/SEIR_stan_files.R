create_SEIR_files <- function(i, j, inv_sigma, inv_gamma, N, meas_model, 
                              stan_fldr, gamma_unk = FALSE,
                              alt_param = FALSE) {
  
  orders <- cross2(i, j)
  
  lapply(orders, SEIR_file, inv_sigma, inv_gamma, N, meas_model, stan_fldr,
         gamma_unk, alt_param)
}

SEIR_file <- function(ord_obj, inv_sigma, inv_gamma,N, meas_model,
                      stan_fldr, gamma_unk, alt_param) {
  
  E_ord <- ord_obj[[1]]
  I_ord <- ord_obj[[2]]
  
  sigma_val    <- 1 / inv_sigma
  gamma_val    <- 1 / inv_gamma
  
  output <- list(E_ord         = E_ord,
                 I_ord         = I_ord,
                 x0            = NULL,
                 n_stocks      = NULL,
                 stan_filepath = NULL)
  
  fp              <- str_glue("./models/SEIR/SE{E_ord}I{I_ord}R.stmx")
  
  if(alt_param) fp <- str_glue("./models/SEIR/SE{E_ord}I{I_ord}R_alt.stmx")
  
  stock_inits     <- get_stock_inits(E_ord, I_ord, N)
  
  mdl             <- read_xmile(fp, stock_list = stock_inits)
  stocks          <- sd_stocks(mdl)
  output$x0       <- stocks$init_value
  output$n_stocks <- nrow(stocks)
  
  const_list <- list(par_sigma = sigma_val,
                     N         = N)
  
  if(meas_model == "pois") meas_mdl <- list("y ~ poisson(net_flow(C))")
  
  if(meas_model == "nbin") meas_mdl <- list("y ~ neg_binomial_2(net_flow(C), phi)")
  
  if(meas_model == "norm") meas_mdl <- list("y ~ normal(net_flow(C), epsilon)")
  
  if(!alt_param) {
    
    data_params <- c("N", "par_sigma", "par_gamma")
    
    estimated_params <- list(
      sd_prior("par_beta", "lognormal", c(0, 1)),
      sd_prior("par_rho", "beta", c(2, 2)),
      sd_prior("I0", "lognormal", c(0, 1), "init"))
    
    if(gamma_unk) {
      
      data_params      <- data_params[1:2]      
      gamma_prior      <- list(sd_prior("par_gamma", "beta", c(2, 2)))
      estimated_params <- c(estimated_params, gamma_prior)
    }
    
    stan_text   <- sd_Bayes(fp, meas_mdl, estimated_params, 
                            const_list  = const_list,
                            data_params = data_params,
                            data_inits  = "xi")
  }
  
  if(alt_param) {
    
    estimated_params <- list(sd_prior("par_inv_R0", "beta", c(2, 2)),
                             sd_prior("par_rho", "beta", c(2,2)),
                             sd_prior("I0", "lognormal", c(0,1), "init"))
    
    if(meas_model == "norm") {
      estimated_params <- c(estimated_params, 
                            list(sd_prior("epsilon", "exponential",
                                          0.2, "meas_par")))
    }
    
    stan_text   <- sd_Bayes(fp, meas_mdl, estimated_params, 
                            const_list = const_list, 
                            data_params = c("N", "par_tau", "par_sigma"),
                            data_inits  = "xi")
  }

  stan_fldr     <- file.path(stan_fldr, meas_model)
  dir.create(stan_fldr, showWarnings = FALSE, recursive = TRUE)  
  stan_filepath <- file.path(stan_fldr, str_glue("SE{E_ord}I{I_ord}R.stan"))
  
  output$stan_filepath <- stan_filepath
  
  create_stan_file(stan_text, stan_filepath)
  
  output
}

