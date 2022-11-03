#' Generate deterministic incidence from a SE_iI_jR model
#'
#' @param i Number of compartments for the Exposed class
#' @param j Number of compartments for the Infectious class
#' @param N Population size
#' @param beta_val Effective contact rate
#' @param inv_sigma Latent period
#' @param inv_gamma Infectious period
#' @param rho Reporting rate
#' @param tf Stop time
SEIR_actual_incidence <- function(i, j, N, beta_val,inv_sigma, inv_gamma, rho, tf) {
  
  orders <- cross2(i, j) 
  
  map_df(orders, function(ord) {
    
    E_ord <- ord[[1]]
    I_ord <- ord[[2]]
    
    fp <- str_glue("./models/SEIR/SE{E_ord}I{I_ord}R.stmx")
    
    E_inits        <- rep(0, E_ord)
    names(E_inits) <- paste0("E", 1:E_ord)
    
    I_inits        <- c(1, rep(0, I_ord - 1))
    names(I_inits) <- paste0("I", 1:I_ord)
    
    C_init         <- list(C = sum(I_inits))
    S_init         <- list(S = N - sum(I_inits))
    
    stock_inits    <- c(S_init, E_inits, I_inits, C_init)
    
    sigma_val      <- 1 / inv_sigma
    gamma_val      <- 1 / inv_gamma
    
    const_vals     <- list(
      par_beta  = beta_val,
      par_gamma = gamma_val,
      par_sigma = sigma_val,
      par_rho   = rho,
      N         = N)
    
    mdl <- read_xmile(fp, 
                      stock_list = stock_inits,
                      const_list = const_vals)
    
    ds_inputs <- mdl$deSolve_components
    
    output    <- sd_simulate(ds_inputs, timestep = 1/8, integ_method = "rk4", 
                             start_time = 0, stop_time = tf)
    

    inc_df <- sd_net_change(output, "C") |> 
      mutate(x = value) |> 
      select(-var, -value) |> mutate(E_order = as.character(E_ord),
                                     I_order = as.character(I_ord))
  })
}


SEIR_measured_incidence <- function(m, n, N, beta_val, inv_sigma, 
                                    inv_gamma, rho, tf, inv_phi_val) {
  
  actual_incidence <- SEIR_actual_incidence(m, n, N, beta_val, inv_sigma, 
                                            inv_gamma, rho, tf)
  
  actual_incidence |> 
      mutate(y = rnbinom(n(), mu = x, size = 1 / inv_phi_val))
}