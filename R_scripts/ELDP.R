estimate_ELDP <- function(fits, M_n_list, L, L_min, flu_data, n_samples, meas_mdl, 
                          L_max) {
  
  D_m <- unique(flu_data$E_order)
  D_n <- unique(flu_data$I_order)
  ds  <- unique(flu_data$dataset)
  
  imap_dfr(M_n_list, function(M_n, j) {
    
    fit_order <- fits[[j]]
       
    loglik_exact     <- pred_scores(fit_order, flu_data, L, n_samples, meas_mdl,
                                    L_min, L_max)
    exact_elpds_1sap <- apply(loglik_exact, 2, log_mean_exp)
    exact_elpd_1sap  <- c(ELPD = sum(exact_elpds_1sap[-(1:L)]))
    
    data.frame(D_m     = D_m,
               D_n     = D_n,
               M_m     = 1,
               M_n     = M_n, 
               dataset = ds,
               ELPD    = exact_elpd_1sap,
               L       = L)
  })
}

pred_scores <- function(fit_order, flu_data, L, n_samples, meas_mdl, L_min, L_max) {
  
  N <- min(nrow(flu_data), L_max)
  
  loglik_exact <- matrix(nrow = n_samples, ncol = N)
  
  for(pred_i in (L + 1):N) {
    
    posterior_df <- fit_order[[pred_i - L_min]]
    
    if(meas_mdl == "nbin") {
      
      loglik <- dnbinom(x    = flu_data$y[pred_i], 
                        mu   = posterior_df$pred, 
                        size = 1 / posterior_df$phi,
                        log  = TRUE)
    }
    
    if(meas_mdl == "pois") {
      
      loglik <- dpois(x      = flu_data$y[pred_i], 
                      lambda = posterior_df$pred,
                      log  = TRUE)
    }
    
    loglik_exact[, pred_i] <- loglik
  }
  
  loglik_exact
}