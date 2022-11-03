format_posterior <- function(fit_list, info_files, D_i, D_j, ds) {
  
  map2_dfr(fit_list, info_files, function(sf, file_obj) {
    
    M_i           <- file_obj$E_ord
    M_j           <- file_obj$I_ord
    
    sf |> 
      mutate(D_i = D_i, D_j = D_j, M_i = M_i, M_j = M_j, dataset = ds)
  })
}

make_x_summary <- function(fits_by_dataset) {
  
  map_dfr(fits_by_dataset, \(df) {
    
    rnd <- 2
    
    df_list <- split(df, df$M_n)
    
    map_dfr(df_list, \(df2) {
      
      extract_timeseries_var("delta_x_1", df2) |> 
        group_by(time) |> 
        summarise(mean    = mean(value),
                  q_val   = quantile(value, c(0.025, 0.25, 0.5, 0.75, 0.975)),
                  q_type  = c("q2.5", "q25", "q50", "q75", "q97.5"), 
                  .groups = "drop") |> 
        mutate(q_val = round(q_val, rnd),
               mean  = round(mean, rnd)) |> 
        pivot_wider(names_from = q_type, values_from = q_val) |>
        mutate(M_n       = unique(df2$M_n),
               dataset   = unique(df2$dataset))
    })
  })
}

extract_incidence <- function(posterior_df, n_samples) {
  
  posterior_df <- posterior_df |> mutate(M_ij = paste0(M_i, M_j))
  
  df_list      <- split(posterior_df, posterior_df$M_ij)
  
  map_df(df_list, \(df) {
    
    delta_df <- extract_timeseries_var("delta_x_1", df)
    
    max_iter   <- max(delta_df$iter)
    sample_ids <- sample.int(max_iter, n_samples)
    
    delta_df |> filter(iter %in% sample_ids) |> 
      mutate(M_i = unique(df$M_i),
             M_j = unique(df$M_j))
  })
}

posterior_pred_checks <- function(j, n_samples, alt_param) {
  
  estimated_params <- ifelse(alt_param,
                             c("par_inv_R0", "par_rho"),
                             c("par_beta", "par_rho"))
  
  filepath <- ifelse(alt_param,
                     str_glue("./models/SEIR/SE1I{j}R_alt.stmx"),
                     str_glue("./models/SEIR/SE1I{j}R.stmx"))
  N_val    <- 5234
  
  
  if(alt_param) {
    const_list <- list(N         = N_val,
                       par_tau   = 2.85)
  } else {
    const_list <- list(N = N_val,
                       par_gamma = 0.5)
  }
  mdl      <- read_xmile(filepath, const_list = const_list,
                         stock_list = list(R = 0.3 * N_val))
  
  ds_inputs <- mdl$deSolve_components
  
  posterior_df <- fit_df |> filter(M_j == j) |> slice_sample(n = n_samples)
  
  stocks_df <- data.frame(I1 = posterior_df$I0) |> 
    mutate(S = N_val * 0.7 - I1,
           C = I1)
  
  consts_df <- select(posterior_df, all_of(estimated_params))
  
  sens_out <- sd_sensitivity_run(ds_inputs    = ds_inputs,
                                 consts_df    = consts_df,
                                 stocks_df    = stocks_df,
                                 start_time   = 0, 
                                 stop_time    = 91,
                                 timestep     = 1/16,
                                 integ_method = "rk4",multicore = TRUE,
                                 n_cores = 7)
  
  head(sens_out)
  
  split_list <- split(sens_out, sens_out$iter)
  
  map_df(split_list, function(sim_df) {
    
    iter_val <- unique(sim_df$iter)
    
    phi_val <- posterior_df |> slice(iter_val) |> pull(phi)
    sd_net_change(sim_df, "C") |> 
      mutate(iter = iter_val,
             y    = rnbinom(n(), mu = value, size = phi_val),
             M_j  = j)
  }) -> incidence_df
}

calculate_coverage <- function(posterior_list, actual_values, 
                               inv_gamma = NULL) {
  
  map_df(posterior_list, \(posterior_df) {
    
    first_cond <- "par_inv_R0" %in% colnames(posterior_df)
    
    if(first_cond) posterior_df <- mutate(posterior_df, R0 = 1 / par_inv_R0)
    
    second_cond <- "par_gamma" %in% colnames(posterior_df)
    
    if(second_cond) posterior_df <- mutate(posterior_df, 
                                           R0 = par_beta / par_gamma)
    
    
    if(!first_cond & !second_cond) {
      posterior_df<- posterior_df |> mutate(R0 = par_beta * inv_gamma)
    }
    
    par_names <- c("R0", "par_rho", "I0")
    
    c_names <- colnames(posterior_df)
    
    if("par_gamma" %in% c_names) {
      
      par_names     <- c(par_names, "inv_gamma")
      posterior_df  <- posterior_df |> mutate(inv_gamma = 1 / par_gamma)
    }
    
    if("inv_phi" %in% c_names) par_names <- c(par_names, "inv_phi")
    
    pars_df <- select(posterior_df, dataset, D_i, D_j,M_i, M_j, 
                      all_of(par_names)) |> 
      mutate(D_ij = paste0(D_i, D_j),
             M_ij = paste0(M_i, M_j)) |> 
      select(-D_i, -D_j, -M_i, -M_j)
    
    tidy_df <- pars_df |> pivot_longer(c(-dataset, -D_ij, -M_ij),
                                       names_to = "par")
    
    tidy_df |> group_by(dataset, D_ij, M_ij, par) |> 
      summarise(lb = quantile(value, 0.025),
                ub = quantile(value, 0.975),
                .groups = "drop") |> 
      left_join(actual_values, by = "par") |> 
      mutate(caught = ifelse(actual_val >= lb & actual_val <= ub, 1, 0))
    
  }) |> 
    group_by(D_ij, M_ij, par) |> summarise(pct_right = sum(caught) / 20,
                                           .groups = "drop")
}