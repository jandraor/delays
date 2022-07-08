CV_datasets <- function(n_data, y_df, D_n, M_n, L_min, n_chains, 
                        samples_per_chain, root_fldr, meas_mdl) {
  
  lapply(1:n_data, function(ds) {
    
    flu_data <- y_df |> filter(I_order == D_n, dataset == ds)
    lt       <- find_last_time(flu_data)
    #print(lt)
    flu_data <- flu_data |> filter(time < (lt + 1))
    
    fits <- run_LFO_CV_SEIR(info_files[M_n], flu_data, L_min, n_chains, 
                            samples_per_chain, root_fldr, meas_mdl)
  })
}

run_LFO_CV_SEIR <- function(info_files, flu_data, L_min, n_chains, 
                           samples_per_chain, root_folder, meas_mdl) {
  
  lapply(info_files, function(info) {
    
    folder_id       <- paste0(sample(c(0:9, LETTERS), 8, T), collapse = '')
    stan_fldr       <- info$stan_filepath |> dirname()
    stan_file_name  <- info$stan_filepath |> basename()
    new_folder_path <- file.path(stan_fldr, folder_id)
    dir.create(new_folder_path, recursive = TRUE)
    
    new_file_path   <- file.path(new_folder_path, stan_file_name)
    
    file.copy(info$stan_filepath, new_file_path)
    
    info$stan_filepath <- new_file_path
    
    fit <- pointwise_pred(flu_data, info, root_folder, meas_mdl, L_min)
    
    unlink(new_folder_path, recursive = TRUE)
    
    fit
  }) -> fits
}

pointwise_pred <- function(flu_data, info, root_folder, meas_mdl, L_min) {
  
  M_m           <- info$E_ord
  M_n           <- info$I_ord
  n_stocks      <- info$n_stocks
  stan_filepath <- info$stan_filepath 
  
  N    <- nrow(flu_data)
  D_m  <- unique(flu_data$E_order)
  D_n  <- unique(flu_data$I_order)
  ds   <- unique(flu_data$dataset)
  
  backup_folder <- str_glue("{root_folder}/D{D_m}{D_n}/dataset_{ds}/M{M_m}{M_n}")
  
  lapply(L_min:(N - 1), function(i) {
    
    n_obs <- i
    
    fn  <- file.path(backup_folder, paste0("fit_", i)) 
    
    if(!file.exists(fn)) {
      
      stan_d <- list(n_obs  = n_obs,
                     y        = flu_data$y[1:n_obs],
                     n_params = 2,
                     n_difeq  = n_stocks,
                     t0       = 0,
                     ts       = 1:(n_obs + 1))
      
      mod <- cmdstan_model(stan_filepath)
      
      fit <- mod$sample(data            = stan_d,
                        chains          = n_chains,
                        parallel_chains = n_chains,
                        iter_warmup     = samples_per_chain,
                        iter_sampling   = samples_per_chain,
                        refresh         = 100,
                        save_warmup     = FALSE)
      
      diag_path <- file.path(backup_folder, str_glue("fit_{i}.txt"))
      diagnosis <- fit$cmdstan_diagnose()
      writeLines(diagnosis$stdout, diag_path)
      
      sf  <- rstan::read_stan_csv(fit$output_files())
      
      posterior_df <- as.data.frame(sf)
      
      select_vars <- c("pred", "beta", "rho", "I0")
      
      if(meas_mdl == "nbin") select_vars <- c(select_vars, "phi")
      
      subset_df <- posterior_df[, select_vars]
      
      saveRDS(subset_df, fn)
    } else {
      subset_df  <- readRDS(fn)
    }
    subset_df 
  }) -> result
  
  result
}
