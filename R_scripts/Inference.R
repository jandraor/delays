fit_SEIR <- function(info_files, D_m, D_n, ds, epi_df, root_fldr, n_chains, 
                    samples_per_chain) {
  
  lapply(info_files, function (info) {
    
    M_m           <- info$E_ord
    M_n           <- info$I_ord
    n_stocks      <- info$n_stocks
    stan_filepath <- info$stan_filepath
    
    backup_folder <- str_glue("{root_fldr}/D{D_m}{D_n}/dataset_{ds}/M{M_m}{M_n}")
    fn            <- file.path(backup_folder, "fit.rds") 
    
    if(!file.exists(fn)) {
      
      N <- nrow(epi_df)
      
      stan_d <- list(n_obs    = N,
                     y        = epi_df$y,
                     n_params = 2,
                     n_difeq  = n_stocks,
                     t0       = 0,
                     ts       = 1:(N + 1))
      
      mod <- cmdstan_model(stan_filepath)
      
      tic.clearlog()
      tic()
      

      fit <- mod$sample(data            = stan_d,
                        chains          = n_chains,
                        parallel_chains = 4,
                        iter_warmup     = samples_per_chain,
                        iter_sampling   = samples_per_chain,
                        adapt_delta     = 0.9,
                        refresh         = 5,
                        save_warmup     = FALSE)
      
      toc(quiet = FALSE, log = TRUE)
      log.lst <- tic.log(format = FALSE)
      
      time_file    <- file.path(backup_folder, "time.csv")
      time_elapsed <- calculate_time(log.lst)
      
      record_df    <- data.frame(D_m = D_m,
                                 D_n = D_n,
                                 M_m = M_m,
                                 M_n = M_n,
                                 dataset = ds,
                                 time    = time_elapsed)
      
      write_csv(record_df, time_file)
      
      diag_path <- file.path(backup_folder, str_glue("diag.txt"))
      diagnosis <- fit$cmdstan_diagnose()
      writeLines(diagnosis$stdout, diag_path)
      
      sf           <- rstan::read_stan_csv(fit$output_files())
      posterior_df <- as.data.frame(sf)
      
      subset_df    <- select(posterior_df, beta, I0, rho, phi, log_lik, pred, 
                             starts_with("x[")) 
        
      saveRDS(subset_df , fn)
      
    } else {
      subset_df  <- readRDS(fn)
    }
    
    subset_df 
  })
}

