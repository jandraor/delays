fit_SEIR <- function(info_files, D_i, D_j, ds, epi_df, root_fldr, n_chains, 
                    samples_per_chain, n_params, alt_param, inv_gamma,
                    inv_sigma, set_gamma = TRUE) {
  
  lapply(info_files, function (info) {
    
    M_i           <- info$E_ord
    M_j           <- info$I_ord
    n_stocks      <- info$n_stocks
    backup_folder <- str_glue("{root_fldr}/D{D_i}{D_j}/dataset_{ds}/M{M_i}{M_j}")
    fn            <- file.path(backup_folder, "fit.rds") 
    
    if(!file.exists(fn)) {
      
      cat("\n M_j: ", M_j, "\n")
      cat("\n ds: ", ds, "\n")
      
      stan_filepath   <- info$stan_filepath
      folder_id       <- paste0(sample(c(0:9, LETTERS), 8, T), collapse = '')
      stan_fldr       <- stan_filepath |> dirname()
      stan_file_name  <- stan_filepath |> basename()
      new_folder_path <- file.path(stan_fldr, folder_id)
      dir.create(new_folder_path, recursive = TRUE)
      
      new_file_path   <- file.path(new_folder_path, stan_file_name)
      
      file.copy(stan_filepath, new_file_path)
      
      stan_filepath <- new_file_path
      
      N <- nrow(epi_df)
      
      stan_d <- list(n_obs     = N,
                     y         = epi_df$y,
                     n_params  = n_params,
                     n_difeq   = n_stocks,
                     t0        = 0,
                     ts        = 1:N,
                     N         = 1e4,
                     xi        = 0,
                     par_sigma = 1 / inv_sigma)
      
      if(set_gamma) stan_d$par_gamma <- 1 / inv_gamma 
      
      if(alt_param) stan_d$par_tau <- mean_generation_time(j         = D_j, 
                                                           sigma_inv = inv_sigma, 
                                                           gamma_inv = inv_gamma)
      
      mod <- cmdstan_model(stan_filepath)
      
      tic.clearlog()
      tic()

      fit <- mod$sample(data            = stan_d,
                        chains          = n_chains,
                        parallel_chains = 4,
                        iter_warmup     = samples_per_chain,
                        iter_sampling   = samples_per_chain,
                        refresh         = 100,
                        save_warmup     = FALSE)
      
      toc(quiet = FALSE, log = TRUE)
      log.lst <- tic.log(format = FALSE)
      
      time_file    <- file.path(backup_folder, "time.csv")
      time_elapsed <- calculate_time(log.lst)
      
      record_df    <- data.frame(D_i = D_i,
                                 D_j = D_j,
                                 M_i = M_i,
                                 M_j = M_j,
                                 dataset = ds,
                                 time    = time_elapsed)
      
      write_csv(record_df, time_file)
      
      diag_path <- file.path(backup_folder, str_glue("diag.txt"))
      diagnosis <- fit$cmdstan_diagnose()
      writeLines(diagnosis$stdout, diag_path)
      
      unlink(new_folder_path, recursive = TRUE)
      
      posterior_df <- posterior::as_draws_df(fit$draws())
      
      subset_df    <- select(posterior_df, contains("par_beta"), par_rho, 
                             log_lik, starts_with("delta_x_1["), 
                             contains("phi"), contains("I0"), 
                             contains("par_gamma"), contains("par_inv_R0"),
                             contains("epsilon")) 
     
      saveRDS(subset_df, fn)
      
    } else {
      subset_df  <- readRDS(fn)
    }
    
    subset_df 
  })
}

fit_datasets <- function(D_i, D_j, n_data, y_df, info_files, root_fldr,
                         n_chains, samples_per_chain, n_params, 
                         alt_param = FALSE, inv_gamma = 2,
                         inv_sigma = 2) {
  
  
  lapply(1:n_data, function(ds) {

    flu_data <- y_df |> filter(E_order == D_i,
                               I_order == D_j,
                               dataset == ds)
    
    lt       <- find_last_time(flu_data)
    
    flu_data <- flu_data |> filter(time < (lt + 1))
    
    fit_list <- fit_SEIR(info_files, D_i, D_j, ds, flu_data, root_fldr, 
                         n_chains, samples_per_chain, n_params, alt_param,
                         inv_gamma, inv_sigma)
    
    format_posterior(fit_list, info_files, D_i, D_j, ds)
  })
}

fit_single_dataset <- function(info_files, epi_df, root_fldr, n_chains, 
                               samples_per_chain, n_params, alt_param, N = 1e4,
                               xi = 0, fixed_tau = NULL) {
  
  map_df(info_files, function (info) {
    
    n_stocks      <- info$n_stocks
    M_i           <- info$E_ord
    M_j           <- info$I_ord
    backup_folder <- str_glue("{root_fldr}/M{M_i}{M_j}")
    if(!dir.exists(backup_folder)) dir.create(backup_folder, recursive = TRUE)
    fn            <- file.path(backup_folder, "fit.rds") 
    
    if(!file.exists(fn)) {
      
      stan_filepath <- create_temp_folder(info)
      
      n_obs <- nrow(epi_df)
      
      stan_d <- list(n_obs    = n_obs,
                     y        = epi_df$y,
                     n_params = n_params,
                     n_difeq  = n_stocks,
                     t0       = 0,
                     ts       = 1:n_obs,
                     N        = N,
                     xi       = xi)
      
      if(alt_param) stan_d$par_tau <- fixed_tau
      
      mod <- cmdstan_model(stan_filepath)
      
      fit <- mod$sample(data            = stan_d,
                        chains          = n_chains,
                        parallel_chains = 4,
                        iter_warmup     = samples_per_chain,
                        iter_sampling   = samples_per_chain,
                        refresh         = 100,
                        save_warmup     = FALSE)
      
      diag_path <- file.path(backup_folder, str_glue("diag.txt"))
      diagnosis <- fit$cmdstan_diagnose()
      writeLines(diagnosis$stdout, diag_path)
      
      unlink(dirname(stan_filepath), recursive = TRUE)
      
      posterior_df <- posterior::as_draws_df(fit$draws())
      subset_df    <- get_subset_posterior(posterior_df)
      
      saveRDS(subset_df, fn)
      } else {
      subset_df  <- readRDS(fn)
    }
    
    subset_df |> mutate(M_i = M_i, M_j = M_j)
  })
}

create_temp_folder <- function(info) {
  
  stan_filepath   <- info$stan_filepath
  folder_id       <- paste0(sample(c(0:9, LETTERS), 8, T), collapse = '')
  stan_fldr       <- stan_filepath |> dirname()
  stan_file_name  <- stan_filepath |> basename()
  new_folder_path <- file.path(stan_fldr, folder_id)
  dir.create(new_folder_path, recursive = TRUE)
  
  new_file_path   <- file.path(new_folder_path, stan_file_name)
  
  file.copy(stan_filepath, new_file_path)
  
  new_file_path
}

get_subset_posterior <- function(posterior_df) {
  
  select(posterior_df, contains("par_beta"), par_rho, 
         log_lik, starts_with("delta_x_1["), 
         contains("phi"), contains("I0"), 
         contains("par_gamma"), contains("par_inv_R0"),
         contains("epsilon")) 
}

