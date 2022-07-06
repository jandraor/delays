create_stan_file <- function(stan_text, filename) {
  file_connection <- file(filename)
  writeLines(stan_text, file_connection)
  close(file_connection)
}

create_backup_fldrs <- function(D_m, D_n, n_data, M_m, M_n, root_fldr) {
  
  data_ids <- 1:n_data
  
  id_list  <- cross3(data_ids,  M_m, M_n)
  
  walk(id_list, function(id_obj) {
    
    ds   <- id_obj[[1]] # dataset id
    M_m  <- id_obj[[2]]
    M_n  <- id_obj[[3]]
    
    fldr <- str_glue("{root_fldr}/D{D_m}{D_n}/dataset_{ds}/M{M_m}{M_n}")
    dir.create(fldr, recursive = TRUE, showWarnings = FALSE)
  })
}

calculate_time <- function(t_list) {
  t_obj <- t_list[[1]]
  (t_obj$toc - t_obj$tic) / 60
}

find_last_time <- function(df) {
  
  df <- df |> arrange(desc(time)) |> 
    mutate(cml_y = cumsum(y)) |> 
    mutate(flag = ifelse(cml_y == 0, TRUE, FALSE)) |> 
    filter(flag == FALSE)
  
  df |> slice(1) |> pull(time)
}

find_first_time <- function(df) {
  
  df <- df |> 
    mutate(cml_y = cumsum(y)) |> 
    mutate(flag = ifelse(cml_y == 0, TRUE, FALSE)) |> 
    filter(flag == FALSE)
  
  df |> slice(1) |> pull(time)
}

save_MLE_choice <- function(MLE_df, root_fldr) {
  
  MLE_df |> group_by(D_n, M_n) |> 
    summarise(n = n(), .groups = "drop") |> 
    mutate(window = "MLE") |> 
    write_csv(file.path(root_fldr, "MLE_choice.csv"))
}

estimate_MLE_score <- function(MLE_results) {
  
  MLE_results |> 
    mutate(sqr_error = ((D_n - M_n) ** 2) * n,
           is_correct = D_n == M_n) |> 
    group_by(window, D_n) |> 
    mutate(n_correct = n * is_correct) |> 
    summarise(sqrt_error  = sum(sqr_error),
              pct_success = n_correct / sum(n),
              .groups     = "drop") |> 
    group_by(window, D_n, sqrt_error) |> 
    summarise(pct_success = sum(pct_success), .groups = "drop")
}