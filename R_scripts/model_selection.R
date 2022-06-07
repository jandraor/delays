bayesian_choice <- function(ll_df, data_orders, model_orders, n_data,
                            cutoff) {

  combs <- cross2(data_orders, 1:n_data)
  
  map_df(combs, function(id_object) {
    
    do_id <- id_object[[1]]
    ds_id <- id_object[[2]]
    
    orders <- model_orders[-length(model_orders)]
    
    for(i in orders) {
      
      df       <- ll_df |> filter(D_n == do_id, dataset == ds_id)
      current  <- i
      proposal <- current + 1
      
      df_current  <- df |> filter(M_n == current)
      df_proposal <- df |> filter(M_n == proposal)
      
      diffs    <- df_proposal$log_lik - df_current$log_lik
      n_higher <- sum(diffs > 0)
      
      prop_higher <- n_higher / nrow(df_current)
      
      if(prop_higher > cutoff) {
        current <- proposal
      } else{
        break
      }
      
    }
    
    data.frame(D_n = do_id, M_n = current, dataset = ds_id) |> 
      mutate(error = abs(M_n - D_n))
  })
}

MLE_choice <- function(ll_df) {
  
  ll_df |> 
    select(D_n, M_n, dataset, log_lik) |> 
    group_by(D_n, dataset) |> 
    filter(log_lik == max(log_lik)) |> 
    filter(M_n == min(M_n)) |> 
    ungroup() |> 
    unique() |> 
    mutate(error      = abs(M_n - D_n),
           sqrt_error = (M_n - D_n) ** 2)
}

ELDP_choice <- function(n_data, y_df, D_n, M_n_list, all_fits, L, L_min, meas_mdl, L_max = Inf) {
  
  map_df(1:n_data, function(ds) {
    
    flu_data  <- y_df |> filter(I_order == D_n, dataset == ds)
    lt        <- find_last_time(flu_data)
    flu_data  <- flu_data |> filter(time < (lt + 1))
    
    estimate_ELDP(all_fits[[ds]], M_n_list, L, L_min, flu_data, n_samples, 
                  meas_mdl, L_max) |> 
      filter(ELPD == max(ELPD))
    
    
  }) -> choice_df
  
  choice_df |> mutate(error      = abs(D_n - M_n),
                      sqrt_error = (D_n - M_n) ** 2,
                      is_correct = ifelse(D_n == M_n, TRUE, FALSE))
}

evaluate_L <- function(n_data, y_df, D_n, M_n_list,all_fits, L_min, L_max,
                       meas_mdl, fp) {
  
  if(!file.exists(fp)) {
    
    L_mins <- L_min:L_max
    L_maxs <- L_min:L_max
    
    combs <- cross2(L_mins, L_maxs) |> 
      
      lapply(function(obj) {
        
        diff <- obj[[2]] - obj[[1]]
        if(diff < 1) return(NULL)
        
        obj
      }) |> compact()
    
    print(length(combs))
    
    plan(multisession, workers = parallel::detectCores() - 1)
    
    with_progress({
      
      p <- progressor(steps = length(combs))
      
      future_map_dfr(combs, function(comb_obj) {
        
        p()
        
        L     <- comb_obj[[1]]
        L_max <- comb_obj[[2]]
        
        choice_df <- ELDP_choice(n_data, y_df, D_n, M_n_list, all_fits, L, 
                                 L_min, meas_mdl, L_max)
        
        data.frame(L_min       = L, 
                   L_max       = L_max,
                   error       = sum(choice_df$error),
                   sqrt_error  = sum(choice_df$sqrt_error),
                   n_successes = sum(choice_df$is_correct))
        
      }) -> result
    })
    
    saveRDS(result, fp)
    
    
    
  } else {
    result <- readRDS(fp)
  }
   
  result
}


  
 