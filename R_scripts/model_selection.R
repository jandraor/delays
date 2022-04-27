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
    mutate(error = abs(M_n - D_n))
}


  
 