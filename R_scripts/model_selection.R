MLE_choice <- function(ll_df) {
  
  ll_df |> 
    select(D_j, M_i, M_j, dataset, log_lik) |> 
    group_by(D_j, dataset) |> 
    filter(log_lik == max(log_lik)) |> 
    filter(M_j == min(M_j)) |> 
    ungroup() |> 
    unique() |> 
    mutate(error      = abs(M_j - D_j),
           sqr_error = (M_j - D_j) ** 2)
}
