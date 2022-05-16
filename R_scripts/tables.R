library(gt)

table_success <- function(df, n_data) {
  
  df$D_n <- as.factor(df$D_n)
  
  error_df <- df |> 
    group_by(D_n) |> summarise(total_error = sum(error))
  
  success_df <- df |> filter(error == 0) |> 
    count(D_n, name = "n_successes", .drop = FALSE) |> 
    mutate(pct_success = n_successes / n_data)
  
  summary_df <- full_join(success_df, error_df, by = "D_n")
  
  summary_df |> gt() |>  
    fmt_percent(
      columns  = pct_success,
      decimals = 0 
    ) |> 
    cols_label(
      D_n         = "Dn",
      n_successes = "# of successes", 
      pct_success = "Success",
      total_error = "Error")
}