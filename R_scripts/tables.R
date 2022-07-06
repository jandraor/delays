library(gt)

table_success <- function(df, n_data) {
  
  df$D_n <- as.factor(df$D_n)
  
  error_df <- df |> 
    group_by(D_n) |> summarise(total_error = sum(error),
                               total_sqr_error = sum(sqr_error))
  
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
      D_n              = "Dn",
      n_successes      = "# of successes", 
      pct_success      = "Success",
      total_error      = "Error",
      total_sqr_error = "Squared error")
}

table_comparison_scores <- function(scores_df) {
  
  scores_df |>
    gt() |> 
    tab_spanner(label   = "Squared error",
                columns = vars(`sqrt_error_LFO (Optimum)`,
                               `sqrt_error_LFO (Default)`, 
                               sqrt_error_MLE)) |> 
    tab_spanner(label   = "Success",
                columns = vars(`pct_success_LFO (Optimum)`,
                               `pct_success_LFO (Default)`, 
                               pct_success_MLE)) |> 
    fmt_percent(
      columns  = vars(pct_success_MLE, 
                      `pct_success_LFO (Default)`,
                      `pct_success_LFO (Optimum)`),
      decimals = 0 
    ) |> 
    cols_label(
      D_n             = html("I<sup>n</sup>"),
      pct_success_MLE = "MLE",
      sqrt_error_MLE  = "MLE",
      `sqrt_error_LFO (Optimum)` = "LFO (Optimum)",
      `pct_success_LFO (Optimum)` = "LFO (Optimum)",
      `sqrt_error_LFO (Default)` = "LFO (Default)",
      `pct_success_LFO (Default)` = "LFO (Default)"
    )
}

table_optim_comparison <- function(wide_optim) {
  
  wide_optim |> 
    gt() |> 
    tab_spanner(label   = "Squared error",
                columns = str_glue("sqrt_error_{1:4}")) |> 
    tab_spanner(label = "Success",
                columns = str_glue("pct_success_{1:4}")) |> 
    fmt_percent(
      columns  = str_glue("pct_success_{1:4}"),
      decimals = 0) |> 
    cols_label(
      sqrt_error_1 = html("I<sup>1</sup>"),
      sqrt_error_2 = html("I<sup>2</sup>"),
      sqrt_error_3 = html("I<sup>3</sup>"),
      sqrt_error_4 = html("I<sup>4</sup>"),
      case         = "Case",
      pct_success_1 = html("I<sup>1</sup>"),
      pct_success_2 = html("I<sup>2</sup>"),
      pct_success_3 = html("I<sup>3</sup>"),
      pct_success_4 = html("I<sup>4</sup>"),)
  
}