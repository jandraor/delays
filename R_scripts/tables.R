library(gt)

table_coverage <- function(df) {
  
  df <- relocate(df, R0, .after = M_ij)
  
  if("inv_phi" %in% colnames(df)) df <- relocate(df, inv_phi, 
                                                 .after = par_rho)
  
  names_col <- colnames(df)
  n_cols    <- ncol(df)
  par_names <- names_col[3:n_cols]
  
  df |> gt() |> 
    fmt_percent(
      columns = par_names,
      decimals = 0
    ) |> 
    cols_label(M_ij    = html("M<sup>ij</sup>"),
               D_ij    = html("D<sup>ij</sup>"),
               I0      = html("I<sub>0</sub>"),
               R0      = html("\u211c<sub>0</sub>",),
               par_rho = html("&rho;")) |> 
    data_color(
      columns = par_names,
      colors = scales::col_numeric(
        palette = c("#8EDACE", "#50BEAD", "#2AA996", "#216176"),
        domain = c(0, 1)
      )
    ) |> 
    tab_style(
      cell_borders(
        sides = c("bottom"),
        color = "#7f7f7f",
        weight = px(2)
      ),
      locations = list(
        cells_body(
          rows = M_ij == "14")
    )) -> gt_table
  
  if("inv_gamma" %in% colnames(df)) {
    
    gt_table <- gt_table |> cols_label(inv_gamma = html("&gamma;<sup>-1</sup>"))
  }
  
  if("inv_phi" %in% colnames(df)) {
    
    gt_table <- gt_table |> cols_label(inv_phi = html("&phi;<sup>-1</sup>"))
  }
  
  gt_table
}

table_success <- function(df, n_data) {
  
  df <- df |> 
    mutate(success       = ifelse(D_j == M_j, 1, 0),
           success_dist  = ifelse(D_dist == M_dist, 1, 0),
           D_j           = as.factor(D_j))
  
  summary_df <- df |> group_by(D_j) |> 
    summarise(n_successes      = sum(success),
              n_successes_dist = sum(success_dist),
              .groups          = "drop") |> 
    mutate(pct_success      = n_successes / n_data,
           pct_success_dist = n_successes_dist / n_data)
  

  
  summary_df |> gt() |>  
    fmt_percent(
      columns  = c(pct_success, pct_success_dist),
      decimals = 0 
    ) |> 
    cols_label(
      D_j              = "Dj",
      n_successes      = "# of successes", 
      n_successes_dist = "# of successes (dist)",
      pct_success      = "Success",
      pct_success_dist = "Success (dist)"
     )
}
 
table_scenarios <- function(df) {
  
  df |> 
    gt()  |> 
    cols_label(R0      = html("\u211c<sub>0</sub>",),
               tau = html("&tau;"))
  
}