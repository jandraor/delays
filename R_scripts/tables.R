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

table_scenarios <- function(df) {
  
  df |> 
    gt()  |> 
    cols_label(R0      = html("\u211c<sub>0</sub>",),
               tau = html("&tau;<sub>e</sub>"))
  
}

table_parameters <- function(df) {
  
  df |> gt() |> 
    tab_style(
      locations = cells_column_labels(columns = everything()),
      style     = list(
        #Give a thick border below
        cell_borders(sides = "bottom", weight = px(3)),
        #Make text bold
        cell_text(weight = "bold")
      ))
}

table_inits <- function(df) {
  
  df |> gt() |> 
    tab_style(
      locations = cells_column_labels(columns = everything()),
      style     = list(
        #Give a thick border below
        cell_borders(sides = "bottom", weight = px(3)),
        #Make text bold
        cell_text(weight = "bold")
      )) |> 
    cols_label(init = "Init value")
}