format_posterior <- function(fit_list, info_files, D_m, D_n, ds) {
  
  map2_dfr(fit_list, info_files, function(sf, file_obj) {
    
    M_m           <- file_obj$E_ord
    M_n           <- file_obj$I_ord
    
    sf |> 
      mutate(D_m = D_m, D_n = D_n, M_m = M_m, M_n = M_n, dataset = ds)
  })
}

make_x_summary <- function(fits_by_dataset) {
  
  map_dfr(fits_by_dataset, \(df) {
    
    rnd <- 2
    
    df_list <- split(df, df$M_n)
    
    map_dfr(df_list, \(df2) {
      
      extract_timeseries_var("x", df2) |> 
        group_by(time) |> 
        summarise(mean   = mean(value),
                  q_val  = quantile(value, c(0.025, 0.25, 0.5, 0.75, 0.975)),
                  q_type = c("q2.5", "q25", "q50", "q75", "q97.5")) |> 
        ungroup() |> 
        mutate(q_val = round(q_val, rnd),
               mean  = round(mean, rnd)) |> 
        pivot_wider(names_from = q_type, values_from = q_val) |>
        mutate(M_n       = unique(df2$M_n),
               dataset   = unique(df2$dataset))
    })
  })
}