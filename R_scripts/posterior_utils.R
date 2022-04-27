format_posterior <- function(fit_list, info_files, D_m, D_n, ds) {
  
  map2_dfr(fit_list, info_files, function(sf, file_obj) {
    
    M_m           <- file_obj$E_ord
    M_n           <- file_obj$I_ord
    
    sf |> 
      mutate(D_m = D_m, D_n = D_n, M_m = M_m, M_n = M_n, dataset = ds)
  })
}