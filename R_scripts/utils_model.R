get_stock_inits <- function(E_ord, I_ord, N) {
  
  E_inits        <- rep(0, E_ord)
  names(E_inits) <- paste0("E", 1:E_ord)
  
  I_inits        <- c(1, rep(0, I_ord - 1))
  names(I_inits) <- paste0("I", 1:I_ord)
  
  C_init         <- list(C = sum(I_inits))
  S_init         <- list(S = N - sum(I_inits))
  
  stock_inits    <- c(S_init, E_inits, I_inits, C_init)
  
}