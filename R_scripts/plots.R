library(ggplot2)
library(ggpubr)

plot_ts_by_I_order <- function(y_df) {
  
  h_df <- y_df |> mutate(highlight = "1", group = I_order)
  
  orders <- unique(y_df$I_order)
  
  b_df <- map_df(orders, function(o) {
    y_df |> mutate(group   = o,
                   highlight = "0")
  })
  
  
  df <- bind_rows(h_df, b_df) |> 
    mutate(id = paste(I_order, group, highlight, sep = "_"))
  
  
  ggplot(df, aes(x = time, y =  y)) +
    geom_line(aes(colour = highlight, group = id, alpha = highlight)) +
    scale_colour_manual(values = c("grey90", "steelblue")) +
    scale_alpha_manual(values = c(0.5, 1)) +
    facet_wrap(~group) +
    labs(y = "Incidence [Cases/day]", x = "Day") +
    theme_classic() +
    theme(legend.position = "none")
  
  
}

plot_hist_errors <- function(df) {
  
  hist_df <- df |> group_by(error) |> 
    summarise(count = n()) |> 
    mutate(pct = count/sum(count))
  
  
  brks <- seq(0, 0.4, 0.1)
  
  ggplot(hist_df, aes(x = error, y = pct)) +
    geom_col(colour = "steelblue",
             fill   = "white") +
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_x_continuous(breaks = 0:8, labels = 0:8, limits = c(NA, 8)) +
    labs(y = "", x = "Order error") +
    theme_pubr()
  
}

plot_freq_order <- function(df) {
  
  count_df <- df |> group_by(D_n, M_n) |> 
    summarise(count = n())
  
  brks <- 1:9
  
  ggplot(count_df, aes(x = D_n, y = M_n)) +
    geom_point(aes(size = count), colour = "steelblue") +
    scale_y_continuous(breaks = brks, limits = c(1,9)) +
    labs(x = "Data order", y = "Model order") +
    theme_pubr() +
    theme(legend.position = "none")
}

plot_nbin_pars <- function(df, inv_gamma_val, beta_val, phi_val, rho_val) {
  
  R0_val <- beta_val * inv_gamma_val
  
  df <- df |> 
    mutate(`R[0]` = beta * inv_gamma_val,
             highlight = ifelse(D_n == M_n, TRUE, FALSE)) |> 
    select(`R[0]`, phi, rho, highlight, M_n) 
  
  tidy_df <- df |> pivot_longer(-c(M_n, highlight))
  
  yi_df <- data.frame(name = c("R[0]", "phi", "rho"), 
                      y_val = c(R0_val, phi_val, rho_val))
  
  ggplot(tidy_df, aes(x = M_n, y = value)) +
    geom_boxplot(aes(group = M_n, colour = highlight)) +
    geom_hline(data = yi_df, aes(yintercept = y_val), linetype = "dashed", 
               colour = "grey65", alpha = 0.5) +
    scale_colour_manual(values = c("grey50", "steelblue")) +
    scale_x_continuous(breaks = 1:9) +
    facet_wrap(~name, scales = "free", labeller = label_parsed, ncol = 1) +
    labs(x = "Model order", y = "Value") +
    theme_pubr() +
    theme(legend.position = "none")
}

plot_log_lik_by_order <- function(posterior_df) {
  
  posterior_df <- posterior_df |> 
    mutate(data_order = D_n,
              model_order = M_n)
  
  MLE_df <- posterior_df |> 
    select(data_order, model_order, dataset, log_lik) |> 
    group_by(model_order, dataset,) |> 
    filter(log_lik == max(log_lik)) |> unique() |> ungroup()
  
  n_data <- c(1, 2, 5, 10, 15, 20)
  
  map_df(n_data, function(d) {
    MLE_df |> filter(dataset <= d) |> 
      group_by(model_order) |> 
      summarise(ll = sum(log_lik)) |> 
      mutate(d = d)
  }) -> comparison_MLE
  
  comparison_MLE |> group_by(d) |> 
    filter(ll == max(ll)) |> ungroup() -> MLE_per_d
  
  
  ggplot(comparison_MLE, aes(x = model_order, y = ll)) +
    geom_line() +
    facet_wrap(~d, scales = "free") +
    geom_vline(data = MLE_per_d, aes(xintercept = model_order),
               alpha = 0.5, colour = "grey50", linetype = "dotted") +
    scale_x_continuous(breaks = 1:9) +
    labs(y = "Log-likelihood") +
    theme_pubr()
}