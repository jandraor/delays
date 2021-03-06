library(ggplot2)
library(ggpubr)
library(patchwork)
library(viridis)

plot_ts_by_I_order <- function(y_df) {
  
  h_df <- y_df |> mutate(highlight = "1", group = I_order)
  
  orders <- unique(y_df$I_order)
  
  b_df <- map_df(orders, function(o) {
    y_df |> mutate(group   = o,
                   highlight = "0")
  })
  
  
  df <- bind_rows(h_df, b_df) |> 
    mutate(id = paste(I_order, group, highlight, sep = "_"))
  
  
  ggplot(df, aes(x = time, y =  x)) +
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

plot_error_by_L <- function(df, subtitle = "") {
  
  ggplot(df, aes(x = L_min, y = L_max)) +
    geom_tile(aes(fill = sqrt(error))) +
    scale_fill_viridis(discrete = FALSE) +
    labs(subtitle = subtitle) +
    theme_pubr()
      
}

plot_min_error <- function(df, tol = 0, error_var = "error") {
  
  if(error_var == "error") {
    df_min <- df |> filter(error <= min(error) + tol) |> 
      mutate(id = row_number())
  }
  
  if(error_var == "sqr_error") {
    df_min <- df |> filter(sqr_error <= min(sqr_error) + tol) |> 
      mutate(id = row_number())
  }
  
  ggplot(df_min) +
    geom_errorbar(aes(x = id, ymin = L_min, ymax = L_max), colour = "steelblue") +
    scale_y_continuous(limits = c(0, 60)) +
    coord_flip() +
    labs(y = "Time range", x= "") +
    theme_pubr()
}

plot_incidence <- function(y_list, dist) {
  
  y_all <- do.call(rbind, y_list)
  
  ggplot(y_all, aes(time, y)) +
    geom_line(aes(group = dataset), colour = "grey80", alpha = 0.10) +
    facet_wrap("I_order") +
    labs(x        = "Days", 
         y        = "Incidence [Cases/day]",
         subtitle = dist) +
    theme_pubr()
}

plot_x_prediction_facets <- function(x_summary, y_subset) {
  
  ggplot(x_summary, aes(time, q50)) +
    geom_line(colour = "steelblue") +
    geom_ribbon(alpha = 0.1, aes(ymin = q2.5, ymax = q97.5),
                fill = "steelblue") +
    geom_ribbon(alpha = 0.2, aes(ymin = q25, ymax = q75),
                fill = "steelblue") +
    geom_point(data = y_subset, aes(time, x), size = 0.5, colour = "grey50") +
    facet_grid(M_n ~ dataset) +
    theme_pubr()
}

plot_x_prediction <- function(x_summary, y_subset) {
  
  ggplot(x_summary |> mutate(M_n = as.factor(M_n)), 
         aes(time, q50)) +
    geom_line(aes(group = M_n, colour = M_n)) +
    geom_point(data = y_subset, aes(time, x), size = 0.5, colour = "grey50") +
    facet_wrap(~ dataset) +
    theme_pubr()
}

plot_CV_boundaries <- function(y_df, boundaries_df) {
  
  demo_df <- y_df |> filter(dataset == 1) |> 
    select(time, x, I_order) |> 
    left_join(boundaries_df, by = "I_order") |> 
    mutate(fill  = ifelse(time >= ll & time <= ul, TRUE, FALSE),
           point = ifelse(time == ll | time == ul, TRUE, FALSE))
  
  point_df <- demo_df |> filter(point == TRUE)
  
  ggplot(demo_df, aes(time, x)) +
    geom_area(aes(fill = fill), show.legend = FALSE) +
    geom_line(colour = "steelblue", size = 0.1) +
    geom_point(data = point_df, size = 0.1, colour = "orange") +
    scale_fill_manual(values = c("white", "grey95")) +
    facet_wrap("I_order", ncol = 1) +
    theme_pubr()
}

plot_CV_boundary <- function(y_df, I_order, ll, ul) {
  
  demo_df <- y_df |> filter(dataset == 1, I_order == !!I_order ) |> 
    select(time, x, I_order) |> 
    mutate(fill  = ifelse(time >= ll & time <= ul, TRUE, FALSE),
           point = ifelse(time == ll | time == ul, TRUE, FALSE))
  
  point_df <- demo_df |> filter(point == TRUE)
  
  ggplot(demo_df, aes(time, x)) +
    geom_area(aes(fill = fill)) +
    geom_line(colour = "steelblue", size = 0.1) +
    geom_point(data = point_df, size = 0.1, colour = "orange") +
    scale_fill_manual(values = c("white", "grey95")) +
    theme_pubr() +
    theme(legend.position = "none")
}

plot_par_comparison <- function(df, par_name, actual_value, D_m, D_n, x_label) {
  
  capt <- paste(str_glue("True E order: {D_m}"),
                str_glue("True I order: {D_n}"),
                sep = "\n")
  
  df <- df |> 
    mutate(M_n = str_glue("I^{M_n}"))
  
  ggplot(df, aes(.data[[par_name]])) +
    geom_histogram(aes(fill = as.factor(M_m)), colour = "white", alpha = 0.6) +
    facet_grid(M_n ~ dataset, labeller = label_parsed) +
    geom_vline(xintercept = actual_value, colour = "grey85", linetype = "dotdash") +
    scale_fill_manual(values = c("steelblue", "grey60"), name = "E order") +
    labs(y = "", caption = capt, x = parse(text = x_label)) +
    theme_pubr() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.line.y  = element_blank())
}

plot_prediction_comparison <- function(ELDP_default, ELDP_best, MLE_results, D_n) {
  
  MLE_predictions <- MLE_results |> 
    filter(D_n == !!D_n) |> 
    select(M_n, n, window)
  
  ELDP_default <- ELDP_default |> mutate(window = "LFO-CV (all points)")
  ELDP_best    <- ELDP_best |> mutate(window = "LFO-CV (optim)")
  
  ELDP_df <- bind_rows(ELDP_default, ELDP_best)
  
  summary_predictions <- ELDP_df |> 
    group_by(M_n, window) |> 
    summarise(n = n(), .groups = "drop") |> 
    bind_rows(MLE_predictions) |> 
    mutate(is_correct = ifelse(D_n == M_n, TRUE, FALSE))
  
  ggplot(summary_predictions, aes(M_n, n)) +
    geom_col(aes(fill = is_correct), colour = "white") +
    scale_fill_manual(values = c("grey65", "steelblue")) +
    facet_wrap(~window) +
    scale_x_continuous(limits = c(0, 5)) +
    labs(y = "Count", x = bquote(I^n)) +
    theme_pubr()
}