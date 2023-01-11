library(GGally)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(scales)
library(viridis)

colours_df <- data.frame(order = 1:4, 
                         colour = c("#C7D9F4", "#304460", "#BCA7D2", "#4C808A"))


plot_latent_incidence <- function(x_df) {
  
  x_df <- x_df |> mutate(label_i = paste0("i = ", E_order),
                         label_j = paste0("j = ", I_order)) 
  
  comparison_df <- x_df |> 
    mutate(label_i = case_when(label_i == "i = 1" ~ "i = 3",
                               label_i == "i = 3" ~ "i = 1"))
 
  ggplot(x_df, aes(time, x)) +
    geom_line(data = comparison_df, aes(linetype = E_order), colour = "grey90",
              alpha = 0.75) +
    geom_line(aes(colour = I_order, linetype = E_order)) +
    facet_grid(label_j ~ label_i) +
    scale_colour_manual(values = colours_df$colour) +
    scale_linetype_manual(values = c("solid", "longdash")) +
    labs(y = "Incidence rate [New cases/day]", x = "Day") +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_rect(colour = "grey80"),
          axis.title = element_text(colour = "grey40"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),)
}

plot_meas_incidence_by_overdispersion <- function(df) {
  
  M_i <- unique(df$E_order)
  
  subtitle_txt <- parse(text = paste0("'Time-series:' ~ D^{'", M_i, "j'}"))
  
  df <- df |> mutate(label_j = paste0("'j = '~", I_order)) 
  
  ggplot(df, aes(time, y)) +
    geom_line(aes(group = dataset, colour = I_order, alpha = highlight)) +
    scale_colour_manual(values = colours_df$colour) +
    scale_alpha_manual(values = c(0.05, 1)) +
    facet_grid(label_j~disp, labeller = label_parsed) +
    labs(y = "Measured incidence [New cases/day]", x = "Day",
         subtitle = parse(text = subtitle_txt)) +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_rect(colour = "grey80"),
          axis.title = element_text(colour = "grey40"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          plot.subtitle = element_text(colour = "grey45"))
}

plot_hist_mase <- function(posterior_sample, y_df, lims = c(0, 0.5)) {
  
  actual_E <- unique(posterior_sample$D_i)
  actual_I <- unique(posterior_sample$D_j)
  ds       <- unique(posterior_sample$dataset)
  
  act_x <- y_df |> filter(E_order == actual_E, I_order == actual_I, 
                          dataset == ds) |> pull(x)
  
  posterior_mase_x <- extract_incidence(posterior_sample, 4000) |> 
    group_by(iter, M_j) |> summarise(mase = mase(act_x, value),
                                     .groups = "drop") |> 
    mutate(label_Mj = paste0("j = ", M_j),
           data = "x")
  
  act_y <- y_df |> filter(E_order == actual_E, I_order == actual_I, 
                        dataset == ds) |> pull(y)
  
  posterior_mase_y <- extract_incidence(posterior_sample, 4000) |> 
    group_by(iter, M_j) |> summarise(mase = mase(act_y, value),
                                     .groups = "drop") |> 
    mutate(label_Mj = paste0("j = ", M_j),
           data = "y")
  
  posterior_mase <- bind_rows(posterior_mase_x, posterior_mase_y)
  
  M_i <- unique(posterior_sample$M_i)
  
  sub_txt <- paste0("'Histograms:' ~ M^{'", M_i, "j'}")
  
  mean_df <- posterior_mase |> group_by(label_Mj, data) |> 
    summarise(mean_mase = mean(mase))
  
  ggplot(posterior_mase, aes(mase)) +
    scale_x_continuous(limits = lims) +
    geom_histogram(aes(fill = as.factor(M_j)), alpha = 0.9) +
    scale_fill_manual(values = colours_df$colour) +
    facet_grid(label_Mj~data, scales = "free_x") +
    geom_vline(data = mean_df, aes(xintercept = mean_mase, colour = label_Mj)) +
    scale_colour_manual(values = colours_df$colour) +
    labs(x        = "Mean absolute scale error (MASE)", y = "",
         caption  = paste0("Dataset: ", ds, "\n", "Line: Mean value"),
         subtitle = parse(text = sub_txt)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.line  = element_line(colour = "grey80"),
          strip.background = element_rect(colour = "grey80"),
          axis.ticks = element_line(colour = "grey60"),
          axis.text  = element_text(colour = "grey60"))
}

plot_hist_mase2 <- function(posterior_sample, y_df, lims = c(0, 0.5)) {
  
  actual_E <- unique(posterior_sample$D_i)
  actual_I <- unique(posterior_sample$D_j)
  ds       <- unique(posterior_sample$dataset)
  
  act_x <- y_df |> filter(E_order == actual_E, I_order == actual_I, 
                          dataset == ds) |> pull(x)
  
  lbls_ordered <- str_glue("i = {rep(c(1, 3), 4)}, j = {rep(1:4, each = 2)}")
  
  posterior_mase_x <- extract_incidence(posterior_sample, 4000) |> 
    group_by(iter, M_i, M_j) |> summarise(mase = mase(act_x, value),
                                     .groups = "drop") |> 
    mutate(label_Mij = paste0("i = ", M_i,", ","j = ", M_j),
           data = "x",
           label_Mij = factor(label_Mij, levels = lbls_ordered))
  
  act_y <- y_df |> filter(E_order == actual_E, I_order == actual_I, 
                          dataset == ds) |> pull(y)
  
  posterior_mase_y <- extract_incidence(posterior_sample, 4000) |> 
    group_by(iter, M_i, M_j) |> summarise(mase = mase(act_y, value),
                                     .groups = "drop") |> 
    mutate(label_Mij = paste0("i = ", M_i,", ","j = ", M_j),
           data = "y",
           label_Mij = factor(label_Mij, levels = lbls_ordered))
  
  
  posterior_mase <- bind_rows(posterior_mase_x, posterior_mase_y)
  
  sub_txt <- paste0("'Histograms:' ~ M^{'ij'}")
  
  mean_df <- posterior_mase |> group_by(label_Mij, M_i, M_j, data) |> 
    summarise(mean_mase = mean(mase))
  
  ggplot(posterior_mase, aes(mase)) +
    scale_x_continuous(limits = lims) +
    scale_y_continuous(n.breaks = 3) +
    geom_histogram(aes(fill = as.factor(M_j)), alpha = 0.9) +
    scale_fill_manual(values = colours_df$colour) +
    facet_grid(label_Mij~data, scales = "free_x") +
    geom_vline(data = mean_df, aes(xintercept = mean_mase, 
                                   colour = as.factor(M_j))) +
    scale_colour_manual(values = colours_df$colour) +
    labs(x        = "Mean absolute scale error (MASE)", y = "",
         caption  = paste0("Dataset: ", ds, "\n", "Line: Mean value"),
         subtitle = parse(text = sub_txt)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.line  = element_line(colour = "grey80"),
          strip.background = element_rect(colour = "grey80"),
          axis.ticks = element_line(colour = "grey60"),
          axis.text  = element_text(colour = "grey60"),
          strip.text.y = element_text(angle = 0, vjust = 0.5))
}

plot_incidence <- function(y_list, dist) {
  
  y_all <- do.call(rbind, y_list) |> 
    mutate(label_Dj = paste0("j = ", I_order))
  
  disp_val <- ifelse(dist == "Poisson", "0", "1/3")
  
  ggplot(y_all, aes(time, y)) +
    geom_line(aes(group = dataset, colour = I_order), alpha = 0.10) +
    scale_colour_manual(values = colours_df$colour) +
    facet_wrap("label_Dj") +
    labs(x        = "Days", 
         y        = "Measured incidence [Cases/day]",
         subtitle = parse(text = paste("phi^{-1}~'=", disp_val, "'"))) +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_rect(colour = "grey80"),
          axis.title = element_text(colour = "grey40"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"))
}

plot_measured_incidence <- function(y_list, disp_val) {
  
  y_all <- do.call(rbind, y_list) |> 
    mutate(label_Dj = paste0("j = ", I_order))
  
  D_i <- unique(y_all$E_order)
  
  subtitle_txt <- paste0("'Lines:' ~ D^{'", D_i, "j'}")
  
  ggplot(y_all, aes(time, y)) +
    geom_line(aes(group = dataset, colour = I_order), alpha = 0.10) +
    scale_colour_manual(values = colours_df$colour) +
    facet_wrap("label_Dj") +
    labs(x        = "Days", 
         y        = "Incidence [Cases/day]",
         subtitle = parse(text = subtitle_txt),
         caption  = parse(text = paste("phi^-1:", disp_val))) +
    theme_pubr() +
    theme(legend.position = "none",
          axis.title = element_text(colour = "grey40"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"))
    
}

plot_x_prediction <- function(x_summary, y_subset) {
  
  ggplot(x_summary |> mutate(M_n = as.factor(M_n)), 
         aes(time, q50)) +
    geom_line(aes(group = M_n, colour = M_n)) +
    geom_point(data = y_subset, aes(time, x), size = 0.5, colour = "grey50") +
    facet_wrap(~ dataset) +
    theme_pubr()
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

summarise_R0_posterior <- function(posterior_df) {
  
  first_cond <- "par_inv_R0" %in% colnames(posterior_df)
  
  if(first_cond)posterior_df <- mutate(posterior_df, R0 = 1 / par_inv_R0)
  
  second_cond <- "par_gamma" %in% colnames(posterior_df)
  
  if(second_cond) posterior_df <- mutate(posterior_df, 
                                         R0 = par_beta / par_gamma)
  
  
  if(!first_cond & !second_cond) {
    
    posterior_df<- posterior_df |> mutate(R0 = par_beta * inv_gamma_val)
  }
  
  ds           <- posterior_df |> slice(1) |> pull(dataset)
  
  posterior_df |> group_by(M_i, M_j) |> 
    summarise(mean_R0     = mean(R0),
              lower_bound = quantile(R0, 0.025),
              upper_bound = quantile(R0, 0.975),
              .groups = "drop") |> 
    mutate(dataset = ds)
}

plot_R0_comparison_by_dataset <- function(fits_list, actual_val, x_min, x_max,
                                          D_i, D_j, text_size = 4) {
  
  summary_df <- map_df(fits_list, summarise_R0_posterior)
  
  M_i_size <- length(unique(summary_df$M_i))
  
  summary_df <- summary_df |> mutate(label_M_i = paste0("i = ", M_i),
                                     label_ds  = paste0("dataset = ", dataset)) 
  
  g <- ggplot(summary_df, aes(x = mean_R0, y = as.factor(M_j))) +
    scale_y_discrete(limits = rev) +
    geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound, 
                      colour = as.factor(M_j)), alpha = 0.75,
                  linetype = "dotted")
  
  if(M_i_size == 1) {
    g         <- g + facet_wrap(~label_ds)
    title_txt <- paste0("'Error bars:' ~ M^{'", M_i, "j'}")
  }
  
  if(M_i_size > 1) {
    g         <- g + facet_grid(label_ds ~ label_M_i)
    title_txt <- paste0("'Error bars:' ~ M^ij")
  }
      
  g + scale_colour_manual(values = colours_df$colour) +
    geom_text(aes(label = percent(mean_R0/actual_val - 1, accuracy = 1),
                  colour = as.factor(M_j)),
              size = text_size, show.legend = FALSE) +
    scale_x_continuous(limits = c(x_min, x_max)) +
    geom_vline(xintercept = actual_val, linetype = "dashed", 
               colour = colours_df[D_j, "colour"]) +
    labs(x = parse(text = "\u211c[0]"),
         y = "j", 
         title    = parse(text = title_txt),
         subtitle = parse(text = paste0("'Vertical line:' ~ D^{'", D_i, D_j,"'}")),
         colour   = "j") +
    theme_classic() +
    theme(legend.position = "none",
          text = element_text(family = "Arial Unicode MS"),
          axis.title.x = element_text(colour = "grey40"),
          axis.title.y  = element_text(colour = "grey40", angle = 0, 
                                       vjust = 0.5), 
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          strip.text = element_text(colour = "grey50"),
          strip.background = element_rect(colour = "grey80"))
}

plot_R0_posterior <- function(posterior_df, actual_val, D_i, D_j) {
  
  summary_df <- summarise_R0_posterior(posterior_df)
  
  g <- ggplot(summary_df, aes(x = mean_R0, y = as.factor(M_j))) +
    scale_x_continuous(limits = c(2, 3), breaks = c(2, 2.5, 3)) +
    geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound, 
                      colour = as.factor(M_j)), size = 0.25, width = 0.4) +
    scale_colour_manual(values = colours_df$colour) +
    geom_vline(xintercept = actual_val, linetype = "dotted",
               colour = colours_df[D_j, "colour"], size = 0.3) +
    labs(x = bquote(R[0]),
         y = "j",
         subtitle = parse(text = paste0("'Vertical line:' ~ D^{'", D_i, D_j,"'}"))) +
    theme_classic() +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = 6),
          axis.title.x  = element_text(colour = "grey40"),
          axis.title.y  = element_text(colour = "grey40", angle = 0, 
                                       vjust = 0.5), 
          axis.line     = element_line(colour = "grey80"),
          axis.text.x   = element_text(colour = "grey60", size = 6.5),
          axis.text.y   = element_text(colour = "grey60", size = 9),
          axis.ticks    = element_line(colour = "grey60"),
          plot.margin   = margin(0.1, 0.08, 0, 0.02, "cm"))
}

plot_par_estimates_by_dataset <- function(fits_list, var_name, actual_val, 
                                          x_min, x_max, x_label, D_i, D_j) {
  
  summary_df <- summarise_par(fits_list, var_name) |> 
    mutate(label_M_i = paste0("i = ", M_i),
           label_ds  = paste0("dataset = ", dataset)) 
    
  M_i_size   <- length(unique(summary_df$M_i))
  
  g <- ggplot(summary_df, aes(x = mean, y = as.factor(M_j))) +
    scale_y_discrete(limits = rev)
  
  if(M_i_size == 1) {
    g         <- g + facet_wrap(~label_ds)
    title_txt <- paste0("'Error bars:' ~ M^{'", M_i, "j'}")
  }
  
  if(M_i_size > 1) {
    g         <- g + facet_grid(label_ds ~ label_M_i)
    title_txt <- paste0("'Error bars:' ~ M^ij")
  }
  
    g +
      geom_vline(xintercept = actual_val, linetype = "dashed", 
                 colour = colours_df[D_j, "colour"]) +
    geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound, 
                      colour = as.factor(M_j)), alpha = 0.6,
                  position = position_dodge(width = 1),
                  linetype = "dotted") +
    scale_colour_manual(values = colours_df$colour) +
    geom_text(aes(label = percent(mean/actual_val - 1, accuracy = 1),
                  colour = as.factor(M_j)),
              position = position_dodge(width = 1.1),
              size = 4, show.legend = FALSE) +
    scale_x_continuous(limits = c(x_min, x_max)) +
    labs(x        = parse(text = x_label),
         y        = "j", 
         title    = parse(text = title_txt),
         subtitle = parse(text = paste0("'Vertical line:' ~ D^{'", D_i, D_j,"'}")),
         colour   = "j") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_text(colour = "grey40"),
          axis.title.y = element_text(colour = "grey40", 
                                      angle = 0, vjust = 0.5), 
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          strip.text = element_text(colour = "grey50"),
          strip.background = element_rect(colour = "grey80"))
}

summarise_par <- function(fits_list, var_name) {
  
  map_df(fits_list, \(posterior_df) {
    
    ds           <- posterior_df |> slice(1) |> pull(dataset)
    
    posterior_df |> group_by(M_i, M_j) |> 
      summarise(mean        = mean(.data[[var_name]]),
                lower_bound = quantile(.data[[var_name]], 0.025),
                upper_bound = quantile(.data[[var_name]], 0.975),
                .groups = "drop") |> 
      mutate(dataset = ds)
  }) -> df
  
  df
}

# ===================Pairs======================================================

dens_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) +
    stat_density2d(aes(fill=..density..), geom = "tile", contour = FALSE) +
    scale_fill_viridis_c()
  p
}

cor_fun <- function(data, mapping, method = "pearson", ndp = 2, sz = 5, 
                    text_size = 3, stars = TRUE, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor.test(x, y, method=method)
  est <- corr$estimate
  lb.size <- sz* abs(est)
  
  palette <- gradient_n_pal(c("lightgrey", "black"))
  
  if(stars){
    stars <- c("***", "**", "*", "")[findInterval(corr$p.value,
                                                  c(0, 0.001, 0.01, 0.05, 1))]
    
    lbl  <- paste0(round(est, ndp), stars)
    cor_colour <-  palette(abs(est))
  }else{
    lbl <- round(est, ndp)
  }
  
  ggplot(data = data, mapping = mapping) +
    annotate("text", x = mean(x, na.rm = TRUE), y = mean(y, na.rm=TRUE),
             label=lbl, size = text_size, colour = cor_colour,...)+
    theme(panel.grid = element_blank())
}

pairs_posterior <- function(posterior, strip_text = 3, axis_text_size = 4,
                            bins = 30, text_size = 3, M_j) {
  
  hist_col <- colours_df$colour[[M_j]]
  
  col_names           <- colnames(posterior)
  col_names           <- str_replace_all(col_names, "par_", "")
  colnames(posterior) <- col_names
  
  if("I0" %in% colnames(posterior)) {
    posterior <- rename(posterior, `I[0]` = I0)
  }
  
  if("R0" %in% colnames(posterior)) {
    posterior <- rename(posterior, `R(0)` = R0)
  }
  
  if("inv_phi" %in% colnames(posterior)) {
    posterior <- rename(posterior, `phi^{-1}` = inv_phi)
  }
  
  if("inv_R0" %in% colnames(posterior)) {
    posterior <- rename(posterior, '\u211c^{-1}' = inv_R0)
  }
  
  ggpairs(posterior, 
          lower = list(continuous = dens_fn),
          diag = list(continuous = wrap("barDiag", colour = hist_col, 
                                        fill = "white", bins = bins,
                                        size = 0.3)), 
          upper = list(continuous = wrap(cor_fun, text_size = text_size)),
          labeller = label_parsed) +
    theme_pubr() +
    theme(text = element_text(family = "Arial Unicode MS"),
          axis.text = element_text(size = axis_text_size),
          strip.text = element_text(size = strip_text))
}

#===============================================================================

plot_incidence_prediction <- function(incidence_df, y_df, D_i, D_j, ds, 
                                      n_col      = 2, 
                                      y_lims     = NULL,
                                      line_size  = 1,
                                      point_size = 1,
                                      axis_labels_size = 10,
                                      axis_x_size = 11,
                                      y_lab = "Incidence [Cases/day]",
                                      strip_text_size = 10,
                                      sub_size = 11,
                                      caption = TRUE,
                                      line_alpha = 0.05) {
  
  epi_df <- y_df |> filter(I_order == D_j, dataset == ds)
  
  sorted_i <- sort(unique(incidence_df$M_i), decreasing = TRUE)
  sorted_j <- sort(unique(incidence_df$M_j), decreasing = TRUE)
  
  ordered_ij <- expand.grid(sorted_i, sorted_j) |>
    mutate(ij = paste0(Var1, Var2)) |> pull(ij)
  
  incidence_df <- incidence_df |> mutate(fit_model = paste0("M^", M_i, M_j))
  
  incidence_df$fit_model <- factor(incidence_df$fit_model, 
                                   levels = paste0("M^", rev(ordered_ij)))
  
  g <- ggplot(incidence_df, aes(time, value)) +
    geom_line(aes(group = iter, colour = as.factor(M_j)), alpha = line_alpha,
              size = line_size) +
    facet_wrap(~fit_model, labeller = label_parsed, ncol = n_col)
  
  if(!is.null(y_lims)) g <- g + scale_y_continuous(limits = y_lims)
  
  g<- g +
    scale_colour_manual(values = colours_df$colour) +
    geom_point(data = epi_df, aes(y = x), size = point_size, 
               colour = colours_df[D_j, "colour"]) +
    labs(x        = "Days", 
         y        = y_lab,
        subtitle = parse(text = str_glue("Points: D^{D_i}{D_j}")))
  
  if(caption) g +  labs(caption  = str_glue("Dataset: {ds}"))
    
    
  g + theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_text(colour = "grey40", size = axis_x_size),
          axis.title.y = element_text(colour = "grey40"),
          axis.line    = element_line(colour = "grey80"),
          axis.text    = element_text(colour = "grey60", 
                                      size = axis_labels_size),
          axis.ticks = element_line(colour = "grey60"),
          strip.background = element_rect(colour = "grey80"),
          strip.text = element_text(margin = margin(0.1, 0, 0.1, 0, "cm"),
                                    size   = strip_text_size),
          plot.subtitle = element_text(size = sub_size))
}

plot_incidence_prediction_by_E_dist <- function(incidence_df, y_df, D_i, D_j, ds) {
  
  incidence_df <- incidence_df |> 
    mutate(fit_model = paste0("M^", M_i, M_j))
  
  M_j <- unique(incidence_df$M_j)
  
  time_points <- c(1, seq(4, 20, 5), seq(22, 40, 2), seq(45, 55, 5))
  
  epi_df <- y_df |> filter(E_order == D_i,
                           I_order == D_j,
                           time %in% time_points)
  
  ggplot(incidence_df, aes(time, value)) +
    geom_line(aes(group = iter, linetype = as.factor(M_i)), 
              colour = colours_df[M_j, "colour"],
              alpha = 0.25, size = 0.25) +
    geom_point(data = epi_df, aes(y = x), size = 1.25, 
               colour = colours_df[D_j, "colour"]) + 
    scale_linetype_manual(values = c("solid", "twodash")) +
    facet_wrap(~fit_model, labeller = label_parsed, ncol = 2) +
    labs(x        = "Days", 
         y        = "x [Cases / day]",
         subtitle = parse(text = str_glue("Points: D^{D_i}{D_j}"))) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_text(colour = "grey40", size = 9),
          axis.title.y = element_text(colour = "grey40", size = 9),
          axis.line    = element_line(colour = "grey80"),
          axis.text    = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          strip.background = element_rect(colour = "grey80"),
          strip.text = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")))
}

plot_incidence_prediction_by_fitting_model <- function(posterior_supralist, 
                                                       y_df, ds, M_i) {
  
  incidence_df <- map_df(posterior_supralist, \(posterior_list) {
    
    posterior_sample <- posterior_list[[ds]]
    D_i              <- unique(posterior_sample$D_i)
    D_j              <- unique(posterior_sample$D_j)
    extract_incidence(posterior_sample, 20) |> 
      mutate(D_i = D_i, D_j = D_j)
  }) |> mutate(label_Dij = paste0("'Points:'~D^",D_i, D_j),
               label_Mij = paste0("'j = ", M_j,"'"))
  
  subtitle_txt <- paste0("'Lines:' ~ M^{'", M_i, "j'}")
  
  epi_df <- y_df |> filter(dataset == ds,
                           time %in% c(1, seq(4, 48, 4))) |> 
    mutate(label_Dij = paste0("'Points:'~D^", E_order, I_order))
  
  ggplot(incidence_df, aes(time, value)) +
    geom_line(aes(group = iter, colour = as.factor(M_j)), alpha = 0.1) +
    geom_point(data = epi_df, aes(y = x, colour = label_Dij), size = 1) +
    scale_colour_manual(values = c(colours_df$colour, colours_df$colour)) +
    facet_grid(label_Mij ~label_Dij, labeller = label_parsed) +
    labs(y = "x [Cases/day]", x = "Days",
         subtitle = parse(text = subtitle_txt)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_text(colour = "grey40"),
          axis.title.y = element_text(colour = "grey40"),
          axis.line    = element_line(colour = "grey80"),
          axis.text    = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          strip.background.x  = element_rect(colour = "white"),
          strip.background.y  = element_rect(colour = "grey80"),
          strip.text.x = element_text(hjust = 0, colour = "grey25"))
}

prior_posterior_comparison <- function(posterior_df, prior_df, par_name,
                                       sub_txt = "Grey histogram (prior): beta(2, 2)") {
  
  x_lab <- str_remove(par_name, "par_")
  
  M_i <- unique(posterior_df$M_i)
  title_txt <- title_txt <- paste0("'Coloured histograms (posterior): ' ~ M^{'", M_i, "j'}")
  
  ds <- unique(posterior_df$dataset)
  caption_txt <- str_glue("Dataset: {ds}")
  
  posterior_df <- posterior_df |> 
    mutate(label_j = paste0("j = ", M_j))
  
  ggplot(posterior_df, aes(.data[[par_name]])) +
    geom_histogram(data = prior_df, aes(sims), fill = "grey70", alpha = 0.5, 
                   colour = "white") +
    geom_histogram(aes(fill = as.factor(M_j)), colour = "white", alpha = 0.55) +
    scale_fill_manual(values = colours_df$colour) +
    labs(x = parse(text = x_lab),
         title    = parse(text = title_txt),
         subtitle = sub_txt,
         y = "Count",
         caption = caption_txt) +
    facet_wrap(~label_j) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(colour = "grey40"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          strip.background = element_rect(colour = "grey80"),
          plot.title = element_text(colour = "grey40"),
          plot.subtitle = element_text(colour = "grey40"),
          plot.caption  = element_text(colour = "grey40"))
}

plot_R0_vs_tau <- function(df, ds, x_lims, y_lims) {
  
  df <- df |> filter(dataset == ds, D_j == 2) |> 
    mutate(label_Dij = paste0("D^", D_i, D_j))
  
  M_i          <- unique(df$M_i)
  subtitle_txt <- title_txt <- paste0("'Points:' ~ M^{'", M_i, "j'}")
  
  ggplot(df, aes(tau, R0)) +
    geom_point(aes(group = id, colour = as.factor(M_j)),
               fill = "white", size = 1, shape = 1, alpha = 0.5) +
    scale_colour_manual(values = colours_df$colour, name = "j") +
    facet_wrap(~label_Dij, labeller = label_parsed, ncol = 4) +
    geom_vline(aes(xintercept = mean_generation_time(D_j, 2, 2)),
               colour = "grey50", linetype = "dashed") +
    geom_hline(yintercept = R0_val, colour = "grey50", linetype = "dotted") +
    scale_x_continuous(limits = x_lims, breaks = c(0, 5, 10)) +
    scale_y_continuous(limits = y_lims) +
    labs(x = bquote(tau), y = parse(text = "\u211c[0]"),
         subtitle = parse(text = subtitle_txt)) +
    guides(alpha = "none") +
    theme_classic() +
    theme(text = element_text(family = "Arial Unicode MS"),
          legend.position = "none",
          axis.line  = element_line(colour = "grey80"),
          axis.ticks = element_line(colour = "grey60"),
          axis.title.x  = element_text(colour = "grey40"),
          axis.title.y  = element_text(colour = "grey40", angle = 0, 
                                       vjust = 0.5),
          strip.background = element_rect(colour = "grey95"),
          axis.text.x   = element_text(colour = "grey60"),
          axis.text.y   = element_text(colour = "grey60"))
}

plot_scatterplot_R0_vs_tau <- function(df, x_lims, y_lims) {
  
  df <- df |> mutate(label_Dij = paste0("D^", D_i, D_j),
                     label_ds  = paste0("'dataset = ", dataset,"'"))
  
  M_i <- unique(df$M_i)
  
  subtitle_txt <- title_txt <- paste0("'Points:' ~ M^{'", M_i, "j'}")
  
  ggplot(df, aes(tau, R0)) +
    geom_point(aes(group = id, colour = as.factor(M_j)),
               fill = "white", size = 1, shape = 1) +
    scale_colour_manual(values = colours_df$colour, name = "j") +
    facet_grid(label_ds~label_Dij, labeller = label_parsed) +
    geom_vline(aes(xintercept = mean_generation_time(D_j, 2, 2)),
               colour = "grey50", linetype = "dashed") +
    geom_hline(yintercept = R0_val, colour = "grey50", linetype = "dotted") +
    scale_x_continuous(limits = x_lims) +
    scale_y_continuous(limits = y_lims) +
    labs(x = bquote(tau), 
         y = parse(text = "\u211c[0]"),
         caption = paste("Dashed line: Data's generation time",
                         "Dotted line: Data's basic reproduction number",
                         sep = "\n"),
         subtitle = parse(text = subtitle_txt)) +
    guides(alpha = FALSE) +
    theme_classic() +
    theme(text = element_text(family = "Arial Unicode MS"),
          legend.position = "bottom",
          plot.caption = element_text(hjust = 0),
          axis.title.x = element_text(colour = "grey40"),
          axis.title.y = element_text(colour = "grey40", angle = 0, 
                                       vjust = 0.5),
          axis.text.x   = element_text(colour = "grey60"),
          axis.text.y   = element_text(colour = "grey60"),
          axis.line  = element_line(colour = "grey80"),
          axis.ticks = element_line(colour = "grey60"))
}

plot_R0_vs_tau_by_fitting_model <- function(df, ds, x_lims, y_lims) {
  
  subset_df <- df |> filter(dataset == ds) |> 
    mutate(label_Mj = paste0("'j = ", M_j, "'"),
           label_Dij = paste0("D^", D_i, D_j)) |> 
    group_by(M_j, D_j) |> 
    slice_sample(n = 250)
  
  subtitle_txt <- title_txt <- paste0("'Points:' ~ M^{'", M_i, "j'}")
  
  ggplot(subset_df, aes(tau, R0)) +
    geom_point(aes(group = id, colour = as.factor(M_j)),
               fill = "white", size = 1.5, shape = 1, alpha = 0.4) +
    scale_colour_manual(values = colours_df$colour, name = "j") +
    facet_grid(label_Mj~label_Dij, labeller = label_parsed) +
    geom_vline(aes(xintercept = mean_generation_time(D_j, 2, 2)),
               colour = "grey50", linetype = "dashed") +
    geom_hline(yintercept = R0_val, colour = "grey50", linetype = "dotted") +
    scale_x_continuous(limits = x_lims, breaks = c(0, 5, 10)) +
    scale_y_continuous(limits = y_lims, breaks = c(0, 3, 6)) +
    labs(x = parse(text = "'Mean generation time ('~tau~')'"), 
         y = parse(text = "\u211c[0]"),
         caption = paste("Dashed line: Data's generation time",
                         "Dotted line: Data's basic reproduction number",
                         sep = "\n"),
         subtitle = parse(text = subtitle_txt)) +
    theme_classic() +
    theme(text = element_text(family = "Arial Unicode MS"),
          legend.position = "none",
          plot.caption = element_text(hjust = 0, size = 6, colour = "grey35"),
          axis.title.x  = element_text(colour = "grey40"),
          axis.title.y  = element_text(colour = "grey40", angle = 0, 
                                       vjust = 0.5), 
          axis.line     = element_line(colour = "grey80"),
          axis.text.x   = element_text(colour = "grey60", size = 6.5),
          axis.text.y   = element_text(colour = "grey60", size = 9),
          axis.ticks    = element_line(colour = "grey60"),
          strip.background = element_rect(colour = "grey95"))
}

plot_R0_vs_Mj <- function(df, R0_val, inv_sigma, inv_gamma, y_lims) {
  
  M_i <- unique(df$M_i)
  
  subtitle_txt <- title_txt <- paste0("'Violin plots:' ~ M^{'", M_i, "j'}")
  
  df <- df |>
    mutate(label_Dij = paste0("D^", D_i, D_j, "~'('~tau~'= ",
                              round(mean_generation_time(D_j, inv_sigma,
                                                         inv_gamma), 1),
                              ")'"))
  
  ggplot(df, aes(as.factor(M_j), R0)) +
    geom_violin(aes(group = as.factor(M_j), fill = as.factor(M_j)),
                colour = "white", alpha = 0.8) +
    scale_fill_manual(values = colours_df$colour) +
    geom_hline(yintercept = R0_val, colour = "grey60", linetype = "dotted",
               alpha = 0.75) +
    scale_y_continuous(limits = y_lims) +
    labs(y = bquote(R[0]), x = "j",
         subtitle = parse(text = subtitle_txt)) +
    facet_grid(dataset~label_Dij, labeller = label_parsed) +
    theme_pubr() +
    theme(legend.position = "none",
          axis.title = element_text(colour = "grey40"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          strip.background = element_rect(colour = "white"))
}

plot_R0_by_fitting_model <- function(df, ds, actual_val, n_row = 1, 
                                     lims, breaks = waiver(),
                                     axis_text_x = 8,
                                     strip_text_size = 6, bar_width = 0.8) {
  
  summary_df <- df |> filter(dataset == ds) |> 
    mutate(D_ij = paste0(D_i, D_j)) |> group_by(D_ij, M_i, M_j) |> 
    summarise(mean_R0     = mean(R0),
              lower_bound = quantile(R0, 0.025),
              upper_bound = quantile(R0, 0.975),
              .groups = "drop") |> 
    mutate(label_Dij = paste0("'Vertical line:'~D^", D_ij))
  
  title_txt <- paste0("'Error bars:' ~ M^{'", 1, "j'}")
  
  vline_df <- data.frame(actual_val = actual_val,
                         label_Dij  = unique(summary_df$label_Dij))
  
  
  ggplot(summary_df, aes(x = mean_R0, y = as.factor(M_j))) +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(breaks = breaks, limits = lims) +
    geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound, 
                      colour = as.factor(M_j)), alpha = 0.75, 
                  width = bar_width) + 
    facet_wrap(~label_Dij, nrow = n_row, labeller = label_parsed) +
    geom_vline(data = vline_df, aes(xintercept = actual_val, 
                                    colour = as.factor(label_Dij)), 
               linetype = "dashed") +
    scale_colour_manual(values = c(colours_df$colour, colours_df$colour)) +
    labs(x = parse(text = "\u211c[0]"),
         y = "j", 
         subtitle    = parse(text = title_txt),
         colour   = "j") +
    theme_classic() +
    theme(text = element_text(family = "Arial Unicode MS"),
          legend.position = "none",
          axis.title.x = element_text(colour = "grey40"),
          axis.title.y = element_text(colour = "grey40", angle = 0, 
                                       vjust = 0.5),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60", size = axis_text_x),
          axis.ticks = element_line(colour = "grey60"),
          strip.background = element_rect(colour = "white"),
          strip.text = element_text(hjust = 0, colour = "grey35",
                                    size = strip_text_size ))
}

plot_vars_by_fitting_model <- function(df, ds, var_names) {
  
  summary_df <- summarise_vars(df, ds, var_names)
  
  title_txt <- paste0("'Error bars:' ~ M^ij")
  
    ggplot(summary_df, aes(x = mean, y = as.factor(M_j))) +
      scale_y_discrete(limits = rev) +
      geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound, 
                        colour = as.factor(M_j)), alpha = 0.75) +
      scale_colour_manual(values = c(colours_df$colour)) +
      facet_grid(label_Dij ~ var_name, labeller = label_parsed, 
                 scales = "free_x") +
      labs(x = "Value",
           y = "j", 
           subtitle    = parse(text = title_txt),
           colour   = "j") +
      theme_classic() +
      theme(legend.position = "none",
            axis.title.x = element_text(colour = "grey40"),
            axis.title.y = element_text(colour = "grey40", angle = 0, 
                                        vjust = 0.5),
            axis.line  = element_line(colour = "grey80"),
            axis.text  = element_text(colour = "grey60", size = 8),
            axis.ticks = element_line(colour = "grey60"),
            strip.background = element_rect(colour = "grey90"),
            strip.text = element_text(colour = "grey35", size = 6))
}

summarise_vars <- function(df, ds, var_names) {
  
  df <- df |> filter(dataset == ds) |> 
    mutate(D_ij = paste0(D_i, D_j)) 
  
  map_df(var_names, \(var_name) {
    
    var_df <- df |> select(D_ij, M_i, M_j, all_of(var_name)) |> 
      rename(var_name = !!var_name) |> 
      group_by(D_ij, M_i, M_j) |> 
      summarise(mean = mean(var_name),
                lower_bound = quantile(var_name, 0.025),
                upper_bound = quantile(var_name, 0.975),
                .groups = "drop") |> 
      mutate(label_Dij = paste0("'Vertical line:'~D^", D_ij),
             var_name = var_name)
  })
}

plot_prior <- function(df, var_name, dist, lims, bins) {
  
  ggplot(prior_df, aes(x)) +
    geom_histogram(colour = "white", bins = bins, fill = "grey75") +
    scale_x_continuous(limits = lims) +
    labs(y = "", x = parse(text = var_name),
         subtitle = parse(text = dist)) +
    theme_classic() +
    theme(text = element_text(family = "Arial Unicode MS"),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          axis.line.x  = element_line(colour = "grey85"))
}

plot_MLE_summary <- function(df) {
  
  df <- df |> 
    mutate(label_Dij = paste0("D^",D_i, D_j),
           match = ifelse(D_j == M_j, TRUE, FALSE))
  
  ggplot(df, aes(x = M_j)) +
    geom_bar(aes(fill = as.factor(M_j), group = M_j)) +
    geom_text(aes(label = ..count.., colour = as.factor(M_j)), stat = "count",
              nudge_y = 1.5) +
    scale_y_continuous(limits = c(0, 20)) +
    scale_fill_manual(values = colours_df$colour) +
    scale_colour_manual(values = colours_df$colour) +
    facet_wrap(~label_Dij, labeller = label_parsed) +
    scale_x_continuous(breaks = 1:9) +
    labs(y = "Number of hits", x = "j",
         subtitle = parse(text = paste0("'Bars: ' ~ M^{'", M_i, "j'}"))) +
    theme_classic() +
    theme(legend.position = "none",
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          strip.text = element_text(colour = "grey20"),
          strip.background = element_rect(colour = "grey80"))
}