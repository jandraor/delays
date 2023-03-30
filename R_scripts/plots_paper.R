fig_3A <- function(df, ds, actual_val, n_row = 1, lims, axis_text_x = 8,
                   strip_text_size = 6, bar_width = 0.8) {
  
  
    
    summary_df <- summarise_R0_for_figs(df, ds) |> 
      mutate(label_Dij = paste0("'Vertical line:'~D^", D_ij))

    vline_df <- data.frame(actual_val = actual_val,
                           label_Dij  = unique(summary_df$label_Dij),
                           D_j        = unique(summary_df$D_j))

    
    label_df <- data.frame(label_Dij = rep(str_glue("'Vertical line:'~D^1{1:4}"),each = 4),
                           M_j = 1:4, x = c(3   , 2.9, 2.85, 2.8,
                                            3   , 2.9, 2.85, 2.8,
                                            3.05, 2.95, 2.90, 2.85,
                                            3.05, 2.95, 2.90, 2.85), y = 4:1) |> 
      mutate(txt = str_glue("M^1{M_j}"))
    
    plot_error_bar(summary_df,
                   limits   = c(2.1, 3.1),
                   breaks   = c(2.2, 2.6, 3.0),
                   vline_df = vline_df,
                   label_df = label_df)
}


fig_4A <- function(df, ds, actual_val, n_row = 1, lims, axis_text_x = 8,
                   strip_text_size = 6, bar_width = 0.8) {
  
  summary_df <- summarise_R0_for_figs(df, ds) |> 
    mutate(label_Dij = paste0("'Vertical line:'~D^", D_ij))

  vline_df <- data.frame(actual_val = actual_val,
                         label_Dij  = unique(summary_df$label_Dij),
                         D_j        = unique(summary_df$D_j))
  

  label_df <- data.frame(label_Dij = rep(str_glue("'Vertical line:'~D^1{1:4}"),each = 4),
                         M_j = 1:4, x = c(1.50, 1.35, 1.30, 1.25,
                                          1.75, 1.65, 1.6, 1.55,
                                          1.75, 1.65, 1.6, 1.55,
                                          1.75, 1.65, 1.6, 1.55), y = 4:1) |> 
    mutate(txt = str_glue("M^1{M_j}"))
  
  
  plot_error_bar(summary_df,
                 limits   = c(1, 6),
                 breaks   = c(2, 4, 6),
                 vline_df = vline_df,
                 label_df = label_df)
  

}

plot_error_bar <- function(df, limits, breaks, vline_df, label_df,
                           legend_pos = "top") {
  
  title_txt <- paste0("'Error bars: Estimates from' ~ M^{'", 1, "j'}~'fitting instances'")
  
  cap_txt <- paste("Vertical line: Actual value \t \t j: # of stages in the infectious class (I)", 
                   "D\U2071\U02B2: Data generator's distribution \t M\U2071\U02B2: Fitting model's distribution",
                   sep = "\n")
  
  ggplot(df, aes(x = mean_R0, y = as.factor(M_j))) +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(breaks = breaks, limits = limits) +
    geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound, 
                      colour = as.factor(M_j)), alpha = 0.6, 
                  width = 0.8) +
    geom_text(data = label_df, aes(x = x, y = y, label = txt,
                                   colour = as.factor(M_j)), parse = TRUE,
              size = 2.5, show.legend = FALSE, alpha = 0.65) +
    facet_wrap(~label_Dij, nrow = 1, labeller = label_parsed) +
    geom_vline(data = vline_df, aes(xintercept = actual_val, 
                                    colour = as.factor(D_j)),
               show.legend = FALSE, linetype = "dashed") +
    scale_colour_manual(values = c(colours_df$colour, colours_df$colour)) +
    labs(x = parse(text = "'Basic reproduction number ('~\u211c[0]~')'"),
         y = "", 
         subtitle    = parse(text = title_txt),
         colour   = "j",
         caption = cap_txt) +
    theme_classic() +
    theme(text = element_text(family = "Arial Unicode MS"),
          legend.position = legend_pos,
          legend.margin = margin(0, 0, 0, 0, "cm"),
          axis.title.x = element_text(colour = "grey65"),
          axis.line.x  = element_line(colour = "grey80"),
          axis.line.y  = element_blank(),
          axis.text.x  = element_text(colour = "grey70", size = 8),
          axis.text.y  = element_blank(),
          axis.ticks.x = element_line(colour = "grey70"),
          axis.ticks.y = element_blank(),
          strip.background = element_rect(colour = "white"),
          strip.text = element_text(hjust = 0, colour = "grey35",
                                    size = 6),
          plot.caption = element_text(hjust = 0, colour = "grey65", size = 8.5),
          plot.subtitle = element_text(colour = "grey55", size = 9))
}


summarise_R0_for_figs <- function(df, ds) {
  
  summary_df <- df |> filter(dataset == ds) |> 
    mutate(D_ij = paste0(D_i, D_j)) |> group_by(D_ij, D_j, M_i, M_j) |> 
    summarise(mean_R0     = mean(R0),
              lower_bound = quantile(R0, 0.025),
              upper_bound = quantile(R0, 0.975),
              .groups = "drop") 
}


fig_4B <- function(df, ds, x_lims, y_lims) {
  
  subset_df <- df |> filter(dataset == ds) |> 
    mutate(label_Mj = paste0("M^", M_i, M_j),
           label_Dij = paste0("D^", D_i, D_j)) |> 
    group_by(M_j, D_j) |> 
    slice_sample(n = 250)
  
  vline_df <- data.frame(label_Dij = unique(subset_df$label_Dij),
                         D_j       = unique(subset_df$D_j)) |> 
    mutate(x_intr = mean_generation_time(D_j, 2, 2))
  
  h_line_df <- data.frame(M_i = 1,
                          M_j = 1:4,
                          y_intr = R0_val) |> 
    mutate(label_Dij = unique(subset_df$label_Dij))
  
  sub_txt <- paste0("'Points: Estimates from' ~ M^{'", 1, "j'}~'fitting instances'")
  
  ggplot(subset_df, aes(tau, R0)) +
    geom_point(aes(group = id, colour = as.factor(M_j)),
               fill = "white", size = 1.5, shape = 1, alpha = 0.4) +
    scale_colour_manual(values = c(colours_df$colour,
                                   colours_df$colour), name = "j") +
    facet_grid(label_Mj~label_Dij, labeller = label_parsed) +
    geom_vline(data = vline_df, aes(xintercept = x_intr, colour = label_Dij),
               linetype = "dashed", alpha = 0.75, linewidth = 0.25) +
    geom_hline(data = h_line_df, aes(yintercept = y_intr, 
                                     colour     = label_Dij),  
               linetype = "dotted", alpha = 0.75, linewidth = 0.5) +
    scale_x_continuous(limits = x_lims, breaks = c(0, 5, 10)) +
    scale_y_continuous(limits = y_lims, breaks = c(0, 3, 6)) +
    labs(x = parse(text = "'Mean generation time ('~tau~')'"), 
         y = parse(text = "'Basic reproduction number ('~\u211c[0]~')'"),
         caption = paste("Dashed line: Data's generation time",
                         "Dotted line: Data's basic reproduction number",
                         sep = "\n"),
         subtitle = parse(text = sub_txt)) +
    theme_classic() +
    theme(text = element_text(family = "Arial Unicode MS"),
          legend.position = "none",
          plot.caption = element_text(hjust = 0, size = 6, colour = "grey65"),
          plot.subtitle = element_text(colour = "grey55", size = 9),
          axis.title.x  = element_text(colour = "grey65"),
          axis.title.y  = element_text(colour = "grey65", size = 8), 
          axis.line     = element_line(colour = "grey90"),
          axis.text.x   = element_text(colour = "grey70", size = 6.5),
          axis.text.y   = element_text(colour = "grey70", size = 9),
          axis.ticks    = element_line(colour = "grey90"),
          strip.background = element_rect(colour = "grey95"),
          strip.text.y  = element_text(angle = 0))
}

fig_4C <- function(df, ds, x_lims, y_lims) {
  
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
               colour = colours_df$colour[[2]], linetype = "dashed") +
    geom_hline(yintercept = R0_val, colour = colours_df$colour[[2]],
               linetype = "dotted") +
    scale_x_continuous(limits = x_lims, breaks = c(0, 5, 10)) +
    scale_y_continuous(limits = y_lims) +
    labs(x = bquote(tau), y = parse(text = "\u211c[0]"),
         subtitle = parse(text = subtitle_txt)) +
    guides(alpha = "none") +
    theme_classic() +
    theme(text = element_text(family = "Arial Unicode MS"),
          legend.position = "none",
          axis.line  = element_line(colour = "grey90"),
          axis.ticks = element_line(colour = "grey90"),
          axis.title.x  = element_text(colour = "grey65"),
          axis.title.y  = element_text(colour = "grey65", angle = 0, 
                                       vjust = 0.5),
          strip.background = element_rect(colour = "grey95"),
          axis.text.x   = element_text(colour = "grey70"),
          axis.text.y   = element_text(colour = "grey70"))
}

fig_5 <- function(df, ds, actual_val) {
  
  summary_df <- summarise_R0_for_figs(df, ds) |> 
    mutate(label_Dij = paste0("'Vertical line:'~D^", D_ij))
  
  vline_df <- data.frame(actual_val = actual_val,
                         label_Dij  = unique(summary_df$label_Dij),
                         D_j        = unique(summary_df$D_j))
  
  label_df <- data.frame(label_Dij = rep(str_glue("'Vertical line:'~D^1{1:4}"),each = 4),
                         M_j = 1:4, x = 2.42 , y = 4:1) |> 
    mutate(txt = str_glue("M^1{M_j}"))
  
  plot_error_bar(summary_df,
                 limits     = c(2.4, 2.6),
                 breaks     = c(2.4, 2.5, 2.6),
                 vline_df   = vline_df,
                 label_df   = label_df,
                 legend_pos = "none")
}

fig_6A <- function(incidence_df, y_df, D_i, D_j, ds) {
  
  incidence_df <- incidence_df |> 
    mutate(fit_model = paste0("M^", M_i, M_j),
           dist      = ifelse(M_i == D_i, "right", "wrong"),
           fit_model = str_glue("M^{M_i}{M_j}~'({dist} latent period distribution)'"))
  
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
         y        = "Latent incidence (x)",
         title    = parse(text = str_glue("'Lines : Simulations from two'~ M^ij~'instances fitted to one D\U00B3\U00B3 dataset'")),
         subtitle = parse(text = str_glue("'Points: Dataset obtained from SE\U00B3I\U00B3R (D\U00B3\U00B3)'"))) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_text(colour = "grey40", size = 9),
          axis.title.y = element_text(colour = "grey40", size = 9),
          axis.line    = element_line(colour = "grey80"),
          axis.text    = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          plot.title    = element_text(colour = "grey65", size = 9,
                                       margin = margin(0.0, 0, 0.05, 0, "cm")),
          plot.subtitle = element_text(colour = "grey65", size = 9,
                                       margin = margin(0., 0, 0., 0, "cm")),
          strip.background = element_rect(colour = "grey80"),
          strip.text = element_text(colour = "grey30",
                                    size = 8,
                                    margin = margin(0.1, 0, 0.1, 0, "cm")))
}

fig_6B <- function(summary_df, R0_val) {
  
  label_df <- data.frame(M_i  = rep(c(rep(1, 4), rep(3, 4)), 2),
                         M_j  = rep(1:4, 4),
                         name = c(rep("I[0]", 8), rep("\u211c[0]", 8)),
                         x    = c(rep(1.3, 4), rep(2.1, 4), 
                                  rep(c(3.13, 2.83, 2.72, 2.68), 2))) |> 
    mutate(txt = str_glue("M^{M_i}{M_j}")) |> 
    mutate(M_j = as.factor(M_j))
  
  vline_df <- data.frame(name = c("I[0]", "\u211c[0]"), val = c(1, R0_val))
                         
  
  ggplot(summary_df, aes(x = mean, y = as.factor(M_j))) +
    geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound, 
                      colour = as.factor(M_j), linetype = as.factor(M_i)),
                  width = 0.8) +
    scale_y_discrete(limits = rev) +
    scale_colour_manual(values = colours_df$colour) +
    scale_linetype_manual(values = c("solid", "twodash")) +
    geom_vline(data = vline_df, aes(xintercept = val), linetype = "twodash", 
               colour = colours_df[D_j, "colour"]) +
    geom_text(data = label_df, aes(x = x, y = M_j, label = txt, 
                                   colour = as.factor(M_j)), parse = TRUE,
              size = 2.5) +
    facet_grid(M_i~name, scales = "free", labeller = label_parsed) +
    geom_blank(data = data.frame(name = c("I[0]", "\u211c[0]"), 
                                 x = c(2.3, 3.2), y = 1), aes(x = x, y = y)) + 
    labs(x = "Value", y = "",
         title = parse(text = str_glue("'Error bars: Estimates from '~ M^ij~'instances fitted to one D\U00B3\U00B3 dataset'")),
         subtitle = parse(text = paste0("'Vertical lines: Values used to configure the SE\U00B3I\U00B3R (D\U00B3\U00B3) data generator'"))) +
    theme_classic() +
    theme(text = element_text(family = "Arial Unicode MS"),
          legend.position = "none",
          strip.background.y = element_blank(),
          strip.text.y       = element_blank(),
          strip.background.x = element_rect(colour = "grey80"),
          axis.text.x   = element_text(colour = "grey60"),
          axis.text.y   = element_blank(),
          axis.ticks.x  = element_line(colour = "grey60"),
          axis.ticks.y  = element_blank(),
          axis.line.x   = element_line(colour = "grey80"),
          axis.line.y   = element_blank(),
          plot.title    = element_text(colour = "grey65", size = 9,
                                       margin = margin(0.0, 0, 0.05, 0, "cm")),
          plot.subtitle = element_text(colour = "grey65", size = 9,
                                       margin = margin(0,0,0,0, "cm")),
          axis.title.y = element_text(colour = "grey40", angle = 0, 
                                      vjust = 0.5),
          axis.title.x = element_text(colour = "grey40"))
}

fig_7_main <- function(incidence_df, data_df) {
  
  incidence_summary <- incidence_df |> group_by(time, M_j, parameterisation) |> 
    summarise(mean = mean(y),
              lb   = quantile(y, 0.025),
              ub   = quantile(y, 0.975),
              .groups = "drop") |> 
    filter(M_j %in% c(1, 4))
  
  title_txt <-"'Line, ribbon, error bars: Estimates from '~ M^{1*j}~'instances fitted to Cumberland\\'s incidence'"
  
  label_df <- data.frame(M_j = c(1, 4),
                         x   = 70,
                         y   = 80) |> 
    mutate(txt = str_glue("M^1{M_j}"))
  
  ggplot(incidence_summary, aes(time, mean)) +
    geom_ribbon(data = incidence_summary, alpha = 0.1,
                aes(y = mean, ymin = lb, ymax = ub, fill = as.factor(M_j))) +
    geom_line(aes(colour = as.factor(M_j))) +
    geom_text(data = label_df, aes(label = txt, x = x, y = y, 
                                   colour = as.factor(M_j)), parse = TRUE) +
    facet_grid(parameterisation~M_j) +
    scale_colour_manual(values = colours_df$colour[c(1, 4)]) +
    scale_fill_manual(values = colours_df$colour[c(1, 4)]) +
    geom_point(data = data_df, aes(y = y), size = 1, colour = Cmb_colour,
               shape = 18) +
    labs(x        = "Days", 
         y        = "Incidence [Cases/day]",
         title    = parse(text = title_txt),
         subtitle = "Points: Cumberland's data",
         caption  = "Fitted model instances stem from the alternative parameterisation") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(colour = "grey40"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          plot.title    = element_text(colour = "grey65", size = 8,
                                       margin = margin(0.0, 0, 0.05, 0, "cm")),
          plot.caption = element_text(colour = "grey60"),
          plot.subtitle = element_text(colour = "grey65", size = 8,
                                       margin = margin(0, 0, 0.05, 0, "cm")))
}

fig_7_inset <- function(summary_df) {
  
  label_df <- data.frame(x = 1.92,
                         M_j = as.factor(1:4)) |> 
    mutate(txt = str_glue("M^1{M_j}"))
  
  ggplot(summary_df, aes(x = mean_R0, y = as.factor(M_j))) +
    scale_x_continuous(limits = c(1.8, 2.2), n.breaks = 3) +
    scale_y_discrete(limits = rev) +
    geom_errorbar(aes(xmin = lower_bound, xmax = upper_bound, 
                      colour = as.factor(M_j)), alpha = 0.75, 
                  width = 0.5) +
    geom_text(data = label_df, aes(x = x, y = M_j, colour = as.factor(M_j),
                                   label = txt), size = 2.2,
              parse = TRUE) + 
    scale_colour_manual(values = c(colours_df$colour, colours_df$colour)) +
    labs(x = parse(text = "\u211c[0]"),
         y = "") +
    theme_classic() +
    theme(text = element_text(family = "Arial Unicode MS"),
          legend.position = "none",
          axis.text.x  = element_text(colour = "grey25", size = 4.5),
          axis.text.y  = element_blank(),
          axis.line.x  = element_line(colour = "grey80"),
          axis.line.y  = element_blank(),
          axis.title.x = element_text(colour = "grey40", size = 6),
          axis.ticks.x = element_line(colour = "grey60"),
          axis.ticks.y = element_blank(),
          plot.background = element_rect(colour = "grey70"))
}
