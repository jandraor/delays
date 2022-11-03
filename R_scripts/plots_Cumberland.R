Cmb_colour <- "#B84C7A"

plot_Cmb_latent_incidence <- function(sim_df, data_df) {
  
  sim_df <- sim_df |> 
    mutate(label_M_j = paste0("j = ", M_j))
  
  ggplot(sim_df, aes(time, value)) +
    geom_line(aes(group = iter, colour = as.factor(M_j)), alpha = 0.1) +
    facet_grid(parameterisation~label_M_j) +
    scale_colour_manual(values = colours_df$colour) +
    geom_point(data = data_df, aes(y =y), size = 1, colour = Cmb_colour,
               shape = 18) +
    labs(x        = "Days", 
         y        = "Incidence [Cases/day]") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(colour = "grey40"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          strip.background = element_rect(colour = "grey80"))
}

plot_epi_bars <- function(df){
  
  ggplot(df, aes(time, y)) +
    geom_col(fill = Cmb_colour, colour = "white") +
    theme_pubclean() +
    labs(x = "Days since the first case reported",
         y = "Incidence (People per day)") 
}

plot_posterior_pred_checks <- function(incidence_df, data_df) {
  
  incidence_summary <- incidence_df |> group_by(time, M_j, parameterisation) |> 
    summarise(mean = mean(y),
              lb   = quantile(y, 0.025),
              ub   = quantile(y, 0.975),
              .groups = "drop")
  
  ggplot(incidence_summary, aes(time, mean)) +
    geom_ribbon(data = incidence_summary, alpha = 0.1,
                aes(y = mean, ymin = lb, ymax = ub, fill = as.factor(M_j))) +
    geom_line(aes(colour = as.factor(M_j))) +
    facet_grid(parameterisation~M_j) +
    scale_colour_manual(values = colours_df$colour) +
    scale_fill_manual(values = colours_df$colour) +
    geom_point(data = data_df, aes(y = y), size = 1, colour = Cmb_colour,
               shape = 18) +
    labs(x        = "Days", 
         y        = "Incidence [Cases/day]") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(colour = "grey40"),
          axis.line  = element_line(colour = "grey80"),
          axis.text  = element_text(colour = "grey60"),
          axis.ticks = element_line(colour = "grey60"),
          strip.background = element_rect(colour = "grey80"))
}

plot_hist_by_M_j <- function(posteriors_df, var_name, x_label, lims, 
                             n_bins = 30, fixed_tau = NULL) {
  
  caption <- ifelse(is.null(fixed_tau), "",
                    parse(text = "'Generation time' ~ (tau)~': 2.85'"))
  
  posteriors_df <- posteriors_df |> mutate(label_M_j = paste0("j = ", M_j))
  
  subtitle_txt <- paste0("'Histogram:' ~ M^{'1j'}")
  
  
  ggplot(posteriors_df, aes(.data[[var_name]])) +
    scale_x_continuous(limits = lims) +
    geom_histogram(aes(fill = as.factor(M_j)), colour = "white",
                   bins = n_bins) +
    facet_grid(label_M_j ~ parameterisation) +
    scale_fill_manual(values = colours_df$colour) +
    labs(x = parse(text = x_label), y = "",
         caption = caption,
         fill = "j",
         subtitle = parse(text = subtitle_txt)) +
    theme_pubr() +
    theme(legend.position = "none",
          axis.title = element_text(colour = "grey40"),
          axis.line.x  = element_line(colour = "grey80"),
          axis.line.y  = element_blank(),
          axis.text.x  = element_text(colour = "grey60"),
          axis.text.y  = element_blank(),
          axis.ticks.x = element_line(colour = "grey60"),
          axis.ticks.y = element_blank(),
          strip.background = element_rect(colour = "white"),
          plot.subtitle = element_text(colour = "grey40"))
  
}