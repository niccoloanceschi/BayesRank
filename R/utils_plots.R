
#' Plot posterior distributions of incremental treatment effects
#'
#' Plots histograms of posterior draws for the "Low" and "High vs. Low"
#' incremental treatment coefficients, faceted by Generation and Rechallenge.
#'
#' @param draws_mcmc Postprocessed MCMC output from `postprocess_BayesRank`
#' @param show Display the plot to the active graphics device? Default TRUE.
#' @param save Save the plot to a PDF file? Default FALSE.
#' @param filename Output file name (without extension). Default 'coeff_incremental'.
#' @param dir Output directory. Default '~/Desktop'.
#'
#' @return Invisibly returns the ggplot object.
#'
#' @export
#' 
plot_coeff_incremental <- function(draws_mcmc, show=T, save=F, filename='coeff_incremental', dir='~/Desktop'){
  
  coef_mcmc <- draws_mcmc$alpha
  
  plot_data <- as_tibble(coef_mcmc) %>%
    pivot_longer(everything()) %>%
    separate(name, into = c("Effect", "Generation", "Retreatment"), sep = ", ") %>%
    mutate(
      Generation = factor(Generation, levels = c("P0", "F1")),
      Retreatment = factor(Retreatment, levels = c("Un-Rechallenged", "Rechallenged")),
      Effect = factor(Effect, levels = c("Low", "High vs. Low")))
  
  my_labeller <- labeller(
    Generation  = c("P0" = "Generation : P0", "F1" = "Generation : F1"),
    Retreatment = c("Un-Rechallenged" = "Un-Rechallenged", "Rechallenged" = "Rechallenged"),
    Effect      = c("Low" = "P0 Early Dose : Low", "High vs. Low" = "P0 Early Dose : High"))
  
  coef_hist <- ggplot(plot_data, aes(x = value, fill = Generation)) +
    # geom_histogram(bins = 40, fill = "steelblue", color = "white", linewidth = 0.1) +
    geom_histogram(bins = 40, color = "white", linewidth = 0.1) +
    scale_fill_manual(values = c("P0" = "#56B4E9", "F1" = "#E69F00")) +
    facet_nested(Effect ~ Generation + Retreatment, labeller = my_labeller,
                 scales = "free_y", independent = "y") +    # makes y truly independent per panel
    theme_minimal() +
    labs(title = "Posterior Distributions of Incremental Treatment Effects",
         x = "Latent Damage Score Increase", y = "Count") +
    theme(
      ggh4x.facet.nestline = element_line(color = "black", linewidth = 0.25),
      strip.background = element_rect(fill = "gray95", color = NA), 
      strip.text = element_text(size = 9, face = "bold"),
      panel.spacing = unit(1, "lines"), legend.position = "none")
  
  if(show) print(coef_hist)
  if(save) ggsave(filename=file.path(dir,paste0(filename,'.pdf')), coef_hist, width=12, height=5)
  invisible(coef_hist)
}


#' Plot posterior distributions of cumulative treatment effects
#'
#' Plots histograms of posterior draws for the "Low" and "High"
#' cumulative treatment coefficients, faceted by Generation and Rechallenge.
#'
#' @param draws_mcmc Postprocessed MCMC output from `postprocess_BayesRank`
#' @param show Display the plot to the active graphics device? Default TRUE.
#' @param save Save the plot to a PDF file? Default FALSE.
#' @param filename Output file name (without extension). Default 'coeff_cumulative'.
#' @param dir Output directory. Default '~/Desktop'.
#'
#' @return Invisibly returns the ggplot object.
#'
#' @export
#' 
plot_coeff_cumulative <- function(draws_mcmc, show=T, save=F, filename='coeff_cumulative', dir='~/Desktop'){
  
  coef_mcmc <- draws_mcmc$beta
  
  plot_data_cumulative <- as_tibble(coef_mcmc) %>%
    pivot_longer(everything()) %>%
    separate(name, into = c("Effect", "Generation", "Retreatment"), sep = ", ") %>%
    mutate(
      Generation  = factor(Generation,  levels = c("P0", "F1")),
      Retreatment = factor(Retreatment, levels = c("Un-Rechallenged", "Rechallenged")),
      Effect      = factor(Effect,      levels = c("Low", "High")))
  
  my_labeller <- labeller(
    Generation  = c("P0" = "Generation : P0", "F1" = "Generation : F1"),
    Retreatment = c("Un-Rechallenged" = "Un-Rechallenged", "Rechallenged" = "Rechallenged"),
    Effect      = c("Low" = "P0 Early Dose : Low", "High" = "P0 Early Dose : High"))
  
  coef_hist_cumulative <- ggplot(plot_data_cumulative, aes(x = value, fill = Generation)) +
    geom_histogram(bins = 40, color = "white", linewidth = 0.1) +
    scale_fill_manual(values = c("P0" = "#56B4E9", "F1" = "#E69F00")) +
    facet_nested(Effect ~ Generation + Retreatment, labeller = my_labeller,
                 scales = "free_y", independent = "y") +
    theme_minimal() +
    labs(title = "Posterior Distributions of Cumulative Treatment Effects",
         x = "Latent Damage Score Increase", y = "Count") +
    theme(ggh4x.facet.nestline = element_line(color = "black", linewidth = 0.25),
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text = element_text(size = 9, face = "bold"),
      panel.spacing = unit(1, "lines"), legend.position = "none")

  if(show) print(coef_hist_cumulative)
  if(save) ggsave(filename=file.path(dir,paste0(filename,'.pdf')), coef_hist_cumulative, width=12, height=5)
  invisible(coef_hist_cumulative)
}

#' Plot posterior distributions of batch effects and intraclass correlation
#'
#' Plots posterior distributions of replicate-level intercept shifts
#' alongside the intraclass correlation (ICC), combined into a single
#' side-by-side figure.
#'
#' @param draws_mcmc Postprocessed MCMC output from `postprocess_BayesRank`
#' @param show Display the plot to the active graphics device? Default TRUE.
#' @param save Save the plot to a PDF file? Default FALSE.
#' @param filename Output file name (without extension). Default 'replicates_icc'.
#' @param dir Output directory. Default '~/Desktop'.
#'
#' @return Invisibly returns the ggplot object.
#'
#' @export
#' 
plot_replicates_ICC <- function(draws_mcmc, show=T, save=F, filename='replicates_icc', dir='~/Desktop'){
  
  rep_int  <- draws_mcmc$mu
  icc_samp <- draws_mcmc$icc
  
  df_icc <- tibble(Value = icc_samp, Parameter = "ICC")
  df_rep <- as_tibble(rep_int) %>%
    pivot_longer(everything(), names_to = "Replicate", values_to = "Value") %>%
    mutate(Replicate = factor(Replicate, levels = c("Replicate 3", "Replicate 2", "Replicate 1")))
  
  plot_replicates <- ggplot(df_rep, aes(x = Value, y = Replicate)) +
    stat_halfeye(.width = c(0.95), fill = "#9b59b6",  alpha = 0.4,             
      color = "black", slab_color = "#7d3c98", slab_linewidth = 0.8) +
    scale_y_discrete(expand = expansion(mult = c(0.05, 0))) + 
    theme_minimal() +
    labs(title = "Posterior Distributions of Batch Effects Shifts", x = NULL, y = NULL) +
    theme(axis.text.y = element_text(face = "bold", size = 11),
          axis.text.x = element_text(size = 11), panel.grid.minor = element_blank())
  
  plot_icc <- ggplot(df_icc, aes(x = Value, y = Parameter)) +
    stat_halfeye(.width = c(0.95), fill = "#66c2a5", alpha = 0.4,
      color = "black", slab_color = "#1b9e77", slab_linewidth = 0.8) +
    scale_y_discrete(expand = expansion(mult = c(0.05, 0))) + 
    theme_minimal() +
    labs(title = "Posterior Distribution of ICC", x = NULL, y = NULL) +
    theme(axis.text.y = element_text(face = "bold", size = 11),
          axis.text.x = element_text(size = 11), panel.grid.minor = element_blank())
  
  combined_plot <- plot_replicates + plot_icc + plot_layout(widths = c(1, 1))
  
  if(show) print(combined_plot)
  if(save) ggsave(filename=file.path(dir,paste0(filename,'.pdf')), combined_plot, width=12, height=5)
  invisible(combined_plot)
}

#' Plot category-wise posterior distributions of latent scores
#'
#' Plots overlapping histograms of the latent utilities, groupe and colored by
#' observed ordinal category, to visualize separation between adjacent groups.
#'
#' @param draws_mcmc Postprocessed MCMC output from `postprocess_BayesRank`
#' @param show Display the plot to the active graphics device? Default TRUE.
#' @param save Save the plot to a PDF file? Default FALSE.
#' @param filename Output file name (without extension). Default 'hist_Z_scores'.
#' @param dir Output directory. Default '~/Desktop'.
#'
#' @return Invisibly returns the ggplot object.
#'
#' @export
#' 
plot_latent_scores <- function(draws_mcmc, show=T, save=F, filename='hist_Z_scores', dir='~/Desktop'){

  Zy = draws_mcmc$Zy
  
  custom_colors <- c(rev(viridis(length(Zy), option = 'viridis')))
  plot_colors <- custom_colors[1:length(Zy)]
  
  df_hist <- lapply(1:length(Zy), function(i) {
    data.frame(value = as.vector(Zy[[i]]), category = factor(i))
  }) %>% bind_rows()
  
  combined_hist <- ggplot(df_hist, aes(x = value, fill = category, color = category)) +
    geom_histogram(aes(y = after_stat(count / sum(count))),
      position="identity", alpha=0.6, binwidth=0.1, linewidth=0.1, boundary=0) +
    scale_fill_manual(values=plot_colors, labels=0:(length(Zy)-1), name="y",
      guide=guide_legend(reverse=TRUE)) +
    scale_color_manual(values=plot_colors, labels=0:(length(Zy)-1), name="y",
      guide=guide_legend(reverse=TRUE)) +
    labs(title="Latent Scores: Category-wise Posterior Distributions",
      x=expression(Z[ij]), y="Proportion") +
    theme_bw() + theme(legend.position="right", panel.grid.minor=element_blank())
  
  if(show) print(combined_hist)
  
  out_file = file.path(dir,paste0(filename,'.pdf'))
  if(save) ggsave(filename=out_file, combined_hist, width=6, height=4, device = cairo_pdf)
  invisible(combined_hist)
}

#' Plot posterior medians of random effects against category boundaries
#'
#' Plots an histogram of posterior medians of the worm-level random effects.
#' Overlays vertical reference lines and colored markers, showing the  
#' estimated boundaries on the latent scale among ordinal categories.
#'
#' @param draws_mcmc Postprocessed MCMC output from `postprocess_BayesRank`
#' @param show Display the plot to the active graphics device? Default TRUE.
#' @param save Save the plot to a PDF file? Default FALSE.
#' @param filename Output file name (without extension). Default 'hist_RE_medians'.
#' @param dir Output directory. Default '~/Desktop'.
#'
#' @return Invisibly returns the ggplot object.
#'
#' @export
#' 
plot_random_effects <- function(draws_mcmc, show=T, save=F, filename='hist_RE_medians', dir='~/Desktop'){
  
  Zmax = draws_mcmc$Zmax
  Zmin = draws_mcmc$Zmin
  
  custom_colors <- c(rev(viridis(ncol(Zmax), option = 'viridis')))
  plot_colors <- custom_colors[1:ncol(Zmax)]
  
  df_median_RE <- tibble(Value = draws_mcmc$median_RE, Parameter = "RE")
  
  y_sep_mean <- 0.5*(colMeans(Zmax)[-ncol(Zmax)]+colMeans(Zmin)[-1])
  inner_centers <- (y_sep_mean[-length(y_sep_mean)] + y_sep_mean[-1]) / 2
  zone_centers <- c(
    2*y_sep_mean[1] - inner_centers[1],
    inner_centers,
    2*y_sep_mean[(ncol(Zmax)-1)] - inner_centers[length(inner_centers)]
  )
  df_bullets <- tibble(x = zone_centers, y = 1.8, color = plot_colors)
  
  plot_RE_extra <- ggplot(df_median_RE, aes(x = Value, y = Parameter)) +
    geom_segment(data = tibble(x = y_sep_mean),
                 aes(x = x, xend = x, y = 1.1, yend = 1.95),
                 linewidth = 0.5, linetype = 1, inherit.aes = FALSE,
                 color = alpha('black', alpha = 0.5)) +
    stat_halfeye(slab_color="#DC267F", slab_linewidth=0.8, .width=c(0.95), 
                 fill="#DC267F", alpha=0.3, interval_size=-100, point_size=-100) +
    geom_point(data = df_bullets,  aes(x = x, y = y), color = df_bullets$color,
               size = 6, inherit.aes = FALSE) +
    scale_x_continuous(breaks = scales::breaks_width(1)) +
    scale_y_discrete(expand = expansion(mult = c(0.05, 0.15))) +
    theme_minimal() + labs(title = "Posterior Medians of RE", x = NULL, y = NULL) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 11),
          axis.text.y = element_text(face = "bold", size = 11))
  
  if(show) print(plot_RE_extra)
  if(save) ggsave(filename=file.path(dir,paste0(filename,'.pdf')), plot_RE_extra, width=6, height=4)
  invisible(plot_RE_extra)
}

#' Plot a confusion matrix of observed vs. predicted neuron scores
#'
#' Plots a heatmap-style confusion matrix comparing observed and predicted ordinal categories
#'
#' @param Observed Vector of observed ordinal categories.
#' @param Predicted Vector of predicted ordinal categories.
#' @param method Optional label identifying the prediction method, appended to the plot title. Default NULL.
#' @param colors Viridis color palette option used for cell fill. Default 'mako'.
#' @param show Display the plot to the active graphics device? Default TRUE.
#' @param save Save the plot to a PDF file? Default FALSE.
#' @param filename Output file name (without extension). Default 'confusion_matrix'.
#' @param dir Output directory. Default '~/Desktop'.
#'
#' @return Invisibly returns the ggplot object.
#'
#' @export
#' 
plot_predictions <- function(Observed,Predicted,method=NULL,colors='mako',show=T,
                             save=F, filename='confusion_matrix', dir='~/Desktop'){
  
  full_title <- "Neuroscore Predictions"
  if(!is.null(method)) full_title <- paste0(full_title," (",method,")")
  
  conf_mat <- table(Observed = Observed, Predicted = Predicted)
  conf_df <- as.data.frame(conf_mat)
  conf_df$Freq[conf_df$Freq == 0] <- NA
  
  plot_conf_counts <- ggplot(conf_df, aes(y = Predicted, x = Observed, fill = Freq)) +
    geom_tile() + geom_text(aes(label = ifelse(is.na(Freq), "0", Freq)), color = "white") +
    scale_fill_viridis_c(trans = "log10", na.value = "white",option = colors,direction = -1) +
    theme_minimal() + theme(plot.title = element_text(size = 12)) + labs(title = full_title)
  
  if(show) print(plot_conf_counts)
  if(save) ggsave(filename=file.path(dir,paste0(filename,'.pdf')), plot_conf_counts, width=4.2, height=4)
  invisible(plot_conf_counts)
}
