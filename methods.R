
# Repeatable plotting, with standard error bars (weighted standard error calculated using diagis package)
plot_hsas <- function(df, ylab) {
  colnames(df) <- c("HSA", "Year", "DCI.Quantile", "Enrollees", "Measure")
  
  default_colors <- rev(hue_pal()(5))
  num_colors <- nlevels(df$DCI.Quantile)
  if (num_colors == 2) {
    color_values <- default_colors[c(1,5)]
  } else {
    color_values <- default_colors
  }
  
  theme_set(theme_bw())
  ggplot(df, aes(Year, Measure, group=as.factor(HSA), color=as.factor(DCI.Quantile))) + 
    geom_line() +
    scale_color_manual(values = color_values) + 
    facet_grid(. ~ DCI.Quantile) +
    labs(y = ylab) +
    theme(legend.position='none') +
    theme(axis.title.x = element_blank()) + 
    theme(axis.ticks.x = element_blank()) + 
    theme(axis.text.x = element_blank())
}

plot_weighted <- function(df, ylab, legend_title="", legend_values = NULL, linetype_values, linesize_values, y_lims = NULL, y_step = NULL, absolute_percent_dollars, errors = FALSE) {
  colnames(df) <- c("HSA", "Year", "DCI.Quantile", "Enrollees", "Measure")
  if (!is.null(legend_values)) {
    levels(df$DCI.Quantile) <- legend_values
  }
  
  df_by_dci <- df %>%
    group_by(DCI.Quantile, Year) %>%
    summarise(n = dplyr::n(),
              wtd.Measure = weighted.mean(Measure, Enrollees),
              se.Measure = weighted_se(Measure, Enrollees),
              .groups = 'drop')
  
  theme_set(theme_bw())
  p <- ggplot(df_by_dci, aes(Year, wtd.Measure, group=DCI.Quantile)) +
    geom_line(aes(linetype=DCI.Quantile, size=DCI.Quantile)) +
    # scale_color_manual(values = color_values) + 
    scale_linetype_manual(name = legend_title, values = linetype_values) +
    scale_size_manual(name = legend_title, values = linesize_values) +
    labs(y = ylab) +
    theme(legend.position = 'top') +
    theme(legend.title = element_text(size=12)) + 
    theme(legend.text = element_text(size=12)) + 
    theme(axis.text = element_text(size=12)) +
    theme(axis.title.y = element_text(size=14, margin = margin(t = 0, r = 8, b = 0, l = 0))) +
    theme(axis.title.x = element_blank()) +
    theme(axis.ticks = element_line(colour='black', size=0.2)) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    theme(panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.y = element_line(colour='black', size=0.2)) +
    theme(panel.background = element_rect(colour='black', size=0.8)) + 
    theme(plot.margin = margin(0, 20, 5, 10)) +
    scale_x_continuous(limits = c(2003, 2015),
                       breaks = 2003:2015,
                       expand = c(0.01, 0.01),
                       labels = 2003:2015) 
  
  if (errors) {
    p <- p + geom_errorbar(aes(ymin=wtd.Measure - se.Measure, ymax=wtd.Measure + se.Measure), width=.1)
  }
  
  # absolute
  if (absolute_percent_dollars == 0 & !is.null(y_lims) & !is.null(y_step)) {
    p + scale_y_continuous(limits = y_lims,
                           breaks = seq(y_lims[1], y_lims[2], by = y_step),
                           expand = c(0.01, 0.01))
  }
  # percent
  else if (absolute_percent_dollars == 1) {
    if (!is.null(y_lims) & !is.null(y_step)) {
      p + scale_y_continuous(limits = y_lims,
                             breaks = seq(y_lims[1], y_lims[2], by = y_step),
                             expand = c(0.01, 0.01),
                             labels=scales::percent_format(scale = 1, accuracy = 1))
    }
    else {
      p + scale_y_continuous(expand = c(0.01, 0.01),
                             labels=scales::percent_format(scale = 1, accuracy = 1))
    }
  }
  # dollars
  else if (absolute_percent_dollars == 2) {
    if (!is.null(y_lims) & !is.null(y_step)) {
      p + scale_y_continuous(limits = y_lims,
                             breaks = seq(y_lims[1], y_lims[2], by = y_step),
                             expand = c(0.01, 0.01),
                             labels=scales::dollar_format())
    }
    else {
      p + scale_y_continuous(expand = c(0.01, 0.01),
                             labels=scales::dollar_format())
    }
    
  } else {
    p
  }
}

library(weights)


# (Welch t-test)
t_by_year <- function(df, quantile_a, quantile_b, compname_suffix = "") {
  colnames(df) <- c("HSA", "Year", "DCI.Quantile", "Enrollees", "Measure")

  comparison_name <- paste(quantile_a, quantile_b, sep=" - ")
  all_p_vals <- c(paste(comparison_name, compname_suffix, sep=""))
  years <- unique(df$Year)

  for (year in years) {
    sub_df = df[df$Year==year,]
    sub_a = sub_df[sub_df$DCI.Quantile == quantile_a,]
    sub_b = sub_df[sub_df$DCI.Quantile == quantile_b,]
    
    ### Welch t-test (unequal variance)
    t_result <- wtd.t.test(sub_a$Measure, sub_b$Measure, weight=sub_a$Enrollees, weighty=sub_b$Enrollees, samedata=FALSE,
                           alternative="two.tailed", drops="pairwise")$coefficients
    all_p_vals <- c(all_p_vals, t_result["p.value"])
  }

  # uncorrected--to be corrected considering all t-tests done
  all_p_vals <- data.frame(t(all_p_vals))
  colnames(all_p_vals) <- c("Tukey.Pair", years)
  all_p_vals
}

### Repeatable statistical testing (one-way ANOVA at each year)
### NOT TWO-WAY ANOVA WITH 'YEAR' AS THE SECOND CATEGORICAL VARIABLE--WOULD VIOLATE ASSUMPTION OF INDEPENDENCE BETWEEN MEASUREMENTS
anova_by_year <- function(df) {
  colnames(df) <- c("HSA", "Year", "DCI.Quantile", "Enrollees", "Measure")

  all_p_vals <- data.frame()
  years <- unique(df$Year)

  for (year in years) {

    sub_df = df[df$Year==year,]

    ### ANOVA (**not assuming equal variances**) for groups in the given year
    anova_p_val <- oneway.test(Measure ~ DCI.Quantile, data=sub_df, na.action=na.omit, var.equal=FALSE)$p.value
    anova_p_val <- data.frame(anova_p_val)
    rownames(anova_p_val) <- c("ANOVA")
    colnames(anova_p_val) <- c("p_vals")
    
    anova_p_val <- cbind(rownames(anova_p_val), anova_p_val)
    rownames(anova_p_val) <- NULL
    colnames(anova_p_val) <- c("Tukey.Pair", year)

    if (year == years[1]) {
      all_p_vals <- anova_p_val
    } else {
      all_p_vals <- left_join(all_p_vals, anova_p_val, by="Tukey.Pair")
    }
  }

  all_p_vals
}
