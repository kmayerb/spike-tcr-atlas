# =====================================================================
# figure_helpers.R
# ---------------------------------------------------------------------
# Shared setup + plotting functions for the ID159 (Voysey) manuscript
# Figure 1 (main) and Figure S3 (supplement).
# =====================================================================

# ---- 1. Packages ----------------------------------------------------
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(ggpubr)
  library(ggbeeswarm)
  library(survival)
  library(survminer)
  library(rms)
})

# # ---- 2. Paths -------------------------------------------------------
# # Resolve the project root robustly whether the script is run with
# # Rscript, sourced, or run interactively from the project directory.
# .find_project_root <- function() {
#   # (a) Rscript --file=
#   args <- commandArgs(trailingOnly = FALSE)
#   file_arg <- grep("^--file=", args, value = TRUE)
#   if (length(file_arg) == 1) {
#     return(normalizePath(dirname(dirname(sub("^--file=", "", file_arg)))))
#   }
#   # (b) sourced file (has ofile)
#   for (fr in sys.frames()) {
#     if (!is.null(fr$ofile)) {
#       return(normalizePath(dirname(dirname(fr$ofile))))
#     }
#   }
#   # (c) {here} if available
#   if (requireNamespace("here", quietly = TRUE)) {
#     return(here::here())
#   }
#   # (d) fall back to working directory
#   normalizePath(getwd())
# }
# 
# PROJECT_ROOT <- .find_project_root()
# DATA_DIR     <- file.path(PROJECT_ROOT, "data")
# OUTPUT_DIR   <- file.path(PROJECT_ROOT, "output")
# if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
# 
# # Look for an input file in data/ first, then the project root.
# resolve_input <- function(filename) {
#   candidates <- c(file.path(DATA_DIR, filename), file.path(PROJECT_ROOT, filename))
#   hit <- candidates[file.exists(candidates)]
#   if (length(hit) == 0) {
#     stop(sprintf(
#       "Could not find '%s'. Looked in:\n  %s\nPlace the file in one of these locations.",
#       filename, paste(candidates, collapse = "\n  ")
#     ))
#   }
#   hit[1]
# }

# ---- 3. Load + clean data ------------------------------------------
# Replace any -Inf (from log2 of zero) with the column minimum finite value.
replace_neginf <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col) && any(is.infinite(col) & col < 0)) {
      col_min <- min(col[is.finite(col)], na.rm = TRUE)
      col[is.infinite(col) & col < 0] <- col_min
    }
    col
  })
  df
}


# ---- 4. Term labels + colour palette (forest plot) -----------------
# NB: in 001_cox_ph_models_dfall.csv the `term` column is ALREADY recoded
# to these label values, so the palette is keyed on the labels.
term_labels <- c(
  "age10"                    = "Age10Y",
  "maleTRUE"                 = "Male",
  "log2nt50"                 = "Log[2]~NT50",
  "log2elispot"              = "Log[2]~ELISpot",
  "log2tcell"                = "Log[2]~T[s]",
  "log2nt50_pre"             = "Log[2]~NT50pre",
  "log2elispot_pre"          = "Log[2]~ELISpotpre",
  "log2tcell_pre"            = "Log[2]~T[s]~pre",
  "delta_log2nt50"           = "delta~Log[2]~NT50",
  "delta_log2elispot"        = "delta~Log[2]~ELISpot",
  "delta_log2tcell"          = "delta~Log[2]~T[s]",
  "high_log2nt50TRUE"        = "High~NT50",
  "high_log2elispotTRUE"     = "High~ELISpot",
  "high_log2tcellTRUE"       = "High~T[s]",
  "high_log2nt50_preTRUE"    = "High~NT50pre",
  "high_log2elispot_preTRUE" = "High~ELISpotpre",
  "high_log2tcell_preTRUE"   = "High~T[s]~pre"
)

colors <- c(
  "Age10Y"               = "darkgray",
  "Male"                 = "gray",
  "Log[2]~NT50"          = "red",
  "Log[2]~ELISpot"       = "orange",
  "Log[2]~T[s]"          = "navy",
  "Log[2]~NT50pre"       = "salmon",
  "Log[2]~ELISpotpre"    = "gold",
  "Log[2]~T[s]~pre"      = "skyblue",
  "delta~Log[2]~NT50"    = "violet",
  "delta~Log[2]~ELISpot" = "green",
  "delta~Log[2]~T[s]"    = "cyan",
  "High~NT50"            = "darkred",
  "High~ELISpot"         = "darkorange",
  "High~T[s]"            = "navy",
  "High~NT50pre"         = "tomato",
  "High~ELISpotpre"      = "darkgoldenrod",
  "High~T[s]~pre"        = "powderblue"
)

# ---- 5. Marker boxplot (C vs NC, Wilcoxon) -------------------------
# Boxplot + beeswarm of a marker, split C vs NC, coloured by median split
# of the plotted marker, with a Wilcoxon p-value bracket.
plot_marker_box <- function(df,
                            y_var,
                            x_var        = "case_abr",
                            split_var    = NULL,   # default: median split on y_var
                            ylab         = NULL, xlab = "",
                            color_values = c("darkgray", "black"),
                            comps        = list(c("C", "NC")),
                            ylim         = NULL,
                            p_label_y    = NULL,
                            p_label      = "p.format",
                            pt_size      = 0.3) {
  df$value <- df[[y_var]]
  # explicit x grouping column (avoids .data pronoun issues with stat_compare_means)
  df$case_group <- factor(df[[x_var]], levels = c("C", "NC"))
  
  # median split (on split_var if supplied, else the plotted value)
  split_col    <- if (!is.null(split_var)) df[[split_var]] else df$value
  median_value <- median(split_col, na.rm = TRUE)
  df$high      <- split_col > median_value
  
  # bracket position for the p-value label
  bracket_y <- if (!is.null(p_label_y)) {
    p_label_y
  } else if (!is.null(ylim)) {
    step <- diff(ylim) * 0.08
    top  <- ylim[2] - diff(ylim) * 0.02
    seq(from = top, by = -step, length.out = length(comps))
  } else {
    NULL
  }
  
  gg <- ggplot(df, aes(x = case_group, y = value)) +
    geom_boxplot(fill = NA, width = 0.75, size = 0.35, outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(width = 0.3, size = pt_size, aes(col = high)) +
    ggpubr::stat_compare_means(
      tip.length  = 0,
      comparisons = comps,
      method      = "wilcox",
      label       = p_label ,
      size        = 2.5,
      label.y     = bracket_y
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_color_manual(values = color_values) +
    theme_classic(base_size = 10) +
    theme(
      strip.background.x = element_blank(),
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line          = element_blank(),
      plot.margin        = margin(t = 2, r = 2, b = 4, l = 2, unit = "pt"),
      axis.text.x        = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 8),
      axis.title.y       = element_text(angle = 90, size = 8),
      legend.position    = "none"
    ) +
    ylab(ylab) +
    xlab(xlab)
  
  if (!is.null(ylim)) gg <- gg + coord_cartesian(ylim = ylim)
  gg
}

# ---- 6. Correlation scatter (Spearman/Pearson) --------------------
# ---- 6. Correlation scatter (Spearman/Pearson) --------------------
plot_cross <- function(df,
                       x_var, y_var,
                       xlab = NULL, ylab = NULL,
                       color = "black",
                       xlim = NULL, ylim = NULL,
                       method = "spearman",
                       stat_above = FALSE,   # place R/p above all points
                       label_x = NULL,       # optional manual x for the label
                       label_y = NULL) {     # optional manual y for the label
  df$x_val <- df[[x_var]]
  df$y_val <- df[[y_var]]
  if (is.null(xlab)) xlab <- x_var
  if (is.null(ylab)) ylab <- y_var
  
  # position the correlation label above the top-most point when requested
  top_expand <- 0.15
  if (stat_above) {
    fin <- is.finite(df$x_val) & is.finite(df$y_val)
    xr  <- range(df$x_val[fin])
    yr  <- range(df$y_val[fin])
    if (is.null(label_x)) label_x <- xr[1]                    # left edge
    if (is.null(label_y)) label_y <- yr[2] + diff(yr) * 0.02  # just above max
    top_expand <- 0.22                                        # headroom so it fits
  }
  
  gg <- ggplot(df, aes(x = x_val, y = y_val)) +
    geom_point(alpha = 0.8, size = 0.25, color = color) +
    ggpubr::stat_cor(method = method, size = 2.5, label.sep = "\n",
                     label.x = label_x, label.y = label_y) +
    theme_classic(base_size = 10) +
    theme(
      strip.background.x = element_blank(),
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line          = element_blank(),
      plot.margin        = margin(t = 2, r = 2, b = 4, l = 2, unit = "pt"),
      axis.text          = element_text(size = 8),
      axis.title         = element_text(size = 8),
      legend.position    = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, top_expand))) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    ylab(ylab) + xlab(xlab)
  
  if (!is.null(ylim) || !is.null(xlim)) gg <- gg + coord_cartesian(xlim = xlim, ylim = ylim)
  gg
}
# ---- 7. Kaplan-Meier cumulative incidence (median split) ----------
# 1 - S(t), on a median (or forced) split of `y_var`, within one stratum.
plot_cumulative_incidence <- function(df,
                                      y_var,
                                      time_var     = "time",
                                      event_var    = "event",
                                      strata_var   = "strata",
                                      strata_val   = "0",
                                      color_values = c("steelblue", "tomato"),
                                      xlab         = "Days post boost",
                                      ylab         = "Cumulative incidence",
                                      force_median = NULL, # used so we get median from the study as a whole not just strata
                                      legend_labs  = NULL,
                                      xlim         = NULL,
                                      ylim         = c(0, 0.5),
                                      title        = NULL,
                                      pval         = TRUE,
                                      conf_int     = FALSE) {
  
  df <- df[df[[strata_var]] == strata_val, ]
  if (nrow(df) == 0)
    stop(sprintf("No rows left after filtering %s == '%s'.", strata_var, strata_val))
  
  df$value     <- df[[y_var]]
  df$time_var  <- df[[time_var]]
  df$event_var <- df[[event_var]]
  
  median_value <- median(df$value, na.rm = TRUE)
  if (!is.null(force_median)) median_value <- force_median
  df$high <- factor(ifelse(df$value > median_value, "High", "Low"),
                    levels = c("Low", "High"))
  
  km <- survival::survfit(survival::Surv(time_var, event_var) ~ high,
                          data = df, weights = wt, robust = TRUE)
  
  lr_p <- if (pval) {
    p_raw <- survival::survdiff(survival::Surv(time_var, event_var) ~ high, data = df)$chisq |>
      (\(chi) pchisq(chi, df = 1, lower.tail = FALSE))()
    if (p_raw < 0.001) "p < 0.001" else sprintf("p = %.3f", p_raw)
  } else NULL
  
  # lr_p <- if (pval) {
  #   fit   <- survival::coxph(survival::Surv(time_var, event_var) ~ high,
  #                            data = df, weights = wt, robust = TRUE)
  #   p_raw <- summary(fit)$robscore["pvalue"]
  #   if (p_raw < 0.001) "p < 0.001" else sprintf("p = %.3f", p_raw)
  # } else NULL
  
  labs <- if (!is.null(legend_labs)) legend_labs else {
    sprintf(c("Low %s\n(≤ %.2f)", "High %s\n(> %.2f)"),
            y_var, median_value, y_var, median_value)
  }
  
  km_tidy <- survminer::surv_summary(km, data = df)
  km_tidy$group <- factor(ifelse(grepl("High", km_tidy$strata), "High", "Low"),
                          levels = c("Low", "High"))
  km_tidy$cum_inc <- 1 - km_tidy$surv
  km_tidy$ci_low  <- 1 - km_tidy$upper
  km_tidy$ci_high <- 1 - km_tidy$lower
  
  ggplot(km_tidy, aes(x = time, y = cum_inc, colour = group, fill = group)) +
    geom_step(linewidth = 0.6) +
    {if (conf_int) geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, colour = NA)} +
    scale_color_manual(values = color_values, labels = labs) +
    scale_fill_manual(values  = color_values, labels = labs) +
    guides(colour = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) +
    scale_y_continuous(limits = ylim, expand = expansion(mult = c(0.01, 0.05))) +
    {if (!is.null(xlim)) coord_cartesian(xlim = xlim)} +
    {if (!is.null(title)) ggtitle(title)} +
    {if (!is.null(lr_p))
      annotate("text",
               x = ifelse(is.null(xlim), max(km_tidy$time) * 0.05, xlim[1] + diff(xlim) * 0.05),
               y = ylim[2] - diff(ylim) * 0.05,
               label = lr_p, hjust = 0, size = 3)} +
    xlab(xlab) + ylab(ylab) +
    theme_classic(base_size = 10) +
    theme(
      panel.border          = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line             = element_blank(),
      plot.margin           = margin(t = 2, r = 2, b = 4, l = 2, unit = "pt"),
      axis.text             = element_text(size = 8),
      axis.title            = element_text(size = 8),
      legend.title          = element_blank(),
      legend.text           = element_text(size = 8),
      legend.key.size       = unit(10, "pt"),
      legend.position       = c(0.05, 0.5),
      legend.justification  = c("left", "top"),
      legend.box.background = element_rect(colour = "black", linewidth = 0.3, fill = "white"),
      legend.margin         = margin(2, 4, 2, 4, "pt"),
      plot.title            = element_text(size = 8, hjust = 0.5)
    )
}

# ---- 8. Restricted-cubic-spline hazard-ratio curve ----------------
# Fits Surv(time, event) ~ rcs(<var>, knots) + age10 + strat(strata) via
# rms::cph, then plots HR vs the marker with a marginal histogram.
# NOTE: the original code adjusts for age10 (not sex) - see README.
setup_rms_datadist <- function(df) {
  dd <- rms::datadist(df)
  options(datadist = dd)
  invisible(dd)
}

plot_spline_hr <- function(df, var, color, xlab,
                           knots = 3,
                           adjust = "age10 + strat(strata)",
                           time_var = "time", event_var = "event",
                           ylim = c(0, 3), hist_bins = 30) {
  fml <- as.formula(sprintf("Surv(%s, %s) ~ rcs(%s, %d) + %s",
                            time_var, event_var, var, knots, adjust))
  # Register a datadist for THIS model's predictors so rms::Predict has
  # variable limits regardless of call order or any stale global datadist.
  df <- as.data.frame(df)
  rhs_vars <- setdiff(all.vars(fml), c(time_var, event_var))
  dd <- rms::datadist(df[, rhs_vars, drop = FALSE])
  old_opts <- options(datadist = dd)
  on.exit(options(old_opts), add = TRUE)
  
  fit <- rms::cph(fml, data = df, x = TRUE, y = TRUE)
  
  preds <- as.data.frame(rms::Predict(fit, name = var, fun = exp))
  y_max <- max(preds$upper, na.rm = TRUE)
  
  ggplot(preds, aes(x = .data[[var]], y = yhat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray90", alpha = 0.5) +
    geom_line(color = color, linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "darkgray") +
    geom_histogram(
      data = df,
      aes(x = .data[[var]],
          y = after_stat(count) / max(after_stat(count)) * (y_max * 0.2)),
      inherit.aes = FALSE,
      fill = "gray80", color = "white", alpha = 0.6, bins = hist_bins
    ) +
    labs(y = "Hazard Ratio", x = xlab) +
    theme_classic(base_size = 10) +
    theme(
      strip.background = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line        = element_blank(),
      plot.margin      = margin(t = 2, r = 2, b = 4, l = 2, unit = "pt"),
      axis.text        = element_text(size = 8),
      axis.title       = element_text(size = 8),
      legend.position  = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    coord_cartesian(ylim = ylim)
}

# ---- 9. Forest / "rope" panel (from dfall model estimates) --------
# Ported from the enhanced make_forest_panel in 004_Small_Figures.R.
make_forest_panel <- function(data, include_forms, term_order,
                              panel_title = NULL,
                              add_pval = TRUE, add_hr = TRUE,
                              add_hr_right = FALSE,
                              hr_x_pos = 0.1, pval_x_pos = 1.75,
                              hr_right_margin = 110) {   # right margin in pt when add_hr_right = TRUE
  gg <- data %>%
    filter(formula %in% include_forms) %>%
    mutate(formula = factor(formula, levels = include_forms)) %>%
    ggplot(aes(y = fct_relevel(term, term_order), x = hr)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(size = 1.5, aes(col = fct_relevel(term, term_labels))) +
    geom_segment(aes(x = hr_ci95_lo, xend = hr_ci95_hi, yend = term,
                     col = fct_relevel(term, term_labels)),
                 linewidth = 0.75, alpha = 0.75) +
    facet_grid(fct_relevel(formula, include_forms) ~ "",
               scales = "free_y", space = "free_y") +
    scale_color_manual("", values = colors) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2)) +
    scale_y_discrete(labels = \(x) parse(text = x), position = "left") +
    labs(x = "HR", y = "") +
    theme_classic(base_size = 10) +
    theme(panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          axis.line        = element_blank(),
          plot.margin      = margin(t = 2, r = 2, b = 4, l = 2, unit = "pt"),
          axis.text        = element_text(size = 8),
          axis.title       = element_text(size = 8),
          axis.text.y      = element_text(size = 8),
          panel.grid       = element_blank(),
          legend.position  = "none",
          strip.background = element_blank(),
          strip.text       = element_blank())

  if (!is.null(panel_title)) gg <- gg + ggtitle(panel_title)
  
  hr_label <- \(hr, lo, hi) paste0("HR = ", round(hr, 2),
                                   "  (", round(lo, 2), ",", round(hi, 2), ")")
  
  if (add_hr_right) {
    # Text anchored at the right edge of the panel, drawn outside it
    gg <- gg +
      geom_text(aes(label = paste0("   ", hr_label(hr, hr_ci95_lo, hr_ci95_hi)),
                    x = Inf),
                size = 1.5, hjust = 0) +
      coord_cartesian(clip = "off") +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
            plot.margin  = margin(2, hr_right_margin, 4, 2, unit = "pt"))
  } else if (add_hr) {
    gg <- gg + geom_text(
      aes(label = hr_label(hr, hr_ci95_lo, hr_ci95_hi), x = hr_x_pos),
      size = 1.5, hjust = 0)
  }
  
  if (add_pval) {
    gg <- gg + geom_text(
      aes(label = paste0("p = ", round(p.value, 3)), x = pval_x_pos),
      size = 1.5)
  }
  gg
}
