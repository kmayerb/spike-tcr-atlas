# =====================================================================
# v1_Fig_S5.R 
#v1_Fig_Sens.R — Sensitivity analysis figure (ID159 / Voysey)
#
# Question: do the day-14 TCR breadth results hold up when some samples
# were sequenced too shallowly (< 100,000 productive templates)?
#
# Three strategies, each fit with and without repertoire-richness and
# NT50 adjustment (9 Cox models total, panel G):
#   1 — all samples, observed breadth               (n = 120)
#   2 — restrict to well-sequenced samples only     (n =  98)
#   3 — impute day-14 breadth for shallow samples
#       from a day-0 -> day-14 model                (n = 120)
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(patchwork); library(scales)
  library(survival); library(ggbeeswarm); library(ggpubr)
})

# ---- 0. Paths -------------------------------------------------------
r        <- '/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey'
adata    <- file.path(r, 'Manuscript', 'adata')
raw_file <- '/Volumes/fh/fast/gilbert_p/fg_data/voysey_covid/2025-03-21_TCR/Oxford_TCR_SampleOverview_03-12-2025_11-56-32_PM.tsv'

df  <- read_csv(file.path(adata, '001_COV_BOOST_ADATA_WIDE.csv'),         show_col_types = FALSE)
d   <- read_csv(file.path(adata, '001_COV_BOOST_ADATA_WIDE_BREADTH.csv'), show_col_types = FALSE)
raw <- read_tsv(raw_file, show_col_types = FALSE)

# ---- 1. Constants ---------------------------------------------------
DEPTH_CUTOFF  <- 1e5     # productive templates required for "good" recovery
BREADTH_FLOOR <- 1e-5    # zero-breadth samples clamped here before log2 (= -16.6)
TRAIN_MIN_D0  <- -15     # day-0 breadth below this is at/near the floor; excluded from the imputation fit

# ---- 2. Analysis frame ---------------------------------------------
# ---- 2. Analysis frame ---------------------------------------------
# NOTE: left_join(d) deliberately joins on ALL shared columns (natural join).
# Do not pin this to a single key — the multi-column match is intentional.
data <- df %>%
  select(event, time, strata, wt, age10, male,
         log2nt50_pre, log2nt50, log2elispot,
         sc_repeat_pid, name_1_day14 = name_1, case, boost_product) %>%
  left_join(d) %>%
  select(event, time, strata, wt, age10, male,
         log2nt50_pre, log2nt50, log2elispot,
         sc_repeat_pid, name_1_day14, name_1_day0,
         breadth_day14, breadth_day0, case, boost_product,
         total_unique_1_day14, total_unique_1_day0) %>%
  # Clamp zero/near-zero breadth so log2 is finite
  mutate(breadth_day14 = ifelse(breadth_day14 < BREADTH_FLOOR, BREADTH_FLOOR, breadth_day14)) %>%
  mutate(breadth_day0  = ifelse(breadth_day0  < BREADTH_FLOOR, BREADTH_FLOOR, breadth_day0)) %>%
  mutate(log2_breadth_day0  = log2(breadth_day0)) %>%
  mutate(log2_breadth_day14 = log2(breadth_day14)) %>%
  select(event, time, strata, wt, age10, male,
         log2nt50_pre, log2nt50, log2elispot,
         sc_repeat_pid, name_1_day14, name_1_day0,
         log2_breadth_day14, log2_breadth_day0,
         case, boost_product,
         total_unique_1_day14, total_unique_1_day0) %>%
  # Sequencing depth per timepoint, from the raw sample overview
  left_join(select(raw, name_1_day14 = sample_name,
                   productive_templates_day14 = productive_templates)) %>%
  left_join(select(raw, name_1_day0 = sample_name,
                   productive_templates_day0 = productive_templates))

low_recovery <- data %>% filter(productive_templates_day14 < DEPTH_CUTOFF) %>% pull(name_1_day14)
data <- data %>% mutate(low_recovery_day14 = ifelse(name_1_day14 %in% low_recovery, 1, 0))

# Repertoire richness, inverted so larger = *less* diverse.
# NOTE: HR < 1 on this term therefore means higher hazard with MORE diversity.
data$log2_richness <- log2(1 / data$total_unique_1_day14)

# Well-sequenced subset: strategy-2 models and the imputation model fit
data_train <- data %>% filter(low_recovery_day14 == 0)
data_train$log2_richness <- log2(1 / data_train$total_unique_1_day14)
# Well-sequenced subset: used for the strategy-2 models and to fit the imputation model
data_train <- data %>% filter(low_recovery_day14 == 0)
data_train$log2_richness <- log2(1 / data_train$total_unique_1_day14)

cat("Low-recovery day-14 samples:", sum(data$low_recovery_day14), "of", nrow(data), "\n")

# ---- 3. Shared theme + scales --------------------------------------
base_pt   <- 8
theme_8pt <- theme_classic(base_size = base_pt) +
  theme(
    plot.title      = element_text(size = base_pt),
    plot.subtitle   = element_text(size = base_pt - 1),
    axis.title      = element_text(size = base_pt),
    axis.text       = element_text(size = base_pt - 1),
    legend.title    = element_text(size = base_pt),
    legend.text     = element_text(size = base_pt - 1),
    legend.key.size = unit(0.35, "cm"),
    strip.text      = element_text(size = base_pt - 1)
  )

shape_scale <- scale_shape_manual(
  "Day-14 sequencing", values = c(`0` = 20, `1` = 15),
  labels = c(`0` = "OK (\u2265100k templates)", `1` = "Low (<100k)"))

recovery_color <- scale_color_manual(
  "Day-14 sequencing", values = c(`0` = "black", `1` = "#D55E00"),
  labels = c(`0` = "OK (\u2265100k)", `1` = "Low (<100k)"))

# ---- 4. Panels A-C: sequencing depth distributions ------------------
# Per-sample depth, joined back to case status, at each timepoint.
raw2 <- data %>%
  select(sample_name = name_1_day14, case, total_unique_1_day14) %>%
  inner_join(raw, by = "sample_name") %>%
  mutate(low_recovery = productive_templates < DEPTH_CUTOFF)

# Density + histogram of a depth measure, split by case status, log-x.
# `n_breaks` sets the histogram resolution (bin edges are log-spaced).
make_plot <- function(var, label, df = raw2, n_breaks = 31) {
  log_breaks <- 10^seq(4, 6, length.out = n_breaks)   # 1e4 .. 1e6
  ggplot(df, aes(x = .data[[var]], col = case, fill = case)) +
    geom_histogram(aes(y = after_stat(density)), breaks = log_breaks,
                   alpha = 0.3, position = "identity", color = NA) +
    geom_density(fill = NA) +
    scale_fill_manual("",  values = c("red", "black")) +
    scale_color_manual("", values = c("red", "black")) +
    scale_x_log10() +
    coord_cartesian(xlim = c(1e4, 1e6)) +
    annotation_logticks(sides = "b", outside = FALSE) +
    xlab(label) +
    theme_8pt
}

k0 <- make_plot("productive_templates",
                "Total productive templates\n(Raw data)", n_breaks = 51) +
  # Dashed line marks the depth cutoff that defines "low recovery"
  geom_vline(xintercept = DEPTH_CUTOFF, col = "gray", linetype = "dashed")

k2 <- make_plot("total_unique_1_day14",
                "Total unique TCRs\n(Resampled to 100,000 templates)\n(All samples, n = 120)",
                n_breaks = 31)

k3 <- make_plot("total_unique_1_day14",
                "Total unique TCRs\n(Resampled to 100,000 templates)\n(Samples > 100,000 templates, n = 98)",
                df = raw2 %>% filter(!low_recovery), n_breaks = 21)

# ---- 5. Panel F: breadth and richness by case status ----------------
# Two measures on one axis: breadth (upper cloud) and inverse richness
# (lower cloud). Wilcoxon tests are computed here, not by stat_compare_means,
# so the bracket labels are guaranteed to match the reported p-values.
make_dual_violin <- function(data_set, xlab = "Breakthrough status") {
  d <- data_set %>% mutate(inv_total = log2(1 / total_unique_1_day14))
  
  p_breadth <- wilcox.test(log2_breadth_day14 ~ case, data = d)$p.value
  p_total   <- wilcox.test(inv_total          ~ case, data = d)$p.value
  
  to_stars <- function(p) as.character(
    cut(p, c(-Inf, 1e-4, 1e-3, 0.01, 0.05, Inf),
        labels = c("****", "***", "**", "*", "ns")))
  
  # Bracket sits just above the top of each cloud
  fmax <- function(x) max(x[is.finite(x)], na.rm = TRUE)
  brackets <- tibble(y     = c(fmax(d$log2_breadth_day14), fmax(d$inv_total)) + 0.4,
                     label = c(to_stars(p_breadth), to_stars(p_total)))
  tip <- 0.2
  
  # Right-hand axis restates the inverse-richness scale as 1-in-N fractions
  frac_labels <- expression(frac(1, "66,000"), frac(1, "16,000"),
                            frac(1, "4000"),   frac(1, "1,000"), frac(1, "250"))
  
  ggplot(d, aes(x = case)) +
    geom_violin(aes(y = log2_breadth_day14, fill = case), col = NA, alpha = 0.25) +
    geom_violin(aes(y = inv_total,          fill = case), col = NA, alpha = 0.05) +
    geom_quasirandom(aes(y = log2_breadth_day14), size = 0.1) +
    geom_quasirandom(aes(y = inv_total),          size = 0.1) +
    # significance bracket: bar, two tips, label
    geom_segment(data = brackets, aes(x = 1, xend = 2, y = y, yend = y),       inherit.aes = FALSE) +
    geom_segment(data = brackets, aes(x = 1, xend = 1, y = y, yend = y - tip), inherit.aes = FALSE) +
    geom_segment(data = brackets, aes(x = 2, xend = 2, y = y, yend = y - tip), inherit.aes = FALSE) +
    geom_text(data = brackets, aes(x = 1.5, y = y + 0.3, label = label),
              inherit.aes = FALSE, size = 3) +
    coord_cartesian(ylim = c(-16.5, -8), clip = "off") +
    scale_y_continuous(
      name   = expression(log[2](frac(1, Total_Unique)) ~ "    " ~
                            log[2](frac("SARS-CoV-2 S hits", Total_Unique))),
      breaks = c(-16, -14, -12, -10, -8),
      sec.axis = sec_axis(~ ., breaks = log2(1 / c(66000, 16000, 4000, 1000, 250)),
                          labels = frac_labels, name = NULL)) +
    scale_fill_manual("", values = c("red", "black")) +
    xlab(xlab) +
    theme_8pt
}

k4 <- make_dual_violin(data_train,
                       "Breakthrough status\n(Samples > 100,000 templates, n = 98)")

# ---- 6. Imputation model (feeds panels D-E and strategy 3) ----------
# Day-0 breadth predicts day-14 breadth. Fit on well-sequenced samples
# only, and drop day-0 values at the clamp floor (they are not real
# measurements and would anchor the slope).
train_df <- data_train %>% filter(log2_breadth_day0 > TRAIN_MIN_D0)
mod <- lm(log2_breadth_day14 ~ log2_breadth_day0, data = train_df)
summary(mod)

# How far outside the training range are we predicting?
d0_range <- range(train_df$log2_breadth_day0)
n_extrap <- data %>% filter(low_recovery_day14 == 1,
                            log2_breadth_day0 < d0_range[1] |
                              log2_breadth_day0 > d0_range[2]) %>% nrow()
cat("Imputed samples extrapolating beyond the day-0 training range:", n_extrap, "\n")

data <- data %>%
  mutate(log2_breadth_day14_predicted = as.numeric(predict(mod, newdata = .)),
         # Recovered variable: imputed where sequencing was shallow, observed otherwise
         log2_breadth_day14_recover = if_else(low_recovery_day14 == 1,
                                              log2_breadth_day14_predicted,
                                              log2_breadth_day14))

# Fitted line for the panel; `mod` has no product term, so one line only
pred_grid <- tibble(log2_breadth_day0 = seq(d0_range[1], d0_range[2], length.out = 100))
pred_grid$log2_breadth_day14 <- predict(mod, newdata = pred_grid)

# ---- 7. Panel D: the imputation model --------------------------------
# Arrows run from each shallow sample's observed value to its imputed
# value, stopping just short so the head does not sit on the fitted line.
gap <- 0.1
arrow_df <- data %>%
  filter(low_recovery_day14 == 1,
         abs(log2_breadth_day14_predicted - log2_breadth_day14) > gap) %>%
  mutate(yend_short = log2_breadth_day14_predicted -
           sign(log2_breadth_day14_predicted - log2_breadth_day14) * gap)

p0 <- ggplot(train_df, aes(x = log2_breadth_day0, y = log2_breadth_day14)) +
  geom_point(aes(shape = factor(low_recovery_day14)), size = 1) +
  shape_scale +
  geom_line(data = pred_grid, linewidth = 0.4) +
  geom_segment(data = arrow_df,
               aes(x = log2_breadth_day0, xend = log2_breadth_day0,
                   y = log2_breadth_day14, yend = yend_short),
               col = "#D55E00", linewidth = 0.3,
               arrow = arrow(length = unit(0.05, "in"), type = "closed")) +
  ggpubr::stat_cor(label.x.npc = .25, method = "pearson")+
  labs(subtitle = "Fit imputation model Day-0 \u2192 Day-14",
       x = expression(log[2] ~ "TCR breadth, pre-boost (Day 0)"),
       y = expression(log[2] ~ "TCR breadth, post-boost (Day 14)")) +
  # coord_cartesian, not xlim(): xlim() would drop rows and change stat_cor
  coord_cartesian(xlim = c(-15, -10)) +
  theme_8pt + theme(legend.position = "bottom")
p0
# ---- 8. Panel E: recovered breadth by case status --------------------
p3 <- ggplot(data, aes(y = log2_breadth_day14_recover, x = case)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  geom_quasirandom(aes(shape = factor(low_recovery_day14),
                       color = factor(low_recovery_day14)), size = 1) +
  shape_scale + recovery_color +
  stat_compare_means(size = base_pt / .pt) +
  labs(subtitle = "Imputed Day-14 breadth for samples with < 100,000 templates",
       x = "Breakthrough status",
       y = expression(log[2] ~ "Day-14 breadth")) +
  theme_8pt
p3
# ---- 9. Panel G: Cox models ------------------------------------------
source(file.path("Manuscript", "v1_Fig_1_helper_functions.R"))  # term_labels, colors, make_forest_panel()

# Weighted fits; robust = TRUE makes the sandwich SE explicit (coxph would
# use it anyway with non-integer weights). All CIs/p-values below are robust.
fit_primary_model                       <- coxph(Surv(time, event) ~ log2_breadth_day14 + age10 + male + strata(strata),                            data = data,       weights = wt, robust = TRUE)
fit_primary_model_rich                  <- coxph(Surv(time, event) ~ log2_breadth_day14 + log2_richness + age10 + male + strata(strata),            data = data,       weights = wt, robust = TRUE)
fit_primary_model_rich_nt50             <- coxph(Surv(time, event) ~ log2_breadth_day14 + log2_richness + log2nt50 + age10 + male + strata(strata), data = data,       weights = wt, robust = TRUE)

fit_sensitivity_high_recovery           <- coxph(Surv(time, event) ~ log2_breadth_day14 + age10 + male + strata(strata),                            data = data_train, weights = wt, robust = TRUE)
fit_sensitivity_high_recovery_rich      <- coxph(Surv(time, event) ~ log2_breadth_day14 + log2_richness + age10 + male + strata(strata),            data = data_train, weights = wt, robust = TRUE)
fit_sensitivity_high_recovery_rich_nt50 <- coxph(Surv(time, event) ~ log2_breadth_day14 + log2_richness + log2nt50 + age10 + male + strata(strata), data = data_train, weights = wt, robust = TRUE)

fit_sensitivity_imputed_low             <- coxph(Surv(time, event) ~ log2_breadth_day14_recover + age10 + male + strata(strata),                            data = data, weights = wt, robust = TRUE)
fit_sensitivity_imputed_low_rich        <- coxph(Surv(time, event) ~ log2_breadth_day14_recover + log2_richness + age10 + male + strata(strata),            data = data, weights = wt, robust = TRUE)
fit_sensitivity_imputed_low_rich_nt50   <- coxph(Surv(time, event) ~ log2_breadth_day14_recover + log2_richness + log2nt50 + age10 + male + strata(strata), data = data, weights = wt, robust = TRUE)

# Pull HR / robust CI / p into the column names make_forest_panel() expects
extract_cox <- function(fit, model_name) {
  s <- summary(fit)
  data.frame(
    formula    = model_name,                 # facet variable
    term       = rownames(s$coefficients),
    hr         = s$conf.int[, "exp(coef)"],
    hr_ci95_lo = s$conf.int[, "lower .95"],
    hr_ci95_hi = s$conf.int[, "upper .95"],
    p.value    = s$coefficients[, "Pr(>|z|)"],
    n          = fit$n,
    n_events   = fit$nevent,
    row.names  = NULL)
}

model_labels <- c("1  (n=120)", "1A (n=120)", "1B (n=120)",
                  "2  (n=98)",  "2A (n=98)",  "2B (n=98)",
                  "3  (n=120)", "3A (n=120)", "3B (n=120)")

results <- rbind(
  extract_cox(fit_primary_model,                       model_labels[1]),
  extract_cox(fit_primary_model_rich,                  model_labels[2]),
  extract_cox(fit_primary_model_rich_nt50,             model_labels[3]),
  extract_cox(fit_sensitivity_high_recovery,           model_labels[4]),
  extract_cox(fit_sensitivity_high_recovery_rich,      model_labels[5]),
  extract_cox(fit_sensitivity_high_recovery_rich_nt50, model_labels[6]),
  extract_cox(fit_sensitivity_imputed_low,             model_labels[7]),
  extract_cox(fit_sensitivity_imputed_low_rich,        model_labels[8]),
  extract_cox(fit_sensitivity_imputed_low_rich_nt50,   model_labels[9]))

# Confirm the n in each label matches the fit (NT50 models lose rows if log2nt50 has NAs)
print(distinct(results, formula, n, n_events))

# Observed and imputed breadth share one row/colour so the models line up
term_labels <- c(term_labels,
                 "log2_breadth_day14"         = "Log[2]~T[s]",
                 "log2_breadth_day14_recover" = "Log[2]~T[s]",
                 "log2_richness"              = "Log[2]~(1/R)",
                 "male"                       = "Male")
colors <- c(colors, "Log[2]~(1/R)" = "purple")   # log2nt50 already mapped in the helper

stopifnot(length(setdiff(unique(results$term), names(term_labels))) == 0)
results$term <- recode(results$term, !!!term_labels)

rope_term_order <- c("Log[2]~NT50", "Log[2]~(1/R)", "Log[2]~T[s]", "Age10Y", "Male")

p_rope <- make_forest_panel(results,
                            include_forms = model_labels,
                            term_order    = rope_term_order,
                            hr_x_pos      = 0.06,
                            pval_x_pos    = 5.5) +
  scale_x_log10(breaks = c(0.5, 1, 2)) +
  # coord_cartesian, not scale limits: limits would turn out-of-range CI
  # endpoints into NA and delete whole segments
  coord_cartesian(xlim = c(0.05, 7), clip = "off") +
  theme(strip.text   = element_text(size = 7),
        strip.text.y = element_text(size = 8, angle = -90, hjust = 0.5),
        axis.text.y = element_text(size = 7))

# ---- 10. Assemble ----------------------------------------------------
# Tags run column-wise: A-C depth, D-F imputation, G models
final <- ((k0 / k2 / k3) | p_rope | (k4 / p0 / p3 ) ) +
  plot_layout(guides = "collect", widths = c(1, 1.3, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right",
        plot.tag = element_text(size = 9, face = "plain"))

final
r = "/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey"
ggsave(file.path(r, "Manuscript", "figures_v1", "v1_Fig_S5.pdf"), final, width = 11, height = 9)
ggsave(file.path(r, "Manuscript", "figures_v1", "v1_Fig_S5.png"), final, width = 11, height = 9)



