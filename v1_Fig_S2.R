# File: v1_Fig_S2.R
# Note: Adapted from 001B_Compare_pre_post includes ts_breadth_figure.pdf  

# ---- Output configuration ----------------------------------------------
# Fits a US-Letter page with 1-inch margins (8.5 - 1 = 7.5 in usable width;
# up to 9 in usable height).
out_dir    <- "Manuscript/figures_v1"   # <-- set destination folder
out_name   <- "v1_Fig_S2"            # filename stem (no extension)
fig_width  <- 7.5                    # inches
fig_height <- 3.0                    # inches (must stay <= 9 for letter)
fig_dpi    <- 300


## Load ADATA with metadata 
require(dplyr)
require(ggplot2)
project_folder = '/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey'
adata   = readr::read_csv(file.path(project_folder , 'Manuscript', 'adata', '001_COV_BOOST_ADATA_WIDE.csv'))
breadth = readr::read_csv(file.path(project_folder , 'Manuscript', 'adata', '001_COV_BOOST_ADATA_WIDE_BREADTH.csv'))
d = readr::read_csv(file.path(project_folder , 'Manuscript', 'adata', '001_COV_BOOST_ADATA.csv'))

r = "/Volumes/fh/fast/gilbert_p/kmayerbl/spike_atlas_data/estimate_20260527/"
scores_fp   = file.path(r, "cov_boost_day0_day14_genotyped_spike_breadth_estimate_x1w2.csv")
d2 = readr::read_csv(scores_fp ) 
d3 = d %>% 
  left_join(d2, by = "name_1")
# --- 2. Reusable plot builder -------------------------------------------
theme_x <- theme_bw(base_size = 10, base_family = "Helvetica") +
  theme(
    # Eliminate all grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Minimal box around data with exact 0.5 linewidth
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    
    # Ensure axis lines and ticks also adhere strictly to 0.5 linewidth
    axis.line = element_blank(), # Covered by panel.border
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(2, "pt"),
    
    # Nature Typography Specs (7pt labels, 8pt titles)
    axis.text = element_text(size = 6, colour = "black"),
    axis.title = element_text(size = 6, colour = "black"),
    axis.title.y = element_text(size =6, margin = margin(r = 2)),
    axis.title.x = element_text(margin = margin(t = 2)),
    
    # Clean facet strip appearance
    strip.background = element_blank(),
    strip.text = element_text(size = 6, face = "bold", colour = "black"),
    
    # Minimize whitespace around each panel
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
    plot.title = element_text(margin = margin(b = 2))
  )


fs = list.files(r)
scores_fp   = file.path(r, "cov_boost_day0_day14_genotyped_spike_breadth_estimate_x1w2.csv")
hits_fp     = file.path(r, "cov_boost_day0_day14_genotyped_spike_breadth_all_hits.csv") 
scores      = readr::read_csv(scores_fp) 
hits        = readr::read_csv(hits_fp)
adata = readr::read_csv(file.path(project_folder , 'Manuscript', 'adata', '001_COV_BOOST_ADATA.csv'))

require(dplyr)
# CHECK COMPLETENESS FOR 240 SAMPLES 
stopifnot(all(adata %>% left_join(scores, by = 'name_1') %>% pull(breadth) %>% is.na() == FALSE))
stopifnot(adata %>% left_join(scores, by = 'name_1') %>% dim() %>% .[1] == 240)
# CASE VS CONTROL 
d3 = adata %>% left_join(scores, by = 'name_1')

# COMPUTE CD4 AND CD8 SEPARATELY
tcell_scores = hits %>%
  filter(hla_detect == 1) %>% 
  select(name_1, total_unique_1, binary, locus = locus_x, x1w2) %>%
  mutate(tcellsub = if_else(locus %in% c("A","B","C"), "CD8", "CD4")) %>% 
  group_by(name_1, total_unique_1, tcellsub) %>% 
  summarise(sum_x1w2 = sum(x1w2)) %>% 
  mutate(breadth = sum_x1w2/total_unique_1 ) %>%
  ungroup()%>%
  select(name_1, tcellsub, breadth) %>%
  mutate(breadth = if_else(is.na(breadth), 1E-5, breadth)) %>%
  tidyr::spread(key = tcellsub , value = breadth) 

d3 %>% left_join(tcell_scores, by = 'name_1') %>%
  ggplot(aes(x = CD8, y = CD4)) + 
  geom_point() + 
  scale_y_log10()+
  scale_x_log10() + 
  annotation_logticks(side = "lb")+
  ggpubr::stat_cor(method = "spearman", cor.coef.name = 'rho') +
  theme_x

d4 = d3 %>% left_join(tcell_scores, by = 'name_1') %>% 
  arrange(sc_repeat_pid, timepoint)

# CONFIRMATORY CHECK -- CD8+CD4 equal total breadth
# ggplot(d4 , aes(x = breadth, y = CD8+CD4)) + 
#   geom_point() + 
#   scale_y_log10()+
#   scale_x_log10() + 
#   annotation_logticks(side = "lb")+
#   ggpubr::stat_cor(method = "spearman", cor.coef.name = 'rho') +
#   theme_x



iqr_labels <- d4 %>%
  group_by(timepoint) %>%
  summarize(
    n = n(),
    median = median(breadth, na.rm = TRUE),
    # Extract the 25th (Q1) and 75th (Q3) percentiles
    q1 = quantile(breadth, 0.25, na.rm = TRUE),
    q3 = quantile(breadth, 0.75, na.rm = TRUE)
  ) %>%
  mutate(
    # Format them together with a dash (e.g., IQR: 5.0e-05 - 2.0e-04)
    label = sprintf("N = %.0f\n%.1e\nIQR:\n[%.1e - %.1e]", n, median, q1, q3),
    y_pos = 1.2e-5 
  )

iqr_labels_cd4 <- d4 %>%
  group_by(timepoint) %>%
  summarize(
    n = n(),
    median = median(CD4, na.rm = TRUE),
    # Extract the 25th (Q1) and 75th (Q3) percentiles
    q1 = quantile(CD4, 0.25, na.rm = TRUE),
    q3 = quantile(CD4, 0.75, na.rm = TRUE)
  ) %>%
  mutate(
    # Format them together with a dash (e.g., IQR: 5.0e-05 - 2.0e-04)
    label = sprintf("N = %.0f\n%.1e\nIQR:\n[%.1e - %.1e]", n, median, q1, q3),
    y_pos = 1.2e-5 
  )
iqr_labels_cd8 <- d4 %>%
  group_by(timepoint) %>%
  summarize(
    n = n(),
    median = median(CD8, na.rm = TRUE),
    # Extract the 25th (Q1) and 75th (Q3) percentiles
    q1 = quantile(CD8, 0.25, na.rm = TRUE),
    q3 = quantile(CD8, 0.75, na.rm = TRUE)
  ) %>%
  mutate(
    # Format them together with a dash (e.g., IQR: 5.0e-05 - 2.0e-04)
    label = sprintf("N = %.0f\n%.1e\nIQR:\n[%.1e - %.1e]", n, median, q1, q3),
    y_pos = 1.2e-5 
  )


library(dplyr)
library(patchwork)
# library(ggplot2) # Assuming already loaded
# library(ggpubr)  # Assuming already loaded

# 1. Create a helper function for the repeated boxplot/quasirandom layers
make_ts_panel <- function(data, y_var, y_label) {
  ggplot(data, aes(x = timepoint, y = .data[[y_var]])) + 
    geom_boxplot(outlier.shape = NA, size = 0.5, width = .5) +
    ggbeeswarm::geom_quasirandom(size = 0.25, width = .25) +
    scale_y_log10(
      breaks = 10^(-10:10), # Forces strictly whole-integer powers of 10
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    annotation_logticks(side = "l",linewidth = .25) +
    coord_cartesian(ylim = c(1E-5, 5E-3)) +
    ggpubr::stat_compare_means(
      label = 'p.signif', paired = TRUE, label.y = -2.75, tip.length = 0,
      comparisons = list(c('Day 0', "Day 14"))
    ) +
    theme_x + 
    labs(x = "", y = y_label) +
    theme(axis.title = element_text(size = 8), 
          axis.text.x = element_text(angle = 0) )
}

# 2. Generate the first three panels using the helper function
panel_both <- make_ts_panel(d4, "breadth", expression("HA-"~T[S]~"breadth"))
panel_cd8  <- make_ts_panel(d4, "CD8", expression(CD8^"+" ~" HA-"~ T[S] ~ "breadth"))
panel_cd4  <- make_ts_panel(d4, "CD4", expression(CD4^"+" ~" HA-"~ T[S] ~ "breadth"))

require(scales)
# 3. Process data and generate the delta scatter plot
panel_delta <- d4 %>% 
  select(sc_repeat_pid, timepoint, CD4, CD8, breadth) %>% 
  tidyr::pivot_wider(names_from = timepoint, values_from = c(CD4, CD8, breadth)) %>% 
  mutate(
    delta_cd4 = `CD4_Day 14` - `CD4_Day 0`,
    delta_cd8 = `CD8_Day 14` - `CD8_Day 0`,
    delta     = `breadth_Day 14` - `breadth_Day 0`
  ) %>% 
  ggplot(aes(x = delta_cd4, y = delta_cd8)) + 
  geom_point(size = 0.25) + 
  ggpubr::stat_cor(cor.coef.name = "rho", size = 3) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  # scale_y_log10(
  #   breaks = 10^(-10:10),
  #   labels = scales::trans_format("log10", scales::math_format(10^.x))
  # ) + 
  # scale_x_log10(
  #   breaks = 10^(-10:10),
  #   labels = scales::trans_format("log10", scales::math_format(10^.x))
  # ) +
  #annotation_logticks(side = "lb") +
  labs(
    x = expression(Delta ~ CD4^"+" ~ " HA-"~  T[S] ~ "breadth"),
    y = expression(Delta ~ CD8^"+" ~ " HA-"~ T[S] ~ "breadth")
  ) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed')+
  geom_vline(aes(xintercept = 0), linetype = 'dashed')+
  theme_x + 
  theme(axis.title = element_text(size = 8)) #+ 
#coord_cartesian(xlim = c(1E-5, 5E-3), ylim = c(1E-5, 5E-3))
panel_delta


# 4. Combine with patchwork
final_plot_S2 <- (panel_both | panel_cd4 | panel_cd8 | panel_delta) +
  plot_layout(widths = c(1, 1, 1, 2)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")) &# Adjust values to tighten
  theme(plot.tag = element_text(size = 10))

final_plot_S2

# ---- Save ---------------------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".png")),
  plot = final_plot_S2 , width = fig_width, height = fig_height, units = "in", dpi = fig_dpi
)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".pdf")),
  plot = final_plot_S2 , width = fig_width, height = fig_height, units = "in"
)


