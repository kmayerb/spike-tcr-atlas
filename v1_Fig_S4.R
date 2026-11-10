## File: v1_Fig_S4.R
# Cross-correlation of TCR depth/breadth metrics vs ELISpot (Day 14 and delta),
# split by CD8 (Class I), CD4 (Class II), and combined.

# ---- Output configuration ----------------------------------------------
# Fits a US-Letter page with 1-inch margins (8.5 - 1 = 7.5 in usable width;
# up to 9 in usable height).
out_dir    <- "Manuscript/figures_v1"   # <-- set destination folder
out_name   <- "v1_Fig_S4"            # filename stem (no extension)
fig_width  <- 7.5                    # inches
fig_height <- 5.0                    # inches (must stay <= 9 for letter)
fig_dpi    <- 300

## Load ADATA with metadata
require(dplyr)
require(ggplot2)
require(patchwork)

# ---- Load data ---------------------------------------------------------
project_folder <- "/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey"
estimate_dir   <- "/Volumes/fh/fast/gilbert_p/kmayerbl/spike_atlas_data/estimate_20260527"

scores_fp <- file.path(estimate_dir, "cov_boost_day0_day14_genotyped_spike_breadth_estimate_x1w2.csv")
hits_fp   <- file.path(estimate_dir, "cov_boost_day0_day14_genotyped_spike_breadth_all_hits.csv")

scores <- readr::read_csv(scores_fp)
hits   <- readr::read_csv(hits_fp)
adata  <- readr::read_csv(
  file.path(project_folder, "Manuscript", "adata", "001_COV_BOOST_ADATA_LONG.csv")
)

locus_binary <- readRDS("locus_binary.rds")
locus_binary <- locus_binary %>%
  mutate(tcellsub = if_else(stringr::str_starts(locus, pattern = "[ABC]"), "CD8", "CD4"))

# ---- Compute Ts metrics (breadth / depth), overall and by T-cell subset -
df_both <- hits %>%
  select(name_1, total_unique_1, binary, hla_detect, templates,
         productive_frequency, x1w2, yule_q_odds) %>%
  left_join(locus_binary) %>%
  filter(hla_detect == 1) %>%
  mutate(depth_yule = productive_frequency * yule_q_odds) %>%
  group_by(name_1, total_unique_1) %>%
  summarise(n = n(), depth = sum(productive_frequency),
            depth2 = sum(depth_yule), s_x1w2 = sum(x1w2), .groups = "drop") %>%
  mutate(breadth = s_x1w2 / total_unique_1)

df_cd4_cd8 <- hits %>%
  select(name_1, total_unique_1, binary, hla_detect, templates,
         productive_frequency, x1w2, yule_q_odds) %>%
  left_join(locus_binary) %>%
  filter(hla_detect == 1) %>%
  mutate(depth_yule = productive_frequency * yule_q_odds) %>%
  group_by(name_1, total_unique_1, tcellsub) %>%
  summarise(n = n(), depth = sum(productive_frequency),
            depth2 = sum(depth_yule), s_x1w2 = sum(x1w2), .groups = "drop") %>%
  mutate(breadth = s_x1w2 / total_unique_1)

df_cd4 <- df_cd4_cd8 %>% filter(tcellsub == "CD4")
df_cd8 <- df_cd4_cd8 %>% filter(tcellsub == "CD8")

df_cd4_cd8_joined <- left_join(df_cd4, df_cd8, suffix = c("_cd4", "_cd8"),
                               by = c("name_1", "total_unique_1"))
df_all <- df_both %>% left_join(df_cd4_cd8_joined, by = c("name_1", "total_unique_1"))

# Attach case / assay metadata and reshape to wide (Day 0 vs Day 14)
df_all <- df_all %>%
  left_join(adata %>% select(name_1, timepoint, sc_repeat_pid, elispot, nt50))
df_all_post <- df_all %>% filter(timepoint == "Day 14")
df_all_pre  <- df_all %>% filter(timepoint == "Day 0")
df_wide <- df_all_post %>%
  left_join(df_all_pre, by = "sc_repeat_pid", suffix = c("_day14", "_day0"))

# ---- Theme -------------------------------------------------------------
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
    axis.ticks.length = unit(3, "pt"),
    
    # Nature Typography Specs (consistent 6pt axis text / titles across figures)
    axis.text = element_text(size = 6, colour = "black"),
    axis.title = element_text(size = 6, colour = "black"),
    axis.title.y = element_text(size = 6, margin = margin(r = 2)),
    axis.title.x = element_text(size = 6, margin = margin(t = 2)),
    
    # Clean facet strip appearance
    strip.background = element_blank(),
    strip.text = element_text(size = 6, face = "bold", colour = "black"),
    
    # Minimize whitespace around each panel
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
    plot.title = element_text(margin = margin(b = 2))
  )

# ---- Reusable scatter builder ------------------------------------------
make_scatter <- function(data, xvar, yvar, xlab_expr = NULL, ylab_expr = NULL) {
  p <- ggplot(data, aes(x = {{ xvar }}, y = {{ yvar }})) +
    geom_point(size = .25) +
    scale_x_log10() +
    scale_y_log10() +
    ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho", size = 2) +
    theme_x
  
  if (!is.null(xlab_expr)) p <- p + xlab(xlab_expr)
  if (!is.null(ylab_expr)) p <- p + ylab(ylab_expr)
  
  p
}

# ---- Combined (CD4 + CD8) ----------------------------------------------
both_depth_d14 <- make_scatter(
  df_wide, depth2_day14, elispot_day14,
  xlab_expr = expression("Depth " * "HA-" * T[S] * " (Day 14)"),
  ylab_expr = expression("SFU/M (ELISpot)")
)
both_breadth_d14 <- make_scatter(
  df_wide, breadth_day14, elispot_day14,
  xlab_expr = expression("Breadth " * "HA-" * T[S] * " (Day 14)"),
  ylab_expr = expression("SFU/M (ELISpot)")
)
both_depth_delta <- make_scatter(
  df_wide, depth2_day14 - depth2_day0, elispot_day14 - elispot_day0,
  xlab_expr = expression(Delta ~ "Depth (" * "HA-" * T[S] - "HA-" * T[S0] * ")"),
  ylab_expr = expression(Delta ~ "SFU/M (ELISpot)")
)
both_breadth_delta <- make_scatter(
  df_wide, breadth_day14 - breadth_day0, elispot_day14 - elispot_day0,
  xlab_expr = expression(Delta ~ "Breadth (" * "HA-" * T[S] - "HA-" * T[S0] * ")"),
  ylab_expr = expression(Delta ~ "SFU/M (ELISpot)")
)

# ---- Class II (CD4) -----------------------------------------------------
cd4_depth_d14 <- make_scatter(
  df_wide, depth2_cd4_day14, elispot_day14,
  xlab_expr = expression("CD4+ Depth " * "HA-" * T[S] * " (Day 14)"),
  ylab_expr = expression("SFU/M (ELISpot)")
)
cd4_breadth_d14 <- make_scatter(
  df_wide, breadth_cd4_day14, elispot_day14,
  xlab_expr = expression("CD4+ Breadth " * "HA-" * T[S] * " (Day 14)"),
  ylab_expr = expression("SFU/M (ELISpot)")
)
cd4_depth_delta <- make_scatter(
  df_wide, depth2_cd4_day14 - depth2_cd4_day0, elispot_day14 - elispot_day0,
  xlab_expr = expression(Delta ~ "CD4+ Depth (" * "HA-" * T[S] - "HA-" *T[S0] * ")"),
  ylab_expr = expression(Delta ~ "SFU/M (ELISpot)")
)
cd4_breadth_delta <- make_scatter(
  df_wide, breadth_cd4_day14 - breadth_cd4_day0, elispot_day14 - elispot_day0,
  xlab_expr = expression(Delta ~ "CD4+ Breadth (" * "HA-" * T[S] - "HA-" *T[S0] * ")"),
  ylab_expr = expression(Delta ~ "SFU/M (ELISpot)")
)

# ---- Class I (CD8) ------------------------------------------------------
cd8_depth_d14 <- make_scatter(
  df_wide, depth2_cd8_day14, elispot_day14,
  xlab_expr = expression("CD8+ Depth " * "HA-" * T[S] * " (Day 14)"),
  ylab_expr = expression("SFU/M (ELISpot)")
)
cd8_breadth_d14 <- make_scatter(
  df_wide, breadth_cd8_day14, elispot_day14,
  xlab_expr = expression("CD8+ Breadth " * "HA-" * T[S] * " (Day 14)"),
  ylab_expr = expression("SFU/M (ELISpot)")
)
cd8_depth_delta <- make_scatter(
  df_wide, depth2_cd8_day14 - depth2_cd8_day0, elispot_day14 - elispot_day0,
  xlab_expr = expression(Delta ~ "CD8+ Depth (" * "HA-" * T[S] - "HA-" *T[S0] * ")"),
  ylab_expr = expression(Delta ~ "SFU/M (ELISpot)")
)
cd8_breadth_delta <- make_scatter(
  df_wide, breadth_cd8_day14 - breadth_cd8_day0, elispot_day14 - elispot_day0,
  xlab_expr = expression(Delta ~ "CD8+ Breadth (" * "HA-" * T[S] - "HA-" *T[S0] * ")"),
  ylab_expr = expression(Delta ~ "SFU/M (ELISpot)")
)

# ---- Assemble grid ------------------------------------------------------
plots <- list(
  cd8_breadth_d14, cd8_depth_d14, cd8_breadth_delta, cd8_depth_delta,
  cd4_breadth_d14, cd4_depth_d14, cd4_breadth_delta, cd4_depth_delta,
  both_breadth_d14, both_depth_d14, both_breadth_delta, both_depth_delta
)

final_plot_S4 <- wrap_plots(plots, nrow = 3, ncol = 4, byrow = TRUE) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8))
final_plot_S4
# ---- Save ---------------------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".png")),
  plot = final_plot_S4, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi
)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".pdf")),
  plot = final_plot_S4, width = fig_width, height = fig_height, units = "in"
)



