# File: Manuscript/v1_Fig_S7.R
# Supplementary Fig. 7: Difference in Antibody and T Cell immune response across prime and booster vaccine modalities

# Adapted from earlier code 04_ID159_Voysey/001A_COV_BOOST_ADATA_CHECK.R
# This version additionally splits homologous 3x-mRNA recipients as "mRNA(3)".

require(ggplot2)
require(dplyr)

# ---- Output configuration ----------------------------------------------
# Fits a US-Letter page with 1-inch margins (8.5 - 2 = 6.5 in usable width;
# up to 9 in usable height).
out_dir    <- "Manuscript/figures_v1"   # <-- set destination folder
out_name   <- "v1_Fig_S7"            # filename stem (no extension)
fig_width  <- 7.5                    # inches
fig_height <- 6.5                    # inches (must stay <= 9 for letter)
fig_dpi    <- 300

# ---- Helpers ------------------------------------------------------------
# Collapse MOD/BNT into a single "mRNA" category, then split off the
# homologous 3x-mRNA recipients (mRNA prime + mRNA boost) as "mRNA(3)",
# mirroring the existing "ChAdOx(3)" homologous group.
add_boost_product2 <- function(df) {
  df %>%
    mutate(
      boost_product2 = ifelse(boost_product1 %in% c("MOD", "BNT"), "mRNA", boost_product1),
      boost_product2 = if_else(prime == "mRNA" & boost_product1 %in% c("MOD", "BNT"),
                               "mRNA(3)", boost_product2)
    )
}

# ---- Data ---------------------------------------------------------------
adata <- readRDS('001_COV_BOOST_ADATA.rds') %>%
  mutate(
    fold_change_nt50    = nt50 / nt50_pre,
    fold_change_elispot = elispot / elispot_pre,
    boost_product1      = if_else(chadvec_only == 1, "ChAdOx(3)", boost_product)
  ) %>%
  add_boost_product2()

# ---- Shared plotting settings ------------------------------------------
vac_levels2 <- c("ChAdOx(3)", "ChAdOx", "Ad26", "NVX", "mRNA(3)","mRNA" )

boost_comparisons <- list(
  c("ChAdOx(3)", "ChAdOx"),
  c("ChAdOx(3)", "Ad26"),
  c("ChAdOx(3)", "NVX"),
  c("ChAdOx(3)", "mRNA(3)"),
  c("ChAdOx(3)", "mRNA"),
  c("mRNA", "mRNA(3)"))
  
prime_comparison <- list(c("ChAdOx", "mRNA"))

prime_colors <- c("#A1C9F4", "#FFB482")

theme_panel <- theme_bw() +
  theme(
    axis.text.y     = element_text(size = 6, hjust = .5),
    axis.text.x     = element_text(angle = 0, size = 6, hjust = .5),
    axis.title      = element_text(size = 8, hjust = .5),
    panel.grid      = element_blank(),
    legend.position = "none"
  )

# Wrapper for the p.signif brackets used on every panel.
signif_brackets <- function(comparisons, label.y = NULL, step.increase = 0.05) {
  ggpubr::stat_compare_means(
    size          = 2,
    comparisons   = comparisons,
    step.increase = step.increase,
    tip.length    = 0,
    label         = "p.signif",
    label.y       = label.y
  )
}

# =========================================================================
# Row 1: NT50 (neutralising antibody)
# =========================================================================
gg_nt50_day0 <- ggplot(
  adata %>% filter(timepoint == "Day 14"),
  aes(x = factor(prime), y = nt50_pre)
) +
  geom_boxplot(aes(fill = prime), outlier.shape = NA, size = .5) +
  ggbeeswarm::geom_quasirandom(size = .5, width = .25) +
  scale_y_log10() +
  scale_fill_manual("", values = prime_colors) +
  signif_brackets(prime_comparison, label.y = 3, step.increase = 0.04) +
  coord_cartesian(ylim = c(10, 10000)) +
  ylab("NT50 (Day 0)") + xlab("Primary vaccine") +
  theme_panel

gg_nt50_day14 <- ggplot(
  adata %>% filter(timepoint == "Day 14"),
  aes(x = factor(boost_product2, levels = vac_levels2), y = nt50)
) +
  geom_boxplot(outlier.shape = NA, size = .5) +
  ggbeeswarm::geom_quasirandom(size = .75, aes(col = prime), width = .25) +
  scale_y_log10() +
  scale_color_manual("", values = prime_colors) +
  signif_brackets(boost_comparisons, label.y = 3.0) +
  coord_cartesian(ylim = c(10, 10000)) +
  ylab("NT50 (Day 28)") + xlab("Booster product") +
  theme_panel

gg_nt50_fc <- ggplot(
  adata %>% filter(timepoint == "Day 14"),
  aes(x = factor(boost_product2, levels = vac_levels2), y = fold_change_nt50)
) +
  geom_boxplot(outlier.shape = NA, size = .5) +
  ggbeeswarm::geom_quasirandom(size = .75, aes(col = prime), width = .25) +
  scale_y_log10() +
  scale_color_manual("", values = prime_colors) +
  signif_brackets(boost_comparisons, label.y = 1.5) +
  ylab("NT50 (Fold change)") + xlab("Booster product") +
  theme_panel

# =========================================================================
# Row 2: ELISpot
# =========================================================================
gg_elispot_day0 <- ggplot(
  adata %>% filter(timepoint == "Day 14"),
  aes(x = factor(prime), y = elispot_pre)
) +
  geom_boxplot(aes(fill = prime), outlier.shape = NA, size = .5) +
  ggbeeswarm::geom_quasirandom(size = .5, width = .25) +
  scale_y_log10() +
  scale_fill_manual("", values = prime_colors) +
  signif_brackets(prime_comparison, label.y = 3, step.increase = 0.04) +
  coord_cartesian(ylim = c(1, 10000)) +
  ylab("ELISpot (Day 0)") + xlab("Primary vaccine") +
  theme_panel

gg_elispot_day14 <- ggplot(
  adata %>% filter(timepoint == "Day 14"),
  aes(x = factor(boost_product2, levels = vac_levels2), y = elispot)
) +
  geom_boxplot(outlier.shape = NA, size = .5) +
  ggbeeswarm::geom_quasirandom(size = .75, aes(col = prime), width = .25) +
  scale_y_log10() +
  scale_color_manual("", values = prime_colors) +
  signif_brackets(boost_comparisons, label.y = 3) +
  coord_cartesian(ylim = c(1, 10000)) +
  ylab("ELISpot (Day 14)") + xlab("Booster product") +
  theme_panel

gg_elispot_fc <- ggplot(
  adata %>% filter(timepoint == "Day 14"),
  aes(x = factor(boost_product2, levels = vac_levels2), y = fold_change_elispot)
) +
  geom_boxplot(outlier.shape = NA, size = .5) +
  ggbeeswarm::geom_quasirandom(size = .75, aes(col = prime), width = .25) +
  scale_y_log10() +
  scale_color_manual("", values = prime_colors) +
  signif_brackets(boost_comparisons, label.y = .75) +
  ylab("ELISpot (Fold change)") + xlab("Booster product") +
  theme_panel

# =========================================================================
# Row 3: T-cell breadth (Ts)
# =========================================================================
adata14_by_prime <- readRDS('adata14_by_prime.rds') %>%
  add_boost_product2()

gg_breadth_day0 <- ggplot(
  adata14_by_prime,
  aes(x = factor(prime), y = breadth_pre)
) +
  geom_boxplot(aes(fill = factor(prime)), outlier.shape = NA, size = .5) +
  ggbeeswarm::geom_quasirandom(size = .5, width = .25) +
  scale_y_log10() +
  geom_abline(linetype = "dashed") +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho", label.y.npc = 1) +
  scale_fill_manual(values = prime_colors) +
  signif_brackets(prime_comparison) +
  coord_cartesian(ylim = c(5e-5, 5e-3)) +
  ylab(expression("HA-"~T[s] ~ "(Day 0)")) + xlab("Primary vaccine") +
  theme_panel

gg_breadth_day14 <- ggplot(
  adata14_by_prime,
  aes(x = factor(boost_product2, levels = vac_levels2), y = breadth)
) +
  geom_boxplot(outlier.shape = NA, size = .5) +
  ggbeeswarm::geom_quasirandom(size = .75, width = .25, aes(col = prime)) +
  scale_y_log10() +
  geom_abline(linetype = "dashed") +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho", label.y.npc = 1) +
  scale_color_manual(values = prime_colors) +
  signif_brackets(boost_comparisons, label.y = -3) +
  coord_cartesian(ylim = c(5e-5, 5e-3)) +
  ylab(expression("HA-"~T[s] ~ "(Day 14)")) + xlab("Booster product") +
  theme_panel

gg_breadth_fc <- adata14_by_prime %>%
  mutate(fc_breadth = breadth / breadth_pre) %>%
  mutate(fc_breadth > 15, 15, fc_breadth) %>%
  ggplot(aes(x = factor(boost_product2, levels = vac_levels2), y = fc_breadth)) +
  geom_boxplot(outlier.shape = NA, size = .5) +
  ggbeeswarm::geom_quasirandom(size = .75, width = .25, aes(col = prime)) +
  scale_y_log10() +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho", label.y.npc = 1) +
  scale_color_manual(values = prime_colors) +
  signif_brackets(boost_comparisons, label.y = 1) +
  coord_cartesian(ylim = c(1, 100)) +
  ylab(expression("HA-"~T[s] ~ "(Fold change)")) + xlab("Booster product") +
  theme_panel

# =========================================================================
# Assemble Figure S6 (3 assays x 3 timepoints) with panel labels A-C
# =========================================================================
fig_s6 <- ggpubr::ggarrange(
  ggpubr::ggarrange(gg_breadth_day0, gg_breadth_day14, gg_breadth_fc, nrow = 1,
                    widths = c(0.5, 1, 1),
                    labels = c("A", "", ""), font.label = list(size = 10, face = "plain" )),
  ggpubr::ggarrange(gg_nt50_day0,    gg_nt50_day14,    gg_nt50_fc,    nrow = 1,
                    widths = c(0.5, 1, 1),
                    labels = c("B", "", ""), font.label = list(size = 10 ,  face = "plain" )),
  ggpubr::ggarrange(gg_elispot_day0, gg_elispot_day14, gg_elispot_fc, nrow = 1,
                    widths = c(0.5, 1, 1),
                    labels = c("C", "", ""), font.label = list(size = 10,  face = "plain"  )),
  ncol = 1
)

# ---- Save ---------------------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".png")),
  plot = fig_s6, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi
)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".pdf")),
  plot = fig_s6, width = fig_width, height = fig_height, units = "in"
)

fig_s6   # also draw to the active device






# =========================================================================
# Panel source data (journal "Source Data" export)
# =========================================================================
# Reconstructs the data plotted in each of the nine panels using only the
# variables mapped in that panel: the x-axis category, the plotted value,
# and the `prime` grouping used for point colour. The log10 axes only render
# finite positive values, so those are what is retained.

src_dir <- file.path(out_dir, "v1_Fig_S6_alt_source_data")
if (!dir.exists(src_dir)) dir.create(src_dir, recursive = TRUE)

# Day-14 subset shared by the NT50 and ELISpot panels
adata14 <- adata %>% filter(timepoint == "Day 14")

# Fold-change breadth: same transform + < 100 filter used in gg_breadth_fc
breadth_fc_df <- adata14_by_prime %>%
  mutate(fc_breadth = breadth / breadth_pre) %>%
  mutate(fc_breadth >15, 15, fc_breadth) 

# Pull x-category / prime / value for one panel into a tidy table.
extract_panel <- function(df, panel, description, x_var, y_var) {
  df %>%
    transmute(
      panel       = panel,
      description = description,
      x_axis      = as.character(.data[[x_var]]),
      prime       = as.character(prime),
      value       = .data[[y_var]]
    ) #%>%
    #filter(is.finite(value) & value > 0)   # log10 axis shows finite positives only
}

panel_data_list <- list(
  # Row A: T-cell breadth (HA-Ts)
  extract_panel(adata14_by_prime, "A1", "HA-Ts (Day 0)",         "prime",          "breadth_pre"),
  extract_panel(adata14_by_prime, "A2", "HA-Ts (Day 14)",        "boost_product2", "breadth"),
  extract_panel(breadth_fc_df,    "A3", "HA-Ts (Fold change)",   "boost_product2", "fc_breadth"),
  # Row B: NT50
  extract_panel(adata14,          "B1", "NT50 (Day 0)",          "prime",          "nt50_pre"),
  extract_panel(adata14,          "B2", "NT50 (Day 28)",         "boost_product2", "nt50"),
  extract_panel(adata14,          "B3", "NT50 (Fold change)",    "boost_product2", "fold_change_nt50"),
  # Row C: ELISpot
  extract_panel(adata14,          "C1", "ELISpot (Day 0)",       "prime",          "elispot_pre"),
  extract_panel(adata14,          "C2", "ELISpot (Day 14)",      "boost_product2", "elispot"),
  extract_panel(adata14,          "C3", "ELISpot (Fold change)", "boost_product2", "fold_change_elispot")
)
names(panel_data_list) <- vapply(panel_data_list, function(d) d$panel[1], character(1))

# One CSV per panel ...
for (nm in names(panel_data_list)) {
  write.csv(panel_data_list[[nm]],
            file.path(src_dir, paste0(out_name, "_", nm, ".csv")),
            row.names = FALSE)
}

# ... plus a single combined long table
panel_data_long <- dplyr::bind_rows(panel_data_list)
write.csv(panel_data_long,
          file.path(src_dir, paste0(out_name, "_source_data.csv")),
          row.names = FALSE)

panel_data_wide = left_join(
  adata %>% filter(timepoint == "Day 14") %>% 
    select(name_1, prime, boost_product2, elispot, elispot_pre, nt50, nt50_pre), 
  adata14_by_prime %>% filter(timepoint == "Day 14") %>% 
    select(name_1, prime, boost_product2, breadth, breadth_pre), 
  by = c("name_1", "prime","boost_product2"))

write.csv(panel_data_wide,
          file.path(src_dir, paste0(out_name, "_source_data_wide.csv")),
          row.names = FALSE)

