## File: v1_Fig_S9.R
# Supplementary Fig. 9:  HLA Genotypes in the COV-BOOST cohort by case and noncase status.
# HLA allele coverage heatmap by case status and locus

# ---- Output configuration ----------------------------------------------
# Fits a US-Letter page with 1-inch margins (8.5 - 1 = 7.5 in usable width;
# up to 9 in usable height).
out_dir    <- "Manuscript/figures_v1"   # <-- set destination folder
out_name   <- "v1_Fig_S9"            # filename stem (no extension)
fig_width  <- 7.5                    # inches
fig_height <- 9.0                    # inches (must stay <= 9 for letter)
fig_dpi    <- 300

## Load ADATA with metadata
require(dplyr)
require(ggplot2)
require(readr)
require(stringr)
require(Polychrome)

# ---- Load data ---------------------------------------------------------
project_folder <- "/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey"

day14        <- readr::read_csv(file.path(project_folder, "Manuscript", "adata", "case_status.csv"))
hla          <- readr::read_csv("/Volumes/fh/fast/gilbert_p/kmayerbl/spike_atlas_data/cov_boost_hla_genotypes.csv")
locus_binary <- readRDS("locus_binary.rds")

data <- hla %>%
  filter(sample_id %in% day14$name_1) %>%
  mutate(hla_value = if_else(is.na(hla_value), 0, hla_value)) %>%
  left_join(locus_binary) %>%
  left_join(day14, by = c("sample_id" = "name_1"))

# ---- Clean allele names for display ------------------------------------
# A: '__' (paired alleles) -> '-'
# B: remaining '_'         -> '*'
# C: split 4-digit groups into 'NN:NN'
raw_alleles <- data$binary %>% unique()

clean_names <- raw_alleles %>%
  gsub("__", "-", .) %>%
  gsub("_", "*", .) %>%
  gsub("([0-9]{2})([0-9]{2})(?![0-9])", "\\1:\\2", ., perl = TRUE)

allele_dict <- setNames(clean_names, raw_alleles)
data <- data %>% mutate(binary = allele_dict[binary])

# ---- Palette -----------------------------------------------------------
poly_60 <- as.vector(createPalette(60, c("#D62728", "#3CB44B", "#9467BD")))

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
    axis.ticks.length = unit(2, "pt"),
    
    # Nature Typography Specs (7pt labels, 8pt titles)
    axis.text = element_text(size = 6, colour = "black"),
    axis.title = element_text(size = 6, colour = "black"),
    axis.title.y = element_text(size = 6, margin = margin(r = 2)),
    axis.title.x = element_text(margin = margin(t = 2)),
    
    # Clean facet strip appearance
    strip.background = element_blank(),
    strip.text = element_text(size = 6, face = "plain", colour = "black"),
    
    # Minimize whitespace around each panel
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
    plot.title = element_text(margin = margin(b = 2))
  )

# ---- Plot --------------------------------------------------------------
final_plot_S8 <- data %>%
  filter(hla_value == 1) %>%
  ggplot(aes(x = factor(binary), y = sample_id)) +
  geom_tile(aes(fill = factor(binary)), col = "white") +
  geom_tile(data = data %>% filter(hla_value == 0),
            fill = "#b0b0b0", col = "white") +
  facet_grid(case ~ locus, scales = "free", space = "free") +
  ylab("") + xlab("") +
  scale_fill_manual(
    "",
    values = poly_60,
    guide  = guide_legend(ncol = 10, keywidth = 0.5, keyheight = 0.5)
  ) +
  theme_x +
  theme(
    axis.text.x     = element_text(size = 7, angle = 90, hjust = 1),
    axis.text.y     = element_blank(),#text(size = 5, hjust = 1),
    legend.position = "bottom",
    legend.text     = element_text(size = 5),
    legend.title    = element_text(size = 5),
    legend.key.size = unit(0.2, "cm"),
    legend.margin   = margin(t = -0.5, unit = "cm")
  )

# ---- Save ---------------------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".png")),
  plot = final_plot_S8, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi
)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".pdf")),
  plot = final_plot_S8, width = fig_width, height = fig_height, units = "in"
)
final_plot_S8