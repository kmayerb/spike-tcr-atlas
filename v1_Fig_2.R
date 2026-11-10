# File: v1_Fig_2.R
# Note: refactored from 002_COV_BOOST_BREADTH.R

# ---- Output configuration ----------------------------------------------
# Fits a US-Letter page with 1-inch margins (8.5 - 1 = 7.5 in usable width;
# up to 9 in usable height).
out_dir    <- "Manuscript/figures_v1"   # <-- set destination folder
out_name   <- "v1_Fig_2"        # filename stem (no extension)
fig_width  <- 7.5                    # inches
fig_height <- 2.5                    # inches (must stay <= 9 for letter)
fig_dpi    <- 300
# ---- Load data ---------------------------------------------------------
project_folder <- "/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey"
estimate_dir   <- "/Volumes/fh/fast/gilbert_p/kmayerbl/spike_atlas_data/estimate_20260527"
r = estimate_dir 
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

require(readr)
require(ggplot2)
require(dplyr)
library(patchwork)

# ---- Shared y-axis: power-of-ten labels + matching log breaks ----------
log_y <- scale_y_log10(
  breaks = 10^(-5:-2),
  labels = scales::trans_format("log10", scales::math_format(10^.x))
)

# Reusable left-side log ticks (consistent 0.5-style weighting)
log_ticks_l <- annotation_logticks(
  sides = "l", size = 0.3,
  short = unit(1, "pt"), mid = unit(1.5, "pt"), long = unit(2.5, "pt")
)

fs = list.files('/Volumes/fh/fast/gilbert_p/kmayerbl/spike_atlas_data/estimate_20260527/')
metadata = readr::read_csv('/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_ID159_tools/data/20260520_breadth_across_13_studies_with_imputed_hla_alleles_draft_metadata.csv')
delmonte = read_csv(file.path(r, "delmonte_spike_breadth_estimate_x1w2.csv"))
delmonte_hits = read_csv(file.path(r, "delmonte_spike_breadth_all_hits.csv"))
hits = delmonte_hits

####
delmonte = read_csv(file.path(r, "delmonte_spike_breadth_estimate_x1w2.csv"))
dm = read_tsv('/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_ID159_tools/data/delmonte_metadata.csv')
table(dm$ICU.admit, dm$Severity)
dm = dm %>% select(name_1= sample_name, 
                   days_since_onset = Days_since_onset, 
                   severity = Severity)
delmonte = delmonte %>% left_join(dm, by = "name_1") 
delmonte = delmonte %>% mutate(severity = ifelse(severity == "Moderate_severe", "Moderate", severity))

# Show T[S] breadth captures time since infection
p1 = delmonte  %>%
  filter(days_since_onset < 45) %>%
  mutate(period = cut(days_since_onset, breaks = c(0, 7,14, 21,45))) %>%#pull(period)
  ggplot(aes(x = period , y = breadth))  + 
  geom_boxplot(outlier.shape= NA)+
  ggbeeswarm::geom_quasirandom(size = .1)+
  log_y +
  log_ticks_l +
  geom_smooth(span = 1, se = F) + 
  theme_x +
  xlab("Days")+
  ylab(expression("HA-"*T[S]~"Breadth"))
  ggpubr::stat_compare_means(label = "p.signif", tip.length = 0, size =2, step.increase = 0.05 ,
                             comparisons = list(c('(0,7]','(7,14]') ,c('(7,14]' , '(14,21]'), c('(14,21]' , '(21,45]')))

# Show difference in severity
p2 = delmonte  %>%
  filter(days_since_onset < 45) %>%
  mutate(period = cut(days_since_onset, breaks = c(0, 7,14, 45))) %>%
  ggplot(aes(x = factor(severity, levels = c("Mild", "Moderate","Critical")) , y = breadth))  + 
  geom_boxplot(aes(fill = severity), outlier.shape = NA)+
  ggbeeswarm::geom_quasirandom(size = .1)+
  log_y +
  log_ticks_l +
  geom_smooth(span = 1, se = F) + 
  scale_fill_manual( values = c("Mild"="#FA8072","Critical"= "red", "Moderate"= "#C00000")) + 
  scale_color_manual( values = c("Mild"="#FA8072","Critical"= "red", "Moderate"= "#C00000")) + 
  ggpubr::stat_compare_means(label = "p.signif", tip.length = 0, size = 2,
                             comparisons = list(c("Mild", "Moderate"),      # shortest, drawn lowest
                                                c("Moderate", "Critical"),  # adjacent
                                                c("Mild", "Critical")),     # longest, drawn highest
                             label.y = c(-2.05, -1.85, -1.65)) +            # log10 positions: low -> high
  theme_x + 
  theme(legend.position = "None") +
  xlab("Severity")+
  ylab(expression("HA-"*T[S]~"Breadth"))

# show timeseries
p3 = delmonte  %>%
  filter(days_since_onset < 45) %>%
  ggplot(aes(x = days_since_onset , y = breadth, col = severity))  + 
  geom_point(size = .1)+
  log_y +
  log_ticks_l +
  geom_smooth(span = 1, se = F) + 
  scale_color_manual( values = c("Mild"="#FA8072","Critical"= "red", "Moderate"= "#C00000")) + 
  facet_wrap(~factor(severity, levels = c("Mild", "Moderate","Critical"))) + 
  theme(strip.background = element_blank())+
  xlab("days")+
  ylab(expression("HA-"*T[S]~"Breadth"))+
  theme_x +
  theme(legend.position = "None") 


# Keep the tight ceiling (no added headroom) across all three panels
p1 <- p1 + coord_cartesian(ylim = c(1E-5, 4E-2))
p2 <- p2 + coord_cartesian(ylim = c(1E-5, 4E-2))
p3 <- p3 + coord_cartesian(ylim = c(1E-5, 4E-2))

final_plot <- p1 + p2 + p3 +
  plot_layout(widths = c(1, 1, 3)) +
  plot_annotation(tag_levels = list(c("B", "C", "D"))) &
  theme(plot.tag = element_text(face = "plain", size = 8))
final_plot


# ---- Save ---------------------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".png")),
  plot = final_plot, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi
)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".pdf")),
  plot = final_plot, width = fig_width, height = fig_height, units = "in"
)