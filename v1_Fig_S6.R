# File: Manuscript/v1_Fig_S6.R
# Supplementary Fig. 6: HLA-aware T cell breadth distinguishes breakthrough cases from non-cases more than HLA-agnostic TCRseq-based metrics



# Note: 9/2 this code adapted from formalized/001A_COV_BOOST_ADATA_CHECK.R

require(ggplot2)
require(dplyr)

# ---- Output configuration ----------------------------------------------
# Fits a US-Letter page with 1-inch margins (8.5 - 1 = 7.5 in usable width;
# up to 9 in usable height).
out_dir    <- "Manuscript/figures_v1"   # <-- set destination folder
out_name   <- "v1_Fig_S6"            # filename stem (no extension)
fig_width  <- 7.5                    # inches
fig_height <- 6.0                    # inches (must stay <= 9 for letter)
fig_dpi    <- 300


## Load ADATA with metadata 
require(dplyr)
require(ggplot2)
project_folder = '/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey'
adata   = readr::read_csv(file.path(project_folder , 'Manuscript', 'adata', '001_COV_BOOST_ADATA_WIDE.csv'))
breadth = readr::read_csv(file.path(project_folder , 'Manuscript', 'adata', '001_COV_BOOST_ADATA_WIDE_BREADTH.csv'))

# <adata> contains adaptive generated metrics
# We confer the partipant ID on breadth
name_1_pid_day14 = adata %>% select(name_1_day14 = name_1, sc_repeat_pid)
breadth = breadth %>% left_join(name_1_pid_day14, by = c("name_1_day14"))
stopifnot(dim(breadth)[1] == 120)

### NOTE THAT WE NEED TO DEAL WITH NON_DETECTS IN ADAPTIVES SET
replace_zeros_with_min <- function(df, col_name) {
  # Extract the specific column
  col_data <- df[[col_name]]
  # Calculate the minimum value, excluding zeros and ignoring NAs
  min_val <- min(col_data[col_data != 0], na.rm = TRUE)
    # Replace the zeros with the calculated minimum
  df[[col_name]][col_data == 0] <- min_val
  return(df)
}

replace_zeros_with_fixed_min <- function(df, col_name, fixed_min = 1E-5) {
  # Extract the specific column
  col_data <- df[[col_name]]
  # Replace the zeros with the calculated minimum
  df[[col_name]][col_data <= 0] <-  fixed_min
  return(df)
}


# --- 1. One join, all columns, done once --------------------------------
adata_joined <- adata %>%
  left_join(
    select(breadth, name_1_day14,
           breadth_day14, breadth_day0,
           breadth_cd4_day14, breadth_cd4_day0,
           breadth_cd8_day14, breadth_cd8_day0),
    by = c("name_1" = "name_1_day14")
  ) %>%   mutate(case2 = if_else(case == "case", "C", "NC"))


y_vars <- c(
  "breadth_day0", "breadth_day14",
  "breadth_cd4_day14", "breadth_cd4_day0",
  "breadth_cd8_day14", "breadth_cd8_day0",
  "adpt_breadth_spike_pre", "adpt_breadth_spike",
  "adpt_COVID_CD4_spike_breadth_pre",
  "adpt_COVID_CD8_spike_breadth_pre",
  "adpt_COVID_CD4_spike_breadth",
  "adpt_COVID_CD8_spike_breadth",
  # --- newly added ---
  "adpt_breadth_MIRA",
  "adpt_score_v3",
  "adpt_score_spike",
  "adpt_breadth_MIRA_pre",
  "adpt_score_v3_pre",
  "adpt_score_spike_pre"
)
# HERE WE ENSURE WE KEEP ZEROS AFTER LOG TRANSFROMATION
adata_joined <- Reduce(
  function(df, col) replace_zeros_with_fixed_min(df, col),
  y_vars,
  init = adata_joined
)

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

log10_nd_labels <- function(x, nd_value = 1e-5) {
  labs <- scales::trans_format("log10", scales::math_format(10^.x))(x)
  labs <- as.list(labs)
  is_nd <- vapply(x, function(v) isTRUE(all.equal(v, nd_value)), logical(1))
  labs[is_nd] <- "ND"
  do.call(expression, labs)
}
# --- Reusable plot builder -------------------------------------------
plot_case_comparison <- function(data, y_var, ylim = c(1e-5, 5e-3), ylab = NULL,
                                 title = NULL,
                                 label_nd = TRUE, nd_value = 1e-5,
                                 violin_fill = "lightgray" ,
                                 comparisons = list(c("C", "NC"))) {
  
  gg <- ggplot(data, aes(x = case2, y = .data[[y_var]])) +
    geom_violin(fill = violin_fill, col = NA) +
    geom_boxplot(width = .1, fill = "white", col = "black", outlier.shape = NA) +
    ggpubr::stat_compare_means(
      tip.length = 0, label = "p.label", size = 2,
      comparisons = comparisons
    ) +
    scale_y_log10(
      labels = if (label_nd) {
        function(x) log10_nd_labels(x, nd_value = nd_value)
      } else {
        scales::trans_format("log10", scales::math_format(10^.x))
      }
    ) +
    theme_x +
    coord_cartesian(ylim = ylim) +
    annotation_logticks(side = "l") +
    xlab("")
  
  if (!is.null(ylab)) {
    gg <- gg + ylab(ylab)
  }
  
  if (!is.null(title)) {
    gg <- gg + ggtitle(title) +
      theme(plot.title = element_text(size = 8, hjust = 0.5))
  }
  
  gg
}

# --- Fill color for each y-variable --------------------------------------
fill_map <- list(
  breadth_day0                       = "lightgray",
  breadth_day14                      = "lightgray",
  breadth_cd4_day0                   = "lightgray",
  breadth_cd4_day14                  = "lightgray",
  breadth_cd8_day0                   = "lightgray",
  breadth_cd8_day14                  = "lightgray",
  adpt_breadth_spike_pre             = "#cba69e",
  adpt_breadth_spike                 = "#cba69e",
  adpt_COVID_CD4_spike_breadth_pre   = "#cba69e",
  adpt_COVID_CD8_spike_breadth_pre   = "#cba69e",
  adpt_COVID_CD4_spike_breadth       = "#cba69e",
  adpt_COVID_CD8_spike_breadth       = "#cba69e",
  adpt_breadth_MIRA           = "#cba69e",
  adpt_score_v3                = "#cba69e",
  adpt_score_spike            = "#cba69e",
  adpt_breadth_MIRA_pre       = "#cba69e",
  adpt_score_v3_pre            = "#cba69e",
  adpt_score_spike_pre        = "#cba69e"
)

# Clean, human-readable labels for each y-variable
ylab_map <- list(
  breadth_day0                       = expression("HA-"*T[S]~"Breadth"),
  breadth_day14                      = expression("HA-"*T[S]~"Breadth"),
  breadth_cd4_day0                   = expression("HA-"*T[S]~"CD4+ Breadth"),
  breadth_cd4_day14                  = expression("HA-"*T[S]~"CD4+ Breadth"),
  breadth_cd8_day0                   = expression("HA-"*T[S]~"CD8+ Breadth"),
  breadth_cd8_day14                  = expression("HA-"*T[S]~"CD8+ Breadth"),
  adpt_breadth_spike_pre             = "Adaptive Spike Breadth",
  adpt_breadth_spike                 = "Adaptive Spike Breadth",
  adpt_COVID_CD4_spike_breadth_pre   = "Adaptive CD4+ Spike Breadth",
  adpt_COVID_CD8_spike_breadth_pre   = "Adaptive CD8+ Spike Breadth",
  adpt_COVID_CD4_spike_breadth       = "Adaptive CD4+ Spike Breadth",
  adpt_COVID_CD8_spike_breadth       = "Adaptive CD8+ Spike Breadth"
)

# Day 0 / Day 14 title for each panel
title_map <- list(
  breadth_day0                       = "Day 0",
  breadth_day14                      = "Day 14",
  breadth_cd4_day0                   = "Day 0",
  breadth_cd4_day14                  = "Day 14",
  breadth_cd8_day0                   = "Day 0",
  breadth_cd8_day14                  = "Day 14",
  adpt_breadth_spike_pre             = "Day 0",
  adpt_breadth_spike                 = "Day 14",
  adpt_COVID_CD4_spike_breadth_pre   = "Day 0",
  adpt_COVID_CD8_spike_breadth_pre   = "Day 0",
  adpt_COVID_CD4_spike_breadth       = "Day 14",
  adpt_COVID_CD8_spike_breadth       = "Day 14"
)

pt <- purrr::map(y_vars, ~ plot_case_comparison(
  adata_joined, .x,
  ylab  = ylab_map[[.x]],
  title = title_map[[.x]],
  violin_fill = fill_map[[.x]] %||% "lightgray"
)) %>%
  setNames(y_vars)

row1 <- pt[["breadth_day0"]] | pt[["breadth_day14"]] |
  pt[["breadth_cd4_day0"]] | pt[["breadth_cd4_day14"]] |
  pt[["breadth_cd8_day0"]] | pt[["breadth_cd8_day14"]]
row2 <- pt[["adpt_breadth_spike_pre"]] | pt[["adpt_breadth_spike"]] |
  pt[["adpt_COVID_CD4_spike_breadth_pre"]] | pt[["adpt_COVID_CD4_spike_breadth"]] |
  pt[["adpt_COVID_CD8_spike_breadth_pre"]] | pt[["adpt_COVID_CD8_spike_breadth"]]







# --- ylab / title entries for the newly added variables -----------------
ylab_map <- c(ylab_map, list(
  adpt_breadth_MIRA           = "Adaptive MIRA Breadth, Day 14",
  adpt_score_v3                = "Adaptive Score v3, Day 14",
  adpt_score_spike            = "Adaptive Spike Score, Day 14",
  adpt_breadth_MIRA_pre       = "Adaptive MIRA Breadth, Day 0",
  adpt_score_v3_pre            = "Adaptive Score v3, Day 0",
  adpt_score_spike_pre        = "Adaptive Spike Score, Day 0"
))

title_map <- c(title_map, list(
  adpt_breadth_MIRA           = "Day 14",
  adpt_score_v3                = "Day 14",
  adpt_score_spike            = "Day 14",
  adpt_breadth_MIRA_pre       = "Day 0",
  adpt_score_v3_pre            = "Day 0",
  adpt_score_spike_pre        = "Day 0"
))


# --- Regenerate pt with the full y_vars (now includes new variables) ----
pt <- purrr::map(y_vars, ~ plot_case_comparison(
  adata_joined, .x,
  ylab  = ylab_map[[.x]],
  title = title_map[[.x]],
  violin_fill = fill_map[[.x]] %||% "lightgray"
)) %>%
  setNames(y_vars)
# --- New comparison figure: newly added variables only -------------------
remove_logticks <- function(gg) {
  gg$layers <- gg$layers[!vapply(gg$layers, function(l) inherits(l$geom, "GeomLogticks"), logical(1))]
  gg
}
new_vars_row1 <-
  pt[["adpt_breadth_MIRA_pre"]] |
  pt[["adpt_breadth_MIRA"]] |
  (remove_logticks(pt[["adpt_score_v3_pre"]])   + scale_y_continuous() + coord_cartesian(ylim = c(1e-5, 50))) |
  (remove_logticks(pt[["adpt_score_v3"]])       + scale_y_continuous() + coord_cartesian(ylim = c(1e-5, 50))) |
  (remove_logticks(pt[["adpt_score_spike_pre"]]) + scale_y_continuous() + coord_cartesian(ylim = c(1e-5, 50))) |
  (remove_logticks(pt[["adpt_score_spike"]])     + scale_y_continuous() + coord_cartesian(ylim = c(1e-5, 50)))


# 
# cell_subsets_comparison <- (row1 / row2 / new_vars_row1) +
#   patchwork::plot_annotation(tag_levels = "A")
# 
# 



# --- Analogous pairs: primary (y) vs adaptive (x) ------------------------
pair_map <- list(
  list(y = "breadth_day0",     x = "adpt_breadth_spike_pre",           title = "Day 0"),
  list(y = "breadth_day14",    x = "adpt_breadth_spike",               title = "Day 14"),
  list(y = "breadth_cd4_day0", x = "adpt_COVID_CD4_spike_breadth_pre", title = "Day 0"),
  list(y = "breadth_cd4_day14",x = "adpt_COVID_CD4_spike_breadth",     title = "Day 14"),
  list(y = "breadth_cd8_day0", x = "adpt_COVID_CD8_spike_breadth_pre", title = "Day 0"),
  list(y = "breadth_cd8_day14",x = "adpt_COVID_CD8_spike_breadth",     title = "Day 14")
)
plot_case_scatter <- function(data, x_var, y_var,
                              xlim = c(1e-5, 5e-3), ylim = c(1e-5, 5e-3),
                              xlab = NULL, ylab = NULL, title = NULL,
                              label_nd = TRUE, nd_value = 1e-5,
                              x_axis_color = "#cba69e") {
  
  gg <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(alpha = 0.95, size = .25)+ #aes(shape = case2)) +
    ggpubr::stat_cor(
      method = "spearman", size = 2, label.sep = "\n",
      label.x.npc = 1.0, label.y.npc = 0.1,
      hjust = 0, vjust = 0,
      show.legend = FALSE
    ) +
    scale_x_log10(
      labels = if (label_nd) function(x) log10_nd_labels(x, nd_value = nd_value)
      else scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      labels = if (label_nd) function(x) log10_nd_labels(x, nd_value = nd_value)
      else scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_shape_manual(values = c("C" = 19, "NC" = 1)) +
    theme_x +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    annotation_logticks(sides = "bl") +
    theme(
      axis.title.x = element_text(size = 6, colour = x_axis_color, margin = margin(t = 2)),
      legend.position = "none"
    ) + 
    geom_abline(linetype = "dashed")
  
  if (!is.null(xlab)) gg <- gg + xlab(xlab)
  if (!is.null(ylab)) gg <- gg + ylab(ylab)
  if (!is.null(title)) {
    gg <- gg + ggtitle(title) +
      theme(plot.title = element_text(size = 8, hjust = 0.5))
  }
  
  gg
}

# --- Build the six scatter panels ----------------------------------------
scatter_pt <- purrr::map(pair_map, ~ plot_case_scatter(
  adata_joined,
  x_var = .x$x, y_var = .x$y,
  xlab  = ylab_map[[.x$x]],
  ylab  = ylab_map[[.x$y]],
  title = .x$title
))

scatter_row1 <- scatter_pt[[1]] | scatter_pt[[2]] | scatter_pt[[3]] |
  scatter_pt[[4]] | scatter_pt[[5]] | scatter_pt[[6]]

cell_subsets_comparison = (row1 / row2 / scatter_row1 / new_vars_row1) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10))



# pdf("/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey/Manuscript/figures/Fig_S2_HA_Ts_v_Adaptive.pdf",
#     width = 12, height= 8)
# cell_subsets_comparison
# dev.off()
cell_subsets_comparison   # also draw to the active device
# ---- Save ---------------------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".png")),
  plot = cell_subsets_comparison, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi
)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".pdf")),
  plot = cell_subsets_comparison, width = fig_width, height = fig_height, units = "in"
)




