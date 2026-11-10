# File: v1_Fig3.R
# Note: previously see 002_COV_BOOST_BREADTH

# ---- Output configuration ----------------------------------------------
# Fits a US-Letter page with 1-inch margins (8.5 - 1 = 7.5 in usable width;
# up to 9 in usable height).
out_dir    <- "Manuscript/figures_v1"   # <-- set destination folder
out_name   <- "v1_Fig_3"        # filename stem (no extension)
fig_width  <- 7.5                    # inches
fig_height <- 4.5                    # inches (must stay <= 9 for letter)
fig_dpi    <- 300

## Load ADATA with metadata
require(dplyr)
require(ggplot2)
require(patchwork)

# ---- Load data ---------------------------------------------------------
project_folder <- "/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey"
estimate_dir   <- "/Volumes/fh/fast/gilbert_p/kmayerbl/spike_atlas_data/estimate_20260527"

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

# ---- Canonical axis-title building blocks ------------------------------
# One source of truth for every axis label, so the figure reads
# consistently: subject term is always "HA-"*T[s] (lowercase subscript),
# the colon lives inside the quoted string, and E[...] uses literal
# brackets to denote the expectation.
lab_est   <- expression("HLA-intrinsic breadth estimate: " ~ E*"["*"HA-"*T[s]*"]")
lab_obs   <- expression("Observed breadth: HA-"*T[s])
lab_ratio <- expression("Ratio: HA-"*T[s]/E*"["*"HA-"*T[s]*"]")
lab_quant <- expression("Quantile: HA-"*T[s])

## Functions
# LET'S LOOK AT ENDOGENOUS
get_breadth_components <- function(HITS,
                                   HLA){
  
  locus_binary = readRDS('locus_binary.rds')
  alleles = c("A_0101", "A_0201", "A_0301", "A_1101", "A_2301", "A_2402", 
              "A_2601", "A_2902", "A_3201", "A_6801", "B_0702", "B_0801", "B_1402", 
              "B_1501", "B_1801", "B_2705", "B_3501", "B_4001", "B_4402", "B_4403", 
              "B_5101", "B_5701", "C_0202", "C_0304", "C_0401", "C_0501", "C_0602", 
              "C_0701", "C_0702", "C_0802", "C_1203", "C_1601", "DPA1_0103__DPB1_0101", 
              "DPA1_0103__DPB1_0201", "DPA1_0103__DPB1_0301", "DPA1_0103__DPB1_0401", 
              "DPA1_0103__DPB1_0402", "DPA1_0201__DPB1_0101", "DPA1_0201__DPB1_1101", 
              "DPA1_0202__DPB1_0101", "DQA1_0101__DQB1_0501", "DQA1_0102__DQB1_0602", 
              "DQA1_0102__DQB1_0604", "DQA1_0103__DQB1_0603", "DQA1_0201__DQB1_0202", 
              "DQA1_0301__DQB1_0302", "DQA1_0401__DQB1_0402", "DQA1_0501__DQB1_0201", 
              "DQA1_0505__DQB1_0301", "DRB1_0101", "DRB1_0301", "DRB1_0401", 
              "DRB1_0402", "DRB1_0404", "DRB1_0701", "DRB1_0801", "DRB1_1101", 
              "DRB1_1301", "DRB1_1302", "DRB1_1501")
  
  hits = HITS %>% filter(hla_detect == 1) %>% 
    left_join(locus_binary) %>%
    mutate(tcellsub = if_else(locus %in% c("A","B","C"), "CD8", "CD4")) %>%
    select(name_1, binary, total_unique_1, x1w2, locus, tcellsub, vfamcdr3) %>% 
    mutate(br = x1w2/total_unique_1) %>% 
    ungroup()%>%
    select(name_1, vfamcdr3, binary, locus, tcellsub, br)
  print(hits)
  # WE NEED TO ENSURE COMPLETENESS IF WE WANT TO COUNT ZEROS IN THE DISTRIBUTION
  hla        =  HLA %>% 
    filter(binary %in% locus_binary$binary) %>%
    mutate(hla_value = pred) %>%
    filter(hla_value == 1) %>%  
    select(name_1 = sample_id, binary, hla_value) %>% 
    arrange(name_1) %>% 
    left_join(locus_binary) %>%
    mutate(tcellsub = if_else(locus %in% c("A","B","C"), "CD8", "CD4"))
  print(hla)
  # THIS IS THE COMPLETE SET OF HLA
  zero_merge_hits = hla %>% 
    left_join(hits, by = c("name_1", "binary", "tcellsub","locus")) %>% 
    mutate(br = if_else(is.na(br), 1E-6, br)) %>% 
    mutate(br = if_else(br < 1E-6, 1E-6, br))
  
  # NOW COMPUTE A PER INDIVIDUAL SCORE
  zero_merge_hits_per_person = zero_merge_hits %>% 
    group_by(name_1, binary, locus, tcellsub) %>%
    summarise(breadth = sum(br))
  df =zero_merge_hits_per_person 
  pools <- split(df$breadth, df$binary)
  return(list('pools'=pools, 'per_allele_df'= df))
}


## TA/E[TA]
simulate_hla_intrinsic_breadth <- function(HLA, POOLS, N = 500, SEED =609){
  locus_binary = readRDS('locus_binary.rds')
  alleles = c("A_0101", "A_0201", "A_0301", "A_1101", "A_2301", "A_2402", 
              "A_2601", "A_2902", "A_3201", "A_6801", "B_0702", "B_0801", "B_1402", 
              "B_1501", "B_1801", "B_2705", "B_3501", "B_4001", "B_4402", "B_4403", 
              "B_5101", "B_5701", "C_0202", "C_0304", "C_0401", "C_0501", "C_0602", 
              "C_0701", "C_0702", "C_0802", "C_1203", "C_1601", "DPA1_0103__DPB1_0101", 
              "DPA1_0103__DPB1_0201", "DPA1_0103__DPB1_0301", "DPA1_0103__DPB1_0401", 
              "DPA1_0103__DPB1_0402", "DPA1_0201__DPB1_0101", "DPA1_0201__DPB1_1101", 
              "DPA1_0202__DPB1_0101", "DQA1_0101__DQB1_0501", "DQA1_0102__DQB1_0602", 
              "DQA1_0102__DQB1_0604", "DQA1_0103__DQB1_0603", "DQA1_0201__DQB1_0202", 
              "DQA1_0301__DQB1_0302", "DQA1_0401__DQB1_0402", "DQA1_0501__DQB1_0201", 
              "DQA1_0505__DQB1_0301", "DRB1_0101", "DRB1_0301", "DRB1_0401", 
              "DRB1_0402", "DRB1_0404", "DRB1_0701", "DRB1_0801", "DRB1_1101", 
              "DRB1_1301", "DRB1_1302", "DRB1_1501")
  
  hla        =  HLA %>% 
    filter(binary %in% locus_binary$binary) %>%
    mutate(hla_value = pred) %>%
    filter(hla_value == 1) %>%  
    select(name_1 = sample_id, binary, hla_value) %>% 
    arrange(name_1) %>% 
    left_join(locus_binary) %>%
    mutate(tcellsub = if_else(locus %in% c("A","B","C"), "CD8", "CD4"))
  
  # Now we draw based on each person
  cat(sprintf("Completeness: %.2f%% of alleles are present in the POOLS dictionary.\n", 
              mean(alleles %in% names(POOLS)) * 100))
  #print(POOLS)
  persons <-hla %>%
    ungroup() %>%
    group_by(name_1) %>%
    summarise(pool_set = list(unique(binary))) %>%
    rename(person_id = name_1)
  
  sample_person <- function(pool_set, N = N, pools, seed = SEED) {
    set.seed(seed)
    replicate(N, {
      sum(sapply(pool_set, function(p) sample(pools[[p]], size = 1)))
    })
  }
  
  persons2 <- persons %>%
    rowwise() %>%
    mutate(
      sums = list(sample_person(pool_set, N, pools= POOLS))
    ) %>%
    ungroup()
  return(persons2)
}


score_observed_quantiles <- function(pool_set, pool_values, pools) {
  # Safety check to ensure inputs line up
  if (length(pool_set) != length(pool_values)) {
    stop("Error: 'pool_set' and 'pool_values' must be the same length.")
  }
  
  # Calculate the quantile for each observed value
  quantiles <- sapply(seq_along(pool_set), function(i) {
    p <- pool_set[i]             # The specific pool name/index
    obs_val <- pool_values[i]    # The observed value for that pool
    pool <- pools[[p]]           # The background distribution for that pool
    
    # Compute the quantile score
    mean(pool <= obs_val)
  })
  
  # Return the composite score (e.g., the average of their quantiles)
  # You can remove mean() if you just want the raw vector of scores back
  mean(quantiles, na.rm = TRUE) 
}

### FIRST WE DEFINE REFERENCE FOR THE EXPECTATION FROM PREVIOUSL mRNA VACCINATION DATA (Fiore-Gartland et al. 2026)
### LETS GENERATE COV3006
HLA_FOLDER = '/Volumes/fh/fast/gilbert_p/kmayerbl/spike_atlas_data/'
HLA_FILE   =  'cov3006_hla_predictions.csv'
HITS_FP     = file.path(estimate_dir, "cov3006_spike_breadth_all_hits.csv") 
hla = readr::read_csv(file.path(HLA_FOLDER, HLA_FILE))
hits        = readr::read_csv(HITS_FP )
# WE NEED TO SUBSET TO 3006 V Exposure
metadata = readr::read_csv('/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_ID159_tools/data/20260520_breadth_across_13_studies_with_imputed_hla_alleles_draft_metadata.csv')
cov3006_vaccine_name1 = metadata %>% filter(archive_name_x == 'cov3006') %>% 
  filter(exposure_level %in% c('V','H')) %>% pull(name_1)
hits = hits %>% filter(name_1 %in% cov3006_vaccine_name1 )
hla = hla %>% filter(sample_id %in% cov3006_vaccine_name1) %>% 
  mutate(hla_value = as.numeric(pred)) %>% 
  mutate(pred = as.numeric(pred)) %>%
  select(sample_id, binary, hla_value, pred)
hits$name_1 %>% unique() %>% length()
hla$sample_id %>% unique() %>% length()
px_cov3006 = get_breadth_components(HITS = hits, 
                                    HLA = hla)
saveRDS(px_cov3006$pools, "pools_cov3006_mRNA.rds")
# Load reference breadth values
pools_3006_mRNA = readRDS("pools_cov3006_mRNA.rds")

# NOW WE WILL COMPUTE EXPECTATON FOR THE COV-BOOST
### ENSURE WE ONLY LOOK AT POST VAC
HLA_FOLDER = '/Volumes/fh/fast/gilbert_p/kmayerbl/spike_atlas_data/'
HLA_FILE   =  'cov_boost_hla_genotypes.csv'
r = '/Volumes/fh/fast/gilbert_p/kmayerbl/spike_atlas_data/estimate_20260527/'
HITS_FP     = file.path(r, "cov_boost_day0_day14_genotyped_spike_breadth_all_hits.csv") 
hla = readr::read_csv(file.path(HLA_FOLDER, HLA_FILE))
scores = readr::read_csv(file.path(r, "cov_boost_day0_day14_genotyped_spike_breadth_estimate_x1w2.csv" ))
hits        = readr::read_csv(HITS_FP )
adata       = readr::read_csv(file.path(project_folder, "Manuscript","adata","001_COV_BOOST_ADATA.csv"))
day14_name_1 = adata %>% filter(timepoint == "Day 14") %>% pull(name_1)
stopifnot(length(day14_name_1) == 120)
hla = hla %>% filter(sample_id %in% day14_name_1)
hits = hits %>% filter(sample_id %in% day14_name_1)
# SIMULATION STEP
require(purrr)
sim_covb_3006 = simulate_hla_intrinsic_breadth(HLA = hla, POOLS=pools_3006_mRNA , N = 1000)
quantile_overall = sim_covb_3006 %>% mutate(name_1 = person_id) %>% left_join(scores, by ="name_1") %>% 
  mutate(
    # map2_dbl iterates row-by-row through breadth and sums simultaneously
    breadth_quantile = map2_dbl(breadth, sums, function(b, s) {
      # mean(s <= b) calculates what proportion of simulated sums are less than or equal to breadth
      mean(s <= b, na.rm = TRUE)
    })
  ) %>% 
  select(name_1, breadth_quantile)

predicted = sim_covb_3006 %>%
  select(person_id, sums) %>%
  tidyr::unnest(sums) %>% 
  group_by(person_id) %>% 
  summarise(estimated_breadth = median(sums))  %>% 
  ungroup() %>% 
  select(name_1 = person_id, estimated_breadth)

comp = scores %>% filter(name_1 %in% day14_name_1) %>% 
  left_join(predicted, by = "name_1") %>%
  left_join(adata, by = 'name_1') %>% 
  mutate(ratio = breadth / estimated_breadth) %>%
  mutate(logratio = log2(breadth) - log2(estimated_breadth)) %>%
  mutate(abs = breadth - estimated_breadth) %>% 
  left_join(quantile_overall)
saveRDS(comp, "comp_covboost.rds")
# 1. Main Scatter Plot
p_scatter <- comp %>%
  ggplot(aes(x = estimated_breadth, y = breadth)) + 
  geom_point(aes(col = case), size = .25) + 
  geom_abline() + 
  scale_y_log10() + 
  scale_x_log10() +
  coord_cartesian(ylim = c(1E-4, 5E-3), xlim = c(1E-4, 5E-3) ) +
  theme_x + 
  theme(legend.position = "None") + 
  scale_color_manual("",values = c("red","darkgray"))+
  scale_fill_manual("",values = c("red","darkgray"))+ 
  xlab(lab_est)+ 
  ylab(lab_obs)

p_scatter

# 2. Top Density Plot (x-axis margin)
p_density_top <- comp %>%
  ggplot(aes(x = estimated_breadth, fill = case, col = case)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +
  coord_cartesian(xlim = c(1E-4, 5E-3)) +
  theme_void() + 
  theme(legend.position = "None") + 
  scale_color_manual("",values = c("red","darkgray"))+
  scale_fill_manual("",values = c("red","darkgray"))
# 3. Right Density Plot (y-axis margin)
# Notice we map to `y` instead of `x` to rotate the density curve
p_density_right <- comp %>%
  ggplot(aes(y = breadth, fill = case, col = case)) +
  geom_density(alpha = 0.4) +
  scale_y_log10() +
  coord_cartesian(ylim = c(1E-4, 5E-3)) +
  theme_void() + 
  theme(legend.position = "None")+
  scale_color_manual(values = c("red","darkgray"))+
  scale_fill_manual(values = c("red","darkgray"))

# 4. Combine them using patchwork
# We add an empty plot_spacer() to the top-right corner to keep things aligned
pall = p_density_top + plot_spacer() + 
  p_scatter     + p_density_right + 
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(10, 1),   # Makes the scatter plot 4x wider than the right density
    heights = c(1, 10)   # Makes the scatter plot 4x taller than the top density
  ) 


# --- Four summary panels as named objects ------------------------------
p_est <- comp %>%
  ggplot(aes(x = case, y = estimated_breadth)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(aes(col = case), size = .25) +
  scale_color_manual(values = c("red","darkgray")) +
  scale_y_log10() +
  ggpubr::stat_compare_means(label = "p.signif", comparison = list(c("case","noncase")),size = 2.5) +
  facet_wrap(~timepoint) +
  coord_cartesian(ylim = c(1E-4, 5E-3)) + theme_x + theme(legend.position = "none") +
  ylab(lab_est) + xlab("")

p_obs <- comp %>%
  ggplot(aes(x = case, y = breadth)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(aes(col = case), size = .25) +
  scale_color_manual(values = c("red","darkgray")) +
  scale_y_log10() +
  ggpubr::stat_compare_means(label = "p.signif", comparison = list(c("case","noncase")),size = 2.5) +
  facet_wrap(~timepoint) +
  coord_cartesian(ylim = c(1E-4, 5E-3)) + theme_x + theme(legend.position = "none") +
  ylab(lab_obs) + xlab("")

p_ratio <- comp %>%
  ggplot(aes(x = case, y = ratio)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(aes(col = case), size = .25) +
  scale_color_manual(values = c("red","darkgray")) +
  ggpubr::stat_compare_means(label = "p.signif", comparison = list(c("case","noncase")), size = 2) +
  facet_wrap(~timepoint) + theme_x + theme(legend.position = "none") +
  ylab(lab_ratio) + xlab("")

p_quant <- comp %>%
  ggplot(aes(x = case, y = breadth_quantile)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(aes(col = case),size = .25) +
  scale_color_manual(values = c("red","darkgray")) +
  ggpubr::stat_compare_means(label = "p.signif", comparison = list(c("case","noncase")), size = 2) +
  facet_wrap(~timepoint) + theme_x + theme(legend.position = "none") +
  ylab(lab_quant) + xlab("")

# Final Plot --- 
# Use patchwork to assemble the final plot
# combine pall with final plot components using patchwork where width 
# of pall should be 2x the size of next panels 

# --- Final 5-column layout: pall (wide) + four summary panels ----------
row1 <- wrap_elements(full = pall) +
  p_est + p_obs + p_ratio + p_quant +
  plot_layout(
    ncol   = 5,
    widths = c(3, 1, 1, 1, 1)   # first column 3x the width of each other
  ) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "plain", size = 8))


#  Load the elispot data
df = readr::read_csv(file.path(project_folder, 'Manuscript','adata','001_COV_BOOST_ADATA_WIDE.csv'))
elispot_data = df %>% select(name_1, log2elispot) %>% # NOTE WE DO THIS TO GET THE 1 IMPUTED VALUE (IMPUTATION WAS DONE IN LOG2 SPACE)
  mutate(elispot = 2^log2elispot)
adata14 = comp %>% left_join(elispot_data)


# Define the color scale once so all panels share an identical legend
# (lowercase T[s] subscript to match every axis title)
ratio_scale <- scale_color_viridis_c(expression(frac("HA-"*T[s], E*"["*"HA-"*T[s]*"]")))
p1 <- ggplot(adata14, aes(x = estimated_breadth, y = breadth, col = ratio)) +
  geom_point(size = .25) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(linetype = "dashed") +
  coord_cartesian(xlim = c(1E-4, 5E-3), ylim = c(1E-4, 5E-3)) +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho", label.y.npc = 1,size = 2) +
  ratio_scale +
  ylab(expression("Day 14 observed breadth: HA-"*T[s])) +
  xlab(lab_est) +
  theme_x

p2 <- ggplot(adata14, aes(x = estimated_breadth, y = elispot, col = ratio)) +
  geom_point(size = .25) +
  scale_y_log10() +
  scale_x_log10() +
  ratio_scale +
  ylab("Day 14 ELISpot (SFU/M)") +
  xlab(lab_est) +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho", label.y.npc = 1,size = 2) +
  theme_x

p3 <- ggplot(adata14, aes(x = breadth, y = elispot, col = ratio)) +
  geom_point(size = .25) +
  scale_y_log10() +
  scale_x_log10() +
  ratio_scale +
  ylab("Day 14 ELISpot (SFU/M)") +
  xlab(expression("Day 14 observed breadth: HA-"*T[s])) +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho", label.y.npc = 1,size = 2) +
  theme_x

p4 <- ggplot(adata14, aes(x = ratio, y = elispot, col = ratio)) +
  geom_point(size = .25) +
  scale_y_log10() +
  ratio_scale +
  ylab("Day 14 ELISpot (SFU/M)") +
  xlab(lab_ratio) +
  geom_vline(aes(xintercept = 1), linetype = "dashed") +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho", label.y.npc = 1, size = 2) +
  theme_x

# Assemble: one row, collect guides into a single shared legend
row2 <- (p1 | p2 | p3 | p4) +
  plot_layout(guides = "collect") &
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 5),
    legend.title     = element_text(size = 6),
    legend.key.width  = unit(12, "pt"),   # colourbar length
    legend.key.height = unit(4,  "pt"),   # colourbar thickness (thinner = less prominent)
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, -4, 0),  # pull the plots up toward the bar
    legend.ticks      = element_blank()       # ggplot >= 3.5; drop tick marks
  )



# --- Final 5-column layout: pall (wide) + four summary panels ----------
row1 <- wrap_elements(full = pall) +
  p_est + p_obs + p_ratio + p_quant +
  plot_layout(
    ncol   = 5,
    widths = c(3, 1, 1, 1, 1)   # first column 3x the width of each other
  ) 

final_plot <- (row1 / row2) +
  plot_layout(heights = c(1, 0.75)) +
  plot_annotation(tag_levels = "A") &
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