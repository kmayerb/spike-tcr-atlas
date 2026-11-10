# File: v1_Fig_4.R
# Note: refactored from 002_COV_BOOST_BREADTH.R

# ---- Output configuration ----------------------------------------------
# Fits a US-Letter page with 1-inch margins (8.5 - 1 = 7.5 in usable width;
# up to 9 in usable height).
out_dir    <- "Manuscript/figures_v1"   # <-- set destination folder
out_name   <- "v1_Fig_4"        # filename stem (no extension)
fig_width  <- 7.5                    # inches
fig_height <- 2.5                    # inches (must stay <= 9 for letter)
fig_dpi    <- 300
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

# ---- Canonical axis-title building blocks (match Fig 3) ----------------
lab_est   <- expression("HLA-intrinsic breadth estimate: " ~ E*"["*"HA-"*T[s]*"]")
lab_ratio <- expression("Ratio: HA-"*T[s]/E*"["*"HA-"*T[s]*"]")

## Load ADATA with metadata
require(dplyr)
require(ggplot2)
require(patchwork)
library(ggplot2)
library(ggpubr)
library(purrr)
library(ggbeeswarm) # For geom_quasirandom

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

######################### #################### ########################
# Load Delmonte data
r = estimate_dir
HLA_FOLDER = '/Volumes/fh/fast/gilbert_p/kmayerbl/spike_atlas_data/'
HLA_FILE   =  'delmonte_hla_predictions.csv'
HITS_FP     = file.path(r, "delmonte_spike_breadth_all_hits.csv") 
SCORES_FP   = file.path(r, "delmonte_spike_breadth_x1w2.csv") 
hla = readr::read_csv(file.path(HLA_FOLDER, HLA_FILE))
scores_fp   = file.path(r, "delmonte_spike_breadth_estimate_x1w2.csv")
scores      = readr::read_csv(scores_fp) 
hits        = readr::read_csv(HITS_FP )

# SIMULATE DELMONTED WITH 3006 Pools
pools_3006_mRNA = readRDS("pools_cov3006_mRNA.rds")
# Here we simulate what breadth is expected
sim_delmonte_3006 = simulate_hla_intrinsic_breadth(HLA = hla, POOLS=pools_3006_mRNA , N = 1000)
# Get quantile score
quantile_overall = sim_delmonte_3006  %>% mutate(name_1 = person_id) %>% left_join(scores, by ="name_1") %>% 
  mutate(
    # map2_dbl iterates row-by-row through breadth and sums simultaneously
    breadth_quantile = map2_dbl(breadth, sums, function(b, s) {
      # mean(s <= b) calculates what proportion of simulated sums are less than or equal to breadth
      mean(s <= b, na.rm = TRUE)
    })
  ) %>% 
  select(name_1, breadth_quantile)

predicted = sim_delmonte_3006  %>%
  select(person_id, sums) %>%
  tidyr::unnest(sums) %>% 
  group_by(person_id) %>% 
  summarise(estimated_breadth = median(sums))  %>% 
  ungroup() %>% 
  select(name_1 = person_id, estimated_breadth)


metadata = readr::read_csv('/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_ID159_tools/data/20260520_breadth_across_13_studies_with_imputed_hla_alleles_draft_metadata.csv')
names(metadata)
ddata = scores %>% 
  inner_join(metadata %>% filter(archive_name_x == "delmonte") %>% select(name_1, severity, exposure_level, icu, sex,death, days_since_onset) %>% 
               filter(exposure_level== "I") , by = 'name_1') %>% 
  filter(days_since_onset > 7) %>% 
  left_join(predicted) %>% 
  filter(estimated_breadth > 1E-5) %>% 
  mutate(severity = factor(severity, levels = c("Mild", "Moderate_severe","Critical")))



# # Plot 2: Severity vs. Estimated Breadth Boxplot
# p2 <- ddata %>% ggplot(aes(x = severity, y = estimated_breadth, fill = severity)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_quasirandom(size = .25) +
#   scale_y_log10() + 
#   scale_fill_manual(values = c("salmon", "red","brown")) +
#   xlab("") + # Removed the duplicate xlab from original code here
#   ylab(lab_est) + 
#   theme_x + 
#   theme(legend.position = "None") + 
#   stat_compare_means(size = 2.5,label = "p.signif",
#                      comparisons = list(c("Mild", "Moderate_severe"),c("Mild", "Critical"),c("Moderate_severe", "Critical")))
# 
# # Plot 3: Severity vs. Ratio Boxplot
# p3 <- ddata %>% ggplot(aes(x = severity, y = breadth/estimated_breadth, fill = severity)) +
#   geom_boxplot(outlier.shape = NA) + 
#   geom_quasirandom(size = .25) +
#   scale_y_sqrt(breaks = c(0.1,1,2,4,8)) + 
#   geom_hline(aes(yintercept = 1), linetype = "dashed") +
#   scale_fill_manual(values = c("salmon", "red","brown")) +
#   ylab(lab_ratio) + 
#   xlab("") +
#   theme_x + 
#   theme(legend.position = "None") + 
#   stat_compare_means(size = 2.5,label = "p.signif",
#                      comparisons = list(c("Mild", "Moderate_severe"),c("Mild", "Critical"),c("Moderate_severe", "Critical")))
# 
# 
# #cp "/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_ID159_tools/data/delmonte_metadata2.csv" \
# #"${PROJECT}/Manuscript/adata/"
# 
death = readr::read_csv(file.path(project_folder, "Manuscript", "adata", 'delmonte_metadata2.csv'))
death = death %>% select(name_1 = sample_name, Alive, ICU.admit, Gender, Age)

ddata2 = scores %>%
  inner_join(metadata %>% filter(archive_name_x == "delmonte") %>% select(name_1, severity, exposure_level, icu, sex,death, days_since_onset) %>%
               filter(exposure_level== "I") , by = 'name_1') %>%
  #filter(days_since_onset > 7) %>%  # we don't care about eventual breadth as we are just looking at predicted
  left_join(predicted) %>%
  filter(estimated_breadth > 1E-5) %>%
  mutate(severity = factor(severity, levels = c("Mild", "Moderate_severe","Critical")))
# 
# r1 = ddata2 %>% left_join(death) %>%
#   filter(icu == 'ICU') %>% 
#   filter(!is.na(Alive)) %>% 
#   filter(Alive != "not known") %>%
#   ggplot(aes(x = Alive, y = estimated_breadth, fill = Alive)) +
#   geom_boxplot(outlier.shape = NA) + 
#   geom_quasirandom(size = .5) +
#   scale_y_log10()+
#   scale_fill_manual(values = c("limegreen", "orange")) +
#   ylab(lab_est) + 
#   #stat_compare_means() + 
#   stat_compare_means(comparison = list(c("Alive", "Deceased")), size = 2.5,label = "p.signif") + 
#   facet_wrap(~sex) + 
#   theme_x+ 
#   theme(legend.position = "None")

# Plot 2: Severity vs. Estimated Breadth Boxplot
p2 <- ddata %>% ggplot(aes(x = severity, y = estimated_breadth, fill = severity)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size = .25) +
  #scale_y_log10(expand = expansion(mult = c(0.05, 0.15))) +   # extra headroom on top
  scale_fill_manual(values = c("salmon", "red","brown")) +
  xlab("") +
  ylab(lab_est) + 
  theme_x + 
  theme(legend.position = "None") + 
  stat_compare_means(size = 2.5, label = "p.signif",
                     step.increase = 0.12,   # more vertical gap between stacked brackets
                     comparisons = list(c("Mild", "Moderate_severe"),
                                        c("Mild", "Critical"),
                                        c("Moderate_severe", "Critical")))+scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                                                                         labels = scales::trans_format("log10", scales::math_format(10^.x)),
                                                                                         expand = expansion(mult = c(0.05, 0.15)))

# Plot 3: Severity vs. Ratio Boxplot
p3 <- ddata %>% ggplot(aes(x = severity, y = breadth/estimated_breadth, fill = severity)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(size = .25) +
  scale_y_sqrt(breaks = c(0.1,1,2,4,8),
               expand = expansion(mult = c(0.05, 0.15))) + 
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  scale_fill_manual(values = c("salmon", "red","brown")) +
  ylab(lab_ratio) + 
  xlab("") +
  theme_x + 
  theme(legend.position = "None") + 
  stat_compare_means(size = 2.5, label = "p.signif",
                     step.increase = 0.12,
                     comparisons = list(c("Mild", "Moderate_severe"),
                                        c("Mild", "Critical"),
                                        c("Moderate_severe", "Critical")))

# r1
r1 = ddata2 %>% left_join(death) %>%
  filter(icu == 'ICU') %>% 
  filter(!is.na(Alive)) %>% 
  filter(Alive != "not known") %>%
  ggplot(aes(x = Alive, y = estimated_breadth, fill = Alive)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(size = .5) +
  #scale_y_log10(expand = expansion(mult = c(0.05, 0.15))) +
  scale_fill_manual(values = c("limegreen", "orange")) +
  ylab(lab_est) + 
  stat_compare_means(comparisons = list(c("Alive", "Deceased")), size = 2.5, label = "p.signif") + 
  facet_wrap(~sex) + 
  theme_x + 
  theme(legend.position = "None")+xlab("")+scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                                          labels = scales::trans_format("log10", scales::math_format(10^.x)),
                                                          expand = expansion(mult = c(0.05, 0.15)))


# ---- Assemble: one row, tagged A–C ------------------------------------
final_plot <- (p2 + p3 + r1) +
  plot_layout(nrow = 1) +
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
