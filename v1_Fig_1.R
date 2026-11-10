

# Because this figure requires output of Cox PH models we need to
# Run some preliminary files first (only needs to be done once)
  # 1. Rscript ./Manuscript/v1_Fig_1_prep_cox_data.R
   # ----> Writes: /Manuscript/adata/id159_prop_hazard_adata.csv 
  # 2. Rscript ./Manuscript/v1_Fig_1_run_cox_models.R
    # ----> Writes: /Manuscript/adata/001a_cox_ph_models_dfall.csv
  # Once those inputs are created we can proceed to generate figure 1.



# ---- data -----------------------------------------------------------

# COX PROPORTIONAL HAZARDS 
project_folder = '/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey/'
source(file.path(project_folder, 'Manuscript', "v1_Fig_1_helper_functions.R"))

OUTPUT_DIR    <- "Manuscript/figures_v1"   # <-- set destination folder
out_name   <- "v1_Fig_1"            # filename stem (no extension)

stopifnot(file.exists(file.path(project_folder, 'Manuscript','adata', '001a_cox_ph_models_dfall.csv')))
stopifnot(file.exists(file.path(project_folder, 'Manuscript','adata' , "id159_prop_hazard_adata.csv")))
# <dfall>  contains model coefficient estimates
dfall <- read_csv(file.path(project_folder, 'Manuscript','adata', '001a_cox_ph_models_dfall.csv'))
# <df>  contains model coefficient estimates
df    <- read_csv(file.path(project_folder, 'Manuscript','adata' , "id159_prop_hazard_adata.csv"))


# ---- shared settings ------------------------------------------------
TS_MEDIAN  <- df$log2tcell %>% median()         # forced T_S median split used in the manuscript
PT_SIZE    <- 0.75

# =====================================================================
# Row 1
# =====================================================================

# Panel C -- cumulative incidence, T_S, full follow-up (strata 0)
C <- plot_cumulative_incidence(
  df, y_var = "log2tcell",
  strata_var   = "strata", strata_val = 0,
  legend_labs  = c(expression("Low HA-"*T[S]), expression("High HA-"*T[S] ) ),
  conf_int     = TRUE,
  force_median = TS_MEDIAN,
  color_values = c("darkgray", "navy")) + 
  guides(
    colour = guide_legend(reverse = FALSE),
    fill   = guide_legend(reverse = FALSE)
  ) 

C
# Panel D -- cumulative incidence, T_S, strata 1
D <- plot_cumulative_incidence(
  df, y_var = "log2tcell",
  strata_var   = "strata", strata_val = 1,
  legend_labs  =  c(expression("Low HA-"*T[S]), expression("High HA-"*T[S] ) ),
  conf_int     = TRUE,
  force_median = TS_MEDIAN,
  color_values = c("darkgray", "navy"))+
  guides(
    colour = guide_legend(reverse = FALSE),
    fill   = guide_legend(reverse = FALSE)
  ) 

D
# Panel E -- six-model forest / rope plot
forest_models <- c(
  "Surv(time, status) ~ age10 + male + high_log2tcell + strata(strata)",
  "Surv(time, status) ~ age10 + male + log2tcell + strata(strata)",
  "Surv(time, status) ~ age10 + male + log2elispot + strata(strata)",
  "Surv(time, status) ~ age10 + male + log2nt50 + strata(strata)",
  "Surv(time, status) ~ age10 + male + log2nt50 + log2elispot + log2tcell + strata(strata)",
  "Surv(time, status) ~ age10 + male + log2nt50_pre + log2tcell + log2elispot + strata(strata)"
)
forest_term_order <- c("High~T[s]", "High~ELISpot", "High~NT50",
                       "Log[2]~T[s]", "Log[2]~ELISpot", "Log[2]~NT50",
                       "Log[2]~NT50pre", "Age10Y", "Male")

E <- make_forest_panel(
  dfall, include_forms = forest_models, term_order = forest_term_order, add_hr = FALSE, add_hr_right = TRUE,
) +
  theme(axis.title.x = element_text(size = 8),
        axis.text    = element_text(size = 8),
        axis.text.y  = element_text(size = 6)) + 
  theme(strip.text.y = element_blank()) 


E

strip_labels <- c(
  "Surv(time, status) ~ age10 + male + high_log2tcell + strata(strata)"                         = "1",
  "Surv(time, status) ~ age10 + male + log2tcell + strata(strata)"                              = "2",
  "Surv(time, status) ~ age10 + male + log2elispot + strata(strata)"                            = "3",
  "Surv(time, status) ~ age10 + male + log2nt50 + strata(strata)"                               = "4",
  "Surv(time, status) ~ age10 + male + log2nt50 + log2elispot + log2tcell + strata(strata)"     = "5",
  "Surv(time, status) ~ age10 + male + log2nt50_pre + log2tcell + log2elispot + strata(strata)" = "6"
)
library(ggforce)
E <- E +
  facet_col(vars(fct_relevel(formula, forest_models)),
            scales = "free_y", space = "free",
            labeller = as_labeller(strip_labels)) +

  theme(strip.text       = element_text(size = 6, hjust = 0,
                                        margin = margin(2, 0, 1, 0)),
        strip.background = element_blank()) 
E

# E <- E +
#   facet_grid(fct_relevel(formula, forest_models) ~ "",
#              scales = "free_y", space = "free_y",
#              labeller = as_labeller(strip_labels)) +
#   theme(strip.text.y     = element_text(angle = 0, size = 7, vjust = 1, hjust = 0),
#         strip.background = element_blank())
# E
new_labs <- c(
  "High~T[s]"      = "Day~14~HA-T[S] > Median",
  "Log[2]~T[s]"    = "Day~14~Log[2]~(HA-T[S])",
  "Log[2]~ELISpot" = "Day~14~Log[2]~(SFU/M)",
  "Log[2]~NT50"    = "Day~28~Log[2]~(NT50)",
  "Log[2]~NT50pre" = "Day~0~Log[2]~(NT50)",
  "Age10Y"         = "Age~(per~10~yr)",
  "Male"           = "Male"
)

relabel <- function(x) {
  out <- ifelse(x %in% names(new_labs), new_labs[x], x)  # unlisted terms keep current label
  parse(text = out)                                       # keep plotmath (subscripts etc.)
}

E <- E + scale_y_discrete(labels = relabel, position = "left") +
  theme(plot.margin  = margin(0, 50, 0, 0),
        panel.spacing.y = unit(1, "pt"), 
        axis.text.y = element_text(size = 6))

E


# =====================================================================
# Row 2 -- marker boxplots (C vs NC)
# =====================================================================
F_p <- plot_marker_box(df, y_var = "log2tcell_pre",   color_values = c("darkgray", "darkgray"),
                       ylab = expression("Day 0  " ~ Log[2]("HA-"*T[S])),   ylim = c(-13, -9), pt_size = PT_SIZE, p_label_y = -9.5, p_label = "p.signif")
G_p <- plot_marker_box(df, y_var = "log2tcell",       color_values = c("darkgray", "navy"),
                       ylab = expression("Day 14  " ~ Log[2]("HA-"*T[S])),  ylim = c(-13, -9), pt_size = PT_SIZE, p_label = "p.signif")
H_p <- plot_marker_box(df, y_var = "log2elispot_pre", color_values = c("darkgray", "darkgray"),
                       ylab = expression("Day 0  " ~ Log[2](SFU/M)),  ylim = c(2, 10),   pt_size = PT_SIZE, p_label = "p.signif")
I_p <- plot_marker_box(df, y_var = "log2elispot",     color_values = c("darkgray", "orange"),
                       ylab = expression("Day 14  " ~ Log[2](SFU/M)), ylim = c(2, 10),   pt_size = PT_SIZE, p_label = "p.signif")
J_p <- plot_marker_box(df, y_var = "log2nt50_pre",    color_values = c("darkgray", "darkgray"),
                       ylab = expression("Day 0  " ~ Log[2](NT50)),   ylim = c(4, 13),   pt_size = PT_SIZE, p_label = "p.signif")
K_p <- plot_marker_box(df, y_var = "log2nt50",        color_values = c("darkgray", "red"),
                       ylab = expression("Day 28  " ~ Log[2](NT50)),  ylim = c(4, 13),   pt_size = PT_SIZE, p_label = "p.signif")

# =====================================================================
# Assemble
# =====================================================================
row1 <- ggpubr::ggarrange(
  C, D, E, ncol = 3, widths = c(1.5, 1.25, 2.0),
  labels = c("C", "D", "E"),
  font.label = list(size = 8, face = "plain")
)
row2 <- ggpubr::ggarrange(
  F_p, G_p, H_p, I_p, J_p, K_p, ncol = 6,
  labels = c("F", "G", "H", "I", "J", "K"),
  font.label = list(size = 8, face = "plain")
)
figure1 <- ggpubr::ggarrange(row1, row2, ncol = 1, heights = c(1.6, 1.0))
figure1

#OUTPUT_DIR = "/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey/Manuscript/figures/"
outfile <- file.path(OUTPUT_DIR,  paste0(out_name, ".pdf"))
pdf(outfile, width = 7.5, height = 5)
print(figure1)
dev.off()
message("Wrote ", outfile)

# ---- also write a PNG -----------------------------------------------
outfile_png <- file.path(OUTPUT_DIR, paste0(out_name, ".png"))
png(outfile_png, width = 7.5, height = 5, units = "in", res = 300)
print(figure1)
dev.off()
message("Wrote ", outfile_png)
