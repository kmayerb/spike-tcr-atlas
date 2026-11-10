# 001_COV_BOOST_ADATA.R
# Script -- ASSEMBLES COV-BOOST ADATA FRAMEWORK BEFORE APPENDING DERIVED
# TCR METRICS
require(dplyr)
require(ggplot2)
project_folder = '/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey'
write_files = TRUE# SOURCE() will overwrite files if TRUE
# Oxford metadata file
# Get the original metadata file 
# NOTE: -- RECALL THAT SOME DATA WAS FLAGGED DUE TO MATERIAL TRANSFER EVENT OR
# FAILUER TO MEET MINIMAL CLONAL OVERLAP WITH PAIRED LONGITUDINAL SAMPLE
# NOTE: See flag (included.in.analysis == 1), if either Day 0 or Day 14 samples are compromised Ptid, is excluded
# NOTE: Only samples with complete Day 0, Day 14 TCR Sample Pairing and HLA-Genotypes are retained
metadata = readr::read_tsv('/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey/metadata/metadata_30092025.tsv')
cx_boost = metadata %>% filter(included.in.analysis == 1) %>% # Remove samples that were excluded by Adaptive QC. Material transfer notes, sample quality, Lack of overlapping serial clones raise risk of sample mislabeling
  filter(Timepoint %in% c("Day 0","Day 14")) %>% 
  filter(Trial == "Cov-Boost") %>%  
  filter(A1 != "XXXX") %>% # Remove sample with missing HLA data (see OXFMV151_TCRB,1130049)
  select(`PID Group`, sc_repeat_pid, `PID Group` ,`Pseudo-ID`, Timepoint, Trial, Status) %>% 
  group_by(`PID Group` , Timepoint) %>% 
  slice(1) %>% 
  group_by(`PID Group`,sc_repeat_pid ) %>% 
  tally() %>% filter(n >=2) %>% 
  pull(`PID Group`) %>% 
  unique()

# Check against the expected full valid list
cohort_boost = c("PID001", "PID002", "PID003", "PID004", "PID005", "PID006", 
                 "PID007", "PID008", "PID009", "PID011", "PID012", "PID013", "PID014", 
                 "PID015", "PID016", "PID017", "PID018", "PID019", "PID020", "PID021", 
                 "PID022", "PID023", "PID024", "PID025", "PID026", "PID027", "PID028", 
                 "PID029", "PID030", "PID031", "PID032", "PID033", "PID034", "PID035", 
                 "PID036", "PID037", "PID038", "PID039", "PID040", "PID041", "PID042", 
                 "PID043", "PID044", "PID045", "PID046", "PID047", "PID048", "PID049", 
                 "PID050", "PID051", "PID052", "PID054", "PID055", "PID056", "PID057", 
                 "PID059", "PID060", "PID061", "PID062", "PID063", "PID066", "PID067", 
                 "PID068", "PID069", "PID070", "PID071", "PID072", "PID073", "PID074", 
                 "PID075", "PID076", "PID077", "PID078", "PID079", "PID080", "PID081", 
                 "PID082", "PID083", "PID084", "PID085", "PID087", "PID088", "PID089", 
                 "PID090", "PID091", "PID092", "PID093", "PID094", "PID095", "PID096", 
                 "PID097", "PID098", "PID099", "PID100", "PID101", "PID102", "PID103", 
                 "PID104", "PID105", "PID106", "PID107", "PID108", "PID109", "PID111", 
                 "PID113", "PID114", "PID115", "PID116", "PID117", "PID118", "PID119", 
                 "PID120", "PID121", "PID122", "PID123", "PID125", "PID126", "PID127", 
                 "PID128", "PID129", "PID148")



cx_boost[!cx_boost %in% cohort_boost]     # 1090058 -- Day 0 (OXFMV195_TCRB) was missing from .zip file so was not paired in version 1, found in Cov2 zip file but with only 15 templates, Day 14 (OXFMV182_TCRB)
cohort_boost[!cohort_boost %in% cx_boost] # 1130049 -- Day 14 missing HLA information, (OXFMV151_TCRB) full dataset is defined as experimentally genotyped
analysis_cohort = intersect(cx_boost, cohort_boost)
stopifnot(length(analysis_cohort) == 120)
unique(analysis_cohort) %>% length()

# 1. Dictionary for the First Component
dict_first_component <- c(
  "2*ChAd+Ad26"     = "ChAdOx",
  "2*ChAd+BNT"      = "ChAdOx",
  "2*BNT+MOD-half"  = "mRNA",
  "2*ChAd+MOD-half" = "ChAdOx",
  "3*BNT"           = "mRNA",
  "2*BNT+Ad26"      = "mRNA",
  "2*BNT+BNT-half"  = "mRNA",
  "2*ChAd+BNT-half" = "ChAdOx",
  "2*BNT+NVX-half"  = "mRNA",
  "2*BNT+NVX"       = "mRNA",
  "2*BNT+ChAd"      = "mRNA",
  "3*ChAd"          = "ChAdOx",
  "2*ChAd+NVX"      = "ChAdOx",
  "2*ChAd+NVX-half" = "ChAdOx",
  "2*BNT+MOD"       = "mRNA",
  "2*ChAd+MOD"      = "ChAdOx"
)

# 2. Dictionary for the Second Component
dict_second_component <- c(
  "2*ChAd+Ad26"     = "Ad26",
  "2*ChAd+BNT"      = "mRNA",
  "2*BNT+MOD-half"  = "mRNA",
  "2*ChAd+MOD-half" = "mRNA",
  "3*BNT"           = "mRNA",
  "2*BNT+Ad26"      = "Ad26",
  "2*BNT+BNT-half"  = "mRNA",
  "2*ChAd+BNT-half" = "mRNA",
  "2*BNT+NVX-half"  = "Protein subunit",
  "2*BNT+NVX"       = "Protein subunit",
  "2*BNT+ChAd"      = "ChAdOx",
  "3*ChAd"          = "ChAdOx",
  "2*ChAd+NVX"      = "Protein subunit",
  "2*ChAd+NVX-half" = "Protein subunit",
  "2*BNT+MOD"       = "mRNA",
  "2*ChAd+MOD"      = "mRNA"
)
dict_second_product <- c(
  "2*ChAd+Ad26"     = "Ad26",
  "2*ChAd+BNT"      = "BNT",
  "2*BNT+MOD-half"  = "MOD",
  "2*ChAd+MOD-half" = "MOD",
  "3*BNT"           = "BNT",
  "2*BNT+Ad26"      = "Ad26",
  "2*BNT+BNT-half"  = "BNT",
  "2*ChAd+BNT-half" = "BNT",
  "2*BNT+NVX-half"  = "NVX",
  "2*BNT+NVX"       = "NVX",
  "2*BNT+ChAd"      = "ChAdOx",
  "3*ChAd"          = "ChAdOx",
  "2*ChAd+NVX"      = "NVX",
  "2*ChAd+NVX-half" = "NVX",
  "2*BNT+MOD"       = "MOD",
  "2*ChAd+MOD"      = "MOD"
)

adata_part_1 = metadata %>% 
  mutate(sc_repeat_pid = `PID Group`) %>% # Note that we will use `PID Group` as sc_repeat_id for publication!
  filter(Timepoint %in% c("Day 0", "Day 14")) %>% 
  filter(Trial == "Cov-Boost") %>%
  mutate(name_1 = paste0(`Pseudo-ID`, "_TCRB")) %>% 
  filter(sc_repeat_pid %in% analysis_cohort) %>% 
  group_by(name_1) %>% slice(1) %>%
  filter(name_1 != 'OXFMV122_TCRB') %>% # note OXFMV178_TCRB and OXFMV122_TCR are both designated the same but OXFM178_TCRB is in repertoire .zip  #  filter(sc_repeat_pid ==  1090023 )  %>%  filter(name_1 %in% c('OXFMV122_TCRB','OXFMV178_TCRB')) %>% View()
  arrange(sc_repeat_pid, name_1) %>% 
  mutate(case = if_else(Status == "Noncase", "noncase", "case")) %>%
  mutate(event = if_else(case == "case", 1, 0)) %>%
  mutate(strata = if_else( `Week of vaccination`< 30, 0,1) ) %>%
  mutate(age10 = Age/10) %>%
  mutate(male = if_else(Sex == "Male",1,0)) %>% 
  mutate(prime = dict_first_component[Arm2]) %>% 
  mutate(boost = dict_second_component[Arm2]) %>% 
  mutate(boost_product = dict_second_product[Arm2]) %>% 
  mutate(chadvec_only = if_else(Arm2 == "3*ChAd", 1, 0)) %>%
  mutate(nt50na = is.na(`PNA(IU/ml)`) ) %>%
  mutate(elispotna = is.na(ELISpot) ) %>%
  mutate(nt50 = as.numeric(`PNA(IU/ml)`)) %>% 
  mutate(nt50 = if_else(is.na(nt50), 1, nt50)) %>% # Note "<6" is now being set to 1, NA will be imputed
  mutate(nt50 = if_else(nt50na, NA, nt50 )) %>%
  #mutate(strata = if_else())
  select(name_1, sc_repeat_pid, 
         timepoint = Timepoint,
         case, 
         event,
         time = follow_up_days,
         strata,
         sex = Sex, 
         male, 
         age = Age,
         age10,
         Arm2,
         platform,
         prime,
         boost,
         boost_product,
         chadvec_only,
         elispot = ELISpot,
         nt50,
         nt50na,
         elispotna,
         `PNA(IU/ml)`,
         week= `Week of vaccination`,
         vac.to.infection.day,
         `included.in.analysis`,
         Status,
         last_vac = last.vac.before.sample) 


# ROUND AGE FOR PUBLIC DATA -- shared in 'adata_v1'
adata_part_1 = adata_part_1 %>% mutate(age = round(age, -1), age10 = round(age10))

if (write_files){
  adata_part_1 %>% 
    write.csv(file.path(project_folder , 'Manuscript', 'adata_v1', '001_COV_BOOST_ADATA.csv'), row.names= F)
  saveRDS(adata_part_1, file.path(project_folder , '001_COV_BOOST_ADATA.rds'))
}


# We need to replace any zeros with 1 before log-transformation 
adata_part_1  = adata_part_1 %>% 
  mutate(elispot = ifelse(elispot == 0, 1, elispot)) %>% # replace 0 with 1 for log transform
  mutate(nt50 = ifelse(nt50 == 0, 1, nt50)) %>% # replace 0 with 1 for log transform
  mutate(log2elispot = log2(elispot)) %>%
  mutate(log2nt50    = log2(nt50))

# Imputation of missing values (PRE DATA)
# We first select only the pre data
adata_pre_raw = adata_part_1 %>% 
  filter(timepoint == "Day 0") 

adata_pre_raw %>% filter(is.na(nt50))
adata_pre_raw %>% filter(is.na(log2nt50))
adata_pre_raw %>% filter(is.na(elispot))
adata_pre_raw %>% filter(is.na(log2elispot))

library(VIM)
require(gower)
# Is there any missing nt50 data
set.seed(1)
k_impute = 9
adata_pre_raw$prime = factor(as.vector(adata_pre_raw$prime))
imputed_df_pre_nt50 <- VIM::kNN(
  adata_pre_raw, 
  variable = "log2nt50", 
  dist_var = c("age","prime"),
  weights = c(1, 1),
  k = k_impute
)
imputed_df_pre_nt50 %>% filter(log2nt50_imp == TRUE)

# Imputation of missing values (POST DATA)
adata_raw = adata_part_1 %>% 
  filter(timepoint == "Day 14") 
adata_raw %>% filter(is.na(nt50))
adata_raw %>% filter(is.na(log2nt50))
adata_raw %>% filter(is.na(elispot))
adata_raw %>% filter(is.na(log2elispot))

library(VIM)
require(gower)
adata_raw$platform = factor(adata_raw$platform)
k_impute = 9
imputed_df <- VIM::kNN(
  adata_raw, 
  variable = "log2nt50", 
  dist_var = c("age", "platform"),
  weights = c(1, 1),
  k = k_impute
)

target_row <- adata_raw[52, c("age", "platform")]
# 3. Identify the "Candidates" (the rows that ARE NOT missing data)
candidates <- adata_raw[!is.na(adata_raw$log2nt50), c("age", "platform")]
# 4. Find the top 9 neighbors
neighbors <- gower_topn(target_row, candidates, n = 9,  weights = c(1, 1))
candidates[neighbors$index[,1],]
# NOW ADD ELISPOT
imputed_df <- VIM::kNN(
  imputed_df, # we just reuse imputed_df and add more
  variable = "log2elispot", 
  dist_var = c("age", "platform"), 
  k = k_impute
)



# We now replace with the impute values
adata_pre_raw$name_1 == imputed_df_pre_nt50$name_1
adata_pre_raw$log2nt50 = imputed_df_pre_nt50$log2nt50
adata_pre_raw$log2nt50_imp =  imputed_df_pre_nt50$log2nt50_imp
adata_raw$log2elispot_imp =  FALSE

# We now replace with the impute values
stopifnot(all(adata_raw$name_1 == imputed_df$name_1))
adata_raw$log2nt50 = imputed_df$log2nt50
adata_raw$log2nt50_imp =  imputed_df$log2nt50_imp
adata_raw$log2elispot = imputed_df$log2elispot
adata_raw$log2elispot_imp =  imputed_df$log2elispot_imp

# Check combine
bind_rows(adata_pre_raw, adata_raw) %>% filter(log2nt50_imp == TRUE | log2elispot_imp == TRUE) %>%
  select(name_1, timepoint, nt50, elispot, log2elispot, log2elispot_imp, log2nt50, log2nt50_imp)

# Combine pre and post into a long form 
adata_long = bind_rows(adata_pre_raw, adata_raw)
# Combine pre and post into a wide
adata_wide = adata_raw %>% left_join(adata_pre_raw, 
                                      by = c("sc_repeat_pid", "platform", "case", "event", "strata", "time"), 
                                      suffix = c("","_pre"))


weight_df = readr::read_csv('/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/02_ID159_Voysey/000_weights/20260317_full_cov_boost_ipw_weights.csv')
ipw_weights = weight_df %>% select(platform, weight_case, weight_noncase) %>% 
  tidyr::gather(case,wt,-platform) %>% 
  mutate(case = stringr::str_replace(pattern = "weight_", replacement = "", string = case)) 
# Add weights 
adata_long_ipw = adata_long %>% left_join(ipw_weights , by = c("platform","case"))
adata_wide_ipw = adata_wide  %>% left_join(ipw_weights , by = c("platform","case"))

# Next consider adaptive metrics
A = readxl::read_xlsx('/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/01_ID159_voysey/UniversityofOxford-Lambe-P01_02_FU2_downsampled_COVID_metrics_v2.xlsx')
names(A)
#Ax = A %>% select(any_of(c("sample_name","breadth_spike","breadth_MIRA","COVID_CD8_spike_breadth","COVID_CD4_spike_breadth","score_v3","score_spike")))
Ax = A %>% select(sample_name, 
                  adpt_breadth_spike     = breadth_spike,
                  adpt_breadth_MIRA      = breadth_MIRA,
                  adpt_COVID_CD8_spike_breadth = COVID_CD8_spike_breadth,
                  adpt_COVID_CD4_spike_breadth = COVID_CD4_spike_breadth, 
                  adpt_score_v3 = score_v3, 
                  adpt_score_spike = score_spike)
Ax
adata_long_ipw = adata_long_ipw %>% left_join(Ax, by = c("name_1"="sample_name"))
adpt_vars = c(
  "adpt_breadth_spike",
  "adpt_breadth_MIRA",
  "adpt_COVID_CD8_spike_breadth",
  "adpt_COVID_CD4_spike_breadth",
  "adpt_score_v3",
  "adpt_score_spike")

for (v in adpt_vars) {
  print(v)
  adata_long_ipw[[paste0("log2", v)]] = log2(adata_long_ipw[[v]])
}


adata_wide_ipw = adata_wide_ipw %>% left_join(Ax, by = c("name_1"= "sample_name"))%>% 
  left_join(Ax, by = c("name_1_pre"="sample_name"),suffix= c("","_pre"))

# Add log2 version]
adpt_vars = c(
  "adpt_breadth_spike",
  "adpt_breadth_MIRA",
  "adpt_COVID_CD8_spike_breadth",
  "adpt_COVID_CD4_spike_breadth",
  "adpt_score_v3",
  "adpt_score_spike",
  "adpt_breadth_spike_pre",
  "adpt_breadth_MIRA_pre",
  "adpt_COVID_CD8_spike_breadth_pre",
  "adpt_COVID_CD4_spike_breadth_pre",
  "adpt_score_v3_pre",
  "adpt_score_spike_pre"
)

for (v in adpt_vars) {
  adata_wide_ipw[[paste0("log2", v)]] = log2(adata_wide_ipw[[v]])
}

# PRELIMINARY FILE
if (write_files){
  adata_part_1 %>% 
    write.csv(file.path(project_folder , '001_COV_BOOST_ADATA.csv'), row.names= F)
  saveRDS(adata_part_1, file.path(project_folder , '001_COV_BOOST_ADATA.rds'))
}

# ADATA WIDE
if (write_files){
  adata_wide_ipw %>% 
    write.csv(file.path(project_folder , 'Manuscript', 'adata_v1', '001_COV_BOOST_ADATA_WIDE.csv'), row.names= F)
  saveRDS(adata_wide_ipw, 
          file.path(project_folder , 'Manuscript', 'adata_v1', '001_COV_BOOST_ADATA_WIDE.rds'))
}

# ADATA LONG
if (write_files){
  adata_long_ipw %>% 
    write.csv(file.path(project_folder , 'Manuscript', 'adata_v1', '001_COV_BOOST_ADATA_LONG.csv'), row.names= F)
  saveRDS(adata_long_ipw, 
          file.path(project_folder , 'Manuscript', 'adata_v1', '001_COV_BOOST_ADATA_LONG.rds'))
}

# ADDITIONAL VISUAL CHECKS OF IMPUTATION 
library(ggplot2)
plot_imputation_post_period <- function(){
  p1 = ggplot(imputed_df, aes(x = age, y = log2nt50, color = log2nt50_imp)) +
    geom_point(alpha = 0.7, aes(size = log2nt50_imp)) +
    scale_size_manual(values = c(1,3))+
    scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "red"),
      labels = c("Observed", "Imputed"),
      name = "Data Status"
    ) +
    labs(
      title = paste0("KNN Imputation Results (k=", k_impute,")"),
      subtitle = "Red points represent the values filled by the KNN algorithm",
      x = "Age",
      y = "log2nt50 Value"
    ) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90))
  p2 = ggplot(imputed_df, aes(y = platform, x = log2nt50, color = log2nt50_imp)) +
    geom_point(alpha = 0.7, aes(size = log2nt50_imp)) +
    scale_size_manual(values = c(1,3))+
    scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "red"),
      labels = c("Observed", "Imputed"),
      name = "Data Status"
    ) +
    labs(
      title = paste0("KNN Imputation Results (k=", k_impute,")"),
      subtitle = "Red points represent the values filled by the KNN algorithm",
      y = "",
      x = "log2elispot Value"
    ) +
    theme_minimal()+ theme(axis.text.x = element_text(angle = 90, size = 3))
  p3 = ggplot(imputed_df, aes(x = age, y = log2elispot, color = log2elispot_imp)) +
    geom_point(alpha = 0.7, aes(size = log2elispot_imp)) +
    scale_size_manual(values = c(1,3))+
    scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "blue"),
      labels = c("Observed", "Imputed"),
      name = "Data Status"
    ) +
    labs(
      title = paste0("KNN Imputation Results (k=", k_impute,")"),
      subtitle = "Blue points represent the values filled by the KNN algorithm",
      x = "Age",
      y = "log2elispot Value"
    ) +
    theme_minimal()
  p4 = ggplot(imputed_df, aes(y = platform, x = log2elispot, color = log2elispot_imp)) +
    geom_point(alpha = 0.7, aes(size = log2elispot_imp)) +
    scale_size_manual(values = c(1,3))+
    scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "blue"),
      labels = c("Observed", "Imputed"),
      name = "Data Status"
    ) +
    labs(
      title = paste0("KNN Imputation Results (k=", k_impute,")"),
      subtitle = "Blue points represent the values filled by the KNN algorithm",
      y = "",
      x = "log2elispot Value"
    ) +
    theme_minimal()+ theme(axis.text.x = element_text(angle = 90, size = 3))
  return(gridExtra::grid.arrange(p1,p2,p3,p4, ncol = 2))
}
plot_imputation_pre_period <- function(){
  p1x = ggplot(imputed_df_pre_nt50, aes(x = age, y = log2nt50, color = log2nt50_imp)) +
    geom_point(alpha = 0.7, aes(size = log2nt50_imp)) +
    scale_size_manual(values = c(1,3))+
    scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "red"),
      labels = c("Observed", "Imputed"),
      name = "Data Status"
    ) +
    labs(
      title = paste0("KNN Imputation Results (k=", k_impute,")"),
      subtitle = "Red points represent the values filled by the KNN algorithm",
      x = "Age",
      y = "log2nt50 Value"
    ) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90))
  p2x = ggplot(imputed_df_pre_nt50, aes(y = prime, x = log2nt50, color = log2nt50_imp)) +
    geom_point(alpha = 0.7, aes(size = log2nt50_imp)) +
    scale_size_manual(values = c(1,3))+
    scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "red"),
      labels = c("Observed", "Imputed"),
      name = "Data Status"
    ) +
    labs(
      title = paste0("KNN Imputation Results (k=", k_impute,")"),
      subtitle = "Red points represent the values filled by the KNN algorithm",
      y = "",
      x = "log2elispot Value"
    ) +
    theme_minimal()+ theme(axis.text.x = element_text(angle = 90, size = 3))
  gridExtra::grid.arrange(p1x,p2x, ncol = 2)
}
pdf('/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey/Manuscript/imputation_notes/post_data_imputation_check.pdf', 
    height = 8, width = 8)
plot_imputation_post_period()
dev.off()
pdf('/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey/Manuscript/imputation_notes/pre_data_imputation_check.pdf', 
    height = 4, width = 8)
plot_imputation_pre_period()
dev.off()
### --- ### 
### --- ### 
### --- ### 
### --- ### 
### Above this point we have combined all the metadata
### Below we load the official breadth estimates 
require(dplyr)
require(ggplot2)
r = "/Volumes/fh/fast/gilbert_p/kmayerbl/spike_atlas_data/estimate_20260527/"
fs = list.files(r)
scores_fp   = file.path(r, "cov_boost_day0_day14_genotyped_spike_breadth_estimate_x1w2.csv")
hits_fp     = file.path(r, "cov_boost_day0_day14_genotyped_spike_breadth_all_hits.csv") 
scores      = readr::read_csv(scores_fp) 
hits        = readr::read_csv(hits_fp)
adata       = readr::read_csv("001_COV_BOOST_ADATA.csv")
hits
names(hits)
locus_binary = readRDS("locus_binary.rds")
locus_binary <- locus_binary %>% mutate(tcellsub = if_else(stringr::str_starts(string = locus, pattern = '[ABC]'),"CD8","CD4"))

# Compute Ts
df_both = hits %>% 
  select(name_1, total_unique_1, binary, hla_detect, templates, productive_frequency, x1w2,yule_q_odds) %>% 
  left_join(locus_binary) %>% 
  filter(hla_detect==1) %>%
  mutate(depth_yule = productive_frequency * yule_q_odds) %>% 
  group_by(name_1,total_unique_1) %>% 
  summarise(n = n(), depth = sum(productive_frequency),depth2 = sum(depth_yule), s_x1w2= sum(x1w2)) %>%
  mutate(breadth = s_x1w2/total_unique_1) %>% ungroup()

df_cd4_cd8 = hits %>% 
  select(name_1, total_unique_1, binary, hla_detect, templates, productive_frequency, x1w2,yule_q_odds) %>% 
  left_join(locus_binary) %>% 
  filter(hla_detect==1) %>%
  mutate(depth_yule = productive_frequency * yule_q_odds) %>% 
  group_by(name_1,total_unique_1, tcellsub) %>% 
  summarise(n = n(), depth = sum(productive_frequency),depth2 = sum(depth_yule), s_x1w2= sum(x1w2)) %>%
  mutate(breadth = s_x1w2/total_unique_1) %>% ungroup()

df_cd4 = df_cd4_cd8 %>% filter(tcellsub == "CD4") 
df_cd8 = df_cd4_cd8 %>% filter(tcellsub == "CD8")

df_cd4_cd8_joined = left_join(df_cd4, df_cd8, suffix = c("_cd4", "_cd8"), by = c("name_1","total_unique_1")) 
df_all = df_both %>% left_join(df_cd4_cd8_joined, by = c("name_1","total_unique_1") )
# Run a check 
df_all = df_all %>% left_join(adata %>% select(name_1, timepoint, sc_repeat_pid, elispot, nt50))
df_all_post = df_all %>% filter(timepoint == "Day 14")
df_all_pre = df_all %>% filter(timepoint == "Day 0")
df_wide = df_all_post %>% left_join(df_all_pre, by = "sc_repeat_pid", suffix = c("_day14", "_day0"))

df_long = df_wide %>% tidyr::gather(key, value, -sc_repeat_pid) %>% 
  left_join(adata %>% filter(timepoint == "Day 14") %>% select( sc_repeat_pid , case, prime, boost, time, strata, male, age10,platform), by = "sc_repeat_pid")


variables = c('depth2_day14', 'depth2_cd4_day14', 'depth2_cd8_day14',
              'depth2_day0', 'depth2_cd4_day0', 'depth2_cd8_day0',
              'breadth_day0', 'breadth_cd4_day0', 'breadth_cd8_day0',
              'breadth_day14', 'breadth_cd4_day14', 'breadth_cd8_day14')

variables = c('depth2_day14', 'depth2_cd4_day14', 'depth2_cd8_day14',
              'breadth_day14', 'breadth_cd4_day14', 'breadth_cd8_day14',
              'elispot_day14', 'nt50_day14')



####

# Compute Ts Unmatched
df_both_UNMATCHED = hits %>% 
  select(name_1, total_unique_1, binary, hla_detect, templates, productive_frequency, x1w2,yule_q_odds) %>% 
  left_join(locus_binary) %>% 
  #filter(hla_detect==1) %>%
  mutate(depth_yule = productive_frequency * yule_q_odds) %>% 
  group_by(name_1,total_unique_1) %>% 
  summarise(n = n(), depth = sum(productive_frequency),depth2 = sum(depth_yule), s_x1w2= sum(x1w2)) %>%
  mutate(breadth = s_x1w2/total_unique_1) %>% ungroup()

df_cd4_cd8_UNMATCHED = hits %>% 
  select(name_1, total_unique_1, binary, hla_detect, templates, productive_frequency, x1w2,yule_q_odds) %>% 
  left_join(locus_binary) %>% 
  #filter(hla_detect==1) %>%
  mutate(depth_yule = productive_frequency * yule_q_odds) %>% 
  group_by(name_1,total_unique_1, tcellsub) %>% 
  summarise(n = n(), depth = sum(productive_frequency),depth2 = sum(depth_yule), s_x1w2= sum(x1w2)) %>%
  mutate(breadth = s_x1w2/total_unique_1) %>% ungroup()

df_cd4_UNMATCHED = df_cd4_cd8_UNMATCHED %>% filter(tcellsub == "CD4") 
df_cd8_UNMATCHED = df_cd4_cd8_UNMATCHED %>% filter(tcellsub == "CD8")

df_cd4_cd8_joined_UNMATCHED = left_join(df_cd4_UNMATCHED, df_cd8_UNMATCHED, suffix = c("_cd4", "_cd8"), by = c("name_1","total_unique_1")) 
df_all_UNMATCHED = df_both_UNMATCHED %>% left_join(df_cd4_cd8_joined_UNMATCHED, by = c("name_1","total_unique_1") )
# Run a check 
df_all_UNMATCHED = df_all_UNMATCHED %>% left_join(adata %>% select(name_1, timepoint, sc_repeat_pid, elispot, nt50))
df_all_post_UNMATCHED = df_all_UNMATCHED %>% filter(timepoint == "Day 14")
df_all_pre_UNMATCHED = df_all_UNMATCHED %>% filter(timepoint == "Day 0")
df_wide_UNMATCHED = df_all_post_UNMATCHED %>% left_join(df_all_pre_UNMATCHED, by = "sc_repeat_pid", suffix = c("_day14", "_day0"))

df_long_UNMATCHED = df_wide_UNMATCHED %>% tidyr::gather(key, value, -sc_repeat_pid) %>% 
  left_join(adata %>% filter(timepoint == "Day 14") %>% select( sc_repeat_pid , case, prime, boost, time, strata, male, age10,platform), by = "sc_repeat_pid")



# Compute Ts MISMATCHED
df_both_MISMATCHED = hits %>% 
  select(name_1, total_unique_1, binary, hla_detect, templates, productive_frequency, x1w2,yule_q_odds) %>% 
  left_join(locus_binary) %>% 
  filter(hla_detect==0) %>%
  mutate(depth_yule = productive_frequency * yule_q_odds) %>% 
  group_by(name_1,total_unique_1) %>% 
  summarise(n = n(), depth = sum(productive_frequency),depth2 = sum(depth_yule), s_x1w2= sum(x1w2)) %>%
  mutate(breadth = s_x1w2/total_unique_1) %>% ungroup()

df_cd4_cd8_MISMATCHED = hits %>% 
  select(name_1, total_unique_1, binary, hla_detect, templates, productive_frequency, x1w2,yule_q_odds) %>% 
  left_join(locus_binary) %>% 
  filter(hla_detect==0) %>%
  mutate(depth_yule = productive_frequency * yule_q_odds) %>% 
  group_by(name_1,total_unique_1, tcellsub) %>% 
  summarise(n = n(), depth = sum(productive_frequency),depth2 = sum(depth_yule), s_x1w2= sum(x1w2)) %>%
  mutate(breadth = s_x1w2/total_unique_1) %>% ungroup()

df_cd4_MISMATCHED = df_cd4_cd8_MISMATCHED %>% filter(tcellsub == "CD4") 
df_cd8_MISMATCHED = df_cd4_cd8_MISMATCHED %>% filter(tcellsub == "CD8")

df_cd4_cd8_joined_MISMATCHED = left_join(df_cd4_MISMATCHED, df_cd8_MISMATCHED, suffix = c("_cd4", "_cd8"), by = c("name_1","total_unique_1")) 
df_all_MISMATCHED = df_both_MISMATCHED %>% left_join(df_cd4_cd8_joined_MISMATCHED, by = c("name_1","total_unique_1") )
# Run a check 
df_all_MISMATCHED = df_all_MISMATCHED %>% left_join(adata %>% select(name_1, timepoint, sc_repeat_pid, elispot, nt50))
df_all_post_MISMATCHED = df_all_MISMATCHED %>% filter(timepoint == "Day 14")
df_all_pre_MISMATCHED = df_all_MISMATCHED %>% filter(timepoint == "Day 0")
df_wide_MISMATCHED = df_all_post_MISMATCHED %>% left_join(df_all_pre_MISMATCHED, by = "sc_repeat_pid", suffix = c("_day14", "_day0"))

df_long_MISMATCHED = df_wide_MISMATCHED %>% tidyr::gather(key, value, -sc_repeat_pid) %>% 
  left_join(adata %>% filter(timepoint == "Day 14") %>% select( sc_repeat_pid , case, prime, boost, time, strata, male, age10,platform), by = "sc_repeat_pid")




# BREADTH ADATA LONG
if (write_files){
  df_long %>% 
    write.csv(file.path(project_folder , 'Manuscript', 'adata_v1', '001_COV_BOOST_ADATA_LONG_BREADTH.csv'), row.names= F)
  saveRDS(df_long, 
          file.path(project_folder , 'Manuscript', 'adata_v1', '001_COV_BOOST_ADATA_LONG_BREADTH.rds'))
}


# BREADTH WIDE
if (write_files){
  df_wide %>% 
    write.csv(file.path(project_folder , 'Manuscript', 'adata_v1', '001_COV_BOOST_ADATA_WIDE_BREADTH.csv'), row.names= F)
  saveRDS( df_wide, 
          file.path(project_folder , 'Manuscript', 'adata_v1', '001_COV_BOOST_ADATA_WIDE_BREADTH.rds'))
}

