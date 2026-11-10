# 003_COMPARATIVE_BREADTH.R

# We need a clean figure that shows breadth estimates across key studies 
# Emerson CMV+ CMV- # 
# May 26 You did this previously
require(ggplot2)
require(dplyr)
require(readr)
library(stringr)

r = "/Volumes/fh/fast/gilbert_p/kmayerbl/spike_atlas_data/estimate_20260527/"
fs = list.files(r)
natgen_information = read_tsv( '/Volumes/fh/fast/gilbert_p/kmayerbl/large_cohorts/emerson_official_natgen_v2/emerson_natgen_SampleOverview_06-11-2026_7-56-05_PM.tsv')
natgen = natgen_information %>% 
  select(name_1 = sample_name, sample_tags) %>% 
  mutate(cmv = case_when(
    str_detect(sample_tags, fixed("Cytomegalovirus +")) ~ "CMV+",
    str_detect(sample_tags, fixed("Cytomegalovirus -")) ~ "CMV-",
    TRUE ~ NA_character_
  )) %>% 
  mutate(sex = case_when(
    str_detect(sample_tags, fixed("Male")) ~ "Male",
    str_detect(sample_tags, fixed("Female")) ~ "Female",
    TRUE ~ NA_character_
  )) %>% 
  mutate(age = as.numeric(str_extract(sample_tags, "\\d+(?=\\s*[Yy]ears?)"))) %>% 
  select(name_1, sex, age, cmv) %>% 
  mutate(exposure_level = "N") %>% 
  mutate(exposure_level_detail = paste0(exposure_level,"_natgen_",cmv)) %>% 
  mutate(group = "Naive")
natgen %>% write.csv('/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_CMV_classifier/data/natgen_serostatus.csv', row.names = F)


# THIS CONTAINS ALL OF THE EXPOSURE DATA WE NEED
metadata = readr::read_csv('/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_ID159_tools/data/20260520_breadth_across_13_studies_with_imputed_hla_alleles_draft_metadata.csv')
# WORKFLOW NOTE: 
# CONSTRUCT <m_*> metadata from each study, 
# COMBIEN into <meta> for plotting purposes

m_natgen = natgen_cmv = natgen %>%
  filter(!is.na(cmv))


m_v4b = readr::read_tsv('/Volumes/fh/fast/gilbert_p/kmayerbl/large_cohorts/2026_L1/hsTCRb_v4b/hsTCRB-V4b_SampleOverview_05-12-2026_6-23-16_PM.tsv') %>% 
  select(name_1= sample_name, sample_tags) %>% 
  mutate(stype = case_when(
    str_detect(sample_tags, fixed("Blood")) ~ "Blood",
    str_detect(sample_tags, fixed("PBMC")) ~ "PBMC",
    #str_detect(sample_tags, fixed("FFPE")) ~ "FFPE",
    TRUE ~ NA_character_
  )) %>% 
  mutate(cmv = case_when(
    str_detect(sample_tags, fixed("Cytomegalovirus +")) ~ "CMV+",
    str_detect(sample_tags, fixed("Cytomegalovirus -")) ~ "CMV-",
    TRUE ~ NA_character_
  )) %>% 
  mutate(sex = case_when(
    str_detect(sample_tags, fixed("Male")) ~ "Male",
    str_detect(sample_tags, fixed("Female")) ~ "Female",
    TRUE ~ NA_character_
  )) %>% 
  mutate(age = as.numeric(str_extract(sample_tags, "\\d+(?=\\s*[Yy]ears?)")))  %>%
  mutate(exposure_level = "N") %>%
  mutate(exposure_level_detail = paste0(exposure_level,"_","v4b", "_", stype)) %>% 
  filter(!is.na(stype)) %>%
  mutate(group = "Naive")


adata = readr::read_csv('001_COV_BOOST_ADATA.csv')

m_covboost = adata %>% 
  select(name_1, sex, age, timepoint) %>% 
  mutate(exposure_level = timepoint) %>% 
  mutate(exposure_level_detail = paste0("cov_boost_", exposure_level)) %>%
  mutate(group = "Vacciation")

m_delmonte_i = metadata %>% filter(archive_name_x == "delmonte") %>%
  select(name_1, sex, age, severity, icu, exposure_level) %>% 
  filter(exposure_level == "I")  %>% 
  mutate(exposure_level_detail = paste0("I_delmonte","_", severity)) %>% 
  mutate(group = if_else(exposure_level == "I", 'Infection COVID-19',NA))

m_cov3006 =  metadata %>% filter(archive_name_x == "cov3006") %>%
  select(name_1, exposure_level) %>% 
  filter(exposure_level %in% c('I', 'V', 'H')) %>% 
  mutate(exposure_level = if_else(exposure_level == "H", "V", exposure_level)) %>%
  mutate(exposure_level_detail = paste0(exposure_level, "_cov3006")) %>% 
  mutate(group = if_else(exposure_level == "I", 'Infection COVID-19', "Vaccination"))


m_elyanow =  metadata %>% filter(archive_name_x == "elyanow") %>%
  filter(exposure_level %in% c("I","I45"))%>%
  select(name_1,exposure_level = exposure_level2, exposure_level1 = exposure_level) %>% 
  filter(!is.na(exposure_level)) %>% 
  mutate(exposure_level_detail = paste0("C_", exposure_level, "_elyanow")) %>% 
  mutate(group = 'Recovered COVID-19')# %>% 
  #mutate(id = str_remove(name_1, "_.*"))


m_alter = metadata %>% filter(archive_name_x == "ad26_alter") %>% 
  filter(exposure_level == "V1") %>%
  select(name_1, exposure_level) %>% 
  mutate(exposure_level_detail = "V_Ad26_dose_1")

m_cov2= metadata %>% filter(archive_name_x == "cov2") %>% 
  filter(exposure_level %in% "V") %>%
  select(name_1, exposure_level) %>% 
  mutate(exposure_level_detail = "V_ChAdOx")


meta = bind_rows(m_natgen, 
                 #m_v4b, 
                 m_covboost,
                 m_delmonte_i,
                 m_cov3006,
                 m_elyanow, 
                 m_cov2,
                 m_alter)

meta$exposure_level_detail %>%table()
# WORKFLOW NOTE: 
  # LOAD breadth from various runs, named as b1,b2 ... bn
  # COMBINE into <br>
  # CONSTRUCT data from LEFT_JOIN <meta> to <br>, 
  # such that only data with available metadata is used.
b1 = read_csv(file.path(r,'emerson_natgen_v2_spike_breadth_estimate_x1w2.csv'))
b2 = read_csv(file.path(r,'v4b_ctrl_spike_breadth_estimate_x1w2.csv'))
b3 = read_csv(file.path(r,'cov_boost_day0_day14_genotyped_spike_breadth_estimate_x1w2.csv'))
b4 = read_csv(file.path(r,'cov3006_spike_breadth_estimate_x1w2.csv'))
b4a= read_csv(file.path(r,'cov2_spike_breadth_estimate_x1w2.csv'))
b4b= read_csv(file.path(r,'ad26_alter_spike_breadth_estimate_x1w2.csv'))
b5 = read_csv(file.path(r,'delmonte_spike_breadth_estimate_x1w2.csv'))
b6= read_csv(file.path(r,'elyanow_spike_breadth_estimate_x1w2.csv'))

br = bind_rows(b1,b3,b4,b4a,b4b, b5, b6)

chars = unique(meta$exposure_level_detail)
two_part_combinations <- combn(chars, 2, simplify = FALSE)
dput(two_part_combinations)

data = meta %>% left_join(br, by = "name_1")
table(data$exposure_level_detail)

thx = readRDS("theme_x.rds")
data$exposure_level_detail %>% unique() %>% dput()


e_levels = c(#"N_v4b_PBMC", "N_v4b_Blood",
             "N_natgen_CMV+", "N_natgen_CMV-",  
  "cov_boost_Day 0", "cov_boost_Day 14",  "V_cov3006", "V_ChAdOx", "V_Ad26_dose_1", "I_cov3006", "I_delmonte_Mild", "I_Mild","I_delmonte_Moderate_severe", 
 "I_delmonte_Critical" , "C_I_elyanow","C_I+fever_elyanow")

theme_x = theme_bw(base_size = 10, base_family = "Helvetica") + 
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
    
    # Nature Typography Specs (7pt labels, 8pt titles)
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 8, colour = "black"),
    
    # Clean facet strip appearance
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold", colour = "black")
  )


data$exposure_level_detail %>% table()
data %>% filter(exposure_level_detail == 'V_Ad26_dose_1')
data%>% mutate(exposure_level_detail  = if_else(exposure_level_detail %in% c( "I_cov3006", "I_delmonte_Mild") , "I_Mild", exposure_level_detail)) %>% 
  ggplot(aes(x = factor(exposure_level_detail, levels = e_levels), y = breadth)) + 
  geom_boxplot(aes(fill = factor(exposure_level_detail, levels = e_levels)), col = "black", outlier.shape = NA, width = .75) + 
  # Replaced geom_boxplot with a blue median line
  ggbeeswarm::geom_quasirandom(width = .35, size = .2) + 
  #stat_summary(fun = median, fun.min = median, fun.max = median, 
  #             geom = "crossbar", width = 0.5, color = "white") +
  coord_cartesian(ylim = c(1E-5, 8E-3)) + 
  theme_x + 
  scale_y_log10()+
  annotation_logticks(side = "l", linewidth = .5)+
  ggpubr::stat_compare_means(
    size = 4, 
    fontface = "bold", # This will work once the lines actually render!
    comparisons = 
      list(
        c("N_natgen_CMV+", "N_natgen_CMV-"), 
        c("N_v4b_PBMC", "N_v4b_Blood"),
        c("N_natgen_CMV-", "cov_boost_Day 0"),
        c("cov_boost_Day 0", "cov_boost_Day 14"),
        c("cov_boost_Day 14","V_cov3006"),
        c("cov_boost_Day 14","V_ChAdOx"),
        c("cov_boost_Day 14","V_Ad26_dose_1"),
        c("cov_boost_Day 14","I_Mild"),
        c("cov_boost_Day 14","I_delmonte_Moderate_severe"),
        c("cov_boost_Day 14","I_delmonte_Critical"),
        c("I_Mild","I_delmonte_Moderate_severe"),
        c("I_delmonte_Moderate_severe","I_delmonte_Critical"),
        c("C_I_elyanow","C_I+fever_elyanow")),
    #step.increase = 0.04,  
    tip.length = 0.00,vjust = .3,    
    
    label.y = c(-3.5, -3.5, -3.3, -2.6, -5.3,-4.8 ,-4.9, -5.0, -5.1, -5.2, -2.6, -2.4, -2.9),
    label = "p.signif") + scale_fill_manual(values = c("gray","gray","gray","gray","steelblue", "#0033CC", "#0066FF","#0099FF", "skyblue","salmon", "red","brown","darkseagreen", "forestgreen")) + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1))


figure_comparative_breadth = data%>% mutate(exposure_level_detail  = if_else(exposure_level_detail %in% c( "I_cov3006", "I_delmonte_Mild") , "I_Mild", exposure_level_detail)) %>% 
  ggplot(aes(x = factor(exposure_level_detail, levels = e_levels), y = breadth)) + 
  geom_boxplot(aes(fill = factor(exposure_level_detail, levels = e_levels)), col = "black", outlier.shape = NA, width = .75) + 
  # Replaced geom_boxplot with a blue median line
  ggbeeswarm::geom_quasirandom(width = .35, size = .2) + 
  #stat_summary(fun = median, fun.min = median, fun.max = median, 
  #             geom = "crossbar", width = 0.5, color = "white") +
  coord_cartesian(ylim = c(1E-5, 8E-3)) + 
  theme_x + 
  scale_y_log10()+
  annotation_logticks(side = "l", linewidth = .5)+
  ggpubr::stat_compare_means(
    size = 3, 
    fontface = "bold", # This will work once the lines actually render!
    comparisons = 
      list(
        c("N_natgen_CMV+", "N_natgen_CMV-"), 
        #c("N_v4b_PBMC", "N_v4b_Blood"),
        c("N_natgen_CMV-", "cov_boost_Day 0"),
        c("cov_boost_Day 0", "cov_boost_Day 14"),
        c("cov_boost_Day 14","V_cov3006"),
        c("cov_boost_Day 14","V_ChAdOx"),
        c("cov_boost_Day 14","V_Ad26_dose_1"),
        c("cov_boost_Day 14","I_Mild"),
        c("cov_boost_Day 14","I_delmonte_Moderate_severe"),
        c("cov_boost_Day 14","I_delmonte_Critical"),
        c("I_Mild","I_delmonte_Moderate_severe"),
        c("I_delmonte_Moderate_severe","I_delmonte_Critical"),
        c("C_I_elyanow","C_I+fever_elyanow")),
    #step.increase = 0.04,  
    tip.length = 0.00,vjust = .3,    
    
    label.y = c(-3.5, -3.3, -2.6, -5.3,-4.8 ,-4.9, -5.0, -5.1, -5.2, -2.6, -2.4, -2.9),
    ) + scale_fill_manual(values = c("gray","gray","steelblue", "#0033CC", "#0066FF","#0099FF", "skyblue","salmon", "red","brown","darkseagreen", "forestgreen")) + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1))


pdf("/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey/Manuscript/figures_v1/Fig_3_comparative_breadth.pdf", 
    width = 12, height = 6)
figure_comparative_breadth 
dev.off()
