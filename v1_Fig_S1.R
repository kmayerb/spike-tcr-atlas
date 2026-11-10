
theme_clean <- function(base_size = 10, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(color = "black")
    )
}

hla_colors = c(
  'A*01:01'='#e6194b',
  'A*02:01'='#3cb44b',
  'A*03:01'='#ffe119',
  'A*11:01'='#0082c8',
  'A*23:01'='#f58231',
  'A*24:02'='#911eb4',
  'A*26:01'='#46f0f0',
  'A*29:02'='#f032e6',
  'A*68:01'='#d2f53c',
  'B*07:02'='#fabebe',
  'B*08:01'='#008080',
  'B*14:02'='#e6beff',
  'B*15:01'='#aa6e28',
  'B*18:01'='#fffac8',
  'B*27:05'='#800000',
  'B*35:01'='#aaffc3',
  'B*40:01'='#808000',
  'B*44:02'='#ffd8b1',
  'B*44:03'='#000080',
  'B*51:01'='black',
  'B*57:01'='#ffffff',
  'C*02:02'='#808080',
  'C*03:04'='#a9a9a9',
  'C*04:01'='#1aff1a',
  'C*05:01'='#ff1a1a',
  'C*06:02'='#1a1aff',
  'C*07:01'='#ff8c00',
  'C*07:02'='#9932cc',
  'C*08:02'='#ff1493',
  'C*12:03'='#00ced1',
  'C*16:01'='#ff6347',
  'A*32:01'='#2e8b57')
hla_colors_dr = c(
  'DRB1*01:01'='#1f77b4',  # Dark Blue
  'DRB1*03:01'='#aec7e8',  # Light Blue
  'DRB1*04:01'='#ff7f0e',  # Dark Orange
  'DRB1*04:02'='#ffbb78',  # Light Orange
  'DRB1*04:04'='#2ca02c',  # Dark Green
  'DRB1*07:01'='#98df8a',  # Light Green
  'DRB1*08:01'='#d62728',  # Red
  'DRB1*11:01'='#ff9896',  # Pink / Light Red
  'DRB1*13:01'='#9467bd',  # Dark Purple
  'DRB1*13:02'='#c5b0d5',  # Light Purple
  'DRB1*15:01'='#8c564b')  # Brown


require(ggplot2)
require(dplyr)

data0 = tibble::tribble(
  ~binary, ~index,
  "DRB1*15:01", 550,
  "DRB1*13:02", 281,
  "DRB1*13:01", 220,
  "DRB1*07:01", 151,
  "DRB1*03:01", 146,
  "DRB1*01:01", 118,
  "DRB1*04:01", 75,
  "DRB1*04:04", 68,
  "DRB1*11:01", 47,
  "DRB1*08:01", 28,
  "DRB1*04:02", 3)

data = tibble::tribble(
  ~binary, ~index,
  "DRB1*15:01", 242,
  "DRB1*07:01", 99,
  "DRB1*03:01", 91,
  "DRB1*13:02", 69,
  "DRB1*13:01", 66,
  "DRB1*01:01", 65,
  "DRB1*04:01", 45,
  "DRB1*04:04", 37,
  "DRB1*08:01", 18,
  "DRB1*11:01", 18
)

data_public = tibble::tribble(
  ~binary, ~index,
  "DRB1*15:01", 105,
  "DRB1*07:01", 50,
  "DRB1*03:01", 45,
  "DRB1*01:01", 32,
  "DRB1*13:01", 23,
  "DRB1*04:01", 20,
  "DRB1*13:02", 17,
  "DRB1*04:04", 12,
  "DRB1*08:01", 6,
  "DRB1*11:01", 6)



plt_dr = ggplot(data=data0, aes(y = forcats::fct_reorder(binary, index), x = index)) + 
  geom_bar(data = data0,  stat = 'identity',
           aes(fill = forcats::fct_relevel(binary,data0$binary)), alpha = .5)+
  geom_bar(data = data0,  
           fill = 'lightgray',stat = 'identity')+
  geom_bar(data = data, stat = 'identity', 
           aes(fill = forcats::fct_relevel(binary,data0$binary)), alpha = .5)+
  geom_bar(data = data_public, stat = 'identity', 
           aes(fill =binary))+
  scale_fill_manual("",values = hla_colors_dr) + 
  theme_clean() + 
  ylab("") + xlab("") + 
  geom_text(aes(label = index), data = data, nudge_x = -.5, hjust = 1, size = 2) + 
  geom_text(aes(label = index, x = 0), data = data_public, nudge_x = 1, hjust = 0, col = "white", size = 2) + 
  theme(axis.text.x = element_text(angle = 0)) + 
  scale_x_sqrt(expand = c(0,0), breaks = c(5, 20, 50, 100, 200, 300, 400, 500)) 


# ---- Save ---------------------------------------------------------------
fig_width = 4.5
fig_height = 3
out_dir    <- "Manuscript/figures_v1"   # <-- set destination folder
out_name   <- "v1_Fig_S1_bar_plot_hla_dr"            # filename stem (no extension)
fig_dpi = 600
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".png")),
  plot = plt_dr, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi
)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".pdf")),
  plot = plt_dr, width = fig_width, height = fig_height, units = "in"
)


#######################
### --- class I --- ###
#######################
require(tibble)
data0 = tribble(
  ~binary, ~index,
  "A*24:02", 509,
  "A*01:01", 404,
  "A*03:01", 404,
  "A*02:01", 212,
  "B*51:01", 206,
  "C*07:02", 137,
  "C*08:02", 126,
  "B*07:02", 86,
  "B*14:02", 79,
  "C*04:01", 78,
  "B*08:01", 68,
  "A*26:01", 52,
  "A*68:01", 42,
  "C*06:02", 38,
  "C*07:01", 37,
  "B*40:01", 35,
  "B*44:02", 35,
  "C*03:04", 33,
  "C*12:03", 32,
  "B*35:01", 30,
  "B*44:03", 29,
  "C*05:01", 29,
  "A*29:02", 26,
  "B*15:01", 24,
  "C*16:01", 21,
  "A*11:01", 18,
  "B*27:05", 17,
  "B*18:01", 12,
  "C*02:02", 7,
  "A*32:01", 6,
  "B*57:01", 4,
  "A*23:01", 3
)
data = tribble(
  ~binary, ~index,
  "A*24:02", 199,
  "A*01:01", 128,
  "A*02:01", 114,
  "A*03:01", 107,
  "C*07:02", 64,
  "B*51:01", 59,
  "B*07:02", 55,
  "C*08:02", 50,
  "B*08:01", 47,
  "B*14:02", 30,
  "C*04:01", 25,
  "C*03:04", 23,
  "C*07:01", 22,
  "B*40:01", 22,
  "B*44:02", 21,
  "C*06:02", 20,
  "B*44:03", 18,
  "C*05:01", 17,
  "B*35:01", 16,
  "A*26:01", 15,
  "B*15:01", 11,
  "A*29:02", 11,
  "C*12:03", 10,
  "C*16:01", 10,
  "B*27:05", 9,
  "A*11:01", 7,
  "B*18:01", 6,
  "C*02:02", 5,
  "B*57:01", 3,
  "A*68:01", 3,
  "A*23:01", 2,
  "A*32:01", 1
)
data_public = tribble(
  ~binary, ~index,
  "A*24:02", 70,
  "A*02:01", 60,
  "A*03:01", 36,
  "C*07:02", 30,
  "A*01:01", 27,
  "B*07:02", 26,
  "B*08:01", 26,
  "B*51:01", 21,
  "C*07:01", 14,
  "C*06:02", 13,
  "C*03:04", 12,
  "B*44:02", 11,
  "B*44:03", 10,
  "B*40:01", 9,
  "C*05:01", 8,
  "C*16:01", 8,
  "C*04:01", 7,
  "B*14:02", 6,
  "B*15:01", 5,
  "C*08:02", 5,
  "C*12:03", 5,
  "A*29:02", 4,
  "B*35:01", 4,
  "B*27:05", 2,
  "A*26:01", 2,
  "B*18:01", 1,
  "B*57:01", 1,
  "C*02:02", 1,
  "A*32:01", 1,
  "A*23:01", 1,
  "A*11:01", 1
)

hla_colors = c(
  'A*01:01'='#e6194b',
  'A*02:01'='#3cb44b',
  'A*03:01'='#ffe119',
  'A*11:01'='#0082c8',
  'A*23:01'='#f58231',
  'A*24:02'='#911eb4',
  'A*26:01'='#46f0f0',
  'A*29:02'='#f032e6',
  'A*68:01'='#d2f53c',
  'B*07:02'='#fabebe',
  'B*08:01'='#008080',
  'B*14:02'='#e6beff',
  'B*15:01'='#aa6e28',
  'B*18:01'='#fffac8',
  'B*27:05'='#800000',
  'B*35:01'='#aaffc3',
  'B*40:01'='#808000',
  'B*44:02'='#ffd8b1',
  'B*44:03'='#000080',
  'B*51:01'='black',
  'B*57:01'='#ffffff',
  'C*02:02'='#808080',
  'C*03:04'='#a9a9a9',
  'C*04:01'='#1aff1a',
  'C*05:01'='#ff1a1a',
  'C*06:02'='#1a1aff',
  'C*07:01'='#ff8c00',
  'C*07:02'='#9932cc',
  'C*08:02'='#ff1493',
  'C*12:03'='#00ced1',
  'C*16:01'='#ff6347',
  'A*32:01'='#2e8b57')
require(ggplot2)
data0
plt_ABC = ggplot(data=data0, aes(y = forcats::fct_reorder(binary, index), x = index)) + 
  geom_bar(data = data0,  stat = 'identity',
           aes(fill = forcats::fct_relevel(binary,data0$binary)), alpha = .5)+
  geom_bar(data = data0,  
           fill = 'lightgray',stat = 'identity')+
  geom_bar(data = data, stat = 'identity', 
           aes(fill = forcats::fct_relevel(binary,data0$binary)), alpha = .5)+
  geom_bar(data = data_public, stat = 'identity', 
           aes(fill =binary))+
  scale_fill_manual("",values = hla_colors) + 
  theme_clean() + 
  ylab("") + xlab("") + 
  geom_text(aes(label = index), data = data, nudge_x = -.5, hjust = 1, size = 2) + 
  geom_text(aes(label = index, x = 0), data = data_public, nudge_x = 1, hjust = 0, col = "white", size = 2) + 
  theme(axis.text.x = element_text(angle = 0)) + 
  scale_x_sqrt(expand = c(0,0), breaks = c(5, 20, 50, 100, 200, 300, 400, 500)) 


# ---- Save ---------------------------------------------------------------
fig_width = 5.6
fig_height = 6
out_dir    <- "Manuscript/figures_v1"   # <-- set destination folder
out_name   <- "v1_Fig_S1_bar_plot_hla_abc"            # filename stem (no extension)
fig_dpi = 300
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".png")),
  plot = plt_ABC, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi
)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(out_name, ".pdf")),
  plot = plt_ABC, width = fig_width, height = fig_height, units = "in"
)

#save_pdf(plt, fig_folder = archive$archive_folder, fig_name = "HLA_A_B_C_catalog_count", width = 5.5, height = 6)

