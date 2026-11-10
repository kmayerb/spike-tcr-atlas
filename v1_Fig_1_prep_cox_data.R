#v1_Fig_1_prep_cox_data.R

# New inputs (place in data/):
#   001_COV_BOOST_ADATA_WIDE.csv          -> df1  (main wide adata)
#   001_COV_BOOST_ADATA_WIDE_BREADTH.csv  -> df2  (T-cell breadth)
#
# What this does:
#   * joins df2's T-cell breadth onto df1 by `sc_repeat_pid`
#   * derives the T_S variables the old pipeline called log2tcell/_pre
#       log2tcell     = log2(breadth_day14)
#       log2tcell_pre = log2(breadth_day0)
#   * derives case_abr, coerces male -> logical, maps status <- event
#   * writes out 
#   * prints a validation report (required columns, NA counts, ranges,
#     and an event-vs-Status cross-tab so you can confirm status==event)
#
# =====================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ---- paths ----------------------------------------------------------
project_folder = '/Volumes/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/04_ID159_Voysey/'
f_df1 = file.path(project_folder, 'Manuscript','adata' ,'001_COV_BOOST_ADATA_WIDE.csv')
f_df2 = file.path(project_folder, 'Manuscript','adata' ,'001_COV_BOOST_ADATA_WIDE_BREADTH.csv')
stopifnot(file.exists(f_df1), file.exists(f_df2))

outfile <-file.path(project_folder, 'Manuscript','adata' , "id159_prop_hazard_adata.csv")

# ---- read -----------------------------------------------------------
df1 <- readr::read_csv(f_df1, show_col_types = FALSE)
df2 <- readr::read_csv(f_df2, show_col_types = FALSE)

# ---- T_S breadth from df2 (confirmed mapping) -----------------------
# Overall spike-reactive T-cell breadth: day14 = post, day0 = pre.
ts_cols <- c("sc_repeat_pid", "breadth_day14", "breadth_day0")
missing_ts <- setdiff(ts_cols, names(df2))
if (length(missing_ts))
  stop("df2 is missing expected column(s): ", paste(missing_ts, collapse = ", "))

df2_ts <- df2 %>%
  select(all_of(ts_cols)) %>%
  distinct(sc_repeat_pid, .keep_all = TRUE)

# ---- join + derive canonical columns --------------------------------
adata <- df1 %>%
  left_join(df2_ts, by = "sc_repeat_pid") %>%
  mutate(
    # --- T_S (T-cell breadth) on log2 scale, as the old `log2tcell` ---
    tcell         = breadth_day14,
    tcell_pre     = breadth_day0,
    log2tcell     = log2(breadth_day14),
    log2tcell_pre = log2(breadth_day0),
    
    # --- x-axis label for the C-vs-NC boxplots -----------------------
    case_abr = if_else(case == "case", "C", "NC"),
    
    # --- male must be logical so the Cox term reads `maleTRUE` --------
    male = as.logical(male),
    
    # --- survival status: old pipeline used Surv(time, status) --------
    #     new data ships `event`; assume identical (VERIFY via report).
    status = event
  )

# ---- write ----------------------------------------------------------
readr::write_csv(adata, outfile)
message("Wrote ", outfile, "  (", nrow(adata), " rows, ", ncol(adata), " cols)")

# ---- validation report ---------------------------------------------
cat("\n================ VALIDATION ================\n")

required <- c("case", "case_abr", "time", "event", "status", "strata", "wt",
              "age10", "male",
              "log2nt50", "log2nt50_pre",
              "log2elispot", "log2elispot_pre",
              "log2tcell", "log2tcell_pre")
present <- required %in% names(adata)
cat("\nRequired columns:\n")
print(data.frame(column = required, present = present))
if (any(!present))
  cat("\n** MISSING:", paste(required[!present], collapse = ", "), "**\n")

cat("\nNA counts (required numeric markers):\n")
num_req <- intersect(c("log2nt50","log2nt50_pre","log2elispot","log2elispot_pre",
                       "log2tcell","log2tcell_pre","time","event","wt","age10"),
                     names(adata))
print(sapply(adata[num_req], function(x) sum(is.na(x))))

cat("\nlog2tcell range (expect roughly -13 .. -9 to match prior figures):\n")
print(round(range(adata$log2tcell, finite = TRUE), 2))
cat("log2tcell_pre range:\n")
print(round(range(adata$log2tcell_pre, finite = TRUE), 2))

# Does breadth contain zeros (-> -Inf after log2)? replace_neginf handles it downstream.
n_neginf <- sum(is.infinite(adata$log2tcell) | is.infinite(adata$log2tcell_pre))
cat("\n-Inf values in log2tcell/_pre (from breadth == 0):", n_neginf, "\n")

# ---- STATUS vs EVENT check -----------------------------------------
cat("\n---- Confirm status == event ----\n")
cat("Cross-tab event vs numeric status:\n")
print(table(event = adata$event, status = adata$status, useNA = "ifany"))
if ("Status" %in% names(df1)) {
  cat("\nCross-tab numeric event vs text Status (from df1):\n")
  print(table(event = adata$event, Status = df1$Status, useNA = "ifany"))
  cat("\n^ Check that this text/numeric correspondence is what you expect.\n")
}
cat("============================================\n")
message("Wrote ", outfile, "  (", nrow(adata), " rows, ", ncol(adata), " cols)")
