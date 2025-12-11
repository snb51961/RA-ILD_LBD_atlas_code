############################################################
# 03) White Map pre-analysis core  — ORIGINAL
# ----------------------------------------------------------
# This script corresponds to the original "White map 全解析コード.txt"
# It computes:
#   - AC co-occurrence (AC_known)
#   - AB/BC evidence summaries (best/max/sum)
#   - priority scoring
#   - white_map_candidates_*.csv
#   - top_white_candidates_*.csv
#
# These outputs are REQUIRED by:
#   - 07_whiteMap_Bbreakdown_fallback2.R  (Figure 5B–C)
#
# Path is generalized:
ROOT <- "/path/to/RA-ILD_textanalysis"
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
  library(tidyr); library(purrr); library(ggplot2)
})

DIR_OUT <- file.path(ROOT, "output")
DIR_FIG <- file.path(ROOT, "fig")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

STAMP <- format(Sys.time(), "%Y%m%d_%H%M")

############################################################
# 1) Load hits_matrix and ABC rankings
############################################################

# latest hits_matrix
hits_file <- list.files(
  file.path(ROOT, "data_proc"),
  "^hits_matrix_.*\\.csv$", full.names = TRUE
) |> sort() |> tail(1)
stopifnot(length(hits_file)==1)
hits <- readr::read_csv(hits_file, show_col_types = FALSE)

# latest ABC
abc_file <- list.files(
  file.path(ROOT,"output"),
  "^abc_rankings_.*\\.csv$", full.names = TRUE
) |> sort() |> tail(1)
stopifnot(length(abc_file)==1)
abc <- readr::read_csv(abc_file, show_col_types = FALSE)

############################################################
# 2) Identify A / C terms
############################################################
hit_cols  <- names(hits)[startsWith(names(hits),"hit__")]
terms_all <- sub("^hit__", "", hit_cols)

C_terms <- unique(abc$C)
A_terms <- intersect(unique(abc$A), terms_all)   # keep only A appearing in hits_matrix

############################################################
# 3) AC_known calculation
############################################################

binv <- function(v) as.integer(ifelse(is.na(v),0L, v>0L))
get_hit <- function(df, t){
  col <- paste0("hit__",t)
  if (!(col %in% names(df))) return(integer(nrow(df)))
  binv(df[[col]])
}

ac_known <- lapply(A_terms, function(a){
  tibble(
    A = a,
    !!!setNames(lapply(C_terms, function(c){
      sum(get_hit(hits,a) & get_hit(hits,c))
    }), paste0("AC_n11__", C_terms))
  )
}) |> bind_rows()

ac_long <- tidyr::pivot_longer(
  ac_known,
  starts_with("AC_n11__"),
  names_to="C",
  values_to="AC_n11"
) %>%
  mutate(
    C = sub("^AC_n11__","",C),
    AC_known = AC_n11 > 0
  )

############################################################
# 4) AB/BC evidence summary by B (max/sum)
############################################################

abc_summ <- abc %>%
  group_by(A,C) %>%
  summarise(
    best_score_q = max(score_q, na.rm=TRUE),
    AB_evi_sum   = sum(AB_n11, na.rm=TRUE),
    BC_evi_sum   = sum(BC_n11, na.rm=TRUE),
    AB_evi_max   = max(AB_n11, na.rm=TRUE),
    BC_evi_max   = max(BC_n11, na.rm=TRUE),
    n_B          = n(),
    .groups="drop"
  )

############################################################
# 5) Combine → WhiteMap table
############################################################

white_map <- ac_long %>%
  left_join(abc_summ, by=c("A","C")) %>%
  mutate(
    ABBC_strength = pmax(AB_evi_max,0) + pmax(BC_evi_max,0),
    priority      = dplyr::coalesce(best_score_q, 0) + 0.1*ABBC_strength
  ) %>%
  arrange(desc(!AC_known), desc(priority))

out_white <- file.path(DIR_OUT, paste0("white_map_candidates_",STAMP,".csv"))
readr::write_csv(white_map, out_white)
message("Saved: ", out_white)

############################################################
# 6) Top 10 extraction (priority-based)
############################################################

cand_unknown <- white_map %>% filter(!AC_known)
if (nrow(cand_unknown)==0) cand_unknown <- white_map

cand_ranked <- cand_unknown %>%
  mutate(ABBC_sum = AB_evi_sum + BC_evi_sum) %>%
  arrange(
    desc(priority),
    desc(best_score_q),
    desc(ABBC_strength),
    desc(ABBC_sum),
    desc(n_B)
  )

topN <- min(10, nrow(cand_ranked))
top10 <- cand_ranked %>% slice_head(n=topN)

out_top10 <- file.path(DIR_OUT, paste0("top_white_candidates_",STAMP,".csv"))
readr::write_csv(
  top10 %>% select(A,C,priority,best_score_q,ABBC_strength,
                   AB_evi_max,BC_evi_max,AB_evi_sum,BC_evi_sum,n_B),
  out_top10
)
message("Saved: ", out_top10)

############################################################
# (Optional) Quick preview tile plot (not used in main Figures)
############################################################

plot_df <- white_map %>%
  mutate(status = ifelse(AC_known,"known A–C","WHITE (A–C unknown)"))

p_white <- ggplot(plot_df, aes(x=C, y=A, fill=best_score_q)) +
  geom_tile(color="grey85") +
  geom_point(aes(size=ABBC_strength, shape=status), alpha=0.85) +
  scale_shape_manual(values=c("WHITE (A–C unknown)"=16, "known A–C"=1)) +
  labs(
    title="White map (pre-analysis)",
    x="Outcome (C)", y="Factor (A)"
  ) +
  theme_minimal(base_size=11) +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggsave(
  file.path(DIR_FIG, paste0("white_map_preanalysis_",STAMP,".pdf")),
  p_white, width=10, height=6
)

############################################################
# END
############################################################
