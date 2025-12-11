## ============================================================
## S1_white_map_core.R   (Supplementary Code)
## Core logic for RA-ILD text-mining — White Map (A–C known/unknown)
## This script provides a minimal reproducible example of how 
## A–C known relationships are detected and how prioritization 
## of unknown A–C pairs (“White Map”) is computed.
##
## NOTE:
## - This is a conceptual, minimal version for Supplementary Material.
## - Full implementation (including canonicalization, dictionary 
##   normalization, fallback rules, and plotting adjustments) is 
##   available on the GitHub repository referenced in the manuscript.
## - Readers should modify file paths according to their environment.
## ============================================================

## ---- Load required packages ----
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)

## ---- Project root path (edit as needed) ----
## Examples:
## ROOT <- "C:/Users/username/projects/RA-ILD_textanalysis"
## ROOT <- "/home/username/projects/RA-ILD_textanalysis"
ROOT <- "/path/to/RA-ILD_textanalysis"   # ★ readers must modify

## ---- Input files (examples) ----
hits_file <- file.path(ROOT, "data_proc", "hits_matrix_YYYYMMDD.csv")
abc_file  <- file.path(ROOT, "output",    "abc_rankings_YYYYMMDD_HHMM.csv")

hits <- readr::read_csv(hits_file, show_col_types = FALSE)
abc  <- readr::read_csv(abc_file,  show_col_types = FALSE)

## ---- Define A and C term candidates ----
## (In the full pipeline these are derived from dictionaries + canonicalization)
hit_cols  <- names(hits)[startsWith(names(hits), "hit__")]
terms_all <- sub("^hit__", "", hit_cols)

C_terms <- unique(abc$C)               # outcomes (AE-ILD, progression, etc.)
A_terms <- setdiff(terms_all, C_terms) # all remaining terms considered A-side

## ---- Utility for binary hit vectors ----
binv <- function(v) as.integer(ifelse(is.na(v), 0L, v > 0L))

get_hit <- function(df, term){
  col <- paste0("hit__", term)
  if (!col %in% names(df)) return(integer(nrow(df)))
  binv(df[[col]])
}

## ---- Determine whether A–C has ever co-occurred (AC_known) ----
ac_tbl <- lapply(A_terms, function(a){
  tibble::tibble(
    A = a,
    !!!setNames(lapply(C_terms, function(c){
      sum(get_hit(hits, a) & get_hit(hits, c))
    }),
    paste0("AC_n11__", C_terms))
  )
}) |> bind_rows()

ac_long <- ac_tbl |>
  tidyr::pivot_longer(starts_with("AC_n11__"),
                      names_to = "C", values_to = "AC_n11") |>
  mutate(
    C        = sub("^AC_n11__", "", C),
    AC_known = AC_n11 > 0     # A–C pair is “known” if they co-occur in ≥1 abstract
  )

## ---- Summarize AB/BC evidence strength and score_q by B ----
abc_summ <- abc |>
  group_by(A, C) |>
  summarise(
    best_score_q = max(score_q, na.rm = TRUE),
    AB_evi_max   = max(AB_n11,   na.rm = TRUE),
    BC_evi_max   = max(BC_n11,   na.rm = TRUE),
    AB_evi_sum   = sum(AB_n11,   na.rm = TRUE),
    BC_evi_sum   = sum(BC_n11,   na.rm = TRUE),
    n_B          = n(),
    .groups      = "drop"
  )

## ---- Construct White Map prioritization table ----
white_map <- ac_long |>
  left_join(abc_summ, by = c("A","C")) |>
  mutate(
    ABBC_strength = pmax(AB_evi_max, 0) + pmax(BC_evi_max, 0),
    priority      = dplyr::coalesce(best_score_q, 0) + 0.1 * ABBC_strength
  ) |>
  arrange(desc(!AC_known), desc(priority))

## ---- Simplified visualization (conceptual) ----
plot_df <- white_map |>
  mutate(status = ifelse(AC_known, "known A–C", "WHITE (A–C unknown)"))

p_white <- ggplot(plot_df, aes(x = C, y = A, fill = best_score_q)) +
  geom_tile(color = "grey85") +
  geom_point(aes(size = ABBC_strength, shape = status), alpha = 0.85) +
  scale_shape_manual(values = c("WHITE (A–C unknown)" = 16,
                                "known A–C"          = 1)) +
  labs(
    title = "White Map (RA-ILD): A–C known/unknown and AB/BC evidence",
    x     = "Outcome (C)",
    y     = "Candidate factor (A)",
    fill  = "best score_q",
    size  = "AB/BC evidence",
    shape = "A–C status"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(ROOT, "fig", "white_map_AC_tiles_supplement.pdf"),
       p_white,
       width = 8,
       height = 0.25 * length(unique(plot_df$A)) + 3)
