## ============================================================
## S1_white_top10.R   (Supplementary Code)
## Extract top-ranked unknown A–C pairs (“White Map candidates”)
## and generate a simple bar plot.
##
## NOTE:
## - Assumes that `white_map` has already been computed 
##   (see S1_white_map_core.R).
## - This is a simplified representation; full logic (e.g., fallback,
##   canonicalization) is available on GitHub.
## ============================================================

library(dplyr)
library(ggplot2)

## ---- Select unknown A–C pairs first ----
cand_unknown <- white_map %>% filter(!AC_known)
if (nrow(cand_unknown) == 0L) {
  cand_unknown <- white_map
}

## ---- Rank candidates using the same priority logic ----
cand_ranked <- cand_unknown %>%
  mutate(ABBC_sum = AB_evi_sum + BC_evi_sum) %>%
  arrange(
    desc(priority),
    desc(best_score_q),
    desc(ABBC_strength),
    desc(ABBC_sum),
    desc(n_B)
  )

## ---- Extract top 10 entries ----
topN  <- min(10L, nrow(cand_ranked))
top10 <- cand_ranked %>%
  slice_head(n = topN) %>%
  select(A, C, priority, best_score_q, ABBC_strength,
         AB_evi_max, BC_evi_max, AB_evi_sum, BC_evi_sum, n_B)

## ---- Simple bar plot (conceptual) ----
p_top10 <- ggplot(top10,
                  aes(x = priority,
                      y = reorder(paste0(A, " → ", C), priority))) +
  geom_col(width = 0.65) +
  labs(
    title = "Top White Map Candidates (A–C unknown)",
    x     = "Priority score",
    y     = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(ROOT, "fig", "top_white_candidates_bar_supplement.pdf"),
       p_top10,
       width = 7,
       height = max(3, 0.4 * nrow(top10) + 1))
