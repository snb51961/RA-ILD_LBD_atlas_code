## ============================================================
## S1_ac_coherence.R   (Supplementary Code)
## A↔C coherence visualization using NPMI and log2(Lift)
##
## NOTE:
## - This is a simplified version of the full RA-ILD pipeline.
## - In the full code, canonicalization, thresholding, and 
##   visualization refinements are implemented.
## ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)

ROOT <- "/path/to/RA-ILD_textanalysis"  # ★ readers must modify
DIR_OUT <- file.path(ROOT, "output")

## ---- Input file (example) ----
npmi_file <- file.path(DIR_OUT, "cooc_npmi_lift_YYYYMMDD.csv")
cooc <- readr::read_csv(npmi_file, show_col_types = FALSE)

## ---- Ensure required columns exist ----
stopifnot(all(c("A","C","pA","pC","pAC","lift","npmi") %in% names(cooc)))

## ---- Basic transformation ----
D <- cooc %>%
  mutate(
    A         = as.character(A),
    C         = as.character(C),
    log2_lift = log2(ifelse(lift > 0, lift, NA_real_))
  ) %>%
  filter(!is.na(npmi), !is.na(log2_lift))

## ---- Select top A per C for visibility ----
TOP_PER_C <- 50L

D <- D %>%
  group_by(C) %>%
  arrange(desc(npmi), desc(log2_lift), .by_group = TRUE) %>%
  mutate(rank_c = row_number()) %>%
  ungroup() %>%
  filter(rank_c <= TOP_PER_C)

## ---- Order terms by mean rank ----
ord_A <- D %>%
  group_by(A) %>%
  summarise(mr = mean(rank_c, na.rm = TRUE), .groups = "drop") %>%
  arrange(mr) %>% pull(A)

ord_C <- D %>%
  count(C, name = "n") %>%
  arrange(desc(n)) %>% pull(C)

D <- D %>%
  mutate(
    A = factor(A, levels = ord_A),
    C = factor(C, levels = ord_C)
  )

## ============================================================
## 1) Heatmap for NPMI & log2(Lift)
## ============================================================

clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

D_heat_npmi <- D %>%
  mutate(
    npmi_clip = clip(npmi, -0.2, 1.0),
    metric    = factor("NPMI", levels = c("NPMI","log2(Lift)")),
    value     = npmi_clip
  )

D_heat_l2l <- D %>%
  mutate(
    l2l_clip  = clip(log2_lift, -2, 2),
    metric    = factor("log2(Lift)", levels = c("NPMI","log2(Lift)")),
    value     = l2l_clip
  )

HH <- bind_rows(D_heat_npmi, D_heat_l2l)

p_heat <- ggplot(HH, aes(x = A, y = C, fill = value)) +
  geom_tile() +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_fill_gradient2(
    low  = "#4575b4",
    mid  = "#ffffbf",
    high = "#d73027",
    midpoint = 0,
    name = "Value"
  ) +
  labs(
    title    = "A–C Coherence Heatmap (NPMI & log2(Lift))",
    subtitle = basename(npmi_file),
    x        = "A term",
    y        = "Outcome (C)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
    panel.grid  = element_blank(),
    strip.text  = element_text(face = "bold")
  )

ggsave(file.path(ROOT, "fig_pub",
                 "ac_npmi_lift_heatmap_supplement.pdf"),
       p_heat, width = 9, height = 7)

## ============================================================
## 2) Scatter: NPMI vs log2(Lift)
## ============================================================

th_npmi <- 0.05
th_l2l  <- log2(1.20)

LABEL_PER_C <- 8L

lab_df <- D %>%
  group_by(C) %>%
  arrange(desc(npmi + 0.5 * log2_lift), .by_group = TRUE) %>%
  slice_head(n = LABEL_PER_C) %>%
  ungroup()

p_sc <- ggplot(D, aes(x = npmi, y = log2_lift)) +
  geom_hline(yintercept = th_l2l, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = th_npmi, linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(size = pAC),
             alpha = 0.8) +
  geom_text(
    data = lab_df,
    aes(label = A),
    size = 3,
    nudge_y = 0.03,
    check_overlap = TRUE
  ) +
  facet_wrap(~ C, scales = "free") +
  scale_size_continuous(
    name  = "p(A,C) (co-occurrence rate)",
    range = c(1, 5)
  ) +
  labs(
    title    = "A–C Coherence Scatter (NPMI vs log2(Lift))",
    subtitle = sprintf("Guides: NPMI ≥ %.2f, Lift ≥ 1.2", th_npmi),
    x        = "NPMI",
    y        = "log2(Lift)"
  ) +
  coord_cartesian(ylim = c(-0.5, 2.5)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(ROOT, "fig_pub",
                 "ac_npmi_lift_scatter_supplement.pdf"),
       p_sc, width = 11, height = 7)
