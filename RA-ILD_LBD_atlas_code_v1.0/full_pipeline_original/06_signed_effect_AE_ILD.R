############################################################
# 06) Signed-effect analysis for AE-ILD — ORIGINAL
# ----------------------------------------------------------
# Reproduces Figure 5A (forest plot) and Figure 5B (bubble plot)
# from the RA-ILD manuscript.
#
# Inputs:
#   output/signed_effects_summary_withCI_{STAMP}.csv
#
# Outputs:
#   fig_pub/signed_effects_AE-ILD_{STAMP}.pdf
#   fig_pub/signed_effects_AE-ILD_bubble_{STAMP}.pdf
#
# This script is based on:
#   - signed_effects_AE-ILD.txt (旧版)
#   - 最終版 signed-effect 部分（RTF）
#
# Adjustments:
#   - ROOT generalized
#   - jpfont = "sans"
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"), pkgType="binary")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
  library(ggplot2); library(scales); library(forcats)
})

ROOT <- "/path/to/RA-ILD_textanalysis"
DIR_OUT <- file.path(ROOT, "output")
DIR_FIG <- file.path(ROOT, "fig_pub")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

jpfont <- "sans"
cairo_pdf_device <- grDevices::cairo_pdf

############################################################
# 1) Load latest signed-effect summary
############################################################

files <- list.files(
  DIR_OUT,
  pattern="^signed_effects_summary_withCI_\\d{8}_\\d{4}\\.csv$",
  full.names = TRUE
)
stopifnot(length(files) > 0)

# Pick the most recent one
f_csv <- files[which.max(file.info(files)$mtime)]

STAMP <- sub(
  "^.*signed_effects_summary_withCI_(\\d{8}_\\d{4})\\.csv$",
  "\\1", basename(f_csv)
)

DF <- readr::read_csv(f_csv, show_col_types = FALSE)
stopifnot(nrow(DF) > 0)

############################################################
# 2) Filter for AE-ILD (C == "AE-ILD")
############################################################

target_C <- "AE-ILD"

dat <- DF %>%
  filter(C == target_C) %>%
  mutate(
    articles     = coalesce(articles, 0),
    balance      = coalesce(balance, 0),
    balance_low  = coalesce(balance_low, NA_real_),
    balance_high = coalesce(balance_high, NA_real_),
    pos_ratio    = coalesce(pos_ratio, NA_real_),
    A_f          = fct_inorder(A)
  )

stopifnot(nrow(dat) > 0)

# Order by article count then by balance
dat <- dat %>%
  arrange(desc(articles), desc(balance)) %>%
  mutate(A_f = fct_inorder(A))

############################################################
# 3) Color scale (balance)
############################################################

col_fun <- scales::col_numeric(
  palette = c(
    "#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7",
    "#FDDBC7","#F4A582","#D6604D","#B2182B"
  ),
  domain = c(-1,1),
  na.color = "grey80"
)

############################################################
# 4) Figure 5A — Forest Plot
############################################################

p_forest <- ggplot(dat, aes(x = balance, y = A_f)) +
  geom_errorbarh(
    aes(xmin = balance_low, xmax = balance_high),
    height = 0, alpha = 0.7, na.rm = TRUE
  ) +
  geom_point(
    aes(size = articles, fill = balance),
    shape = 21, stroke = 0.3, color = "black"
  ) +
  scale_fill_gradientn(
    colours = col_fun(seq(-1,1,length.out=9)),
    limits  = c(-1,1), oob = scales::squish
  ) +
  scale_size_continuous(
    name="Articles", range=c(2.3,7),
    breaks = pretty_breaks(3)
  ) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.4) +
  coord_cartesian(xlim=c(-1,1)) +
  labs(
    title = "Signed effect vs AE-ILD (forest)",
    subtitle = paste0("Balance with CI | ", target_C, " | ", STAMP),
    x = "Balance (risk_decrease ← 0 → risk_increase)",
    y = NULL, fill="Balance"
  ) +
  theme_minimal(base_size = 12, base_family = jpfont) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

outfile_forest <- file.path(DIR_FIG, paste0("signed_effects_AE-ILD_", STAMP, ".pdf"))
ggsave(outfile_forest, p_forest, device=cairo_pdf_device, width=8.5, height=11)
message("Saved: ", outfile_forest)

############################################################
# 5) Figure 5B — Bubble Plot
############################################################

dat_bub <- dat %>% filter(!is.na(pos_ratio))

p_bubble <- ggplot(dat_bub, aes(x = pos_ratio, y = A_f)) +
  geom_point(
    aes(size = articles, fill = balance),
    shape = 21, color = "black", stroke = 0.3, alpha = 0.9
  ) +
  scale_fill_gradientn(
    colours = col_fun(seq(-1,1,length.out=9)),
    limits = c(-1,1), oob = scales::squish
  ) +
  scale_size_continuous(
    name = "Articles",
    range = c(2.3,7),
    breaks = pretty_breaks(3)
  ) +
  scale_x_continuous(
    limits = c(0,1),
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0,1,0.2)
  ) +
  labs(
    title = "Signed effect vs AE-ILD (bubble)",
    subtitle = paste0("Positive-article ratio | ", target_C, " | ", STAMP),
    x = "Positive articles ratio",
    y = NULL, fill="Balance"
  ) +
  theme_minimal(base_size = 12, base_family = jpfont) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

outfile_bubble <- file.path(DIR_FIG, paste0("signed_effects_AE-ILD_bubble_",STAMP,".pdf"))
ggsave(outfile_bubble, p_bubble, device=cairo_pdf_device, width=10, height=11)
message("Saved: ", outfile_bubble)

############################################################
# END
############################################################
