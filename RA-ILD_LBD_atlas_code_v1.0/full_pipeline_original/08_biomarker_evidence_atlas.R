############################################################
# 08) Biomarker Evidence Atlas — ORIGINAL
# ----------------------------------------------------------
# Reproduces Figure 6A and Figure 6B of the RA-ILD manuscript.
#
# Inputs:
#   - signed_effects_summary_withCI_*.csv
#   - cooc_npmi_lift_*.csv
#   - raalid_pub.csv / raalid_v07.csv (dictionary)
#   - external_ac_evidence.csv (if present)
#
# Outputs:
#   fig_pub/atlas_biomarker_heatmap_{STAMP}.pdf   (Figure 6A)
#   fig_pub/atlas_biomarker_network_{STAMP}.pdf   (Figure 6B)
#
# Origin:
#   This script is the cleaned integration of the FINAL code block
#   in "Figure調整用追加コード.rtf".
#
# Adjustments:
#   - ROOT generalized
#   - jpfont = "sans"
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
  library(ggplot2); library(forcats); library(glue); library(igraph)
  library(ggraph); library(purrr)
})

ROOT <- "/path/to/RA-ILD_textanalysis"
DIR_OUT <- file.path(ROOT,"output")
DIR_FIG <- file.path(ROOT,"fig_pub")
DIR_DIC <- file.path(ROOT,"dic")

dir.create(DIR_FIG, recursive=TRUE, showWarnings=FALSE)

jpfont <- "sans"
cairo_pdf_device <- grDevices::cairo_pdf
STAMP <- format(Sys.time(), "%Y%m%d_%H%M")

############################################################
# 1) Load dictionary (for selecting biomarker classes)
############################################################

dic_pub  <- read_csv(file.path(DIR_DIC,"raalid_pub.csv"), show_col_types=FALSE)
dic_v07  <- read_csv(file.path(DIR_DIC,"raalid_v07.csv"), show_col_types=FALSE)
dic      <- bind_rows(dic_pub, dic_v07) %>% distinct(term, .keep_all=TRUE)

############################################################
# 2) Load signed-effect summary WITH CI
############################################################

f_se <- list.files(
  DIR_OUT,
  pattern="^signed_effects_summary_withCI_\\d{8}_\\d{4}\\.csv$",
  full.names=TRUE
) |> sort() |> tail(1)

stopifnot(length(f_se)==1)
SE <- read_csv(f_se, show_col_types=FALSE)

############################################################
# 3) Load A–C npmi/lift table
############################################################

f_npmi <- list.files(
  DIR_OUT,
  pattern="^cooc_npmi_lift_\\d{8}_\\d{4}\\.csv$",
  full.names=TRUE
) |> sort() |> tail(1)

stopifnot(length(f_npmi)==1)
NP <- read_csv(f_npmi, show_col_types=FALSE)

############################################################
# 4) Select biomarker-like A terms (role in B or molecular/cell)
############################################################

# Classes that constitute biomarker/family-level entries
BM_CLASSES <- c("biomarker","molecular","microbiome","cell")

# Which A-terms are biomarker-like?
BM_A <- dic %>%
  filter(class %in% BM_CLASSES) %>%
  pull(term) %>% unique()

# Also include family-level A created during dictionary building
FAMILY_A <- dic %>%
  filter(str_detect(class,"family") | str_detect(term,"_pathway")) %>%
  pull(term) %>% unique()

A_KEEP <- unique(c(BM_A, FAMILY_A))

############################################################
# 5) Join signed-effects + npmi/lift for A–C pairs
############################################################

# Signed-effects table:
#   A, C, balance, balance_low, balance_high, articles, pos_ratio, ...
# NP table:
#   A, C, pA, pC, pAC, lift, npmi

AC_MERGE <- SE %>%
  filter(A %in% A_KEEP) %>%
  left_join(
    NP %>% select(A,C,pA,pC,pAC,lift,npmi),
    by=c("A","C")
  )

############################################################
# 6) Define grouped outcomes for atlas (AE-ILD, progression, mortality)
############################################################

group_C <- function(c){
  case_when(
    c == "AE-ILD"                                ~ "AE-ILD",
    c %in% c("progression","FVC_decline","DLCO_decline") ~ "Progression",
    c == "mortality"                              ~ "Mortality",
    TRUE                                          ~ "Other"
  )
}

ACG <- AC_MERGE %>%
  mutate(Cgrp = group_C(C))

############################################################
# 7) Compute atlas score (RTF final logic)
# ----------------------------------------------------------
# scoreraw = (|balance| + 0.20) * (evidence + 0.50) * (max(0,npmi) + 0.05)
# evidence = log10(max(1, articles))
#
# If ext_flag == TRUE → atlas_score = scoreraw * 0.5
############################################################

ACG <- ACG %>%
  mutate(
    evidence = log10(pmax(1, articles)),
    scoreraw = (abs(balance)+0.20) * (evidence+0.50) * (pmax(0,npmi)+0.05)
  )

# External evidence flag（あれば処理）
f_ext <- file.path(DIR_DIC,"external_ac_evidence.csv")
if (file.exists(f_ext)){
  EXT <- read_csv(f_ext, show_col_types=FALSE)
  ACG <- ACG %>% left_join(EXT, by=c("A","Cgrp"))
  ACG <- ACG %>%
    mutate(
      ext_flag = coalesce(ext_flag, FALSE),
      atlas_score = ifelse(ext_flag, scoreraw*0.5, scoreraw)
    )
} else {
  ACG <- ACG %>% mutate(ext_flag = FALSE, atlas_score = scoreraw)
}

############################################################
# 8) Family labeling (RTF final)
############################################################

family_of <- function(a){
  dic_row <- dic %>% filter(term==a)
  if (!nrow(dic_row)) return(a)
  if (!is.na(dic_row$family[1])) return(dic_row$family[1])
  return(a)
}

ACG <- ACG %>%
  mutate(
    family = purrr::map_chr(A, family_of),
    A_label = ifelse(family!=A, paste0(family," (family)"), A)
  )

############################################################
# 9) Prepare Figure 6A bubble heatmap
############################################################

# Row order = atlas_score (descending)
row_order <- ACG %>%
  group_by(A_label) %>%
  summarise(max_score = max(atlas_score, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(max_score)) %>% pull(A_label)

col_order <- c("AE-ILD","Progression","Mortality")  #論文と同順

DF6A <- ACG %>%
  filter(Cgrp %in% col_order) %>%
  mutate(
    A_label = factor(A_label, levels = row_order),
    Cgrp    = factor(Cgrp,    levels = col_order)
  )

# bubble fill = balance, bubble size = evidence
p6A <- ggplot(DF6A, aes(x=Cgrp, y=A_label)) +
  geom_point(
    aes(fill=balance, size=evidence),
    shape=21, color="black", alpha=0.9
  ) +
  scale_fill_gradient2(
    low="#2166AC", mid="white", high="#B2182B",
    midpoint=0, name="Balance"
  ) +
  scale_size(range=c(2,12), name="log10 evidence") +
  labs(
    title="Biomarker evidence atlas (direction vs evidence)",
    x="Outcome group (C)", y="Biomarker / Family"
  ) +
  theme_minimal(base_size=12) +
  theme(
    text=element_text(family=jpfont),
    axis.text.x=element_text(angle=45, hjust=1)
  )

f6A <- file.path(DIR_FIG, paste0("atlas_biomarker_heatmap_",STAMP,".pdf"))
ggsave(f6A, p6A, device=cairo_pdf_device, width=11, height=9)
message("Saved (Fig6A): ", f6A)

############################################################
# 10) Figure 6B — network (biomarker/family → outcome)
############################################################

# Select high-scoring edges
EDGE <- ACG %>%
  filter(Cgrp %in% col_order) %>%
  group_by(A_label, Cgrp) %>%
  summarise(
    score = max(atlas_score, na.rm=TRUE),
    balance = mean(balance, na.rm=TRUE),
    .groups="drop"
  ) %>%
  filter(score > quantile(score, 0.75, na.rm=TRUE)) # top quartile edges

# Build graph
nodes_terms <- unique(c(EDGE$A_label, EDGE$Cgrp))
nodes <- tibble(name = nodes_terms) %>%
  mutate(
    type = ifelse(name %in% col_order, "Outcome", "Biomarker")
  )

g <- igraph::graph_from_data_frame(
  d = EDGE %>% select(A_label, Cgrp, score, balance),
  directed = FALSE,
  vertices = nodes
)

# Node aesthetics
V(g)$color <- ifelse(V(g)$type=="Outcome", "#7B241C", "#1F618D")
V(g)$size  <- ifelse(V(g)$type=="Outcome", 12, 8)
V(g)$label <- V(g)$name

# Edge aesthetics
E(g)$width <- scales::rescale(E(g)$score, c(0.6,4))
E(g)$color <- scales::col_numeric(
  palette=c("#2166AC","#F7F7F7","#B2182B"),
  domain=c(-1,1)
)(E(g)$balance)

# Plot
set.seed(123)
p6B <- ggraph::ggraph(g, layout="fr") +
  ggraph::geom_edge_link(aes(color=I(E(g)$color), width=I(E(g)$width)), alpha=0.9) +
  ggraph::geom_node_point(aes(color=I(V(g)$color), size=I(V(g)$size))) +
  ggraph::geom_node_text(aes(label=label), family=jpfont, repel=TRUE, size=3) +
  labs(
    title="Biomarker/Family – Outcome network",
    subtitle="Edges weighted by atlas score; colored by balance"
  ) +
  theme_void(base_size=12) +
  theme(text=element_text(family=jpfont))

f6B <- file.path(DIR_FIG, paste0("atlas_biomarker_network_",STAMP,".pdf"))
ggsave(f6B, p6B, device=cairo_pdf_device, width=11, height=8)
message("Saved (Fig6B): ", f6B)

############################################################
# END
############################################################
