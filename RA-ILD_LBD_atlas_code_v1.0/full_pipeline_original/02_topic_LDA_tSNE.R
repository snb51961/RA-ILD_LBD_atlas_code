############################################################
# 02) Topic modelling (LDA) + t-SNE mapping  — ORIGINAL
# ----------------------------------------------------------
# Reconstructs Figure 2 of the RA-ILD manuscript:
#   - LDA top terms panel
#   - t-SNE topic map (clusters labeled T1–T8)
#   - topic × year trend plots
#
# Based on the full "図の作成コード topicmap.txt" (final version)
# with these adjustments:
#   - ROOT generalized
#   - font family forced to "sans"
#   - topic labels simplified to "T1"…"T8"
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(tibble); library(ggplot2); library(quanteda)
  library(topicmodels); library(Rtsne); library(tidytext)
  library(ggrepel)
})

# ---------------------------
# Project paths
# ---------------------------
ROOT <- "/path/to/RA-ILD_textanalysis"
DIR_PUB <- file.path(ROOT, "fig_pub")
DIR_OUT <- file.path(ROOT, "output")
dir.create(DIR_PUB, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

STAMP_TIME <- format(Sys.time(), "%Y%m%d_%H%M")

jpfont <- "sans"
cairo_pdf_device <- grDevices::cairo_pdf

# ---------------------------
# Utility — pick latest file
# ---------------------------
pick_latest <- function(pattern, dirs){
  xs <- unlist(lapply(dirs, function(d) Sys.glob(file.path(d, pattern))))
  if (!length(xs)) return(NA_character_)
  xs[which.max(file.info(xs)$mtime)]
}

# ---------------------------
# Load latest articles_*.csv
# ---------------------------
f_articles <- pick_latest(
  "articles_*.csv",
  c(file.path(ROOT,"data_proc"), file.path(ROOT,"data_raw"))
)
stopifnot(!is.na(f_articles), file.exists(f_articles))

articles <- readr::read_csv(f_articles, show_col_types = FALSE) %>%
  mutate(
    text = coalesce(text, paste(title, abstract, sep=" ")),
    text = ifelse(is.na(text), "", text)
  )

if (nrow(articles) < 30)
  stop("Too few documents (<30). Increase corpus or relax filters.")

# ---------------------------
# Tokenization & DFM
# ---------------------------
sw_general <- quanteda::stopwords("en")

sw_domain_base <- c(
  "rheumatoid","arthritis","interstitial","lung","disease","ild","ip","pulmonary","fibrosis",
  "patients","patient","study","studies","review","case","cases","report","reports",
  "introduction","conclusion","methods","background","objective","aim","result","results",
  "purpose","mg","ml","day","days","week","weeks","year","years"
)

sw_domain_extra <- c(
  "clinical","associated","factors","group","groups","levels","activity","use","used","using",
  "evidence","management","systemic","diseases","connective","tissue","months","month","ci",
  "significant","significance","primary","syndrome","pneumonia"
)

sw_domain <- unique(c(sw_general, sw_domain_base, sw_domain_extra))

corp <- quanteda::corpus(articles, text_field = "text")

toks <- corp %>%
  quanteda::tokens(remove_punct=TRUE, remove_numbers=TRUE, remove_symbols=TRUE) %>%
  quanteda::tokens_tolower() %>%
  quanteda::tokens_remove(sw_domain)

dfm_raw <- quanteda::dfm(toks)

# Stepwise trimming
dfm_trimmed <- quanteda::dfm_trim(dfm_raw, min_docfreq = 5, docfreq_type = "count")
dfm_trimmed <- quanteda::dfm_trim(dfm_trimmed, max_docfreq = 0.5, docfreq_type = "prop")

if (quanteda::nfeat(dfm_trimmed) < 100){
  dfm_trimmed <- quanteda::dfm_trim(dfm_raw, min_docfreq = 3, docfreq_type = "count")
  dfm_trimmed <- quanteda::dfm_trim(dfm_trimmed, max_docfreq = 0.8, docfreq_type = "prop")
}

stopifnot(quanteda::nfeat(dfm_trimmed) > 0, quanteda::ndoc(dfm_trimmed) > 0)

# ---------------------------
# LDA (K=8)
# ---------------------------
set.seed(42)
K <- 8
dtm <- quanteda::convert(dfm_trimmed, to = "topicmodels")

lda_fit <- topicmodels::LDA(
  dtm, k = K, method = "Gibbs",
  control = list(burnin=500, iter=2000, thin=50, seed=42)
)

# ---------------------------
# Extract β (top terms)
# ---------------------------
beta <- tidytext::tidy(lda_fit, matrix = "beta")

top_terms <- beta %>%
  group_by(topic) %>%
  slice_max(order_by = beta, n = 12, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(term_plot = tidytext::reorder_within(term, beta, topic))

p_topics <- ggplot(top_terms, aes(x = beta, y = term_plot)) +
  geom_col() +
  tidytext::scale_y_reordered() +
  facet_wrap(~ paste0("Topic ", topic), scales = "free_y", ncol = 2) +
  labs(
    title = paste0("LDA Top Terms (K=", K, ")"),
    x = "β (term | topic)", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = jpfont),
        strip.text = element_text(face = "bold"))

f_lda <- file.path(DIR_PUB, paste0("lda_topics_", STAMP_TIME, ".pdf"))
ggsave(f_lda, p_topics, device = cairo_pdf_device, width = 10, height = 8)
message("Saved: ", f_lda)

# ---------------------------
# γ matrix → t-SNE
# ---------------------------
gamma <- topicmodels::posterior(lda_fit)$topics
stopifnot(nrow(gamma) == quanteda::ndoc(dfm_trimmed))

doc_ids <- quanteda::docnames(dfm_trimmed)
idx <- suppressWarnings(as.integer(doc_ids))
if (any(is.na(idx)) || length(idx) != nrow(gamma))
  idx <- seq_len(nrow(gamma))

doc_meta <- articles[idx, c("pmid","year","journal","title")]

year_bin <- cut(
  doc_meta$year,
  breaks = c(-Inf,2004,2009,2014,2019,Inf),
  labels = c("<=2004","2005-09","2010-14","2015-19","2020-"),
  right = TRUE
)

set.seed(42)
perp <- max(5, floor(nrow(gamma)/50))

tsne_res <- Rtsne::Rtsne(
  gamma, perplexity = perp, theta = 0.5,
  check_duplicates = FALSE, pca = TRUE, max_iter = 1000
)

tsne_df <- tibble(
  x = tsne_res$Y[,1],
  y = tsne_res$Y[,2],
  year = doc_meta$year,
  year_bin = year_bin,
  pmid = doc_meta$pmid,
  top_topic = max.col(gamma, ties.method = "first")
)

# ---------------------------
# Topic labels = "T1" … "T8"
# ---------------------------
label_df <- tibble(
  topic = 1:K,
  label = paste0("T", topic)
)

centers <- tsne_df %>%
  group_by(top_topic) %>%
  summarise(
    cx = median(x), cy = median(y), n = dplyr::n(),
    .groups = "drop"
  ) %>%
  rename(topic = top_topic) %>%
  left_join(label_df, by = "topic")

# ---------------------------
# t-SNE scatter (with T1–T8 labels)
# ---------------------------
p_tsne <- ggplot(tsne_df, aes(x = x, y = y)) +
  geom_point(
    aes(color = factor(top_topic), shape = year_bin),
    size = 2, alpha = 0.9
  ) +
  ggrepel::geom_label_repel(
    data = centers,
    aes(x = cx, y = cy, label = label),
    size = 3.5, fontface = "bold", fill = "white", label.size = 0.2,
    segment.color = "grey50", max.overlaps = Inf
  ) +
  labs(
    title = paste0("t-SNE of LDA Document–Topic Mixtures (K=", K, ")"),
    x = "t-SNE 1", y = "t-SNE 2",
    color = "Dominant topic", shape = "Year bin"
  ) +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = jpfont))

f_tsne <- file.path(DIR_PUB, paste0("topicmap_", STAMP_TIME, ".pdf"))
ggsave(f_tsne, p_tsne, device = cairo_pdf_device, width = 9, height = 7)
message("Saved: ", f_tsne)

# ---------------------------
# Topic × Year: trends / heatmap / dominant topic counts
# ---------------------------
doc_year <- doc_meta$year
topic_year <- as.data.frame(gamma) %>%
  mutate(year = doc_year) %>%
  filter(!is.na(year)) %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "topic", values_to = "gamma") %>%
  mutate(topic = as.integer(stringr::str_extract(topic,"\\d+"))) %>%
  group_by(year, topic) %>%
  summarise(mean_gamma = mean(gamma), n_docs=n(), .groups="drop")

# Trend plot
p_trend <- ggplot(topic_year, aes(x = year, y = mean_gamma)) +
  geom_line() + geom_point(size=.8) +
  facet_wrap(~ paste0("T",topic), ncol=4, scales="free_y") +
  labs(title="Yearly mean topic proportion", x="Year", y="Mean γ") +
  theme_minimal(base_size=12) +
  theme(text = element_text(family=jpfont),
        strip.text = element_text(face="bold"))

ggsave(file.path(DIR_PUB, paste0("topic_year_trends_",STAMP_TIME,".pdf")),
       p_trend, device=cairo_pdf_device, width=12, height=7)

# Heatmap
p_heat <- ggplot(topic_year, aes(x = year, y = factor(topic), fill = mean_gamma)) +
  geom_tile() +
  scale_y_discrete(labels = function(v) paste0("T",v)) +
  labs(title="Topic × Year heatmap", x="Year", y="Topic", fill="Mean γ") +
  theme_minimal(base_size=12) +
  theme(text=element_text(family=jpfont))

ggsave(file.path(DIR_PUB, paste0("topic_year_heatmap_",STAMP_TIME,".pdf")),
       p_heat, device=cairo_pdf_device, width=12, height=6)

# Dominant topic counts
dom_df <- tibble(
  year = doc_year,
  topic = max.col(gamma, ties.method = "first")
) %>%
  filter(!is.na(year)) %>%
  count(year, topic, name="n")

p_dom <- ggplot(dom_df, aes(x = year, y = n, color=factor(topic))) +
  geom_line() + geom_point(size=.9) +
  labs(title="Counts of documents whose dominant topic = Tk",
       x="Year", y="Count", color="Topic") +
  theme_minimal(base_size=12) +
  theme(text=element_text(family=jpfont))

ggsave(file.path(DIR_PUB, paste0("topic_dominant_counts_",STAMP_TIME,".pdf")),
       p_dom, device=cairo_pdf_device, width=12, height=6)

# Save CSV summaries
write_csv(topic_year,
          file.path(DIR_OUT,paste0("topic_year_mean_gamma_",STAMP_TIME,".csv")))
write_csv(dom_df,
          file.path(DIR_OUT,paste0("topic_dominant_counts_",STAMP_TIME,".csv")))

message("=== LDA + t-SNE COMPLETE ===")
