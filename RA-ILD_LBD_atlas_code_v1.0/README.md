# RA-ILD Literature-Based Discovery (LBD) & Biomarker Evidence Atlas  
### Full Reproducible Code for Manuscript Figures (RA-ILD_LBD_atlas_code_v1.0)

This repository provides all R scripts required to reproduce every figure in the RA-ILD text-mining manuscript, including:

- Topic modelling & t-SNE maps (Figure 2)  
- A↔C coherence analysis (Figure 3)  
- ABC triad network (Figure 4)  
- Signed-effect analysis (Figure 5A–B)  
- White-map decomposition & fallback2 (Figure 5C)  
- Biomarker Evidence Atlas (Figure 6A–B)

All scripts are organized into a reproducible pipeline with clear execution order.

---

## 🔧 Directory Structure

RA-ILD_LBD_atlas_code_v1.0/
│
├── README.md
│
├── full_pipeline_original/
│     ├── 01_corpus_textmining_core.R
│     ├── 02_topic_LDA_tSNE.R
│     ├── 03_whiteMap_preanalysis_core.R
│     ├── 04_coherence_AC_npmi_lift.R
│     ├── 05_ABC_triad_score_network.R
│     ├── 06_signed_effect_AE_ILD.R
│     ├── 07_whiteMap_Bbreakdown_fallback.R
│     └── 08_biomarker_evidence_atlas.R
│
├── core/
│     ├── S1_white_map_core.R
│     ├── S1_white_top10.R
│     └── S1_ac_coherence.R
│
└── dic/
  　  ├── raalid_pub.csv
　    ├── raalid_v07.csv
　    └── external_ac_evidence.csv




## ▶️ How to Reproduce Manuscript Figures

1. **Prepare directories**

/path/to/RA-ILD_textanalysis/
data_raw/
data_proc/
output/
fig_pub/



2. **Place input files**

- `articles_*.csv` → `data_proc/`
- Dictionary files → `dic/`

3. **Run scripts in this order (full_pipeline_original/):**

01_corpus_textmining_core.R
02_topic_LDA_tSNE.R
03_whiteMap_preanalysis_core.R
04_coherence_AC_npmi_lift.R
05_ABC_triad_score_network.R
06_signed_effect_AE_ILD.R
07_whiteMap_Bbreakdown_fallback2.R
08_biomarker_evidence_atlas.R


4. **Figures will be saved automatically** in:

fig_pub/


---

## 📌 Notes

- All scripts use `ROOT <- "/path/to/RA-ILD_textanalysis"` → please modify this path for your environment.
- Font is set to `family = "sans"` for maximal cross-platform compatibility.
- The “core” folder contains minimal reproducible code used in Supplementary Material.
- The “extended” folder contains optional visualizations that were helpful for analysis but are not part of the main figures.

---

## 📚 Citation

Please cite the associated manuscript when using these scripts.

A DOI for this repository can be issued through Zenodo once uploaded to GitHub.

---

## ❓ Contact

For questions regarding the analysis pipeline, dictionary, or reproducibility,  
please contact the authors.


