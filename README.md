# Ecological-drivers-of-gut-microbiota-diversity-and-evolution-in-brazilian-stingless-bees
Here, we characterized the gut microbiota of 29 species of Brazilian stingless bees collected across four states (Minas Gerais, Espírito Santo, Amazonas, and Pernambuco), encompassing a broad range of ecological traits, nesting behaviors, and management conditions (wild and managed colonies).
# Gut Microbiota Analyses in Brazilian stingless bees
[![DOI](https://zenodo.org/badge/1157461999.svg)](https://doi.org/10.5281/zenodo.22864094)
This repository contains all scripts and workflows used for the microbiome analyses described in our manuscript.  
All analyses were conducted in **R v4.5.1** and rely on a set of specialized packages for sequence processing, diversity analyses, differential abundance testing, functional prediction, and network inference.

---

## 📂 Repository Contents
- `scripts/` – R scripts for each analysis step
- `data/` – Input files (ASV tables, metadata, taxonomy assignments)
- `results/` – Output files (figures, tables, statistical summaries)
- `README.md` – Documentation (this file)

---

## 🔧 Software and Packages

All analyses were run in **R v4.5.1 (R Core Team, 2025)**.  
Below is the list of packages and versions referenced in the workflow:

| Step | Package | Version | Reference |
|------|---------|---------|-----------|
| Sequence processing | **DADA2** | v1.52.0 | Callahan et al., 2016 |
| Primer removal | **Cutadapt** | (external tool) | Martin, 2011 |
| Data integration | **phyloseq** | v1.52.0 | McMurdie & Holmes, 2012 |
| Sequence alignment | **DECIPHER** | v3.4.0 | Wright, 2016 |
| Phylogenetic tree | **phangorn** | (functions: NJ, optim.pml) | Schliep, 2011 |
| Diversity analysis | **vegan** | v2.7-2 | Oksanen et al., 2025 |
| Post-hoc tests | **emmeans** | v2.0.0 | Lenth & Piaskowski, 2017 |
| Differential abundance | **DESeq2** | v1.48.1 | Love et al., 2014 |
| Differential abundance | **ANCOM-BC2** | v2.10.1 | Lin & Peddada, 2020 |
| Visualization | **ggplot2** | (latest CRAN) | Wickham, 2016 |
| Alluvial plots | **ggalluvial** | v0.12.5 | Brunson, 2020 |
| Functional prediction | **PICRUSt2** | v2.3.0-b | Douglas et al., 2020 |
| KEGG enrichment | **clusterProfiler** | v4.16.0 | Yu et al., 2012 |
| KEGG REST API | **KEGGREST** | v1.48.1 | Tenenbaum, 2025 |
| Network inference | **SpiecEasi** | v1.1.13 | Kurtz et al., 2015 |
| Network inference | **SparCC** | (via SpiecEasi) | Friedman & Alm, 2012 |
| Network visualization | **igraph** | v2.2.1 | Csárdi et al., 2025 |
| Cytoscape interface | **RCy3** | v2.28.1 | Shannon et al., 2003 |
| Statistical tests | **rstatix** | v0.7.3 | Kassambara, 2020 |
| Tree visualization | **ggtree** | (Bioconductor) | Yu et al., 2017 |

📖 Citation
If you use this code, please cite our manuscript and the relevant software packages listed above.

📜 License
This repository is released under the MIT License. You are free to use, modify, and distribute the code with proper attribution.
