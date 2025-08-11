# GBM Plasticity Project

## Overview

Glioblastoma (GBM) is a highly aggressive brain tumor with poor prognosis, largely due to its exceptional capacity for relapse and resistance to therapy. **Cellular plasticity**—the ability of tumor cells to transition between transcriptional states—has emerged as a key driver of this adaptability ([Neftel et al., 2019](https://doi.org/10.1016/j.cell.2019.06.024)).

This project aims to **uncover the molecular mechanisms underlying GBM plasticity** to identify potential therapeutic targets. We apply a **multi-omics approach** combining **single-cell RNA sequencing (scRNA-seq)**, **bulk RNA-seq**, and **proteomics**.

---

## Objectives

1. **Model GBM heterogeneity** using patient-derived primary GBM tumor samples.
2. **Map transcriptional states** to the three clinical subtypes:
   - Proneural (PN)
   - Classical (CL)
   - Mesenchymal (MES)
3. **Identify genes and pathways** driving cellular state transitions.
4. **Integrate multi-omics data** to find regulators consistent at both RNA and protein levels.

---

## Methods

### 1. Single-Cell RNA-seq Analysis
- Performed scRNA-seq on three primary GBM samples.
- Confirmed high intra-tumoral heterogeneity and subtype representation (PN, CL, MES).
- Used **CellRank** for trajectory inference to identify dynamic gene expression changes.
- Discovered **CD44** (MES) and **CD24** (NPC-like) as surface markers for state enrichment.

### 2. Functional Multi-Omics Profiling (ongoing)
- Selected the most heterogeneous sample (based on CD24/CD44 expression) for detailed profiling.
- Performed **RNA-seq** and **proteomics profiling** on four CD24/CD44-sorted subpopulations.
- Integrated datasets to:
  - Identify RNA/protein co-regulated genes.
  - Reveal post-transcriptional regulation.
  - Highlight pathways linked to GBM plasticity.

### 3. Temporal State Transition Analysis
- Isolated **CD44⁻/CD24⁻** cells via FACS—capable of re-establishing parental heterogeneity within days.
- Collected cells at 0, 1, and 3 days post-sorting for scRNA-seq.
- Identified transcriptional programs and candidate regulators driving transitions, including:
  - **DLX gene family**
  - **MYC**

---

## References
1. Neftel C, Laffy J, Filbin MG, et al. *An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma.* Cell. 2019;178(4):835–849.e21. doi:[10.1016/j.cell.2019.06.024](https://doi.org/10.1016/j.cell.2019.06.024)
