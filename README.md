# Viral Metagenomics & Host-Virus Interaction Pipeline 🦠🧬

A comprehensive and highly robust pipeline for viral discovery, quality control, abundance quantification, and host-virus interaction network reconstruction from metagenomic assemblies. Developed for **EliLab**.

This workflow integrates state-of-the-art tools to reliably identify viral contigs, pair them with their bacterial/archaeal hosts using both empirical (CRISPR/BLAST) and machine-learning (iPHoP) evidence, and visualize the complex ecological dynamics of the community.

## 📊 Pipeline Workflow

The pipeline consists of modular steps, going from raw assemblies to publication-ready visualizations.

```mermaid
graph TD
    %% Estilos de nodos
    classDef process fill:#f9f,stroke:#333,stroke-width:2px;
    classDef filter fill:#bbf,stroke:#333,stroke-width:2px;
    classDef output fill:#bfb,stroke:#333,stroke-width:2px;
    classDef rscript fill:#ff9,stroke:#333,stroke-width:2px;

    %% 1. ENTRADA DE DATOS
    A[Metagenomic Assemblies / Raw Data] --> B(Viral Discovery Pipeline)
    A --> C(MAG Binning Pipeline)
    
    %% 2. FILTRADO Y DESCUBRIMIENTO VIRAL
    B -->|geNomad, VirSorter2, VIBRANT| D[Putative Viral Contigs]
    D --> F{CheckV Quality Filter}
    F -->|>= 5kb & Contig Quality| G[Filtered Viral Contigs]
    
    %% INCLUSIÓN DE VRHYME
    G --> V[vRhyme: Viral Binning & vMAG Reconstruction]
    V --> W[Linked vMAGs & Unbinned Contigs]
    W --> CheckV2{CheckV End-to-End Quality Check}
    CheckV2 -->|Completeness / Quality Assessment| H_vOTU[Dereplication & Final vOTUs]

    %% 3. PROCESAMIENTO DE HOST MAGs
    C -->|metaBAT2, VAMB, etc.| E[Raw MAGs]
    E --> H{CheckM2 Filter}
    H -->|Completeness >=50%, Contamination <=10%| I[High-Quality MAGs]
    I --> J(Taxonomic Classification: GTDB-Tk)

    %% 4. ANOTACIÓN FUNCIONAL (DRAM-V)
    H_vOTU --> DRAM[DRAM-v: Functional & AMG Annotation]
    DRAM --> AMGs[Verified Auxiliary Metabolic Genes / AMGs: phoH, Sulfur, Nucleotides]

    %% 5. ABUNDANCIA E INTERACCIONES
    H_vOTU --> K(Abundance Quantification)
    I --> K
    K -->|CoverM TPM (Lane Averaging)| L[Robust Abundance Matrix]

    H_vOTU --> M(Host-Virus Network Discovery)
    I --> M
    M -->|MinCED| N[CRISPR Spacers]
    M -->|BLASTn| O[Provirus Integration / Homology]
    M -->|iPHoP v1.4.2| P[Machine Learning Host Prediction]

    %% 6. INTEGRACIÓN Y ESTADÍSTICAS EN R
    L --> Q((Master R Script))
    J --> Q
    N --> Q
    O --> Q
    P --> Q
    AMGs --> Q

    %% 7. PRODUCTOS Y RESULTADOS VISUALES
    Q --> R[Host-Virus Interaction Networks & VHR Analysis]
    Q --> S[Taxonomic & Genomic Architecture Panels]
    Q --> T[Ecological Statistics: PERMANOVA, PERMDISP, NMDS & PCoA]
    Q --> U[AMG Functional Profiles & Bubble Heatmaps]

    %% Aplicación de clases
    class F,H,CheckV2 filter;
    class H_vOTU,L,AMGs output;
    class Q rscript;
``` 

## 1. Viral Discovery & Quality Control
Run the consensus discovery tools on your assembled contigs. Merge the outputs and evaluate viral quality using CheckV.
The master R scripts apply a strict filter to retain only Medium-quality, High-quality, and Complete viral genomes that are >= 5,000 bp (5 kb), adhering to MIUViG standards to prevent false positives from small fragments.

## 2. MAG Quality & Taxonomy
Evaluate MAGs using CheckM2. The pipeline automatically filters MAGs retaining those with >=50% completeness and <=10% contamination. Ensure taxonomy is assigned using GTDB-Tk (the R scripts are pre-configured to parse KBase GTDB CSV outputs).

## 3. Abundance Calculation (CoverM)
Generate a unified reference catalog containing both validated MAGs and vOTUs. Map clean reads back to this catalog using CoverM to calculate Trimmed Mean (TPM) values for the entire community simultaneously.

Crucial Step: Ensure MAG contigs are prepended with their MAG ID and the ___ separator (e.g., MAG_001___NODE_1) before mapping so the R script can accurately aggregate contig abundances into total host abundance.

## 4. Host-Virus Interactions
Run the three pillars of interaction prediction:

CRISPR (MinCED): Extract spacers from MAGs and map them against the viral catalog.

Provirus (BLAST): Identify viral sequences integrated directly into host contigs.

Machine Learning (iPHoP): Predict host-virus pairs based on genomic signatures.

## 5. Visualization & Statistics (Master R Scripts)
The repository contains standalone R scripts that merge all the above outputs, handle naming inconsistencies (automatically translating contigs to MAGs using a dictionary), and generate publication-ready plots:

Network_Master.R: Builds the purified Host-Virus interaction network (nodes scaled by TPM, edges colored by detection tool, hosts colored by taxonomic order).

Viral_Catalog_Viz.R: Generates CheckV quality profiles, predicted lifestyles (Lytic vs. Lysogenic), and non-proportional Sankey diagrams tracing viral taxonomy (Class to Order), including the assignment of Unclassified viruses to "Viral Dark Matter".

Ecological_Stats.R: Computes alpha diversity, beta diversity (NMDS & PERMANOVA), and Indicator Species Analysis (IndVal) with corresponding heatmaps and volcano plots.

# Notes on Contamination Filtering
The pipeline includes a dual-filtering mechanism using external blacklists (contaminacion.xlsx and contaminacion_sylph.xlsx) to rigorously remove known physical contaminants from the network based on their GTDB order.

Maintained by the EliLab team.
