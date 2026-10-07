# Viral Metagenomics & Host-Virus Interaction Pipeline 🦠🧬

A comprehensive and highly robust pipeline for viral discovery, quality control, abundance quantification, and host-virus interaction network reconstruction from metagenomic assemblies. Developed for **EliLab**.

This workflow integrates state-of-the-art tools to reliably identify viral contigs, pair them with their bacterial/archaeal hosts using both empirical (CRISPR/BLAST) and machine-learning (iPHoP) evidence, and visualize the complex ecological dynamics of the community.

## 📊 Pipeline Workflow

The pipeline consists of modular steps, going from raw assemblies to publication-ready visualizations.

```mermaid
graph TD
    A[Metagenomic Assemblies] --> B(Viral Discovery)
    A --> C(MAG Binning)
    
    B -->|geNomad, VirSorter2, VIBRANT| D[Putative Viral Contigs]
    C -->|metaBAT2, VAMB, etc.| E[Raw MAGs]
    
    D --> F{CheckV Quality Filter}
    F -->|>= 5kb & Med/High/Complete| G[Validated vOTUs]
    
    E --> H{CheckM2 Filter}
    H -->|Completeness >=50, Contamination <=10| I[High-Quality MAGs]
    I --> J(Taxonomy: GTDB-Tk via KBase)
    
    G --> K(Abundance)
    I --> K
    K -->|CoverM TPM| L[Abundance Matrix]
    
    G --> M(Interactions)
    I --> M
    M -->|MinCED| N[CRISPR Spacers]
    M -->|BLAST| O[Provirus Integration]
    M -->|iPHoP| P[ML Host Prediction]
    
    L --> Q((Master R Script))
    J --> Q
    N --> Q
    O --> Q
    P --> Q
    
    Q --> R[Host-Virus Network Graphs]
    Q --> S[Taxonomic Sankey Diagrams]
    Q --> T[Ecological Statistics NMDS/IndVal]
