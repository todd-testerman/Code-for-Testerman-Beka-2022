# A Large-Scale, Multi-Year Microbial Community Survey of a Freshwater Trout Aquaculture Facility

[![DOI](https://img.shields.io/badge/DOI-10.1093%2Ffemsec%2Ffiac101-blue)](https://doi.org/10.1093/femsec/fiac101)

Code used to process, analyze, and visualize microbiome data from a trout aquaculture facility. This repository contains the complete bioinformatics and statistical analysis pipeline for 16S rRNA gene sequencing data published in [FEMS Microbiology Ecology](https://doi.org/10.1093/femsec/fiac101).

## Project Overview

This project examines bacterial community composition in trout aquaculture facilities to understand:
- How microbial communities differ between sample types (wall swabs, water filters, etc.)
- Alpha and beta diversity patterns across sampling locations
- Relationships between microbiome composition and fish disease status (*Flavobacterium columnare*)
- Temporal and spatial variation in bacterial communities across multiple hatchery sites

## Repository Structure

```
.
├── README.md
├── scripts/
│   └── Qiime2_script_for_USDA_paper.sh   # Bioinformatics pipeline (QIIME2)
├── analysis/
│   └── USDA_Analysis_rmarkdown_Github.md # Statistical analysis (R Markdown)
└── figures/
    └── *.png                              # Generated visualizations
```

## Analysis Workflow

### 1. Bioinformatics Processing (QIIME2)

The shell script (`scripts/Qiime2_script_for_USDA_paper.sh`) performs:
- Import of demultiplexed FASTQ files
- Quality filtering and denoising using DADA2
- Sequence alignment with MAFFT
- Phylogenetic tree construction
- Alpha and beta diversity calculations
- Taxonomic classification using GreenGenes and Silva databases

### 2. Statistical Analysis (R)

The R Markdown document (`analysis/USDA_Analysis_rmarkdown_Github.md`) includes:
- **Data Processing**: Import QIIME2 artifacts, taxonomic filtering
- **Alpha Diversity**: Faith's Phylogenetic Diversity, species richness
- **Beta Diversity**: Bray-Curtis, Jaccard, Weighted/Unweighted UniFrac ordinations
- **Statistical Testing**: PERMANOVA, pairwise comparisons, normality tests
- **Differential Abundance**: ANCOM-BC analysis for compositional data
- **Visualizations**: Barplots, heatmaps, ordination plots, Venn diagrams

## Requirements

### QIIME2 Pipeline
- [QIIME2](https://qiime2.org/) (version 2019.10 or compatible)
- DADA2 plugin
- GreenGenes and/or Silva reference databases

### R Analysis
Key R packages required:
- `phyloseq` - Microbiome data handling
- `qiime2R` - Import QIIME2 artifacts
- `vegan` - Diversity statistics and PERMANOVA
- `ggplot2`, `patchwork` - Visualization
- `ANCOMBC` - Compositional analysis
- `microViz`, `microbiome` - Microbiome-specific tools
- `DESeq2`, `corncob` - Differential abundance

## Data Summary

- **Total Samples**: 134 samples analyzed
- **ASVs (Amplicon Sequence Variants)**: 25,790 initial; 12,554 after filtering
- **Sample Types**: Wall swabs, Baffle, Sterivex filters, Tailscreen, Prefilter
- **Time Period**: 2017-2019 (3 years of data)
- **Rarefaction Depth**: 10,000 reads per sample

## Usage

### Running the QIIME2 Pipeline

```bash
# Activate QIIME2 environment
conda activate qiime2-2019.10

# Run the bioinformatics pipeline
bash scripts/Qiime2_script_for_USDA_paper.sh
```

**Note**: Modify file paths in the script to match your data location.

### Running the R Analysis

Open the R Markdown file in RStudio and execute the code chunks sequentially. Ensure the required QIIME2 output artifacts are available:
- `filtered-table.qza` - Feature table
- `taxonomy.qza` - Taxonomic assignments
- `rooted-tree.qza` - Phylogenetic tree
- `Mapping_File_filtered_FINAL.txt` - Sample metadata

## Key Findings

The analysis addresses several research questions:
1. Bacterial community composition differences between surface and water samples
2. Impact of fish maturity stage on microbiome structure
3. Association between microbial communities and *F. columnare* disease status
4. Comparison of communities across three aquaculture sites

## Authors

- **Todd Testerman** - University of Connecticut, Department of Molecular and Cell Biology
- **Lidia Beka** - University of Connecticut, Department of Molecular and Cell Biology
- **Stephen R. Reichley** - Riverence Provisions LLC
- **Stacy King** - Riverence Provisions LLC
- **Timothy J. Welch** - USDA/ARS National Center for Cool and Cold Water Aquaculture
- **Gregory D. Wiens** - USDA/ARS National Center for Cool and Cold Water Aquaculture
- **Joerg Graf** - University of Connecticut, Department of Molecular and Cell Biology

## Citation

If you use this code or data, please cite:

> Testerman T, Beka L, Reichley SR, King S, Welch TJ, Wiens GD, Graf J. **A large-scale, multi-year microbial community survey of a freshwater trout aquaculture facility.** *FEMS Microbiology Ecology*. 2022;98(10):fiac101. doi: [10.1093/femsec/fiac101](https://doi.org/10.1093/femsec/fiac101)

### BibTeX

```bibtex
@article{testerman2022large,
  title={A large-scale, multi-year microbial community survey of a freshwater trout aquaculture facility},
  author={Testerman, Todd and Beka, Lidia and Reichley, Stephen R and King, Stacy and Welch, Timothy J and Wiens, Gregory D and Graf, Joerg},
  journal={FEMS Microbiology Ecology},
  volume={98},
  number={10},
  pages={fiac101},
  year={2022},
  publisher={Oxford University Press},
  doi={10.1093/femsec/fiac101}
}
```

## License

This project is shared for academic and research purposes.
