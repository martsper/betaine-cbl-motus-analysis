# Betaine metabolism and cobalamin biosynthesis analysis in mOTUs-db genomes

This repository contains the R workflow used to analyze genomic potential for cobalamin biosynthesis, cobalamin-dependent glycine betaine demethylation, and cobalamin-independent glycine betaine oxidation in high-quality Ocean Microbiomics Database genomes represented in mOTUs-db.

The analysis starts from PyHMMER/Pfam count tables, filters for high-quality OMDB genomes, calculates order-level prevalence estimates, and exports iTOL-compatible files for visualization on the GTDB R220 bacterial tree.

## Repository structure

```text
.
├── scripts/
│   └── betaine_cbl_analysis.R
├── input_tables/
├── input_GTDB/
├── output_tables/
├── output_itol/
├── environment.yml
└── README.md
```

## Requirements

The analysis can be reproduced using the provided conda environment:

```bash
conda env create -f environment.yml
conda activate betaine-cbl-analysis
```

## Input files

The script expects the following input files:

```text
input_tables/query_pfams_in_motus4.genomes.tsv.gz
input_tables/motus4.representatives.genome_export.csv
input_tables/omdb.genome_export.csv
input_tables/gene_accession_PFAM.xlsx
input_GTDB/bac120_r220.tree
input_GTDB/bac120_taxonomy_r220.tsv.gz
```

The mOTUs-db representative genome metadata was downloaded from:

```text
https://motus-db.org/genome-cols
```

The OMDB genome metadata was downloaded from:

```text
https://omdb.microbiomics.io/repository/ocean/genome-cols
```

GTDB R220 files were obtained from:

```text
https://data.gtdb.ecogenomic.org/releases/release220/220.0/
```

## Running the analysis

Run the workflow from the repository root:

```bash
Rscript scripts/betaine_cbl_analysis.R
```

## Outputs

The script writes genome-level and order-level summary tables to:

```text
output_tables/
```

It also writes iTOL-compatible files to:

```text
output_itol/
```

Main output files:

```text
output_tables/betaine_PFAMs_in_filtered_mOTUs_genomes.tsv
output_tables/betaine_PFAMs_in_filtered_mOTUs_genomes_abundance.tsv
output_itol/itol_pruned_tree.newick
output_itol/itol_tree_colors.txt
output_itol/itol_heatmap_genes.txt
output_itol/itol_colorstrip_groups.txt
output_itol/itol_bar_n_genomes.txt
```

## Notes on input data

The mOTUs-db representative genome metadata was exported from the mOTUs-db web interface after filtering for genomes with representative status.

The OMDB genome metadata was exported from the Ocean Microbiomics Database web interface.

The GTDB tree and taxonomy files correspond to GTDB release R220.

## Citation

If using this workflow, please cite the associated manuscript.

