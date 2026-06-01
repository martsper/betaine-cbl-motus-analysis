# Betaine degradation and cobalamin biosynthesis analysis in mOTUs-db genomes

This repository contains the R workflow used to analyze genomic potential for cobalamin-dependent glycine betaine demethylation, cobalamin biosynthesis, and cobalamin-independent glycine betaine oxidation in mOTUs-db genomes.

The analysis starts from PyHMMER/Pfam count tables, filters for high-quality genomes (>= 90% completeness, <= 5% contamination), removes genomes that are not part of the Ocean Microbiomics Database (OMDB), calculates order-level prevalences, and exports iTOL-compatible files for visualization on the GTDB R220 bacterial tree.

## Repository structure

```text
.
├── scripts/
│   └── betaine-cbl-motus-analysis.R
├── input_tables/
├── input_GTDB/
├── output_tables/
├── output_itol/
├── environment_short.yml
├── environment_long.yml
└── README.md
```

## Requirements

The analysis can be reproduced using the provided conda environment:

```bash
conda env create -f environment_short.yml
conda activate betaine-cbl-motus-analysis
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
(also accessible from microbiomics.io landing page)
```

The OMDB genome metadata was downloaded from:

```text
https://omdb.microbiomics.io/repository/ocean/genome-cols
(also accessible from microbiomics.io landing page)
```

GTDB R220 files were obtained from:

```text
https://data.gtdb.ecogenomic.org/releases/release220/220.0/
```

## Running the analysis

Run the workflow from the repository root:

```bash
Rscript scripts/betaine-cbl-motus-analysis.R
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

## References

The analysis is part of the following pre-print:

> The glycine betaine-cobalamin feedback loop drives cross-feeding between marine bacteria and algae<br>
Jonathan Hammer, Myriel Staack, Martin Sperfeld, Tom Haufschild, Nicolai Kallscheuer, Delia A. Narváez-Barragán, Carl-Eric Wegner, Kirsten Küsel, Shinichi Sunagawa, Georg Pohnert, Torsten Schubert, Einat Segev, Christian Jogler.<br>
bioRxiv 2025.12.19.695462; [https://doi.org/10.64898/2025.12.19.695462](https://doi.org/10.64898/2025.12.19.695462)

The mOTUs-db genomes were described here:

> The mOTUs online database provides web-accessible genomic context to taxonomic profiling of microbial communities<br>
Marija Dmitrijeva, Hans-Joachim Ruscheweyh, Lilith Feer, Kang Li, Samuel Miravet-Verde, Anna Sintsova, Daniel R Mende, Georg Zeller, Shinichi Sunagawa.<br>
Nucleic Acids Research, Volume 53, Issue D1, 6 January 2025, Pages D797–D805, [https://doi.org/10.1093/nar/gkae1004](https://doi.org/10.1093/nar/gkae1004)


