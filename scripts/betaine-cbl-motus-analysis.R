# Run from the repository root:
# Rscript scripts/betaine-cbl-motus-analysis.R

################################################################################
## Betaine metabolism and cobalamin biosynthesis analysis in mOTUs-db genomes
################################################################################

# This script analyzes the genomic potential for:
#   1. Cobalamin (Cbl) biosynthesis
#   2. Cbl-dependent glycine betaine (GB) demethylation
#   3. Cbl-independent GB oxidation
#
# The script starts from PyHMMER/Pfam count tables for mOTUs-db representative
# genomes, filters for high-quality Ocean Microbiomics Database (OMDB) genomes,
# calculates order-level prevalence estimates, and exports iTOL-compatible files
# for GTDB R220 tree visualization.

################################################################################
## Setup
################################################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(data.table)
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
  library(tidyr)
})

input_dir <- "input_tables"
gtdb_dir <- "input_GTDB"
output_table_dir <- "output_tables"
output_itol_dir <- "output_itol"

dir.create(output_table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_itol_dir, showWarnings = FALSE, recursive = TRUE)

input_files <- list(
  pfam_counts = file.path(input_dir, "query_pfams_in_motus4.genomes.tsv.gz"),
  motus_reps = file.path(input_dir, "motus4.representatives.genome_export.csv"),
  omdb_genomes = file.path(input_dir, "omdb.genome_export.csv"),
  gene_metadata = file.path(input_dir, "gene_accession_PFAM.xlsx"),
  gtdb_tree = file.path(gtdb_dir, "bac120_r220.tree"),
  gtdb_taxonomy = file.path(gtdb_dir, "bac120_taxonomy_r220.tsv.gz")
)

output_files <- list(
  genome_level_table = file.path(
    output_table_dir,
    "betaine_PFAMs_in_filtered_mOTUs_genomes.tsv"
  ),
  order_level_abundance = file.path(
    output_table_dir,
    "betaine_PFAMs_in_filtered_mOTUs_genomes_abundance.tsv"
  ),
  itol_pruned_tree = file.path(output_itol_dir, "itol_pruned_tree.newick"),
  itol_tree_colors = file.path(output_itol_dir, "itol_tree_colors.txt"),
  itol_heatmap = file.path(output_itol_dir, "itol_heatmap_genes.txt"),
  itol_colorstrip = file.path(output_itol_dir, "itol_colorstrip_groups.txt"),
  itol_barplot = file.path(output_itol_dir, "itol_bar_n_genomes.txt")
)

min_genomes_per_order <- 10

################################################################################
## Helper functions
################################################################################

check_required_columns <- function(data, required_columns, data_name) {
  missing_columns <- setdiff(required_columns, names(data))
  
  if (length(missing_columns) > 0) {
    stop(
      data_name, " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

check_required_files <- function(paths) {
  missing_files <- unlist(paths[!file.exists(unlist(paths))], use.names = FALSE)
  
  if (length(missing_files) > 0) {
    stop(
      "Missing required input file(s):\n",
      paste(missing_files, collapse = "\n"),
      call. = FALSE
    )
  }
}

strip_pfam_version <- function(pfam) {
  if_else(
    is.na(pfam),
    NA_character_,
    sub("\\..*$", "", as.character(pfam))
  )
}

normalize_gtdb_accession <- function(accession) {
  accession %>%
    toupper() %>%
    sub("^(RS_|GB_)", "", x = .) %>%
    sub("\\.\\d+$", "", x = .) %>%
    str_extract("GC[AF]_\\d+")
}

extract_accession_from_tip <- function(tip_label) {
  tip_label %>%
    toupper() %>%
    sub("\\|.*$", "", x = .) %>%
    normalize_gtdb_accession()
}

read_gtdb_taxonomy <- function(path) {
  taxonomy <- data.table::fread(path, header = FALSE, sep = "\t")
  
  has_header <- ncol(taxonomy) >= 2 &&
    tolower(as.character(taxonomy[[1]][1])) %in% c("accession", "user_genome") &&
    tolower(as.character(taxonomy[[2]][1])) %in% c("taxonomy", "classification")
  
  if (ncol(taxonomy) != 2 || has_header) {
    taxonomy <- data.table::fread(path, header = TRUE, sep = "\t")
    taxonomy <- taxonomy[, 1:2]
  }
  
  data.table::setnames(taxonomy, c("accession", "taxonomy"))
  
  taxonomy %>%
    as_tibble() %>%
    separate_wider_delim(
      taxonomy,
      delim = ";",
      names = c(
        "domain", "phylum", "class", "order",
        "family", "genus", "species"
      ),
      too_few = "align_start",
      too_many = "drop"
    ) %>%
    mutate(accession_core = normalize_gtdb_accession(accession))
}

node_label_from_ape_node <- function(tree, node_number) {
  node_index <- as.integer(node_number) - length(tree$tip.label)
  
  if (node_index < 1 || node_index > length(tree$node.label)) {
    return(NA_character_)
  }
  
  tree$node.label[node_index]
}

write_itol_file <- function(lines, path) {
  writeLines(lines, con = path)
  message("Wrote: ", path)
}

assert_columns_exist <- function(data, columns, data_name) {
  missing_columns <- setdiff(columns, names(data))
  
  if (length(missing_columns) > 0) {
    stop(
      data_name, " is missing expected column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

################################################################################
## Input data
################################################################################

check_required_files(input_files)

# PyHMMER results: Pfam hit counts per mOTUs-db genome.
pfam_counts <- read_tsv(input_files$pfam_counts, show_col_types = FALSE)

# Metadata for 124,295 non-redundant representative genomes was downloaded
# from the mOTUs-db web interface:
#
#   https://motus-db.org/genome-cols
#   (also accessible from microbiomics.io landing page)
#
# The total genome collection was filtered by selecting entries with
# "representative" in the "mOTUS4 Status" column. All filtered entries were then
# selected and exported using "Export Metadata".
motus_reps <- read_csv(input_files$motus_reps, show_col_types = FALSE)

# Metadata for marine genomes was downloaded from the OMDB web interface:
#
#   https://omdb.microbiomics.io/repository/ocean/genome-cols
#   (also accessible from microbiomics.io landing page)
#
# All entries were selected and exported using "Export Metadata".
omdb <- read_csv(input_files$omdb_genomes, show_col_types = FALSE)

check_required_columns(pfam_counts, "GENOME", "PyHMMER count table")
check_required_columns(
  motus_reps,
  c(
    "genome", "completeness", "contamination", "domain", "phylum",
    "class", "order", "family", "genus", "species"
  ),
  "mOTUs-db representative genome metadata"
)
check_required_columns(omdb, "Genome", "OMDB genome metadata")

################################################################################
## Filter to high-quality OMDB mOTUs-db representative genomes
################################################################################

motus_reps_filtered <- motus_reps %>%
  mutate(omdb_subset = genome %in% omdb$Genome) %>%
  filter(
    completeness >= 90,
    contamination <= 5,
    omdb_subset
  )

pfam_counts_filtered <- pfam_counts %>%
  filter(GENOME %in% motus_reps_filtered$genome)

################################################################################
## Gene and Pfam metadata
################################################################################

# Manually curated metadata table linking queried Pfam families to gene names.
# This enables downstream grouping of Pfams by putative function.
#
# Notes on gene synonyms:
#   mttB  = mtgB  (trimethylamine methyltransferase or GB--Cbl methyltransferase)
#   cbp   = mtgC  (Cbl-binding protein)
#   metH' = mtgD  (Cbl--tetrahydrofolate methyltransferase)
#   bmt   = mtgE  (GB--homocysteine S-methyltransferase or
#                  Cbl--homocysteine methyltransferase module of split
#                  methionine synthase)
#   gbcA          (GB monooxygenase)
#   gbcB          (GB monooxygenase)

gene_metadata <- read_excel(input_files$gene_metadata) %>%
  mutate(PFAM_stripped = strip_pfam_version(PFAM))

check_required_columns(
  gene_metadata,
  c("Gene name", "PFAM_stripped"),
  "Gene/Pfam metadata"
)

cobalamin_genes <- c(
  "cobA1", "cbiM", "cobV", "cobJ/cbiH", "cobD", "cobN", "cobQ/cbiP",
  "cobB/cbiA", "cbiD", "cbiG", "cobH/cbiA", "cobP/cobU", "cobK/cbiJ",
  "pduO", "cobL/cbiET", "cobI/cbiL"
)

gbc_genes <- c("gbcA", "gbcB")

cobalamin_pfams <- gene_metadata %>%
  filter(`Gene name` %in% cobalamin_genes) %>%
  pull(PFAM_stripped) %>%
  unique() %>%
  na.omit() %>%
  as.character()

gbc_pfams <- gene_metadata %>%
  filter(`Gene name` %in% gbc_genes) %>%
  pull(PFAM_stripped) %>%
  unique() %>%
  na.omit() %>%
  as.character()

single_gene_pfams <- c(
  bmt = "PF02574",
  cbp = "PF02310",
  `metH'` = "PF00809",
  mttB = "PF06253"
)

assert_columns_exist(
  pfam_counts_filtered,
  c(cobalamin_pfams, gbc_pfams, unname(single_gene_pfams)),
  "Filtered PyHMMER count table"
)

################################################################################
## Convert Pfam counts to gene/module presence-absence
################################################################################

pfam_presence <- pfam_counts_filtered %>%
  mutate(across(where(is.numeric), ~ as.integer(.x >= 1)))

pfam_presence_modules <- pfam_presence %>%
  mutate(
    gbcAB = as.integer(rowSums(across(all_of(gbc_pfams), ~ .x == 1)) == length(gbc_pfams)),
    Cob_syn_75_perc = as.integer(
      rowSums(across(all_of(cobalamin_pfams), ~ .x == 1)) >=
        0.75 * length(cobalamin_pfams)
    )
  ) %>%
  rename(
    bmt = all_of(single_gene_pfams[["bmt"]]),
    cbp = all_of(single_gene_pfams[["cbp"]]),
    `metH'` = all_of(single_gene_pfams[["metH'"]]),
    mttB = all_of(single_gene_pfams[["mttB"]])
  ) %>%
  mutate(
    `mttB-cbp-metH'` = as.integer(mttB == 1 & cbp == 1 & `metH'` == 1)
  ) %>%
  select(-starts_with("PF"))

################################################################################
## Add taxonomy and save genome-level table
################################################################################

motus_taxonomy <- motus_reps_filtered %>%
  select(genome, domain, phylum, class, order, family, genus, species) %>%
  rename(GENOME = genome)

pfam_presence_taxonomy <- pfam_presence_modules %>%
  left_join(motus_taxonomy, by = "GENOME") %>%
  filter(domain != "Archaea")

write_tsv(pfam_presence_taxonomy, output_files$genome_level_table)
message("Wrote: ", output_files$genome_level_table)

################################################################################
## Calculate order-level prevalence
################################################################################

gene_module_columns <- c("mttB-cbp-metH'", "gbcAB", "Cob_syn_75_perc")

order_abundance <- pfam_presence_taxonomy %>%
  mutate(
    order = paste0("o__", sub("^o__", "", order)),
    phylum = paste0("p__", sub("^p__", "", phylum)),
    class = paste0("c__", sub("^c__", "", class))
  ) %>%
  group_by(order, phylum, class) %>%
  mutate(n_genomes = n()) %>%
  filter(n_genomes >= min_genomes_per_order) %>%
  summarise(
    n_genomes = first(n_genomes),
    across(all_of(gene_module_columns), ~ mean(.x == 1, na.rm = TRUE) * 100),
    .groups = "drop"
  ) %>%
  arrange(desc(n_genomes))

write_tsv(order_abundance, output_files$order_level_abundance)
message("Wrote: ", output_files$order_level_abundance)

################################################################################
## Export iTOL input files for GTDB R220 order-level tree visualization
################################################################################

# This section links order-level prevalence estimates to the GTDB R220 bacterial
# species tree. It prunes the tree to one representative genome per target order
# and exports iTOL-compatible datasets:
#   1. Pruned Newick tree with internal node labels
#   2. TREE_COLORS dataset for highlighted taxonomic groups
#   3. DATASET_HEATMAP dataset for gene/prevalence rings
#   4. DATASET_COLORSTRIP dataset used as a legend for highlighted clades
#   5. DATASET_SIMPLEBAR dataset showing the number of mOTUs representatives
#      per order
#
# Required GTDB input files:
#   input_GTDB/bac120_r220.tree
#   input_GTDB/bac120_taxonomy_r220.tsv.gz
#
# GTDB R220 download source:
#   https://data.gtdb.ecogenomic.org/releases/release220/220.0/
#
# Note: mOTUs-db genomes were annotated with GTDB R220 using GTDB-Tk v2.4.

highlight_phyla <- c(
  "p__Desulfobacterota",
  "p__Actinomycetota",
  "p__Chloroflexota",
  "p__Acidobacteriota",
  "p__Cyanobacteriota",
  "p__Bacillota_C",
  "p__Bacillota_A"
)

highlight_classes <- c(
  "c__Alphaproteobacteria",
  "c__Gammaproteobacteria"
)

highlight_palette <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#84D9D4"
)

other_color <- "#D0D0D0"

required_abundance_columns <- c(
  "order", "phylum", "class", "n_genomes",
  "mttB-cbp-metH'", "gbcAB", "Cob_syn_75_perc"
)

check_required_columns(
  order_abundance,
  required_abundance_columns,
  "Order-level abundance table"
)

order_abundance <- order_abundance %>%
  filter(
    str_starts(order, "o__"),
    str_starts(phylum, "p__"),
    str_starts(class, "c__")
  ) %>%
  mutate(
    n_genomes = as.numeric(n_genomes),
    `mttB-cbp-metH'` = as.numeric(`mttB-cbp-metH'`),
    gbcAB = as.numeric(gbcAB),
    Cob_syn_75_perc = as.numeric(Cob_syn_75_perc)
  ) %>%
  distinct(order, .keep_all = TRUE)

if (nrow(order_abundance) == 0) {
  stop(
    "No valid order-level rows remained after filtering the abundance table.",
    call. = FALSE
  )
}

gtdb_taxonomy <- read_gtdb_taxonomy(input_files$gtdb_taxonomy)
check_required_columns(
  gtdb_taxonomy,
  c("accession_core", "phylum", "class", "order"),
  "GTDB taxonomy"
)

gtdb_tree <- ape::read.tree(input_files$gtdb_tree)

################################################################################
## Prune GTDB tree to one representative genome per target order
################################################################################

tip_taxonomy <- tibble(
  tip_label = gtdb_tree$tip.label,
  accession_core = extract_accession_from_tip(gtdb_tree$tip.label)
) %>%
  filter(!is.na(accession_core)) %>%
  left_join(
    gtdb_taxonomy %>% select(accession_core, phylum, class, order),
    by = "accession_core"
  ) %>%
  filter(!is.na(order), str_starts(order, "o__"))

representative_tips <- tip_taxonomy %>%
  filter(order %in% order_abundance$order) %>%
  group_by(order) %>%
  slice(1) %>%
  ungroup()

if (nrow(representative_tips) == 0) {
  stop(
    "No tree tips matched the target orders in the abundance table.",
    call. = FALSE
  )
}

itol_tree <- ape::keep.tip(gtdb_tree, representative_tips$tip_label)
itol_tree <- ape::ladderize(itol_tree, right = TRUE)

tip_metadata <- representative_tips %>%
  select(tip_label, order, phylum, class) %>%
  distinct(tip_label, .keep_all = TRUE) %>%
  slice(match(itol_tree$tip.label, tip_label)) %>%
  mutate(
    label = sub("^o__", "", order),
    phylum_clean = sub("^p__", "", phylum),
    class_clean = sub("^c__", "", class)
  )

stopifnot(identical(tip_metadata$tip_label, itol_tree$tip.label))

itol_tree$tip.label <- tip_metadata$label
itol_tree$node.label <- paste0("N", seq_len(itol_tree$Nnode))

tip_metadata <- tip_metadata %>%
  select(label, order, phylum, class, phylum_clean, class_clean) %>%
  mutate(
    color_key = case_when(
      phylum %in% highlight_phyla ~ phylum_clean,
      class %in% highlight_classes ~ class_clean,
      TRUE ~ "Other"
    )
  )

################################################################################
## Prepare iTOL data tables
################################################################################

heatmap_data <- order_abundance %>%
  mutate(order_clean = sub("^o__", "", order)) %>%
  select(order_clean, `mttB-cbp-metH'`, gbcAB, Cob_syn_75_perc) %>%
  filter(order_clean %in% itol_tree$tip.label) %>%
  pivot_longer(
    cols = c(`mttB-cbp-metH'`, gbcAB, Cob_syn_75_perc),
    names_to = "gene",
    values_to = "percent"
  ) %>%
  mutate(
    label = order_clean,
    gene = recode(
      gene,
      `mttB-cbp-metH'` = "mtgBCD",
      gbcAB = "gbcAB",
      Cob_syn_75_perc = "Cbl"
    )
  ) %>%
  select(label, gene, percent)

heatmap_wide <- heatmap_data %>%
  pivot_wider(names_from = gene, values_from = percent) %>%
  right_join(tibble(label = itol_tree$tip.label), by = "label") %>%
  arrange(match(label, itol_tree$tip.label)) %>%
  replace_na(list(mtgBCD = 0, gbcAB = 0, Cbl = 0))

stopifnot(identical(heatmap_wide$label, itol_tree$tip.label))

bar_data <- order_abundance %>%
  mutate(order_clean = sub("^o__", "", order)) %>%
  distinct(order_clean, .keep_all = TRUE) %>%
  select(order_clean, n_genomes) %>%
  right_join(tibble(label = itol_tree$tip.label), by = c("order_clean" = "label")) %>%
  transmute(
    label = order_clean,
    n_genomes = coalesce(as.numeric(n_genomes), 0)
  ) %>%
  arrange(match(label, itol_tree$tip.label))

stopifnot(identical(bar_data$label, itol_tree$tip.label))

highlight_keys <- sort(setdiff(unique(tip_metadata$color_key), "Other"))

if (length(highlight_keys) > length(highlight_palette)) {
  warning("More highlighted groups than palette colors; colors will be recycled.")
}

highlight_colors <- setNames(
  rep(highlight_palette, length.out = length(highlight_keys)),
  highlight_keys
)
highlight_colors <- c(highlight_colors, Other = other_color)

################################################################################
## Export A: pruned Newick tree
################################################################################

ape::write.tree(itol_tree, file = output_files$itol_pruned_tree)
message("Wrote: ", output_files$itol_pruned_tree)

################################################################################
## Export B: iTOL TREE_COLORS dataset
################################################################################

tree_color_lines <- c(
  "TREE_COLORS",
  "SEPARATOR TAB",
  "DATA"
)

for (tip_label in itol_tree$tip.label) {
  color_key <- tip_metadata$color_key[match(tip_label, tip_metadata$label)]
  color_key <- ifelse(is.na(color_key), "Other", color_key)
  
  color_hex <- unname(highlight_colors[[color_key]])
  color_hex <- ifelse(is.na(color_hex), other_color, color_hex)
  
  tree_color_lines <- c(
    tree_color_lines,
    sprintf("%s\tlabel\t%s\tnormal\t1", tip_label, color_hex),
    sprintf("%s\tbranch\t%s\tnormal\t2", tip_label, color_hex)
  )
}

for (color_key in highlight_keys) {
  group_tips <- tip_metadata %>%
    filter(color_key == .env$color_key) %>%
    pull(label) %>%
    intersect(itol_tree$tip.label)
  
  if (length(group_tips) < 2 || !ape::is.monophyletic(itol_tree, group_tips)) {
    next
  }
  
  mrca_node <- ape::getMRCA(itol_tree, group_tips)
  node_label <- node_label_from_ape_node(itol_tree, mrca_node)
  
  if (!is.na(node_label)) {
    tree_color_lines <- c(
      tree_color_lines,
      sprintf("%s\tclade\t%s\tnormal\t2", node_label, highlight_colors[[color_key]])
    )
  }
}

write_itol_file(tree_color_lines, output_files$itol_tree_colors)

################################################################################
## Export C: iTOL DATASET_HEATMAP dataset
################################################################################

heatmap_header <- c(
  "DATASET_HEATMAP",
  "SEPARATOR TAB",
  "DATASET_LABEL\t% prevalence",
  "COLOR\t#555555",
  "FIELD_LABELS\tmtgBCD\tgbcAB\tCbl",
  "FIELD_COLORS\t#440154\t#21918C\t#FDE725",
  "COLOR_MIN\t#440154",
  "COLOR_MID\t#21918C",
  "COLOR_MAX\t#FDE725",
  "LEGEND_TITLE\t% prevalence",
  "DATA"
)

heatmap_lines <- sprintf(
  "%s\t%.6g\t%.6g\t%.6g",
  heatmap_wide$label,
  heatmap_wide$mtgBCD,
  heatmap_wide$gbcAB,
  heatmap_wide$Cbl
)

write_itol_file(
  c(heatmap_header, heatmap_lines),
  output_files$itol_heatmap
)

################################################################################
## Export D: iTOL DATASET_COLORSTRIP legend for highlighted clades
################################################################################

legend_keys <- c(sort(highlight_keys), "Other")
legend_colors <- unname(highlight_colors[legend_keys])

colorstrip_lines <- c(
  "DATASET_COLORSTRIP",
  "SEPARATOR TAB",
  "DATASET_LABEL\tHighlighted clades",
  "COLOR\t#000000",
  "STRIP_WIDTH\t0",
  "SHOW_LABELS\t0",
  "SHOW_STRIP_LABELS\t0",
  "BORDER_WIDTH\t0",
  "MARGIN\t-2",
  "LEGEND_TITLE\tHighlighted clades",
  sprintf("LEGEND_SHAPES\t%s", paste(rep(1, length(legend_keys)), collapse = "\t")),
  sprintf("LEGEND_COLORS\t%s", paste(legend_colors, collapse = "\t")),
  sprintf("LEGEND_LABELS\t%s", paste(legend_keys, collapse = "\t")),
  "DATA",
  sprintf("%s\trgba(0,0,0,0)", itol_tree$tip.label)
)

write_itol_file(
  colorstrip_lines,
  output_files$itol_colorstrip
)

################################################################################
## Export E: iTOL DATASET_SIMPLEBAR dataset for clade sizes
################################################################################

bar_header <- c(
  "DATASET_SIMPLEBAR",
  "SEPARATOR TAB",
  "DATASET_LABEL\tClade sizes",
  "COLOR\t#555555",
  "WIDTH\t80",
  "MARGIN\t2",
  "SHOW_VALUES\t1",
  "LEGEND_TITLE\tClade sizes",
  "LEGEND_SHAPES\t1",
  "LEGEND_COLORS\t#555555",
  "LEGEND_LABELS\t# mOTUs Reps. per Order",
  "DATA"
)

bar_lines <- sprintf("%s\t%.6g", bar_data$label, bar_data$n_genomes)

write_itol_file(
  c(bar_header, bar_lines),
  output_files$itol_barplot
)

message("Done. Output tables are available in: ", output_table_dir)
message("Done. iTOL files are available in: ", output_itol_dir)

sessionInfo()
