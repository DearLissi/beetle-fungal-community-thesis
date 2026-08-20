# ============================================================
# H3 — HOST PHYLOGENETIC DISTANCE AND FUNGAL-COMMUNITY DISSIMILARITY
#
# This script reconstructs the final H3 analytical workflow.
#
# Transparency note:
# The original working H3 script was not retained. This version was
# reconstructed from the documented workflow and validated against
# the saved final H3 analysis object. Validation reproduced the same:
#   - retained host tips (135 Scolytinae; 52 Platypodinae)
#   - fungal OTU sets (2,463; 1,630)
#   - host patristic-distance matrices
#   - fungal binary-Jaccard distance matrices
#   - Mantel statistics and permutation P values
#
# Input data are not included in this repository because they were
# supplied by the research group.
# ============================================================

library(tidyverse)
library(readxl)
library(ape)
library(vegan)
library(dendextend)

data_dir <- Sys.getenv("BEETLE_DATA_DIR", unset = "data")
results_dir <- Sys.getenv("BEETLE_RESULTS_DIR", unset = "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

clean_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_remove(x, "\\*+$")
  x <- stringr::str_trim(x)
  x <- toupper(x)
  x[x %in% c("", "0", "NA", "N/A", "#N/A")] <- NA_character_
  x
}

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "#N/A", "0", "NA")] <- NA_character_
  x
}

clean_taxonomy <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x[
    is.na(x) |
      tolower(x) %in% c(
        "", "na", "n/a", "#n/a",
        "unknown", "unidentified",
        "not identified", "none"
      )
  ] <- NA_character_
  x
}

first_non_missing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & stringr::str_trim(x) != ""]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

# ------------------------------------------------------------
# 1. Host phylogenies
# ------------------------------------------------------------

scoly_tree <- read.tree(
  file.path(data_dir, "scoly_new_metabarcode_samples.tre")
)

platy_tree <- read.tree(
  file.path(data_dir, "platy_new_metabarcode_samples.tre")
)

scoly_tree$tip.label <- clean_id(scoly_tree$tip.label)
platy_tree$tip.label <- clean_id(platy_tree$tip.label)

stopifnot(
  anyDuplicated(scoly_tree$tip.label) == 0,
  anyDuplicated(platy_tree$tip.label) == 0
)

all_tree_ids <- union(
  scoly_tree$tip.label,
  platy_tree$tip.label
)

# ------------------------------------------------------------
# 2. Conservative specimen crosswalk
#
# Matching follows the same logic used in the validated H1 workflow:
# project_asv_id -> metadata asv_id where available; otherwise
# mt_id -> metadata mt_id.
# ------------------------------------------------------------

metadata_csv <- read_csv(
  file.path(data_dir, "Sco_Pla_FG_Borneo_metadata.csv"),
  show_col_types = FALSE
)

metadata_xlsx <- read_excel(
  file.path(data_dir, "beetles id.xlsx")
)

beetle_crosswalk <- metadata_xlsx %>%
  transmute(
    xlsx_mt_id = clean_id(mt_id),
    xlsx_asv_id = clean_id(project_asv_id),
    TaxSubfamily = clean_taxonomy(subfamily),
    TaxTribe = clean_taxonomy(tribe),
    TaxGenus = clean_taxonomy(genus),
    TaxSpecies = clean_taxonomy(species)
  ) %>%
  mutate(
    beetle = coalesce(xlsx_mt_id, xlsx_asv_id),
    join_type = case_when(
      !is.na(xlsx_asv_id) ~ "ASV",
      !is.na(xlsx_mt_id) ~ "MT",
      TRUE ~ NA_character_
    ),
    join_id = case_when(
      join_type == "ASV" ~ xlsx_asv_id,
      join_type == "MT" ~ xlsx_mt_id,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    beetle %in% all_tree_ids,
    !is.na(join_id)
  ) %>%
  group_by(beetle, join_type, join_id) %>%
  summarise(
    TaxSubfamily = first_non_missing(TaxSubfamily),
    TaxTribe = first_non_missing(TaxTribe),
    TaxGenus = first_non_missing(TaxGenus),
    TaxSpecies = first_non_missing(TaxSpecies),
    .groups = "drop"
  )

crosswalk_conflicts <- beetle_crosswalk %>%
  distinct(join_type, join_id, beetle) %>%
  count(join_type, join_id, name = "n_tree_tips") %>%
  filter(n_tree_tips > 1)

if (nrow(crosswalk_conflicts) > 0) {
  print(crosswalk_conflicts, n = Inf)
  stop("One lookup ID maps to more than one tree tip.")
}

sample_metadata_clean <- metadata_csv %>%
  mutate(
    sample_row_id = row_number(),
    fungal_sample = clean_chr(fungal_metabarcode_ID),
    csv_mt_id = clean_id(mt_id),
    csv_asv_id = clean_id(asv_id),
    Country = case_when(
      tolower(trimws(as.character(country))) %in%
        c("malaysia", "borneo") ~ "Borneo",
      tolower(trimws(as.character(country))) ==
        "french guiana" ~ "French Guiana",
      TRUE ~ trimws(as.character(country))
    )
  )

matches_by_asv <- sample_metadata_clean %>%
  inner_join(
    beetle_crosswalk %>%
      filter(join_type == "ASV") %>%
      select(
        join_id, beetle,
        TaxSubfamily, TaxTribe, TaxGenus, TaxSpecies
      ),
    by = c("csv_asv_id" = "join_id")
  ) %>%
  mutate(match_rule = "xlsx project_asv_id -> CSV asv_id")

matches_by_mt <- sample_metadata_clean %>%
  inner_join(
    beetle_crosswalk %>%
      filter(join_type == "MT") %>%
      select(
        join_id, beetle,
        TaxSubfamily, TaxTribe, TaxGenus, TaxSpecies
      ),
    by = c("csv_mt_id" = "join_id")
  ) %>%
  mutate(match_rule = "xlsx mt_id -> CSV mt_id")

all_match_candidates <- bind_rows(
  matches_by_asv,
  matches_by_mt
) %>%
  filter(!is.na(fungal_sample)) %>%
  distinct(sample_row_id, beetle, .keep_all = TRUE)

ambiguous_sample_matches <- all_match_candidates %>%
  distinct(sample_row_id, fungal_sample, beetle, match_rule) %>%
  group_by(sample_row_id) %>%
  summarise(
    n_tree_tips = n_distinct(beetle),
    candidate_tree_tips =
      paste(sort(unique(beetle)), collapse = " | "),
    .groups = "drop"
  ) %>%
  filter(n_tree_tips > 1)

if (nrow(ambiguous_sample_matches) > 0) {
  print(ambiguous_sample_matches, n = Inf)
  stop("Some fungal samples map to more than one tree tip.")
}

sample_to_beetle <- all_match_candidates %>%
  arrange(sample_row_id) %>%
  distinct(sample_row_id, .keep_all = TRUE)

# ------------------------------------------------------------
# 3. Fungal OTU filtering
# ------------------------------------------------------------

OTU_counts <- read_tsv(
  file.path(data_dir, "OTUsotu_table.tsv"),
  show_col_types = FALSE
) %>%
  slice(-1) %>%
  column_to_rownames(var = "OTU_ID")

OTU_counts <- OTU_counts[
  rowSums(OTU_counts) >= 10,
  ,
  drop = FALSE
]

OTU_counts <- OTU_counts[
  rowSums(OTU_counts > 0) >= 3,
  ,
  drop = FALSE
]

OTU_by_sample <- OTU_counts %>%
  t() %>%
  as.data.frame()

sample_to_beetle <- sample_to_beetle %>%
  filter(fungal_sample %in% rownames(OTU_by_sample))

if (anyDuplicated(sample_to_beetle$fungal_sample)) {
  stop(
    "A fungal metabarcoding sample appears more than once after matching."
  )
}

otu_samples <- OTU_by_sample[
  sample_to_beetle$fungal_sample,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(otu_samples),
    sample_to_beetle$fungal_sample
  )
)

otu_samples$beetle <- sample_to_beetle$beetle

# ------------------------------------------------------------
# 4. Pool replicate samples by beetle-tree tip
# ------------------------------------------------------------

otu_by_beetle <- otu_samples %>%
  as_tibble() %>%
  group_by(beetle) %>%
  summarise(
    across(everything(), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

beetle_ids <- otu_by_beetle$beetle

otu_by_beetle_matrix <- otu_by_beetle %>%
  select(-beetle) %>%
  as.data.frame()

rownames(otu_by_beetle_matrix) <- beetle_ids

# ------------------------------------------------------------
# 5. H3 analysis
# ------------------------------------------------------------

run_H3 <- function(host_tree, subfamily_name) {

  valid_beetles <- intersect(
    host_tree$tip.label,
    rownames(otu_by_beetle_matrix)
  )

  host_tree_sub <- keep.tip(
    host_tree,
    valid_beetles
  )

  fungal_counts_sub <- otu_by_beetle_matrix[
    host_tree_sub$tip.label,
    ,
    drop = FALSE
  ]

  keep_nonempty_hosts <- rowSums(fungal_counts_sub) > 0

  if (any(!keep_nonempty_hosts)) {
    host_tree_sub <- drop.tip(
      host_tree_sub,
      rownames(fungal_counts_sub)[!keep_nonempty_hosts]
    )

    fungal_counts_sub <- fungal_counts_sub[
      keep_nonempty_hosts,
      ,
      drop = FALSE
    ]
  }

  fungal_counts_sub <- fungal_counts_sub[
    ,
    colSums(fungal_counts_sub) > 0,
    drop = FALSE
  ]

  fungal_counts_sub <- fungal_counts_sub[
    host_tree_sub$tip.label,
    ,
    drop = FALSE
  ]

  otu_sub <- as.data.frame(
    (fungal_counts_sub > 0) * 1L
  )

  fungal_dist <- vegan::vegdist(
    otu_sub,
    method = "jaccard",
    binary = TRUE
  )

  beetle_phylo_matrix <- ape::cophenetic.phylo(
    host_tree_sub
  )

  beetle_phylo_matrix <- beetle_phylo_matrix[
    host_tree_sub$tip.label,
    host_tree_sub$tip.label,
    drop = FALSE
  ]

  beetle_dist <- as.dist(beetle_phylo_matrix)

  set.seed(123)

  mantel_result <- vegan::mantel(
    xdis = beetle_dist,
    ydis = fungal_dist,
    method = "spearman",
    permutations = 9999
  )

  # Average-linkage clustering of the host-distance and fungal
  # community-distance matrices.
  beetle_dendrogram <- as.dendrogram(
    hclust(beetle_dist, method = "average")
  )

  fungal_dendrogram <- as.dendrogram(
    hclust(fungal_dist, method = "average")
  )

  dend_pair <- dendextend::dendlist(
    beetle_dendrogram,
    fungal_dendrogram
  )

  dend_pair_plot <- dendextend::untangle(
    dend_pair,
    method = "step2side"
  )

  beetle_dendrogram_plot <- dend_pair_plot[[1]]
  fungal_dendrogram_plot <- dend_pair_plot[[2]]

  cophenetic_correlation <- dendextend::cor_cophenetic(
    beetle_dendrogram,
    fungal_dendrogram
  )

  bakers_gamma <- dendextend::cor_bakers_gamma(
    beetle_dendrogram,
    fungal_dendrogram
  )

  entanglement_value <- dendextend::entanglement(
    dend_pair_plot
  )

  summary_row <- tibble(
    Subfamily = subfamily_name,
    Host_tips = length(host_tree_sub$tip.label),
    Fungal_OTUs = ncol(otu_sub),
    Mantel_r = unname(mantel_result$statistic),
    Mantel_P = mantel_result$signif,
    Cophenetic_correlation = cophenetic_correlation,
    Bakers_gamma = bakers_gamma,
    Entanglement = entanglement_value
  )

  list(
    tree = host_tree_sub,
    otu = otu_sub,
    beetle_distance = beetle_dist,
    fungal_distance = fungal_dist,
    beetle_dendrogram = beetle_dendrogram,
    fungal_dendrogram = fungal_dendrogram,
    beetle_dendrogram_plot = beetle_dendrogram_plot,
    fungal_dendrogram_plot = fungal_dendrogram_plot,
    mantel = mantel_result,
    cophenetic_correlation = cophenetic_correlation,
    bakers_gamma = bakers_gamma,
    entanglement = entanglement_value,
    summary = summary_row
  )
}

scoly_H3 <- run_H3(
  scoly_tree,
  "Scolytinae"
)

platy_H3 <- run_H3(
  platy_tree,
  "Platypodinae"
)

# ------------------------------------------------------------
# 6. Host metadata used by the final tanglegram plotting script
# ------------------------------------------------------------

final_beetles <- union(
  scoly_H3$tree$tip.label,
  platy_H3$tree$tip.label
)

beetle_metadata <- sample_to_beetle %>%
  filter(beetle %in% final_beetles) %>%
  group_by(beetle) %>%
  summarise(
    Subfamily = first_non_missing(TaxSubfamily),
    Tribe = first_non_missing(TaxTribe),
    Genus = first_non_missing(TaxGenus),
    Species = first_non_missing(TaxSpecies),
    Country = first_non_missing(Country),
    n_fungal_samples = n_distinct(fungal_sample),
    .groups = "drop"
  )

final_tanglegram_analysis <- list(
  Scolytinae = scoly_H3,
  Platypodinae = platy_H3,
  beetle_metadata = beetle_metadata
)

saveRDS(
  final_tanglegram_analysis,
  file.path(results_dir, "H3_analysis_objects.rds")
)

validation_summary <- bind_rows(
  scoly_H3$summary,
  platy_H3$summary
)

write_csv(
  validation_summary,
  file.path(results_dir, "H3_summary.csv")
)

print(validation_summary)
