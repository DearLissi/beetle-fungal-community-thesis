library(ape)
library(readr)
library(dplyr)
library(readxl)

data_dir <- Sys.getenv("BEETLE_DATA_DIR", unset = "data")
setwd(data_dir)



metabarcode_file <- "beetles id.xlsx"

platy_new <- "platy_dry_MMG_BOLD_COX1_constrained_rooted.treefile"

scoly_new <- "scoly_ent_MMG_BOLD_COX1_constrained_rooted.treefile"


meta <- read_xlsx(metabarcode_file)

platy_new <- read.tree(platy_new)
scoly_new <- read.tree(scoly_new)

meta <- meta %>%
  mutate(mt_id = coalesce(mt_id, project_asv_id))

platy_sub <- meta %>%
  filter(subfamily== "Platypodinae")
scoly_sub <- meta %>%
  filter(subfamily=="Scolytinae")
platy_ids <- platy_sub$mt_id
platy_ids <- platy_ids[!is.na(platy_ids)]
platy_ids <- platy_ids[platy_ids != ""]
platy_ids <- unique(platy_ids)

scoly_ids <- scoly_sub$mt_id
scoly_ids <- scoly_ids[!is.na(scoly_ids)]
scoly_ids <- scoly_ids[scoly_ids != ""]
scoly_ids <- unique(scoly_ids)

cat(
  "Number of unique Platypodinae samples:",
  length(platy_ids),
  "\n"
)

cat(
  "Number of unique Scolytinae samples:",
  length(scoly_ids),
  "\n\n"
)

# -----------------------------
# Function: find tips from tree
# -----------------------------

find_matching_tips <- function(tree, ids) {
  matched_tips <- unique(unlist(lapply(ids, function(id) {
    grep(id, tree$tip.label, value = TRUE, fixed = TRUE)
  })))
  
  return(matched_tips)
}

# -----------------------------
# Function: prune one tree
# -----------------------------

prune_tree <- function(tree, subfamily_name, ids, output_file) {
  

  sample_tips <- find_matching_tips(tree, ids)
  
  cat("Tree:", output_file, "\n")
  cat("Tips matching metabarcode table:", length(sample_tips), "\n")
  
  if (length(sample_tips) == 0) {
    stop(paste("No metabarcode sample tips found in tree for", output_file))
  }
  

  subtree <- keep.tip(tree, sample_tips)
  

  write.tree(subtree, output_file)
  
  cat("Written:", output_file, "\n\n")
  
  return(subtree)
}



platy_new <- prune_tree(
  tree = platy_new,
  subfamily_name = "Platypodinae",
  ids = platy_ids,
  output_file = "platy_new_metabarcode_samples.tre"
)

scoly_new <- prune_tree(
  tree = scoly_new,
  subfamily_name = "Scolytinae",
  ids = scoly_ids,
  output_file = "scoly_new_metabarcode_samples.tre"
)
