#Alpha diversity multivariate analysis
library(tidyverse)
library(vegan)
library(cowplot)
library(nlme)
library(lme4)
library(flexplot)
library(ggplot2)
library(readxl)
library(ape)
library(phyr)

data_dir <- Sys.getenv("BEETLE_DATA_DIR", unset = "data")
setwd(data_dir)

scoly_tree <- read.tree("scoly_new_metabarcode_samples.tre")
platy_tree <- read.tree("platy_new_metabarcode_samples.tre")
metadata_xlsx <- read_excel("beetles id.xlsx")

#read in data with OTUs and their presence in each sample
OTU_table <- read_tsv('OTUsotu_table.tsv') %>% 
  slice(-1) %>%
  column_to_rownames(var = "OTU_ID") %>% 
  t() %>%
  as.data.frame()

#visualise sequencing depth
OTU_table %>%
  mutate(total = rowSums(.,)) %>%
  ggplot(aes(x = total)) + 
  geom_histogram(binwidth = 500) +
  scale_y_continuous(labels = scales::comma)

#see where bulk of data sits
OTU_table %>%
  mutate(total = rowSums(.,)) %>%
  ggplot(aes(x=1,y=total)) +
  geom_jitter() +
  scale_y_log10()

#sort by the total number of sequences
OTU_table %>%
  mutate(total = rowSums(.,)) %>%
  arrange(total) %>%
  select(total) %>%
  head(50)

#plot rarefaction
rarecurve(OTU_table, step=1000, label=FALSE)
#zoom in
rarecurve(OTU_table, step=1000, label=FALSE, xlim = c(0,10000))

#rarefy to 1000 seqs
rarefied_richness <- OTU_table %>%
  rarefy(., 1000) %>%
  as_tibble(rownames ='name') %>%
  select(name, rarefied_richness = value)

#non-rarefied richness
alpha_stats <- tibble(
  name = rownames(OTU_table), 
  richness = rowSums(OTU_table > 0)
  ) %>%
  inner_join(rarefied_richness, by = 'name')

#read phylogenetic diversity indices
MPD_MNTD <- read_csv('MPD_MNTD_results.csv')
#merge the richness and phylogenetic diversity
alpha_stats <- left_join(alpha_stats, 
                         MPD_MNTD, 
                         by = c('name' = 'rowname'))

#read in sample data
sample_data <- read_csv('Sco_Pla_FG_Borneo_metadata.csv') %>%
  #remove samples not in dissimilarity matrix
  filter(fungal_metabarcode_ID %in% alpha_stats$name)



alpha_stats_metadata_base <- alpha_stats %>%
  inner_join(
    sample_data,
    by = c("name" = "fungal_metabarcode_ID")
  )

# ============================================================
# Match beetle IDs using beetles id.xlsx as the authority
# ============================================================


clean_id <- function(x) {
  
  x <- as.character(x)
  x <- stringr::str_trim(x)
  
  # 删除末尾备注星号
  x <- stringr::str_remove(
    x,
    "\\*+$"
  )
  
  x <- stringr::str_trim(x)
  x <- toupper(x)
  
  x[
    x %in% c(
      "",
      "0",
      "NA",
      "N/A",
      "#N/A"
    )
  ] <- NA_character_
  
  x
}

clean_taxonomy <- function(x) {
  
  x <- as.character(x)
  x <- stringr::str_squish(x)
  
  x[
    is.na(x) |
      tolower(x) %in% c(
        "",
        "na",
        "n/a",
        "#n/a",
        "unknown",
        "unidentified",
        "not identified",
        "none"
      )
  ] <- NA_character_
  
  x
}

first_non_missing <- function(x) {
  
  x <- as.character(x)
  x <- x[
    !is.na(x) &
      stringr::str_trim(x) != ""
  ]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  x[1]
}

safe_mean <- function(x) {
  
  # 避免数值列被意外读成字符
  x <- suppressWarnings(
    as.numeric(x)
  )
  
  # 整组都没有有效值时返回 NA，而不是 NaN
  if (
    length(x) == 0L ||
    all(is.na(x))
  ) {
    return(NA_real_)
  }
  
  mean(
    x,
    na.rm = TRUE
  )
}
scoly_tree$tip.label <- clean_id(
  scoly_tree$tip.label
)

platy_tree$tip.label <- clean_id(
  platy_tree$tip.label
)

stopifnot(
  anyDuplicated(scoly_tree$tip.label) == 0,
  anyDuplicated(platy_tree$tip.label) == 0
)

# ============================================================
# 1. IDs actually present in the two pruned trees
# ============================================================

all_tree_ids <- union(
  clean_id(scoly_tree$tip.label),
  clean_id(platy_tree$tip.label)
)

# ============================================================
# 2. Build the Excel-to-tree crosswalk
# ============================================================

beetle_crosswalk <- metadata_xlsx %>%
  transmute(
    xlsx_mt_id =
      clean_id(mt_id),
    
    xlsx_asv_id =
      clean_id(project_asv_id),
    
    TaxSubfamily =
      clean_taxonomy(subfamily),
    
    TaxTribe =
      clean_taxonomy(tribe),
    
    TaxGenus =
      clean_taxonomy(genus),
    
    TaxSpecies =
      clean_taxonomy(species)
  ) %>%
  
  mutate(
    # This is exactly the rule used when pruning the trees
    beetle = coalesce(
      xlsx_mt_id,
      xlsx_asv_id
    ),
    
    # project_asv_id connects to CSV asv_id whenever available
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
  
  group_by(
    beetle,
    join_type,
    join_id
  ) %>%
  
  summarise(
    TaxSubfamily =
      first_non_missing(TaxSubfamily),
    
    TaxTribe =
      first_non_missing(TaxTribe),
    
    TaxGenus =
      first_non_missing(TaxGenus),
    
    TaxSpecies =
      first_non_missing(TaxSpecies),
    
    .groups = "drop"
  )

crosswalk_conflicts <- beetle_crosswalk %>%
  distinct(
    join_type,
    join_id,
    beetle
  ) %>%
  
  count(
    join_type,
    join_id,
    name = "n_tree_tips"
  ) %>%
  
  filter(
    n_tree_tips > 1
  )


print(
  crosswalk_conflicts,
  n = Inf
)


if (nrow(crosswalk_conflicts) > 0) {
  stop(
    "One Excel lookup ID maps to more than one tree tip."
  )
}

sample_metadata_clean <- alpha_stats_metadata_base %>%
  mutate(
    sample_row_id = row_number(),
    
    csv_mt_id =
      clean_id(mt_id),
    
    csv_asv_id =
      clean_id(asv_id),
    
    Country = case_when(
      tolower(
        trimws(as.character(Country))
      ) %in% c(
        "malaysia",
        "borneo"
      ) ~
        "Borneo",
      
      tolower(
        trimws(as.character(Country))
      ) == "french guiana" ~
        "French Guiana",
      
      TRUE ~
        trimws(as.character(Country))
    )
  )

matches_by_asv <- sample_metadata_clean %>%
  inner_join(
    beetle_crosswalk %>%
      filter(
        join_type == "ASV"
      ) %>%
      select(
        join_id,
        beetle,
        TaxSubfamily,
        TaxTribe,
        TaxGenus,
        TaxSpecies
      ),
    
    by = c(
      "csv_asv_id" = "join_id"
    )
  ) %>%
  
  mutate(
    match_rule =
      "xlsx project_asv_id -> CSV asv_id"
  )

matches_by_mt <- sample_metadata_clean %>%
  inner_join(
    beetle_crosswalk %>%
      filter(
        join_type == "MT"
      ) %>%
      select(
        join_id,
        beetle,
        TaxSubfamily,
        TaxTribe,
        TaxGenus,
        TaxSpecies
      ),
    
    by = c(
      "csv_mt_id" = "join_id"
    )
  ) %>%
  
  mutate(
    match_rule =
      "xlsx mt_id -> CSV mt_id"
  )

all_match_candidates <- bind_rows(
  matches_by_asv,
  matches_by_mt
) %>%
  
  distinct(
    sample_row_id,
    beetle,
    .keep_all = TRUE
  )

ambiguous_sample_matches <- all_match_candidates %>%
  distinct(
    sample_row_id,
    name,
    csv_mt_id,
    csv_asv_id,
    beetle,
    match_rule
  ) %>%
  
  group_by(
    sample_row_id
  ) %>%
  
  summarise(
    name = first(name),
    csv_mt_id = first(csv_mt_id),
    csv_asv_id = first(csv_asv_id),
    
    n_tree_tips =
      n_distinct(beetle),
    
    candidate_tree_tips =
      paste(
        sort(unique(beetle)),
        collapse = " | "
      ),
    
    candidate_rules =
      paste(
        sort(unique(match_rule)),
        collapse = " | "
      ),
    
    .groups = "drop"
  ) %>%
  
  filter(
    n_tree_tips > 1
  )


print(
  ambiguous_sample_matches,
  n = Inf,
  width = Inf
)


if (nrow(ambiguous_sample_matches) > 0) {
  stop(
    "Some fungal samples map to more than one tree tip."
  )
}

alpha_stats_metadata <- all_match_candidates %>%
  arrange(
    sample_row_id
  ) %>%
  
  distinct(
    sample_row_id,
    .keep_all = TRUE
  ) %>%
  
  select(
    -sample_row_id
  )

alpha_phylo_new <- alpha_stats_metadata %>%
  filter(
    !is.na(beetle),
    !is.na(Country),
    Country %in% c(
      "French Guiana",
      "Borneo"
    )
  ) %>%
  
  group_by(
    beetle,
    Country
  ) %>%
  
  summarise(
    Subfamily =
      first_non_missing(TaxSubfamily),
    
    Tribe =
      first_non_missing(TaxTribe),
    
    Genus =
      first_non_missing(TaxGenus),
    
    Species =
      first_non_missing(TaxSpecies),
    
    rarefied_richness =
      safe_mean(rarefied_richness),
    
    mpd.obs =
      safe_mean(mpd.obs),
    
    mntd.obs =
      safe_mean(mntd.obs),
    
    .groups = "drop"
  )

beetle_country_conflicts <- alpha_phylo_new %>%
  distinct(
    beetle,
    Country
  ) %>%
  
  count(
    beetle,
    name = "n_countries"
  ) %>%
  
  filter(
    n_countries > 1
  )


print(
  beetle_country_conflicts,
  n = Inf
)

tree_matching_summary <- tibble(
  Tree = c(
    "Scolytinae",
    "Platypodinae"
  ),
  
  Total_tree_tips = c(
    length(scoly_tree$tip.label),
    length(platy_tree$tip.label)
  ),
  
  Tips_in_alpha_phylo = c(
    sum(
      clean_id(scoly_tree$tip.label) %in%
        alpha_phylo_new$beetle
    ),
    
    sum(
      clean_id(platy_tree$tip.label) %in%
        alpha_phylo_new$beetle
    )
  )
) %>%
  
  mutate(
    Missing_tips =
      Total_tree_tips -
      Tips_in_alpha_phylo
  )


print(
  tree_matching_summary,
  n = Inf
)

saveRDS(
  alpha_phylo_new,
  "alpha_phylo_new_validated_2026-08-04.rds"
)

saveRDS(
  alpha_stats_metadata,
  "alpha_stats_metadata_validated_2026-08-04.rds"
)

saveRDS(
  beetle_crosswalk,
  "beetle_crosswalk_validated_2026-08-04.rds"
)

saveRDS(
  scoly_tree,
  "scoly_tree_validated_2026-08-04.rds"
)

saveRDS(
  platy_tree,
  "platy_tree_validated_2026-08-04.rds"
)

write_csv(
  beetle_crosswalk,
  "beetle_crosswalk_validated_2026-08-04.csv"
)

write_csv(
  alpha_phylo_new,
  "alpha_phylo_new_validated_2026-08-04.csv"
)

write_csv(
  tree_matching_summary,
  "tree_matching_summary_2026-08-04.csv"
)

unmatched_tree_tips <- tibble(
  beetle = c(
    "6MSL1388",
    "6MSL2104",
    "6MSL2701",
    "6MSL5685"
  ),
  Tree = "Platypodinae",
  status = "No matching asv_id in current fungal metadata CSV"
)

write_csv(
  unmatched_tree_tips,
  "unmatched_platypodinae_tips_2026-08-04.csv"
)

model_sample_size <- bind_rows(
  tibble(
    Tree = "Scolytinae",
    response = "Rarefied richness",
    n = sum(
      alpha_phylo_new$beetle %in% scoly_tree$tip.label &
        is.finite(alpha_phylo_new$rarefied_richness)
    )
  ),
  
  tibble(
    Tree = "Scolytinae",
    response = "MPD",
    n = sum(
      alpha_phylo_new$beetle %in% scoly_tree$tip.label &
        is.finite(alpha_phylo_new$mpd.obs)
    )
  ),
  
  tibble(
    Tree = "Scolytinae",
    response = "MNTD",
    n = sum(
      alpha_phylo_new$beetle %in% scoly_tree$tip.label &
        is.finite(alpha_phylo_new$mntd.obs)
    )
  ),
  
  tibble(
    Tree = "Platypodinae",
    response = "Rarefied richness",
    n = sum(
      alpha_phylo_new$beetle %in% platy_tree$tip.label &
        is.finite(alpha_phylo_new$rarefied_richness)
    )
  ),
  
  tibble(
    Tree = "Platypodinae",
    response = "MPD",
    n = sum(
      alpha_phylo_new$beetle %in% platy_tree$tip.label &
        is.finite(alpha_phylo_new$mpd.obs)
    )
  ),
  
  tibble(
    Tree = "Platypodinae",
    response = "MNTD",
    n = sum(
      alpha_phylo_new$beetle %in% platy_tree$tip.label &
        is.finite(alpha_phylo_new$mntd.obs)
    )
  )
)

print(
  model_sample_size,
  n = Inf
)

prepare_pglmm_data <- function(
    tree,
    data,
    response
) {
  
  dat <- data %>%
    filter(
      beetle %in% tree$tip.label,
      Country %in% c(
        "French Guiana",
        "Borneo"
      ),
      is.finite(.data[[response]])
    ) %>%
    mutate(
      Country = factor(
        Country,
        levels = c(
          "French Guiana",
          "Borneo"
        )
      )
    ) %>%
    distinct(
      beetle,
      .keep_all = TRUE
    )
  
  
  shared_ids <- intersect(
    tree$tip.label,
    dat$beetle
  )
  
  
  tree_use <- ape::keep.tip(
    tree,
    shared_ids
  )
  
  
  dat <- dat %>%
    filter(
      beetle %in% tree_use$tip.label
    ) %>%
    arrange(
      match(
        beetle,
        tree_use$tip.label
      )
    )
  
  
  rownames(dat) <- dat$beetle
  
  
  stopifnot(
    identical(
      rownames(dat),
      tree_use$tip.label
    )
  )
  
  
  list(
    data = dat,
    tree = tree_use,
    response = response
  )
}

scoly_richness_input <- prepare_pglmm_data(
  tree = scoly_tree,
  data = alpha_phylo_new,
  response = "rarefied_richness"
)

scoly_mpd_input <- prepare_pglmm_data(
  tree = scoly_tree,
  data = alpha_phylo_new,
  response = "mpd.obs"
)

scoly_mntd_input <- prepare_pglmm_data(
  tree = scoly_tree,
  data = alpha_phylo_new,
  response = "mntd.obs"
)


platy_richness_input <- prepare_pglmm_data(
  tree = platy_tree,
  data = alpha_phylo_new,
  response = "rarefied_richness"
)

platy_mpd_input <- prepare_pglmm_data(
  tree = platy_tree,
  data = alpha_phylo_new,
  response = "mpd.obs"
)

platy_mntd_input <- prepare_pglmm_data(
  tree = platy_tree,
  data = alpha_phylo_new,
  response = "mntd.obs"
)

c(
  scoly_richness = nrow(scoly_richness_input$data),
  scoly_MPD = nrow(scoly_mpd_input$data),
  scoly_MNTD = nrow(scoly_mntd_input$data),
  platy_richness = nrow(platy_richness_input$data),
  platy_MPD = nrow(platy_mpd_input$data),
  platy_MNTD = nrow(platy_mntd_input$data)
)

scoly_richness_model <- phyr::pglmm_compare(
  rarefied_richness ~ Country,
  data = scoly_richness_input$data,
  phy = scoly_richness_input$tree,
  family = "gaussian",
  REML = TRUE
)


scoly_mpd_model <- phyr::pglmm_compare(
  mpd.obs ~ Country,
  data = scoly_mpd_input$data,
  phy = scoly_mpd_input$tree,
  family = "gaussian",
  REML = TRUE
)


scoly_mntd_model <- phyr::pglmm_compare(
  mntd.obs ~ Country,
  data = scoly_mntd_input$data,
  phy = scoly_mntd_input$tree,
  family = "gaussian",
  REML = TRUE
)


platy_richness_model <- phyr::pglmm_compare(
  rarefied_richness ~ Country,
  data = platy_richness_input$data,
  phy = platy_richness_input$tree,
  family = "gaussian",
  REML = TRUE
)


platy_mpd_model <- phyr::pglmm_compare(
  mpd.obs ~ Country,
  data = platy_mpd_input$data,
  phy = platy_mpd_input$tree,
  family = "gaussian",
  REML = TRUE
)


platy_mntd_model <- phyr::pglmm_compare(
  mntd.obs ~ Country,
  data = platy_mntd_input$data,
  phy = platy_mntd_input$tree,
  family = "gaussian",
  REML = TRUE
)

pglmm_models_validated <- list(
  Scolytinae_richness = scoly_richness_model,
  Scolytinae_MPD = scoly_mpd_model,
  Scolytinae_MNTD = scoly_mntd_model,
  
  Platypodinae_richness = platy_richness_model,
  Platypodinae_MPD = platy_mpd_model,
  Platypodinae_MNTD = platy_mntd_model
)

saveRDS(
  pglmm_models_validated,
  "pglmm_models_validated_2026-08-04.rds"
)

capture.output(
  {
    cat("Validated PGLMM models\n")
    cat("======================\n")
    cat("Generated:", format(Sys.time()), "\n\n")
    
    for (model_name in names(pglmm_models_validated)) {
      
      cat("\n")
      cat("============================================================\n")
      cat("MODEL:", model_name, "\n")
      cat("============================================================\n\n")
      
      current_model <- pglmm_models_validated[[model_name]]
      
      cat("MODEL OUTPUT\n")
      cat("------------\n")
      print(current_model)
      
      cat("\nMODEL SUMMARY\n")
      cat("-------------\n")
      print(summary(current_model))
      
      cat("\n\n")
    }
  },
  
  file = "pglmm_models_validated_2026-08-04.txt"
)

