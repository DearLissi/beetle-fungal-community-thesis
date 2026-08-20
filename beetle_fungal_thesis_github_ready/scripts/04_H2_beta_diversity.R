library(readxl)
library(tidyverse)
library(vegan)
library(cowplot)

data_dir <- Sys.getenv("BEETLE_DATA_DIR", unset = "data")

# By default, results are saved with the source data. For a different results
# folder, set the BETA_RESULTS_DIR environment variable before running.
results_dir <- Sys.getenv("BETA_RESULTS_DIR", unset = data_dir)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 1. Read and filter the data
# -------------------------------------------------------------------------

metadata_csv <- read_csv(file.path(data_dir, "Sco_Pla_FG_Borneo_metadata.csv"))
metadata_xlsx <- read_excel(file.path(data_dir, "beetles id.xlsx"))

# OTUs are rows and samples are columns at this stage.
OTU_counts <- read_tsv(file.path(data_dir, "OTUsotu_table.tsv")) %>%
  slice(-1) %>%
  column_to_rownames(var = "OTU_ID")

# Retain OTUs with at least 10 reads in total and occurring in >= 3 samples.
OTU_counts <- OTU_counts[rowSums(OTU_counts) >= 10, , drop = FALSE]
OTU_counts <- OTU_counts[rowSums(OTU_counts > 0) >= 3, , drop = FALSE]

# Samples become rows and OTUs become columns.
OTU_counts <- OTU_counts %>%
  t() %>%
  as.data.frame()

# Do not apply a Hellinger transformation for presence/absence Jaccard.

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "#N/A", "0", "NA")] <- NA_character_
  x
}

metadata_csv_clean <- metadata_csv %>%
  transmute(
    index = clean_chr(fungal_metabarcode_ID),
    country = clean_chr(country),
    principal_id = clean_chr(principal_id)
  )

metadata_xlsx_clean <- metadata_xlsx %>%
  transmute(
    index = clean_chr(fungal_metabarcode_ID),
    mt_id = clean_chr(mt_id),
    project_asv_id = clean_chr(project_asv_id),
    subfamily = clean_chr(subfamily),
    BeetleTribe = clean_chr(tribe),
    BeetleSubtribe = clean_chr(subtribe),
    BeetleGenus = clean_chr(genus),
    BeetleSpecies = clean_chr(species)
  )

sample_data_model <- metadata_csv_clean %>%
  left_join(metadata_xlsx_clean, by = "index") %>%
  filter(
    index %in% rownames(OTU_counts),
    !is.na(country),
    !is.na(BeetleGenus),
    !is.na(BeetleTribe)
  ) %>%
  arrange(index) %>%
  mutate(
    country = factor(country),
    BeetleGenus = factor(BeetleGenus),
    BeetleTribe = factor(BeetleTribe)
  )

# Stop early if the metadata join created repeated sample IDs.
if (anyDuplicated(sample_data_model$index)) {
  stop("Duplicated sample IDs were found after joining the metadata.")
}

OTU_counts_model <- OTU_counts[
  sample_data_model$index,
  ,
  drop = FALSE
]

stopifnot(all(sample_data_model$index == rownames(OTU_counts_model)))

# Empty samples do not have a meaningful Jaccard composition.
keep_nonempty <- rowSums(OTU_counts_model) > 0
OTU_counts_model <- OTU_counts_model[keep_nonempty, , drop = FALSE]
sample_data_model <- sample_data_model[keep_nonempty, , drop = FALSE]

# Remove OTUs that are absent from all retained samples.
OTU_counts_model <- OTU_counts_model[
  ,
  colSums(OTU_counts_model) > 0,
  drop = FALSE
]

stopifnot(all(sample_data_model$index == rownames(OTU_counts_model)))

# Optional checks of the final sample sizes.
sample_data_model %>% count(country)
sample_data_model %>% count(BeetleTribe, sort = TRUE)
sample_data_model %>% count(BeetleGenus, sort = TRUE)

# -------------------------------------------------------------------------
# 2. Presence/absence Jaccard distance and PERMANOVA
# -------------------------------------------------------------------------

# binary = TRUE converts every positive abundance to presence (1).
jaccard_dist <- vegdist(
  OTU_counts_model,
  method = "jaccard",
  binary = TRUE
)

set.seed(123)

# Genus is nested within tribe. A sequential nested model gives tribe its own
# test, followed by a test of differences among genera within tribes.
adonis_jaccard <- adonis2(
  jaccard_dist ~ country + BeetleTribe / BeetleGenus,
  data = sample_data_model,
  permutations = 999,
  by = "terms"
)

print(adonis_jaccard)

# -------------------------------------------------------------------------
# 3. PERMDISP for each grouping variable
# -------------------------------------------------------------------------

run_permdisp <- function(dist_obj, group_vec, permutations = 999) {
  group_vec <- as.factor(group_vec)
  dist_mat <- as.matrix(dist_obj)
  
  if (length(group_vec) != nrow(dist_mat)) {
    stop("The group vector and distance matrix have different sample counts.")
  }
  
  keep <- !is.na(group_vec)
  group_counts <- table(group_vec[keep])
  valid_groups <- names(group_counts[group_counts >= 2])
  keep <- keep & group_vec %in% valid_groups
  
  group_vec2 <- droplevels(group_vec[keep])
  dist_obj2 <- as.dist(dist_mat[keep, keep, drop = FALSE])
  
  if (nlevels(group_vec2) < 2) {
    stop("Fewer than two valid groups remain after filtering.")
  }
  
  # Spatial median is preferred for the permutation test.
  bd <- betadisper(
    dist_obj2,
    group_vec2,
    type = "median",
    bias.adjust = TRUE
  )
  
  pt <- permutest(
    bd,
    permutations = permutations
  )
  
  list(
    betadisper = bd,
    permutest = pt,
    group_counts = table(group_vec2)
  )
}

set.seed(123)

permdisp_jaccard_country <- run_permdisp(
  jaccard_dist,
  sample_data_model$country
)

permdisp_jaccard_genus <- run_permdisp(
  jaccard_dist,
  sample_data_model$BeetleGenus
)

permdisp_jaccard_tribe <- run_permdisp(
  jaccard_dist,
  sample_data_model$BeetleTribe
)

print(permdisp_jaccard_country$permutest)
print(permdisp_jaccard_genus$permutest)
print(permdisp_jaccard_tribe$permutest)

# Save all statistical results in one text file.
capture.output(
  {
    cat("Presence/absence Jaccard PERMANOVA\n")
    print(adonis_jaccard)
    
    cat("\nJaccard PERMDISP: country\n")
    print(permdisp_jaccard_country$permutest)
    
    cat("\nJaccard PERMDISP: genus\n")
    print(permdisp_jaccard_genus$permutest)
    
    cat("\nJaccard PERMDISP: tribe\n")
    print(permdisp_jaccard_tribe$permutest)
  },
  file = file.path(results_dir, "beta_diversity_Jaccard_results.txt")
)
# =========================================================================
# 4. JACCARD PCoA
# =========================================================================


# -------------------------------------------------------------------------
# 4.1 Colour palettes
#
# These are the same fixed colours used in the fruit-tree figures.
# Defining them here also allows this PCoA script to run independently.
# -------------------------------------------------------------------------

tribe_cols <- c(
  
  "Corthylini" =
    "#9BC2F1",
  
  "Cryphalini" =
    "#F7C283",
  
  "Dryocoetini" =
    "#F4E085",
  
  "Platypodini" =
    "#B4E197",
  
  "Scolytoplatypodini" =
    "#8CD3C5",
  
  "Tesserocerini" =
    "#F28E8E",
  
  "Trypophloeini" =
    "#C7B2E3",
  
  "Xyleborini" =
    "#F3B2C8",
  
  "Xyloterini" =
    "#D9B8A2",
  
  "Unknown" =
    "white"
)


region_cols <- c(
  
  "French Guiana" =
    "#F9A882",
  
  "Borneo" =
    "#87C0E8",
  
  "No linked fungal data" =
    "white"
)


cat(
  "\nTribe colours used for PCoA:\n"
)

print(
  tribe_cols
)


cat(
  "\nRegion colours used for PCoA:\n"
)


print(
  region_cols
)



# -------------------------------------------------------------------------
# 4.2 Create slightly darker colours for thin point borders
#
# IMPORTANT:
# This does NOT change the original Tribe colours.
#
# tribe_cols = fill
# tribe_border_cols = thin outline only
# -------------------------------------------------------------------------

make_darker <- function(
    colours,
    factor = 0.90
) {
  
  rgb_values <- grDevices::col2rgb(
    colours
  )
  
  darker_values <- round(
    rgb_values * factor
  )
  
  darker_colours <- grDevices::rgb(
    red = darker_values[1, ],
    green = darker_values[2, ],
    blue = darker_values[3, ],
    maxColorValue = 255
  )
  
  names(
    darker_colours
  ) <- names(
    colours
  )
  
  darker_colours
}


region_border_cols <- make_darker(
  region_cols,
  factor = 0.90
)


tribe_border_cols <- make_darker(
  tribe_cols,
  factor = 0.90
)


region_border_cols[
  "No linked fungal data"
] <- "#C7C7C7"


if ("Unknown" %in% names(tribe_border_cols)) {
  
  tribe_border_cols[
    "Unknown"
  ] <- "#C7C7C7"
}



# -------------------------------------------------------------------------
# 4.3 PCoA function
# -------------------------------------------------------------------------

make_pcoa_df <- function(
    dist_obj,
    metadata
) {
  
  # Lingoes correction avoids problems caused by negative eigenvalues.
  
  pco <- wcmdscale(
    dist_obj,
    eig = TRUE,
    add = "lingoes"
  )
  
  
  # Extract first two PCoA axes.
  
  pco_df <- as.data.frame(
    pco$points[
      ,
      1:2,
      drop = FALSE
    ]
  ) %>%
    rownames_to_column(
      "index"
    )
  
  
  colnames(
    pco_df
  )[2:3] <- c(
    "PCoA1",
    "PCoA2"
  )
  
  
  # Percentage variation explained.
  
  positive_eig <- pco$eig[
    pco$eig > 0
  ]
  
  
  axis1_var <- round(
    100 *
      positive_eig[1] /
      sum(positive_eig),
    1
  )
  
  
  axis2_var <- round(
    100 *
      positive_eig[2] /
      sum(positive_eig),
    1
  )
  
  
  # Join metadata.
  
  pco_df <- pco_df %>%
    left_join(
      metadata,
      by = "index"
    )
  
  
  list(
    data = pco_df,
    axis1_var = axis1_var,
    axis2_var = axis2_var,
    model = pco
  )
}



# -------------------------------------------------------------------------
# 4.4 Run Jaccard PCoA
#
# Uses exactly the same Jaccard distance matrix as PERMANOVA.
# -------------------------------------------------------------------------

pcoa_jaccard <- make_pcoa_df(
  jaccard_dist,
  sample_data_model
)



# -------------------------------------------------------------------------
# 4.5 Prepare plotting metadata
#
# This only changes labels used for plotting.
# It does NOT recalculate the PCoA.
# -------------------------------------------------------------------------

pcoa_jaccard$data <- pcoa_jaccard$data %>%
  mutate(
    
    Region = case_when(
      
      as.character(country) %in% c(
        "Malaysia",
        "Malayisia",
        "Borneo"
      ) ~
        "Borneo",
      
      as.character(country) ==
        "French Guiana" ~
        "French Guiana",
      
      is.na(country) ~
        "No linked fungal data",
      
      TRUE ~
        as.character(country)
    ),
    
    
    Region = factor(
      Region,
      levels = names(
        region_cols
      )
    ),
    
    
    BeetleTribe = case_when(
      
      is.na(BeetleTribe) ~
        "Unknown",
      
      trimws(
        as.character(BeetleTribe)
      ) == "" ~
        "Unknown",
      
      TRUE ~
        as.character(BeetleTribe)
    ),
    
    
    # IMPORTANT:
    # Factor order follows the EXISTING fruit-tree tribe_cols object.
    
    BeetleTribe = factor(
      BeetleTribe,
      levels = names(
        tribe_cols
      )
    )
  )



# -------------------------------------------------------------------------
# 4.6 Confirm sample count has not changed
# -------------------------------------------------------------------------

stopifnot(
  nrow(
    pcoa_jaccard$data
  ) ==
    nrow(
      sample_data_model
    )
)



# -------------------------------------------------------------------------
# 4.7 Check colour matching
# -------------------------------------------------------------------------

observed_regions <- unique(
  na.omit(
    as.character(
      pcoa_jaccard$data$Region
    )
  )
)


observed_tribes <- unique(
  na.omit(
    as.character(
      pcoa_jaccard$data$BeetleTribe
    )
  )
)



missing_region_colours <- setdiff(
  observed_regions,
  names(
    region_cols
  )
)


missing_tribe_colours <- setdiff(
  observed_tribes,
  names(
    tribe_cols
  )
)



if (length(missing_region_colours) > 0) {
  
  stop(
    paste0(
      "Missing colours for these Region values: ",
      paste(
        missing_region_colours,
        collapse = ", "
      )
    )
  )
}



if (length(missing_tribe_colours) > 0) {
  
  stop(
    paste0(
      "These tribes are not present in the existing fruit-tree tribe_cols: ",
      paste(
        missing_tribe_colours,
        collapse = ", "
      ),
      ". No new Tribe colours will be generated here."
    )
  )
}



# -------------------------------------------------------------------------
# 4.8 Keep only groups present in this dataset
# -------------------------------------------------------------------------

region_levels_present <- intersect(
  names(
    region_cols
  ),
  observed_regions
)


tribe_levels_present <- intersect(
  names(
    tribe_cols
  ),
  observed_tribes
)



region_cols_present <- region_cols[
  region_levels_present
]


region_border_cols_present <- region_border_cols[
  region_levels_present
]


# IMPORTANT:
# This is only a subset of the EXISTING fruit-tree palette.

tribe_cols_present <- tribe_cols[
  tribe_levels_present
]


tribe_border_cols_present <- tribe_border_cols[
  tribe_levels_present
]



# -------------------------------------------------------------------------
# 4.9 Genus colours
#
# Genus colours remain automatically generated.
#
# This does NOT affect Tribe colours.
# -------------------------------------------------------------------------

genus_levels_present <- pcoa_jaccard$data %>%
  pull(
    BeetleGenus
  ) %>%
  as.character() %>%
  na.omit() %>%
  unique() %>%
  sort()



genus_cols <- scales::hue_pal(
  l = 78,
  c = 65
)(
  length(
    genus_levels_present
  )
)


names(
  genus_cols
) <- genus_levels_present



genus_border_cols <- make_darker(
  genus_cols,
  factor = 0.88
)



pcoa_jaccard$data <- pcoa_jaccard$data %>%
  mutate(
    
    BeetleGenus = factor(
      as.character(
        BeetleGenus
      ),
      levels = genus_levels_present
    )
    
  )



# -------------------------------------------------------------------------
# 4.10 Axis labels
# -------------------------------------------------------------------------

pcoa_x_label <- paste0(
  "PCoA 1 (",
  pcoa_jaccard$axis1_var,
  "% of variation)"
)


pcoa_y_label <- paste0(
  "PCoA 2 (",
  pcoa_jaccard$axis2_var,
  "% of variation)"
)



# -------------------------------------------------------------------------
# 4.11 Shared PCoA theme
#
# Each A / B / C panel keeps its original complete headline.
# -------------------------------------------------------------------------

pcoa_theme <- theme_classic(
  base_size = 13
) +
  theme(
    
    # ---------------------------------------------------------------------
    # Individual panel headline
    # ---------------------------------------------------------------------
    
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5,
      colour = "black",
      margin = margin(
        b = 10
      )
    ),
    
    plot.title.position = "plot",
    
    
    # ---------------------------------------------------------------------
    # Axis
    # ---------------------------------------------------------------------
    
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text = element_text(
      colour = "black",
      size = 10
    ),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.65
    ),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.55
    ),
    
    
    # ---------------------------------------------------------------------
    # Legend
    # ---------------------------------------------------------------------
    
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    
    legend.text = element_text(
      colour = "black",
      size = 9.5
    ),
    
    legend.position = "right",
    
    legend.key.height = grid::unit(
      0.48,
      "cm"
    ),
    
    legend.spacing.y = grid::unit(
      0.04,
      "cm"
    ),
    
    
    # ---------------------------------------------------------------------
    # Background
    # ---------------------------------------------------------------------
    
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    
    # ---------------------------------------------------------------------
    # Margin
    # ---------------------------------------------------------------------
    
    plot.margin = margin(
      t = 12,
      r = 12,
      b = 12,
      l = 12
    )
  )



# -------------------------------------------------------------------------
# 4.12 Region ellipse data
#
# Ellipses require at least three samples in each group.
# -------------------------------------------------------------------------

region_ellipse_data <- pcoa_jaccard$data %>%
  filter(
    !is.na(
      Region
    )
  ) %>%
  group_by(
    Region
  ) %>%
  filter(
    n() >= 3
  ) %>%
  ungroup()



# =========================================================================
# 4.13 PANEL A — REGION
# =========================================================================

p_jaccard_region <- ggplot(
  pcoa_jaccard$data,
  aes(
    x = PCoA1,
    y = PCoA2
  )
) +
  
  
  stat_ellipse(
    
    data = region_ellipse_data,
    
    aes(
      colour = Region,
      group = Region
    ),
    
    type = "norm",
    level = 0.95,
    linewidth = 0.75,
    show.legend = FALSE
  ) +
  
  
  geom_point(
    
    aes(
      fill = Region,
      colour = Region
    ),
    
    shape = 21,
    size = 3,
    stroke = 0.16,
    alpha = 1
  ) +
  
  
  scale_fill_manual(
    
    values =
      region_cols_present,
    
    breaks =
      region_levels_present,
    
    drop = TRUE,
    na.translate = FALSE
  ) +
  
  
  scale_colour_manual(
    
    values =
      region_border_cols_present,
    
    breaks =
      region_levels_present,
    
    drop = TRUE,
    na.translate = FALSE
  ) +
  
  
  coord_equal() +
  
  
  labs(
    
    title =
      "Beetle-associated fungal communities coloured by sampling region",
    
    x = pcoa_x_label,
    y = pcoa_y_label,
    
    fill = "Sampling region",
    colour = "Sampling region"
  ) +
  
  
  pcoa_theme +
  
  
  guides(
    
    colour = "none",
    
    fill = guide_legend(
      
      override.aes = list(
        shape = 21,
        size = 4,
        colour = unname(
          region_border_cols_present
        ),
        stroke = 0.16,
        alpha = 1
      )
      
    )
    
  )



# =========================================================================
# 4.14 PANEL B — BEETLE TRIBE
#
# Uses the EXISTING fruit-tree tribe_cols.
# =========================================================================

p_jaccard_tribe <- ggplot(
  pcoa_jaccard$data,
  aes(
    x = PCoA1,
    y = PCoA2
  )
) +
  
  
  geom_point(
    
    aes(
      fill = BeetleTribe,
      colour = BeetleTribe
    ),
    
    shape = 21,
    size = 3,
    stroke = 0.16,
    alpha = 1
  ) +
  
  
  scale_fill_manual(
    
    values =
      tribe_cols_present,
    
    breaks =
      tribe_levels_present,
    
    drop = TRUE,
    na.translate = FALSE
  ) +
  
  
  scale_colour_manual(
    
    values =
      tribe_border_cols_present,
    
    breaks =
      tribe_levels_present,
    
    drop = TRUE,
    na.translate = FALSE
  ) +
  
  
  coord_equal() +
  
  
  labs(
    
    title =
      "Beetle-associated fungal communities coloured by beetle tribe",
    
    x = pcoa_x_label,
    y = pcoa_y_label,
    
    fill = "Beetle tribe",
    colour = "Beetle tribe"
  ) +
  
  
  pcoa_theme +
  
  
  guides(
    
    colour = "none",
    
    fill = guide_legend(
      
      override.aes = list(
        shape = 21,
        size = 4,
        colour = unname(
          tribe_border_cols_present
        ),
        stroke = 0.16,
        alpha = 1
      )
      
    )
    
  )



# =========================================================================
# 4.15 PANEL C — BEETLE GENUS
# =========================================================================

p_jaccard_genus <- ggplot(
  pcoa_jaccard$data,
  aes(
    x = PCoA1,
    y = PCoA2
  )
) +
  
  
  geom_point(
    
    aes(
      fill = BeetleGenus,
      colour = BeetleGenus
    ),
    
    shape = 21,
    size = 2.8,
    stroke = 0.14,
    alpha = 0.88
  ) +
  
  
  scale_fill_manual(
    
    values =
      genus_cols,
    
    breaks =
      genus_levels_present,
    
    drop = TRUE,
    na.translate = FALSE
  ) +
  
  
  scale_colour_manual(
    
    values =
      genus_border_cols,
    
    breaks =
      genus_levels_present,
    
    drop = TRUE,
    na.translate = FALSE
  ) +
  
  
  coord_equal() +
  
  
  labs(
    
    title =
      "Beetle-associated fungal communities coloured by beetle genus",
    
    x = pcoa_x_label,
    y = pcoa_y_label,
    
    fill = "Beetle genus",
    colour = "Beetle genus"
  ) +
  
  
  pcoa_theme +
  
  
  theme(
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.key.height = grid::unit(
      0.38,
      "cm"
    )
    
  ) +
  
  
  guides(
    
    colour = "none",
    
    fill = guide_legend(
      
      override.aes = list(
        shape = 21,
        size = 3.5,
        colour = "#D0D0D0",
        stroke = 0.14,
        alpha = 1
      )
      
    )
    
  )



# =========================================================================
# 4.16 SHARED FIGURE HEADLINE
#
# REQUIREMENTS:
#
# 1. Main headline centred
# 2. Sub-headline centred
# 3. Sub-headline BLACK
# =========================================================================

pcoa_headline <- ggdraw() +
  
  
  draw_label(
    
    "Community composition of beetle-associated fungi",
    
    x = 0.5,
    y = 0.64,
    
    hjust = 0.5,
    vjust = 0.5,
    
    fontface = "bold",
    size = 17,
    colour = "black"
  ) +
  
  
  draw_label(
    
    "Samples from French Guiana and Borneo analysed using presence–absence Jaccard PCoA",
    
    x = 0.5,
    y = 0.28,
    
    hjust = 0.5,
    vjust = 0.5,
    
    fontface = "plain",
    size = 11.5,
    colour = "black"
  )



# Small gap between the shared headline and the figure body.

pcoa_headline_gap <- ggdraw()



# =========================================================================
# 4.17 DISPLAY INDIVIDUAL PANELS
# =========================================================================

print(
  p_jaccard_region
)


print(
  p_jaccard_tribe
)


print(
  p_jaccard_genus
)



# -------------------------------------------------------------------------
# 4.18 Save PCoA coordinates
# -------------------------------------------------------------------------

write_csv(
  
  pcoa_jaccard$data,
  
  file.path(
    results_dir,
    "Jaccard_PCoA_coordinates.csv"
  )
)



# =========================================================================
# 4.19 INDIVIDUAL FINAL FIGURES
#
# Each figure keeps its own full original headline.
#
# The shared headline is NOT added to these separate files,
# because otherwise the figure would contain two very similar headlines.
# =========================================================================


PCoA_Jaccard_region_FINAL <- p_jaccard_region


PCoA_Jaccard_tribe_FINAL <- p_jaccard_tribe


PCoA_Jaccard_genus_FINAL <- p_jaccard_genus



# -------------------------------------------------------------------------
# Display
# -------------------------------------------------------------------------

print(
  PCoA_Jaccard_region_FINAL
)


print(
  PCoA_Jaccard_tribe_FINAL
)


print(
  PCoA_Jaccard_genus_FINAL
)



# =========================================================================
# 4.20 SAVE INDIVIDUAL FIGURES
# =========================================================================


# -------------------------------------------------------------------------
# Region
# -------------------------------------------------------------------------

ggsave(
  
  filename = file.path(
    results_dir,
    "PCoA_Jaccard_region_FINAL.png"
  ),
  
  plot =
    PCoA_Jaccard_region_FINAL,
  
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


ggsave(
  
  filename = file.path(
    results_dir,
    "PCoA_Jaccard_region_FINAL.pdf"
  ),
  
  plot =
    PCoA_Jaccard_region_FINAL,
  
  width = 8,
  height = 6,
  
  device = cairo_pdf,
  bg = "white"
)



# -------------------------------------------------------------------------
# Tribe
# -------------------------------------------------------------------------

ggsave(
  
  filename = file.path(
    results_dir,
    "PCoA_Jaccard_tribe_FINAL.png"
  ),
  
  plot =
    PCoA_Jaccard_tribe_FINAL,
  
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)


ggsave(
  
  filename = file.path(
    results_dir,
    "PCoA_Jaccard_tribe_FINAL.pdf"
  ),
  
  plot =
    PCoA_Jaccard_tribe_FINAL,
  
  width = 10,
  height = 7,
  
  device = cairo_pdf,
  bg = "white"
)



# -------------------------------------------------------------------------
# Genus
# -------------------------------------------------------------------------

ggsave(
  
  filename = file.path(
    results_dir,
    "PCoA_Jaccard_genus_FINAL.png"
  ),
  
  plot =
    PCoA_Jaccard_genus_FINAL,
  
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)


ggsave(
  
  filename = file.path(
    results_dir,
    "PCoA_Jaccard_genus_FINAL.pdf"
  ),
  
  plot =
    PCoA_Jaccard_genus_FINAL,
  
  width = 12,
  height = 8,
  
  device = cairo_pdf,
  bg = "white"
)



# =========================================================================
# 4.21 COMBINED REGION + TRIBE
#
# Combined figures use short panel headings because the shared
# headline and subtitle already define the biological context.
# =========================================================================

p_jaccard_region_combined <- p_jaccard_region +
  labs(
    title = "Sampling region"
  )


p_jaccard_tribe_combined <- p_jaccard_tribe +
  labs(
    title = "Beetle tribe"
  )


p_jaccard_genus_combined <- p_jaccard_genus +
  labs(
    title = "Beetle genus"
  )

PCoA_Jaccard_region_tribe_body <- plot_grid(
  
  p_jaccard_region_combined,
  p_jaccard_tribe_combined,
  
  labels = c(
    "A",
    "B"
  ),
  
  label_fontface = "bold",
  label_size = 16,
  
  label_x = 0.01,
  label_y = 0.99,
  
  hjust = 0,
  vjust = 1,
  
  ncol = 2,
  
  align = "hv",
  axis = "tblr",
  
  rel_widths = c(
    1,
    1.20
  )
)



PCoA_Jaccard_region_tribe_FINAL <- plot_grid(
  
  pcoa_headline,
  pcoa_headline_gap,
  PCoA_Jaccard_region_tribe_body,
  
  ncol = 1,
  
  rel_heights = c(
    0.11,
    0.005,
    1
  )
)



print(
  PCoA_Jaccard_region_tribe_FINAL
)



ggsave(
  
  filename = file.path(
    results_dir,
    "PCoA_Jaccard_region_and_tribe_FINAL.png"
  ),
  
  plot =
    PCoA_Jaccard_region_tribe_FINAL,
  
  width = 14,
  height = 7.5,
  dpi = 300,
  bg = "white"
)


ggsave(
  
  filename = file.path(
    results_dir,
    "PCoA_Jaccard_region_and_tribe_FINAL.pdf"
  ),
  
  plot =
    PCoA_Jaccard_region_tribe_FINAL,
  
  width = 14,
  height = 7.5,
  
  device = cairo_pdf,
  bg = "white"
)



# =========================================================================
# 4.22 COMBINED REGION + TRIBE + GENUS
#
# MAIN FULL FIGURE
#
# Shared headline:
#   Fungal community composition
#
# Shared sub-headline:
#   Jaccard dissimilarity among beetle-associated fungal communities
#
# A:
#   Fungal community composition by region
#
# B:
#   Fungal community composition by beetle tribe
#
# C:
#   Fungal community composition by beetle genus
# =========================================================================

PCoA_Jaccard_all_body <- plot_grid(
  
  p_jaccard_region_combined,
  p_jaccard_tribe_combined,
  p_jaccard_genus_combined,
  
  labels = c(
    "A",
    "B",
    "C"
  ),
  
  label_fontface = "bold",
  label_size = 16,
  
  label_x = 0.01,
  label_y = 0.99,
  
  hjust = 0,
  vjust = 1,
  
  ncol = 1,
  
  align = "v",
  axis = "lr",
  
  rel_heights = c(
    1,
    1,
    1.05
  )
)



# -------------------------------------------------------------------------
# Complete figure
#
# Header is a full-width object.
# x = 0.5 and hjust = 0.5 above therefore place BOTH headline lines
# exactly at the centre of the entire exported figure.
# -------------------------------------------------------------------------

PCoA_Jaccard_all_FINAL <- plot_grid(
  
  pcoa_headline,
  
  pcoa_headline_gap,
  
  PCoA_Jaccard_all_body,
  
  ncol = 1,
  
  rel_heights = c(
    0.04,
    0.002,
    1
  )
)



# -------------------------------------------------------------------------
# Display main figure
# -------------------------------------------------------------------------

print(
  PCoA_Jaccard_all_FINAL
)



# =========================================================================
# 4.23 SAVE MAIN FULL FIGURE
# =========================================================================

ggsave(
  
  filename = file.path(
    results_dir,
    "PCoA_Jaccard_all_groups_FINAL.png"
  ),
  
  plot =
    PCoA_Jaccard_all_FINAL,
  
  width = 12,
  height = 20,
  dpi = 300,
  bg = "white"
)


ggsave(
  
  filename = file.path(
    results_dir,
    "PCoA_Jaccard_all_groups_FINAL.pdf"
  ),
  
  plot =
    PCoA_Jaccard_all_FINAL,
  
  width = 12,
  height = 20,
  
  device = cairo_pdf,
  bg = "white"
)



# =========================================================================
# 4.24 FINAL CHECKS
# =========================================================================

cat(
  "\n----------------------------------------\n"
)


cat(
  "Final PCoA sample count:",
  nrow(
    pcoa_jaccard$data
  ),
  "\n"
)


cat(
  "Jaccard distance sample count:",
  attr(
    jaccard_dist,
    "Size"
  ),
  "\n"
)


cat(
  "PCoA axis 1:",
  pcoa_jaccard$axis1_var,
  "%\n"
)


cat(
  "PCoA axis 2:",
  pcoa_jaccard$axis2_var,
  "%\n"
)


cat(
  "\nExisting fruit-tree Tribe colours used:\n"
)


print(
  tribe_cols
)


cat(
  "\nTribe colours actually present in this PCoA:\n"
)


print(
  tribe_cols_present
)


cat(
  "\nMain combined figure saved as:\n"
)


cat(
  file.path(
    results_dir,
    "PCoA_Jaccard_all_groups_FINAL.png"
  ),
  "\n"
)


cat(
  "----------------------------------------\n"
)



stopifnot(
  
  nrow(
    pcoa_jaccard$data
  ) ==
    attr(
      jaccard_dist,
      "Size"
    )
)
