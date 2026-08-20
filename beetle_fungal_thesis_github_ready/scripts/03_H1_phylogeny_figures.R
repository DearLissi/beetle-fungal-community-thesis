# NOTE:
# This script uses analysis objects created by 02_H1_alpha_diversity.R.
# Run scripts 02 and 03 in the same R session.
#
# ============================================================
# 1. Packages
# ============================================================

library(ape)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(scales)
library(grid)


# ============================================================
# 2. Plotting helper functions
# ============================================================

clean_plot_id <- function(x) {
  
  x <- as.character(x)
  x <- stringr::str_trim(x)
  
  # Remove trailing annotation stars
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


clean_plot_text <- function(x) {
  
  x <- as.character(x)
  x <- stringr::str_squish(x)
  
  x[
    toupper(x) %in% c(
      "",
      "NA",
      "N/A",
      "#N/A"
    )
  ] <- NA_character_
  
  x
}


first_non_missing_plot <- function(x) {
  
  x <- clean_plot_text(x)
  x <- x[!is.na(x)]
  
  if (length(x) == 0L) {
    return(NA_character_)
  }
  
  x[1]
}


get_metric_limits <- function(x) {
  
  x <- suppressWarnings(
    as.numeric(x)
  )
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    
    stop(
      "No finite values found for one of the plotting metrics."
    )
  }
  
  limits <- range(
    x,
    na.rm = TRUE
  )
  
  if (limits[1] == limits[2]) {
    
    limits <- limits + c(
      -0.5,
      0.5
    )
  }
  
  limits
}


# ============================================================
# 3. Prepare fungal diversity values
# ============================================================

alpha_plot_data <- alpha_phylo_new %>%
  
  transmute(
    
    tip_id =
      clean_plot_id(beetle),
    
    Country_plot = case_when(
      
      Country == "French Guiana" ~
        "French Guiana",
      
      Country %in% c(
        "Borneo",
        "Malaysia"
      ) ~
        "Borneo",
      
      TRUE ~
        NA_character_
    ),
    
    Richness =
      suppressWarnings(
        as.numeric(
          rarefied_richness
        )
      ),
    
    MPD =
      suppressWarnings(
        as.numeric(
          mpd.obs
        )
      ),
    
    MNTD =
      suppressWarnings(
        as.numeric(
          mntd.obs
        )
      )
  ) %>%
  
  filter(
    !is.na(tip_id)
  )


# ------------------------------------------------------------
# There must be one row per beetle tip
# ------------------------------------------------------------

duplicated_alpha_tips <- alpha_plot_data %>%
  
  count(
    tip_id,
    name = "n_rows"
  ) %>%
  
  filter(
    n_rows > 1
  )


if (nrow(duplicated_alpha_tips) > 0) {
  
  print(
    duplicated_alpha_tips,
    n = Inf
  )
  
  stop(
    "Duplicated beetle tips are present in alpha_plot_data."
  )
}


# ============================================================
# 4. Prepare taxonomy from alpha_phylo_new
# ============================================================

taxonomy_from_alpha <- alpha_phylo_new %>%
  
  transmute(
    
    tip_id =
      clean_plot_id(beetle),
    
    Subfamily =
      clean_plot_text(Subfamily),
    
    Tribe =
      clean_plot_text(Tribe),
    
    Genus =
      clean_plot_text(Genus),
    
    Species =
      clean_plot_text(Species)
  ) %>%
  
  filter(
    !is.na(tip_id)
  ) %>%
  
  group_by(
    tip_id
  ) %>%
  
  summarise(
    
    Subfamily =
      first_non_missing_plot(Subfamily),
    
    Tribe =
      first_non_missing_plot(Tribe),
    
    Genus =
      first_non_missing_plot(Genus),
    
    Species =
      first_non_missing_plot(Species),
    
    .groups = "drop"
  )


# ============================================================
# 5. Prepare taxonomy from beetles id.xlsx
#
# Tree-tip ID follows:
# coalesce(mt_id, project_asv_id)
# ============================================================

taxonomy_from_xlsx <- metadata_xlsx %>%
  
  transmute(
    
    xlsx_mt_id =
      clean_plot_id(mt_id),
    
    xlsx_asv_id =
      clean_plot_id(project_asv_id),
    
    tip_id =
      coalesce(
        xlsx_mt_id,
        xlsx_asv_id
      ),
    
    Subfamily =
      clean_plot_text(subfamily),
    
    Tribe =
      clean_plot_text(tribe),
    
    Genus =
      clean_plot_text(genus),
    
    Species =
      clean_plot_text(species)
  ) %>%
  
  filter(
    !is.na(tip_id)
  ) %>%
  
  group_by(
    tip_id
  ) %>%
  
  summarise(
    
    Subfamily =
      first_non_missing_plot(Subfamily),
    
    Tribe =
      first_non_missing_plot(Tribe),
    
    Genus =
      first_non_missing_plot(Genus),
    
    Species =
      first_non_missing_plot(Species),
    
    .groups = "drop"
  )


# ============================================================
# 6. Combine taxonomy sources
#
# Prefer taxonomy already present in alpha_phylo_new.
# Use Excel taxonomy for tree tips without fungal data.
# ============================================================

taxonomy_lookup <- full_join(
  
  taxonomy_from_xlsx,
  taxonomy_from_alpha,
  
  by = "tip_id",
  
  suffix = c(
    "_xlsx",
    "_alpha"
  )
  
) %>%
  
  transmute(
    
    tip_id,
    
    Subfamily =
      coalesce(
        Subfamily_alpha,
        Subfamily_xlsx
      ),
    
    Tribe =
      coalesce(
        Tribe_alpha,
        Tribe_xlsx
      ),
    
    Genus =
      coalesce(
        Genus_alpha,
        Genus_xlsx
      ),
    
    Species =
      coalesce(
        Species_alpha,
        Species_xlsx
      )
  )


# ============================================================
# 7. Create plotting copies of the complete trees
#
# No tips are removed here.
# ============================================================

scoly_plot_tree <- scoly_tree
platy_plot_tree <- platy_tree


scoly_plot_tree$tip.label <- clean_plot_id(
  scoly_plot_tree$tip.label
)


platy_plot_tree$tip.label <- clean_plot_id(
  platy_plot_tree$tip.label
)


if (
  anyDuplicated(
    scoly_plot_tree$tip.label
  ) > 0
) {
  
  stop(
    "Duplicated Scolytinae tree-tip labels after cleaning."
  )
}


if (
  anyDuplicated(
    platy_plot_tree$tip.label
  ) > 0
) {
  
  stop(
    "Duplicated Platypodinae tree-tip labels after cleaning."
  )
}


cat(
  "Scolytinae tree tips:",
  length(
    scoly_plot_tree$tip.label
  ),
  "\n"
)


cat(
  "Platypodinae tree tips:",
  length(
    platy_plot_tree$tip.label
  ),
  "\n"
)


# ============================================================
# 8. Identify tips without linked fungal data
# ============================================================

scoly_unlinked_tips <- setdiff(
  
  scoly_plot_tree$tip.label,
  alpha_plot_data$tip_id
)


platy_unlinked_tips <- setdiff(
  
  platy_plot_tree$tip.label,
  alpha_plot_data$tip_id
)


cat(
  "\nScolytinae tips without linked fungal data:\n"
)

print(
  scoly_unlinked_tips
)


cat(
  "\nPlatypodinae tips without linked fungal data:\n"
)

print(
  platy_unlinked_tips
)


cat(
  "\nNumber of unlinked Scolytinae tips:",
  length(scoly_unlinked_tips),
  "\n"
)


cat(
  "Number of unlinked Platypodinae tips:",
  length(platy_unlinked_tips),
  "\n"
)


# ============================================================
# 9. Global Tribe colour palette
#
# A fixed palette is defined once and reused for both trees,
# so every Tribe keeps exactly the same colour in all figures.
#
# Unknown is always white.
# ============================================================

all_tree_ids <- union(
  
  scoly_plot_tree$tip.label,
  platy_plot_tree$tip.label
)


tribes_present <- taxonomy_lookup %>%
  
  filter(
    tip_id %in% all_tree_ids
  ) %>%
  
  pull(
    Tribe
  ) %>%
  
  clean_plot_text() %>%
  
  unique() %>%
  
  na.omit() %>%
  
  sort()


known_tribes <- setdiff(
  tribes_present,
  "Unknown"
)


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


missing_tribe_colours <- setdiff(
  known_tribes,
  names(tribe_cols)
)


if (length(missing_tribe_colours) > 0L) {
  
  stop(
    "No fixed colour supplied for Tribe(s): ",
    paste(
      missing_tribe_colours,
      collapse = ", "
    )
  )
}


cat(
  "\nTribe colour palette:\n"
)

print(
  tribe_cols
)


# ============================================================
# 10. Region colours
# ============================================================

region_cols <- c(
  
  "French Guiana" =
    "#F9A882",
  
  "Borneo" =
    "#87C0E8",
  
  "No linked fungal data" =
    "white"
)


print(
  region_cols
)


# ============================================================
# 11. Shared continuous colour ranges
#
# The same colour corresponds to the same value in both trees.
# ============================================================

richness_limits <- get_metric_limits(
  alpha_plot_data$Richness
)


mpd_limits <- get_metric_limits(
  alpha_plot_data$MPD
)


mntd_limits <- get_metric_limits(
  alpha_plot_data$MNTD
)


cat(
  "\nRarefied richness range:\n"
)

print(
  richness_limits
)


cat(
  "\nMPD range:\n"
)

print(
  mpd_limits
)


cat(
  "\nMNTD range:\n"
)

print(
  mntd_limits
)


# ============================================================
# 12. Tip-label function
#
# Species available:
#   Tribe-Species-ID
#
# Genus available:
#   Tribe-Genus-ID
#
# Tribe only:
#   Tribe-ID
#
# Tribe unavailable:
#   Subfamily-ID
# ============================================================

make_tip_text <- function(
    subfamily,
    tribe,
    genus,
    species,
    tip_id
) {
  
  case_when(
    
    !is.na(tribe) &
      tribe != "Unknown" &
      !is.na(species) &
      species != "" ~
      
      paste(
        tribe,
        species,
        tip_id,
        sep = "-"
      ),
    
    
    !is.na(tribe) &
      tribe != "Unknown" &
      !is.na(genus) &
      genus != "" ~
      
      paste(
        tribe,
        genus,
        tip_id,
        sep = "-"
      ),
    
    
    !is.na(tribe) &
      tribe != "Unknown" ~
      
      paste(
        tribe,
        tip_id,
        sep = "-"
      ),
    
    
    !is.na(subfamily) &
      subfamily != "" ~
      
      paste(
        subfamily,
        tip_id,
        sep = "-"
      ),
    
    
    TRUE ~
      
      tip_id
  )
}


# ============================================================
# 13. Complete fruit-tree function
# ============================================================

make_fruit_tree <- function(
    
  tree,
  alpha_plot_data,
  taxonomy_lookup,
  
  subfamily_name,
  
  title_text,
  subtitle_text,
  
  output_prefix,
  
  tribe_palette,
  region_palette,
  
  richness_scale_limits,
  mpd_scale_limits,
  mntd_scale_limits,
  
  open_angle = 10,
  
  # Separate value can be supplied for each tree.
  first_ring_offset = 0.75,
  
  tip_size = 2.05,
  tip_offset = 0.03,
  
  # Separate point size can also be supplied for each tree.
  tip_point_size = 0.50,
  
  output_width_cm = 25,
  output_height_cm = 18
) {
  
  
  # ----------------------------------------------------------
  # Start with every tip in the complete tree
  # ----------------------------------------------------------
  
  dat <- tibble(
    
    tip_id =
      tree$tip.label
    
  ) %>%
    
    left_join(
      
      taxonomy_lookup,
      
      by =
        "tip_id"
    ) %>%
    
    left_join(
      
      alpha_plot_data,
      
      by =
        "tip_id"
    ) %>%
    
    mutate(
      
      Subfamily =
        coalesce(
          Subfamily,
          subfamily_name
        ),
      
      Tribe_plot =
        coalesce(
          Tribe,
          "Unknown"
        ),
      
      Country_plot =
        coalesce(
          Country_plot,
          "No linked fungal data"
        ),
      
      tip_text =
        make_tip_text(
          
          subfamily =
            Subfamily,
          
          tribe =
            Tribe_plot,
          
          genus =
            Genus,
          
          species =
            Species,
          
          tip_id =
            tip_id
        )
    ) %>%
    
    arrange(
      
      match(
        tip_id,
        tree$tip.label
      )
    )
  
  
  # ----------------------------------------------------------
  # Safety checks
  # ----------------------------------------------------------
  
  if (
    nrow(dat) !=
    length(tree$tip.label)
  ) {
    
    stop(
      title_text,
      ": plotting rows do not equal the number of tree tips."
    )
  }
  
  
  if (
    anyDuplicated(
      dat$tip_id
    ) > 0
  ) {
    
    stop(
      title_text,
      ": duplicated plotting tip IDs."
    )
  }
  
  
  if (
    !setequal(
      dat$tip_id,
      tree$tip.label
    )
  ) {
    
    stop(
      title_text,
      ": tree tips and plotting data do not match."
    )
  }
  
  
  missing_tribe_colours <- setdiff(
    
    unique(
      dat$Tribe_plot
    ),
    
    names(
      tribe_palette
    )
  )
  
  
  missing_region_colours <- setdiff(
    
    unique(
      dat$Country_plot
    ),
    
    names(
      region_palette
    )
  )
  
  
  if (
    length(
      missing_tribe_colours
    ) > 0
  ) {
    
    stop(
      
      "Missing Tribe colours for: ",
      
      paste(
        missing_tribe_colours,
        collapse = ", "
      )
    )
  }
  
  
  if (
    length(
      missing_region_colours
    ) > 0
  ) {
    
    stop(
      
      "Missing Region colours for: ",
      
      paste(
        missing_region_colours,
        collapse = ", "
      )
    )
  }
  
  
  n_unlinked <- sum(
    
    dat$Country_plot ==
      "No linked fungal data"
  )
  
  
  # ----------------------------------------------------------
  # Data attached directly to the tree
  #
  # IMPORTANT:
  # Tip_region is used ONLY for the small circles at tree tips.
  #
  # Country_plot remains separate and is used ONLY for the
  # Region fruit ring.
  # ----------------------------------------------------------
  
  label_dat <- dat %>%
    
    transmute(
      
      label =
        tip_id,
      
      tip_text =
        tip_text,
      
      Tip_region =
        Country_plot
    )
  
  
  # ----------------------------------------------------------
  # Ring datasets
  # ----------------------------------------------------------
  
  tribe_dat <- dat %>%
    
    transmute(
      
      tip_id,
      
      ring_value =
        1,
      
      Tribe_plot
    )
  
  
  region_dat <- dat %>%
    
    transmute(
      
      tip_id,
      
      ring_value =
        1,
      
      Country_plot
    )
  
  
  richness_dat <- dat %>%
    
    transmute(
      
      tip_id,
      
      ring_value =
        1,
      
      Richness
    )
  
  
  mpd_dat <- dat %>%
    
    transmute(
      
      tip_id,
      
      ring_value =
        1,
      
      MPD
    )
  
  
  mntd_dat <- dat %>%
    
    transmute(
      
      tip_id,
      
      ring_value =
        1,
      
      MNTD
    )
  
  
  # ----------------------------------------------------------
  # Legend entries actually present in the current tree
  # ----------------------------------------------------------
  
  tribe_legend_breaks <- intersect(
    names(tribe_palette),
    unique(tribe_dat$Tribe_plot)
  )
  
  
  region_legend_breaks <- intersect(
    names(region_palette),
    unique(region_dat$Country_plot)
  )
  
  
  # ----------------------------------------------------------
  # Figure caption
  # ----------------------------------------------------------
  
  caption_text <- paste(
    
    "Rings from inner to outer represent beetle tribe,",
    "sampling region, fungal rarefied richness, MPD, and MNTD."
  )
  
  
  if (n_unlinked > 0) {
    
    caption_text <- paste0(
      
      caption_text,
      
      " \nWhite cells indicate ",
      
      n_unlinked,
      
      " tree tip",
      
      ifelse(
        n_unlinked == 1,
        "",
        "s"
      ),
      
      " for which no linked fungal diversity data were available."
    )
  }
  
  
  # ==========================================================
  # Build plot
  # ==========================================================
  
  Tree <- ggtree(
    
    tree,
    
    layout =
      "fan",
    
    open.angle =
      open_angle,
    
    linewidth =
      0.25
    
  ) %<+%
    
    label_dat +
    
    
    # ========================================================
  # Ring 1: Tribe
  #
  # ORIGINAL geometry restored.
  # ========================================================
  
  ggtreeExtra::geom_fruit(
    
    data =
      tribe_dat,
    
    geom =
      geom_tile,
    
    mapping =
      aes(
        y = tip_id,
        x = ring_value,
        fill = Tribe_plot
      ),
    
    inherit.aes =
      FALSE,
    
    pwidth =
      0.030,
    
    offset =
      first_ring_offset,
    
    width =
      0.120,
    
    height =
      0.95,
    
    colour =
      "white",
    
    linewidth =
      0.05
  ) +
    
    
    scale_fill_manual(
      
      values =
        tribe_palette,
      
      breaks =
        tribe_legend_breaks,
      
      name =
        "Tribe",
      
      drop =
        FALSE,
      
      guide =
        guide_legend(
          
          order =
            1,
          
          ncol =
            1,
          
          byrow =
            TRUE,
          
          # Only the white "Unknown" key receives an outline.
          override.aes =
            list(
              
              colour =
                ifelse(
                  tribe_legend_breaks == "Unknown",
                  "#B8B8B8",
                  NA_character_
                ),
              
              linewidth =
                ifelse(
                  tribe_legend_breaks == "Unknown",
                  0.4,
                  0
                )
            )
        )
    ) +
    
    
    ggnewscale::new_scale_fill() +
    
    
    # ========================================================
  # Ring 2: Region
  #
  # ORIGINAL geometry restored.
  # Colours unchanged.
  # ========================================================
  
  ggtreeExtra::geom_fruit(
    
    data =
      region_dat,
    
    geom =
      geom_tile,
    
    mapping =
      aes(
        y = tip_id,
        x = ring_value,
        fill = Country_plot
      ),
    
    inherit.aes =
      FALSE,
    
    pwidth =
      0.030,
    
    offset =
      0.010,
    
    width =
      0.120,
    
    height =
      0.95,
    
    colour =
      "white",
    
    linewidth =
      0.05
  ) +
    
    
    scale_fill_manual(
      
      values =
        region_palette,
      
      breaks =
        region_legend_breaks,
      
      name =
        "Region",
      
      drop =
        FALSE,
      
      guide =
        guide_legend(
          
          order = 2,
          
          # Only the white missing-data key receives an outline.
          override.aes =
            list(
              
              colour =
                ifelse(
                  region_legend_breaks ==
                    "No linked fungal data",
                  "#B8B8B8",
                  NA_character_
                ),
              
              linewidth =
                ifelse(
                  region_legend_breaks ==
                    "No linked fungal data",
                  0.4,
                  0
                )
            )
        )
    ) +
    
    
    ggnewscale::new_scale_fill() +
    
    
    # ========================================================
  # Ring 3: Rarefied richness
  #
  # ORIGINAL geometry restored.
  # Colours unchanged.
  # ========================================================
  
  ggtreeExtra::geom_fruit(
    
    data =
      richness_dat,
    
    geom =
      geom_tile,
    
    mapping =
      aes(
        y = tip_id,
        x = ring_value,
        fill = Richness
      ),
    
    inherit.aes =
      FALSE,
    
    pwidth =
      0.032,
    
    offset =
      0.010,
    
    width =
      0.100,
    
    height =
      0.95,
    
    colour =
      "white",
    
    linewidth =
      0.05
  ) +
    
    
    scale_fill_gradient(
      
      low =
        "#FFF1BF",
      
      high =
        "#FFD23F",
      
      limits =
        richness_scale_limits,
      
      oob =
        scales::squish,
      
      na.value =
        "white",
      
      name =
        "Rarefied\nrichness",
      
      guide =
        guide_colorbar(
          
          order =
            3,
          
          barheight =
            grid::unit(
              18,
              "mm"
            )
        )
    ) +
    
    
    ggnewscale::new_scale_fill() +
    
    
    # ========================================================
  # Ring 4: MPD
  #
  # ORIGINAL geometry restored.
  # Colours unchanged.
  # ========================================================
  
  ggtreeExtra::geom_fruit(
    
    data =
      mpd_dat,
    
    geom =
      geom_tile,
    
    mapping =
      aes(
        y = tip_id,
        x = ring_value,
        fill = MPD
      ),
    
    inherit.aes =
      FALSE,
    
    pwidth =
      0.028,
    
    offset =
      0.010,
    
    width =
      0.087,
    
    height =
      0.95,
    
    colour =
      "white",
    
    linewidth =
      0.05
  ) +
    
    
    scale_fill_gradient(
      
      low =
        "#E2F0FA",
      
      high =
        "#66B2FF",
      
      limits =
        mpd_scale_limits,
      
      oob =
        scales::squish,
      
      na.value =
        "white",
      
      name =
        "MPD",
      
      guide =
        guide_colorbar(
          
          order =
            4,
          
          barheight =
            grid::unit(
              18,
              "mm"
            )
        )
    ) +
    
    
    ggnewscale::new_scale_fill() +
    
    
    # ========================================================
  # Ring 5: MNTD
  #
  # ORIGINAL geometry restored.
  # Colours unchanged.
  # ========================================================
  
  ggtreeExtra::geom_fruit(
    
    data =
      mntd_dat,
    
    geom =
      geom_tile,
    
    mapping =
      aes(
        y = tip_id,
        x = ring_value,
        fill = MNTD
      ),
    
    inherit.aes =
      FALSE,
    
    pwidth =
      0.014,
    
    offset =
      0.010,
    
    width =
      0.0447,
    
    height =
      0.95,
    
    colour =
      "white",
    
    linewidth =
      0.05
  ) +
    
    
    scale_fill_gradient(
      
      low =
        "#FCE4EC",
      
      high =
        "#FF8FB1",
      
      limits =
        mntd_scale_limits,
      
      oob =
        scales::squish,
      
      na.value =
        "white",
      
      name =
        "MNTD",
      
      guide =
        guide_colorbar(
          
          order =
            5,
          
          barheight =
            grid::unit(
              18,
              "mm"
            )
        )
    ) +
    
    
    # ========================================================
  # Redraw tree branches above the fruit layers
  # ========================================================
  
  geom_tree(
    linewidth =
      0.25
  ) +
    
    
    # ========================================================
  # Small Region points at tree tips
  #
  # Uses the SAME region colours as before.
  #
  # Point size is controlled separately for Scolytinae
  # and Platypodinae.
  # ========================================================
  
  geom_tippoint(
    
    aes(
      colour =
        Tip_region
    ),
    
    size =
      tip_point_size,
    
    alpha =
      0.90
  ) +
    
    
    scale_colour_manual(
      
      values =
        region_palette,
      
      na.value =
        "white",
      
      guide =
        "none"
    ) +
    
    
    # ========================================================
  # Tip labels
  #
  # Original label sizes and offsets can be supplied
  # separately for each tree.
  # ========================================================
  
  geom_tiplab2(
    
    aes(
      label =
        tip_text
    ),
    
    size =
      tip_size,
    
    offset =
      tip_offset,
    
    align =
      FALSE
  ) +
    
    
    # ========================================================
  # Titles and caption
  # ========================================================
  
  labs(
    
    title =
      title_text,
    
    subtitle =
      subtitle_text,
    
    caption =
      caption_text
  ) +
    
    
    # ========================================================
  # Theme
  # ========================================================
  
  theme(
    
    legend.position =
      "right",
    
    legend.box =
      "vertical",
    
    plot.title.position =
      "plot",
    
    plot.caption.position =
      "plot",
    
    plot.title =
      element_text(
        
        hjust =
          0.5,
        
        size =
          12.5,
        
        face =
          "bold",
        
        margin =
          margin(
            b = 3
          )
      ),
    
    plot.subtitle =
      element_text(
        
        hjust =
          0.5,
        
        size =
          10,
        
        margin =
          margin(
            b = 6
          )
      ),
    
    plot.caption =
      element_text(
        
        hjust =
          0,
        
        size =
          7.5,
        
        colour =
          "#444444",
        
        margin =
          margin(
            t = 6
          )
      ),
    
    legend.title =
      element_text(
        
        size =
          8.5,
        
        face =
          "bold"
      ),
    
    legend.text =
      element_text(
        size = 7.5
      ),
    
    legend.key.height =
      grid::unit(
        4.5,
        "mm"
      ),
    
    legend.spacing.y =
      grid::unit(
        2,
        "mm"
      ),
    
    plot.margin =
      margin(
        
        t = 8,
        r = 10,
        b = 8,
        l = 8,
        
        unit =
          "mm"
      )
  )
  
  
  # ----------------------------------------------------------
  # Display plot
  # ----------------------------------------------------------
  
  print(
    Tree
  )
  
  
  # ----------------------------------------------------------
  # Vector PDF
  # ----------------------------------------------------------
  
  ggsave(
    
    filename =
      paste0(
        output_prefix,
        ".pdf"
      ),
    
    plot =
      Tree,
    
    width =
      output_width_cm,
    
    height =
      output_height_cm,
    
    units =
      "cm",
    
    bg =
      "white",
    
    limitsize =
      FALSE
  )
  
  
  # ----------------------------------------------------------
  # 600-dpi PNG
  # ----------------------------------------------------------
  
  ggsave(
    
    filename =
      paste0(
        output_prefix,
        ".png"
      ),
    
    plot =
      Tree,
    
    width =
      output_width_cm,
    
    height =
      output_height_cm,
    
    units =
      "cm",
    
    dpi =
      600,
    
    bg =
      "white",
    
    limitsize =
      FALSE
  )
  
  
  # ----------------------------------------------------------
  # SVG
  # ----------------------------------------------------------
  
  if (
    requireNamespace(
      "svglite",
      quietly = TRUE
    )
  ) {
    
    ggsave(
      
      filename =
        paste0(
          output_prefix,
          ".svg"
        ),
      
      plot =
        Tree,
      
      width =
        output_width_cm,
      
      height =
        output_height_cm,
      
      units =
        "cm",
      
      device =
        svglite::svglite,
      
      bg =
        "white",
      
      limitsize =
        FALSE
    )
    
  } else {
    
    message(
      "SVG was not saved because package 'svglite' is not installed."
    )
  }
  
  
  return(
    Tree
  )
}


# ============================================================
# 14. Draw Scolytinae fruit tree
#
# Important:
# Scolytinae has many more tips, therefore:
#
# - first ring remains closer than Platypodinae
# - Region tip points are smaller
# ============================================================

scoly_fruit_tree <- make_fruit_tree(
  
  tree =
    scoly_plot_tree,
  
  alpha_plot_data =
    alpha_plot_data,
  
  taxonomy_lookup =
    taxonomy_lookup,
  
  subfamily_name =
    "Scolytinae",
  
  title_text =
    "Phylogenetic distribution of beetle-associated fungal diversity",
  
  subtitle_text =
    "Scolytinae sampled in French Guiana and Borneo",
  
  output_prefix =
    "Scolytinae_fruit_tree_macaron",
  
  tribe_palette =
    tribe_cols,
  
  region_palette =
    region_cols,
  
  richness_scale_limits =
    richness_limits,
  
  mpd_scale_limits =
    mpd_limits,
  
  mntd_scale_limits =
    mntd_limits,
  
  open_angle =
    10,
  
  # Original Scolytinae ring distance
  first_ring_offset =
    0.75,
  
  # Original Scolytinae label settings
  tip_size =
    1.15,
  
  tip_offset =
    0.015,
  
  # Small Region dots because the tree contains many tips
  tip_point_size =
    0.45,
  
  output_width_cm =
    24,
  
  output_height_cm =
    19
)


# ============================================================
# 15. Draw Platypodinae fruit tree
#
# Platypodinae has fewer tips, therefore:
#
# - the first ring remains further away
# - Region tip points can be slightly larger
# ============================================================

platy_fruit_tree <- make_fruit_tree(
  
  tree =
    platy_plot_tree,
  
  alpha_plot_data =
    alpha_plot_data,
  
  taxonomy_lookup =
    taxonomy_lookup,
  
  subfamily_name =
    "Platypodinae",
  
  title_text =
    "Phylogenetic distribution of beetle-associated fungal diversity",
  
  subtitle_text =
    "Platypodinae sampled in French Guiana and Borneo",
  
  output_prefix =
    "Platypodinae_fruit_tree_macaron",
  
  tribe_palette =
    tribe_cols,
  
  region_palette =
    region_cols,
  
  richness_scale_limits =
    richness_limits,
  
  mpd_scale_limits =
    mpd_limits,
  
  mntd_scale_limits =
    mntd_limits,
  
  open_angle =
    10,
  
  # Original Platypodinae ring distance
  first_ring_offset =
    1.5,
  
  # Original Platypodinae label settings
  tip_size =
    1.55,
  
  tip_offset =
    0.020,
  
  # Slightly larger dots because this tree has fewer tips
  tip_point_size =
    0.65,
  
  output_width_cm =
    24,
  
  output_height_cm =
    19
)

