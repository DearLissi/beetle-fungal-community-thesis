# ============================================================
# THESIS COPHYLOGENY V11
#
# OPPOSING SEMICIRCULAR PHYLOGENIES + TAXONOMIC BEETLE LABELS
#
# Four panels:
#   A. Scolytinae – Microascales
#   B. Scolytinae – Ophiostomatales
#   C. Platypodinae – Microascales
#   D. Platypodinae – Ophiostomatales
#
# Display rules copied from the final fruit-tree/tanglegram figures:
#   Tribe-Species-ID
#   Tribe-Genus-ID
#   Tribe-ID
#   Subfamily-ID
#
# Original beetle IDs and fungal OTU IDs remain the internal keys.
# Tree pruning, association construction, cophylo rotation, and all
# biological relationships therefore remain unchanged.
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(readxl)
  library(tidyverse)
})


# ============================================================
# 0. PATHS
# ============================================================

data_dir <- Sys.getenv("BEETLE_DATA_DIR", unset = "data")

default_plot_dir <- file.path(
  data_dir,
  "Cophylogeny_THESIS_v11_semicircle"
)

plot_dir <- Sys.getenv(
  "COPHYLOGENY_OUTPUT_DIR",
  unset = default_plot_dir
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. FIXED PROJECT COLOURS AND GEOMETRY
# ============================================================

tribe_cols <- c(
  "Corthylini" = "#9BC2F1",
  "Cryphalini" = "#F7C283",
  "Dryocoetini" = "#F4E085",
  "Platypodini" = "#B4E197",
  "Scolytoplatypodini" = "#8CD3C5",
  "Tesserocerini" = "#F28E8E",
  "Trypophloeini" = "#C7B2E3",
  "Xyleborini" = "#F3B2C8",
  "Xyloterini" = "#D9B8A2",
  "Unknown" = "white"
)

TREE_MIXED_COLOUR <- "#A8A8A8"
LABEL_COLOUR <- "#202020"
TIP_BORDER_COLOUR <- "#444444"

TIP_RADIUS <- 1.00
ROOT_RADIUS <- 1.34
LABEL_RADIUS <- 0.965


# ============================================================
# 2. CLEANING HELPERS
# ============================================================

clean_tax <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c(
    "", "NA", "N/A", "#N/A", "0",
    "na", "n/a", "NaN", "nan"
  )] <- NA_character_
  x
}

clean_id <- function(x) {
  toupper(clean_tax(x))
}

has_text <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}


# ============================================================
# 3. READ BEETLE AND FUNGAL TREES
# ============================================================

scoly_tree <- read.tree(
  file.path(data_dir, "scoly_new_metabarcode_samples.tre")
)

platy_tree <- read.tree(
  file.path(data_dir, "platy_new_metabarcode_samples.tre")
)

scoly_tree$tip.label <- clean_id(scoly_tree$tip.label)
platy_tree$tip.label <- clean_id(platy_tree$tip.label)

ophio_tree_full <- read.tree(
  file.path(
    data_dir,
    "Ophiostomatales_OTUs_FG_Borneo_constr_taxonomised.treefile"
  )
)

micro_tree_full <- read.tree(
  file.path(
    data_dir,
    "Microascales_OTUs_FG_Borneo_constr_taxonomised.treefile"
  )
)

ophio_tree <- keep.tip(
  ophio_tree_full,
  grep("A00", ophio_tree_full$tip.label, value = TRUE)
)

micro_tree <- keep.tip(
  micro_tree_full,
  grep("A00", micro_tree_full$tip.label, value = TRUE)
)


# ============================================================
# 4. READ OTU TABLE
# ============================================================

OTU_raw <- read_tsv(
  file.path(data_dir, "OTUsotu_table.tsv"),
  show_col_types = FALSE
)

OTU_counts <- OTU_raw %>%
  slice(-1) %>%
  mutate(
    OTU_ID = gsub(":", "_", OTU_ID)
  ) %>%
  mutate(
    across(-OTU_ID, as.numeric)
  )

OTU_by_sample <- OTU_counts %>%
  column_to_rownames("OTU_ID") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("fungal_metabarcode_ID")


# ============================================================
# 5. READ BEETLE TAXONOMY AND BUILD DISPLAY LABELS
# ============================================================

beetle_id <- read_xlsx(
  file.path(data_dir, "beetles id.xlsx")
) %>%
  mutate(
    fungal_metabarcode_ID = clean_tax(fungal_metabarcode_ID),
    mt_id = clean_id(mt_id),
    project_asv_id = clean_id(project_asv_id),
    subfamily = clean_tax(subfamily),
    tribe = clean_tax(tribe),
    genus = clean_tax(genus),
    species = clean_tax(species)
  )

make_tip_text <- function(
    subfamily,
    tribe,
    genus,
    species,
    tip_id
) {
  subfamily <- as.character(subfamily)
  tribe <- as.character(tribe)
  genus <- as.character(genus)
  species <- as.character(species)
  tip_id <- as.character(tip_id)
  
  output <- tip_id
  valid_tribe <- has_text(tribe) & tribe != "Unknown"
  valid_subfamily <- has_text(subfamily)
  valid_genus <- has_text(genus)
  valid_species <- has_text(species)
  
  use_subfamily <- !valid_tribe & valid_subfamily
  output[use_subfamily] <- paste(
    subfamily[use_subfamily],
    tip_id[use_subfamily],
    sep = "-"
  )
  
  use_tribe <- valid_tribe
  output[use_tribe] <- paste(
    tribe[use_tribe],
    tip_id[use_tribe],
    sep = "-"
  )
  
  use_genus <- valid_tribe & valid_genus
  output[use_genus] <- paste(
    tribe[use_genus],
    genus[use_genus],
    tip_id[use_genus],
    sep = "-"
  )
  
  use_species <- valid_tribe & valid_species
  output[use_species] <- paste(
    tribe[use_species],
    species[use_species],
    tip_id[use_species],
    sep = "-"
  )
  
  output
}

# Either mt_id or project_asv_id can occur in a beetle tree.
# When an ID occurs more than once, retain the row with the most
# resolved taxonomy (species > genus > tribe > subfamily).
beetle_taxonomy <- bind_rows(
  beetle_id %>%
    transmute(
      beetle = mt_id,
      subfamily,
      tribe,
      genus,
      species
    ),
  beetle_id %>%
    transmute(
      beetle = project_asv_id,
      subfamily,
      tribe,
      genus,
      species
    )
) %>%
  filter(!is.na(beetle)) %>%
  mutate(
    taxonomy_score =
      8L * as.integer(has_text(species)) +
      4L * as.integer(has_text(genus)) +
      2L * as.integer(has_text(tribe) & tribe != "Unknown") +
      as.integer(has_text(subfamily))
  ) %>%
  arrange(beetle, desc(taxonomy_score)) %>%
  group_by(beetle) %>%
  slice(1L) %>%
  ungroup() %>%
  mutate(
    display_label = make_tip_text(
      subfamily,
      tribe,
      genus,
      species,
      beetle
    ),
    label_rule = case_when(
      has_text(tribe) & tribe != "Unknown" & has_text(species) ~
        "Tribe-Species-ID",
      has_text(tribe) & tribe != "Unknown" & has_text(genus) ~
        "Tribe-Genus-ID",
      has_text(tribe) & tribe != "Unknown" ~
        "Tribe-ID",
      has_text(subfamily) ~
        "Subfamily-ID",
      TRUE ~
        "ID only"
    )
  )

if (anyDuplicated(beetle_taxonomy$beetle) > 0L) {
  stop("The beetle taxonomy lookup contains duplicated IDs.")
}

if (anyDuplicated(beetle_taxonomy$display_label) > 0L) {
  stop("The taxonomic display labels are not unique.")
}

display_lookup <- setNames(
  beetle_taxonomy$display_label,
  beetle_taxonomy$beetle
)

tribe_lookup <- setNames(
  beetle_taxonomy$tribe,
  beetle_taxonomy$beetle
)

get_beetle_display_label <- function(ids) {
  ids <- as.character(ids)
  output <- unname(display_lookup[ids])
  output[is.na(output)] <- ids[is.na(output)]
  output
}

get_beetle_tribe <- function(ids) {
  output <- unname(tribe_lookup[as.character(ids)])
  output[!has_text(output)] <- "Unknown"
  output
}


# ============================================================
# 6. BUILD ONE COPHYLOGENY DATASET
# ============================================================

make_cophylo_data <- function(
    beetle_tree,
    fungal_tree,
    fungal_group
) {
  cat("\nBuilding ", fungal_group, "\n", sep = "")
  
  beetle_map <- beetle_id %>%
    mutate(
      join_key = case_when(
        !is.na(mt_id) & mt_id %in% beetle_tree$tip.label ~ mt_id,
        !is.na(project_asv_id) &
          project_asv_id %in% beetle_tree$tip.label ~ project_asv_id,
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(fungal_metabarcode_ID),
      !is.na(join_key)
    ) %>%
    select(fungal_metabarcode_ID, join_key) %>%
    distinct()
  
  joined <- beetle_map %>%
    left_join(
      OTU_by_sample,
      by = "fungal_metabarcode_ID"
    )
  
  fungal_OTUs <- intersect(
    fungal_tree$tip.label,
    colnames(joined)
  )
  
  if (length(fungal_OTUs) == 0L) {
    stop(paste0("No fungal OTUs matched for ", fungal_group, "."))
  }
  
  beetle_OTU <- joined %>%
    select(join_key, all_of(fungal_OTUs)) %>%
    group_by(join_key) %>%
    summarise(
      across(everything(), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    filter(
      rowSums(across(all_of(fungal_OTUs))) > 0
    )
  
  beetle_OTU_mat <- beetle_OTU %>%
    column_to_rownames("join_key") %>%
    as.matrix()
  
  beetle_OTU_mat <- beetle_OTU_mat[
    ,
    colSums(beetle_OTU_mat) > 0,
    drop = FALSE
  ]
  
  association_index <- which(
    beetle_OTU_mat > 0,
    arr.ind = TRUE
  )
  
  assoc <- cbind(
    rownames(beetle_OTU_mat)[association_index[, 1]],
    colnames(beetle_OTU_mat)[association_index[, 2]]
  )
  assoc <- as.matrix(assoc)
  colnames(assoc) <- c("beetle", "fungus")
  
  output <- list(
    beetle_tree = keep.tip(
      beetle_tree,
      rownames(beetle_OTU_mat)
    ),
    fungal_tree = keep.tip(
      fungal_tree,
      colnames(beetle_OTU_mat)
    ),
    association = assoc,
    matrix = beetle_OTU_mat
  )
  
  cat(
    "  beetles: ", Ntip(output$beetle_tree),
    "; fungal OTUs: ", Ntip(output$fungal_tree),
    "; associations: ", nrow(output$association),
    "\n",
    sep = ""
  )
  
  output
}


# ============================================================
# 7. CREATE ALL FOUR DATASETS
# ============================================================

platy_ophio <- make_cophylo_data(
  platy_tree,
  ophio_tree,
  "Platypodinae – Ophiostomatales"
)

platy_micro <- make_cophylo_data(
  platy_tree,
  micro_tree,
  "Platypodinae – Microascales"
)

scoly_ophio <- make_cophylo_data(
  scoly_tree,
  ophio_tree,
  "Scolytinae – Ophiostomatales"
)

scoly_micro <- make_cophylo_data(
  scoly_tree,
  micro_tree,
  "Scolytinae – Microascales"
)


# ============================================================
# 8. COLOUR HELPERS
# ============================================================

state_colour <- function(states) {
  states <- as.character(states)
  colours <- unname(tribe_cols[states])
  bad <- is.na(colours) | is.na(states) | states == "Mixed"
  colours[bad] <- TREE_MIXED_COLOUR
  colours
}

make_edge_colours <- function(tree, tip_states) {
  n_tip <- Ntip(tree)
  tip_states <- tip_states[tree$tip.label]
  edge_colours <- rep(TREE_MIXED_COLOUR, nrow(tree$edge))
  
  for (i in seq_len(nrow(tree$edge))) {
    child <- tree$edge[i, 2]
    
    if (child <= n_tip) {
      state <- tip_states[tree$tip.label[child]]
    } else {
      descendants <- phytools::getDescendants(tree, child)
      descendant_tips <- descendants[descendants <= n_tip]
      states <- unique(
        tip_states[tree$tip.label[descendant_tips]]
      )
      states <- states[!is.na(states)]
      state <- if (length(states) == 1L) states else "Mixed"
    }
    
    edge_colours[i] <- state_colour(state)
  }
  
  edge_colours
}


# ============================================================
# 9. ROTATE AND STYLE ONE COPHYLOGENY
# ============================================================

prepare_cophylo_panel <- function(x) {
  obj <- cophylo(
    x$beetle_tree,
    x$fungal_tree,
    assoc = x$association,
    rotate = TRUE
  )
  
  left_tip_states <- get_beetle_tribe(
    obj$trees[[1]]$tip.label
  )
  names(left_tip_states) <- obj$trees[[1]]$tip.label
  
  assoc_df <- tibble(
    beetle = obj$assoc[, 1],
    fungus = obj$assoc[, 2],
    tribe = get_beetle_tribe(obj$assoc[, 1])
  )
  
  fungal_states <- assoc_df %>%
    group_by(fungus) %>%
    summarise(
      fungal_state = {
        known <- unique(tribe[tribe != "Unknown"])
        if (length(known) == 1L) {
          known
        } else if (length(known) == 0L) {
          "Unknown"
        } else {
          "Mixed"
        }
      },
      .groups = "drop"
    )
  
  right_tip_states <- setNames(
    fungal_states$fungal_state,
    fungal_states$fungus
  )[obj$trees[[2]]$tip.label]
  right_tip_states[is.na(right_tip_states)] <- "Unknown"
  names(right_tip_states) <- obj$trees[[2]]$tip.label
  
  used_tribes <- unique(left_tip_states)
  legend_tribes <- names(tribe_cols)[
    names(tribe_cols) %in% used_tribes &
      names(tribe_cols) != "Unknown"
  ]
  
  list(
    obj = obj,
    left_edge_colours = make_edge_colours(
      obj$trees[[1]],
      left_tip_states
    ),
    right_edge_colours = make_edge_colours(
      obj$trees[[2]],
      right_tip_states
    ),
    left_tip_colours = state_colour(left_tip_states),
    right_tip_colours = state_colour(right_tip_states),
    link_colours = state_colour(
      get_beetle_tribe(obj$assoc[, 1])
    ),
    legend_tribes = legend_tribes
  )
}


# ============================================================
# 10. SEMICIRCULAR PHYLOGENY COORDINATES
# ============================================================

build_phylo_semicircle_layout <- function(
    tree,
    edge_colours,
    start_angle,
    end_angle,
    tip_radius = TIP_RADIUS,
    root_radius = ROOT_RADIUS
) {
  n_tip <- Ntip(tree)
  n_node <- tree$Nnode
  total_nodes <- n_tip + n_node
  
  # The plotting order reflects rotations made by cophylo().
  plot_coordinates <- ape:::plotPhyloCoor(
    tree,
    type = "phylogram",
    use.edge.length = !is.null(tree$edge.length),
    direction = "rightwards"
  )
  
  tip_order <- tree$tip.label[
    order(plot_coordinates[seq_len(n_tip), 2])
  ]
  
  tip_angles <- seq(
    start_angle,
    end_angle,
    length.out = n_tip
  )
  names(tip_angles) <- tip_order
  
  theta <- rep(NA_real_, total_nodes)
  theta[seq_len(n_tip)] <- unname(
    tip_angles[tree$tip.label]
  )
  
  child_list <- split(
    tree$edge[, 2],
    tree$edge[, 1]
  )
  
  descendant_tip_cache <- vector("list", total_nodes)
  
  descendant_tips <- function(node) {
    if (!is.null(descendant_tip_cache[[node]])) {
      return(descendant_tip_cache[[node]])
    }
    
    if (node <= n_tip) {
      output <- node
    } else {
      children <- child_list[[as.character(node)]]
      output <- unlist(
        lapply(children, descendant_tips),
        use.names = FALSE
      )
    }
    
    descendant_tip_cache[[node]] <<- output
    output
  }
  
  internal_nodes <- seq.int(n_tip + 1L, total_nodes)
  for (node in internal_nodes) {
    tips <- descendant_tips(node)
    theta[node] <- mean(range(theta[tips]))
  }
  
  if (!is.null(tree$edge.length)) {
    node_depth <- node.depth.edgelength(tree)
  } else {
    raw_depth <- node.depth(tree, method = 1)
    node_depth <- max(raw_depth) - raw_depth
  }
  
  root_node <- setdiff(
    unique(tree$edge[, 1]),
    unique(tree$edge[, 2])
  )[1]
  node_depth <- node_depth - node_depth[root_node]
  maximum_tip_depth <- max(node_depth[seq_len(n_tip)], na.rm = TRUE)
  
  if (!is.finite(maximum_tip_depth) || maximum_tip_depth <= 0) {
    maximum_tip_depth <- 1
  }
  
  radius <- root_radius -
    (node_depth / maximum_tip_depth) *
    (root_radius - tip_radius)
  
  radius[seq_len(n_tip)] <- tip_radius
  radius[root_node] <- root_radius
  radius <- pmax(tip_radius, pmin(root_radius, radius))
  
  nodes <- data.frame(
    node = seq_len(total_nodes),
    theta = theta,
    radius = radius,
    is_tip = seq_len(total_nodes) <= n_tip,
    tip_id = c(
      tree$tip.label,
      rep(NA_character_, n_node)
    ),
    stringsAsFactors = FALSE
  )
  
  edges <- data.frame(
    parent = tree$edge[, 1],
    child = tree$edge[, 2],
    colour = edge_colours,
    stringsAsFactors = FALSE
  )
  
  tips <- nodes[nodes$is_tip, , drop = FALSE]
  tips <- tips[match(tip_order, tips$tip_id), , drop = FALSE]
  rownames(tips) <- NULL
  
  list(
    nodes = nodes,
    edges = edges,
    tips = tips,
    tip_order = tip_order
  )
}


# ============================================================
# 11. DRAWING HELPERS
# ============================================================

polar_x <- function(radius, theta) radius * cos(theta)
polar_y <- function(radius, theta) radius * sin(theta)

draw_arc <- function(radius, theta_1, theta_2, colour, lwd) {
  n_points <- max(3L, ceiling(abs(theta_2 - theta_1) * 90))
  theta <- seq(theta_1, theta_2, length.out = n_points)
  lines(
    polar_x(radius, theta),
    polar_y(radius, theta),
    col = colour,
    lwd = lwd
  )
}

draw_semicircle_tree <- function(
    layout,
    tip_colours,
    branch_lwd,
    tip_cex
) {
  nodes <- layout$nodes
  
  for (i in seq_len(nrow(layout$edges))) {
    edge <- layout$edges[i, ]
    parent <- nodes[nodes$node == edge$parent, , drop = FALSE]
    child <- nodes[nodes$node == edge$child, , drop = FALSE]
    
    draw_arc(
      parent$radius,
      parent$theta,
      child$theta,
      edge$colour,
      branch_lwd
    )
    
    segments(
      polar_x(parent$radius, child$theta),
      polar_y(parent$radius, child$theta),
      polar_x(child$radius, child$theta),
      polar_y(child$radius, child$theta),
      col = edge$colour,
      lwd = branch_lwd
    )
  }
  
  tip_colour_lookup <- setNames(
    tip_colours,
    layout$nodes$tip_id[layout$nodes$is_tip]
  )
  
  points(
    polar_x(layout$tips$radius, layout$tips$theta),
    polar_y(layout$tips$radius, layout$tips$theta),
    pch = 21,
    bg = unname(tip_colour_lookup[layout$tips$tip_id]),
    col = TIP_BORDER_COLOUR,
    cex = tip_cex,
    lwd = 0.55
  )
}

cubic_bezier <- function(p0, p1, p2, p3, n = 80L) {
  t <- seq(0, 1, length.out = n)
  one_minus_t <- 1 - t
  
  list(
    x = one_minus_t^3 * p0[1] +
      3 * one_minus_t^2 * t * p1[1] +
      3 * one_minus_t * t^2 * p2[1] +
      t^3 * p3[1],
    y = one_minus_t^3 * p0[2] +
      3 * one_minus_t^2 * t * p1[2] +
      3 * one_minus_t * t^2 * p2[2] +
      t^3 * p3[2]
  )
}

label_rotation <- function(theta) {
  degrees <- (theta * 180 / pi) %% 360
  left_half <- cos(theta) < 0
  degrees[left_half] <- degrees[left_half] + 180
  degrees
}

label_adjustment <- function(theta) {
  ifelse(cos(theta) >= 0, 1, 0)
}

draw_tip_labels <- function(layout, display_labels, label_cex) {
  display_lookup_local <- setNames(
    display_labels,
    layout$nodes$tip_id[layout$nodes$is_tip]
  )
  ordered_labels <- unname(
    display_lookup_local[layout$tips$tip_id]
  )
  
  for (i in seq_len(nrow(layout$tips))) {
    theta <- layout$tips$theta[i]
    text(
      polar_x(LABEL_RADIUS, theta),
      polar_y(LABEL_RADIUS, theta),
      labels = ordered_labels[i],
      srt = label_rotation(theta),
      adj = c(label_adjustment(theta), 0.5),
      cex = label_cex,
      col = LABEL_COLOUR,
      family = "sans"
    )
  }
}

draw_association_links <- function(
    upper_layout,
    lower_layout,
    association,
    upper_display_labels,
    lower_display_labels,
    upper_label_cex,
    lower_label_cex,
    link_colours,
    link_lwd,
    link_alpha
) {
  upper_display_lookup <- setNames(
    upper_display_labels,
    upper_layout$nodes$tip_id[upper_layout$nodes$is_tip]
  )
  lower_display_lookup <- setNames(
    lower_display_labels,
    lower_layout$nodes$tip_id[lower_layout$nodes$is_tip]
  )
  
  upper_widths <- strwidth(
    unname(upper_display_lookup[upper_layout$tips$tip_id]),
    units = "user",
    cex = upper_label_cex,
    family = "sans"
  )
  lower_widths <- strwidth(
    unname(lower_display_lookup[lower_layout$tips$tip_id]),
    units = "user",
    cex = lower_label_cex,
    family = "sans"
  )
  
  upper_endpoint_radius <- setNames(
    pmax(0.12, LABEL_RADIUS - upper_widths + 0.002),
    upper_layout$tips$tip_id
  )
  lower_endpoint_radius <- setNames(
    pmax(0.12, LABEL_RADIUS - lower_widths + 0.002),
    lower_layout$tips$tip_id
  )
  upper_theta <- setNames(
    upper_layout$tips$theta,
    upper_layout$tips$tip_id
  )
  lower_theta <- setNames(
    lower_layout$tips$theta,
    lower_layout$tips$tip_id
  )
  
  visible_colours <- adjustcolor(
    link_colours,
    alpha.f = link_alpha
  )
  
  for (i in seq_len(nrow(association))) {
    beetle <- association[i, 1]
    fungus <- association[i, 2]
    theta_upper <- unname(upper_theta[beetle])
    theta_lower <- unname(lower_theta[fungus])
    radius_upper <- unname(upper_endpoint_radius[beetle])
    radius_lower <- unname(lower_endpoint_radius[fungus])
    
    if (
      any(!is.finite(c(
        theta_upper,
        theta_lower,
        radius_upper,
        radius_lower
      )))
    ) {
      stop("An association endpoint could not be matched to a plotted tip.")
    }
    
    p0 <- c(
      polar_x(radius_upper, theta_upper),
      polar_y(radius_upper, theta_upper)
    )
    p3 <- c(
      polar_x(radius_lower, theta_lower),
      polar_y(radius_lower, theta_lower)
    )
    p1 <- p0 * 0.28
    p2 <- p3 * 0.28
    curve <- cubic_bezier(p0, p1, p2, p3)
    
    lines(
      curve$x,
      curve$y,
      col = visible_colours[i],
      lwd = link_lwd
    )
  }
}


# ============================================================
# 12. BOTTOM LEGEND
# ============================================================

draw_bottom_legend <- function(legend_tribes, legend_cex) {
  par(mar = c(0.2, 1.2, 0.2, 1.2), family = "sans", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  text(
    0.035,
    0.66,
    "Lines",
    adj = c(0, 0.5),
    cex = legend_cex,
    font = 2,
    col = LABEL_COLOUR
  )
  segments(
    0.105,
    0.66,
    0.16,
    0.66,
    col = adjustcolor("#777777", alpha.f = 0.72),
    lwd = 1.7
  )
  text(
    0.18,
    0.66,
    "connect beetle taxa to associated fungal OTUs; colours indicate beetle tribe.",
    adj = c(0, 0.5),
    cex = legend_cex * 0.86,
    col = LABEL_COLOUR
  )
  
  text(
    0.035,
    0.24,
    "Tribe",
    adj = c(0, 0.5),
    cex = legend_cex,
    font = 2,
    col = LABEL_COLOUR
  )
  
  tribe_x <- seq(
    0.17,
    0.92,
    length.out = length(legend_tribes)
  )
  points(
    tribe_x,
    rep(0.24, length(tribe_x)),
    pch = 21,
    bg = unname(tribe_cols[legend_tribes]),
    col = TIP_BORDER_COLOUR,
    cex = legend_cex * 0.95
  )
  text(
    tribe_x + 0.012,
    rep(0.24, length(tribe_x)),
    labels = legend_tribes,
    adj = c(0, 0.5),
    cex = legend_cex * 0.76,
    col = LABEL_COLOUR
  )
}


# ============================================================
# 13. SAVE ONE COMPLETE PANEL
# ============================================================

save_semicircle_cophylogeny <- function(
    x,
    panel_letter,
    panel_title,
    file_stem,
    beetle_label_cex,
    fungal_label_cex,
    link_lwd,
    link_alpha,
    branch_lwd,
    tip_cex,
    legend_cex,
    png_width,
    png_height,
    pdf_width,
    pdf_height
) {
  styled <- prepare_cophylo_panel(x)
  upper_tree <- styled$obj$trees[[1]]
  lower_tree <- styled$obj$trees[[2]]
  association <- styled$obj$assoc
  
  missing_taxonomy <- setdiff(
    upper_tree$tip.label,
    names(display_lookup)
  )
  if (length(missing_taxonomy) > 0L) {
    stop(
      paste0(
        panel_title,
        ": beetle IDs missing from taxonomy: ",
        paste(missing_taxonomy, collapse = ", ")
      )
    )
  }
  
  angle_gap <- 3.5 * pi / 180
  upper_layout <- build_phylo_semicircle_layout(
    upper_tree,
    styled$left_edge_colours,
    start_angle = pi - angle_gap,
    end_angle = angle_gap
  )
  lower_layout <- build_phylo_semicircle_layout(
    lower_tree,
    styled$right_edge_colours,
    start_angle = pi + angle_gap,
    end_angle = 2 * pi - angle_gap
  )
  
  beetle_display_labels <- get_beetle_display_label(
    upper_tree$tip.label
  )
  fungal_display_labels <- lower_tree$tip.label
  
  draw_complete_figure <- function() {
    layout(
      matrix(c(1, 2), nrow = 2),
      heights = c(15.2, 1.8)
    )
    
    par(
      mar = c(0.1, 0.7, 5.5, 0.7),
      family = "sans",
      xpd = NA
    )
    plot.new()
    plot.window(
      xlim = c(-1.58, 1.58),
      ylim = c(-1.58, 1.58),
      asp = 1
    )
    
    mtext(
      "Phylogenetic distribution of fungal associations across bark beetle lineages",
      side = 3,
      line = 3.6,
      cex = 1.40,
      font = 2,
      col = LABEL_COLOUR
    )
    mtext(
      paste0(
        "The upper semicircle shows beetle phylogeny and the lower semicircle shows fungal phylogeny; ",
        "curves indicate observed beetle–fungus associations."
      ),
      side = 3,
      line = 2.0,
      cex = 0.82,
      col = LABEL_COLOUR
    )
    mtext(
      panel_letter,
      side = 3,
      line = 0.25,
      adj = 0.02,
      cex = 1.28,
      font = 2,
      col = LABEL_COLOUR
    )
    mtext(
      panel_title,
      side = 3,
      line = 0.25,
      adj = 0.5,
      cex = 1.08,
      font = 2,
      col = LABEL_COLOUR
    )
    
    text(
      0,
      1.46,
      "Beetle phylogenetic relationships",
      cex = 0.76,
      font = 2,
      col = LABEL_COLOUR
    )
    text(
      0,
      -1.46,
      "Fungal phylogenetic relationships",
      cex = 0.76,
      font = 2,
      col = LABEL_COLOUR
    )
    
    draw_association_links(
      upper_layout = upper_layout,
      lower_layout = lower_layout,
      association = association,
      upper_display_labels = beetle_display_labels,
      lower_display_labels = fungal_display_labels,
      upper_label_cex = beetle_label_cex,
      lower_label_cex = fungal_label_cex,
      link_colours = styled$link_colours,
      link_lwd = link_lwd,
      link_alpha = link_alpha
    )
    
    draw_tip_labels(
      upper_layout,
      beetle_display_labels,
      beetle_label_cex
    )
    draw_tip_labels(
      lower_layout,
      fungal_display_labels,
      fungal_label_cex
    )
    
    draw_semicircle_tree(
      upper_layout,
      styled$left_tip_colours,
      branch_lwd,
      tip_cex
    )
    draw_semicircle_tree(
      lower_layout,
      styled$right_tip_colours,
      branch_lwd,
      tip_cex
    )
    
    draw_bottom_legend(
      styled$legend_tribes,
      legend_cex
    )
  }
  
  png_file <- file.path(plot_dir, paste0(file_stem, ".png"))
  pdf_file <- file.path(plot_dir, paste0(file_stem, ".pdf"))
  
  png(
    png_file,
    width = png_width,
    height = png_height,
    res = 300,
    bg = "white"
  )
  draw_complete_figure()
  dev.off()
  
  pdf(
    pdf_file,
    width = pdf_width,
    height = pdf_height,
    family = "Helvetica",
    useDingbats = FALSE
  )
  draw_complete_figure()
  dev.off()
  
  tibble(
    Panel = panel_letter,
    Comparison = panel_title,
    n_beetles = Ntip(upper_tree),
    n_fungal_OTUs = Ntip(lower_tree),
    n_associations = nrow(association),
    PNG = png_file,
    PDF = pdf_file
  )
}


# ============================================================
# 14. GENERATE THE FOUR FIGURES
#
# Dense panels use slightly lower opacity than sparse panels, but
# every connector retains the original fixed tribe colour and is
# thicker than in the old rectangular cophylogeny figures.
# ============================================================

summary_A <- save_semicircle_cophylogeny(
  x = scoly_micro,
  panel_letter = "A",
  panel_title = "Scolytinae – Microascales",
  file_stem = "A_Scolytinae_Microascales_COPHYLO_v11_semicircle",
  beetle_label_cex = 0.285,
  fungal_label_cex = 0.39,
  link_lwd = 1.00,
  link_alpha = 0.52,
  branch_lwd = 1.00,
  tip_cex = 0.45,
  legend_cex = 0.78,
  png_width = 6500,
  png_height = 6800,
  pdf_width = 19.5,
  pdf_height = 20.4
)

summary_B <- save_semicircle_cophylogeny(
  x = scoly_ophio,
  panel_letter = "B",
  panel_title = "Scolytinae – Ophiostomatales",
  file_stem = "B_Scolytinae_Ophiostomatales_COPHYLO_v11_semicircle",
  beetle_label_cex = 0.275,
  fungal_label_cex = 0.275,
  link_lwd = 0.90,
  link_alpha = 0.46,
  branch_lwd = 0.92,
  tip_cex = 0.41,
  legend_cex = 0.77,
  png_width = 7000,
  png_height = 7300,
  pdf_width = 20.5,
  pdf_height = 21.4
)

summary_C <- save_semicircle_cophylogeny(
  x = platy_micro,
  panel_letter = "C",
  panel_title = "Platypodinae – Microascales",
  file_stem = "C_Platypodinae_Microascales_COPHYLO_v11_semicircle",
  beetle_label_cex = 0.47,
  fungal_label_cex = 0.48,
  link_lwd = 1.65,
  link_alpha = 0.82,
  branch_lwd = 1.22,
  tip_cex = 0.64,
  legend_cex = 0.87,
  png_width = 5700,
  png_height = 6000,
  pdf_width = 18.0,
  pdf_height = 18.8
)

summary_D <- save_semicircle_cophylogeny(
  x = platy_ophio,
  panel_letter = "D",
  panel_title = "Platypodinae – Ophiostomatales",
  file_stem = "D_Platypodinae_Ophiostomatales_COPHYLO_v11_semicircle",
  beetle_label_cex = 0.43,
  fungal_label_cex = 0.36,
  link_lwd = 1.30,
  link_alpha = 0.66,
  branch_lwd = 1.12,
  tip_cex = 0.56,
  legend_cex = 0.84,
  png_width = 6200,
  png_height = 6500,
  pdf_width = 18.8,
  pdf_height = 19.6
)

cophylogeny_v11_summary <- bind_rows(
  summary_A,
  summary_B,
  summary_C,
  summary_D
)

write_csv(
  cophylogeny_v11_summary,
  file.path(
    plot_dir,
    "Cophylogeny_THESIS_v11_semicircle_summary.csv"
  )
)

used_beetle_ids <- unique(c(
  scoly_micro$beetle_tree$tip.label,
  scoly_ophio$beetle_tree$tip.label,
  platy_micro$beetle_tree$tip.label,
  platy_ophio$beetle_tree$tip.label
))

write_csv(
  beetle_taxonomy %>%
    filter(beetle %in% used_beetle_ids) %>%
    select(
      beetle,
      subfamily,
      tribe,
      genus,
      species,
      label_rule,
      display_label
    ),
  file.path(
    plot_dir,
    "Cophylogeny_taxonomic_tip_label_key.csv"
  )
)

cat("\nTHESIS COPHYLOGENY V11 COMPLETE\n")
cat("Saved to:\n", plot_dir, "\n\n", sep = "")
print(cophylogeny_v11_summary, n = Inf, width = Inf)

