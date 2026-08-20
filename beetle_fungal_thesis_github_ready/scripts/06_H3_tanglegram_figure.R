# ============================================================
# THESIS TANGLEGRAM V11
#
# SEMICIRCULAR TANGLEGRAM + TAXONOMIC TIP LABELS
#
# This script loads the H3 analysis objects created by script 05. It does
# not recalculate distances, clustering, untangling, or Mantel tests.
#
# V11 changes only the presentation:
#   1. Beetle and fungal-community dendrograms are drawn as two
#      opposing semicircular trees.
#   2. Identical beetle samples are connected through the centre.
#   3. Tip labels use exactly the same hierarchy as the fruit trees:
#        Tribe-Species-ID
#        Tribe-Genus-ID
#        Tribe-ID
#        Subfamily-ID
#
# IMPORTANT:
# Original sample IDs remain the internal matching keys. Taxonomic
# labels are display text only, so tree topology, tip correspondence,
# connector colours, and all saved statistics remain unchanged.
# ============================================================


# ============================================================
# 0. INPUT AND OUTPUT PATHS
# ============================================================

data_dir <- Sys.getenv("BEETLE_DATA_DIR", unset = "data")
results_dir <- Sys.getenv("BEETLE_RESULTS_DIR", unset = "results")

analysis_object_file <- file.path(
  results_dir,
  "H3_analysis_objects.rds"
)

default_plot_dir <- file.path(
  results_dir,
  "Figure5_tanglegram"
)

# This optional environment variable is useful for testing a preview
# elsewhere. When it is not set, figures go to default_plot_dir.
plot_dir <- Sys.getenv(
  "TANGLEGRAM_OUTPUT_DIR",
  unset = default_plot_dir
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. FIXED PROJECT COLOURS AND VISUAL SETTINGS
# ============================================================

FINAL_tribe_cols <- c(
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
# 2. LOAD SAVED ANALYSIS OBJECTS
# ============================================================

if (!file.exists(analysis_object_file)) {
  stop(
    paste0(
      "The saved tanglegram analysis file was not found:\n",
      analysis_object_file
    )
  )
}

final_tanglegram_analysis <- readRDS(analysis_object_file)

required_saved_components <- c(
  "Scolytinae",
  "Platypodinae",
  "beetle_metadata"
)

missing_saved_components <- setdiff(
  required_saved_components,
  names(final_tanglegram_analysis)
)

if (length(missing_saved_components) > 0L) {
  stop(
    paste0(
      "Saved analysis file is missing: ",
      paste(missing_saved_components, collapse = ", ")
    )
  )
}

scoly_tangle <- final_tanglegram_analysis$Scolytinae
platy_tangle <- final_tanglegram_analysis$Platypodinae
beetle_metadata <- final_tanglegram_analysis$beetle_metadata


# ============================================================
# 3. TAXONOMIC DISPLAY LABELS
#
# This is deliberately the same rule used by make_tip_text() in
# the final fruit-tree script.
# ============================================================

required_metadata_columns <- c(
  "beetle",
  "Subfamily",
  "Tribe",
  "Genus",
  "Species"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  names(beetle_metadata)
)

if (length(missing_metadata_columns) > 0L) {
  stop(
    paste0(
      "beetle_metadata is missing: ",
      paste(missing_metadata_columns, collapse = ", ")
    )
  )
}

if (anyDuplicated(beetle_metadata$beetle) > 0L) {
  stop("beetle_metadata contains duplicated beetle IDs.")
}

has_text <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

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

beetle_metadata$display_label <- make_tip_text(
  subfamily = beetle_metadata$Subfamily,
  tribe = beetle_metadata$Tribe,
  genus = beetle_metadata$Genus,
  species = beetle_metadata$Species,
  tip_id = beetle_metadata$beetle
)

valid_tribe <- has_text(beetle_metadata$Tribe) &
  beetle_metadata$Tribe != "Unknown"

beetle_metadata$label_rule <- ifelse(
  valid_tribe & has_text(beetle_metadata$Species),
  "Tribe-Species-ID",
  ifelse(
    valid_tribe & has_text(beetle_metadata$Genus),
    "Tribe-Genus-ID",
    ifelse(
      valid_tribe,
      "Tribe-ID",
      ifelse(
        has_text(beetle_metadata$Subfamily),
        "Subfamily-ID",
        "ID only"
      )
    )
  )
)

if (anyDuplicated(beetle_metadata$display_label) > 0L) {
  stop(
    "Taxonomic display labels are not unique. Keep the sample ID in every label."
  )
}

display_lookup <- stats::setNames(
  beetle_metadata$display_label,
  beetle_metadata$beetle
)

tribe_lookup <- stats::setNames(
  as.character(beetle_metadata$Tribe),
  beetle_metadata$beetle
)

get_display_label <- function(ids) {
  ids <- as.character(ids)
  missing_ids <- setdiff(ids, names(display_lookup))
  
  if (length(missing_ids) > 0L) {
    stop(
      paste0(
        "Tip IDs missing from beetle_metadata: ",
        paste(missing_ids, collapse = ", ")
      )
    )
  }
  
  unname(display_lookup[ids])
}

get_tribe <- function(ids) {
  unname(tribe_lookup[as.character(ids)])
}

get_tribe_colour <- function(tribes) {
  colours <- unname(FINAL_tribe_cols[as.character(tribes)])
  colours[is.na(colours)] <- TREE_MIXED_COLOUR
  colours
}

get_branch_colour <- function(descendant_ids) {
  tribes <- unique(get_tribe(descendant_ids))
  tribes <- tribes[has_text(tribes) & tribes != "Unknown"]
  
  if (length(tribes) == 1L) {
    get_tribe_colour(tribes)
  } else {
    TREE_MIXED_COLOUR
  }
}


# ============================================================
# 4. FORMAT MANTEL P
# ============================================================

format_mantel_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.0001) return("< 0.0001")
  if (p < 0.001) return(formatC(p, format = "f", digits = 4))
  formatC(p, format = "f", digits = 3)
}


# ============================================================
# 5. CONVERT A DENDROGRAM TO SEMICIRCULAR COORDINATES
#
# Tips lie on the inner radius. Increasing dendrogram height moves
# branches outwards, keeping the centre free for correspondence lines.
# ============================================================

build_semicircle_layout <- function(
    dend,
    start_angle,
    end_angle,
    tip_radius = TIP_RADIUS,
    root_radius = ROOT_RADIUS
) {
  tip_order <- as.character(labels(dend))
  
  if (length(tip_order) < 2L) {
    stop("A semicircular dendrogram needs at least two tips.")
  }
  
  tip_angles <- seq(
    start_angle,
    end_angle,
    length.out = length(tip_order)
  )
  names(tip_angles) <- tip_order
  
  root_height <- as.numeric(attr(dend, "height"))
  if (!is.finite(root_height) || root_height <= 0) root_height <- 1
  
  store <- new.env(parent = emptyenv())
  store$next_id <- 0L
  store$nodes <- list()
  store$edges <- list()
  
  walk_node <- function(node) {
    store$next_id <- store$next_id + 1L
    node_id <- store$next_id
    
    if (is.leaf(node)) {
      tip_id <- as.character(attr(node, "label"))
      theta <- unname(tip_angles[tip_id])
      
      if (!is.finite(theta)) {
        stop(paste0("Could not assign an angle to tip: ", tip_id))
      }
      
      store$nodes[[node_id]] <- list(
        id = node_id,
        theta = theta,
        radius = tip_radius,
        is_tip = TRUE,
        tip_id = tip_id,
        descendants = tip_id
      )
      
      return(node_id)
    }
    
    child_ids <- vapply(
      seq_along(node),
      function(i) walk_node(node[[i]]),
      integer(1)
    )
    
    descendants <- unlist(
      lapply(
        child_ids,
        function(child_id) store$nodes[[child_id]]$descendants
      ),
      use.names = FALSE
    )
    
    descendant_angles <- unname(tip_angles[descendants])
    theta <- mean(range(descendant_angles))
    
    node_height <- as.numeric(attr(node, "height"))
    if (!is.finite(node_height)) node_height <- 0
    
    radius <- tip_radius +
      (node_height / root_height) * (root_radius - tip_radius)
    
    store$nodes[[node_id]] <- list(
      id = node_id,
      theta = theta,
      radius = radius,
      is_tip = FALSE,
      tip_id = NA_character_,
      descendants = descendants
    )
    
    for (child_id in child_ids) {
      child_descendants <- store$nodes[[child_id]]$descendants
      store$edges[[length(store$edges) + 1L]] <- list(
        parent = node_id,
        child = child_id,
        colour = get_branch_colour(child_descendants)
      )
    }
    
    node_id
  }
  
  root_id <- walk_node(dend)
  
  nodes <- do.call(
    rbind,
    lapply(
      store$nodes,
      function(x) {
        data.frame(
          id = x$id,
          theta = x$theta,
          radius = x$radius,
          is_tip = x$is_tip,
          tip_id = x$tip_id,
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  edges <- do.call(
    rbind,
    lapply(
      store$edges,
      function(x) {
        data.frame(
          parent = x$parent,
          child = x$child,
          colour = x$colour,
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  tips <- nodes[nodes$is_tip, , drop = FALSE]
  tips <- tips[match(tip_order, tips$tip_id), , drop = FALSE]
  rownames(tips) <- NULL
  
  list(
    nodes = nodes,
    edges = edges,
    tips = tips,
    tip_order = tip_order,
    root_id = root_id
  )
}


# ============================================================
# 6. LOW-LEVEL DRAWING HELPERS
# ============================================================

polar_x <- function(radius, theta) radius * cos(theta)
polar_y <- function(radius, theta) radius * sin(theta)

draw_arc <- function(radius, theta_1, theta_2, colour, lwd) {
  n_points <- max(3L, ceiling(abs(theta_2 - theta_1) * 90))
  theta <- seq(theta_1, theta_2, length.out = n_points)
  graphics::lines(
    polar_x(radius, theta),
    polar_y(radius, theta),
    col = colour,
    lwd = lwd
  )
}

draw_semicircle_tree <- function(layout, branch_lwd, tip_cex) {
  nodes <- layout$nodes
  
  for (i in seq_len(nrow(layout$edges))) {
    edge <- layout$edges[i, ]
    parent <- nodes[nodes$id == edge$parent, , drop = FALSE]
    child <- nodes[nodes$id == edge$child, , drop = FALSE]
    
    draw_arc(
      radius = parent$radius,
      theta_1 = parent$theta,
      theta_2 = child$theta,
      colour = edge$colour,
      lwd = branch_lwd
    )
    
    graphics::segments(
      x0 = polar_x(parent$radius, child$theta),
      y0 = polar_y(parent$radius, child$theta),
      x1 = polar_x(child$radius, child$theta),
      y1 = polar_y(child$radius, child$theta),
      col = edge$colour,
      lwd = branch_lwd
    )
  }
  
  tip_colours <- get_tribe_colour(get_tribe(layout$tips$tip_id))
  
  graphics::points(
    polar_x(layout$tips$radius, layout$tips$theta),
    polar_y(layout$tips$radius, layout$tips$theta),
    pch = 21,
    bg = tip_colours,
    col = TIP_BORDER_COLOUR,
    cex = tip_cex,
    lwd = 0.55
  )
}

cubic_bezier <- function(p0, p1, p2, p3, n = 80L) {
  t <- seq(0, 1, length.out = n)
  one_minus_t <- 1 - t
  
  x <- one_minus_t^3 * p0[1] +
    3 * one_minus_t^2 * t * p1[1] +
    3 * one_minus_t * t^2 * p2[1] +
    t^3 * p3[1]
  
  y <- one_minus_t^3 * p0[2] +
    3 * one_minus_t^2 * t * p1[2] +
    3 * one_minus_t * t^2 * p2[2] +
    t^3 * p3[2]
  
  list(x = x, y = y)
}

draw_connectors <- function(
    upper_layout,
    lower_layout,
    label_cex,
    connector_lwd,
    connector_alpha
) {
  if (!setequal(upper_layout$tip_order, lower_layout$tip_order)) {
    stop("The two dendrograms do not contain the same beetle IDs.")
  }
  
  lower_index <- stats::setNames(
    seq_len(nrow(lower_layout$tips)),
    lower_layout$tips$tip_id
  )
  
  # Each line reaches the inner end of its own taxonomic label.
  # Because label lengths differ, the two endpoint radii are
  # calculated separately for every matched beetle sample.
  upper_label_widths <- graphics::strwidth(
    get_display_label(upper_layout$tips$tip_id),
    units = "user",
    cex = label_cex,
    family = "sans"
  )
  
  lower_label_widths <- graphics::strwidth(
    get_display_label(lower_layout$tips$tip_id),
    units = "user",
    cex = label_cex,
    family = "sans"
  )
  
  upper_endpoint_radii <- pmax(
    0.16,
    LABEL_RADIUS - upper_label_widths + 0.002
  )
  
  lower_endpoint_radii <- pmax(
    0.16,
    LABEL_RADIUS - lower_label_widths + 0.002
  )
  
  for (i in seq_len(nrow(upper_layout$tips))) {
    tip_id <- upper_layout$tips$tip_id[i]
    j <- unname(lower_index[tip_id])
    
    theta_upper <- upper_layout$tips$theta[i]
    theta_lower <- lower_layout$tips$theta[j]
    
    p0 <- c(
      polar_x(upper_endpoint_radii[i], theta_upper),
      polar_y(upper_endpoint_radii[i], theta_upper)
    )
    p3 <- c(
      polar_x(lower_endpoint_radii[j], theta_lower),
      polar_y(lower_endpoint_radii[j], theta_lower)
    )
    
    # Control points pull each connector gently towards the centre.
    p1 <- p0 * 0.30
    p2 <- p3 * 0.30
    
    curve <- cubic_bezier(p0, p1, p2, p3)
    connector_colour <- grDevices::adjustcolor(
      get_tribe_colour(get_tribe(tip_id)),
      alpha.f = connector_alpha
    )
    
    graphics::lines(
      curve$x,
      curve$y,
      col = connector_colour,
      lwd = connector_lwd
    )
  }
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

draw_tip_labels <- function(
    layout,
    label_cex
) {
  display_labels <- get_display_label(layout$tips$tip_id)
  
  for (i in seq_len(nrow(layout$tips))) {
    theta <- layout$tips$theta[i]
    
    graphics::text(
      x = polar_x(LABEL_RADIUS, theta),
      y = polar_y(LABEL_RADIUS, theta),
      labels = display_labels[i],
      srt = label_rotation(theta),
      adj = c(label_adjustment(theta), 0.5),
      cex = label_cex,
      col = LABEL_COLOUR,
      family = "sans"
    )
  }
}


# ============================================================
# 7. BOTTOM LEGEND PANEL
# ============================================================

draw_bottom_panel <- function(mantel_text, present_tribes, legend_cex) {
  graphics::par(mar = c(0.2, 1.2, 0.2, 1.2), family = "sans", xpd = NA)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  graphics::text(
    0.5,
    0.82,
    mantel_text,
    cex = legend_cex * 1.06,
    font = 2,
    col = LABEL_COLOUR
  )
  
  graphics::text(
    0.035,
    0.51,
    "Lines",
    adj = c(0, 0.5),
    cex = legend_cex,
    font = 2,
    col = LABEL_COLOUR
  )
  graphics::segments(0.11, 0.51, 0.16, 0.51, col = "#A8A8A8", lwd = 1.2)
  graphics::text(
    0.18,
    0.51,
    "connect identical beetle taxa; colours indicate beetle tribe.",
    adj = c(0, 0.5),
    cex = legend_cex * 0.88,
    col = LABEL_COLOUR
  )
  
  graphics::text(
    0.035,
    0.19,
    "Tribe",
    adj = c(0, 0.5),
    cex = legend_cex,
    font = 2,
    col = LABEL_COLOUR
  )
  
  legend_x <- seq(0.17, 0.94, length.out = length(present_tribes))
  graphics::points(
    legend_x,
    rep(0.19, length(legend_x)),
    pch = 21,
    bg = get_tribe_colour(present_tribes),
    col = TIP_BORDER_COLOUR,
    cex = legend_cex * 0.93
  )
  graphics::text(
    legend_x + 0.012,
    rep(0.19, length(legend_x)),
    labels = present_tribes,
    adj = c(0, 0.5),
    cex = legend_cex * 0.79,
    col = LABEL_COLOUR
  )
}


# ============================================================
# 8. MAIN FIGURE FUNCTION
# ============================================================

plot_semicircle_tanglegram <- function(
    analysis_object,
    subfamily_name,
    panel_label,
    output_dir = plot_dir
) {
  left_dend <- analysis_object$beetle_dendrogram_plot
  right_dend <- analysis_object$fungal_dendrogram_plot
  
  if (is.null(left_dend) || is.null(right_dend)) {
    stop(paste0(subfamily_name, ": dendrogram object missing."))
  }
  
  original_left_ids <- as.character(labels(left_dend))
  original_right_ids <- as.character(labels(right_dend))
  
  if (!setequal(original_left_ids, original_right_ids)) {
    stop(paste0(subfamily_name, ": left and right tip IDs differ."))
  }
  
  # This strict check guarantees every plotted taxonomic label is valid.
  get_display_label(c(original_left_ids, original_right_ids))
  
  mantel_r <- unname(analysis_object$mantel$statistic)
  mantel_p <- analysis_object$mantel$signif
  mantel_text <- paste0(
    "Mantel test: r = ",
    sprintf("%.3f", mantel_r),
    ", P = ",
    format_mantel_p(mantel_p)
  )
  
  if (identical(subfamily_name, "Scolytinae")) {
    label_cex <- 0.285
    tip_cex <- 0.43
    branch_lwd <- 1.02
    connector_lwd <- 1.55
    connector_alpha <- 0.82
    legend_cex <- 0.79
    png_width <- 6200
    png_height <- 6500
    pdf_width <- 19.2
    pdf_height <- 20.0
  } else {
    label_cex <- 0.46
    tip_cex <- 0.62
    branch_lwd <- 1.24
    connector_lwd <- 1.90
    connector_alpha <- 0.88
    legend_cex <- 0.88
    png_width <- 5600
    png_height <- 5900
    pdf_width <- 18.0
    pdf_height <- 18.8
  }
  
  angle_gap <- 3.5 * pi / 180
  
  upper_layout <- build_semicircle_layout(
    left_dend,
    start_angle = pi - angle_gap,
    end_angle = angle_gap
  )
  
  lower_layout <- build_semicircle_layout(
    right_dend,
    start_angle = pi + angle_gap,
    end_angle = 2 * pi - angle_gap
  )
  
  present_tribes <- names(FINAL_tribe_cols)
  present_tribes <- present_tribes[
    present_tribes %in% unique(get_tribe(original_left_ids)) &
      present_tribes != "Unknown"
  ]
  
  draw_complete_figure <- function() {
    graphics::layout(
      matrix(c(1, 2), nrow = 2),
      heights = c(15.2, 2.1)
    )
    
    graphics::par(
      mar = c(0.1, 0.7, 5.5, 0.7),
      family = "sans",
      xpd = NA
    )
    
    graphics::plot.new()
    graphics::plot.window(
      xlim = c(-1.56, 1.56),
      ylim = c(-1.56, 1.56),
      asp = 1
    )
    
    graphics::mtext(
      "Phylogenetic concordance between bark beetles and associated fungal communities",
      side = 3,
      line = 3.6,
      cex = 1.44,
      font = 2,
      col = LABEL_COLOUR
    )
    graphics::mtext(
      paste0(
        "The upper semicircle shows beetle phylogeny and the lower semicircle shows ",
        "Jaccard-based fungal community similarity; curves link identical beetle taxa."
      ),
      side = 3,
      line = 2.0,
      cex = 0.83,
      col = LABEL_COLOUR
    )
    graphics::mtext(
      panel_label,
      side = 3,
      line = 0.25,
      adj = 0.02,
      cex = 1.28,
      font = 2,
      col = LABEL_COLOUR
    )
    graphics::mtext(
      subfamily_name,
      side = 3,
      line = 0.25,
      adj = 0.5,
      cex = 1.10,
      font = 2,
      col = LABEL_COLOUR
    )
    
    graphics::text(
      0,
      1.455,
      "Beetle phylogenetic relationships",
      cex = 0.78,
      font = 2,
      col = LABEL_COLOUR
    )
    graphics::text(
      0,
      -1.455,
      "Fungal community similarity",
      cex = 0.78,
      font = 2,
      col = LABEL_COLOUR
    )
    
    # Central links are drawn first so labels and tree branches remain crisp.
    draw_connectors(
      upper_layout,
      lower_layout,
      label_cex = label_cex,
      connector_lwd = connector_lwd,
      connector_alpha = connector_alpha
    )
    
    draw_tip_labels(
      upper_layout,
      label_cex = label_cex
    )
    draw_tip_labels(
      lower_layout,
      label_cex = label_cex
    )
    
    draw_semicircle_tree(
      upper_layout,
      branch_lwd = branch_lwd,
      tip_cex = tip_cex
    )
    draw_semicircle_tree(
      lower_layout,
      branch_lwd = branch_lwd,
      tip_cex = tip_cex
    )
    
    draw_bottom_panel(
      mantel_text = mantel_text,
      present_tribes = present_tribes,
      legend_cex = legend_cex
    )
  }
  
  png_file <- file.path(
    output_dir,
    paste0(subfamily_name, "_tanglegram_THESIS_v11_semicircle.png")
  )
  pdf_file <- file.path(
    output_dir,
    paste0(subfamily_name, "_tanglegram_THESIS_v11_semicircle.pdf")
  )
  
  grDevices::png(
    filename = png_file,
    width = png_width,
    height = png_height,
    res = 300,
    bg = "white"
  )
  draw_complete_figure()
  grDevices::dev.off()
  
  grDevices::pdf(
    file = pdf_file,
    width = pdf_width,
    height = pdf_height,
    family = "Helvetica",
    useDingbats = FALSE
  )
  draw_complete_figure()
  grDevices::dev.off()
  
  data.frame(
    Subfamily = subfamily_name,
    n_beetles = length(original_left_ids),
    Mantel_r = mantel_r,
    Mantel_P = mantel_p,
    PNG = png_file,
    PDF = pdf_file,
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 9. DRAW BOTH TANGLEGRAMS
# ============================================================

scoly_thesis_v11 <- plot_semicircle_tanglegram(
  analysis_object = scoly_tangle,
  subfamily_name = "Scolytinae",
  panel_label = "A"
)

platy_thesis_v11 <- plot_semicircle_tanglegram(
  analysis_object = platy_tangle,
  subfamily_name = "Platypodinae",
  panel_label = "B"
)

thesis_v11_summary <- rbind(
  scoly_thesis_v11,
  platy_thesis_v11
)

utils::write.csv(
  thesis_v11_summary,
  file.path(plot_dir, "Tanglegram_THESIS_v11_semicircle_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  beetle_metadata[
    ,
    c(
      "beetle",
      "Subfamily",
      "Tribe",
      "Genus",
      "Species",
      "label_rule",
      "display_label"
    )
  ],
  file.path(plot_dir, "Tanglegram_taxonomic_tip_label_key.csv"),
  row.names = FALSE
)

cat("\nTHESIS TANGLEGRAM V11 COMPLETE\n")
cat("Saved to:\n", plot_dir, "\n\n", sep = "")
print(thesis_v11_summary, row.names = FALSE)
cat("\nTaxonomic label-rule counts:\n")
print(table(beetle_metadata$label_rule))
