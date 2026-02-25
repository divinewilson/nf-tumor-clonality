#!/usr/bin/env Rscript

library(clonevol)
library(ggplot2)

ClonEvol_function <- function(tree_file, variant_file, outpdf) {

  
  ## Load and prepare CloneEvol tree
  y <- import.tree(tree_file, variant_file)

  y <- convert.consensus.tree.clone.to.branch(
    y,
    branch.scale = "sqrt"
  )

  y <- transfer.events.to.consensus.trees(
    y,
    y$variants[y$variants$is.driver, ],
    cluster.col.name = "cluster",
    event.col.name   = "gene"
  )

  
  ## Extract first merged tree
  mt <- y$matched$merged.trees[[1]]

  
  ## Generate germination data (coordinates)
  germ_data <- clonevol:::germinate(
    list(
      trunk.height            = 32,
      branches                 = mt$branches,
      lengths                  = mt$blengths,
      branch.colors            = mt$color,
      branch.border.colors     = mt$branch.border.color,
      branch.border.linetypes  = mt$branch.border.linetype,
      branch.border.widths     = mt$branch.border.width,
      node.colors              = mt$color,
      node.border.colors       = mt$node.border.color,
      node.border.widths       = mt$node.border.width,
      node.labels              = mt$lab,
      node.texts               = mt$samples.with.nonzero.cell.frac,
      branch.texts             = mt$event
    ),
    angle = 30,
    plot  = FALSE
  )

  
  ## Rotate tree 180 degrees
  germ_data$x0 <- -germ_data$x0
  germ_data$x1 <- -germ_data$x1
  germ_data$y0 <- -germ_data$y0
  germ_data$y1 <- -germ_data$y1

  
  ## Node and branch sizes
  germ_data$node.size <- 8

  branch_size_map <- setNames(mt$branch_size, mt$lab)
  germ_data$branch.size <- branch_size_map[
    as.character(germ_data$node.labels)
  ]

  
  ## Format gene labels (horizontal stacking)
  germ_data$branch.texts <- sapply(
  strsplit(as.character(germ_data$branch.texts), ","),
  function(genes) {

    genes <- trimws(genes)           # remove extra spaces
    if (length(genes) == 0) return("")

    # Split into groups of 5
    gene_groups <- split(
      genes,
      ceiling(seq_along(genes) / 5)
    )

    # Join genes within each group by ", "
    lines <- sapply(gene_groups, function(g) {
      paste(g, collapse = ", ")
    })

    # Join groups by newline
    paste(lines, collapse = "\n")
  }
  )

  ## Adjust branch length scale
  germ_data$length <- germ_data$length / 3

  ## Sample label (title)
  sample_info <- paste(
    strsplit(as.character(mt$sample_info[1]), ", ")[[1]],
    collapse = "\n"
  )

  
  ## Plot: tree with gene labels
  p <- ggplot(
    germ_data,
    aes(x = x0, y = y0, xend = x1, yend = y1)
  ) +
    geom_segment(
      aes(color = branch.colors, linewidth = branch.size),
      lineend = "round",
      show.legend = FALSE
    ) +
    scale_color_identity() +

    geom_point(
      aes(x = x1, y = y1, fill = node.colors, size = node.size),
      shape = 21,
      color = "black",
      stroke = 1.2,
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    scale_size_identity() +

    ## Clone number
    geom_text(
      aes(x = x1, y = y1, label = node.labels),
      size  = 3,
      vjust = 0.5,
      hjust = 0.5
    ) +

    ## Sample labels (S1, S2, ...)
    geom_text(
      aes(x = x1, y = y1, label = node.texts),
      size  = 3,
      vjust = 3,
      hjust = 0.5
    ) +

    ## Gene labels on branches
    geom_text(
      aes(
        x = (x0 + x1) / 2,
        y = (y0 + y1) / 2,
        label = branch.texts
      ),
      size  = 2,
      hjust = 0.5,
      vjust = 0.5
    ) +

    labs(title = sample_info) +
    theme_void() +
    theme(plot.title = element_text(size = 8)) +

    expand_limits(
      x = range(germ_data$x0, germ_data$x1) + c(-2, 2),
      y = range(germ_data$y0, germ_data$y1) + c(-2, 2)
    )
  ggsave(outpdf, plot = p, width = 5, height = 10, dpi = 300)
    
    
  ### plot tree without labels
  outpdf_clean <- sub("\\.pdf$", "_clean.pdf", outpdf)
  q <- ggplot(germ_data, aes(x = x0, y = y0, xend = x1, yend = y1)) +
  geom_segment(
    aes(color = branch.colors, linewidth = branch.size),
    lineend = "round",
    show.legend = FALSE
  ) +
  scale_color_identity() +

  geom_point(
    aes(x = x1, y = y1, fill = node.colors, size = node.size),
    shape = 21,
    color = "black",
    stroke = 1.2,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  scale_size_identity() +

  geom_text(
    aes(x = x1, y = y1, label = node.labels),
    size = 3
  ) +

  geom_text(
    aes(x = x1, y = y1, label = node.texts),
    size = 3,
    vjust = 3
  ) +

  labs(title = sample_info) +
  theme_void() +
  theme(plot.title = element_text(size = 8)) +
  expand_limits(
    x = range(germ_data$x0, germ_data$x1) + c(-2, 2),
    y = range(germ_data$y0, germ_data$y1) + c(-2, 2)
  )

  ggsave(outpdf_clean, plot = q, width = 2.5, height = 5.5, dpi = 300)

}



tree_file    <- as.character(commandArgs(trailingOnly = TRUE)[1])
variant_file <- as.character(commandArgs(trailingOnly = TRUE)[2])
outpdf       <- as.character(commandArgs(trailingOnly = TRUE)[3])

ClonEvol_function(tree_file, variant_file, outpdf)
