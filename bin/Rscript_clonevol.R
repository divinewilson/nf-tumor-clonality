#!/usr/bin/env Rscript

library(clonevol)
library(ggplot2)

ClonEvol_function <- function(tree_file, variant_file, outpdf, plot_w, plot_h) {

  
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

  ## Format gene labels (Split into Left and Right columns)
  split_genes <- function(gene_string, side = "left") {
    genes <- trimws(strsplit(as.character(gene_string), ",")[[1]])
    if (length(genes) == 0 || gene_string == "") return("")
    
    n <- length(genes)
    mid <- ceiling(n / 2)
    
    if (side == "left") {
      selected <- genes[1:mid]
    } else {
      # If only 1 gene, put it on left, return empty for right
      if (n == 1) return("") 
      selected <- genes[(mid + 1):n]
    }
    
    # Stack them vertically (1 gene per line)
    # For very long lists, you can do 2 per line by changing this:
    paste(selected, collapse = "\n")
  }
  germ_data$text_left  <- sapply(germ_data$branch.texts, split_genes, side = "left")
  germ_data$text_right <- sapply(germ_data$branch.texts, split_genes, side = "right")
  
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

    ## Gene labels: LEFT SIDE
    geom_text(
      aes(
        x = (x0 + x1) / 2,
        y = (y0 + y1) / 2,
        label = text_left
      ),
      size  = 2,
      hjust = 1.1,      # Right-align the text against the branch
      vjust = 0.5,
      nudge_x = -0.3,   # Push slightly to the left
      lineheight = 0.8  # Tighten vertical space between genes
    ) +

    ## Gene labels: RIGHT SIDE
    geom_text(
      aes(
        x = (x0 + x1) / 2,
        y = (y0 + y1) / 2,
        label = text_right
      ),
      size  = 2,
      hjust = -0.1,     # Left-align the text against the branch
      vjust = 0.5,
      nudge_x = 0.3,    # Push slightly to the right
      lineheight = 0.8
    ) +

    labs(title = sample_info) +
    theme_void() +
    theme(plot.title = element_text(size = 8)) +

    expand_limits(
      x = range(germ_data$x0, germ_data$x1) + c(-2, 2),
      y = range(germ_data$y0, germ_data$y1) + c(-2, 2)
    )
  ggsave(outpdf, plot = p, width = plot_w, height = plot_h, dpi = 400)
    
    
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

  ggsave(outpdf_clean, plot = q, width = 2.5, height = 5.5, dpi = 400)

}



tree_file    <- as.character(commandArgs(trailingOnly = TRUE)[1])
variant_file <- as.character(commandArgs(trailingOnly = TRUE)[2])
outpdf       <- as.character(commandArgs(trailingOnly = TRUE)[3])
plot_w       <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
plot_h       <- as.numeric(commandArgs(trailingOnly = TRUE)[5])

ClonEvol_function(tree_file, variant_file, outpdf, plot_w, plot_h)
