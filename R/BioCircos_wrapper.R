#' Convenience wrapper function for the BioCircos package
#'
#' This function takes as arguments a data.frame and a list of columns
#'
#' @param data A data frame or tibble.
#' @param genome A named list of chromosomes and their total length such as
#'   \code{list(chromsome_1 = 1000)}. If not provided the genome information
#'   is retrieved from gene annotation.
#' @param genome_color (character) (optional) color vector of same length as \code{genome}.
#' @param gene_name (character) column name of the data frame containing gene names.
#' @param gene_start (character) column name of the data frame containing start positions in bp.
#' @param gene_end (character) column name of the data frame containing end positions in pb.
#' @param gene_chromosome (character) column name of the data frame containing the chromosome ID of each gene.
#' @param gene_strand (character) (optional) column name of the data frame containing the strand of each gene ('+' or '-').
#' @param gene_category (character) (optional) column name of the data frame containing gene annotation such 
#'   as 'Respiratory chain' or 'Unknown'.
#' @param gene_abundance (character) (optional) column name of the data frame containing an absolute abundance 
#'   such as protein length or mass fraction.
#' @param gene_foldchange (character) (optional) column name of the data frame containing fold change values for a gene.
#' @param gene_heatmap (character) (optional) column name(s) of the data frame containing abundances plotted as heatmap. 
#'   Can contain more than one element.
#' @param gene_linetrack (character) (optional) column name(s) of the data frame containing abundances plotted as lines.
#'   Can contain more than one element.
#' @param gene_labels (character) (optional) column name with labels that will be plotted as arcs. 
#'   Can contain NAs for labels hat should not be plotted.
#' @param color_heatmap (character) (optional) vector of length 2 indicating colors to use for heatmap gradient
#' @param range_abundance (numeric) (optional) range for plotting gene/protein abundances. Vector of length 2.
#' @param range_foldchange (numeric) (optional) range for plotting gene/protein fold changes. Vector of length 2.
#' @param range_heatmap (numeric) (optional) range for plotting gene/protein abundances. Vector of length 2.
#' @param range_linetrack (numeric) (optional) range for plotting gene/protein abundances. Vector of length 2.
#' @param track_size (numeric) (optional) size of one circular layer added to the plot, a scalar between 0 and 1.
#' @param track_gap (numeric) (optional) size of the gap between circular layers, a scalar between 0 and 1.
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_at
#' @importFrom dplyr summarise
#' @importFrom scales rescale
#' @importFrom colorspace qualitative_hcl
#'
#' @export
BioCircos_wrapper <- function(
  data, 
  genome = NULL,
  genome_color = NULL,
  gene_name,
  gene_start,
  gene_end,
  gene_chromosome,
  gene_strand = NULL,
  gene_category = NULL,
  gene_abundance = NULL,
  gene_foldchange = NULL,
  gene_heatmap = NULL,
  gene_linetrack = NULL,
  gene_labels = NULL,
  color_heatmap = c("#F5F5F5", "#4E962D"),
  range_abundance = NULL, 
  range_foldchange = NULL,
  range_heatmap = NULL, 
  range_linetrack = NULL,
  track_size = 0.15,
  track_gap = 0.02
) {
  
  
  # DEFINE GLOBAL VARIABLES ----------------------------------------------------
  #
  # re-order data frame by chromosome and gene position
  data <- dplyr::arrange(data, get(gene_chromosome), get(gene_start))
  
  # guess genome information from gene list in case it's not provided
  if (is.null(genome)) {
    genome = group_by(data, get(gene_chromosome)) %>% 
      dplyr::summarise(
        chromosome = unique(get(gene_chromosome)), 
        length = max(get(gene_end))
      )
    genome = genome %>% pull(length) %>% as.list %>% 
      setNames(genome$chromosome)
  }
  
  
  # set colors for genome circles
  if (is.null(genome_color)) {
    genome_color <- qualitative_hcl(n = length(genome), h = c(20, 360), c = 100, l = 70)
  } else if (length(genome_color) != length(genome) | !is.character(genome_color)) {
    stop("genome_color must be a character vector of length(genome)")
  }
  
  
  # set ranges from data if not provided
  for (dat_col in c("abundance", "foldchange", "heatmap", "linetrack")) {
    
    gene_var <- paste0("gene_", dat_col)
    range_var <- paste0("range_", dat_col)
    
    if (!is.null(get(gene_var))) {
      # clip abundances to the provided range
      if (!is.null(get(range_var))) {
        data <- dplyr::mutate_at(data, vars(get(gene_var)), 
          function(x) pmax(get(range_var)[1], pmin(x, get(range_var)[2])))
      } else {
      # alternatively set range 
        assign(range_var, c(
          min(data[get(gene_var)], na.rm = TRUE), 
          max(data[get(gene_var)], na.rm = TRUE)
        ))
      }
    }
    
  }
  
  # define color code for genes
  if (is.null(gene_category)) {
    data <- dplyr::mutate(data, gene_color = "#56D69A")
  } else {
    custom_palette =  colorspace::qualitative_hcl(
      n = length(unique(data[[gene_category]])), 
      h = c(20, 360), c = 100, l = 50)
    data$gene_color <- data[[gene_category]] %>% 
      as.factor %>% as.numeric %>% custom_palette[.]
  }
  
  
  # PLOTTING OF DIFFERENT LAYERS -----------------------------------------------
  #
  # define starting radii for inner and outer ring dimension
  maxRad <- 1-track_gap; minRad <- maxRad-track_size
  
  
  # ADD GENES AS ARCS
  if (!is.null(gene_strand)) {
    if (unique(data[[gene_strand]]) %in% c("+", "-") %>% all) {
      data_plus <- dplyr::filter(data, get(gene_strand) == "+")
      data_minus <- dplyr::filter(data, get(gene_strand) == "-")
    } else {
      stop("'gene_strand' index must contain only '+' or '-'")
    }
    
    # need two individual tracks for plus and minus strand
    tracklist = BioCircosArcTrack("track_genes_plus", 
      labels = data_plus[[gene_name]],
      colors = data_plus$gene_color,
      chromosomes = data_plus[[gene_chromosome]],
      starts = data_plus[[gene_start]],
      ends = data_plus[[gene_end]],
      minRadius = 1.2, maxRadius = 1.22)
  
    tracklist = tracklist + BioCircosArcTrack("track_genes_minus", 
      labels = data_minus[[gene_name]],
      colors = data_minus$gene_color,
      chromosomes = data_minus[[gene_chromosome]],
      starts = data_minus[[gene_start]],
      ends = data_minus[[gene_end]],
      minRadius = 1.2, maxRadius = 1.18)
    
  } else {
    
    # or just one plus strand
    tracklist = BioCircosArcTrack("track_genes_plus", 
      labels = data[[gene_name]],
      colors = data$gene_color,
      chromosomes = data[[gene_chromosome]],
      starts = data[[gene_start]],
      ends = data[[gene_end]],
      minRadius = 1.2, maxRadius = 1.22)
  }
  
  
  # ADD GENES AS ARCS
  if (!is.null(gene_labels)) {
    
    # generate label arcs from column
    d <- dplyr::filter(data, !is.na(get(gene_labels))) %>%
      dplyr::group_by(get(gene_labels)) %>%
      dplyr::summarise(
        start = get(gene_start)[1],
        end = tail(get(gene_end), 1),
        chromosome = get(gene_chromosome)[1]
      )

    tracklist = tracklist + BioCircosArcTrack("track_gene_labels", 
      labels = d[['get(gene_labels)']],
      colors = grey(0.6),
      chromosomes = d[["chromosome"]],
      starts = d[["start"]],
      ends = d[["end"]],
      minRadius = 1.28, maxRadius = 1.30)
    
  }
  
  
  # ADD BARS FOR PROTEIN ABUNDANCES
  if (!is.null(gene_abundance)) {
    
    # we have to loop through all chromosomes
    for (chr in unique(data[[gene_chromosome]])) {
      d <- dplyr::filter(data, get(gene_chromosome) == chr, !is.na(get(gene_abundance)))
      tracklist = tracklist + BioCircosBarTrack("track_abundance_bars", 
        chromosomes = chr,
        starts = d[[gene_start]],
        ends = d[[gene_end]],
        values = d[[gene_abundance]],
        color = grey(0.7),
        range = range_abundance,
        maxRadius = maxRad, minRadius = minRad
      )
    }
    
    # abundance bar background
    tracklist = tracklist + BioCircosBackgroundTrack(
      "background", fillColors = grey(0.9), borderSize = 0,
      maxRadius = maxRad, minRadius = minRad)
    
    # adjust radii for next ring
    maxRad <- minRad-track_gap; minRad <- maxRad-track_size
  }
  
  # ADD BARS FOR PROTEIN FOLD CHANGES
  if (!is.null(gene_foldchange)) {
    
    # we have to loop through all chromosomes
    for (chr in unique(data[[gene_chromosome]])) {
      d <- dplyr::filter(data, get(gene_chromosome) == chr, !is.na(get(gene_foldchange)))
      tracklist = tracklist + BioCircosBarTrack("track_fold_change_bars",
        chromosomes = chr,
        starts = d[[gene_start]],
        ends = d[[gene_end]],
        values = d[[gene_foldchange]],
        color = grey(0.7),
        range = c(0, -range_foldchange[1]),
        maxRadius = maxRad, minRadius = minRad
      )
    }
    
    # log FC bar backgrounds
    tracklist = tracklist + BioCircosBackgroundTrack(
      "background", fillColors = "#FFB8B8", borderSize = 0,
      maxRadius = maxRad, minRadius = minRad)
    
    tracklist = tracklist + BioCircosBackgroundTrack(
      "background", fillColors = "#AFD9E0", borderSize = 0,
      maxRadius = minRad, minRadius = minRad-track_size)
    
    # adjust radii for next ring
    maxRad <- minRad-track_size-track_gap; minRad <- maxRad-track_size
    
  }
  
  # ADD HEATMAP (MULTIPLE COLUMNS POSSIBLE)
  if (!is.null(gene_heatmap)) {
    
    # loop through a set of different conditions and make a track
    # with values for each of those
    for (cond in gene_heatmap) {
      
      tracklist = tracklist + BioCircosHeatmapTrack("track_heatmap",
        chromosomes = data[[gene_chromosome]],
        starts = data[[gene_start]],
        ends = data[[gene_end]],
        values = data[[cond]],
        color = color_heatmap,
        range = range_heatmap,
        maxRadius = maxRad, minRadius = maxRad-track_size/length(gene_heatmap))
      
      # adjust radii for next ring
      maxRad <- maxRad-track_size/length(gene_heatmap)
      minRad <- maxRad-track_size/length(gene_heatmap)
      
    }
    
    # adjust radii for next ring
    maxRad <- minRad-track_gap; minRad <- maxRad-track_size
    
  }
  
  # ADD LINE TRACKS (MULTIPLE COLUMNS POSSIBLE)
  if (!is.null(gene_linetrack)) {
    
    # loop through a set of different conditions and make a track
    # with values for each of those
    color_lines <- qualitative_hcl(n = length(gene_linetrack), 
      h = c(20, 360), c = 100, l = 70) %>% setNames(gene_linetrack)
    
    for (cond in gene_linetrack) {
      
      d <- dplyr::filter(data, !is.na(get(cond)))
      tracklist = tracklist + BioCircosLineTrack("track_lines",
        chromosomes = d[[gene_chromosome]],
        positions = d[[gene_start]] + (d[[gene_end]]-d[[gene_start]])/2,
        values = d[[cond]],
        color = color_lines[[cond]],
        range = range_linetrack,
        maxRadius = maxRad, minRadius = minRad
      )
      
    }
    
    # abundance bar background
    tracklist = tracklist + BioCircosBackgroundTrack(
      "background", fillColors = grey(0.9), borderSize = 0,
      maxRadius = maxRad, minRadius = minRad)
    
    # adjust radii for next ring
    maxRad <- minRad-track_gap; minRad <- maxRad-track_size
    
  }
  
  # PLOT COMPLETE TRACKLIST
  BioCircos(tracklist = tracklist, width = "100%", height = "1000px",
    genome = genome, displayGenomeBorder = FALSE, 
    genomeFillColor = genome_color,
    genomeTicksScale = 10000, genomeTicksLen = 3,
    genomeTicksTextSize = "0.2em"
  )
  
}

