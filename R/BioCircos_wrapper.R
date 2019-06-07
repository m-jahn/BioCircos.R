#' Convenience wrapper function for the BioCircos package
#'
#' This function takes as arguments a data.frame and a list of columns
#'
#' @param data A data frame or tibble.
#' @param genome A named list of chromosomes and their total length such as
#'   \code{list(chromsome_1 = 1000)}. If not provided the genome information
#'   is retrieved from gene annotation.
#' @param gene_name (character) column name of the data frame containing gene names.
#' @param gene_start (character) column name of the data frame containing start positions in bp.
#' @param gene_end (character) column name of the data frame containing end positions in pb.
#' @param gene_chromosome (character) column name of the data frame containing the chromosome ID of each gene.
#' @param gene_strand (character) (optional) column name of the data frame containing the strand of each gene ('+' or '-').
#' @param gene_category (character) (optional) column name of the data frame containing gene annotation such as 'Respiratory chain' or 'Unknown'.
#' @param gene_abundance (character) (optional) column name of the data frame containing an absolute abundance such as protein length or mass fraction.
#' @param gene_foldchange (character) (optional) column name of the data frame containing fold change values for a gene.
#' @param gene_heatmap (character) (optional) column name(s) of the data frame containing abundances plotted as heatmap.
#' @param color_heatmap (character) (optional) vector of length 2 indicating colors to use for heatmap gradient
#' @param range_abundance (numeric) (optional) range for plotting gene/protein abundances. Vector of length 2.
#' @param range_foldchange (numeric) (optional) range for plotting gene/protein fold changes. Vector of length 2.
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise
#' @importFrom scales rescale
#' @importFrom colorspace qualitative_hcl
#'
#' @export
BioCircos_wrapper <- function(data, 
  genome = NULL, 
  gene_name,
  gene_start,
  gene_end,
  gene_chromosome,
  gene_strand = NULL,
  gene_category = NULL,
  gene_abundance = NULL,
  gene_foldchange = NULL,
  gene_heatmap = NULL,
  color_heatmap = c("#F5F5F5", "#4E962D"),
  range_abundance = c(
    min(data[[gene_abundance]], na.rm = TRUE), 
    max(data[[gene_abundance]], na.rm = TRUE)),
  range_foldchange = c(
    min(data[[gene_foldchange]], na.rm = TRUE),
    max(data[[gene_foldchange]], na.rm = TRUE)),
  track_size = 0.15,
  track_gap = 0.02
) {
  
  
  # DEFINE GLOBAL VARIABLES ----------------------------------------------------
  #
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
  
  # re-order data frame by chromosome and gene position
  data <- dplyr::arrange(data, get(gene_chromosome), get(gene_start))
  
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
  maxRad = 1-track_gap; minRad = maxRad-track_size
  
  
  # add genes as arcs
  if (!is.null(gene_strand)) {
    if (unique(data[[gene_strand]]) %in% c("+", "-") %>% all) {
      data_plus <- filter(data, get(gene_strand) == "+")
      data_minus <- filter(data, get(gene_strand) == "-")
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
  
  # add bars for absolute protein abundance
  if (!is.null(gene_abundance)) {
    
    tracklist = tracklist + BioCircosBarTrack("track_abundance_bars", 
      chromosome = unique(data[[gene_chromosome]]),
      starts = data[[gene_start]],
      ends = data[[gene_end]],
      values = data[[gene_abundance]],
      color = grey(0.7),
      range = range_abundance,
      maxRadius = maxRad, minRadius = minRad
    )
    
    # abundance bar background
    tracklist = tracklist + BioCircosBackgroundTrack(
      "background", fillColors = grey(0.9), borderSize = 0,
      maxRadius = maxRad, minRadius = minRad)
    
    # adjust radii for next ring
    maxRad = minRad-track_gap; minRad = maxRad-track_size
  }
  
  if (!is.null(gene_foldchange)) {
    
    tracklist = tracklist + BioCircosBarTrack("track_fold_change_bars",
      chromosome = unique(data[[gene_chromosome]]),
      starts = data[[gene_start]],
      ends = data[[gene_end]],
      values = data[[gene_foldchange]],
      color = grey(0.7),
      range = c(0, -range_foldchange[1]),
      maxRadius = maxRad, minRadius = minRad)

    # log FC bar backgrounds
    tracklist = tracklist + BioCircosBackgroundTrack(
      "background", fillColors = "#FFB8B8", borderSize = 0,
      maxRadius = maxRad, minRadius = minRad)

    tracklist = tracklist + BioCircosBackgroundTrack(
      "background", fillColors = "#AFD9E0", borderSize = 0,
      maxRadius = minRad, minRadius = minRad-track_size)
    
    # adjust radii for next ring
    maxRad = minRad-track_size-track_gap; minRad = maxRad-track_size
    
  }
  
  if (!is.null(gene_heatmap)) {
    
    # loop through a set of different conditions and make a track
    # with log FC values for each of those
    
    for (cond in gene_heatmap) {
      
      tracklist = tracklist + BioCircosHeatmapTrack("track_heatmap",
        chromosome = data[[gene_chromosome]],
        starts = data[[gene_start]],
        ends = data[[gene_end]],
        values = data[[cond]],
        color = color_heatmap,
        range = range_foldchange,
        maxRadius = maxRad, minRadius = maxRad-track_size/length(gene_heatmap))
      
      # adjust radii for next ring
      maxRad = maxRad-track_size/length(gene_heatmap)
      minRad = maxRad-track_size/length(gene_heatmap)
      
    }
    
    # adjust radii for next ring
    maxRad = minRad-track_gap; minRad = maxRad-track_size
    
  }
    
  # plot complete circos tracklist
  BioCircos(tracklist = tracklist, width = "100%", height = "1000px",
    genome = genome, genomeFillColor = qualitative_hcl(n = length(genome), 
      h = c(20, 360), c = 100, l = 70), displayGenomeBorder = FALSE, 
    genomeTicksScale = 10000, genomeTicksLen = 3,
    genomeTicksTextSize = "0.2em"
  )
  
}

