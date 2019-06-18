#' BioCircos inspired plot of genomes in linear representation
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
#' @param track_margin (numeric) (optional) the margin added to each side of the genome, a scalar in bp.
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
BioCircos_linear <- function(
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
  gene_linetrack = NULL,
  range_abundance = NULL, 
  range_foldchange = NULL,
  range_linetrack = NULL,
  track_margin = 500
) {
  
  # DEFINE GLOABL VARIABLES
  # --------------------------------
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
  
  # set xlims for genome
  genome_xlim <- range(data[c(gene_start, gene_end)]) + 
    c(-track_margin, track_margin)
  
  # set colors for genome
  if (is.null(genome_color)) {
    genome_color <- qualitative_hcl(n = length(genome), h = c(20, 360), c = 100, l = 70)
  } else if (length(genome_color) != length(genome) | !is.character(genome_color)) {
    stop("genome_color must be a character vector of length(genome)")
  }
  
  # set strand to only "+" if not provided
  if (!is.null(gene_strand)) {
    if (!(unique(data[[gene_strand]]) %in% c("+", "-")) %>% any) {
      stop("'gene_strand' index must contain only '+' or '-'")
    }
  } else {
    data <- mutate(data, gene_strand = "+")
  }
  
  # set ranges from data if not provided
  for (dat_col in c("abundance", "foldchange", "linetrack")) {
    
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
    data <- dplyr::mutate(data, gene_color = grey(0.8))
  } else {
    custom_palette =  colorspace::qualitative_hcl(
      n = length(unique(data[[gene_category]])), 
      h = c(20, 360), c = 100, l = 50)
    data$gene_color <- data[[gene_category]] %>% 
      as.factor %>% as.numeric %>% custom_palette[.]
  }
  
  
  # LAYER 1: GENE POSITIONS
  # --------------------------------
  tracklist <- list()
  tracklist$genes <- xyplot(get(gene_end) ~ get(gene_start), data,
    xlim = genome_xlim,
    ylim = c(-4, 2), 
    xlab = "location [bp]", ylab = "",
    par.settings = custom.lattice,
    scales = list(alternating = FALSE, y = list(rot = 0)),
    panel = function(x, y, ...) {
      strand = ifelse(data[[gene_strand]] == "+", 1, -1)
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.rect(xleft = x, ybottom = 0, xright = y, ytop = strand,
        border = grey(0.3), fill = data[["gene_color"]], lwd = 1.5)
      panel.text(x+(y-x)/2, -1.5, data[[gene_name]], col = grey(0.3), 
        srt = 35, cex = 0.7, adj = c(1.0, 1.0), pos = NULL)
      panel.segments(x0 = x+(y-x)/2, y0 = ifelse(strand == -1, -1, 0), x1 = x+(y-x)/2, y1 = -1.4, 
        col = grey(0.3), lwd = 0.8)
      panel.abline(h = 0, lwd = 1.5, col = grey(0.3))
    }
  )
  
  
  # LAYER 2: ABUNDANCES
  # --------------------------------
  if (!is.null(gene_abundance)) {
    
    form <- paste0(paste(gene_abundance, collapse = " + "),
      " ~ (get(gene_start) + (get(gene_end) - get(gene_start))/2)")
    
    tracklist$abundance <- xyplot(as.formula(form), data,
      par.settings = custom.lattice, 
      xlim = genome_xlim, ylim = range_abundance,
      panel = function(x, y, ...) {
        panel.grid(h = -1, v = -1, col = grey(0.9))
        panel.xyarea(x, y, lwd = 1.5, alpha = 0.5, origin = 0, ...)
      }
    )
    
  }
  
  
  # return trellis plot without printing
  tracklist <- unname(tracklist)
  tracklist$layout = c(1, length(tracklist))
  do.call(c, args = tracklist)
  
}
