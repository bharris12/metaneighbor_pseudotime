





`%>%` <- magrittr::`%>%`
# log<-file('log/test_log.log', open='wt')
# log <- file(snakemake@log[[1]], open="wt")
# sink(log, type="message")

parser <-
  argparser::arg_parser(description = "Arguments for Heatmap generation")
parser <-
  argparser::add_argument(parser, "--fn", help = "Filename with MetaNeighbor Output")
parser <-
  argparser::add_argument(parser, "--outdir", help = "Directory to store figure files in")
parser <-
  argparser::add_argument(parser,
                          "--split-lineages",
                          default = 0,
                          help = "Argument to split heatmaps into each lineage, pass non-zero for T")
parser <-
  argparser::add_argument(parser, '--file-ext', default = 'png', help = 'Extensio to use for saving, png or pdf')

args <- argparser::parse_args(parser)
fn <- args$fn
outdir <- args$outdir
split_lineages <- as.logical(args$split_lineages)
extension <- args$file_ext

if (!extension %in% c('pdf', 'png')) {
  message('File Extension Passed was neither pdf or png, setting extension to png')
  extension <- 'png'
}
build_annotation <- function(meta, side = "column") {
  cols <- colnames(meta)
  color_list <- vector(mode = "list", length = length(cols))
  if ("Lineage" %in% cols) {
    lins <- unique(meta$Lineage)
    pal <- RColorBrewer::brewer.pal(length(lins), "Dark2")
    if (length(pal) > length(lins)) {
      # RColorBrewer won't return less than 3 colors
      pal <- pal[1:length(lins)]
    }
    names(pal) <- lins
    color_list$Lineage <- pal
  }
  if ("Study" %in% cols) {
    studies <- unique(meta$Study)
    pal <- RColorBrewer::brewer.pal(length(studies), "Set3")
    names(pal) <- studies
    color_list$Study <- setNames(pal, studies)
  }
  if ("Pseudotime.Bin" %in% cols) {
    n <- nrow(meta)
    min_v <- min(meta[, "Pseudotime.Bin"])
    max_v <- max(meta[, "Pseudotime.Bin"])
    Var <- circlize::colorRamp2(seq(min_v, max_v, length = n),
                                rev(hcl.colors(n, "Blues")))
    
    color_list$Pseudotime.Bin <- Var
  }
  annot <-
    ComplexHeatmap::HeatmapAnnotation(
      df = meta,
      col = color_list,
      which = side,
      show_legend = side == "column",
      show_annotation_name = F
    )
  annot
}

compute_splits <- function(nbins, nstudies) {
  rep(1:nstudies, each = nbins)
}

build_heatmap <- function(mat,
                          figfile,
                          col_meta = NULL,
                          row_meta = NULL,
                          row_title = NULL,
                          column_title = NULL,
                          split = NULL) {
  message(figfile)
  
  if (!is.null(col_meta)) {
    col_meta <- build_annotation(col_meta)
  }
  if (!is.null(row_meta)) {
    row_meta <- build_annotation(row_meta, side = "row")
  }
  
  if (tools::file_ext(figfile) == 'png'){
    png(file = figfile,
      width = 700,
      height = 500)
  } else {
    pdf(file=figfile,
        width = 7,
        height = 5)
  }
  
  print({
    mat %>% ComplexHeatmap::Heatmap(
      .,
      cluster_rows = F,
      cluster_columns = F,
      show_column_names = F,
      show_row_names = F,
      row_title = row_title,
      column_title = column_title,
      row_order = 1:dim(.)[1],
      column_order = 1:dim(.)[1],
      top_annotation = col_meta,
      width = 1,
      height = 1,
      left_annotation = row_meta,
      col = circlize::colorRamp2(c(0, .5, 1), c("blue", "white", "red")),
      row_split = split,
      column_split = split
    )
  })
  dev.off()
}


df <- read.csv(fn)
results_name <- tools::file_path_sans_ext(fn)
outfile <-
  paste(outdir, basename(results_name), '.', extension, sep = "")
mat <- df %>%
  dplyr::select(!c(Pseudotime.Bin, Study, Lineage)) %>%
  as.matrix()
mat <- mat[,-1]

meta <- df %>% dplyr::select(c(Lineage, Study, Pseudotime.Bin))
n_lins <- meta$Lineage %>%
  unique() %>%
  length()

title <- basename(results_name)

n_bins <-
  max(meta[, "Pseudotime.Bin"]) - min(meta[, "Pseudotime.Bin"]) + 1
n_studies <- length(unique(meta$Study))


if (n_lins > 1) {
  build_heatmap(
    mat,
    outfile,
    col_meta = meta,
    row_meta = meta[, "Lineage", drop = F],
    column_title = title
  )
  if (split_lineages) {
    lins <- unique(meta$Lineage)
    for (lin in lins) {
      slice <- meta$Lineage == lin
      
      outf <- tools::file_path_sans_ext(outfile)
      outf <- paste0(outf, '_', lin, '.', extension)
      
      build_heatmap(
        mat[slice, slice],
        outf,
        col_meta = meta[slice, c('Study', 'Pseudotime.Bin')],
        split = compute_splits(n_bins, n_studies),
        column_title = paste(lin, title)
      )
    }
  }
  
  
} else {
  build_heatmap(
    mat,
    outfile,
    col_meta = meta[, c("Study", "Pseudotime.Bin")],
    column_title = title,
    split = compute_splits(n_bins, n_studies)
  )
}
