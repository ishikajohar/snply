#' Plot Chromosome Painting for One Chromosome
#'
#' Given a table of painted ancestry blocks, draw a Plotly figure showing two
#' haplotypes (rows) along a single chromosome. Each block is a colored
#' segment whose color corresponds to the inferred ancestry for haplotype 1
#' or 2.
#'
#' @param res A \code{data.frame} or \code{data.table} with columns
#'   \describe{
#'     \item{\code{CHR}}{Chromosome identifier (character or integer).}
#'     \item{\code{block}}{Block index (integer).}
#'     \item{\code{start}}{Start position of the block (base‐pairs).}
#'     \item{\code{end}}{End position of the block (base‐pairs).}
#'     \item{\code{hap1}}{Ancestry label for haplotype 1 (must match names of \code{pop_cols}).}
#'     \item{\code{hap2}}{Ancestry label for haplotype 2 (must match names of \code{pop_cols}).}
#'   }
#' @param chrom Chromosome to plot. Must exactly match one of the values in \code{res\$CHR}.
#' @param pop_cols Named character vector of hex‐color codes, where names correspond
#'   to ancestry labels (e.g. \code{c(CEU = "#4363d8", YRI = "#e41a1c", …)}).
#'
#' @return A \code{plotly} object with two horizontal rows of colored segments.
#'
#' @examples
#' \dontrun{
#'   blocks <- paint_mc(user_df(), load_aisnp_ref())
#'   plot_paint_blocks(blocks, chrom = "1", pop_cols)
#' }
#'
#' @import plotly
#' @export
plot_paint_blocks <- function(res, chrom, pop_cols) {
  # sub‐set to the chosen chromosome
  dat <- res[ res$CHR == chrom, ]
  if (nrow(dat) == 0) {
    warning("No blocks found for chromosome ", chrom)
    return(NULL)
  }

  # Initialize empty plotly object
  p <- plot_ly()

  # Add one segment per block per haplotype
  for (i in seq_len(nrow(dat))) {
    seg <- dat[i, ]

    # Hap1 (plotted at y = 1)
    p <- add_segments(
      p,
      x     = seg$start / 1e6,
      xend  = seg$end   / 1e6,
      y     = 1,
      yend  = 1,
      line  = list(color = pop_cols[[ seg$hap1 ]], width = 7),
      hoverinfo = "text",
      text  = paste0(
        "Block ", seg$block,
        "<br>", seg$hap1, " (hap1)",
        "<br>", format(seg$start, big.mark = ","), "–",
        format(seg$end,   big.mark = ",")
      )
    )

    # Hap2 (plotted at y = 0)
    p <- add_segments(
      p,
      x     = seg$start / 1e6,
      xend  = seg$end   / 1e6,
      y     = 0,
      yend  = 0,
      line  = list(color = pop_cols[[ seg$hap2 ]], width = 7),
      hoverinfo = "text",
      text  = paste0(
        "Block ", seg$block,
        "<br>", seg$hap2, " (hap2)",
        "<br>", format(seg$start, big.mark = ","), "–",
        format(seg$end,   big.mark = ",")
      )
    )
  }

  # Final layout adjustments
  layout(
    p,
    yaxis = list(
      title      = "Haplotype",
      tickvals   = c(0, 1),
      ticktext   = c("2",  "1"),
      zeroline   = FALSE
    ),
    xaxis = list(
      title    = "Position (Mbp)",
      zeroline = FALSE
    ),
    showlegend = FALSE
  )
}
