#' Plot genome-wide ancestry painting (all chromosomes)
#'
#' @param all_blocks data.table with columns CHR, block, start, end, hap1, hap2
#' @param pop_cols   named vector of colors for each ancestry
#' @return ggplot
#' @export
plot_all_chromosomes <- function(all_blocks, pop_cols) {
  library(data.table)
  library(ggplot2)

  # 1) Melt hap1/hap2 into long format
  blocks_long <- melt(
    setDT(all_blocks),
    id.vars       = c("CHR","start","end"),
    measure       = list(c("hap1","hap2")),
    variable.name = "haplotype",
    value.name    = "ancestry"
  )
  # map haplotype factor to y-coordinate
  blocks_long[, y := fifelse(haplotype == "hap1", 1, 0)]

  # 2) Draw segments, faceted by chromosome
  ggplot(blocks_long,
         aes(x = start/1e6, xend = end/1e6, y = y, yend = y, color = ancestry)) +
    geom_segment(size = 3) +
    facet_wrap(~CHR, scales = "free_x", ncol = 4) +
    scale_color_manual(values = pop_cols) +
    scale_y_continuous(
      breaks = c(0,1),
      labels = c("Hap2","Hap1"),
      expand = expansion(add = 0.1)
    ) +
    labs(x = "Position (Mbp)", y = NULL, color = "Ancestry") +
    theme_minimal(base_size = 12) +
    theme(
      strip.text       = element_text(face = "bold"),
      axis.ticks.y     = element_blank(),
      panel.grid.major.y = element_blank()
    )
}
