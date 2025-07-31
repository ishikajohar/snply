# plot_chromosomes.R
library(ggplot2)
library(dplyr)

#' Draw chromosome silhouettes (scaled, sorted 1–22, X, Y)
#'
#' @return ggplot
#' @export
draw_chromosome_silhouette <- function() {
  # 1) Setup
  chr_levels <- c(as.character(1:22), "X", "Y")
  lengths_bp <- c(
    249250621,243199373,198022430,191154276,180915260,
    171115067,159138663,146364022,141213431,135534747,
    135006516,133851895,115169878,107349540,102531392,
    90354753, 81195210, 78077248,  59128983,  63025520,
    48129895, 51304566,155270560,  59373566
  )
  cent_s <- c(
    121535434,  92326171,  90504854,  50000000,  44000000,
    58500000,  58000000,  43000000,  49000000,  40200000,
    34850000,  35800000,  17500000,  17800000,  19000000,
    36000000,  25700000,  23000000,  24600000,  26400000,
    13300000,  17600000,  12100000,  10200000
  )
  cent_e <- cent_s + 3e6  # 3 Mb centromere width

  # 2) Assemble master df
  df <- tibble(
    chr    = factor(chr_levels, levels=chr_levels),
    length = lengths_bp,
    cent_s = cent_s,
    cent_e = cent_e
  ) %>%
    mutate(row = as.integer(chr))     # 1..24

  max_len <- max(df$length)

  # 3) Duplicate for paternal/maternal, compute a length‐scaled vertical gap
  df2 <- bind_rows(
    df %>% mutate(copy="paternal", gap = (length/max_len)*0.4, y = row + gap),
    df %>% mutate(copy="maternal", gap = (length/max_len)*0.4, y = row - gap)
  )

  # 4) Plot
  ggplot(df2) +
    # chromosome arm
    geom_segment(aes(x=0, xend=length/1e6, y=y, yend=y),
                 linewidth=4, colour="#DDDDDD") +
    # centromere
    geom_rect(aes(xmin=cent_s/1e6, xmax=cent_e/1e6,
                  ymin=y - gap*0.2, ymax=y + gap*0.2),
              fill="black") +
    scale_y_continuous(breaks = seq_along(chr_levels),
                       labels = chr_levels,
                       expand = expansion(add = .5)) +
    labs(x="Position (Mbp)", y=NULL) +
    theme_minimal(base_size=12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.ticks.y       = element_blank()
    )
}


#' Paint introgressed SNPs on top of a chromosome silhouette
#'
#' @param base_plot Output of draw_chromosome_silhouette()
#' @param snps      Data.frame with columns chromosome, position (bp), archaic_count
#' @param tick_shape "segment" (default) or "rect"
#' @return ggplot
#' @export
paint_snps <- function(base_plot, snps, tick_shape = c("segment","rect")) {
  tick_shape <- match.arg(tick_shape)
  if (nrow(snps)==0) return(base_plot)

  library(dplyr)
  library(ggplot2)

  # 1) Clean & Mbp‐convert
  sn2 <- snps %>%
    mutate(
      chr      = as.character(chromosome),
      position = as.numeric(as.character(position))
    ) %>%
    filter(!is.na(position)) %>%
    mutate(pos_Mbp = position/1e6)

  # 2) Explicitly create paternal + maternal copies
  sn2b <- bind_rows(
    sn2 %>% mutate(copy = "paternal"),
    sn2 %>% mutate(copy = "maternal")
  )
  # 3) Grab silhouette coords (y & gap)
  sil <- base_plot$data %>% select(chr, copy, y, gap)

  # 4) Join so each SNP gets its y/gap
  ov <- sn2b %>%
    left_join(sil, by = c("chr","copy")) %>%
    mutate(
      ymin   = y - gap*0.2,
      ymax   = y + gap*0.2,
      colour = case_when(
        archaic_count==2 ~ "homo",
        archaic_count==1 ~ "hetero",
        TRUE             ~ NA_character_
      )
    )

  # 5) Build the overlay layer
  layer <- if (tick_shape=="segment") {
    geom_segment(aes(x=pos_Mbp, xend=pos_Mbp, y=ymin, yend=ymax, color=colour),
                 data=ov, size=0.5)
  } else {
    geom_rect(aes(xmin=pos_Mbp-0.05, xmax=pos_Mbp+0.05,
                  ymin=ymin, ymax=ymax, fill=colour),
              data=ov, color=NA)
  }

  # 6) Return combined plot
  base_plot +
    layer +
    scale_color_manual(values=c(homo="firebrick", hetero="seagreen"), name="Neandertal") +
    scale_fill_manual(values=c(homo="firebrick", hetero="seagreen"), guide="none") +
    guides(color = guide_legend(override.aes=list(size=3)))
}
