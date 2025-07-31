#' Monte Carlo Ancestry Painting
#'
#' Given cleaned genotype data and an AISNP reference table, splits each
#' chromosome into blocks of adjacent SNPs and assigns ancestry to each block
#' for both haplotypes by Monte Carlo sampling.
#'
#' @param genotype_dt A \code{data.table} or \code{data.frame} with columns:
#'   \describe{
#'     \item{\code{rsid}}{SNP identifier (character).}
#'     \item{\code{chromosome}}{Chromosome (character or factor).}
#'     \item{\code{position}}{Genomic position (numeric or integer).}
#'     \item{\code{genotype}}{Two‐character string (e.g. "AG"), uppercase.}
#'   }
#' @param freq_dt A \code{data.table} from \code{load_aisnp_ref()} containing
#'   columns \code{SNP}, \code{CHR}, \code{A1}, \code{CEU}, \code{CHB},
#'   \code{GIH}, \code{JPT}, \code{YRI}, and optionally \code{FST}.  Frequencies
#'   must already be clipped to (0,1).
#' @param block_size Integer ≥1.  Number of SNPs per block (default 50).
#' @param n_sim Integer ≥1. Number of Monte Carlo samples per block (default 1000).
#'
#' @return A \code{data.table} with one row per chromosome‐block and columns:
#'   \describe{
#'     \item{\code{CHR}}{Chromosome.}
#'     \item{\code{block}}{Block index (within that chromosome).}
#'     \item{\code{start_pos}}{Minimum SNP position in block.}
#'     \item{\code{end_pos}}{Maximum SNP position in block.}
#'     \item{\code{n_snps}}{Number of SNPs in block.}
#'     \item{\code{CEU_votes}, …, \code{YRI_votes}}{Vote counts (out of \code{n_sim}) for each population.}
#'     \item{\code{top_ancestry}}{Most‐voted population for the block.}
#'     \item{\code{top_frac}}{Fraction of votes for \code{top_ancestry}.}
#'   }
#'
#' @examples
#' \dontrun{
#' gens  <- user_df()            # cleaned genotype data.table
#' freqs <- load_aisnp_ref()     # CSV‐loaded AISNP reference
#' blocks <- paint_mc(gens, freqs, block_size = 50, n_sim = 1000)
#' }
#'
#' @import data.table
#' @export
paint_mc <- function(genotype_dt,
                     freq_dt,
                     block_size = 50,
                     n_sim      = 1000) {

  # — defensive checks —
  library(data.table)
  genotype_dt <- as.data.table(genotype_dt)
  freq_dt     <- as.data.table(freq_dt)
  stopifnot(
    all(c("rsid","chromosome","position","genotype") %in% names(genotype_dt)),
    all(c("SNP","CHR","A1","CEU","CHB","GIH","JPT","YRI") %in% names(freq_dt)),
    is.numeric(block_size), block_size >= 1,
    is.numeric(n_sim),      n_sim >= 1
  )
  if (nrow(genotype_dt)==0) stop("genotype_dt is empty")
  if (nrow(freq_dt)==0)     stop("freq_dt is empty")

  # Population labels
  pops <- c("CEU","CHB","GIH","JPT","YRI")

  # Convert frequency cols to numeric (in‐place)
  freq_dt[, (pops) := lapply(.SD, as.numeric), .SDcols = pops]

  # Only work with SNPs present in both
  common <- intersect(genotype_dt$rsid, freq_dt$SNP)
  if (length(common)==0) stop("No overlapping SNPs between genotype_dt and freq_dt")

  # Subset & merge genotype with frequencies
  gt <- genotype_dt[rsid %in% common]
  setnames(gt, "chromosome", "CHR")        # align to freq_dt naming
  setnames(gt, "position",   "POS")
  merged <- merge(gt, freq_dt, by.x = "rsid", by.y = "SNP", all = FALSE)

  # Sort by chromosome & position
  setorder(merged, CHR, POS)

  # Prepare container for results
  results_list <- vector("list", length = 0L)

  # Process each chromosome separately
  for (ch in unique(merged$CHR)) {
    data_ch <- merged[CHR == ch]
    n_snps  <- nrow(data_ch)
    n_blocks <- ceiling(n_snps / block_size)

    # Preallocate vote matrix (blocks × populations)
    vote_mat <- matrix(0L,
                       nrow = n_blocks,
                       ncol = length(pops),
                       dimnames = list(NULL, pops))

    # Monte Carlo sampling
    for (sim in seq_len(n_sim)) {
      # For each block
      for (b in seq_len(n_blocks)) {
        i_start <- (b - 1L) * block_size + 1L
        i_end   <- min(b * block_size, n_snps)
        block_snps <- data_ch[i_start:i_end]

        # Initialize log‐likelihood for this simulated haplotype
        logL <- numeric(length(pops))

        # Sample one allele per SNP
        for (j in seq_len(nrow(block_snps))) {
          snp   <- block_snps[j]
          a1    <- substr(snp$genotype, 1, 1)
          a2    <- substr(snp$genotype, 2, 2)
          pick  <- if (a1 == a2) a1 else if (runif(1) < 0.5) a1 else a2
          # Reference allele:
          refA  <- snp$A1

          # Compute log‐likelihood increment
          p_vec <- as.numeric(snp[, ..pops])
          if (pick == refA) {
            logL <- logL + log(p_vec)
          } else {
            logL <- logL + log(1 - p_vec)
          }
        }

        # Convert to probabilities
        maxL <- max(logL)
        w    <- exp(logL - maxL)    # avoid underflow
        probs <- w / sum(w)

        # Draw one vote
        vote_mat[b, which.max(probs)] <- vote_mat[b, which.max(probs)] + 1L
      }
    }

    # Build per‐block summary for this chromosome
    dt_chr <- data.table(
      CHR   = ch,
      block = seq_len(n_blocks),
      start_pos = vapply(seq_len(n_blocks), function(b) {
        (merged$POS[ (b-1)*block_size + 1 ])
      }, integer(1)),
      end_pos   = vapply(seq_len(n_blocks), function(b) {
        merged$POS[ min(b*block_size, n_snps) ]
      }, integer(1)),
      n_snps = vapply(seq_len(n_blocks), function(b) {
        min(block_size, n_snps - (b-1)*block_size)
      }, integer(1))
    )

    # Attach vote counts
    for (pop in pops) {
      dt_chr[[ paste0(pop, "_votes") ]] <- vote_mat[, pop]
    }

    # Determine top ancestry and fraction
    dt_chr[, top_ancestry := pops[ which.max(.SD) ],
           .SDcols = paste0(pops, "_votes"), by = block]
    dt_chr[, top_frac := do.call(pmax, .SD) / n_sim,
           .SDcols = paste0(pops, "_votes")]

    results_list[[as.character(ch)]] <- dt_chr
  }

  # Combine all chromosomes
  result_dt <- rbindlist(results_list, use.names = TRUE)
  return(result_dt)
}
