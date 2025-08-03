#' Perform Chromosome Painting Analysis
#'
#' Compares user genotypes against pre-calculated reference population allele
#' frequencies to infer local ancestry segments across the genome.
#'
#' @param user_data A data.frame with columns:
#'   - **rsid** (optional): SNP ID
#'   - **chromosome** (character or numeric)
#'   - **position** (numeric)
#'   - **genotype** (two-letter string, e.g. "AG")
#' @param ref_freq_data A data.frame or data.table with columns:
#'   - **chr**, **start_pos**, **end_pos**, **freq_data**
#'     where **freq_data** is a list-column; each element is
#'     a named list of frequency vectors per population
#' @param snp_map_ref A data.frame with columns:
#'   - **snp.name**: "chr:pos:ref:alt"
#'   - **allele.1**, **allele.2**: reference/alternate alleles
#' @param target_pops Character vector of population codes.
#'   Default: c("GBR","FIN","PUR","PJL","GIH","YRI","IBS")
#' @param min_snps_per_window Minimum overlapping SNPs per window. Default: 10.
#' @return A data.frame with columns: chr, start, end, ancestry.
#' @import data.table
#' @importFrom stats setNames
#' @export
paint_chromosomes <- function(user_data,
                              ref_freq_data,
                              snp_map_ref,
                              target_pops = c("GBR","FIN","PUR","PJL","GIH","YRI","IBS"),
                              min_snps_per_window = 10) {
  # DEBUG
  message("▶ Entering paint_chromosomes()")
  message("   user_data rows:       ", nrow(user_data))
  message("   ref_freq_data windows:", nrow(ref_freq_data))
  message("   snp_map_ref rows:     ", nrow(snp_map_ref))

  # 1) build user snp.name = chr:pos:ref:alt
  u <- user_data %>%
    mutate(
      chromosome = as.character(chromosome),
      position   = as.integer(position),
      genotype   = toupper(genotype),
      allele1    = substr(genotype,1,1),
      allele2    = substr(genotype,2,2),
      snp.name   = paste(chromosome, position, allele1, allele2, sep=":")
    ) %>%
    filter(nchar(genotype)==2)
  message("   after building snp.name, u rows: ", nrow(u))

  # 2) join to map
  m <- dplyr::inner_join(u, snp_map_ref, by="snp.name")
  message("   after inner_join, matched rows: ", nrow(m))
  if (nrow(m)==0) {
    warning("No user SNPs matched the reference map.")
    return(data.frame(
      chr      = character(),
      start    = numeric(),
      end      = numeric(),
      ancestry = character(),
      stringsAsFactors = FALSE
    ))
  }

  # 3) convert to numeric dosage of allele2
  m <- m %>% mutate(
    dosage = as.integer(substr(genotype,1,1)==allele2) +
      as.integer(substr(genotype,2,2)==allele2)
  )

  # 4) window-by-window
  segs <- list()
  for (i in seq_len(nrow(ref_freq_data))) {
    win    <- ref_freq_data[i, ]
    freq_l <- win$freq_data[[1]]
    if (length(freq_l)==0) next

    # gather all SNP IDs in this window
    all_ids <- unique(unlist(lapply(freq_l, names)))
    in_win  <- intersect(all_ids, m$snp.name)
    if (length(in_win) < min_snps_per_window) next

    user_vals <- setNames(m$dosage[match(in_win, m$snp.name)], in_win)
    pop_scores <- setNames(rep(NA_real_, length(target_pops)), target_pops)

    for (pop in target_pops) {
      pf <- freq_l[[pop]]
      if (!is.null(pf)) {
        common <- intersect(in_win, names(pf))
        if (length(common) >= min_snps_per_window) {
          pop_scores[pop] <- mean(abs(user_vals[common] - pf[common]), na.rm=TRUE)
        }
      }
    }
    vs <- !is.na(pop_scores)
    if (any(vs)) {
      best <- names(which.min(pop_scores[vs]))
      segs[[length(segs)+1]] <- data.frame(
        chr      = win$chr,
        start    = win$start_pos,
        end      = win$end_pos,
        ancestry = best,
        stringsAsFactors = FALSE
      )
    }
  }
  message("▶ Window processing complete")

  if (length(segs)>0) {
    result <- do.call(rbind, segs)
    message("▶ Painting complete: assigned ", nrow(result), " segments.")
    return(result)
  } else {
    message("▶ Painting complete: no segments assigned.")
    return(data.frame(
      chr      = character(),
      start    = numeric(),
      end      = numeric(),
      ancestry = character(),
      stringsAsFactors = FALSE
    ))
  }
}
