#' Count copies of an alternate allele in an *un-phased* genotype
#'
#' @param geno   Character scalar such as `"AA"`, `"AG"`, `"CT"`.
#'               May contain one or two letters (or `NA`).
#' @param alt_alle Single-letter character vector – the alternate (archaic)
#'                 allele whose dosage you want to count.
#'
#' @return An integer: `0`, `1`, or `2` (or `NA` if the genotype was `NA`
#'         or empty).
#' @examples
#' alt_dosage("AG", "G")   # 1
#' alt_dosage("GG", "G")   # 2
#' alt_dosage("AA", "G")   # 0
#' @export
alt_dosage <- function(geno, alt_alle) {
  if (is.na(geno) || nchar(geno) == 0) return(NA_integer_)
  sum(strsplit(geno, "")[[1]] == alt_alle)
}

#' Log-likelihood of one un-phased genotype given an alternate-allele frequency
#'
#' @param dosage Integer 0, 1, or 2 – number of alternate-allele copies
#'               (output of [alt_dosage()]).
#' @param p      Numeric scalar 0 ≤ *p* ≤ 1 – frequency of the alternate
#'               allele in the reference population.
#'
#' @return Natural-log likelihood `log P(genotype | p)`.
#' @keywords internal
genotype_loglike <- function(dosage, p) {
  switch(as.character(dosage),
         "0" = log((1 - p)^2  + 1e-16),
         "1" = log(2 * p * (1 - p) + 1e-16),
         "2" = log(p^2 + 1e-16),
         log(1e-16))   # fallback for unexpected dosage
}

#' Assign ancestry to every fixed-size SNP window
#'
#' The function scores each 200-SNP window (or any size used to build
#' the `freq_list`) for all reference populations and picks the most likely
#' one.
#'
#' @param user_df   `data.frame` after merging the user SNPs with the
#'                  reference SNP map **and** adding `alt_count`
#'                  (number of alt-allele copies). Must contain the
#'                  columns `window_id`, `rsid`, `alt_count`,
#'                  `chromosome`, `position`.
#' @param freq_list Named list produced in *STEP 0*; each element is a
#'                  matrix of allele frequencies (rows = rsids,
#'                  columns = populations) for one window.
#'
#' @return A data-frame with one row per window and the columns
#'   \describe{
#'     \item{window_id}{integer window index}
#'     \item{chrom}{chromosome}
#'     \item{start_bp,end_bp}{first / last SNP position in the window}
#'     \item{best_pop}{population with the highest log-likelihood}
#'   }
#' @export
assign_windows <- function(user_df, freq_list) {

  by_win <- split(user_df, user_df$window_id)

  window_call <- lapply(by_win, function(df) {

    win  <- unique(df$window_id)
    freq <- freq_list[[as.character(win)]]

    # align rows of the frequency matrix to user SNP order
    freq <- freq[df$rsid, , drop = FALSE]

    # vectorised: one log-likelihood per population
    ll <- colSums(mapply(genotype_loglike,
                         dosage = df$alt_count,
                         p      = as.data.frame(freq),
                         SIMPLIFY = FALSE))

    data.frame(window_id = win,
               chrom      = unique(df$chromosome),
               start_bp   = min(df$position),
               end_bp     = max(df$position),
               best_pop   = names(which.max(ll)),
               stringsAsFactors = FALSE)
  })

  data.table::rbindlist(window_call)
}
