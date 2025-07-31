#' Monte-Carlo ancestry painting
#'
#' @param user_snps data.table with columns rsid, chromosome, position, genotype
#' @param ref_freq  data.table from load_aisnp_ref(), with columns SNP,CHR,A1,CEU,CHB,GIH,JPT,YRI
#' @param block_n   SNPs per block
#' @param iter      number of MC iterations
#' @param thresh    mono-ancestry vote threshold (0–1)
#' @return data.table: CHR, start, end, hap1, hap2
#' @export
paint_mc <- function(user_snps, ref_freq, block_n = 50, iter = 150, thresh = 0.8) {
  library(data.table)
  setDT(user_snps)
  # 1) Rename to match ref_freq
  setnames(user_snps,
           old = c("rsid","chromosome","position"),
           new = c("SNP","CHR","POS"))

  # 2) Prepare ref_freq: ensure pop columns numeric
  setDT(ref_freq)
  pop_cols <- c("CEU","CHB","GIH","JPT","YRI")
  ref_freq[, (pop_cols) := lapply(.SD, as.numeric), .SDcols = pop_cols]

  # 3) Merge on SNP & CHR
  # 3) Merge and validate
  dt <- merge(user_snps, ref_freq, by = c("SNP","CHR"))
  if (nrow(dt) < 100) stop("Too few overlapping AISNPs after merge")

  # 3b) Make sure our positions are numeric (so we can later divide by 1e6)
  dt[, POS := as.numeric(POS)]
  if (any(is.na(dt$POS))) stop("Some POS values could not be turned into numeric – check your input file")

  # 4) Split genotype and count reference allele (A1) copies
  dt[, `:=`(
    a1 = substr(genotype,1,1),
    a2 = substr(genotype,2,2),
    countA1 = (substr(genotype,1,1)==A1) + (substr(genotype,2,2)==A1)
  )]

  # 5) Compute per-population log-likelihoods under HWE
  for (p in pop_cols) {
    dt[, paste0("logL_",p) := {
      p_i <- pmax(get(p), 1e-6)
      log( ifelse(countA1==2,    p_i^2,
                  ifelse(countA1==1, 2*p_i*(1-p_i),
                         (1-p_i)^2)) )
    }]
  }

  # 6) Build numeric log-matrix
  logs_matrix <- as.matrix(dt[, paste0("logL_",pop_cols), with=FALSE])

  # 7) Define SNP blocks
  setorder(dt, CHR, POS)
  dt[, block := ceiling(.I / block_n), by=CHR]

  # 8) Initialize vote tallies
  votes <- dt[, .(
    hap1 = integer(length(pop_cols)),
    hap2 = integer(length(pop_cols))
  ), by=.(CHR, block)]

  # 9) Monte Carlo sampling
  set.seed(42)
  for (b in unique(dt$block)) {
    idx <- dt$block==b
    for (it in seq_len(iter)) {
      # Random phasing
      p1 <- dt$a1[idx]; p2 <- dt$a2[idx]
      het <- p1 != p2
      p1[het] <- ifelse(runif(sum(het))<.5, p1[het], p2[het])
      p2[het] <- ifelse(p1[het]==dt$a1[idx][het],
                        dt$a2[idx][het], dt$a1[idx][het])
      # Block log-sums
      subL <- logs_matrix[idx,,drop=FALSE]
      li1  <- colSums(subL * (p1==dt$A1[idx]) +
                        log(1 - exp(subL)) * (p1!=dt$A1[idx]))
      li2  <- colSums(subL * (p2==dt$A1[idx]) +
                        log(1 - exp(subL)) * (p2!=dt$A1[idx]))
      # Tally votes
      pop1 <- which.max(li1); pop2 <- which.max(li2)
      votes[CHR==dt$CHR[idx][1] & block==b, `:=`(
        hap1 = hap1 + (seq_along(pop_cols)==pop1),
        hap2 = hap2 + (seq_along(pop_cols)==pop2)
      )]
    }
  }

  # 10) Summarize into painted blocks
  res <- votes[, {
    tv  <- hap1 + hap2
    top <- which.max(tv)
    frac<- tv[top] / (2*iter)
    if (frac>=thresh) {
      h1<-h2<-pop_cols[top]
    } else {
      sec <- which.max(replace(tv, top, -Inf))
      h1<-pop_cols[top]; h2<-pop_cols[sec]
    }
    st <- dt[CHR==.BY$CHR & block==.BY$block, min(POS)]
    en <- dt[CHR==.BY$CHR & block==.BY$block, max(POS)]
    .(start=st, end=en, hap1=h1, hap2=h2)
  }, by=.(CHR, block)]

  # 11) (Optional) Smooth flickers
  res <- res[, .(
    start=min(start),
    end  =max(end),
    hap1=hap1[1],
    hap2=hap2[1]
  ), by=.(CHR, run=rleid(hap1, hap2))]
  res[, run:=NULL][]
}
