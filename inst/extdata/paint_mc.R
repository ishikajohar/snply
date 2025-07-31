#' Monte-Carlo ancestry painting
#' @param user_snps data.table (rsid, chromosome, position, genotype)
#' @param ref_freq  data.table from load_aisnp_ref()
#' @param block_n   SNPs per block
#' @param iter      number of MC iterations
#' @param thresh    mono-ancestry vote threshold (0–1)
#' @return data.table: chr, start, end, hap1, hap2
#' @export
paint_mc <- function(user_snps, ref_freq, block_n = 50,
                     iter = 150, thresh = 0.8) {

  library(data.table)
  setDT(user_snps)
  keycols <- c("rsid","chromosome","position")
  setnames(user_snps, keycols, c("SNP","CHR","POS"))   # align names
  setkey(user_snps, SNP)

  ## merge AISNP panel with user genotypes
  dt <- ref_freq[user_snps, nomatch = 0]      # inner join
  if (nrow(dt) < 100) stop("Too few AISNPs found in user data")

  ## split genotype into alleles
  dt[, `:=`(a1 = substr(genotype,1,1),
            a2 = substr(genotype,2,2))]

  ## block index
  setorder(dt, CHR, POS)
  dt[, block := ceiling(seq_len(.N) / block_n), by = CHR]

  pops  <- grep("^[A-Z]{2,}$", names(ref_freq), value = TRUE)
  logs  <- lapply(pops, function(p) log(dt[[p]]))          # log(freq)
  comp  <- lapply(pops, function(p) log(1 - dt[[p]]))      # log(1-freq)

  votes <- dt[, .(cnt  = vector("list", 2)), by = .(CHR, block)][
    , `:=`(hap1 = rep(0,length(pops)), hap2 = rep(0,length(pops)))]

  ## Monte-Carlo loop
  set.seed(42)
  for (b in unique(dt$block)) {
    snp_idx <- dt$block == b
    L1 <- L2 <- rep(0, length(pops))           # log-lik per pop

    for (it in seq_len(iter)) {
      pick1 <- dt$a1[snp_idx]                  # copy early
      pick2 <- dt$a2[snp_idx]
      het   <- pick1 != pick2
      pick1[het] <- ifelse(runif(sum(het)) < .5,
                           dt$a1[snp_idx][het],
                           dt$a2[snp_idx][het])
      pick2[het] <- ifelse(pick1[het] == dt$a1[snp_idx][het],
                           dt$a2[snp_idx][het],
                           dt$a1[snp_idx][het])

      ## log-likelihood vectorised
      li1 <- li2 <- rep(0, length(pops))
      for (i in seq_along(pops)) {
        li1[i] <- sum(ifelse(pick1 == dt$A1[snp_idx], logs[[i]][snp_idx],
                             comp[[i]][snp_idx]))
        li2[i] <- sum(ifelse(pick2 == dt$A1[snp_idx], logs[[i]][snp_idx],
                             comp[[i]][snp_idx]))
      }
      votes[CHR == dt$CHR[snp_idx][1] & block == b,
            hap1 := hap1 + (seq_along(pops) == which.max(li1))]
      votes[CHR == dt$CHR[snp_idx][1] & block == b,
            hap2 := hap2 + (seq_along(pops) == which.max(li2))]
    }
  }

  ## translate votes to ancestry labels
  res <- votes[, {
    tot1 <- hap1; tot2 <- hap2
    top1 <- which.max(tot1); top2 <- which.max(tot2)
    prop1 <- max(tot1)/iter;  prop2 <- max(tot2)/iter
    call1 <- pops[top1]; call2 <- pops[top2]
    if (prop1 >= thresh & prop2 >= thresh & call1 == call2) {
      call2 <- call1
    }
    st <- dt[CHR == .BY[[1]] & block == .BY[[2]],
             min(POS)]; en <- dt[CHR == .BY[[1]] & block == .BY[[2]],
                                 max(POS)]
    .(start = st, end = en, hap1 = call1, hap2 = call2)
  }, by = .(CHR = CHR, block)]
  res[]
}
