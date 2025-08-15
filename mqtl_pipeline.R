#!/usr/bin/env Rscript

# mqtl_pipeline.R
# Perform QTL mapping using R/qtl for each phenotype in a cross file.
# The script scans the genome and fits multiple QTL models.
# Usage: Rscript mqtl_pipeline.R -i cross.csv -o results --permutations 1000
#
# Requirements: R packages qtl and optparse.

suppressPackageStartupMessages({
  library(qtl)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Path to csv cross file (required)"),
  make_option(c("-o", "--outdir"), type="character", default="results",
              help="Directory to save output files [default %default]"),
  make_option(c("-p", "--permutations"), type="integer", default=1000,
              help="Number of permutations for LOD threshold [default %default]"),
  make_option(c("--bc-gen"), type="integer", default=1,
              help="Number of backcross generations [default %default]"),
  make_option(c("--f-gen"), type="integer", default=6,
              help="Number of inbreeding generations [default %default]"),
  make_option(c("--step"), type="double", default=0,
              help="Step size for genotype probability calculation [default %default]"),
  make_option(c("--na"), type="character", default="NA",
              help="Missing value string in input file [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)) {
  stop("Input cross file is required. Use -i or --input to specify the file.")
}

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive=TRUE)

#' Run QTL mapping for a single phenotype.
#' @param geno cross object from R/qtl.
#' @param pheno.col name of the phenotype column.
#' @param perms number of permutations for threshold estimation.
#' @param outdir directory to write results.
run_mapper <- function(geno, pheno.col, perms, outdir) {
  out <- scanone(geno, pheno.col=pheno.col, method="hk")
  out.perms <- scanone(geno, pheno.col=pheno.col, method="hk",
                       n.perm=perms, verbose=TRUE)
  threshold <- summary(out.perms, alpha=.05)[1,1]
  write.table(summary(out.perms, alpha=.05),
              file=file.path(outdir, paste0("LOD_threshold_for_", pheno.col, ".txt")),
              col.names=NA, sep="\t")

  single.summary <- summary(out, threshold=threshold, lodcolumn=1)
  write.csv(out,
            file=file.path(outdir, paste0(pheno.col, "_scanone.csv")),
            quote=FALSE, row.names=TRUE)
  write.csv(single.summary,
            file=file.path(outdir, paste0(pheno.col, "_scanone_summary.csv")),
            quote=FALSE, row.names=TRUE)

  if (nrow(single.summary) == 0) return()

  chr <- single.summary$chr
  pos <- single.summary$pos
  qtl <- makeqtl(geno, chr=chr, pos=pos, what="prob")

  looper <- TRUE
  while (looper) {
    term <- rownames(summary(qtl))
    formula1 <- as.formula(paste("y ~", paste(term, collapse=" + ")))

    qtl <- refineqtl(geno, pheno.col=pheno.col, qtl=qtl,
                     formula=formula1, method="hk")
    out.aq <- addqtl(geno, pheno.col=pheno.col, qtl=qtl,
                     formula=formula1, method="hk")
    additional <- summary(out.aq, threshold=threshold, lodcolumn=1)

    if (nrow(additional) > 0) {
      qtl <- addtoqtl(geno, qtl, additional$chr, additional$pos)
    } else {
      looper <- FALSE
    }
  }

  pdf(file.path(outdir, paste0(pheno.col, "_mqm_qtl.pdf")), width=6, height=4)
  plotLodProfile(qtl, main=pheno.col, ylab="LOD", cex=0.4, lwd=1)
  add.threshold(out, perms=out.perms, alpha=0.05, gap=0)
  dev.off()

  pdf(file.path(outdir, paste0(pheno.col, "_mqm_all.pdf")), width=6, height=4)
  plot(out, lwd=1.5, gap=0, bandcol="gray70", incl.markers=TRUE,
       main=pheno.col,
       xlab=paste("Threshold for alpha=.05 using", perms, "permutations"))
  add.threshold(out, perms=out.perms, alpha=0.05, gap=0)
  dev.off()

  out.mqm <- fitqtl(geno, pheno.col=pheno.col, qtl=qtl,
                    formula=formula1, method="hk", dropone=TRUE,
                    get.ests=TRUE)
  write.table(summary(out.mqm)[[1]],
              file=file.path(outdir, paste0(pheno.col, "_mqm_full.txt")),
              quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
  write.table(summary(out.mqm)[[2]],
              file=file.path(outdir, paste0(pheno.col, "_mqm_dropone.txt")),
              quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
  write.table(summary(out.mqm)[[3]],
              file=file.path(outdir, paste0(pheno.col, "_mqm_effect.txt")),
              quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

  lodprof <- attr(qtl, "lodprofile")
  cat(capture.output(lodprof),
      file=file.path(outdir, paste0(pheno.col, "_mqm_lodprof.txt")),
      sep="\n", append=FALSE)

  cat("", file=file.path(outdir, paste0(pheno.col, "_mqm_qtlint.txt")),
      sep="\n")
  for (i in seq_along(lodprof)) {
    qtlint <- lodint(qtl, qtl.index=i)
    cat(capture.output(qtlint),
        file=file.path(outdir, paste0(pheno.col, "_mqm_qtlint.txt")),
        sep="\n", append=TRUE)
  }
}

geno1 <- read.cross("csv", file=opt$input, na.strings=opt$na,
                    genotypes=c("AA", "BB"),
                    alleles=c("A", "B"), estimate.map=FALSE)
attr(geno1, "class")[1] <- "f2"
geno <- convert2bcsft(geno1, BC.gen=opt$bc_gen,
                      F.gen=opt$f_gen, estimate.map=FALSE)
geno <- jittermap(geno)
geno <- calc.genoprob(geno, step=opt$step, off.end=0,
                      error.prob=0.0001, map.function="kosambi",
                      stepwidth="fixed")

for (pheno in names(geno$pheno)) {
  run_mapper(geno, pheno, opt$permutations, opt$outdir)
}
