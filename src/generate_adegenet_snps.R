#!/usr/bin/env Rscript

library(adegenet)

###########
# GLOBALS #
###########

populations_genepop <- snakemake@input[["genepop"]]
adegenet_output <- snakemake@output[["adegenet_snps"]]
log_file <- snakemake@log[["log"]]

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read SNPs
snps <- read.genepop(populations_genepop)
pop(snps) <- gsub("^([[:alpha:]]+).*", "\\1", indNames(snps))

# write output
saveRDS(snps, adegenet_output)

# write session info
sessionInfo()
