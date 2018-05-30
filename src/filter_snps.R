#!/usr/bin/env Rscript

library(SNPRelate)

###########
# GLOBALS #
###########

gds_file <- snakemake@input[["gds"]]
pedfn <- snakemake@params[["pedfn"]]
maf <- snakemake@params[["maf"]]
missing_rate <- snakemake@params[["missing_rate"]]
sample_missing_quantile <- snakemake@params[["sample_missing_quantile"]]
log_file <- snakemake@log[["log"]]

# dev 
# gds_file <- "output/11_stacks-populations/r0/populations.snps.gds"
# maf <- 0.05
# missing_rate <- 0.2
# sample_missing_quantile <- 0.8

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read the GDS
gds <- snpgdsOpen(gds_file)

# subset the SNPs
snpset <- snpgdsSelectSNP(
    gds,
    maf = maf,
    autosome.only = FALSE,
    verbose = TRUE,
    missing.rate = missing_rate)
snpset_ids <- unlist(snpset)

# get a vector of individuals missing < sample_missing_quantile missing_rate
sample_missing_rates <- snpgdsSampMissRate(gdsobj = gds,
                   snp.id = snpset_ids,
                   with.id = TRUE)
kept_individuals <- sample_missing_rates[
    sample_missing_rates < quantile(sample_missing_rates, 
                                    sample_missing_quantile)]

# write the ped
snpgdsGDS2PED(gds,
              ped.fn = pedfn,
              snp.id = snpset_ids,
              sample.id = names(kept_individuals),
              verbose = TRUE)

# log session info
sessionInfo()
