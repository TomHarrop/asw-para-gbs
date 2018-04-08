
library(adegenet)
library(data.table)
library(ggplot2)

############
# FUNCTION #
############

###########
# GLOBALS #
###########

key_file <- 'data/SQ0003.txt'
adegenet_snps <- "output/11_stacks-populations/r0.8/snps.Rds"

########
# MAIN #
########

snps <- readRDS(adegenet_snps)

# set up individual populations
ind_pop <- data.table(individual = indNames(snps))
ind_pop[, population := gsub("^([[:alpha:]]+).*", "\\1", individual)]
ind_pop[endsWith(individual, "p"), parasitised := "parasitised"]
ind_pop[!endsWith(individual, "p"), parasitised := "unparasitised"]
ind_pop[, pop_para := paste(population, parasitised, sep = "_")]
pop(snps) <- ind_pop[, pop_para]
n_pop <- ind_pop[, length(unique(pop_para))]

# calculate minor allele frequency
maf <- minorAllele(snps)
keep_snps <- unique(names(maf[maf > 0.05]))
filtered_snps <- snps[loc = keep_snps, drop = TRUE]

# run dapc
dapc_results <- dapc(filtered_snps,
                     pop(filtered_snps),
                     n.pca = 100,
                     n.da = n_pop)

# what's driving LD2?
loadingplot(dapc_results$var.contr, axis=2, thres=.07, lab.jitter=1)
vc <- data.table(dapc_results$var.contr, keep.rownames = TRUE)
setnames(vc, "rn" , "id")
setorder(vc, -LD2)


x <- vc[which.max(LD1), id]


sl <- seploc(filtered_snps)
snp1 <- data.table(truenames(sl[["13343_89"]]), keep.rownames = TRUE)


# extract components
dapc_dt <- data.table(dapc_results$ind.coord, keep.rownames = TRUE)
setnames(dapc_dt, "rn", "individual")
dapc_dt[, population := gsub("^([[:alpha:]]+).*", "\\1", individual)]
dapc_long <- melt(dapc_dt,
                  id.vars = c("individual", "population"))
dapc_long[endsWith(individual, "p"), parasitised := "parasitised"]
dapc_long[!endsWith(individual, "p"), parasitised := "unparasitised"]
dapc_long[, pop_para := paste(population, parasitised, sep = "_")]

# add contributions
pct <- data.table(
    variable = colnames(dapc_results$ind.coord),
    pct = round(100 * dapc_results$eig / sum(dapc_results$eig), 1))
pct_var <- pct[, structure(paste0(variable, " (", pct, "%)"),
                           .Names = variable)]
dapc_long[, pct_var := plyr::revalue(variable, pct_var)]

# 2d plots
pd_2d <- dcast(dapc_long,
               individual + population + parasitised + pop_para ~ pct_var,
               value.var = "value")
ggplot(pd_2d, aes(x = get(dapc_long[, levels(pct_var)][1]),
                  y = get(dapc_long[, levels(pct_var)][2]),
                  colour = pop_para)) +
    xlab(pd[, levels(pct_var)][1]) + ylab(pd[, levels(pct_var)][2]) +
    theme_minimal(base_size = 12) +
    scale_colour_brewer(palette = "Paired",
                        guide = guide_legend(title = "Population")) +
    geom_point(alpha = 0.8,
               shape = 16,
               size = 2)

# components vs. population
ggplot(dapc_long, aes(x = population, y = value, colour = pop_para)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_blank(),
          strip.placement = "outside") +
    facet_wrap(~ pct_var, ncol = 3, strip.position = "left") +
    ggtitle("Discriminants vs. populations") +
    xlab(NULL) + ylab(NULL) +
    scale_colour_brewer(palette = "Paired",
                        guide = guide_legend(title = NULL)) +
    geom_point(position = position_jitter(width = 0.25),
               alpha = 0.6,
               shape = 16,
               size = 2)

