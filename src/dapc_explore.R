
library(adegenet)
library(data.table)
library(ggplot2)

############
# FUNCTION #
############

###########
# GLOBALS #
###########

adegenet_snps <- "output/11_stacks-populations/r0/genlight.Rds"

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

# impute?
imputed_snps <- tab(snps, NA.method = "mean")
imputed_genlight <- as.genlight(imputed_snps)
imputed_genlight$pop <- snps$pop


# run the PCA
snps_pca <- dudi.pca(imputed_genlight,
                     center = TRUE,
                     scale = FALSE,
                     scannf = FALSE,
                     nf = 9)
s.class(snps_pca$li, pop(snps), col=rainbow(nPop(snps)))
add.scatter.eig(snps_pca$eig[1:10], xax=1, yax=2)

# try a ggplot of the results
var_exp <- snps_pca$eig / sum(snps_pca$eig)
pv_labs <- paste0("PC", 1:length(var_exp), " (", round(var_exp*100, 1), "%)")
facet_labeller <- function(x) {
    list(pv_labs[x[, 1]])
}

pd <- data.table(snps_pca$li, keep.rownames = TRUE)
pd_long <- melt(pd, id.vars = "rn")
pd_long[, PC := as.numeric(gsub("[[:alpha:]]+", "", variable))]
pd_long[, population := snps$pop]


ggplot(pd_long, aes(x = population, y = value, colour = population)) +
    scale_colour_brewer(palette = "Paired") +
    facet_wrap(~ PC, ncol = 3, labeller = facet_labeller) +
    geom_point(position = position_jitter(width = 0.2),
               alpha = 0.6,
               shape = 16,
               size = 2)

# run dapc
dapc_results <- dapc(imputed_genlight,
                     pop(imputed_genlight),
                     n.pca = 100,
                     n.da = n_pop)

# what's driving LD3?
loadingplot(dapc_results$var.contr, axis=3, thres=.07, lab.jitter=1)

vc <- data.table(dapc_results$var.contr, id = snps$loc.names)
setorder(vc, -LD3)

x <- data.table(sapply(snps[, "3911_82_T"]$gen, as.integer),
           snps$pop)
table(x)

sl <- seploc(snps, block.size = 1)
snp1 <- data.table(truenames(sl[["414"]]), keep.rownames = TRUE)


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
    pct = round(100 * dapc_results$eig / sum(dapc_results$eig), 1)[1:n_pop - 1])
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
    xlab(pd_2d[, levels(pct_var)][1]) +
    ylab(pd_2d[, levels(pct_var)][2]) +
    theme_minimal(base_size = 12) +
    scale_colour_brewer(palette = "Paired",
                        guide = guide_legend(title = "Population")) +
    geom_point(alpha = 0.8,
               shape = 16,
               size = 2)

# components vs. population
ggplot(dapc_long, aes(x = population, y = value, colour = parasitised)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_blank(),
          strip.placement = "outside") +
    facet_wrap(~ pct_var, ncol = 3, strip.position = "left") +
    ggtitle("Discriminants vs. populations") +
    xlab(NULL) + ylab(NULL) +
    scale_colour_brewer(palette = "Set1",
                        guide = guide_legend(title = NULL)) +
    geom_point(position = position_jitter(width = 0.25),
               alpha = 0.6,
               shape = 16,
               size = 2)

