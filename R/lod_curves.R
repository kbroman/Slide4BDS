# part of Fig 1 from Broman et al. (2019) R/qtl2 paper
# (reproduction of Fig 5 from Gatti et al (2014))
#
# "We regressed log neutrophil counts on founder allele dosages at
# each marker using a kinship correction with sex and log white blood cell counts as covariates"

library(qtl2)
set.seed(33003221)

# pseudomarker maps
file <- "cache/maps_n_phe.RData"
if(file.exists(file)) {
    load(file)
} else {
    gmap <- insert_pseudomarkers(do$gmap, stepwidth="max", step=0.2)
    pmap <- interp_map(gmap, do$gmap, do$pmap)

    # phenotypes and covariates
    phe <- log10(do$pheno[,2])
    covar <- cbind(sex=(do$is_female*1),
                   wbc=log10(do$pheno[,1]))

    save(gmap, pmap, phe, covar, file=file)
}

# genome scan with additive model
file <- "cache/out_add.rds"
if(file.exists(file)) {
    out_add <- readRDS(file)
} else {
    out_add <- scan1(apr, phe, k, addcovar=covar, cores=0)
    saveRDS(out_add, file)
}

# functions to query snps and genes
qv <- create_variant_query_func("~/Data/CCdb/cc_variants.sqlite")
qg <- create_gene_query_func("~/Data/CCdb/mouse_genes_mgi.sqlite")

# GWAS at SNPs


file <- "cache/blup_c1.rds"
if(file.exists(file)) {
    blup <- readRDS(file)
} else {
    blup <- scan1blup(apr[,1], phe, k[1], addcovar=covar, cores=0)
    saveRDS(blup, file)
}


# line colors
altcolor <- "green4"
linecolor <- "violetred"



# load permutation results
operm_add <- readRDS("cache/operm_add.rds")

# calculate thresholds
thr_add <- summary(operm_add)

# ylim
ymx_add <- maxlod(out_add)*1.04

# make the plots
panel_lab_adj <- c(0.12, 0.06)
panel_lab_cex <- 1.3

res <- 256
png("../Figs/lod_curves.png", height=7.5*res, width=6*res, pointsize=14, res=res)
layout(rbind(1,1,1,1,1,2,2,2,3,3))
par(mar=c(2.1, 4.1, 1.6, 1.1))

plot(out_add, pmap, xlab="", ylim=c(0, ymx_add), altcol=altcolor)


par(mar=c(0.5,4.1,1.6,1.1))
ymx <- max(abs(blup[,1:8]))*1.04 * 1.6
mgp <- c(2.1, 0.3, 0)
plot_coefCC(blup, pmap, xaxt="n", ylim=c(-ymx, ymx), xlab="", mgp=mgp)
legend("topleft", ncol=4, lwd=2, col=CCcolors, legend=names(CCcolors), bg="gray92")
par(mar=c(3.1,4.1,0,1.1))
plot(out_add, pmap[1], xlab="", xaxt="n", mgp=mgp)
axis(side=1, at=pretty(par("usr")[1:2]), tick=FALSE, mgp=c(0, 0.2, 0))
title(xlab="Chr 1 position (Mbp)", mgp=c(1.8, 0, 0))
dev.off()
