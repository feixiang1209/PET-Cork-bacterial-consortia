# This document contains all statistical analyses conducted for the manuscript.
# Due to the stochastic and iterative nature of certain analyses, some figure
# parameters may vary slightly upon reanalysis; however, the core results remain unchanged.
# For editing purposes, certain figure colors, shapes, and layout arrangements
# have been modified without affecting the underlying results.
# All data required to reproduce the analyses can be found here:
# https://github.com/feixiang1209/PET-Cork-bacterial-consortia




library(compositions)
library(pheatmap)

rpkm <- read.table("coverm_results_all-F_genus.txt", header=TRUE, sep="\t", row.names=1)

rpkm_clr <- clr(rpkm + 1e-6)

taxonomy <- sapply(strsplit(rownames(rpkm_clr), "-"), `[`, 2)

print(taxonomy)

rpkm_clr_sorted <- rpkm_clr[order(taxonomy), ]

annotation_row <- data.frame(Taxonomy = taxonomy[order(taxonomy)])
rownames(annotation_row) <- rownames(rpkm_clr_sorted)

svg("heatmap_CLR_taxonomy_sorted.svg", width = 12, height = 20)
pheatmap(rpkm_clr_sorted,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "average",
         fontsize_row = 6,
         fontsize_col = 8,
         annotation_row = annotation_row,  
         cluster_rows = FALSE)             

dev.off()
