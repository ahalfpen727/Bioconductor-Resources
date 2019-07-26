## ----first_load, echo=TRUE, warning=FALSE, results='hide', message=FALSE----
library(RITANdata)
library(RITAN)

## ----geneset_overlap_d, echo=TRUE, warning=FALSE, fig.width = 7, fig.height = 7, fig.align='center'----
# Show the fraction of genes common to multiple disease-gene relationships.
o <- geneset_overlap( geneset_list$DisGeNet )
plot( density(c( o[upper.tri(o)] )), log='y', ylim = c(1e-3, 1e3),
      main='', xlab='Fraction Shared Genes', ylab='Prob()')

## ----geneset_overlap_h, echo=TRUE, warning=FALSE, fig.width = 7, fig.height = 7, fig.align='center'----
# Show the diseases and their gene-level overlap, when each disease has some overlap of at least 80%.
diag(o) <- NA # ignore self-overlap
i <- which( o > 0.8, arr.ind = TRUE )
mat <- o[ unique(i[,1]), unique(i[,2]) ]
heatmap( mat, col = rev(gray(seq(0,1,length.out = 15))),
         cexRow=.7, cexCol=0.7, margins = c(7,7) )

## ----geneset_overlap2_d, echo=TRUE, warning=FALSE, fig.width = 7, fig.height = 7, fig.align='center'----
# Show the fraction of genes common between disease-gene relationships and GO-Slim term definitions
o <- geneset_overlap( geneset_list$DisGeNet, geneset_list$GO_slim_generic )
o <- o[ , !( colnames(o) %in% c('biological_process', 'molecular_function', 'cellular_component')) ] # remove the root of each sub-ontology
plot( density(c(o)), log='y',
      main='', xlab='Fraction Shared Genes', ylab='Prob()')

## ----geneset_overlap2_h, echo=TRUE, warning=FALSE, fig.width = 7, fig.height = 7, fig.align='center'----
# Show the heatmap for relationships where a disease and term share 95% of genes in common
i <- which( o > 0.95, arr.ind = TRUE )
mat <- o[ unique(i[,1]), unique(i[,2]) ]
heatmap( mat, col = rev(gray(seq(0,1,length.out = 15))),
         cexRow=.7, cexCol=0.7, margins = c(7,7) )

## ----geneset_overlap3, echo=TRUE, warning=FALSE--------------------------
rownames(o)[ o[ , "chromosome_organization" ] > 0.66 ]

## ----geneset_overlap4, echo=TRUE, warning=FALSE--------------------------
d <- rownames(o)[ o[ , "cell_motility" ] > 0.66 ]
str(d)
new_geneset <- intersect( unique(unlist(geneset_list$DisGeNet[d])),
                          unique(unlist(geneset_list$GO_slim_generic$cell_motility)) )
str(new_geneset)

## ----resource_reduce_1, echo=TRUE, warning=FALSE-------------------------
unique_diseases <- resource_reduce( geneset_list$DisGeNet )

## ----resource_reduce_2, echo=TRUE, warning=FALSE-------------------------
unique_disease_slim <- resource_reduce( c(geneset_list$DisGeNet, geneset_list$GO_slim_generic), min_overlap = 0.95 )

