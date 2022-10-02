# section 3 in HavardXbioconductor
BiocManager::install(c("Biobase",
                       "GEOquery",
                       "genomicsclass/GSE5859Subset",
                       "affy",
                       "hgu95acdf",
                       "genefilter",
                       "parathyroidSE",
                       "airway",
                       "pasillaBamSubset",
                       "Rsamtools",
                       "GenomicAlignments",
                       "ArrayExpress",
                       "NGScopyData",
                       "AnnotationDbi"))

library(Biobase)
library(GEOquery)

geoq <- getGEO("GSE9514")    # download a microarray dataset from GEO
names(geoq)    
e <- geoq[[1]]    # extract ExpressionSet
e

# exprs gives matrix of microarray values
dim(e)    # number of features and samples in ExpressionSet
ncol(e)
nrow(e)

exprs(e)[1:3,1:3]
head(exprs(e))[,1]    # first column
exprs(e)[1,]    # first row
exprs(e)["10000_at",]    # can also index by name
rownames(e)[1]    # row names are probe sets
dim(exprs(e))    # rows are features, columns are samples

# pData gives phenotype data (sample information)
pData(e)[1:3,1:6]
names(pData(e))
pData(e)$characteristics_ch1    # column in GEO to describe experimental state/condition
as.numeric(factor(pData(e)$characteristics_ch1))    # help see replicates of each state
dim(pData(e))    # rows of pData correspond to columns of exprs
dim(e)

# fData gives feature data (probe information)
fData(e)[1:3,1:3]
dim(fData(e))    # rows of fData correspond to rows of exprs
names(fData(e))
head(fData(e)$"Gene Symbol")
head(rownames(e))

# additional annotation tied to ExpressionSet
experimentData(e)
annotation(e)
