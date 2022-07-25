library(TCGAWorkflowData)
library(TCGAWorkflow)
library(DT)
library(tidyverse)
library(UCSCXenaTools)

library(TCGAbiolinks)

query.met.gbm <- GDCquery(
  project = "TCGA-KICH", 
  legacy = TRUE,
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450"
)

query.met.gbm

GDCdownload(query.met.gbm)


met.gbm.450 <- GDCprepare(
  query = query.met.gbm,
  save = TRUE, 
  save.filename = "gbmDNAmet450k.rda",
  summarizedExperiment = TRUE
)

query.met.lgg <- GDCquery(
  project = "TCGA-LGG", 
  legacy = TRUE,
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450",
  barcode = c("TCGA-HT-7879-01A-11D-2399-05", "TCGA-HT-8113-01A-11D-2399-05")
)
GDCdownload(query.met.lgg)
met.lgg.450 <- GDCprepare(
  query = query.met.lgg,
  save = TRUE, 
  save.filename = "lggDNAmet450k.rda",
  summarizedExperiment = TRUE
)

met.lgg.450$days_to_death <- NA
met.lgg.450$year_of_death <- NA
met.gbm.lgg <- SummarizedExperiment::cbind(met.lgg.450, met.gbm.450)


query.exp.lgg <- GDCquery(
  project = "TCGA-LGG", 
  legacy = TRUE,
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type = "results",
  sample.type = "Primary solid Tumor"
)
GDCdownload(query.exp.lgg)
exp.lgg <- GDCprepare(query = query.exp.lgg, save = TRUE, save.filename = "lggExp.rda")

query.exp.gbm <- GDCquery(
  project = "TCGA-GBM", 
  legacy = TRUE,
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type = "results",
  sample.type = "Primary solid Tumor"
)
GDCdownload(query.exp.gbm)
exp.gbm <- GDCprepare(query = query.exp.gbm, save = TRUE, save.filename = "gbmExp.rda")
exp.gbm.lgg <- SummarizedExperiment::cbind(exp.lgg, exp.gbm)