library(TCGAWorkflowData)
library(TCGAWorkflow)
library(DT)
library(tidyverse)
library(UCSCXenaTools)
# ---------------------------------------------------------------------------
# 以下是使用TCGAbiolinks包进行数据下载的一个示例，由于速度太慢，故弃用
# 检索 %>% 下载 %>% 处理
# ---------------------------------------------------------------------------
library(TCGAbiolinks)

query.met.gbm <- GDCquery(
  project = "TCGA-GBM", 
  legacy = TRUE,
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450", 
  barcode = c("TCGA-76-4926-01B-01D-1481-05", "TCGA-28-5211-01C-11D-1844-05")
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
# -----------------------------------------------------------------------------
# 以下使用ucscxenatools进行数据下载及处理
# -----------------------------------------------------------------------------

# 检索所有数据库

data(XenaData)
head(XenaData)

XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "copynumber") %>%
  XenaFilter(filterDatasets= "LUNG") -> met_gbm

# 可选的hostname包括
# 可选的filterdataset包括copynumber
met_gbm
options(use_hiplot = TRUE)
XenaQuery(met_gbm) %>%
  XenaDownload() -> xe_download
cli = XenaPrepare(xe_download)
class(cli)

XenaScan(pattern="Blood")
  