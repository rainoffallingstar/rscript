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
datatable(XenaData)

XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "methylation") %>%
  XenaFilter(filterDatasets= "KICH") -> met_kich

# 可选的hostname包括tcgaHub\publicHub\gdcHub\icgcHub\toilHub\pcawgHub\
# atacseqHub\singlecellHub
# 可选的filterdataset包括copynumber、clinical、genomic\mutation\expression\methylation\
# pathway\RPPA\ gene

XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "gene") %>%
  XenaFilter(filterDatasets= "KICH") -> exp_kich
XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "genomic") %>%
  XenaFilter(filterDatasets= "KICH") -> exp_kich2
exp_kich2 #genomicegment
exp_kich
met_kich
options(use_hiplot = TRUE)
XenaQuery(exp_kich) %>%
  XenaDownload() -> xe_download_exp
cli = XenaPrepare(xe_download)
cli_exp = XenaPrepare(xe_download_exp)
str(cli_exp)
problems()
class(cli)
datatable(cli)
str(cli)

XenaScan(pattern="Blood")

# --------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#                   Data.category: Copy number variation aligned to hg38
#-----------------------------------------------------------------------------
query <- GDCquery(
  project = "TCGA-ACC",
  data.category = "Copy Number Variation",
  data.type = "Copy Number Segment",
  barcode = c( "TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01")
)

GDCdownload(query)
data <- GDCprepare(query)

query <- GDCquery(
  project = "TCGA-ACC",
  data.category = "Copy Number Variation",
  data.type = "Masked Copy Number Segment",
  sample.type = c("Primary solid Tumor")
) # see the barcodes with getResults(query)$cases
GDCdownload(query)
data <- GDCprepare(query)
library(SummarizedExperiment)

# Load object from TCGAWorkflowData package
# THis object will be created in the further sections,
data(GBMIllumina_HiSeq) 

# get expression matrix
data <- assay(gbm.exp)
datatable(
  data = data[1:10,], 
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = TRUE
)

# get genes information
genes.info <- rowRanges(gbm.exp)
genes.info

# get sample information
sample.info <- colData(gbm.exp)
datatable(
  data = as.data.frame(sample.info), 
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)

# get indexed clinical patient data for GBM samples
gbm_clin <- GDCquery_clinic(project = "TCGA-GBM", type = "Clinical")

# get indexed clinical patient data for LGG samples
lgg_clin <- GDCquery_clinic(project = "TCGA-LGG", type = "Clinical")

# Bind the results, as the columns might not be the same,
# we will will plyr rbind.fill, to have all columns from both files
clinical <- plyr::rbind.fill(gbm_clin,lgg_clin)

datatable(clinical[1:10,], options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

# Fetch clinical data directly from the clinical XML files.
# if barcode is not set, it will consider all samples.
# We only set it to make the example faster
query <- GDCquery(
  project = "TCGA-GBM",
  file.type = "xml",
  data.category = "Clinical",
  barcode = c("TCGA-08-0516","TCGA-02-0317")
) 
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")
datatable(clinical, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
datatable(clinical.drug, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
datatable(clinical.radiation, options = list(scrollX = TRUE,  keys = TRUE), rownames = FALSE)

clinical.admin <- GDCprepare_clinic(query, clinical.info = "admin")
datatable(clinical.admin, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")
data(mafMutect2LGGGBM)
datatable(LGGmut[1:10,], options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

gbm.subtypes <- TCGAquery_subtype(tumor = "gbm")
datatable(gbm.subtypes[1:10,], options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

library(RTCGAToolbox)
# Get the last run dates
lastRunDate <- getFirehoseRunningDates()[1]

# get DNA methylation data, RNAseq2 and clinical data for GBM
gbm.data <- getFirehoseData(
  dataset = "GBM",
  runDate = lastRunDate, 
  gistic2Date = getFirehoseAnalyzeDates(1),
  Methylation = FALSE,  
  clinical = TRUE,
  RNASeq2GeneNorm  = FALSE, 
  Mutation = TRUE,
  fileSizeLimit = 10000
)

gbm.mut <- getData(gbm.data,"Mutation")
gbm.clin <- getData(gbm.data,"clinical")

# Download GISTIC results
lastanalyzedate <- getFirehoseAnalyzeDates(1)
gistic <- getFirehoseData("GBM",GISTIC = TRUE, gistic2Date = lastanalyzedate)

# get GISTIC results
gistic.allbygene <- getData(gistic, type = "GISTIC", platform = "AllByGene")
gistic.thresholedbygene <- getData(gistic, type = "GISTIC", platform = "ThresholdedByGene")
data(GBMGistic)
gistic.allbygene %>% head() %>% gt::gt()

gistic.thresholedbygene %>% head() %>% gt::gt()

# Genomic analysis

library(TCGAbiolinks)
# Select common CN technology available for GBM and LGG
#############################
## CNV data pre-processing ##
#############################
query.gbm.nocnv <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Copy number variation",
  legacy = TRUE,
  file.type = "nocnv_hg19.seg",
  sample.type = c("Primary solid Tumor")
)
# to reduce time we will select only 20 samples
query.gbm.nocnv$results[[1]] <- query.gbm.nocnv$results[[1]][1:20,]

GDCdownload(query.gbm.nocnv, files.per.chunk = 100)

gbm.nocnv <- GDCprepare(query.gbm.nocnv, save = TRUE, save.filename = "GBMnocnvhg19.rda")

library(maftools)
# Download Mutation Annotation Format (MAF) files
LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")
GBMmut <- GDCquery_Maf(tumor = "GBM", pipelines = "mutect2")

# Merge them 
mut <- plyr::rbind.fill(LGGmut, GBMmut)
save(maf,file ="mafMutect2LGGGBM.rda", compress = "xz")

library(maftools)
# recovering data from TCGAWorkflowData package.
data(mafMutect2LGGGBM)

# To prepare for maftools we will also include clinical data
# For a mutant vs WT survival analysis 
# get indexed clinical patient data for GBM samples
gbm_clin <- GDCquery_clinic(project = "TCGA-GBM", type = "Clinical")
# get indexed clinical patient data for LGG samples
lgg_clin <- GDCquery_clinic(project = "TCGA-LGG", type = "Clinical")
# Bind the results, as the columns might not be the same,
# we will will plyr rbind.fill, to have all columns from both files
clinical <- plyr::rbind.fill(gbm_clin,lgg_clin)
colnames(clinical)[1] <- "Tumor_Sample_Barcode"

# we need to create a binary variable 1 is dead 0 is not dead
plyr::count(clinical$vital_status)
clinical$Overall_Survival_Status <- 1 # dead
clinical$Overall_Survival_Status[which(clinical$vital_status != "Dead")] <- 0

# If patient is not dead we don't have days_to_death (NA)
# we will set it as the last day we know the patient is still alive
clinical$time <- clinical$days_to_death
clinical$time[is.na(clinical$days_to_death)] <- clinical$days_to_last_follow_up[is.na(clinical$days_to_death)]

# Create object to use in maftools
maf <- read.maf(maf = mut, clinicalData = clinical, isTCGA = TRUE)
plotmafSummary(
  maf = maf,
  rmOutlier = TRUE,
  addStat = 'median',
  dashboard = TRUE
)

oncoplot(
  maf = maf,
  top = 20,
  legendFontSize = 8,
  clinicalFeatures = c("tissue_or_organ_of_origin")
)

plot <- mafSurvival(
  maf = maf,
  genes = "TP53",
  time = 'time',
  Status = 'Overall_Survival_Status',
  isTCGA = TRUE
)

# 转录分析

query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "results", 
  sample.type = c("Primary solid Tumor"),
  legacy = TRUE
)
# We will use only 20 samples to make the example faster
query$results[[1]] <-  query$results[[1]][1:20,]                  
GDCdownload(query)
gbm.exp <- GDCprepare(
  query = query, 
  save = TRUE, 
  summarizedExperiment = TRUE, 
  save.filename = "GBMIllumina_HiSeq.rda"
)

query <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "results", 
  sample.type = c("Primary solid Tumor"),
  legacy = TRUE
)
# We will use only 20 samples to make the example faster
query$results[[1]] <-  query$results[[1]][1:20,]
GDCdownload(query)
lgg.exp <- GDCprepare(
  query = query, 
  save = TRUE, 
  summarizedExperiment = TRUE, 
  save.filename = "LGGIllumina_HiSeq.rda"
)

data("LGGIllumina_HiSeq")
data("GBMIllumina_HiSeq")

dataPrep_LGG <- TCGAanalyze_Preprocessing(
  object = lgg.exp,
  cor.cut = 0.6,    
  datatype = "raw_count",
  filename = "LGG_IlluminaHiSeq_RNASeqV2.png"
)

dataPrep_GBM <- TCGAanalyze_Preprocessing(
  object = gbm.exp,
  cor.cut = 0.6, 
  datatype = "raw_count",
  filename = "GBM_IlluminaHiSeq_RNASeqV2.png"
)

dataNorm <- TCGAanalyze_Normalization(
  tabDF = cbind(dataPrep_LGG, dataPrep_GBM),
  geneInfo = TCGAbiolinks::geneInfo,
  method = "gcContent"
) #18323   672

dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile",
  qnt.cut =  0.25
)  # 13742   672

save(dataFilt, file = paste0("LGG_GBM_Norm_IlluminaHiSeq.rda"))

dataFiltLGG <- subset(
  dataFilt, 
  select = substr(colnames(dataFilt),1,12) %in% lgg_clin$bcr_patient_barcode
)

dataFiltGBM <- subset(
  dataFilt, 
  select = substr(colnames(dataFilt),1,12) %in% gbm_clin$bcr_patient_barcode
)

dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFiltLGG,
  mat2 = dataFiltGBM,
  Cond1type = "LGG",
  Cond2type = "GBM",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT"
)
# Number of differentially expressed genes (DEG)
nrow(dataDEGs)

# EA: enrichment analysis

#-------------------  4.2 EA: enrichment analysis             --------------------
ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes LGG Vs GBM", 
  RegulonList = rownames(dataDEGs)
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOBPTab = ansEA$ResBP,
  nRGTab = rownames(dataDEGs),
  nBar = 20)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOCCTab = ansEA$ResCC,
  nRGTab = rownames(dataDEGs),
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOMFTab = ansEA$ResMF,
  nRGTab = rownames(dataDEGs),
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  PathTab = ansEA$ResPat,
  nRGTab = rownames(dataDEGs),
  nBar = 20
)

# PEA: Pathways enrichment analysis

library(SummarizedExperiment)
GenelistComplete <- rownames(assay(lgg.exp,1))

# DEGs TopTable
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = dataDEGs,
  typeCond1 = "LGG",
  typeCond2 = "GBM",
  TableCond1 = dataFilt[,colnames(dataFiltLGG)],
  TableCond2 = dataFilt[,colnames(dataFiltGBM)]
)

dataDEGsFiltLevel$GeneID <- 0

library(clusterProfiler)
# Converting Gene symbol to geneID
eg = as.data.frame(
  bitr(dataDEGsFiltLevel$mRNA,
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = "org.Hs.eg.db")
)
eg <- eg[!duplicated(eg$SYMBOL),]

dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]

dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]

# table(eg$SYMBOL == dataDEGsFiltLevel$mRNA) should be TRUE
all(eg$SYMBOL == dataDEGsFiltLevel$mRNA)
dataDEGsFiltLevel$GeneID <- eg$ENTREZID

dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
library(pathview)
# pathway.id: hsa05214 is the glioma pathway
# limit: sets the limit for gene expression legend and color
hsa05214 <- pathview::pathview(
  gene.data  = genelistDEGs,
  pathway.id = "hsa05214",
  species    = "hsa",
  limit = list(gene = as.integer(max(abs(genelistDEGs))))
)

# Inference of gene regulatory networks
### read biogrid info (available in TCGAWorkflowData as "biogrid")
### Check last version in https://thebiogrid.org/download.php 
file <- "http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.146/BIOGRID-ALL-3.4.146.tab2.zip"
if(!file.exists(gsub("zip","txt",basename(file)))){
  downloader::download(file,basename(file))
  unzip(basename(file),junkpaths =TRUE)
}

tmp.biogrid <- read.csv(gsub("zip","txt",basename(file)), header=TRUE, sep="\t", stringsAsFactors=FALSE)
### plot details (colors & symbols)
mycols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628')

### load network inference libraries
library(minet)
library(c3net)

### deferentially identified genes using TCGAbiolinks
# we will use only a subset (first 50 genes) of it to make the example faster
names.genes.de <- rownames(dataDEGs)[1:30]

data(biogrid)
net.biogrid.de <- getAdjacencyBiogrid(tmp.biogrid, names.genes.de)

mydata <- dataFiltLGG[names.genes.de, ]

### infer networks
t.mydata <- t(mydata)
net.aracne <- minet(t.mydata, method = "aracne")
net.mrnet <- minet(t.mydata)
net.clr <- minet(t.mydata, method = "clr")
net.c3net <- c3net(mydata)

### validate compared to biogrid network
tmp.val <- list(
  validate(net.aracne, net.biogrid.de), 
  validate(net.mrnet, net.biogrid.de),
  validate(net.clr, net.biogrid.de), 
  validate(net.c3net, net.biogrid.de)
)

### plot roc and compute auc for the different networks
dev1 <- show.roc(tmp.val[[1]],cex=0.3,col=mycols[1],type="l")
res.auc <- auc.roc(tmp.val[[1]])
for(count in 2:length(tmp.val)){
  show.roc(tmp.val[[count]],device=dev1,cex=0.3,col=mycols[count],type="l")
  res.auc <- c(res.auc, auc.roc(tmp.val[[count]]))
}

legend("bottomright", legend=paste(c("aracne","mrnet","clr","c3net"), signif(res.auc,4), sep=": "),
       col=mycols[1:length(tmp.val)],lty=1, bty="n" )
# Please, uncomment this line to produce the pdf files.
# dev.copy2pdf(width=8,height=8,device = dev1, file = paste0("roc_biogrid_",cancertype,".pdf"))

# Epigenetic analysis
#----------------------------
# Obtaining DNA methylation
#----------------------------
# Samples
lgg.samples <- matchedMetExp("TCGA-LGG", n = 10)
gbm.samples <- matchedMetExp("TCGA-GBM", n = 10)
samples <- c(lgg.samples,gbm.samples)

#-----------------------------------
# 1 - Methylation
# ----------------------------------
# For methylation it is quicker in this case to download the tar.gz file
# and get the samples we want instead of downloading files by files
query <- GDCquery(
  project = c("TCGA-LGG","TCGA-GBM"),
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450",
  legacy = TRUE, 
  barcode = samples
)
GDCdownload(query)
met <- GDCprepare(query, save = FALSE)

# We will use only chr9 to make the example faster
met <- subset(met,subset = as.character(seqnames(met)) %in% c("chr9"))
# This data is avaliable in the package (object elmerExample)
data(elmerExample)
#----------------------------
# Mean methylation
#----------------------------
# Plot a barplot for the groups in the disease column in the
# summarizedExperiment object

# remove probes with NA (similar to na.omit)
met <- met[rowSums(is.na(assay(met))) == 0,]

df <- data.frame(
  "Sample.mean" = colMeans(assay(met), na.rm = TRUE),
  "groups" = met$project_id
)

library(ggpubr)
ggpubr::ggboxplot(
  data = df,
  y = "Sample.mean",
  x = "groups",
  color = "groups",
  add = "jitter",
  ylab = expression(paste("Mean DNA methylation (", beta, "-values)")),
  xlab = ""
) + stat_compare_means() 
#------- Searching for differentially methylated CpG sites     ----------
dmc <- TCGAanalyze_DMC(
  data = met,
  groupCol = "project_id", # a column in the colData matrix
  group1 = "TCGA-GBM", # a type of the disease type column
  group2 = "TCGA-LGG", # a type of the disease column
  p.cut = 0.05,
  diffmean.cut = 0.15,
  save = FALSE,
  legend = "State",
  plot.filename = "LGG_GBM_metvolcano.png",
  cores = 1 # if set to 1 there will be a progress bar
)

#--------------------------
# DNA Methylation heatmap
#-------------------------
library(ComplexHeatmap)
clinical <- plyr::rbind.fill(gbm_clin,lgg_clin)

# get the probes that are Hypermethylated or Hypomethylated
# met is the same object of the section 'DNA methylation analysis'
status.col <- "status"
probes <- rownames(dmc)[grep("hypo|hyper",dmc$status,ignore.case = TRUE)]
sig.met <- met[probes,]


# top annotation, which sampples are LGG and GBM
# We will add clinical data as annotation of the samples
# we will sort the clinical data to have the same order of the DNA methylation matrix
clinical.order <- clinical[match(substr(colnames(sig.met),1,12),clinical$bcr_patient_barcode),]

ta <- HeatmapAnnotation(
  df = clinical.order[, c("disease", "gender", "vital_status", "race")],
  col = list(
    disease = c("LGG" = "grey", "GBM" = "black"),
    gender = c("male" = "blue", "female" = "pink")
  )
)

# row annotation: add the status for LGG in relation to GBM
# For exmaple: status.gbm.lgg Hypomethyated means that the
# mean DNA methylation of probes for lgg are hypomethylated
# compared to GBM ones.
ra = rowAnnotation(
  df = dmc[probes, status.col],
  col = list(
    "status.TCGA.GBM.TCGA.LGG" =
      c("Hypomethylated" = "orange",
        "Hypermethylated" = "darkgreen")
  ),
  width = unit(1, "cm")
)

heatmap  <- Heatmap(
  matrix = assay(sig.met),
  name = "DNA methylation",
  col = matlab::jet.colors(200),
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  bottom_annotation = ta,
  column_title = "DNA Methylation"
) 
# Save to pdf
png("heatmap.png",width = 600, height = 400)
draw(heatmap, annotation_legend_side =  "bottom")
dev.off()

library(rGADEM)
library(BSgenome.Hsapiens.UCSC.hg19)
library(motifStack)
library(SummarizedExperiment)
library(dplyr)

probes <- rowRanges(met)[rownames(dmc)[grep("hypo|hyper",dmc$status,ignore.case = TRUE)],]

# Get hypo/hyper methylated probes and make a 200bp window 
# surrounding each probe.
sequence <- GRanges(
  seqnames = as.character(seqnames(probes)),
  IRanges(
    start = ranges(probes) %>% as.data.frame() %>% dplyr::pull("start") - 100,
    end = ranges(probes) %>% as.data.frame() %>% dplyr::pull("end") + 100), 
  strand = "*"
)
#look for motifs
gadem <- GADEM(sequence, verbose = FALSE, genome = Hsapiens)

# How many motifs were found?
nMotifs(gadem)

# get the number of occurrences
nOccurrences(gadem)

# view all sequences consensus
consensus(gadem)

# Print motif
pwm <- getPWM(gadem)
pfm  <- new("pfm",mat = pwm[[1]],name = "Novel Site 1")
plotMotifLogo(pfm)

# Number of instances of motif 1?
length(gadem@motifList[[1]]@alignList)

# Integrative (Epigenomic & Transcriptomic) analysis

library(ChIPseeker)
library(pbapply)
library(ggplot2)

#------------------ Working with ChipSeq data ---------------
# Step 1: download histone marks for a brain and non-brain samples.
#------------------------------------------------------------
# loading annotation hub database
library(AnnotationHub)
ah = AnnotationHub()

# Searching for brain consolidated epigenomes in the roadmap database
bpChipEpi_brain <- query(ah , c("EpigenomeRoadMap", "narrowPeak", "chip", "consolidated","brain","E068"))
# Get chip-seq data
histone.marks <- pblapply(names(bpChipEpi_brain), function(x) {ah[[x]]})
names(histone.marks) <- names(bpChipEpi_brain) 
# OBS: histone.marks is available in TCGAWorkflowData package

data(histoneMarks)
# Create a GR object based on the hypo/hypermethylated probes.
probes <- keepStandardChromosomes(rowRanges(met)[rownames(dmc)[dmc$status %in% c("Hypermethylated in TCGA-GBM", "Hypomethylated in TCGA-GBM")],])
# Defining a window of 3kbp - 3kbp_probe_3kbp
# to make it work with ChIPseeker package version "1.31.3.900"
attributes(probes)$type <- "start_site"
attributes(probes)$downstream <- 3000
attributes(probes)$upstream <- 3000
probes <- GenomicRanges::resize(probes,6001,fix = "center") 

### Profile of ChIP peaks binding to TSS regions
# First of all, to calculate the profile of ChIP peaks binding to TSS regions, we should
# prepare the TSS regions, which are defined as the flanking sequence of the TSS sites.
# Then align the peaks that are mapping to these regions and generate the tagMatrix.
tagMatrixList <- pbapply::pblapply(histone.marks, function(x) {
  getTagMatrix(keepStandardChromosomes(x), windows = probes, weightCol = "score")
})

# change names retrieved with the following command: basename(bpChipEpi_brain$title)
names(tagMatrixList) <- c("H3K4me1","H3K4me3", "H3K9ac", "H3K9me3", "H3K27ac",  "H3K27me3", "H3K36me3")

# To plot the enrichment heatmap use the function tagHeatmap
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000),color = NULL)

p <- plotAvgProf(tagMatrixList, xlim = c(-3000,3000), xlab = "Genomic Region (5'->3', centered on CpG)")
# We are centreing in the CpG instead of the TSS. So we'll change the labels manually
p <- p + scale_x_continuous(breaks=c(-3000,-1500,0,1500,3000),labels=c(-3000,-1500,"CpG",1500,3000))
library(ggthemes)
p + theme_few() + scale_colour_few(name="Histone marks") +  guides(colour = guide_legend(override.aes = list(size=4)))

# Enhancer Linking by Methylation/Expression Relationship

#----------- 8.3 Identification of Regulatory Enhancers   -------
library(TCGAbiolinks)
# Samples: primary solid tumor w/ DNA methylation and gene expression
lgg.samples <- matchedMetExp("TCGA-LGG", n = 10)
gbm.samples <- matchedMetExp("TCGA-GBM", n = 10)
samples <- c(lgg.samples,gbm.samples)

#-----------------------------------
# 1 - Methylation
# ----------------------------------
query.met <- GDCquery(
  project = c("TCGA-LGG","TCGA-GBM"),
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450",
  legacy = TRUE, 
  barcode = samples
)
GDCdownload(query.met)
met <- GDCprepare(query.met, save = FALSE)
met <- subset(met,subset = as.character(GenomicRanges::seqnames(met)) %in% c("chr9"))

#-----------------------------------
# 2 - Expression
# ----------------------------------
query.exp <- GDCquery(
  project = c("TCGA-LGG","TCGA-GBM"),
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "results", 
  legacy = TRUE, 
  barcode =  samples
)
GDCdownload(query.exp)
exp <- GDCprepare(query.exp, save = FALSE)
save(exp, met, gbm.samples, lgg.samples, file = "elmer.example.rda", compress = "xz")

library(ELMER)
library(MultiAssayExperiment)
library(GenomicRanges)
distal.probes <- get.feature.probe(genome = "hg19", 
                                   met.platform = "450K")

# Recover the data created in the last step
data(elmerExample)
rownames(exp) <- values(exp)$ensembl_gene_id
mae <- createMAE(
  exp = assay(exp), 
  met = met,
  save = TRUE,
  linearize.exp = TRUE,
  save.filename = "mae.rda",
  filter.probes = distal.probes,
  met.platform = "450K",
  genome = "hg19",
  TCGA = TRUE
)
mae

cores <- 1 # you can use more cores if you want
group.col <- "project_id"
group1 <- "TCGA-GBM"
group2 <- "TCGA-LGG"

# Available directions are hypo and hyper, we will use only hypo
# due to speed constraint
direction <- "hypo"
dir.out <- paste0("elmer/",direction)
dir.create(dir.out, showWarnings = FALSE, recursive = TRUE)

#--------------------------------------
# STEP 3: Analysis                     |
#--------------------------------------
# Step 3.1: Get diff methylated probes |
#--------------------------------------
message("Get diff methylated probes")
Sig.probes <- get.diff.meth(
  data = mae, 
  group.col = group.col,
  group1 = group1,
  group2 =  group2,
  minSubgroupFrac = 1.0, # Use all samples
  sig.dif = 0.2, # defualt is 0.3
  diff.dir = direction, # Search for hypomethylated probes in group 1
  cores = cores, 
  dir.out = dir.out, 
  pvalue = 0.1
)
datatable(
  data = Sig.probes[1:10,], 
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = TRUE
)

# Identification of putative target gene(s) for differentially methylated distal probes
#-------------------------------------------------------------
# Step 3.2: Identify significant probe-gene pairs            |
#-------------------------------------------------------------
# Collect nearby 20 genes for Sig.probes
message("Get nearby genes")
nearGenes <- GetNearGenes(
  data = mae,
  numFlankingGenes = 4, # default is 20 genes
  probes = Sig.probes$probe
)
length(Sig.probes$probe)
dim(nearGenes)
head(nearGenes)

message("Get anti correlated probes-genes")
pair <- get.pair(
  data = mae,
  group.col = group.col,
  group1 = group1,
  group2 =  group2,
  nearGenes = nearGenes,
  mode = "supervised",
  minSubgroupFrac = 1, # % of samples to use in to create groups U/M
  raw.pvalue = 0.1,   # defualt is 0.001
  Pe = 0.5, # Please set to 0.001 to get significant results
  filter.probes = TRUE, # See preAssociationProbeFiltering function
  filter.percentage = 0.05,
  save = FALSE, # Create CVS file
  filter.portion = 0.3,
  dir.out = dir.out,
  diff.dir = direction,
  cores = cores,
  label = direction
)
datatable(
  pair[1:10,], 
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = TRUE
)

# Motif enrichment analysis

#-------------------------------------------------------------
# Step 3.3: Motif enrichment analysis on the selected probes |
#-------------------------------------------------------------
enriched.motif <- get.enriched.motif(
  data = mae,
  probes = pair$Probe, 
  dir.out = dir.out, 
  label = direction,
  pvalue = 1, # default is FDR < 0.05
  min.incidence = 10,
  lower.OR = 1.1
)
# One of the output from the  previous function is a file with the motif, OR and Number of probes
# It will be used for plotting purposes
motif.enrichment <- read.csv(paste0(dir.out,"/getMotif.",direction, ".motif.enrichment.csv"))
head(enriched.motif[names(enriched.motif)[1]]) ## probes in the given set that have the first motif.
motif.enrichment %>% head %>% gt::gt()

#-------------------------------------------------------------
# Step 3.4: Identifying regulatory TFs                        |
#-------------------------------------------------------------
TF <- get.TFs(
  data = mae, 
  group.col = group.col,
  group1 = group1,
  group2 =  group2,
  mode = "supervised",
  enriched.motif = enriched.motif,
  dir.out = dir.out, 
  cores = cores,
  save.plots = FALSE,
  diff.dir = direction,
  label = direction
)

# One of the output from the previous function is a file with the raking of TF,
# for each motif. It will be used for plotting purposes
TF.meth.cor <- get(load(paste0(dir.out,"/getTF.",direction,".TFs.with.motif.pvalue.rda")))
datatable(TF, 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

datatable(TF.meth.cor[1:10,1:6], 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = TRUE)

# ELMER visualization functions

heatmapPairs(
  data = mae, 
  group.col = group.col,
  group1 = group1, 
  group2 = group2, 
  annotation.col = c("gender"),
  pairs = pair,
  filename =  NULL
)

motif.enrichment.plot(
  motif.enrichment = motif.enrichment, 
  save = FALSE,
  significant = list(lowerOR = 1.1)
) # Filter motifs in the plot lowerOR > 1.3

grid:TF.rank.plot(motif.pvalue=TF.meth.cor, motif=TF$motif[1], save=FALSE)

png("TF.png",width = 800, height = 400)
scatter.plot(
  data = mae, 
  category = group.col, 
  save = FALSE, 
  lm_line = TRUE,
  byTF = list(
    TF = unlist(stringr::str_split(TF[1,"top_5percent_TFs"],";"))[1:4], 
    probe = enriched.motif[[TF$motif[1]]]
  )
)
dev.off()









  