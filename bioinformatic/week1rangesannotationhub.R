library(AnnotationHub)
ahub <- AnnotationHub() 
qhs <- ahub %>% 
  subset(species == "Homo sapiend") %>% 
  query(c("H3K4me3","Gm12878"))

qhs
qhs <- subset(ahub,species == "Homo sapiend")

gr1 <- qhs[[2]]
gr2 <- qhs[[4]]
summary(width(gr1))
summary(width(gr2))
table(width(gr2))
peaks = gr2

qhs <- ahub %>% 
  query("RefSeq")

qhs$genome
genes <- qhs[[1]]
table(gene$name)
prom <- promoters(genes)
table(width(prom))
args(promoters)
ov <- findOverlaps(prom, peaks)
length(unique(queryHits(ov)))
length(unique(subjectHits(ov)))
# add a trick
length(subsetByOverlaps(peaks,prom,ignore.strand = TRUE))
length(subsetByOverlaps(peaks,prom,ignore.strand = TRUE))/length(peaks)
length(subsetByOverlaps(prom,peaks,ignore.strand = TRUE))/length(prom)
sum(width(Reduce(peaks,ignore.strand = TRUE)))/10^6
sum(width(Reduce(prom,ignore.strand = TRUE)))/10^6
inter <- intersect(peaks,prom,ignore.strand = TRUE)/10^6 %>% 
  sum(width())
inOut <- matrix(0,ncol=2,nrow=2) %>% 
  colnames() <- c("in","out") %>% 
  rownames() <- c("in","out") 
inOut[1,1] <- inter
inOut[1,2] <- setdiff(peaks,prom,ignore.strand = TRUE) %>% 
  sum(width())

inOut[2,1] <- sum(width(setdiff(prom,peaks,ignore.strand = TRUE)))

colSums(inOut)
rowSums(inOut)
inOut[2,2] <- 3*10^9-sum(inOut)

fisher.test(inOut)$statistic
addsRatio <- inOut[1,1] * inOut[2,2] / inOut[1,2] * inOut[2,1]
 







