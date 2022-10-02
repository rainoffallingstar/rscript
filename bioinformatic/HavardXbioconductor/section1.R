BiocManager::install("gwascat")
BiocManager::install(c("Homo.sapiens",
                       "GenomicFeatures",
                       "genomicsclass/ERBS",
                       "genomicsclass/ph525x"))

library(ERBS)
library(tidyverse)
data(GM12878)
data(HepG2)
class(HepG2)
values(HepG2)
HepG2[1:10,]
chr= seqnames(HepG2)
seqnames(HepG2) %>% as.character() %>% table()
table(chr)
table(chr)[1:24]
x = HepG2[order(HepG2),]
seqnames(x) %>% as.character()

# iranges operations
library(IRanges)
ir <- IRanges(5,10)
ir
# intra-ranges method 
shift(ir,-2)
narrow(ir,)
flank(ir,width=3)
# inter-ranges method
ir <- IRanges(start=c(3,5,17),end=c(10,8,20))
range(ir)
ir
reduce(ir)
gaps(ir)
disjoin(ir)
# granges
library(GenomicRanges)
gr <- GRanges("chrZ",IRanges(start=c(5,19),end=c(35,45)),strand="+",
              seqlengths=c(chrZ=100L))
genome(gr) <- "hg19"

shift(gr,80)
trim(shift(gr,80))
mcols(gr)
mcols(gr)$value <- NULL

gr2 <- GRanges("chr2",IRanges(11:13,51:53))

grl <- GRangesList(gr, gr2)
grl
ir <- IRanges(c(3, 8, 14, 15, 19, 34, 40),
              width = c(12, 6, 6, 15, 6, 2, 7))

par(mfrow=c(4,1), mar=c(4,2,2,2))
plotRanges(ir, xlim=c(0,60))
plotRanges(reduce(ir), xlim=c(0,60))
plotRanges(disjoin(ir), xlim=c(0,60))
plotRanges(gaps(ir), xlim=c(0,60))

# overlaps

gr1 <- GRanges("chrZ",IRanges(c(1,11,21,31,41),width=5),strand="*")

# load packages
library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(ERBS)

# load ESRRA ChIP data
data(HepG2)
data(GM12878)

browseVignettes("GenomicRanges")

# find binding sites common to both HepG2 and GM12878
?findOverlaps
# for each row in query, return overlapping row in subject
res = findOverlaps(HepG2, GM12878)
class(res)
res

# ranges from the query for which we found a hit in the subject
index = queryHits(res)
erbs = HepG2[index,]
erbs

# extract only the ranges
granges(erbs)
erbs

# define human genes
library(Homo.sapiens)
ghs = genes(Homo.sapiens)
ghs

# learn about the precede function (and related functions like follow)
?precede

# for each range in erbs, find the closest preceding range in ghs
index = precede(erbs, ghs)
ghs[index[1:3]]
erbs[1:3]    # note result is strand-aware

# distance between binding sites and nearest preceding genes
distance(erbs, ghs[index])

# find transcription start site nearest to each binding site
tssgr = resize(ghs, 1)
tssgr

# distance between binding site and nearest TSS
d = distanceToNearest(erbs, tssgr)
queryHits(d)
dists = values(d)$distance
hist(dists, nc=100, xlim=c(0,100000))

index = subjectHits(d)[dists < 1000]
index

tssgr[index,]
keytypes(Homo.sapiens)
keys = as.character(values(tssgr[index])$GENEID)
columns(Homo.sapiens)
res = select(Homo.sapiens, keys = keys,
             columns = c("SYMBOL", "GENENAME"), keytype="GENEID")
res[1:2,]

library(Biostrings)

# basics of DNAStrings
dna <- DNAString("TCGAGCAAT")    # define a DNAString
dna
length(dna)    # number of bases in a DNAString
DNAString("JQX")    # error - invalid bases
DNAString("NNNACGCGC-TTA-CGGGCTANN")    # valid sequence with unknowns and gaps
dna[4:6]    # extract a substring
as.character(dna)    # convert DNAString to character

# basics of DNAStringSets
set1 <- DNAStringSet(c("TCA", "AAATCG", "ACGTGCCTA", "CGCGCA", "GTT", "TCA"))    # define a DNAStringSet
set1
set1[2:3]    # extract subset of sequences
set1[[4]]    # extract one sequence as a single DNAString
length(set1)    # number of DNAstrings in set
width(set1)    # size of each DNAString
duplicated(set1)    # detect which sequences are duplicated
unique(set1)    # keep only unique sequences
sort(set1)

# dna
dna_seq <- DNAString("ATCGCGCGCGGCTCTTTTAAAAAAACGCTACTACCATGTGTGTCTATC")

# analyze DNAStrings
letterFrequency(dna_seq, "A")    # count A in sequence
letterFrequency(dna_seq, "GC")    # count G or C in sequence
dinucleotideFrequency(dna_seq)    # frequencies of all dinucleotides
trinucleotideFrequency(dna_seq)    # frequencies of all trinucleotides

# convert DNAStrings
reverseComplement(dna_seq)    # find reverse complement
translate(dna_seq)    # amino acid translation

eco <- DNAString("GGTTTCACCGCCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCGGTGCTGGCTGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGTCAGGTGCCCGATGCGAGGTTGTTGAAGTCGA")
eco


library(ERBS)
data(HepG2)
HepG2

# load and inspect human reference genome
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens

# extract chromosome 17 sequence
c17 = Hsapiens$chr17
c17

?getSeq
class(Hsapiens)
showMethods("getSeq")

# collection of DNA strings with ChIP-seq binding peaks
hepseq = getSeq(Hsapiens, HepG2)
length(HepG2)    # same number of sequences
width(HepG2)[1:5]    # widths match

# collection of shifted DNA strings with no relationship to binding sequences - essentially random
rhepseq = getSeq(Hsapiens, shift(hepG2, 2500))

# count occurrences of a motif in DNA sequences
mot = "TCAAGGTCA"
?vmatchPattern
vcountPattern(mot, hepseq)

# consider both forward matches and reverse complement matches 
sum(vcountPattern(mot, hepseq))    # forward pattern match
sum(vcountPattern(mot, reverseComplement(hepseq)))    # reverse pattern match

## compare motif occurrence in binding peak to random upstream sequences
# count of motifs in binding peaks
sum(vcountPattern(mot, hepseq)) +
  sum(vcountPattern(mot, reverseComplement(hepseq)))
# count of motifs in randomly selected regions of equal length
sum(vcountPattern(mot, rhepseq)) +
  sum(vcountPattern(mot, reverseComplement(rhepseq)))

# for real analysis, use MotifDb package, probabilistic binding packages like MEME and FIMO