#library(crisprDesign)
#library(GenomicRanges)
#library(GenomeInfoDb)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
genome=BSgenome.Hsapiens.UCSC.hg38
seq <- getSeq(genome, "chr12", start=1, end=171330, as.character=TRUE)
seq <- DNAStringSet(seq)
names(seq) <- "chr12"
writeXStringSet(seq, "../fasta/chr12.fa")
