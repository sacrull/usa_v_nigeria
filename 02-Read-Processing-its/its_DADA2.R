#Install packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#BiocManager::install("dada2")
#install.packages("magrittr")
#install.packages("stringr")
#install.packages("data.table")
#install.packages("seqinr")
#devtools::install_github("jolars/qualpalr")
 
#Load packages
library(dada2, verbose = FALSE)
library(stringr, verbose = FALSE)
library(data.table, verbose = FALSE)
library(ShortRead, verbose = FALSE)
library(Biostrings, verbose = FALSE)
library(seqinr, verbose = FALSE) 
library(qualpalr, verbose = FALSE)   
set.seed(12349)

#Set up path to raw files
rawpath <- "~/usa_nigeria/its/all_its" #put in path to files
fnFs <- sort(list.files(rawpath, pattern="_R1_001.fastq.gz", full.names=T))
fnRs <- sort(list.files(rawpath, pattern="_R2_001.fastq.gz", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names, 50)
paste("Number of input samples: ", length(sample.names))

#Make quality plot scores
system("mkdir img") # ignore warning
fwdqual <- plotQualityProfile(fnFs[10:25])
revqual <- plotQualityProfile(fnRs[10:25])
fwdqual
revqual

pdf(paste("img/", "forward_quality_plot.pdf", sep=""))
fwdqual
dev.off()
pdf(paste("img/", "reverse_quality_plot.pdf", sep=""))
revqual
dev.off()

#Filter out uncalled bases
fnFs.filtN <- file.path(rawpath, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(rawpath, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)

#Primer Removal
cutadapt <- as.character(system("which cutadapt", intern=T))
system("cutadapt --version")
path.cut <- file.path(rawpath, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc("CTTGGTCATTTAGAGGAAGTAA")
REV.RC <- dada2:::rc("GCTGCGTTCTTCATCGATGC")
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", "CTTGGTCATTTAGAGGAAGTAA", "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", "GCTGCGTTCTTCATCGATGC", "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("--cores=0", R1.flags, R2.flags, "-n", 2,"-o", fnFs.cut[i], "-p", fnRs.cut[i], fnFs.filtN[i], fnRs.filtN[i]))
}
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE))

#Quality Filter and Trim Reads
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, trimRight=0, truncLen=200, maxN=c(0,0), maxEE=c(4,6), rm.phix=TRUE, matchIDs=TRUE, compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$precentage_retained <- retained$reads.out/retained$reads.in * 100
retained

# Learn error rates
errF <- learnErrors(filtFs, multithread=T, random=T, nbases=1e9)
errR <- learnErrors(filtRs, multithread=T, random=T, nbases=1e9)
err.f.plt <- plotErrors(errF, nominalQ=TRUE) 
err.f.plt
pdf(paste(rawpath,"/img/", "error_plot.pdf", sep=""))
err.f.plt
dev.off()

# Dereplicating
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample infrence
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, verbose = FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, verbose = FALSE)

# Merge Sample reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)

# Generate sequence tables
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Check Length distrubtion 
table(nchar(colnames(seqtab)))

# Filter Sequences by length
seqlens <- nchar(getSequences(seqtab))
seqtab.filt <- seqtab[,seqlens <= 400]
dim(seqtab.filt)

# get histogram of length distribution after filter
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab.filt))))
len.plt <- plot(x=length.histogram[,1], y=length.histogram[,2])
pdf(paste(rawpath, "/img/", "length_hist.pdf", sep=""))
plot(x=length.histogram[,1], y=length.histogram[,2])
dev.off()

# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab.filt, method="pooled", multithread=T, verbose=T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab.filt)

# Processing summary
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochimeras")
rownames(track) <- sample.names
track
# write to file
write.table(data.frame("row_names"=rownames(track),track),"read_retention.txt", row.names=FALSE, quote=F, sep="\t")
uniquesToFasta(seqtab.nochim, "rep_set_its.fa")
# fix ASV names 
system("awk '/^>/{print \">ASV\" ++i; next}{print}' < rep_set.fa > rep_set_fix.fa")
system("mv rep_set_fix.fa rep_set.fa")
# write sequence table to file, fix ASV names
my_otu_table <- t(as.data.frame(seqtab.nochim)) 
ASV.seq <- as.character(unclass(row.names(my_otu_table))) 
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') 
colnames(seqtab.nochim) <- ASV.num 
write.table(data.frame("row_names"=rownames(seqtab.nochim),seqtab.nochim),"sequence_table_its.merged.txt", row.names=FALSE, quote=F, sep="\t")
     