# load the required packages for dada2
library(dada2)
library(ShortRead)
library(Biostrings)
library(magrittr)
library(dplyr)

# Assign the path were the data is
# Assign the path for cutadapt -> install cutadapt: http://cutadapt.readthedocs.io/en/stable/index.html
path <- "/home/winny.thoen/arise-metabarcoding-biodiversity/data/raw_sequences_NovaSeq"
cutadapt <- "/home/winny.thoen/.local/bin/cutadapt"
system2(cutadapt, args = "--version")

# Import all functions
source("/home/winny.thoen/arise-metabarcoding-biodiversity/src/FunctionsDADA2.R")

list.files(path)

# We have to separate the two different reads (read1 and read2)
fnFs <- sort(list.files(path, pattern = "_R1_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_", full.names = TRUE))

# Quality plot tests per two reads. To view different read change [..:..]
# Maybe interesting for in the paper ect
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])

# Identify the primers
# The forward primers are all the same primer, except for the barcodes/adapters
FWD1 <- "CTAGACTCGTCATCGATGAAGAACGCAG"
FWD2 <- "CTAGACTCGTCAACGATGAAGAACGCAG"
FWD3 <- "CTAGACTCGTCACCGATGAAGAACGCAG"
FWD4 <- "CTAGACTCGTCATCGATGAAGAACGTAG"
FWD5 <- "CTAGACTCGTCATCGATGAAGAACGTGG"
REV <- "TCCTSCGCTTATTGATATGC"

# All the primers (in all orders) are saved in .orients
# To view: //.orients
FWD1.orients <- allOrients(FWD1)
FWD2.orients <- allOrients(FWD2)
FWD3.orients <- allOrients(FWD3)
FWD4.orients <- allOrients(FWD4)
FWD5.orients <- allOrients(FWD5)
REV.orients <- allOrients(REV)

# Here we make a new path where we can save the filtered reads.
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
# We filter the reads for other characters than A/G/T/C
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# Here we make a table of the primers found and how many
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD1.ReverseReads = sapply(FWD1.orients, primerHits, fn = fnRs.filtN[[1]]),
      FWD2.ForwardReads = sapply(FWD2.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD2.ReverseReads = sapply(FWD2.orients, primerHits, fn = fnRs.filtN[[1]]),
      FWD3.ForwardReads = sapply(FWD3.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD3.ReverseReads = sapply(FWD3.orients, primerHits, fn = fnRs.filtN[[1]]),
      FWD4.ForwardReads = sapply(FWD4.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD4.ReverseReads = sapply(FWD4.orients, primerHits, fn = fnRs.filtN[[1]]),
      FWD5.ForwardReads = sapply(FWD5.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD5.ReverseReads = sapply(FWD5.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Make a new path where the reads can be saved after primer removal
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD1.RC <- dada2:::rc(FWD1)
FWD2.RC <- dada2:::rc(FWD2)
FWD3.RC <- dada2:::rc(FWD3)
FWD4.RC <- dada2:::rc(FWD4)
FWD5.RC <- dada2:::rc(FWD5)
REV.RC <- dada2::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
# Only FWD1 is used, error rate default of 0.1 takes out all forward reads
R1.flags <- paste("-g", FWD1, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD1.RC) 
# Only FWD1 is used, the default error rate is 0.1 which means that the 
# difference between all different primers are covered

# Removing the primers from the reads
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# Make a table again to see if it is all removed.
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD1.ReverseReads = sapply(FWD1.orients, primerHits, fn = fnRs.cut[[1]]),
      FWD2.ForwardReads = sapply(FWD2.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD2.ReverseReads = sapply(FWD2.orients, primerHits, fn = fnRs.cut[[1]]),
      FWD3.ForwardReads = sapply(FWD3.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD3.ReverseReads = sapply(FWD3.orients, primerHits, fn = fnRs.cut[[1]]),
      FWD4.ForwardReads = sapply(FWD4.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD4.ReverseReads = sapply(FWD4.orients, primerHits, fn = fnRs.cut[[1]]),
      FWD5.ForwardReads = sapply(FWD5.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD5.ReverseReads = sapply(FWD5.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_", full.names = TRUE))

# Extract sample names, assuming filenames have format:
sample.names <- unname(sapply(cutFs, get.sample.name))
sample.namesR <- unname(sapply(cutRs, get.sample.name))
head(sample.namesR)

# Now you could look at some Quality plots:
# plotQualityProfile(cutRs[1:2])

# Filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxEE=c(2,2), truncLen=c(240,200), truncQ=2, maxN=0, rm.phix=TRUE,
                     minLen = 200, compress=TRUE, verbose=TRUE, multithread=TRUE)  # on windows, set multithread = FALSE
head(out)

# Double check for identical sample names
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Controleren op de output 
out %>% 
  data.frame() %>% 
  mutate(Samples = rownames(.),
         percent_kept = 100*(reads.out/reads.in)) %>%
  select(Samples, everything()) %>%
  summarise(min_remaining = paste0(round(min(percent_kept), 2), "%"), 
            median_remaining = paste0(round(median(percent_kept), 2), "%"),
            mean_remaining = paste0(round(mean(percent_kept), 2), "%"), 
            max_remaining = paste0(round(max(percent_kept), 2), "%"))

# Error Rates Default
#errF <- learnErrors(filtFs, multithread = TRUE)
#errR <- learnErrors(filtRs, multithread = TRUE)

#plotErrors(errF, nominalQ = TRUE)
#plotErrors(errR, nominalQ = TRUE)

# Error Rates Outcome

errR_1 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)
  
errF_1 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)

# Option 2
# check what this looks like
errF_2 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e8,
  errorEstimationFunction = loessErrfun_mod2,
  verbose = TRUE
)

errR_2 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e8,
  errorEstimationFunction = loessErrfun_mod2,
  verbose = TRUE
)

# Option 3
# check what this looks like
errF_3 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod3,
  verbose = TRUE
)

# check what this looks like
errR_3 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod3,
  verbose = TRUE
)

# Create all plots
plotErrors(errF_1, nominalQ = TRUE)
plotErrors(errR_1, nominalQ = TRUE)

plotErrors(errF_2, nominalQ = TRUE)
plotErrors(errR_2, nominalQ = TRUE)

plotErrors(errF_3, nominalQ = TRUE)
plotErrors(errR_3, nominalQ = TRUE)

# Dereplicate identical reads
# This reduces the running time

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample Inference 

dadaFs <- dada(derepFs, err = errF_2, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR_2, multithread = TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(seqtab.nochim)))

# Track reads through the pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy
unite.ref <- "/home/winny.thoen/arise-metabarcoding-biodiversity/data/sh_general_release_dynamic_29.11.2022.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Show Results:

# Optinal
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# To view primers: //.orients
# plotQualityProfile(cutRs[1:2])

# fwd_qual_plots_filt
# rev_qual_plots_filt

