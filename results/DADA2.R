# load the required packages for dada2
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

# Assign the path were the data is
path <- "/Users/winnythoen/Desktop/BioInformatica/Afstuderen/TestData"
list.files(path)

# In the files you have read 1 and read 2 (like forward and reversed reads)
# We have to separate the two different reads.
fnFs <- sort(list.files(path, pattern = "_R1_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_", full.names = TRUE))

# Quality plot tests -> takes a long time and only does a view.
# Maybe interesting for in the paper. But little unnecessary.
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Identify the primers
# Primers where given to me by Vincent Merkx
# The forward primers are all the same primer, except for the barcodes/adapters

FWD1 <- "CTAGACTCGTCATCGATGAAGAACGCAG"
FWD2 <- "CTAGACTCGTCAACGATGAAGAACGCAG"
FWD3 <- "CTAGACTCGTCACCGATGAAGAACGCAG"
FWD4 <- "CTAGACTCGTCATCGATGAAGAACGTAG"
FWD5 <- "CTAGACTCGTCATCGATGAAGAACGTGG"
REV <- "TCCTSCGCTTATTGATATGC"

# This function makes a forward, reversed, compliment and reversed compliment primer of every primer given

allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}

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

# This may not work, but it is worth a shot
fnFs.filtN <- file.path("/Users/winnythoen/Desktop/BioInformatica/Afstuderen/TestData/filtN", basename(fnFs))
fnRs.filtN <- file.path("/Users/winnythoen/Desktop/BioInformatica/Afstuderen/TestData/filtN", basename(fnRs))

# Here we detect the primers
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

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


# Install cutadapt for removing the primers: http://cutadapt.readthedocs.io/en/stable/index.html
# Find where u saved cutadapt
cutadapt <- "/Users/winnythoen/opt/anaconda3/bin/cutadapt"
system2(cutadapt, args = "--version")

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
R1.flags <- paste("-g", FWD1, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD1.RC) 

# Function for removing the primer
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
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
plotQualityProfile(cutFs[0:1])

# Make a path for filtered and cut reads.
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 1, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)  # on windows, set multithread = FALSE

head(out)
