library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

#path <- "/Users/winnythoen/Library/CloudStorage/GoogleDrive-winnywinbin@gmail.com/.shortcut-targets-by-id/1E5lAF2Q142BBy6yXw-NM1VRV6_5N3le3/21022-375801_ARISE_Soil_Pilot/BaseClear/145561_21022_NovaseqSoilMetabarcoding _nafiesa/raw_sequences"
path <- "/Users/winnythoen/Desktop/BioInformatica/Afstuderen/TestData"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_", full.names = TRUE))

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Install cutadapt for removing the primers: http://cutadapt.readthedocs.io/en/stable/index.html

errF <- learnErrors(fnFs, multithread = TRUE)

cutadapt <- "/Users/winnythoen/opt/anaconda3/bin/cutadapt"
system2(cutadapt, args = "--version")

# Identify the primers

FWD1 <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTCCTSCGCTTATTGATATGC"
FWD2 <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCTAGACTCGTCAACGATGAAGAACGCAG"
FWD3 <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCTAGACTCGTCACCGATGAAGAACGCAG"
FWD4 <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCTAGACTCGTCATCGATGAAGAACGTAG"
FWD5 <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCTAGACTCGTCATCGATGAAGAACGTGG"

allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}

FWD1.orients <- allOrients(FWD1)
FWD2.orients <- allOrients(FWD2)
FWD3.orients <- allOrients(FWD3)
FWD4.orients <- allOrients(FWD4)
FWD5.orients <- allOrients(FWD5)
FWD1.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
# Hier krijg ik nu een error, misschien opnieuw na het verwijderen van de primers

primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnFs[[1]]),
      FWD1.ReverseReads = sapply(FWD1.orients, primerHits, fn = fnRs[[1]]),
      FWD2.ForwardReads = sapply(FWD2.orients, primerHits, fn = fnFs[[1]]),
      FWD2.ReverseReads = sapply(FWD2.orients, primerHits, fn = fnRs[[1]]),
      FWD3.ForwardReads = sapply(FWD3.orients, primerHits, fn = fnFs[[1]]),
      FWD3.ReverseReads = sapply(FWD3.orients, primerHits, fn = fnRs[[1]]),
      FWD4.ForwardReads = sapply(FWD4.orients, primerHits, fn = fnFs[[1]]),
      FWD4.ReverseReads = sapply(FWD4.orients, primerHits, fn = fnRs[[1]]),
      FWD5.ForwardReads = sapply(FWD5.orients, primerHits, fn = fnFs[[1]]),
      FWD5.ReverseReads = sapply(FWD5.orients, primerHits, fn = fnRs[[1]]))
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnFs[[1]]),
      FWD1.ReverseReads = sapply(FWD1.orients, primerHits, fn = fnRs[[1]]))
