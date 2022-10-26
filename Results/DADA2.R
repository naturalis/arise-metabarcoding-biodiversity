library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

#path <- "/Users/winnythoen/Library/CloudStorage/GoogleDrive-winnywinbin@gmail.com/.shortcut-targets-by-id/1E5lAF2Q142BBy6yXw-NM1VRV6_5N3le3/21022-375801_ARISE_Soil_Pilot/BaseClear/145561_21022_NovaseqSoilMetabarcoding _nafiesa/raw_sequences"
path <- "/Users/winnythoen/Desktop/BioInformatica/Afstuderen/Naturalis/Project/Data/TestData"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_", full.names = TRUE))

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Install cutadapt for removing the primers: http://cutadapt.readthedocs.io/en/stable/index.html

errF <- learnErrors(fnFs, multithread = TRUE)