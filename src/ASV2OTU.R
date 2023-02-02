library(Biostrings)
library(dplyr)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocParallel")
library(tibble)
library(dplyr)
library(DECIPHER)
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

nproc <- 4 # set to number of cpus/processors to use for the clustering
asv_sequences <- colnames(seqtab.nochim)
sample_names <- rownames(seqtab.nochim)
dna <- Biostrings::DNAStringSet(asv_sequences)

## Find clusters of ASVs to form the new OTUs
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
clusters <- DECIPHER::TreeLine(
  myDistMatrix=d,
  method = "complete",
  cutoff = 0.03, # use `cutoff = 0.03` for a 97% OTU
  type = "clusters",
  processors = nproc)

## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
clusters <- clusters %>% add_column(sequence = asv_sequences)

merged_seqtab <- seqtab.nochim %>% 
  t %>%
  rowsum(clusters$cluster) %>%
  t
# Optional renaming of clusters to OTU<cluster #>
colnames(merged_seqtab) <- paste0("OTU", colnames(merged_seqtab))

otu_results <- as.data.frame(t(merged_seqtab))
write.csv(otu_results, "OTU_test.csv", row.names=TRUE)

