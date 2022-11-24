# COMPOSITION FILES

results
- BasicInfo.R > basic information from the datasets
- DADA2.R > dada2 its workflow pipeline NovaSeq data
- RarefactionCurve.R > Rarefaction curve script metabarcoding soil data from ±Leiden
- Testdata3Uitwerkingen.docx > Results from test Novaseq data (50 random samples).  
                              Code plus printout from dada2 its workflow pipeline on test data.
- *.png > results from error models tested on the NovaSeq data.

data
- ARISE_Sample information_logbook (2).xlsx > sample information from the soil data (locations ect)
- ASVtab_raw.csv > csv data from the soil data
- OTU97tab_tax.csv > OTU table from the soil data
- Primers.xlsx > Primers used by BaseClear for the NovaSeq data

# PART 1

DADA2 ITS pipeline workflow for Novaseq fungi data
This pipeline runs the dada2 its workflow for NovaSeq paired illumina data

Basic information from the data:
-> even uit R script halen en hierin samenvatten

R packages:
- BiocManager (: BiocManager::install("dada2", version = "3.8"))
- dada2
- ShortRead
- Biostrings
- magrittr
- dplyr
Other installations:
- cutadapt (anaconda: )

Workflow:
-> stappenplan uitwerken

