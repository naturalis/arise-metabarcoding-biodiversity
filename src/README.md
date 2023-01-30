# Source code

This section contains R scripts that perform various processing and reporting steps:

- [BasicInfo.R](BasicInfo.R) - basic information from the datasets
- [DADA2.R](DADA2.R) - dada2 ITS workflow pipeline NovaSeq data
- [RarefactionCurve.R](RarefactionCurve.R) - rarefaction curve script metabarcoding soil data from ±Leiden
- [FunctionsDADA2.R](FunctionsDADA2.R)

## [DADA2.R](DADA2.R)

### Needed R packages
- dada2
- ShortRead
- Biostrings
- magrittr
- dplyr

Install packages in console using:

    install.package(..)

Last package to install BioParallel:

    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("BiocParallel")

Load all packages before working with the pipeline using `library()`, line 1-11.

### Other packages/data:
Import the 'your_data', that is saved in the 'data' directory, as path:

    path <- "/home/winny.thoen/arise-metabarcoding-biodiversity/data/your_data"

For removing the primers we will need to install cutadapt.
Install cutadapt : http://cutadapt.readthedocs.io/en/stable/index.html

    python3 -m pip install --user --upgrade cutadapt

After installing cutadapt it will show the directory where it is saved, use that path for in your R script. After loading cutadapt using the following command, it is possible to check if the importing cutadapt in R was successful using `system2(cutadapt, args = "--version")` in the console.

    cutadapt <- "/home/winny.thoen/.local/bin/cutadapt"

The DADA2 workflow uses UNITE database as a fungal taxonmic assignment. The general files are saved in 'data/UNITE_database'.

    unite.ref <- "/home/winny.thoen/arise-metabarcoding-biodiversity/data/UNITE_database/sh_general_release_dynamic_29.11.2022.fasta"

The functions of the DADA2 workflow are imported seperately due to the size of the script.

    source("/home/winny.thoen/arise-metabarcoding-biodiversity/src/FunctionsDADA2.R")

# Seperating data

The first thing we do is separate the data from R1 and R2. Change this if needed depending on your file naming [line:23,26-27]. 

# (Optional) Quality check

Thoughout the pipeline there are a couple of moments when you can make a quality plot. This is not a necessary step, but may be something you want to do depending on your research. The 'Quality Check' moments will be presented as an optional step.

    plotQualityProfile(fnFs[1:2])
    plotQualityProfile(fnRs[1:2])

# Identify primers

Identify the primers that were used for your dataset [line:34-41]. The function `allOrients()` [FunctionsDADA2.R](FunctionsDADA2.R) is used to make a list of all possible nominations for each of the primer [line:43-50].

# First Filter

Create new directory for filtered results [line:53-55]. The first filter options are given before removing the primers. Otherwise the primers might not get detected and/or removed. `maxN = 0` states that there not be more then 0 other characters (N's) detected.

    filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
    
Detected primers are shown in a table by using `primerHits` function [FunctionsDADA2.R](FunctionsDADA2.R) [line:59-71].

# Primer removal

Create new directory for results after primer removal [line:75-78]. The next lines define the parameters that cutadapt will use. As u can see we only use forward primer 1, instead of all 5 that were used. This reason is that the default error rate that is used as the primers are removed is 0.1. Which means that there can be a difference of 10% in the primer as it is detected and removed. As all 5 forward primers are almost the same, they all fit exactly in this error rate. After removing forward read one, all five forward reads are removed. We check this as we make a new table to detect the primers using the function `primerHits()` [FunctionsDADA2.R](FunctionsDADA2.R) [line:98-110]. After primer removal seperate the results again by read 1 and read2 [line:112-114].

# (Optional) Quality check Moment

# Filter and trim your data

Create new directory for results after filter and trimming [line:125-126].
The code used to filter my data:

    out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxEE=c(2,2), truncLen=c(240,200), truncQ=2, maxN=0, rm.phix=TRUE, minLen = 200, compress=TRUE, verbose=TRUE, multithread=TRUE)  # on windows, set multithread = FALSE

With my data I used `truncLen=c(240,200)` because my data is ±251bp per sample. 
The `maxE` parameter sets the maximum number of expected errors allowed in a read. `truncQ` filters out low quality reads by comparing if they are less or equel to the truncQ score. When handling your own data use this for extra information: https://rdrr.io/bioc/dada2/man/filterAndTrim.html

After filtering your data you can test your data:

    out %>% 
      data.frame() %>% 
      mutate(Samples = rownames(.),
           percent_kept = 100*(reads.out/reads.in)) %>%
    select(Samples, everything()) %>%
    summarise(min_remaining = paste0(round(min(percent_kept), 2), "%"), 
              median_remaining = paste0(round(median(percent_kept), 2), "%"),
              mean_remaining = paste0(round(mean(percent_kept), 2), "%"), 
              max_remaining = paste0(round(max(percent_kept), 2), "%"))

If scores come out very low, try to change `truncLen`. 

# Error rates

Simply said; NovaSeq data has a lot more data so a lot more errors. This is why we need to use a different error model. There are a lot of people online who have tested multiple different error models. When testing the error models look at the error plots; the closer the black dots are to the black line ánd keep decreacing are the best. The discussions that directed me to try these different error models; https://github.com/benjjneb/dada2/issues/791 and https://github.com/ErnakovichLab/dada2_ernakovichlab#learn-the-error-rates. You can read these for extra information. I have compared all 4 error models and the differences where very small, at the end I decided to use error model 1 as they had the best results. All 4 error models are in [FunctionsDADA2.R](FunctionsDADA2.R). This step takes at least an hour, at the end you can make a plot result from the error models using PlotErrors [line:170-172]. The plots will appear in a Plots.pdf file in your home directory. 
# Qualoity results 

Dereplication is to remove almost identical reads, making the results high quality and reducing the running time [line:184-186]. Sample inference is the last step for higher quality results [line:192-194], before merging the reads together [line:196-197]. After merging we can make an amplicon sequence variant table (ASV) [line:199-201], and lastly remove de chimeras [line:203-204]. To inspect the lengths of the sequences use the following command in console:

    table(nchar(getSequences(seqtab.nochim)))

# Track back

Is is possible to inspect all the results we saved from the pipeline in one table. You can do a trackback [line:206-212].

# Taxonomic assignment










