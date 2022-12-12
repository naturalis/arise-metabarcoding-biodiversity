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

    install.package(..)

And BioParallel:

    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("BiocParallel")

Load all packages before working with the pipeline.

### Other packages:
For removing the primers we will need to install cutadapt.
Install cutadapt : http://cutadapt.readthedocs.io/en/stable/index.html

    python3 -m pip install --user --upgrade cutadapt

After installing cutadapt it will show the directory where it is saved, use that path for in your R script.

# Filter your data

The code used to filter my data:

    out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxEE=c(2,2), truncLen=c(240,200), truncQ=2, maxN=0, rm.phix=TRUE, minLen = 200, compress=TRUE, verbose=TRUE, multithread=TRUE)  # on windows, set multithread = FALSE

With my data I used `truncLen=c(240,200)` because my data is ±251bp per sample. 
The `maxE` parameter sets the maximum number of expected errors allowed in a read.
`truncQ` filters out low quality reads by comparing if they are less or equel to the truncQ score. 

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

### Error models

Simply said; NovaSeq data has a lot more data so a lot more errors. 
This is why we need to use a different error model.. There are a lot of people online who have tested multiple different error models. When testing the error models look at the error plots; the closer the black dots are to the black line ánd keep decreacing are the best.
I have tested 5 different error models on my data, including the default for comparison. The discussions that directed me to try these different error models; https://github.com/benjjneb/dada2/issues/791 and https://github.com/ErnakovichLab/dada2_ernakovichlab#learn-the-error-rates. You can read these for extra information. 
After trying all the different error models and comparing the plots if chose error model ?.

All 4 error models are in [FunctionsDADA2.R](FunctionsDADA2.R)
