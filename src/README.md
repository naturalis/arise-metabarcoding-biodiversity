# Source code

This section contains R scripts that perform various processing and reporting steps:

- [DADA2.R](DADA2.R) - dada2 ITS workflow pipeline NovaSeq data
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

Load all packages before working with the pipeline using `library()`, line 1-12.

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

### Seperating data

The first thing we do is separate the data from R1 and R2. Change this if needed depending on your file naming [line:24,27-28]. 

### Optional) Quality check

Thoughout the pipeline there are a couple of moments when you can make a quality plot. This is not a necessary step, but may be something you want to do depending on your research. The 'Quality Check' moments will be presented as an optional step.

    plotQualityProfile(fnFs[1:2])
    plotQualityProfile(fnRs[1:2])

### Identify primers

Identify the primers that were used for your dataset [line:35-42]. The function `allOrients()` [FunctionsDADA2.R](FunctionsDADA2.R) is used to make a list of all possible nominations for each of the primer [line:44-51].

### First Filter

Create new directory for filtered results [line:54-56]. The first filter options are given before removing the primers. Otherwise the primers might not get detected and/or removed. `maxN = 0` states that there not be more then 0 other characters (N's) detected.

    filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
    
Detected primers are shown in a table by using `primerHits` function [FunctionsDADA2.R](FunctionsDADA2.R) [line:60-72].

### Primer removal

Create new directory for results after primer removal [line:76-79]. The next lines define the parameters that cutadapt will use. As u can see we only use forward primer 1, instead of all 5 that were used. This reason is that the default error rate that is used as the primers are removed is 0.1. Which means that there can be a difference of 10% in the primer as it is detected and removed. As all 5 forward primers are almost the same, they all fit exactly in this error rate. After removing forward read one, all five forward reads are removed. We check this as we make a new table to detect the primers using the function `primerHits()` [FunctionsDADA2.R](FunctionsDADA2.R) [line:99-111]. After primer removal seperate the results again by read 1 and read2 [line:113-115].

### Optional) Quality check Moment

### Filter and trim your data

Create new directory for results after filter and trimming [line:126-127].
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

### Error rates

Simply said; NovaSeq data has a lot more data so a lot more errors. This is why we need to use a different error model. There are a lot of people online who have tested multiple different error models. When testing the error models look at the error plots; the closer the black dots are to the black line ánd keep decreacing are the best. The discussions that directed me to try these different error models; https://github.com/benjjneb/dada2/issues/791 and https://github.com/ErnakovichLab/dada2_ernakovichlab#learn-the-error-rates. You can read these for extra information. I have compared all 4 error models and the differences where very small, at the end I decided to use error model 1 as they had the best results. All 4 error models are in [FunctionsDADA2.R](FunctionsDADA2.R). This step takes at least an hour, at the end you can make a plot result from the error models using PlotErrors [line:171-173]. The plots will appear in a Plots.pdf file in your home directory. 

### Quality results 

Dereplication is to remove almost identical reads, making the results high quality and reducing the running time [line:185-187]. Sample inference is the last step for higher quality results [line:192-194], before merging the reads together [line:197-198]. After merging we can make an amplicon sequence variant table (ASV) [line:200-202], and lastly remove de chimeras [line:204-205]. To inspect the lengths of the sequences use the following command in console:

    table(nchar(getSequences(seqtab.nochim)))

### Track back

Is is possible to inspect all the results we saved from the pipeline in one table. You can do a trackback [line:207-213].

### Taxonomic assignment

UNITE database is used as taxonomic assignment [line:215-220]. The results are done, and can be exported. Use the following command to save from the global environment:

    save(taxa, file = "results.Rdata")
    
 After saving to your home directory, export the file.

### Output files

- ASV table
- OTU table

To creat an ASV table, we need to combine the taxonomic assembly [Taxonomic assignment](Taxonomic assignment) with the sequence table [Quality results](Quality results) we both made earlier in the script [line:224-226].

    seqtab_results <- as.data.frame(t(seqtab.nochim))
    results <- cbind(taxa, seqtab_results)
    write.matrix(results, file = "ASVtab_raw.csv")
    
To create an OTU table, we need to...


## [FunctionsDADA2.R](FunctionsDADA2.R)

### AllOrients

This function makes a forward, reversed, compliment and reversed compliment primer of every primer given. This information is needed when trying to locate all the primers in the data.
Note: To avoid mixups check is the REV primer matches REVcomp exactly -> if so: replace REV <- REV.orient[[“RevComp”]] before proceeding with the rest of the pipeline.

### primerHits
This function detects the primers in your data. Outcome of this function is shown after creating a rbind table [line:59-71] [DADA2.R](DADA2.R).
Expected outcome of this function:
                 Forward Complement Reverse RevComp
FWD.ForwardReads    4000          0       0       0
FWD.ReverseReads       0          0       0    4000
REV.ForwardReads       0          0       0    4000
REV.ReverseReads    4000          0       0       0
*4000 is an estimation -> expected a big number
* 0 expect exactly 0

### Get.sample.name
After the primers are removed, the data is ready to be further analysed. First, we look at the sample names; assuming the filenames have this format; change if needed.

### Error rates
We need to look at the error rates of the data, by visualizing it in a plot. The default mode to do so doesn’t work on this data, because it is NovaSeq Illumina data. NovaSeq generates a lot more data, so therefor there will be a lot more errors. There are four great options as an alternative for the error rates. There will not be one ‘perfect’ plot, that’s just the way it is with NovaSeq data. Best way to check is to pick the one where the dots match the black line and keep increasing. Best way to test is on test data, keep in mind that this step may take ±24 hours to test (even with smaller testdata).
-	loessErrfun_mod1 
Error model option 1
-	loessErrfun_mod2 
Error model option 2
-	loessErrfun_mod3 
Error model option 3
-	loessErrfun_mod4 
Error model option 4

### getN
This function in used around the end of the pipeline to make a table of all the unique reads from every output. This table can be useful examining the outputs.



