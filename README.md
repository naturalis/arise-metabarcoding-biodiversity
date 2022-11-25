# COMPOSITION FILES

results
- Testdata3Uitwerkingen.docx > Results from test Novaseq data (50 random samples).  
                              Code plus printout from dada2 its workflow pipeline on test data.
- *.png > results from error models tested on the NovaSeq data.

# PART 1

DADA2 ITS pipeline workflow for Novaseq fungi data
This pipeline runs the dada2 its workflow for NovaSeq paired illumina data

Basic information from the data:
-> even uit R script halen en hierin samenvatten

R packages:
- BiocManager (: BiocManager::install("dada2", version = "3.8"))
- dada2 (BC)
- ShortRead (BC)
- Biostrings (BC)
- magrittr
- dplyr

Other installations:
- cutadapt (anaconda: )

Workflow:
-> stappenplan uitwerken

# PART 2

Implementing the pipeline on HPC infrastructure. In this part we scale
up the analysis on a [MaaS](https://maas.io/how-it-works) machine. This 
machine is loaded with the SSH public keys of project participants, who
can access the machine as:

    ssh -i id_rsa ubuntu@145.136.253.38
 
where `id_rsa` is the path to the associated private key, which is 
normally located at `~/.ssh/id_rsa`. *Note that this SSH connection can 
only be made from behind an EduVPN connection logged in as Naturalis.*

On this machine is installed:

- [R](https://linuxize.com/post/how-to-install-r-on-ubuntu-20-04/) (version 4.2.2)
- [RStudio server](https://www.how2shout.com/linux/install-rstudio-server-open-source-on-ubuntu-20-04-lts/)
- `sudo apt install git libcurl4-openssl-dev python3 python3-pip`
- `python3 -m pip install --user --upgrade cutadapt` # https://cutadapt.readthedocs.io/en/stable/installation.html
- The above R packages, those flagged with (BC) via `Bioc::Manager::install()`, others via `install.packages()`
