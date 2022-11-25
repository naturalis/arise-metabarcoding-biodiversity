# Introduction

This section describes how a user with an account with which
they can log in on the Rstudio webserver can further customize
their environment for dada2.

## Installing cutadapt

Cutadapt is a python program that can be installed with pip.
As an RStudio server user this can be done through the "terminal"
tab in the browser interface:

    python3 -m pip install --user --upgrade cutadapt

Note that as a non-root user, the program is installed under 
`~/.local`, so the environment variables PYTHONPATH and PATH
need to be updated so that the executable script is found and
it can load its library dependencies.

## Installing R packages

Within the RStudio server interface, users can install packages
as needed through the console:

    install.packages('BiocManager', 'magrittr', 'dplyr')
    library(BiocManager)
    BiocManager::install('dada2')
    BiocManager::install('ShortRead')
    BiocManager::install('Biostrings')

