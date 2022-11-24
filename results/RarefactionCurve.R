library(tidyverse)
library(vegan)
library(ggplot2)
library(readr)

otu <- read_csv("/Users/winnythoen/Desktop/BioInformatica/Afstuderen/Naturalis/arise-metabarcoding-biodiversity/data/OTU97tab_tax.csv")

#rarecurve(t(otu), step=50, cex=0.5)
