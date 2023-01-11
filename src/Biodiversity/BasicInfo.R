library(ShortRead)
packageVersion("ShortRead")
library(ggplot2)

# First we make a list of all the data there is.
path <- "/Users/winnythoen/Desktop/BioInformatica/Afstuderen/Testdata2"
list.files(path)

fq <- sort(list.files(path, pattern = "e11", full.names = TRUE))
fq1 <- readFastq(fq[52])
# Now that we have the data in 1 list, we can look up some basic information
head(fq1)
summary(fq1)
reads = sread(fq1)
widths = reads@ranges@width

# Plotten van de gevonden reads in fq1
widths = as.data.frame(reads@ranges@width)

ggplot(widths) +
  geom_histogram(aes(x=reads@ranges@width))

# Kwaliteit bekijken
quals = quality(fq1)
numqscores = as(quals, 'matrix')
avgscore = rowMeans(numqscores, na.rm = T)
avgscore = as.data.frame(avgscore)

ggplot(avgscore) +
  geom_histogram(aes(x=avgscore))
