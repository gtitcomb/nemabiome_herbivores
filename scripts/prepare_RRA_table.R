
### Prepare RRA data

## 1. Import data and finish filtering

# libraries
library(here)
library(tidyverse)


### 1.1 Setup and visualization

# read in data
OTUtab <- read.csv(here("./data/allMOTUs_sumatra98_GTmethod2020.csv"))
names(OTUtab)%>%tail(3) %>% .[1]

# sequencing depth
median(rowSums(OTUtab[,-c(1,570,571)])) 
mean(rowSums(OTUtab[,-c(1,570,571)]))
hist(rowSums(OTUtab[,-c(1,570,571)]))

# Remove samples with reads = 0
OTUtabsmall <- OTUtab[-which(rowSums(OTUtab[,-c(1,570,571)])==0),]
dim(OTUtabsmall)


### 1.2 Filtering and rarefaction

# Remove samples with reads<1000. 
OTUtabsmall <- OTUtab[-which(rowSums(OTUtab[,-c(1,570,571)])<1000),]

# now Rarefy
set.seed(123)
OTUtabsmall2 = GUniFrac::Rarefy(OTUtabsmall[,-c(1,570,571)], depth=1000)$otu.tab.rff
OTUtabsmall2 = as.data.frame(OTUtabsmall2)
RRAtab = as.data.frame(decostand(OTUtabsmall2, method="total"))
RRAtab$X = OTUtabsmall$X
head(RRAtab)
dim(RRAtab)

# Filter out <1% RRA
RRAtab_filt <- RRAtab %>% mutate_at(vars(MOTU_1:MOTU_568), funs(ifelse(.<0.01,0,.)))

# How many OTUs are now zero?
length(which(colSums(RRAtab_filt[,1:568])==0))

# Remove those columns
RRAtab_filt2 <- RRAtab_filt[,-(which(colSums(RRAtab_filt[,1:568])==0))]
dim(RRAtab_filt2)

# Currently we have 165 MOTUs in 560 samples
RRAtab_filt2 = RRAtab_filt2[,c(167,1:166)]
names(RRAtab_filt2)[1]="Sample_ID"

# save a copy
write.csv(RRAtab_filt2, here("data/RRAtab_filt2_g98.csv"), row.names=F)

# (this is a slightly different version due to seed difference for rarefaction)

# Read in and join datasets
RRA_rare= read.csv(here("data/RRAtab_filt2_g98.csv"))
RRA_hosts=read.csv(here("data/check_summary.csv"))


# Combine mOTU table with host data
RRA2 = full_join(RRA_hosts, RRA_rare)
dim(RRA2)

# Use only the qPCR prevalence
RRA2 = RRA2 %>% 
  filter(rd_2_qpcr_conducted == "Y")
dim(RRA2)

# select columns of interest
head(RRA2)
lastmotu = names(RRA2)%>%tail(1)

# Replace NA values with 0 for samples that were undetected
# This includes samples that were not screened via metabarcoding
RRA2 = RRA2 %>% mutate_at(vars(MOTU_1:all_of(lastmotu)), funs(ifelse(is.na(.),0,.)))

# exclude HIPPO 93 which has the same parasites as a zebra, and none with other hippos
RRA2 = RRA2 %>% 
  filter(Sample_ID != "MRC_17_HIP_93")
