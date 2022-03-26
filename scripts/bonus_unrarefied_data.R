# 2. Unrarefied data comparison
# Georgia Titcomb

# Use Ctrl+Shift+F10 to restart R, clear objects from workspace

# this script depends on script 1
#####################
library(vegan)
library(binom)
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(ape)
library(MCMCglmm)
library(emmeans)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
library(here) # ensure wd is in parent folder before running
#####################

# comparison between rarefied (table_2 from old analysis) and unrarefied (table_2b is new)
# A. read in RRA tables and combine with metadata

table_1 = read.csv(here("data/RRA_table_1b.csv"))
table_1og = read.csv(here("data/RRA_table_1.csv"))

table_2 = read.csv(here("data/RRA_table_2b.csv"))
table_2og = read.csv(here("data/RRA_table_2.csv"))
table_2og==table_2
table_1og$mOTU_1[1:10]
table_1$mOTU_1[1:10]

max(abs(table_1[,1:96]-table_1og[,1:96]))
max(abs(table_2[,1:80]-table_2og[,1:80]))

mean(unlist(abs(table_1[,1:96]-table_1og[,1:96])))
mean(unlist(abs(table_2[,1:80]-table_2og[,1:80])))

cor.test( unlist(table_1[,1:96]), unlist(table_1og[,1:96]) , method="pearson")
cor.test( unlist(table_2[,1:80]), unlist(table_2og[,1:80]) , method="pearson")



head(table_1)[,1:10]

ut2 = read.csv(here("data/unrarefied2pct.csv"))
ut1 = read.csv(here("data/unrarefied02pct.csv"))
ut3 = read.csv(here("data/unrarefied_unfiltered.csv"))
ut4 = read.csv(here("data/rarefied_unfiltered.csv"))
names(ut4)

names(ut2)[82]="sample"
names(ut1)[98]="sample"
names(ut3)[1]="sample"
names(ut4)[215]="sample"

# transpose - ut2
tut2 = as.data.frame(t(ut2[,-c(81:82)]))
names(tut2) = ut2$sample
head(tut2)
tname = as.data.frame(t(ut2[,81:82]))
names(tname)=tname[2,]
tname2 = as.data.frame(t(tname))

# ut1
tut1 = as.data.frame(t(ut1[,-c(97:98)]))
names(tut1) = ut1$sample
head(tut1)
tname = as.data.frame(t(ut1[,97:98]))
names(tname)=tname[2,]
tname1 = as.data.frame(t(tname))

# ut3
tut3 = as.data.frame(t(ut3[,2:215]))
names(tut3)=ut3$sample
tname = as.data.frame(t(ut3[,c(1,216)]))
names(tname)=tname[1,]
tname3 = as.data.frame(t(tname))

# ut4
ut4 = ut4 %>% 
  mutate_at(vars(mOTU_1:mOTU_568), funs(round(.,0)))
which(rowSums(ut4[,-215])==0)

tut4 = as.data.frame(t(ut4[,1:214]))
names(tut4)=ut4$sample
tname = as.data.frame(t(ut4[,c(215)]))
names(tname)=tname[1,]
tname4 = as.data.frame(t(tname))
anyNA(tut4)
tut4=tut4[-which(rowSums(tut4)<200),]

library(phyloseq)
ex1 = phyloseq(sample_data(tname1), otu_table(tut1[-c(97:98),], taxa_are_rows=T))
ex2 = phyloseq(sample_data(tname2), otu_table(tut2[-c(81,82),], taxa_are_rows=T))
ex3 = phyloseq(sample_data(tname3), otu_table(tut3[-c(1,216),], taxa_are_rows=T))
ex4 = phyloseq(sample_data(tname4), otu_table(tut4, taxa_are_rows=T))

library(mirlyn)
rrex1 = rarefy_whole_rep(ex1,rep=10, steps=seq(from=0.0001, to=0.02, by=0.00005))
rrex2 = rarefy_whole_rep(ex2,rep=10, steps=seq(from=0.0001, to=0.02, by=0.00005))
rrex3 = rarefy_whole_rep(ex3,rep=10, steps=seq(from=0.0001, to=0.02, by=0.00005))
rrex4 = rarefy_whole_rep(ex4,rep=10, steps=seq(from=0.0001, to=0.02, by=0.00005))

rcurv = rarecurve(rrex1)

p1 = ggplot(rrex1, aes(x=LibSize, y=ObsASVCount))+
  geom_line(aes(col=Sample), alpha=0.2)+
  guides(col="none")+
  theme_classic(base_size=14)+
  geom_vline(xintercept=1000)+
  labs(x="Rarefied Library Size", y="Observed mOTUs", title="RRA <0.002% excluded")

p2 = ggplot(rrex2, aes(x=LibSize, y=ObsASVCount))+
  geom_line(aes(col=Sample), alpha=0.2)+
  guides(col="none")+
  theme_classic(base_size=14)+
  geom_vline(xintercept=1000)+
  labs(x="Rarefied Library Size", y="Observed mOTUs", title="RRA <0.02% excluded")

p3 = ggplot(rrex3, aes(x=LibSize, y=ObsASVCount))+
  geom_line(aes(col=Sample), alpha=0.2)+
  guides(col="none")+
  theme_classic(base_size=14)+
  geom_vline(xintercept=1000)+
  labs(x="Rarefied Library Size", y="Observed mOTUs", title="No RRA filtering")


gridExtra::grid.arrange(p2, p1, p3, ncol=1)


# rarefaction isn't really necessary after performing the 0.2 and 0.02 exclusion
# it is exactly the same


