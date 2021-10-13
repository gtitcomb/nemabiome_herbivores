
# 3. Construct Nematode Tree
# 6 October 2021
# Georgia Titcomb

# Clear environment
# CTRL+SHIFT+F10

# this script depends on script 1

library(tidyverse)
library(dada2)
library(msa)
library(adegenet)
library(phyloseq)
library(QsRutils)
library(phangorn)
library(seqinr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
library(here)

# read in data
table_1 = read.csv(here("data/RRA_table_1.csv"))
table_2 = read.csv(here("data/RRA_table_2.csv"))
motuinfo = read.csv(here("data/raw_motu_table_mpala_nematodes.csv"))

motu_names = names(table_1)[1:(dim(table_1)[2]-2)]

# select most abundant sequence per cluster
seqsmall = motuinfo %>% 
  arrange(desc(count)) %>% 
  group_by(cluster_sumatra98) %>% 
  dplyr::slice(1)

seqsmall = dplyr::select(seqsmall, id, cluster_sumatra98, sequence) 
seqsmall$MOTU = paste("mOTU_",seqsmall$cluster_sumatra98, sep="")

# all motus
# use only motus in table
seqsmall = seqsmall %>% 
 filter(MOTU %in% motu_names)

# assign taxonomy
# read in latest database
assigned = assignTaxonomy(seqsmall$sequence,
                          here("data/Nematode_ITS2_1.0.0_dada2.fasta"),
                          minBoot = 80,
                          outputBootstraps = T,
                          verbose = T)
nems_complete = as.data.frame(matrix(unlist(assigned$tax), nrow = length(seqsmall$sequence)))
nems_boot = as.data.frame(matrix(unlist(assigned$boot), nrow=length(seqsmall$sequence)))

names(nems_complete) = c("Kingdom", "Other", "Phylum", "Class", "Order", "Family", "Genus", "Species")
names(nems_boot) = paste("boot.",names(nems_complete), sep="")
nems_complete = cbind(nems_complete, nems_boot)
nems_complete$seq_id = seqsmall$id
nems_complete$mOTU = seqsmall$MOTU
head(nems_complete)

# any boot < 100 for order strongylida?
nems_complete = nems_complete %>% 
  filter(mOTU!="mOTU_NA")

# now use only these sequences
# make a subset
myfasta = read.fasta(file = here("data/nematode.fasta"), seqtype = "AA",as.string = TRUE, set.attributes = FALSE)

my_fasta_sub = myfasta[names(myfasta) %in% nems_complete$seq_id]

# save subset
write.fasta(sequences = my_fasta_sub, names = names(my_fasta_sub), nbchar = 80, file.out = here("data/nem_subset.fasta"))



#### CONSTRUCT TREE

# ClustalOmega appeared to make the most parsimonious trees
align = msaClustalOmega(here("data/nem_subset.fasta"),
                        type="dna",
                        auto=F,
                        cluster=100,
                        dealign=F,
                        order="input",
                        useKimura=T) 

aldna = as.DNAbin(align)
dna = fasta2DNAbin(file=here("data/nem_subset.fasta")) # not sure why this takes so long

# Labeling
tips = rownames(dna)
rownames(align) = tips
rownames(aldna) = tips
dnaphydatAll = phyDat(aldna, type="DNA", levels=NULL)

# compare models
mt2 = modelTest(dnaphydatAll,
                model=c("JC69","K80","F81", "K81", "F84","T92" ,"TN93", "GG95", "LOGDET", "BH87" ,"PARALIN", "N" ,"TS" ,"TV", "INDEL", "INDELBLOCK"))

min(mt2$BIC)
which.min(mt2$BIC)
mt2[12,] # K81+G+I is the best out of the options available for dist.dna

#top 10 overall:
mt2[order(mt2$BIC),][1:10,]

env = attr(mt2, "env")
ls(envir=env)
(fit = eval(get("K80+G+I", env), env))

# Use model
dna_distbest = ape::dist.dna(aldna, model="K80", gamma=T) #note, gamma correction not available for K81
# based on the BIC values, K80 with gamma is better than K81 without

DistMatrixK80 = as.matrix(dna_distbest)
max(DistMatrixK80)

dna_distK81 = dist.dna(aldna, model="K81") #note, gamma correction not available for K81
DistMatrixK81 = as.matrix(dna_distK81)
max(DistMatrixK81)

dist_raw = dist.dna(aldna, model="raw")
DistMatrix_raw = as.matrix(dist_raw)
max(DistMatrix_raw)

hist(DistMatrix_raw, breaks=100)
hist(DistMatrixK80, breaks=50)
hist(DistMatrixK81, breaks=100)


treeNJ = NJ(dna_distbest) # K80, gamma T
treeNJ2 = NJ(dna_distK81) # K81, no gamma
treeNJraw = NJ(dist_raw) # K81, no gamma


parsimony(treeNJ, dnaphydatAll)
parsimony(treeNJ2, dnaphydatAll)
parsimony(treeNJraw, dnaphydatAll)

treeNJ = ladderize(treeNJ)
plot(treeNJ)

# rename tips
# first make a final ID column based on the boot level
nems_complete$Best_ID
class_id = nems_complete %>% 
  filter(is.na(Order)) %>% 
  dplyr::select(Class, boot.Class, seq_id, mOTU) %>% 
  mutate(Level = "Class") %>% 
  rename(Class = "Taxa", boot.Class="Boot")
order_id = nems_complete %>% 
  filter(is.na(Family)==T, is.na(Order)==F) %>% 
  dplyr::select(Order, boot.Order, seq_id, mOTU) %>% 
  mutate(Level = "Order") %>% 
  rename(Order = "Taxa", boot.Order="Boot")
family_id = nems_complete %>% 
  filter(is.na(Genus)==T, is.na(Family)==F) %>% 
  dplyr::select(Family, boot.Family, seq_id, mOTU) %>% 
  mutate(Level = "Family") %>% 
  rename(Family = "Taxa", boot.Family="Boot")
genus_id = nems_complete %>% 
  filter(is.na(Species)==T, is.na(Genus)==F) %>% 
  dplyr::select(Genus, boot.Genus, seq_id, mOTU) %>% 
  mutate(Level="Genus") %>% 
  rename(Genus = "Taxa", boot.Genus="Boot")
species_id = nems_complete %>% 
  filter(is.na(Species)==F) %>% 
  dplyr::select(Species, boot.Species, seq_id, mOTU) %>% 
  mutate(Level="Species") %>% 
  rename(Species = "Taxa", boot.Species="Boot")

nems_complete2 = rbind(class_id, order_id, family_id, genus_id, species_id)
# retain all additional information
nems_complete3 = left_join(nems_complete2, nems_complete)

#treeNJ$tip.label = nems_complete2[match(treeNJ$tip.label, nems_complete2$seq_id),]$Taxa
#plot(treeNJ, cex=0.6, show.tip=T)

treeNJ$edge.length[treeNJ$edge.length<0] = 0
plot(treeNJ, cex=0.6, show.tip=T)

write.tree(treeNJ, file=here("data/treeNJ_K80_gamma_table_1.tree")) 
write.csv(nems_complete3, file=here("data/nem_taxa_table_1.csv"), row.names=F)

# prune tree to table 2 taxa
table_2_nems_exclude = nems_complete3[-which(nems_complete3$mOTU %in% names(table_2)[-c(dim(table_2)[2]-1,dim(table_2)[2])]),]
table_2_nems = nems_complete3[-which(nems_complete3$mOTU %in% table_2_nems_exclude$mOTU),]

new_nj = drop.tip(treeNJ,table_2_nems_exclude$seq_id)
plot(new_nj, cex=0.6, show.tip=T)

write.tree(new_nj, file=here("data/treeNJ_K80_gamma_table_2.tree")) 
write.csv(table_2_nems, file=here("data/nem_taxa_table_2.csv"), row.names=F)
