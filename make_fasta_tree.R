
# phangorn stuff

# first get the list of MOTUs of interest
motunames = read.csv("RRAtab_filt2_g98.csv") # this is rarified
# can also try using my filtering method
head(motunames)

#motunames=RRA2

# now get the corresponding sequences

seqs = read.delim("nematodes_2018_assembled_setid_filterE2_nl_uniq_c10_ann_assignEMBL_nematoda_sumatra_sumaclust_assignedR134.tab")
head(seqs)

seqsmall = seqs[,c(1,7,9,18,22,27,28)]
head(seqsmall)
seqsmall$MOTU = paste("MOTU_",seqsmall$cluster_sumatra98, sep="")

relevant_seqs = seqsmall[which(seqsmall$MOTU %in% names(motunames)),]
dim(relevant_seqs)
dim(seqsmall)
# these are very similar, suggesting that the MOTUs we kept had a wider array of subsequences

# now select the most abundant sequence
library(tidyverse)
topseq = relevant_seqs %>% group_by(MOTU) %>% top_n(1, count)
head(topseq)
dim(topseq)
# add taxa info
topseq = left_join(topseq, seqs[,18:28])

# now use only these sequences
list.files()
library(phangorn)
library(seqinr)


# nemf = read.FASTA("nematode.fasta")
# subsetty = nemf[c(which(names(nemf) %in% topseq$id))]
# str(subsetty)
# # we now have a list of 181 sequences, rather than >60K!
# hmm = phyDat(nemf,type="DNA")

# make a subset
myfasta= read.fasta(file = "nematode.fasta", seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
#subsetlist<-read.table("list.txt", header=TRUE)

my_fasta_sub = myfasta[names(myfasta) %in% topseq$id]

# since sequences are so different from the MOTU method, then I'll just take the first 1000 names and see if results look any different
my_fasta_sub2 = myfasta[names(myfasta)[1:1000]]

write.fasta(sequences = my_fasta_sub, names = names(my_fasta_sub), nbchar = 80, file.out = "nem_subset.fasta")


#save another subset
write.fasta(sequences = my_fasta_sub2, names = names(my_fasta_sub), nbchar = 80, file.out = "nem_subset2.fasta")


###########################################
# find the animal that is most common for each motu
# bring in RRA2 from analysis -- note that this needs to have the same motu IDs as what was used to construct the phylo
head(RRA2)

RRAagg = RRA2 %>% group_by(Species) %>% summarize_at(vars(MOTU_1:MOTU_557), funs(mean))
# which species has the highest relative frequency for each MOTU?
head(RRAagg)
#RRAagg = RRAagg[-which(RRAagg$Species %in% c("Pos","Neg", "Blank")),]
specs = unlist(apply(RRAagg[,-1], 2, which.max))
specieslink = RRAagg$Species[specs]

plot(tre4, show.tip=FALSE, edge.width=2)
title("Maximum-likelihood tree")
tiplabels(specieslink, cex=.5, fg="transparent", bg="transparent", offset=0.5, frame="n")
axisPhylo()

tre4
length(specieslink) 
# also not much makes sense

# how about if the tip label is the number of species?

head(RRAagg)
RRA2_bin = RRA2[,-c(1:14)] %>% mutate_at(vars(MOTU_1:MOTU_557), funs(ifelse(.>0,1,0)))
RRA2_bin$Species = RRA2$Species
# use only positives because high prev species overall will bias
RRA2_bin_present = RRA2_bin[-which(rowSums(RRA2_bin[,-169])==0),]
RRA2_binagg = RRA2_bin_present %>% group_by(Species) %>% summarize_at(vars(MOTU_1:MOTU_557), funs(mean))
numberhost = RRAagg[,-1] %>% mutate_at(vars(MOTU_1:MOTU_557), funs(ifelse(.>0,1,0))) %>% colSums()
binnumberhost = RRA2_binagg[,-1] %>% mutate_at(vars(MOTU_1:MOTU_557), funs(ifelse(.>0.0,1,0))) %>% colSums()


#############################
library(msa)#multiple sequence alignment package
#By default, msa() runs ClustalW with default parameters

align<-msaClustalOmega('nem_subset.fasta',type="dna", auto=F, cluster=158, dealign=F, order="input", useKimura=T)#this takes a minute
aldna<-as.DNAbin(align)#all sequences have the same length: 464  (481 for me)
dna = fasta2DNAbin(file="nem_subset.fasta")
# wait....I need to check these are actually the same

tips<-rownames(dna)
rownames(align)<-tips
rownames(aldna)<-tips
dnaphydatAll<-phyDat(aldna, type="DNA", levels=NULL)
mt<-modelTest(dnaphydatAll, model=c("JC", "F81", "K80", "HKY", "TrNe", "TrN", "TPM1", "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u", "TIM1e", "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", "SYM", "GTR"))
mt2 = modelTest(dnaphydatAll, model=c("JC69","K80","F81", "K81", "F84","T92" ,"TN93", "GG95", "LOGDET", "BH87" ,"PARALIN", "N" ,"TS" ,"TV", "INDEL", "INDELBLOCK"))

min(mt2$BIC)
which.min(mt2$BIC)
mt2[12,] # K81+G+I is the best out of the options available for dist.dna

#top 10 overall:
mt[order(mt$BIC),][1:10,]
# TPM2u+G+I

env<-attr(mt2, "env")
ls(envir=env)
(fit<-eval(get("K81+G+I", env), env))

dna_distbest<-dist.dna(aldna, model="K80", gamma=F) #note, gamma correction not available for K81
# based on the BIC values, K80 with gamma is better than K81 without

DistMatrixK80<-as.matrix(dna_distbest)
max(DistMatrixK80) # 0.65 is very high; muscle algorithm brings this to 84! 0.46 with omega

dna_distK81<-dist.dna(aldna, model="K81")#note, gamma correction not available for K81
DistMatrixK81<-as.matrix(dna_distK81)
max(DistMatrixK81) # 0.45

DistMatrix_raw<-as.matrix(dist.dna(aldna, model="raw"))
max(DistMatrix_raw)#0.33

hist(DistMatrix_raw, breaks=100)
hist(DistMatrixK80, breaks=50)
hist(DistMatrixK81, breaks=100)

treeNJ<-NJ(dna_distbest)
treeNJ2 = NJ(dna_distK81)
treeUPGMA<-upgma(dna_distbest)
treeUPGMA2=upgma(dna_distK81)
parsimony(treeNJ, dnaphydatAll)#2361
parsimony(treeUPGMA, dnaphydatAll)#2377
parsimony(treeNJ2, dnaphydatAll)
parsimony(treeUPGMA2, dnaphydatAll)
# NJ slightly more parsimonious

treeNJ<-ladderize(treeNJ)
treeNJ2=ladderize(treeNJ2)
par(mfrow=c(1,2))
plot(treeNJ, cex=0.6, show.tip=T)
plot(treeNJ2, cex=0.6, show.tip=T)


# these look very similar to me


############ Plotting with metadata

splink = data.frame(species=specieslink, richness = binnumberhost, MOTU=names(RRAagg[-1]))
ttlab = data.frame((topseq$MOTU[which(topseq$id %in% treeNJ$tip.label)]))
ttlab =topseq[which(topseq$id %in% treeNJ$tip.label),]
names(ttlab)[2]="MOTU"
ttlab$MOTU = paste("MOTU",ttlab$MOTU, sep="_")
ttlab=as.data.frame(ttlab)
# some clusters have multi IDs?
ttlab = left_join(ttlab[,c(1,2,5,6,7,8)], splink, by=c("MOTU"))

par(mfrow=c(1,1))
plot(treeNJ,show.tip=F)
tiplabels(ttlab$scientific_name_ok, cex=0.5, col=transp(num2col(ttlab$bid_ok*100, col.pal=myPal), 0.7), frame="n", offset=0.02)
temp=pretty(min(ttlab$bid_ok*100):max(ttlab$bid_ok*100),5)
legend("topright",fill=transp(num2col(temp, col.pal=myPal),.7),leg=temp, ncol=1, cex=0.7, title="best identity score")
#### This makes zero sense whatsoever
### Family names are completely randomly distributed across the phylogeny? What?
### Have I completely messed up the tip labels? What is going on?

### using another alignment algorithm improved things somewhat

plot(treeNJ, show.tip=T, cex=0.5)
treeNJ$tip.label


#ttlab$richness[which(ttlab$richness==0)] = NA
#treeNJ$tip.label = ttlab$species
plot(treeNJ, cex=0.6, show.tip=F, align.tip.label=F, label.offset=0.005)
tiplabels(ttlab$species, cex=0.5, frame="n", col=transp(num2col(ttlab$richness, col.pal=myPal), 0.7), offset=0.005)
temp <- c(1,2,3,5,10)
legend("topright", fill=transp(num2col(temp, col.pal=myPal),.7), leg=c(1,2,3,5,10), ncol=1, cex=0.7)

# 






##############################
#### Repeat with johans method

# phangorn stuff

# first get the list of MOTUs of interest
motunames = read.csv("RRAtab_filt2_j98.csv") # this is rarified
# can also try using my filtering method
head(motunames)

# now get the corresponding sequences
seqs = read.delim("nematodes_2018_assembled_setid_filterE2_nl_uniq_c10_ann_assignEMBL_nematoda_sumatra_sumaclust_assignedR134.tab")
head(seqs)

seqsmall = seqs[,c(1,7,9,18,22,27,28)]
head(seqsmall)
seqsmall$MOTU = paste("MOTU_",seqsmall$cluster_sumatra98, sep="")

relevant_seqs = seqsmall[which(seqsmall$MOTU %in% names(motunames)),]
dim(relevant_seqs)
dim(seqsmall)
# these are very similar, suggesting that the MOTUs we kept had a wider array of subsequences

# now select the most abundant sequence
library(tidyverse)
topseq = relevant_seqs %>% group_by(MOTU) %>% top_n(1, count)
head(topseq)
dim(topseq)
# add taxa info
topseq = left_join(topseq, seqs[,18:28])

# now use only these sequences
# make a subset
#myfasta= read.fasta(file = "nematode.fasta", seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
my_fasta_sub = myfasta[names(myfasta) %in% topseq$id]

write.fasta(sequences = my_fasta_sub, names = names(my_fasta_sub), nbchar = 80, file.out = "nem_subset_j98.fasta")

###########################################
# find the animal that is most common for each motu
# bring in RRA2 from analysis -- note that this needs to have the same motu IDs as what was used to construct the phylo
RRA2 = read.csv("j98_RRA2_filt_rare.csv")

RRAagg = RRA2 %>% group_by(Species) %>% summarize_at(vars(MOTU_1:MOTU_557), funs(mean))
# which species has the highest relative frequency for each MOTU?
head(RRAagg)
levels(as.factor(RRAagg$Species))
RRAagg = RRAagg[-which(RRAagg$Species %in% c("EC","NTC")),]
specs = unlist(apply(RRAagg[,-1], 2, which.max))
specieslink = RRAagg$Species[specs]

# how about if the tip label is the number of species?
head(RRAagg)
RRA2_bin = RRA2[,-c(1:14)] %>% mutate_at(vars(MOTU_1:MOTU_557), funs(ifelse(.>0,1,0)))
RRA2_bin$Species = RRA2$Species
# use only positives because high prev species overall will bias
RRA2_bin_present = RRA2_bin[-which(rowSums(RRA2_bin[,-159])==0),]
RRA2_binagg = RRA2_bin_present %>% group_by(Species) %>% summarize_at(vars(MOTU_1:MOTU_557), funs(mean))
numberhost = RRAagg[,-1] %>% mutate_at(vars(MOTU_1:MOTU_557), funs(ifelse(.>0,1,0))) %>% colSums()
binnumberhost = RRA2_binagg[,-1] %>% mutate_at(vars(MOTU_1:MOTU_557), funs(ifelse(.>0.0,1,0))) %>% colSums()


#############################
library(msa)#multiple sequence alignment package
#By default, msa() runs ClustalW with default parameters

align<-msaClustalOmega('nem_subset_j98.fasta',type="dna", auto=F, cluster=158, dealign=F, order="input", useKimura=T)#this takes a minute
aldna<-as.DNAbin(align)#all sequences have the same length: 464  (481 for me)
dna = fasta2DNAbin(file="nem_subset_j98.fasta")
# wait....I need to check these are actually the same

tips<-rownames(dna)
rownames(align)<-tips
rownames(aldna)<-tips
dnaphydatAll<-phyDat(aldna, type="DNA", levels=NULL)
mt<-modelTest(dnaphydatAll, model=c("JC", "F81", "K80", "HKY", "TrNe", "TrN", "TPM1", "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u", "TIM1e", "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", "SYM", "GTR"))
mt2 = modelTest(dnaphydatAll, model=c("JC69","K80","F81", "K81", "F84","T92" ,"TN93", "GG95", "LOGDET", "BH87" ,"PARALIN", "N" ,"TS" ,"TV", "INDEL", "INDELBLOCK"))

min(mt2$BIC)
which.min(mt2$BIC)
mt2[12,] # K81+G+I is the best out of the options available for dist.dna

#top 10 overall:
mt[order(mt$BIC),][1:10,]
# TPM2u+G+I

env<-attr(mt2, "env")
ls(envir=env)
(fit<-eval(get("K81+G+I", env), env))

dna_distbest<-dist.dna(aldna, model="K80", gamma=T) #note, gamma correction not available for K81
# based on the BIC values, K80 with gamma is better than K81 without

DistMatrixK80<-as.matrix(dna_distbest)
max(DistMatrixK80) # 0.65 is very high; muscle algorithm brings this to 84! 0.46 with omega

dna_distK81<-dist.dna(aldna, model="K81")#note, gamma correction not available for K81
DistMatrixK81<-as.matrix(dna_distK81)
max(DistMatrixK81) # 0.45

DistMatrix_raw<-as.matrix(dist.dna(aldna, model="raw"))
max(DistMatrix_raw)#0.33

hist(DistMatrix_raw, breaks=100)
hist(DistMatrixK80, breaks=50)
hist(DistMatrixK81, breaks=100)

treeNJ<-NJ(dna_distbest)
treeNJ2 = NJ(dna_distK81)
treeUPGMA<-upgma(dna_distbest)
treeUPGMA2=upgma(dna_distK81)
parsimony(treeNJ, dnaphydatAll)#2361
parsimony(treeUPGMA, dnaphydatAll)#2377
parsimony(treeNJ2, dnaphydatAll)
parsimony(treeUPGMA2, dnaphydatAll)
# NJ slightly more parsimonious

treeNJ<-ladderize(treeNJ)
treeNJ2=ladderize(treeNJ2)
par(mfrow=c(1,2))
plot(treeNJ, cex=0.6, show.tip=T)
plot(treeNJ2, cex=0.6, show.tip=T)


# these look very similar to me


############ Plotting with metadata

splink = data.frame(species=specieslink, richness = binnumberhost, MOTU=names(RRAagg[-1]))
ttlab = data.frame((topseq$MOTU[which(topseq$id %in% treeNJ$tip.label)]))
ttlab =topseq[which(topseq$id %in% treeNJ$tip.label),]
names(ttlab)[2]="MOTU"
ttlab$MOTU = paste("MOTU",ttlab$MOTU, sep="_")
ttlab=as.data.frame(ttlab)
# some clusters have multi IDs?
ttlab = left_join(ttlab[,c(1,2,5,6,7,8)], splink, by=c("MOTU"))

par(mfrow=c(1,1))
plot(treeNJ,show.tip=F)
tiplabels(ttlab$scientific_name_ok, cex=0.5, col=transp(num2col(ttlab$bid_ok*100, col.pal=myPal), 0.7), frame="n", offset=0.002, adj=c(0,0.5))
temp=pretty(min(ttlab$bid_ok*100):max(ttlab$bid_ok*100),5)
legend("topright",fill=transp(num2col(temp, col.pal=myPal),.7),leg=temp, ncol=1, cex=0.7, title="best identity score")


#ttlab$richness[which(ttlab$richness==0)] = NA
#treeNJ$tip.label = ttlab$species
plot(treeNJ, cex=0.6, show.tip=F, align.tip.label=F, label.offset=0.005)
tiplabels(ttlab$species, cex=0.5, frame="n", col=transp(num2col(ttlab$richness, col.pal=myPal), 0.7), offset=0.003, adj=c(0,0.5))
temp <- c(1,2,3,5,10)
legend("topright", fill=transp(num2col(temp, col.pal=myPal),.7), leg=c(1,2,3,5,10), ncol=1, cex=0.7, title="host richness")

# 







#######################################################################





# next step -- compare the species IDs from Kaia

kaia = read.csv("NemRefSeqs_forMakingRefDBfasta.csv")
head(kaia)
dim(kaia)

# find possible matches from kaia
kaia[which(kaia$F_TrimmedSeq %in% seqs$sequence),]$KinsellaTombak.Assignment
length(seqs$sequence)
# only one!
seqs[which(seqs$sequence%in%kaia$F_TrimmedSeq),]$scientific_name_ok
seqs[which(seqs$sequence%in%kaia$F_TrimmedSeq),]$bid_ok

# Genbank said Coronocyclus with 82% similarity

