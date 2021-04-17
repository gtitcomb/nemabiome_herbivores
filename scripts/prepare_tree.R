library(msa)
library(adegenet)
library(phyloseq)
library(QsRutils)
library(phangorn)
library(seqinr)
library(here)

# ClustalOmega appeared to make the most parsimonious trees
align = msaClustalOmega(here("data/nem_subset.fasta"),
                        type="dna",
                        auto=F,
                        cluster=100,
                        dealign=F,
                        order="input",
                        useKimura=T) 

aldna = as.DNAbin(align)
dna = fasta2DNAbin(file=here("data/nem_subset.fasta"))

# Labelling
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
(fit = eval(get("K81+G+I", env), env))

# Use model
dna_distbest = dist.dna(aldna, model="K80", gamma=T) #note, gamma correction not available for K81
# based on the BIC values, K80 with gamma is better than K81 without

DistMatrixK80 = as.matrix(dna_distbest)
max(DistMatrixK80)

dna_distK81 = dist.dna(aldna, model="K81") #note, gamma correction not available for K81
DistMatrixK81 = as.matrix(dna_distK81)
max(DistMatrixK81)

DistMatrix_raw = as.matrix(dist.dna(aldna, model="raw"))
max(DistMatrix_raw)

hist(DistMatrix_raw, breaks=100)
hist(DistMatrixK80, breaks=50)
hist(DistMatrixK81, breaks=100)


treeNJ = NJ(dna_distbest)
treeNJ2 = NJ(dna_distK81)
treeUPGMA = upgma(dna_distbest)
treeUPGMA2 = upgma(dna_distK81)
parsimony(treeNJ, dnaphydatAll)
parsimony(treeUPGMA, dnaphydatAll)
parsimony(treeNJ2, dnaphydatAll)
parsimony(treeUPGMA2, dnaphydatAll)
# NJ slightly more parsimonious
#
treeNJ = ladderize(treeNJ)
treeNJ2 = ladderize(treeNJ2)
treeUPGMA = ladderize(treeUPGMA)
par(mfrow = c(1,2))
plot(treeNJ, cex=0.6, show.tip=T)
plot(treeNJ2, cex=0.6, show.tip=T)

# Similar
# save best tree

treeNJ$edge.length[treeNJ$edge.length<0] = 0

write.tree(treeNJ, file=here("data/g98_treeNJ_K80_gamma_check.tree")) # open this using dendro program and save as new newick file to convert negative branch lengths
#write.tree(treeUPGMA, file=here("data/g98_treeUPGMA_K80_gamma"))