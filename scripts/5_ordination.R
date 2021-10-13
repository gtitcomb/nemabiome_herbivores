
# 5 Ordinations

library(phyloseq)
library(tidyverse)
library(vegan)
library(ggConvexHull)

# this script depends on 1 and 3

# Import files
# tree files
treeNJ1 = read.tree(here("data/treeNJ_K80_gamma_table_1.tree")) 
tipdata1 = read.csv(here("data/nem_taxa_table_1.csv"))
treeNJ2 = read.tree(here("data/treeNJ_K80_gamma_table_2.tree")) 
tipdata2 = read.csv(here("data/nem_taxa_table_2.csv"))

# RRA files
table_1 = read.csv(here("data/RRA_table_1.csv"))
table_2 = read.csv(here("data/RRA_table_2.csv"))

# host files
hosts = read.csv(here("data/host_metadata.csv"))
tree = read.tree(here("data/new_mammal_tree_pruned.newick"))

# decide which dataframe
treeNJ = treeNJ2
data_table = table_2
tipdata = tipdata2

animal_colors=c("darkorchid4","goldenrod1","blueviolet", "deepskyblue3", "hotpink",  "dodgerblue","green3",  "goldenrod","dodgerblue3", "maroon1",  "deepskyblue2", "greenyellow", "dodgerblue4", "lightskyblue","lightskyblue1", "maroon3", "green4", "cyan2")


# create unifrac matrix
treeNJ$tip.label
tipdata = tipdata[match(treeNJ$tip.label, tipdata$seq_id),]
treeNJ$tip.label = tipdata$mOTU
treetree = phy_tree(treeNJ)

# format for phyloseq
treetree = multi2di(treetree)
is.rooted(treetree)

# ensure matches
treetree$tip.label %in% names(data_table)
names(data_table) %in% treetree$tip.label

# create otu table=
otuotu = otu_table(data_table[,-c(dim(data_table)[2]-2, dim(data_table)[2]-1)], taxa_are_rows = F)

# combine into phyloseq object
comboseq = phyloseq(treetree,otuotu)

# run unifrac
set.seed(123)
ufc = phyloseq::UniFrac(comboseq, normalized=T, weighted = T) 
hist(ufc)

# ufc are the unifrac distances -- save this as a file
ufcmat = as.matrix(ufc)
rownames(ufcmat)=data_table$Sample
colnames(ufcmat)=data_table$Sample

#write.csv(ufcmat, here("data/g98_UW_unifrac_dists_2021_check2.csv"),row.names = T)


# Species level
data_table = data_table %>% 
  left_join(hosts, by=c("Sample"="Sample_ID"))

endmotu = names(data_table)[dim(data_table)[2]-27]

otuagg = data_table %>%
  group_by(Species) %>%
  summarize_at(vars(mOTU_1:endmotu), funs(mean(., na.rm=T))) %>% 
  ungroup()

otuotu_agg = otu_table(otuagg[,-1], taxa_are_rows = F)
comboseq_agg = phyloseq(treetree,otuotu_agg)

# create unifrac
set.seed(123)
ufc_agg = UniFrac(comboseq_agg, normalized=T, weighted =T) 
hist(ufc_agg)

ufcmat_agg = as.matrix(ufc_agg)
rownames(ufcmat_agg)=otuagg$Species
colnames(ufcmat_agg)=otuagg$Species

#write.csv(ufcmat_agg, here("data/g98_UW_unifrac_dists_sp_2021_check2.csv"))



### Ordinate ####

### Individual ####

### mOTU table ###
ordMRC = metaMDS(dplyr::select(data_table, mOTU_1:endmotu), k=3, distance="bray", binary=F, autotransform = F)
ordMRC2 = metaMDS(dplyr::select(data_table, mOTU_1:endmotu), k=2, distance="bray",binary=F, autotransform = F)

stressplot(ordMRC)
stressplot(ordMRC2)

gof3 = goodness(ordMRC)
gof2 = goodness(ordMRC2)
plot(gof2 ~ gof3);abline(0,1)

plot(ordMRC);points(ordMRC, display="sites",cex=gof3*100)
plot(ordMRC2);points(ordMRC2, display="sites",cex=gof2*100) # probably two points benefit

pro3 = procrustes(ordMRC, ordMRC2)
plot(pro3) 

# create data tables to use for ggplot
mdsplot = data.frame(x1 = ordMRC$points[,1], x2=ordMRC$points[,2], x3=ordMRC$points[,3], species = data_table$Species, period=data_table$Period, sample=data_table$Sample)
spmds = as.data.frame(scores(ordMRC,"species"))
spmds$mOTU = row.names(spmds)
spmds = left_join(spmds, tipdata)

# calculate average for the same taxa to simplify plot
spmds = spmds %>% group_by(Taxa) %>% summarize_at(vars(NMDS1:NMDS3), funs(mean))
spmds$Taxa2 = sub("_"," ",spmds$Taxa)

MRCord1 = ggplot(mdsplot, aes(x=x1, y=x2))+
  geom_point(aes(col=species), size=2)+
  geom_convexhull(aes(fill=species), alpha=0.2)+
  geom_point(data=spmds, aes(x=NMDS1, y=NMDS2), col="gray", alpha=0.8)+
  ggrepel::geom_text_repel(data=spmds, aes(x=NMDS1, y=NMDS2, label=Taxa2), size=2.5, fontface="italic")+
  scale_fill_manual(values=animal_colors)+
  scale_color_manual(values=animal_colors)+
  theme_bw()+
  coord_fixed(ratio = 1)+
  guides(fill="none", col="none")+
  labs(x="NMDS 1", y="NMDS 2")+
  annotate(geom="text", x=min(mdsplot$x1)+abs(0.15*min(mdsplot$x1)), y=min(mdsplot$x2), label=paste("Stress = ",round(ordMRC$stress,2)))

MRCord2 = ggplot(mdsplot, aes(x=x2, y=x3))+
  geom_point(aes(col=species), size=2)+
  geom_convexhull(aes(fill=species), alpha=0.2)+
  geom_point(data=spmds, aes(x=NMDS1, y=NMDS3), col="gray", alpha=0.8)+
  ggrepel::geom_text_repel(data=spmds, aes(x=NMDS2, y=NMDS3, label=Taxa2), size=2.5, fontface="italic")+
  scale_fill_manual(values=animal_colors)+
  scale_color_manual(values=animal_colors)+
  theme_bw()+
  coord_fixed(ratio = 1)+
  guides(fill="none", col="none")+
  labs(x="NMDS 2", y="NMDS 3")+
  annotate(geom="text", x=min(mdsplot$x1)+abs(0.5*min(mdsplot$x2)), y=min(mdsplot$x3), label=paste("Stress = ",round(ordMRC$stress,2)))

gridExtra::grid.arrange(MRCord1,MRCord2,ncol=1)



### UniFrac Distances ###
unif_ordMRC = metaMDS(ufcmat, k=3, distance="jaccard", binary=T, autotransform = F)
unif_ordMRC2 = metaMDS(ufcmat,  distance="jaccard", binary=T, autotransform = F)

stressplot(unif_ordMRC)
stressplot(unif_ordMRC2)

gof3 = goodness(unif_ordMRC)
gof2 = goodness(unif_ordMRC2)
plot(gof2 ~ gof3);abline(0,1)

plot(unif_ordMRC);points(unif_ordMRC, display="sites",cex=gof3*100)
plot(unif_ordMRC2);points(unif_ordMRC2, display="sites",cex=gof2*100) # probably two points benefit

pro3 = procrustes(unif_ordMRC, unif_ordMRC2)
plot(pro3) 

# create data table
unif_mdsplot = data.frame(x1 = unif_ordMRC$points[,1], x2=unif_ordMRC$points[,2],x3=unif_ordMRC$points[,3], species = data_table$Species, period=data_table$Period, sample=data_table$Sample)

MRCord = ggplot(unif_mdsplot, aes(x=x1, y=x2))+
  geom_point(aes(col=species), size=2)+
  geom_convexhull(aes(fill=species), alpha=0.2)+
  scale_fill_manual(values=animal_colors)+
  scale_color_manual(values=animal_colors)+
  theme_bw()+
  coord_fixed(ratio = 1)+
  guides(fill="none", col="none")+
  annotate(geom="text", x=min(unif_mdsplot$x1)+abs(0.15*min(unif_mdsplot$x1)), y=min(unif_mdsplot$x2), label=paste("Stress = ",round(unif_ordMRC$stress,2)))

MRCord


#### Test Variation Explained by Species ####

# mOTU table

# create env matrix
envMRC = dplyr::select(data_table, Species, Period)
# create distance matrix
MRCdist = vegdist(dplyr::select(data_table, mOTU_1:endmotu), distance="bray")
# implement anova
MRCdiv = adonis2(MRCdist ~ Species,  data=envMRC, permutations=999)
MRCdiv


# UniFrac

# already a distance matrix, so vegdist does not need to be used
# implement anova
MRCdiv_unif = adonis2(ufcmat ~ Species,  data=envMRC, permutations=999)
MRCdiv_unif




#### Species ####

data_table_aggregated = data_table %>% 
  dplyr::select(Species:MSW93_Binomial,BM_KG:UNDERSTORY_SP_MEAN) %>% 
  unique()
data_table_aggregated = data_table_aggregated[match(rownames(ufcmat_agg), data_table_aggregated$Species),]
data_table_aggregated$Species == rownames(ufcmat_agg)

### mOTU ###

motu_agg = data_table %>% 
  group_by(Species) %>% 
  summarize_at(vars(mOTU_1:endmotu), funs(mean))

ordMRCagg = metaMDS(dplyr::select(motu_agg, mOTU_1:endmotu), k=2, distance="bray", binary=F, autotransform = F)
stressplot(ordMRCagg)
plot(ordMRCagg);points(ordMRCagg, display="sites")

# create data tables to use for ggplot
mdsplot_agg = data.frame(x1 = ordMRCagg$points[,1], x2=ordMRCagg$points[,2], Species = motu_agg$Species)
mdsplot_agg = mdsplot_agg %>% left_join(data_table_aggregated)

MRCord = ggplot(mdsplot_agg, aes(x=x1, y=x2))+
  geom_point(aes(col=MSW93_Order, shape=MSW93_Family), size=3, stroke=2)+
  geom_convexhull(aes(fill=GUT),color="gray",alpha=0.3)+
  guides(linetype="none")+
  scale_shape_manual(name="Family",values=c(1,3,5,7,12, 15, 16))+
  scale_fill_manual(name="Fermentation",values=c("aquamarine4", "deepskyblue4","gray"), labels=c("Foregut Fermenter","Hindgut Fermenter"))+
  scale_color_manual(name="Order",values=c("black","gray70","gray90"))+
  ggrepel::geom_text_repel(aes(label=str_sub(Species,1)), nudge_y=0, max.overlaps = 15)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  guides(fill="none", col="none")+
  labs(x="NMDS 1", y="NMDS 2")+
  annotate(geom="text", x=min(mdsplot_agg$x1)+abs(0.15*min(mdsplot_agg$x1)), y=min(mdsplot_agg$x2), label=paste("Stress = ",round(ordMRCagg$stress,2)))

MRCord


### Unifrac ###
unif_ordMRCagg = metaMDS(ufcmat_agg, k=2, distance="bray", binary=F, autotransform = F)


# create data tables to use for ggplot
mdsplot_agg_unif = data.frame(x1 = unif_ordMRCagg$points[,1], x2=unif_ordMRCagg$points[,2], Species = data_table_aggregated$Species)
mdsplot_agg_unif = mdsplot_agg_unif %>% left_join(data_table_aggregated)

MRCord_agg_unif = ggplot(mdsplot_agg_unif, aes(x=x1, y=x2))+
  geom_point(aes(col=MSW93_Order, shape=MSW93_Family), size=3, stroke=2)+
  geom_convexhull(aes(fill=GUT),color="gray",alpha=0.3)+
  guides(linetype="none")+
  scale_shape_manual(name="Family",values=c(1,3,5,7,12, 15, 16))+
  scale_fill_manual(name="Fermentation",values=c("aquamarine4", "deepskyblue4","gray"), labels=c("Foregut Fermenter","Hindgut Fermenter"))+
  scale_color_manual(name="Order",values=c("black","gray70","gray90"))+
  annotate(geom="text", x=min(mdsplot_agg_unif$x1)+abs(0.15*min(mdsplot_agg_unif$x1)), y=min(mdsplot_agg_unif$x2), label=paste("Stress = ",round(unif_ordMRCagg$stress,2)))+
  ggrepel::geom_text_repel(aes(label=str_sub(Species,1)), nudge_y=0, max.overlaps = 15)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  guides(fill="none", col="none")+
  labs(x="NMDS 1", y="NMDS 2")

MRCord_agg_unif


#### Test Variation Explained by Species Traits ####

# mOTU table

# exclude cattle for this particular analysis
MRCdist
# create env matrix
envMRC = data_table_aggregated %>% filter(Species != "Cattle")
motu_agg2 = motu_agg %>% filter(Species != "Cattle")
# reordinate
ordMRCagg = metaMDS(dplyr::select(motu_agg2, mOTU_1:endmotu), k=2, distance="bray", binary=F, autotransform = F)


# create distance matrix
MRCdist = vegdist(dplyr::select(motu_agg2, mOTU_1:endmotu), distance="bray")
# check env and dist match
envMRC$Species == motu_agg2$Species

# run together
par(mfrow=c(2,2))
ordisurf(ordMRCagg ~ BM_KG, envMRC);
ordisurf(ordMRCagg ~ RS_KM2, envMRC);
ordisurf(ordMRCagg ~ GS, envMRC);
ordisurf(ordMRCagg ~ UNDERSTORY_SP_MEAN, envMRC)
par(mfrow=c(1,1))
#

# implement test
MRCdiv = adonis2(MRCdist ~ BM_KG + RS_KM2 + GS + GUT + UNDERSTORY_SP_MEAN,  data=envMRC, permutations=999, by="margin")
MRCdiv


# UniFrac

# already a distance matrix, so vegdist does not need to be used
# check env and dist match
# exclude cattle
exclude = which(row.names(ufcmat_agg)=="Cattle")
ufcmat_agg2 = ufcmat_agg[-exclude,-exclude]
envMRC$Species == row.names(ufcmat_agg2)

unif_ordMRCagg2 = metaMDS(ufcmat_agg2, k=2)
plot(unif_ordMRCagg)

otu_table(comboseq_agg) = otu_table(comboseq_agg)[-exclude,]
uniford = ordinate(comboseq_agg, method="NMDS", distance="unifrac", weighted=T)
plot(uniford)

# run together
par(mfrow=c(2,2))
ordisurf(uniford~ BM_KG, envMRC);
ordisurf(uniford~ RS_KM2, envMRC);
ordisurf(uniford ~ GS, envMRC);
ordisurf(uniford ~ UNDERSTORY_SP_MEAN, envMRC)
par(mfrow=c(1,1))
#

# implement permanova
MRCdiv_unif = adonis2(ufcmat_agg2 ~ BM_KG + RS_KM2 + GS + GUT + UNDERSTORY_SP_MEAN,  data=envMRC, permutations=999, by="margin")
MRCdiv_unif



