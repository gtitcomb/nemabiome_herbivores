
# 5 Ordinations

#######
library(phyloseq)
library(tidyverse)
library(vegan)
library(ggConvexHull)
library(ape)
library(paco)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
library(here)
######

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
threshold_used = "0.02"

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



### Ordinate ####

### Individual ####

### mOTU table ###

optimizeMDS = function(data_table, k_vector){
  # create empty df
  ord_elbow = data.frame(stress = rep(0,length(k_vector)),
                         dimensions = rep(0,length(k_vector)))
  for(i in 1:length(k_vector)){
    ordi = metaMDS(dplyr::select(data_table, mOTU_1:endmotu), k=k_vector[i], distance="bray", binary=F, autotransform = F, try=5, trymax=5)
    ord_elbow$stress[i]=ordi$stress
    ord_elbow$dimensions[i]=k_vector[i]
  }
  return(ord_elbow)
}

elbow_ord = optimizeMDS(data_table, c(1,2,3,4,5,6,7,8,9,10))

ggplot(elbow_ord, aes(x=dimensions, y=stress))+
  geom_point()+
  geom_line()

# use three 
ordMRC = metaMDS(dplyr::select(data_table, mOTU_1:endmotu), k=3, distance="bray", binary=F, autotransform = F)
gof = goodness(ordMRC2)
mean(gof)
plot(ordMRC);points(ordMRC, display="sites",cex=gof*100)

# create data tables to use for ggplot
mdsplot = data.frame(x1 = ordMRC$points[,1], x2=ordMRC$points[,2], x3=ordMRC$points[,3], species = data_table$Species, period=data_table$Period, sample=data_table$Sample)
spmds = as.data.frame(scores(ordMRC,"species"))
spmds$mOTU = row.names(spmds)
spmds = left_join(spmds, tipdata)

# calculate average for the same genus to simplify plot
spmds = spmds %>% group_by(Genus) %>% summarize_at(vars(NMDS1:NMDS3), funs(mean))
spmds$Taxa2 = sub("_"," ",spmds$Genus)
spmds

MRCord1 = ggplot(mdsplot, aes(x=x1, y=x2))+
  geom_point(aes(col=species), size=2, alpha=0.5)+
  geom_convexhull(aes(fill=species), alpha=0.2)+
  geom_point(data=spmds, aes(x=NMDS1, y=NMDS2), col="gray", alpha=0.8)+
  annotate(geom="text", x=min(mdsplot$x1)+abs(0.15*min(mdsplot$x1)),
           y=min(mdsplot$x2), label=paste("Stress = ",round(ordMRC$stress,2)))+
  ggrepel::geom_text_repel(data=spmds, aes(x=NMDS1, y=NMDS2, label=Taxa2), size=2.5, fontface="italic")+
  scale_fill_manual(values=animal_colors)+
  scale_color_manual(values=animal_colors)+
  theme_bw(base_size=14)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  coord_fixed(ratio = 1)+
  guides(fill="none", col="none")+
  labs(x="NMDS 1", y="NMDS 2")

MRCord1+facet_wrap(~species)

MRCord2 = ggplot(mdsplot, aes(x=x2, y=x3))+
  geom_point(aes(col=species), size=2, alpha=0.5)+
  geom_convexhull(aes(fill=species), alpha=0.2)+
  geom_point(data=spmds, aes(x=NMDS1, y=NMDS3), col="gray", alpha=0.8)+
  annotate(geom="text", x=min(mdsplot$x2)+abs(0.2*min(mdsplot$x2)), y=min(mdsplot$x3), label=paste("Stress = ",round(ordMRC$stress,2)))+
  ggrepel::geom_text_repel(data=spmds, aes(x=NMDS2, y=NMDS3, label=Taxa2), size=2.5, fontface="italic")+
  scale_fill_manual(values=animal_colors)+
  scale_color_manual(values=animal_colors)+
  theme_bw(base_size=14)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  coord_fixed(ratio = 1)+
  guides(fill="none", col="none")+
  labs(x="NMDS 2", y="NMDS 3")
MRCord2+facet_wrap(~species)

mOTU_NMDS = gridExtra::grid.arrange(MRCord1,MRCord2,ncol=2, widths=c(1,1.3))

#ggsave(here(paste("plots/5_mOTU_NMDS",threshold_used,".pdf", sep="")), mOTU_NMDS, width=12, height=5, dpi=300, device="pdf")
ggsave(here(paste("plots/5_mOTU_NMDS",threshold_used,".png", sep="")), mOTU_NMDS, width=12, height=5, dpi=300, device="png")



### UniFrac Distances ###
unif_ordMRC = metaMDS(ufcmat, k=3, distance="jaccard", binary=T, autotransform = F)
#unif_ordMRC2 = metaMDS(ufcmat,  distance="jaccard", binary=T, autotransform = F)

stressplot(unif_ordMRC)
plot(unif_ordMRC);points(unif_ordMRC, display="sites")

 
# create data table
unif_mdsplot = data.frame(x1 = unif_ordMRC$points[,1], x2=unif_ordMRC$points[,2],x3=unif_ordMRC$points[,3], species = data_table$Species, period=data_table$Period, sample=data_table$Sample)

MRCord = ggplot(unif_mdsplot, aes(x=x1, y=x2))+
  geom_point(aes(col=species), size=2)+
  geom_convexhull(aes(fill=species), alpha=0.2)+
  scale_fill_manual(values=animal_colors)+
  scale_color_manual(values=animal_colors)+
  theme_bw()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
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

result_1 = as.data.frame(MRCdiv)


# UniFrac

# already a distance matrix, so vegdist does not need to be used
# implement anova
MRCdiv_unif = adonis2(ufcmat ~ Species,  data=envMRC, permutations=999)
MRCdiv_unif

result_2 = as.data.frame(MRCdiv_unif)



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

elbow_mds_agg = optimizeMDS(motu_agg, c(1:10))
ggplot(elbow_mds_agg, aes(x=dimensions, y=stress))+
  geom_point()+
  geom_line()

ordMRCagg = metaMDS(dplyr::select(motu_agg, mOTU_1:endmotu), k=2, distance="bray", binary=F, autotransform = F)
stressplot(ordMRCagg)
plot(ordMRCagg);points(ordMRCagg, display="sites")

# create data tables to use for ggplot
mdsplot_agg = data.frame(x1 = ordMRCagg$points[,1], x2=ordMRCagg$points[,2], Species = motu_agg$Species)
mdsplot_agg = mdsplot_agg %>% left_join(data_table_aggregated)

MRCord_agg_motu = ggplot(mdsplot_agg, aes(x=x1, y=x2))+
  geom_point(aes(col=MSW93_Order, shape=MSW93_Family), size=3, stroke=2)+
  geom_convexhull(aes(fill=GUT),color="gray",alpha=0.3)+
  guides(linetype="none")+
  scale_shape_manual(name="Family",values=c(1,3,5,7,12, 15, 16))+
  scale_fill_manual(name="Fermentation",values=c("aquamarine4", "deepskyblue4","gray"), labels=c("Foregut Fermenter","Hindgut Fermenter"))+
  scale_color_manual(name="Order",values=c("black","gray70","gray90"))+
  ggrepel::geom_text_repel(aes(label=str_sub(Species,1)), nudge_y=0, max.overlaps = 15)+
  theme_bw(base_size=14)+
  theme(aspect.ratio = 1,
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  guides(fill="none", col="none", shape="none")+
  labs(x="NMDS 1", y="NMDS 2")+
  annotate(geom="text", x=min(mdsplot_agg$x1)+abs(0.25*min(mdsplot_agg$x1)), y=min(mdsplot_agg$x2), label=paste("Stress = ",round(ordMRCagg$stress,2)))

MRCord_agg_motu


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
  annotate(geom="text", x=min(mdsplot_agg_unif$x1)+abs(0.25*min(mdsplot_agg_unif$x1)), y=min(mdsplot_agg_unif$x2), label=paste("Stress = ",round(unif_ordMRCagg$stress,2)))+
  ggrepel::geom_text_repel(aes(label=str_sub(Species,1)), nudge_y=0, max.overlaps = 15)+
  theme_bw(base_size=14)+
  theme(aspect.ratio = 1,
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  guides(fill="none", col="none")+
  labs(x="NMDS 1", y="NMDS 2")

MRCord_agg_unif

sp_ords = gridExtra::grid.arrange(MRCord_agg_motu, MRCord_agg_unif, widths=c(1.75,3), ncol=2)

ggsave(here(paste("plots/5_sp_ords",threshold_used,".pdf",sep="")), sp_ords, width=15, height=5, dpi=300, device="pdf")


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

envMRC = envMRC %>% 
  mutate(GUT = recode(GUT, `FG*`="FG"))

# run together
par(mfrow=c(2,2))
ordisurf(ordMRCagg ~ log(BM_KG), envMRC);
ordisurf(ordMRCagg ~ log(RS_KM2), envMRC);
ordisurf(ordMRCagg ~ GS, envMRC);
ordisurf(ordMRCagg ~ UNDERSTORY_SP_MEAN, envMRC)
par(mfrow=c(1,1))
#

# implement test
MRCdiv = adonis2(MRCdist ~ MSW93_Order +  GUT + log(BM_KG) + GS + UNDERSTORY_SP_MEAN,  data=envMRC, permutations=999, by="margin")
# check range size
MRCdivb = adonis2(MRCdist ~ MSW93_Order +  GUT + log(RS_KM2) + GS + UNDERSTORY_SP_MEAN,  data=envMRC, permutations=999, by="margin")

# compare variance explained by gut and order
MRCgut = adonis2(MRCdist ~ GUT,  data=envMRC, permutations=999, by="margin")
MRCorder = adonis2(MRCdist ~ MSW93_Order,  data=envMRC, permutations=999, by="margin")
MRCboth = adonis2(MRCdist ~ GUT + MSW93_Order,  data=envMRC, permutations=999, by=NULL)


# save results
result_3 = as.data.frame(MRCdiv)
result_3b = data.frame(Predictor = c("Order", "Gut", "Both"), R2 = c(as.data.frame(MRCorder)$R2[1],
                                                                     as.data.frame(MRCgut)$R2[1],
                                                                     as.data.frame(MRCboth)$R2[1]))


# UniFrac

# already a distance matrix, so vegdist does not need to be used
# check env and dist match
# exclude cattle
exclude = which(row.names(ufcmat_agg)=="Cattle")
ufcmat_agg2 = ufcmat_agg[-exclude,-exclude]
envMRC$Species == row.names(ufcmat_agg2)

unif_ordMRCagg2 = metaMDS(ufcmat_agg2, k=2, dist="bray", binary=F, autotransform = F)
plot(unif_ordMRCagg)

# run together
par(mfrow=c(2,2))
ordisurf(unif_ordMRCagg2~ log(BM_KG), envMRC);
ordisurf(unif_ordMRCagg2~ log(RS_KM2), envMRC);
ordisurf(unif_ordMRCagg2 ~ GS, envMRC);
ordisurf(unif_ordMRCagg2~ UNDERSTORY_SP_MEAN, envMRC)
par(mfrow=c(1,1))
#

# implement permanova
MRCdiv_unif = adonis2(ufcmat_agg2 ~ MSW93_Order + GUT + log(BM_KG) + GS + UNDERSTORY_SP_MEAN,  data=envMRC, permutations=999, by="margin")
MRCdiv_unifb = adonis2(ufcmat_agg2 ~ MSW93_Order +  GUT + log(RS_KM2) + GS + UNDERSTORY_SP_MEAN,  data=envMRC, permutations=999, by="margin")

MRCdiv_unif
MRCdiv_unifb # minimal difference

MRCunif_order = adonis2(ufcmat_agg2 ~ MSW93_Order,  data=envMRC, permutations=999, by="margin")
MRCunif_gut = adonis2(ufcmat_agg2 ~ GUT,  data=envMRC, permutations=999, by="margin")
MRCunif_both = adonis2(ufcmat_agg2 ~ GUT+MSW93_Order,  data=envMRC, permutations=999, by=NULL)

result_4 = as.data.frame(MRCdiv_unif)
result_4b = data.frame(Predictor = c("Order", "Gut", "Both"), R2 = c(as.data.frame(MRCunif_order)$R2[1],
                                                        as.data.frame(MRCunif_gut)$R2[1],
                                                        as.data.frame(MRCunif_both)$R2[1]))

## Export the results

species_anova = rbind(result_1,result_2)
write.table(species_anova, here(paste("docs/5_species_anova_ord",threshold_used,".txt",sep="")))

trait_permanova = round(rbind(result_3,result_4),3)
trait_permanovab = rbind(result_3b, result_4b)
write.table(trait_permanova, here(paste("docs/5_trait_permanova_ord",threshold_used,".txt",sep="")))
write.table(trait_permanovab, here(paste("docs/5_trait_permanova_b_ord",threshold_used,".txt",sep="")))


## Cophylogeny

# prune host tree
tree_host = drop.tip(tree, c("Tragelaphus_scriptus","Kobus_ellipsiprymnus"))

# Format correctly
htree = cophenetic(tree_host)
ptree = cophenetic(treeNJ)
gllink = left_join(motu_agg, dplyr::select(data_table_aggregated, Species, MSW93_Binomial))%>%unique()
rownames(gllink)=gllink$MSW93_Binomial
rownames(gllink)[c(5,6,9,16)]=c("Equus_africanus","Tragelaphus_oryx","Nanger_granti","Equus_quagga")
gllink2 = as.data.frame(gllink[,-c(1,dim(gllink)[2])])
rownames(gllink2) = rownames(gllink)

# # Prepare and run analysis
# D = prepare_paco_data(H=htree, P=ptree, HP=gllink2)
# 
# # Commented out due to running time; uncomment to run
# 
# D = add_pcoord(D, correction="cailliez")
# 
# D = PACo(D, nperm=1000, seed=12, method="r0", symmetric=F)
# D2 = PACo(D, nperm=1000, seed=12, method="r0", symmetric=T)
# 
# D = paco_links(D)
# D2 = paco_links(D2)
# 
# res = residuals_paco(D$proc)
# plot(res~D$jackknife)
# 
# sum(res^2)
# print(D$gof)
# print(D2$gof)
# 
# # Additional plotting
# Dres = data.frame(res=res, connection=names(D$jackknife))
# Dres =Dres[order(Dres$res),]
# Dres$sp = rep("",dim(Dres)[1])
# Dres$MOTU = rep("",dim(Dres)[1])
# for(i in 1:length(Dres$res)){
#   Dres$sp[i]=strsplit(Dres$connection, split="-")[[i]][1]
#   Dres$MOTU[i]=strsplit(Dres$connection, split="-")[[i]][2]
# }
# 
# ggplot(Dres, aes(x=reorder(sp, -res), y=res))+
#   geom_point()+
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))

