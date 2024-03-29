# 6. Create Networks
# Georgia Titcomb

#########
library(igraph)
library(MuMIn)
library(fitdistrplus)
library(tidyverse)
library(ape)
library(caper)
library(picante)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
library(here)
#####

# Import files
# tree files
treeNJ1 = read.tree(here("data/treeNJ_K80_gamma_table_1.tree")) 
tipdata1 = read.csv(here("data/nem_taxa_table_1.csv"))
treeNJ2 = read.tree(here("data/treeNJ_K80_gamma_table_2.tree")) 
tipdata2 = read.csv(here("data/nem_taxa_table_2.csv"))

# RRA files
table_1 = read.csv(here("data/RRA_table_1r.csv"))
table_2 = read.csv(here("data/RRA_table_2r.csv"))

# host files
hosts = read.csv(here("data/host_metadata.csv"))
tree = read.tree(here("data/new_mammal_tree_pruned.newick"))

# decide which dataframe
treeNJ = treeNJ2
data_table = table_2
tipdata = tipdata2
threshold_used = "0.02"

animal_colors=c("darkorchid4","goldenrod1","blueviolet", "deepskyblue3", "hotpink",  "dodgerblue","green3",  "goldenrod","dodgerblue3", "maroon1",  "deepskyblue2", "greenyellow", "dodgerblue4", "lightskyblue","lightskyblue1", "maroon3", "green4", "cyan2")
endmotu = names(data_table)[dim(data_table)[2]-27]
hosts = hosts %>% 
  mutate_at(vars(GUT), funs(ifelse(.=="FG*", "FG",.)))

#############

data_table = data_table %>% 
  left_join(hosts, by=c("Sample"="Sample_ID"))


## Characteristics of host/parasite links ##

# Calculate average RRA for each species and mOTU combination
# (this will also be used later)

avgRRA = data_table %>%
  group_by(Species) %>%
  summarize_at(vars(mOTU_1:endmotu), funs(mean)) %>% 
  ungroup()

# Explore thresholding for average
find_threshold = function(avg_otu_table, vector_to_try){
  store = data.frame(threshold = vector_to_try, n_links = rep(0,length(vector_to_try)))
  for(i in 1:length(vector_to_try)){
    store$threshold[i] = vector_to_try[i]
    store$n_links[i] = avg_otu_table %>% 
      mutate_at(vars(mOTU_1:endmotu), funs(ifelse(.< vector_to_try[i],0,.))) %>% 
      mutate_at(vars(mOTU_1:endmotu), funs(ifelse(.>0,1,0))) %>% 
      dplyr::select(mOTU_1:endmotu) %>% 
      sum()
  }
  return(store)
}

find_threshold(avgRRA, seq(from = 0, to = 0.1, by=0.01)) %>% 
  ggplot(aes(x=threshold, y=n_links))+
  geom_point()+
  geom_line()+
  scale_y_continuous(limits=c(0,300))+
  geom_vline(xintercept=c(0.020), linetype="dotted", col="red", size=1)+
  theme_bw()+
  labs(x="Average RRA Threshold", y="Number of Edges")


# set average RRA thresh.
avgRRA2 = avgRRA %>% 
  mutate_at(vars(mOTU_1:endmotu), funs(ifelse(.< as.numeric(threshold_used),0,.)))

# calculate prevalence of mOTU instances
nRRA = data_table %>%
  mutate_at(vars(mOTU_1:endmotu), funs(ifelse(. > 0,1,0))) %>%
  group_by(Species) %>%
  summarize_at(vars(mOTU_1:endmotu), funs(sum(.)))%>%
  mutate_at(vars(mOTU_1:endmotu), funs(ifelse(.>0,1,0))) %>% 
  ungroup()

# calculate binary version of avgRRA
bin_avgRRA = avgRRA2 %>% 
  mutate_at(vars(mOTU_1:endmotu), funs(ifelse(.>0,1,0))) %>% 
  dplyr::select(mOTU_1:endmotu) %>% 
  colSums()

bin_avgRRA[which.max(bin_avgRRA)]
tipdata %>% 
  filter(mOTU=="mOTU_1")

# Create a plot showing the distribution
bin_avgRRA = bin_avgRRA %>% 
  as.data.frame() %>% 
  dplyr::rename("N_Hosts"=".")
bin_avgRRA$mOTU = row.names(bin_avgRRA)
bin_avgRRA = left_join(bin_avgRRA, tipdata)

# exclude zeros
bin_avgRRA = bin_avgRRA %>% filter(N_Hosts>0)

# give higher taxonomic assignment to NA
bin_avgRRA = bin_avgRRA %>% 
  mutate(new_tax = Family) %>% 
  mutate_at(vars(new_tax), funs(ifelse(is.na(.), "Strongylida (Order)", .)))

# reorder family levels
bin_avgRRA$new_tax = factor(bin_avgRRA$new_tax, levels=c("Chabertiidae", "Cooperiidae", "Haemonchidae",
                                     "Strongylidae", "Trichostrongylidae", "Strongylida (Order)"))

# now reorder MOTUs based on Hosts AND name2
bin_avgRRA = bin_avgRRA[order(bin_avgRRA$new_tax, bin_avgRRA$N_Hosts),] # this should be the correct order
MOTUorder = unique(bin_avgRRA$mOTU)
bin_avgRRA$mOTU = factor(bin_avgRRA$mOTU, levels=MOTUorder)

psite_distribution = bin_avgRRA %>% 
  mutate(new_tax = dplyr::recode(new_tax, Chabertiidae = "Chabertiidae (nodular worms)",
         Cooperiidae = "Cooperiidae (cooperia)",
         Haemonchidae = "Haemonchidae (barber pole worms)",
         Trichostrongylidae = "Trichostrongylidae (trichostrongyles)",
         Strongylidae = "Strongylidae (strongyles)",
         `Strongylida (Order)` = "Strongylida (order)")) %>% 
  ggplot(aes(x=N_Hosts, y=mOTU))+
  geom_point(aes(col=new_tax), size=1.6)+
  geom_line(aes(col=new_tax, group=new_tax), size=1)+
  scale_x_continuous(breaks=c(1:max(bin_avgRRA$N_Hosts)))+
  scale_y_discrete(labels=NULL)+
  #theme_classic(base_size = 14)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(color="gray"),
        text= element_text(size=14),
        axis.text.x = element_text(size=14))+
  guides(col="none")+
  facet_wrap(~new_tax, nrow=6, scales="free_y")+
  labs(x="Number of Host Species", y="mOTU ID")
psite_distribution

#ggsave(here(paste("plots/6_psite_distribution",threshold_used,".png", sep="")),psite_distribution, width=3, height=12, dpi=300, device="png")

### Test aggregation ###
# which theoretical distribution has best fit?
descdist(bin_avgRRA$N_Hosts)
fit.lnorm = fitdist(bin_avgRRA$N_Hosts, "lnorm", discrete = T)
plot(fit.lnorm)
fit.lnorm$aic

fit.pois = fitdist(bin_avgRRA$N_Hosts, "pois")
plot(fit.pois)
fit.pois$aic

fit.nb = fitdist(bin_avgRRA$N_Hosts, "nbinom")
fit.nb$aic
plot(fit.nb)

fit.exp = fitdist(bin_avgRRA$N_Hosts, "exp", discrete=T)
plot(fit.exp)
fit.exp$aic

# Combine
distribution_table = data.frame(Distribution = c("Log-Normal", "Exponential", "Negative Binomial", "Poisson"),
           AIC = c(fit.lnorm$aic, fit.exp$aic, fit.nb$aic, fit.pois$aic))

plot(fit.lnorm); plot(fit.exp); plot(fit.nb); plot(fit.pois)
result_2 = distribution_table

#write_delim(result_2, here(paste("docs/6_distribution_fit",threshold_used,".txt",sep="")))

### Alternative aggregation method ###

head(avgRRA2)
head(hosts)
names(avgRRA2)
t_RRA = hosts %>% 
  group_by(Species, MSW93_Binomial) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  left_join(avgRRA2) %>% 
  filter(is.na(mOTU_1)==F) %>%
  pivot_longer(mOTU_1:mOTU_59, names_to="mOTU", values_to="avgRRA") %>% 
  filter(avgRRA>0) %>% 
  mutate(avgRRA = 1) %>% 
  pivot_wider(names_from=MSW93_Binomial, values_from="avgRRA", -c(n,Species),
              values_fill = 0)
t_RRA = t_RRA[which(specnumber(t_RRA[,-1])>1),]


## core and satellite parasites ##
head(data_table)
tall = data_table %>% 
  mutate_at(vars(mOTU_1:mOTU_94), funs(ifelse(.>0,1,0))) %>% 
  group_by(MSW93_Binomial) %>% 
  summarize_at(vars(mOTU_1:mOTU_94), funs(sum))

tall %>% 
  pivot_longer(cols=mOTU_1:mOTU_94) %>% 
  filter(value>0) %>% 
  ggplot(aes(x=value))+
  geom_histogram(aes(fill=MSW93_Binomial), binwidth = 3)+
  facet_wrap(~MSW93_Binomial, scales="free")+
  guides(fill="none")


### Create Network ###

# need to put all in one column
mat_1_long = avgRRA2 %>% 
  pivot_longer(mOTU_1:endmotu)
head(mat_1_long)
names(mat_1_long)=c("Species", "mOTU", "RRA")

# add species information
mat_1_long = mat_1_long %>% left_join(tipdata, by="mOTU")

# format taxa for plotting
# take the motu id and paste to taxa
mat_1_long$Taxa2 = paste(mat_1_long$Taxa, sub("mOTU_","",mat_1_long$mOTU))
mat_1_long$Taxa2 = sub("_"," ", mat_1_long$Taxa2)

# remove zero
mat_1_long = mat_1_long %>% 
  filter(RRA > 0)

# use taxa id
mat_1_long = mat_1_long %>% dplyr::select(Species.x, Taxa2, RRA)

g = graph.data.frame(mat_1_long, directed=FALSE)

V(g)$type = bipartite_mapping(g)$type 
V(g)$color = ifelse(V(g)$type, "lightblue", "salmon")
V(g)$shape = ifelse(V(g)$type, "square", "circle")
E(g)$color = "lightgray"
V(g)$frame.color =  "white"
V(g)$label.dist = 0.5
E(g)$weight = E(g)$RRA

set.seed(123);plot(g,
     vertex.label.cex = 0.8,
     vertex.size=ifelse(V(g)$type,1,7),
     vertex.label.color = ifelse(V(g)$type,"gray","black"),
     edge.width=E(g)$weight,
     layout=layout_with_fr)

# use at least 8 samples
data_table %>% group_by(Species) %>% summarize(n()) 

names(mat_1_long)[1]="Species"
mat_2_long = mat_1_long %>% 
  filter(Species != "Kudu")

#reorder for plotting
mat_2_long = mat_2_long %>% 
  arrange(RRA)

g = graph.data.frame(mat_2_long, directed=FALSE)

V(g)$type = bipartite_mapping(g)$type 
V(g)$color = ifelse(V(g)$type, "red", "gray40")
V(g)$shape = ifelse(V(g)$type, "square", "circle")
edgecolors = round(E(g)$RRA*10,0)
edgecolors = ifelse(edgecolors == 0, "gray90",
                    ifelse(edgecolors == 1, "yellow",
                    ifelse(edgecolors == 2, "orange",
                    ifelse(edgecolors == 3, "red", 
                    "darkred"))))
E(g)$color = edgecolors
V(g)$frame.color =  "gray40"
V(g)$label.dist = 1
E(g)$weight = E(g)$RRA*10

mat_2_long$Species = as.factor(mat_2_long$Species)
plot_names = V(g)$name %in% c(levels(mat_2_long$Species), mat_2_long$Taxa2[which(mat_2_long$RRA>0.1)])


# Save Network Image

pdf(here(paste("plots/6_bipartite_network",threshold_used,".pdf",sep="")))
set.seed(123); plot(g,
     vertex.label = ifelse(plot_names, V(g)$name, ""),
     vertex.label.cex = 0.8,
     vertex.size=ifelse(V(g)$type,1,7),
     vertex.label.color = ifelse(V(g)$type,"gray60","black"),
     edge.width=E(g)$weight,
     layout=layout_with_fr)
dev.off()




#### Bipartite to Unipartite Projection ###
# project host matrix
g2 = bipartite_projection(g, types=E(g)$type)

# retrieve edgelist
V(g2$proj1)$name
comparisons = combn(V(g2$proj1)$name,2, simplify=T) %>% t()

lala = igraph::as_data_frame(g, "edges")
df = lala %>% 
  dplyr::select(from, to, RRA) %>% 
  pivot_wider(names_from = "to", values_from = "RRA", values_fill=0) # takes a second

df = dplyr::rename(df, name2=from)
df_num = df[,-1]
# this has the RRA for every species and mOTU

head(comparisons)
head(df)

# create empty vector
vec=c()

# calculate the sum product
###################
extent = length(comparisons)/2
for(i in 1:extent){
  v = apply(df_num[df$name2%in%(comparisons[i,]),], MARGIN = 2, FUN=prod) %>% sum()
  vec=c(vec,v)
}
####################
g
checklist = paste(comparisons[,1],comparisons[,2], sep="|")
which((checklist == attr(E(g2$proj1),"vnames")) == F) # should be 0 if in order

# does rearranging work? (all should be true)
checklist[match(attr(E(g2$proj1),"vnames"), checklist)] == attr(E(g2$proj1),"vnames")

# use this rearrangement if so
E(g2$proj1)$weight = vec[match(attr(E(g2$proj1),"vnames"),checklist)]
plot(g2$proj1, edge.width = E(g2$proj1)$weight*100, vertex.label.cex = 0.7)

hostgraph = g2$proj1

## Calcualte centrality
c = igraph::closeness(hostgraph, weights=(1/E(hostgraph)$weight))
b = igraph::betweenness(hostgraph, weights=(1/E(hostgraph)$weight), directed=F)
d = igraph::degree(hostgraph)
e = igraph::eigen_centrality(hostgraph, weights=E(hostgraph)$weight)$vector

# combine into one dataframe
netdata = data.frame(species=names(c),c=c, b=b, d=d, e=e)
netlong = netdata %>% 
  pivot_longer(cols=c(c,b,d,e), names_to="Measure", values_to="Value")

ggplot(netlong, aes(x=reorder(species,-Value), y=Value))+
  geom_col(aes(fill=species), alpha=0.6)+
  theme_bw(base_size=14)+
  theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0))+
  scale_fill_manual(values=animal_colors)+guides(fill="none")+
  facet_wrap(~Measure, scales="free_y")


# phylogenetic distance of hosts
phylodists = cophenetic(tree)

# calculate average distance from all others
totaldists = as.data.frame(rowSums(phylodists)/18) # this is Average pairwise distance
totaldists2 = evol.distinct(tree, type="fair.proportion")
# does this correlate with equal splits?
td3 = evol.distinct(tree, type="equal.splits")
cor.test(totaldists2$w, td3$w)

# save in dataframe
totaldists$sciname = rownames(totaldists)
totaldists2$sciname = totaldists2$Species


cp = cov(phylodists) %>% 
  as.data.frame() %>% 
  mutate(sp1 = row.names(.)) %>% 
  pivot_longer(.,Equus_africanus:Loxodonta_africana) %>% 
  rename("sp2"=name)
pp = phylodists %>% 
  as.data.frame() %>% 
  mutate(sp1 = row.names(.)) %>% 
  pivot_longer(.,Equus_africanus:Loxodonta_africana) %>% 
  rename("sp2"=name, "pdis"=value)
left_join(cp, pp) %>% 
  ggplot(., aes(x=value, y=pdis))+
  geom_point()

evol.distinct(tree);totaldists
cor.test(evol.distinct(tree)$w,totaldists$`rowSums(phylodists)/18`,method="pearson")
# both methods are very similar

# correct sci name
species_link = data.frame(sciname = c("Loxodonta_africana",
                    "Bos_taurus",
                    "Madoqua_kirkii",
                    "Equus_grevyi",
                    "Equus_quagga",
                    "Aepyceros_melampus",
                    "Tragelaphus_oryx",
                    "Alcelaphus_buselaphus",
                    "Giraffa_camelopardalis",
                    "Nanger_granti",
                    "Syncerus_caffer",
                    "Phacochoerus_africanus" ,
                    "Equus_africanus",
                    "Hippopotamus_amphibius",
                    "Oryx_gazella", 
                    "Camelus_dromedarius"), 
                    species = c("Elephant", "Cattle", "DikDik", "Grevys zebra", "Plains zebra", 
                                "Impala", "Eland", "Hartebeest", "Giraffe", "Grants gazelle",
                                "Buffalo", "Warthog", "Donkey", "Hippo", "Oryx", "Camel")
)


# combine and plot data
netdata = left_join(netdata,species_link) %>%
  left_join(.,totaldists) %>% 
  left_join(.,totaldists2)
names(netdata)[7]="Distance"

ggplot(netdata, aes(x=w, y=e))+
  geom_point(size=3,aes(col=species))+
  scale_color_manual(values=animal_colors[-14])+
  geom_smooth(method="lm")+
  guides(col="none")+
  ggrepel::geom_text_repel(aes(label=species))+
  theme_bw()


netdata %>% 
  dplyr::select(species, c:e, w) %>% 
  pivot_longer(c:e, names_to="metric", values_to="value") %>%
  mutate(metric=recode(metric, b="Betweenness", c="Closeness", d="Degree", e="Eigenvector")) %>% 
  ggplot(aes(x=w, y=value))+
  geom_point(size=3,aes(col=species))+
  scale_color_manual(values=animal_colors[-14])+ 
  geom_smooth(method="lm", se=F)+
  guides(col=F)+
  ggrepel::geom_text_repel(aes(label=species))+
  facet_wrap(~metric, scales="free_y")+
  labs(x="Evolutionary Distinctiveness", y="Value")+
  theme_bw()


netlong2 = pivot_longer(netdata, cols=c(c,b,d,e),names_to="Measure", values_to="Value")

netdata_all = left_join(netdata, dplyr::select(data_table,
                                         Species,
                                         MSW93_Order:MSW93_Binomial,BM_KG:UNDERSTORY_SP_MEAN),
                           by=c("species"="Species")) %>% unique()

netlong3 = left_join(netlong2, dplyr::select(netdata_all,
                                          sciname,
                                          MSW93_Order:MSW93_Binomial,BM_KG:UNDERSTORY_SP_MEAN))

measuredict = c(b="Betweenness",c="Closeness", d="Degree", e="Eigenvector")

nodedata2 = ggplot(netlong3, aes(x=w, y=Value))+
  theme_bw(base_size=14)+
  geom_point(aes(col=MSW93_Family, shape=GUT),
             alpha=0.6,
             size=4)+
  scale_color_manual(values=c("dodgerblue","goldenrod","green3","maroon3","goldenrod1","greenyellow","green4"))+ # use animalcols2 for values if color
  facet_wrap(~Measure,
             scales="free_y",
             labeller = labeller(Measure=measuredict))+
  labs(x="Mean Pairwise Distance")

nodedata2

# edit graph for plotting
V(hostgraph)$color = ifelse(V(hostgraph)$name %in% c("Hippo"),
                            "greenyellow",
                            ifelse(V(hostgraph)$name %in% c("Elephant"),
                                   "green3",
                                   ifelse(V(hostgraph)$name %in% c("Warthog"),
                                          "green4",
                                          ifelse(V(hostgraph)$name %in% c("Giraffe"),
                                                 "goldenrod1",
                                                 ifelse(V(hostgraph)$name %in% c("Grevys zebra","Plains zebra", "Donkey"),
                                                        "maroon3", "dodgerblue")))))
E(hostgraph)$color = "lightgray"
V(hostgraph)$frame.color =  "white"
V(hostgraph)$label.dist = 0


# save network image
pdf(here(paste("plots/6_hostgraph",threshold_used,".pdf",sep="")), width=8, height=8)
set.seed(123); plot(hostgraph,
     vertex.label.cex = 1.5,
     vertex.size=10,
     vertex.label.color = "gray30",
     edge.width=E(hostgraph)$weight*100)
dev.off()


# Create models with host information

# exclude cattle because of parasite treatment

netdata2 = netdata_all %>% 
  filter(species != "Cattle")

comp.data = comparative.data(tree, netdata2, names.col="sciname", vcv.dim=3, warn.dropped=TRUE)
comp.data$dropped

ccor = cor.test(comp.data$data$c, comp.data$data$w)
dcor = cor.test(comp.data$data$d, comp.data$data$w)
ecor = cor.test(comp.data$data$e, comp.data$data$w)
bcor = cor.test(comp.data$data$b, comp.data$data$w)

correlations = data.frame(t = c(dcor$statistic,ecor$statistic,bcor$statistic,ccor$statistic),
                          r =c(dcor$estimate,ecor$estimate,bcor$estimate,ccor$estimate),
                          P = c(dcor$p.value,ecor$p.value,bcor$p.value,ccor$p.value),
                          test = c("Degree","Eigenvector","Betweenness","Closeness"))
result_1 = correlations %>% 
  pivot_longer(t:P, names_to="Correlation Details", values_to="value") %>% 
  pivot_wider(names_from=test, values_from=value)
result_1
write_delim(result_1, here(paste("docs/6_centrality_cors",threshold_used, ".txt", sep="")))

### Models ###

# rescale closeness 
comp.data$data$c = comp.data$data$c * 100
# don't use log(BM) and log(RS) at the same time
c_mod = pgls(c~ GUT+ UNDERSTORY_SP_MEAN + log(BM_KG) + GS + w, data = comp.data, lambda=1, kappa=1, delta=1)
dredge(c_mod)[1:5,]
c_mod2 = pgls(c~ GUT, data=comp.data, lambda=1, kappa=1, delta=1)
summary(c_mod2)
cmodval = as.data.frame(summary(c_mod)$coefficients %>% round(.,3))
cmodval$Predictor = row.names(cmodval)
cmodval$Metric = "Closeness"
cmodval2 = as.data.frame(summary(c_mod2)$coefficients %>% round(.,3))
cmodval2$Predictor = row.names(cmodval2)
cmodval2$Metric = "Closeness"

# eigenvector
e_mod = pgls(e ~ GUT + UNDERSTORY_SP_MEAN + log(BM_KG) + GS + w, data=comp.data, lambda=1, kappa=1, delta=1)
dredge(e_mod)[1:5,]
e_mod2 = pgls(e ~ GUT, data=comp.data, lambda=1, kappa=1, delta=1)
emodval = as.data.frame(summary(e_mod)$coefficients %>% round(.,3))
emodval$Predictor = row.names(emodval)
emodval$Metric = "Eigenvector"
emodval2 = as.data.frame(summary(e_mod2)$coefficients %>% round(.,3))
emodval2$Predictor = row.names(emodval2)
emodval2$Metric = "Eigenvector"

# betweenness
b_mod = pgls(b ~ GUT+ UNDERSTORY_SP_MEAN + log(BM_KG) +GS +w, data=comp.data, lambda=1, kappa=1, delta=1)
dredge(b_mod)[1:5,]
b_mod2 = pgls(b ~ 1, data=comp.data, lambda=1, kappa=1, delta=1)
bmodval = as.data.frame(summary(b_mod)$coefficients %>% round(.,3))
bmodval$Predictor = row.names(bmodval)
summary(b_mod)
bmodval$Metric = "Betweenness"
bmodval2 = as.data.frame(summary(b_mod2)$coefficients %>% round(.,3))
bmodval2$Predictor = row.names(bmodval2)
summary(b_mod2)
bmodval2$Metric = "Betweenness"

# degree -- not relevant for 0.002 dataset
d_mod = pgls(d ~ GUT+ UNDERSTORY_SP_MEAN + log(BM_KG) + GS +w, data=comp.data, lambda=1,kappa=1,delta=1)
dredge(d_mod)
# Understory and Gut have some evidence, but not substantially 
# more than the intercept
d_mod2 = pgls(d ~ 1, data=comp.data, lambda=1,kappa=1,delta=1)
dmodval = as.data.frame(summary(d_mod)$coefficients %>% round(.,3))
dmodval$Predictor = row.names(dmodval)
dmodval$Metric = "Degree"
summary(d_mod2)
dmodval2 = as.data.frame(summary(d_mod2)$coefficients %>% round(.,3))
dmodval2$Predictor = row.names(dmodval2)
dmodval2$Metric = "Degree"
summary(d_mod2)

# collect values in table
network_table_vals = as.data.frame(rbind(dmodval, emodval, bmodval, cmodval))
network_table_vals2 = as.data.frame(rbind(dmodval2, emodval2, bmodval2, cmodval2))

# export table
write.csv(network_table_vals, here(paste("docs/6_network_table_vals",threshold_used,".csv",sep="")), row.names = F)
write.csv(network_table_vals2, here(paste("docs/6_network_table_vals2_",threshold_used,".csv",sep="")), row.names = F)


#####################################################################
# combine centrality metrics for visualization
netdata$sc = BBmisc::normalize(netdata$c, method = "range", range=c(0,1))
netdata$sd = BBmisc::normalize(netdata$d, method="range", range=c(0,1))
netdata$se = BBmisc::normalize(netdata$e, method="range", range=c(0,1))
netdata$sb = BBmisc::normalize(netdata$b, method="range", range=c(0,1))
netdata$tot = netdata$sc+netdata$se+netdata$sd+netdata$sb
neworder = netdata$species[order((netdata$sc)+(netdata$se)+(netdata$sd)+(netdata$sb))]
neworder = netdata$species[order(netdata$tot)]

netlongsc = tidyr::pivot_longer(netdata, cols=c(sc,sb,sd,se), names_to="Measure", values_to="Value")

netlongsc$species = factor(netlongsc$species, levels = neworder)

netlongsc$Measure=as.factor(netlongsc$Measure)
levels(netlongsc$Measure)=c("Betweenness","Closeness","Degree","Eigenvector")

# format species labels
change = match(c("Grevys zebra","DikDik","Grants gazelle"),as.character(levels(netlongsc$species)))
levels(netlongsc$species)[change]=c("Grevy's zebra","Dik-dik","Grant's gazelle")

# plot
centscore = ggplot(netlongsc, aes(x=species, y=Measure))+
  geom_tile(aes(fill=Value))+
  scale_fill_viridis_c(option="A")+
  labs(x="")+
  theme_bw(base_size=18)+
  theme(aspect.ratio = 1/5, axis.text.x=element_text(angle=90,hjust=1, vjust=0))

centscore
ggsave(here(paste("plots/6_centrality_scores",threshold_used,".png",sep="")), centscore, dpi=300, device="png", width=7, height=5)


#### Sharing Matrix ####
lar = avgRRA2 %>% 
  pivot_longer(mOTU_1:endmotu) %>% 
  filter(value >0) %>% 
  filter(Species!="Kudu")
g3 = graph_from_data_frame(lar, directed=F)
V(g3)$type = bipartite_mapping(g3)$type 

g4 = bipartite_projection(g3, types=E(g3)$type)

proj_host_data = data.frame(get.edgelist(g4$proj1), E(g4$proj1)$weight)
names(proj_host_data)=c("Species1","Species2","N_Shared")
proj_host_data2 = proj_host_data
proj_host_data2$Species1 = proj_host_data$Species2
proj_host_data2$Species2 = proj_host_data$Species1
proj_host_data2$RRA_sum_prod = proj_host_data$RRA_sum_prod

proj_host_data3 = rbind(proj_host_data, proj_host_data2)

addto = get.adjacency(g4$proj1) %>% as.matrix() %>% rowSums() %>% as.data.frame()
names(addto)="N_H_Share"
addto$Species = rownames(addto)
addto = addto[order(addto$N_H_Share),]

proj_host_data3$Species1 = factor(proj_host_data3$Species1, levels= rev(c("Giraffe", "Eland","Camel", "Oryx", "Hartebeest","Buffalo", "Impala", "DikDik","Grants gazelle", "Donkey", "Cattle", "Hippo","Warthog","Plains zebra","Grevys zebra", "Elephant")))
proj_host_data3$Species2 = factor(proj_host_data3$Species2, levels=c("Giraffe", "Eland","Camel", "Oryx", "Hartebeest","Buffalo", "Impala", "DikDik","Grants gazelle", "Donkey", "Cattle", "Hippo", "Warthog","Plains zebra","Grevys zebra", "Elephant"))


# add self
g_ad = get.adjacency(g3, sparse=F) %>% as.data.frame()
self = colSums(g_ad)[1:16]
add_self = data.frame(Species1=names(self), Species2=names(self), N_Shared=self)
proj_host_data3=rbind(proj_host_data3, add_self)

levels(proj_host_data3$Species1)[c(9,8,2)]=c("Dik-dik","Grant's gazelle","Grevy's zebra")
levels(proj_host_data3$Species2)[c(8,9,15)]=c("Dik-dik","Grant's gazelle","Grevy's zebra")

share_matrix = ggplot(proj_host_data3, aes(x=Species1, y=Species2))+
  geom_tile(aes(fill=N_Shared))+
  scale_fill_gradient(low="gray90", high="red")+
  geom_text(aes(label=ifelse(N_Shared >0, round(N_Shared, 2),"")))+
  theme_bw(base_size=14)+
  theme(axis.text.x = element_text(angle=90, hjust=0, vjust=0))+
  theme(axis.text.y = element_text(hjust=0, vjust=1))


share_matrix

ggsave(here(paste("plots/6_share_matrix",threshold_used,".png", sep="")),
       share_matrix, dpi=300, device="png", width=10, height=7)



## Generate a summary 

# Remake long matrix
mat_long = avgRRA2 %>% 
  pivot_longer(mOTU_1:endmotu)
head(mat_long)
names(mat_long)=c("Host_Species", "mOTU", "RRA")
# add taxa information
mat_long = mat_long %>% left_join(tipdata, by="mOTU")
mat_long$Taxa2 = paste(mat_long$Taxa, sub("mOTU_","",mat_long$mOTU))
mat_long$Taxa2 = sub("_"," ", mat_long$Taxa2)
head(mat_long)

parasite_summary = mat_long %>% 
  group_by(Host_Species, Taxa2) %>% 
  summarize_at(vars(RRA, Boot), funs(mean)) %>% 
  pivot_wider(names_from="Host_Species", values_from="RRA") %>% 
  mutate_at(vars(Buffalo:Warthog), funs(round(.,2)))

write.csv(parasite_summary, here(paste("docs/parasite_summary",threshold_used,".csv")), row.names=F)


# Feeding type and family
gut_fam = mat_long %>% 
  left_join(hosts, by=c("Host_Species"="Species")) %>% 
  group_by(Sample_ID, GUT, Family) %>% 
  summarize_at(vars(RRA), funs(sum)) %>% 
  group_by(GUT, Family) %>% 
  summarize_at(vars(RRA), funs(mean)) %>% 
  mutate_at(vars(Family), funs(ifelse(is.na(Family), "No Match", .)))

gut_fam$Family = as.factor(gut_fam$Family)
gut_fam$Family = factor(gut_fam$Family, levels=rev(c("Strongylidae","Cooperiidae","Trichostrongylidae","Haemonchidae","Chabertiidae","No Match")))

gut_fam_plot = ggplot(gut_fam,aes(x=GUT, y=Family))+
  geom_tile(aes(fill=RRA))+
  scale_fill_gradient(low="white", high="darkred")+
  theme_bw(base_size=14)+
  labs(x="Digestive Strategy", y="Parasite Family", fill="Mean RRA")
gut_fam_plot

gut_plot_data = mat_long %>% 
  left_join(dplyr::select(data_table,Sample, Species,GUT), by=c("Host_Species"="Species")) %>% 
  group_by(Sample, GUT, Family) %>% 
  summarize_at(vars(RRA), funs(sum)) %>% 
  group_by(GUT, Family) %>% 
  mutate_at(vars(Family), funs(ifelse(is.na(.), "No match", .))) %>% 
  mutate(Family = dplyr::recode(Family, Chabertiidae = "Chabertiidae (nodular worms)",
                          Cooperiidae = "Cooperiidae (cooperia)",
                          Haemonchidae = "Haemonchidae (barber pole worms)",
                          Trichostrongylidae = "Trichostrongylidae (trichostrongyles)",
                          Strongylidae = "Strongylidae (strongyles)"))
gut_plot_data$Family=as.factor(gut_plot_data$Family)
gut_plot_data$Family = factor(gut_plot_data$Family,
                              levels=levels(gut_plot_data$Family)[c(1,2,3,5,6,4)])

gut_fam = ggplot(gut_plot_data, aes(x=GUT, y=RRA))+
  geom_boxplot(aes(fill=Family), position="dodge")+
  #geom_errorbar(aes(col=Family, ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.5))
  theme_bw(base_size=14)+
  labs(x="Digestive Strategy", y="RRA")
gut_fam

ggsave(here(paste("plots/6_gut_para_families",threshold_used,".png")), gut_fam, device="png", dpi=300)
          