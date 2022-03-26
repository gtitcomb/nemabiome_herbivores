

list.files("data/")

tab1 = read.csv("data/RRA_table_2.csv")
head(tab1)
meta = read.csv("data/host_metadata.csv")

library(dplyr)
names(tab1);names(meta)
tab1 = left_join(tab1, meta, by=c("Sample"="Sample_ID"))

library(ggplot2)
tab1 %>% group_by(Species, Period) %>% 
  summarize(n=n()) %>% 
  ggplot(aes(x=Period, y=n))+
  geom_col()+
  facet_wrap(~Species)

tab1 %>% group_by(Species, Period) %>% 
  summarize(n=n()) %>% 
  filter(n>3) %>% 
  group_by(Species) %>% 
  summarize(n=n()) %>% 
  filter(n>1)

# replicated species = cattle, elephant, giraffe, grants, grevys, impala

period_subset = tab1 %>% 
  filter(Species %in%
           c("Cattle", "Elephant", "Giraffe", "Grants gazelle",
             "Grevys zebra", "Impala")) %>% 
  group_by(Species, Period) %>% 
  summarize(n=n()) %>% 
  filter(n>0) %>% 
  select(Species, Period) %>%
  ungroup() %>% 
  left_join(tab1)

library(vegan)

names(period_subset)
ps = select(period_subset, mOTU_1:mOTU_94)
mm = metaMDS(ps)

df = data.frame(x1 = mm$points[,1], x2=mm$points[,2], Species = period_subset$Species,
           Period = period_subset$Period)

rdf = data.frame(Period = levels(as.factor(df$Period)), Rain = c(228, 231, 170, 104,24.8,228,68.48))
df = left_join(df, rdf)

df %>% 
  mutate(Species=recode(Species, `Grants gazelle`="Grant's gazelle", `Grevys zebra`="Grevy's zebra")) %>% 
  ggplot(aes(x=x1, y=x2))+
  ggConvexHull::geom_convexhull(aes(fill=Species, linetype=Period), col="gray50", alpha=0.32, size=0.1)+
  geom_point(aes(col=Species))+
  guides(linetype="none")+
  scale_color_manual(values=c("blueviolet",
                              "green3", "goldenrod", "dodgerblue3", "maroon1", "dodgerblue4"))+
  scale_fill_manual(values=c("blueviolet","green3", "goldenrod", "dodgerblue3", "maroon1", "dodgerblue4"))+
  theme_bw(base_size=14)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  coord_fixed(ratio = 1)+
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","solid","solid"))+
  labs(x="NMDS 1", y="NMDS 2")+
  annotate(geom="text", label=paste("Stress =",round(mm$stress,2)), x=1.9, y=-2)
  

df %>% 
  mutate(Species=recode(Species, `Grants gazelle`="Grant's gazelle", `Grevys zebra`="Grevy's zebra")) %>% 
  ggplot(aes(x=x1, y=x2))+
  ggConvexHull::geom_convexhull(aes(fill=Rain, group=Period, col=120), alpha=0.3)+
  geom_point(aes(col=Rain), alpha=0.5)+
  guides(linetype="none")+
  scale_fill_gradient(low="darkorange", high="darkgreen")+
  scale_color_gradient(low="darkorange", high="darkgreen")+
 theme_bw(base_size=14)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  coord_fixed(ratio = 1)+
  labs(fill="Prior 90-day Rainfall",col="Prior 90-day Rainfall", x="NMDS 1", y="NMDS 2")+
  facet_wrap(~Species)


# now a test
t1d = period_subset %>% 
  select(mOTU_1:mOTU_94)
dm = vegdist(t1d)
#perms = with(tab1, how(nperm = 1000, blocks = Species))
ad = adonis2(dm ~ Species+Period, data=period_subset, permutations = 999, by="margin")
ad


##### Species by Species

sp_to_test = tab1 %>% 
  group_by(Species, Period) %>% 
  summarize(n=n()) %>% 
  filter(n>4) %>% 
  group_by(Species) %>% 
  summarize(n=n()) %>% 
  filter(n>1) %>% 
  select(Species)

find_period_diff = function(species){
  period_use = tab1 %>% 
    filter(Species==species) %>% 
    group_by(Period) %>% 
    summarize(n=n()) %>% 
    select(Period)
  cattle = tab1 %>%
    filter(Species==species) %>%
    filter(Period %in% period_use$Period) %>% 
    select(Period, Sample, mOTU_1:mOTU_94)
  row.names(cattle)=cattle$Sample
  cm = vegdist(cattle[,-c(1,2)])
  cd = adonis2(cm ~ Period, data=cattle, by="margin", permutations = 10000)
  return(cd)
}

andf = data.frame(Species = sp_to_test$Species, nPeriod = rep(0,length(sp_to_test$Species)), p=rep(0,length(sp_to_test$Species)))
for(i in 1:length(sp_to_test$Species)){
  pf = find_period_diff(sp_to_test$Species[i])
  andf[i,]=c(sp_to_test$Species[i], pf$Df[1]+1, pf$`Pr(>F)`[1])
}

p.adjust(andf$p, method="holm")
andf




### Code to find the right rainfall is below:
tab1 %>% 
  group_by(Period) %>% 
  summarize_at(vars(Rain90), funs(mean(., na.rm=T)))

# 2015JUN is 103.378 using the Mpala data
# 2017JUN is 68.48
# 3/1-3/25:
# 7/3-7/23:
rain %>% filter(Date=="1-Mar-17" |Date =="25-Mar-17")
rain %>% filter(Date=="3-Jul-17" |Date =="23-Jul-17")

# typical cumulative rainfall for may is 231

rain %>% filter((JulianDate >1760 & JulianDate <1784) | (JulianDate >17184 & JulianDate<17204)) %>% 
  summarize_at(vars(Cum30), funs(mean))

tab1 %>% filter(Period=="2014MAY")

animal_colors=c("darkorchid4","goldenrod1","blueviolet", "deepskyblue3", "hotpink",  "dodgerblue","green3",  "goldenrod","dodgerblue3", "maroon1",  "deepskyblue2", "greenyellow", "dodgerblue4", "lightskyblue","lightskyblue1", "maroon3", "green4", "cyan2")
levels(as.factor(tab1$Species))
