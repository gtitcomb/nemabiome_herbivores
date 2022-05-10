
# Phylogenetic Diversity

# this script depends on scripts 1 and 3
# Use CTRL+SHIFT+F10 to clear and clear workspace objects

########

library(picante)
library(car)
library(DHARMa)
library(MCMCglmm)
library(glmmTMB)
library(tidyverse)
library(adegenet)
library(emmeans)
library(ggtree)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
library(here)

########

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
treeNJ = treeNJ1
data_table = table_1
tipdata = tipdata1
threshold_used = "0.002"

# order the tipdata correctly
tipdata = tipdata[match(treeNJ$tip.label, tipdata$seq_id),]

myPal = colorRampPalette(c("red","yellow","green","blue"))

png(here(paste("plots/4_tree",threshold_used,".png",sep="")), res=300, width=8, height=8, units="in")
par(mfrow=c(1,1))
plot(treeNJ, show.tip=F, x.lim = 0.85)
plot(treeNJ, show.tip=F)
tiplabels(tipdata$Taxa, cex=0.5, col=transp(num2col(tipdata$Boot, col.pal=myPal), 0.7), frame="n", offset=0, adj=0)
temp = pretty(min(tipdata$Boot):max(tipdata$Boot),5)
legend("topright",fill=transp(num2col(temp, col.pal=myPal),.7), leg=temp, ncol=1, cex=0.7, title="Boot")
dev.off()


# ladderize
phylo_nems=ladderize(treeNJ)

# rename tips to the mOTU
phylo_nems$tip.label = tipdata$mOTU

# find end label
end = dim(data_table)[2]-2

# calculate PD
set.seed(123)
PDData_sum = data.frame(data_table[,c(end+1,end+2)], ses.pd(data_table[,1:end], tree=phylo_nems, include.root = F, runs=100, null.model="richness"))
# singletons are given NA


# add host data
PDData_sum = left_join(PDData_sum, hosts, by=c("Sample"="Sample_ID"))

PDData_sum %>% 
  filter(is.na(pd.obs)) %>% 
  group_by(Species) %>% 
  summarize(n())

# exclude NA
PDData_sum = PDData_sum %>% 
  filter(is.na(pd.obs)==F)

# Exploratory Plots:
# relationship between SR and PD
ggplot(PDData_sum, aes(x=ntaxa, y=pd.obs))+
  geom_point(aes(col=Species))+
  geom_smooth(aes(col=Species),method="lm", formula=y~log(x), se=F)+
  facet_wrap(~Species)

# reassign gut type
PDData_sum = PDData_sum %>% 
  mutate_at(vars(GUT), funs(ifelse(.=="FG*", "FG", .)))

animal_colors=c("darkorchid4","goldenrod1","blueviolet", "deepskyblue3", "hotpink",  "dodgerblue","green3",  "goldenrod","dodgerblue3", "maroon1",  "deepskyblue2", "greenyellow", "dodgerblue4", "lightskyblue","lightskyblue1", "maroon3", "green4", "cyan2")

pdz = ggplot(PDData_sum, aes(x=GUT, y=pd.obs.z))+
  geom_hline(yintercept=-2, linetype="dotted")+
  geom_hline(yintercept = 2, linetype="dotted")+
  geom_boxplot(fill=NA, col="gray50", size=1)+
  geom_boxplot(data=PDData_sum, aes(x=GUT, y=pd.obs.z, col=Species), size=0.8, position=position_dodge(width=0.5), width=0.25, alpha=0.6)+
  theme_bw(base_size=16)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  labs(x="Gut Type", y="PD Z-Value")+
  scale_color_manual(values=c(animal_colors))+
  guides(col="none")
  pdz

pd = ggplot(PDData_sum, aes(x=GUT, y=pd.obs))+
  geom_boxplot(fill=NA, col="gray50", size=1)+
  geom_boxplot(data=PDData_sum, aes(x=GUT, y=pd.obs, col=Species), size=0.8, position=position_dodge(width=0.5), width=0.25, alpha=0.6)+
  theme_bw(base_size=16)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  labs(x="Diet", y="Faith's PD")+
  scale_color_manual(values=c(animal_colors[-13]))+
  #facet_wrap(~Fermentation, scales="free_x")+
  guides(color="none")
pd

pdplot = gridExtra::grid.arrange(pd, pdz, ncol=2)

ggsave(here(paste("plots/4_pd_plots_sp",threshold_used,".png", sep="")), pdplot, device="png", dpi=300, width=10, height=5, units="in")


# species, ungrouped

pdplot_sp = PDData_sum %>% 
  mutate(Species = recode(Species, `Grevys zebra` = "Grevy's zebra", `Grants gazelle`="Grant's gazelle", `DikDik`="Dik-dik")) %>% 
  ggplot(aes(x=reorder(Species,-pd.obs,na.rm=T), y=pd.obs))+
  geom_jitter(width=0.01,height=0.01, alpha=0.5, aes(color=Species))+
  geom_boxplot(aes(fill=Species), color="gray50",outlier.shape=NA)+
  guides(fill="none",col="none")+
  theme_bw(base_size=16)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  scale_color_manual(values=(animal_colors))+
  scale_fill_manual(values=(animal_colors))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0))+
  labs(x="Host Species",y="PD")
pdplot_sp

# total PD by species
PD_spmod = glm(pd.obs~ Species, data=PDData_sum, family="gaussian")
summary(PD_spmod)
Anova(PD_spmod)

# results
y_hat = predict(PD_spmod, type="response")
y_bar = mean(PDData_sum$pd.obs)
PDtemp = 1-sum((PDData_sum$pd.obs-y_hat)^2)/sum((PDData_sum$pd.obs-y_bar)^2)
PDtemp
performance::r2(PD_spmod)

hist(PDData_sum$pd.obs) # data are bounded at 0
qqnorm(residuals(PD_spmod))
testResiduals(simulateResiduals(PD_spmod))
# some deviation

result_1 = data.frame(car::Anova(PD_spmod), R2 = performance::r2(PD_spmod)[1])

empd = emmeans(PD_spmod, ~Species)

# PD diff by species
PD_spmodz = glm(pd.obs.z ~ Species, data=PDData_sum, family="gaussian")
summary(PD_spmodz)
Anova(PD_spmodz)
DHARMa::testResiduals(simulateResiduals(PD_spmodz))
performance::r2(PD_spmodz)

# results
y_hat = predict(PD_spmodz, type="response")
y_bar = mean(PDData_sum$pd.obs.z)
PDztemp = 1-sum((PDData_sum$pd.obs.z-y_hat)^2)/sum((PDData_sum$pd.obs.z-y_bar)^2)
PDztemp

result_2 = data.frame(car::Anova(PD_spmodz), R2 = PDztemp)

empdz= emmeans(PD_spmodz, ~Species)

pdzplot_sp = PDData_sum %>% 
  mutate(Species = recode(Species, `Grevys zebra` = "Grevy's zebra", `Grants gazelle`="Grant's gazelle", `DikDik`="Dik-dik")) %>% 
  ggplot(aes(x=reorder(Species,-pd.obs.z,na.rm=T), y=pd.obs.z))+
  geom_hline(yintercept = -2, linetype="dotted")+
  geom_jitter(width=0.01,height=0.01, alpha=0.5, aes(color=Species))+
  geom_boxplot(aes(fill=Species), color="gray50",outlier.shape=NA)+
  guides(fill="none",col="none")+
  theme_bw(base_size=16)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  scale_color_manual(values=animal_colors)+
  scale_fill_manual(values=animal_colors)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0))+
  annotate(geom="text", x=4.5, y=-5.75, label="Lower than expected PD")+
  labs(x="Host Species",y="sesPD")


pdspplot = gridExtra::grid.arrange(pdplot_sp,pdzplot_sp,ncol=2)

ggsave(here(paste("plots/4_pdsplots",threshold_used,".png")),pdspplot, device="png", dpi=300, width=10, height=5, units="in")


# Models at the species-level
# Setup
# format mammal tree
tree$tip.label[which(tree$tip.label=="Nanger_granti")]="Gazella_granti"
tree$tip.label[which(tree$tip.label=="Tragelaphus_oryx")]="Taurotragus_oryx"
tree$tip.label[which(tree$tip.label=="Equus_africanus")]="Equus_asinus"
tree$tip.label[which(tree$tip.label=="Equus_quagga")]="Equus_burchellii"

#remove unused tips
tree = drop.tip(tree, c("Tragelaphus_scriptus","Kobus_ellipsiprymnus"))

inv.phylo = inverseA(tree, nodes="TIPS", scale=T)

# using default normal for B, inverse wishart for variance parts
prior_normal = list(G=list(G1=list(V=1,nu=0.002),
                           G2=list(V=1,nu=0.002),
                           G3=list(V=1, nu=0.002)),
                    R=list(V=1,nu=0.002))

prior_normal2 = list(G=list(G1=list(V=1,nu=0.002),
                           G2=list(V=1,nu=0.002)),
                    R=list(V=1,nu=0.002))

PDMCMC = PDData_sum %>% 
  filter(is.na(pd.obs.p)==F)
PDMCMC$phylo = PDMCMC$MSW93_Binomial

# exclude cattle
PDMCMC2 = PDMCMC[-which(PDMCMC$Species=="Cattle"),]

# replace missing values with species average for understory
PDMCMC2[which(is.na(PDMCMC2$UNDERSTORY_PROP)),]$UNDERSTORY_PROP=PDMCMC2[which(is.na(PDMCMC2$UNDERSTORY_PROP)),]$UNDERSTORY_SP_MEAN

# exclude RS because it is collinear with BM

# full model
set.seed(123)
model_simple2pd = MCMCglmm(pd.obs.z ~ GUT+ UNDERSTORY_PROP+ GS+ log(BM_KG), random=~phylo+ MSW93_Binomial+ Period,
                           family = "gaussian",
                           ginverse=list(phylo=inv.phylo$Ainv),
                           prior=prior_normal,
                           data=PDMCMC2,
                           nitt=500000, burnin=10000, thin=250, verbose=F, pl=T)

plot(model_simple2pd)

# plot autocorrelation function
plot.acfs = function(x) {
  n <- dim(x)[2]
  par(mfrow=c(ceiling(n/2),2), mar=c(3,2,3,0))
  for (i in 1:n) {
    acf(x[,i], lag.max=100, main=colnames(x)[i])
    grid()
  }
}

plot.acfs(model_simple2pd$VCV)
plot.acfs(model_simple2pd$Sol)

# view results
summary(model_simple2pd)

# non phylo
set.seed(123)
model_pd_2 = MCMCglmm(pd.obs.z ~ GUT+ UNDERSTORY_PROP+ GS+ log(BM_KG),
                           random=~MSW93_Binomial+ Period,
                           family = "gaussian",
                           prior=prior_normal2,
                           data=PDMCMC2,
                           nitt=500000, burnin=10000, thin=250, verbose=F, pl=T)

plot(model_pd_2)
plot.acfs(model_pd_2$VCV)
plot.acfs(model_pd_2$Sol)

# view results
summary(model_pd_2)


# Calculate lambda
lambda_pdz = model_simple2pd$VCV[,'phylo']/
  (model_simple2pd$VCV[,'phylo']+model_simple2pd$VCV[,'units']+model_simple2pd$VCV[,'MSW93_Binomial']+model_simple2pd$VCV[,'Period'])
plot(density(lambda_pdz))
posterior.mode(lambda_pdz)
mean(lambda_pdz)

MCMCpdz = as.data.frame(summary(model_simple2pd)$solutions)
MCMCpdz_var = as.data.frame(summary(model_simple2pd)$Gcovariances)
lambda_pd_z = c(mean = mean(lambda_pdz), mode = posterior.mode(lambda_pdz))

MCMCpdz2 = as.data.frame(summary(model_pd_2)$solutions)
MCMCpdz_var2 = as.data.frame(summary(model_pd_2)$Gcovariances)


pdz_em = PDMCMC2 %>% 
  group_by(Species, GUT) %>% 
  summarize_at(vars(pd.obs.z), funs(mean(., na.rm=T))) %>% 
  group_by(GUT) %>% 
  summarize_at(vars(pd.obs.z), funs(mean, sd, n())) %>% 
  mutate(se=sd/sqrt(n))
  

###
pdz_sp = PDMCMC2 %>% 
  ggplot(aes(x=GUT, y=pd.obs.z))+
  geom_boxplot(aes(col=Species), position=position_dodge(0.8), alpha=0.2)+
  geom_errorbar(data=pdz_em, aes(x=GUT, y=mean, ymin=mean-1.96*se, ymax=mean+1.96*se), size=1)+
  geom_point(data=pdz_em, aes(x=GUT, y=mean), size=2)+
  scale_color_manual(values=animal_colors[-3])+
  theme_bw(base_size=14)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  labs(x="Gut Type", y="PD Z-Value")+guides(col="none")

pdz_sp


#### Just PD
# slow model
set.seed(123)
pd_mod_1 = MCMCglmm(pd.obs ~GUT+UNDERSTORY_PROP+ GS+log(BM_KG),
                    random=~phylo+MSW93_Binomial+Period,
                    family = "gaussian",
                    ginverse=list(phylo=inv.phylo$Ainv),
                    prior=prior_normal,
                    data=PDMCMC2,
                    nitt=500000, burnin=10000, thin=250,
                    verbose=F, pl=T)

plot(pd_mod_1)
plot.acfs(pd_mod_1$VCV)
plot.acfs(pd_mod_1$Sol)
summary(pd_mod_1)


set.seed(123)
pd_mod_2 = MCMCglmm(pd.obs ~GUT+UNDERSTORY_PROP+ GS+log(BM_KG),
                    random=~MSW93_Binomial+Period,
                    family = "gaussian",
                    prior=prior_normal2,
                    data=PDMCMC2,
                    nitt=500000, burnin=10000, thin=250,
                    verbose=F, pl=T)

plot(pd_mod_2)
plot.acfs(pd_mod_2$VCV)
plot.acfs(pd_mod_2$Sol)
summary(pd_mod_2)


# Calculate lambda
lambda_pdobs = pd_mod_1$VCV[,'phylo']/
  (pd_mod_1$VCV[,'phylo']+pd_mod_1$VCV[,'units']+pd_mod_1$VCV[,'MSW93_Binomial']+pd_mod_1$VCV[,'Period'])
plot(density(lambda_pdobs))
posterior.mode(lambda_pdobs)
mean(lambda_pdobs)

MCMCpd = as.data.frame(summary(pd_mod_1)$solutions)
MCMCpd_var = as.data.frame(summary(pd_mod_1)$Gcovariances)
lambda_pd = c(mean = mean(lambda_pdobs), mode = posterior.mode(lambda_pdobs))
MCMCpd2 = as.data.frame(summary(pd_mod_2)$solutions)
MCMCpd_var2 = as.data.frame(summary(pd_mod_2)$Gcovariances)

pd_em = PDMCMC2 %>% 
  group_by(Species, GUT) %>% 
  summarize_at(vars(pd.obs), funs(mean(., na.rm=T))) %>% 
  group_by(GUT) %>% 
  summarize_at(vars(pd.obs), funs(mean, sd, n())) %>% 
  mutate(se=sd/sqrt(n))

###
pd_sp = PDMCMC2 %>% 
  ggplot(aes(x=GUT, y=pd.obs))+
  geom_boxplot(aes(col=Species), position=position_dodge(0.8), alpha=0.2)+
  geom_errorbar(data=pd_em, aes(x=GUT, y=mean, ymin=mean-1.96*se, ymax=mean+1.96*se), size=1)+
  geom_point(data=pd_em, aes(x=GUT, y=mean), size=2)+
  scale_color_manual(values=animal_colors[-3])+
  theme_bw(base_size=14)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  labs(x="Gut Type", y="Faith's PD")+guides(col="none")

pd_sp

pdspplot = gridExtra::grid.arrange(pd_sp,pdz_sp,ncol=2)

ggsave(here(paste("plots/4_pd_pdz_sp_plots",threshold_used,".png",sep="")),pdspplot, device="png", dpi=300, width=10, height=5, units="in")



##### New plot
# combine datasets
PDMCMC2
head(sp_sum)
plot_data = PDMCMC2 %>% 
  group_by(Species,GUT) %>% 
  summarize_at(vars(pd.obs, pd.obs.z), funs(mean)) %>% 
  ungroup() %>% 
  left_join(sp_sum) %>% 
  dplyr::select(MSW93_Binomial, Species, GUT, Prevalence, Richness, pd.obs, pd.obs.z)

# manually add cattle data
# 0.02
plot_data=rbind(plot_data,data.frame(MSW93_Binomial="Bos_taurus",
           Species="Cattle",
           GUT="FG",
           Prevalence=0.57,
           Richness=4.71,
           pd.obs=0.39,
           pd.obs.z=-0.42))
# 0.002
plot_data=rbind(plot_data,data.frame(MSW93_Binomial="Bos_taurus",
                                     Species="Cattle",
                                     GUT="FG",
                                     Prevalence=0.57,
                                     Richness=11.7,
                                     pd.obs=1.00,
                                     pd.obs.z=-0.43))

p = ggtree(tree)
p2 <- p %<+%
  plot_data +
  geom_tiplab(aes(color=GUT), size=2)+
  scale_color_manual(values=c("aquamarine4", "deepskyblue4"))

pp=p+geom_facet(panel="Prevalence",
             data=plot_data,
             geom=geom_tile,
             mapping=aes(x=1, fill=Prevalence))+
  scale_fill_gradient(low="black",high="#F2300F")
pr=p+geom_facet(panel="Richness",
             data=plot_data,
             geom=geom_tile,
             mapping=aes(x=1, fill=Richness))+
  scale_fill_gradient(low="black",high="#F2AD00")
ppd=p+geom_facet(panel="PD",
             data=plot_data,
             geom=geom_tile,
             mapping=aes(x=1, fill=pd.obs))+
  scale_fill_gradient(low="black", high="#81A88D")
ppdz=p+geom_facet(panel="sesPD",
                  data=plot_data,
                  geom=geom_tile,
                  mapping=aes(x=1, fill=pd.obs.z))+
  scale_fill_gradient(low="black", high="#46ACC8")

all = gridExtra::grid.arrange(pp, pr, ppd, ppdz, ncol=4)

ggsave("plots/fig1_tree0.02.png",p2, height=8, width=4)
ggsave("plots/fig1_metrics0.02.png",all, height=8, width=10)

ggsave("plots/fig1_tree0.002.png",p2, height=8, width=4)
ggsave("plots/fig1_metrics0.002.png",all, height=8, width=10)


#############################################################

# write out information
result_2_3 = rbind(result_1, result_2)
write.table(result_2_3, here(paste("docs/4_anova_tests_pd",threshold_used,".txt", sep="")),sep="\t")


### species summary
# species summary
pd_result = PDData_sum %>% 
  group_by(Species) %>% 
  summarize_at(vars(pd.obs, pd.obs.z), funs(mean(., na.rm=T)))
samp_exluded= PDData_sum %>% 
  filter(is.na(pd.obs)) %>% 
  group_by(Species) %>% 
  summarize(n())

write.csv(pd_result, here(paste("docs/4_pd_summary",threshold_used,".csv", sep="")), row.names=F)
write.table(samp_exluded, here(paste("docs/4_pd_excluded",threshold_used,".txt",sep="")), sep="\t")
