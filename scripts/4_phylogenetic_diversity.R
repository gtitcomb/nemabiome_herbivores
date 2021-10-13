
# Phylogenetic Diversity

# this script depends on scripts 1 and 3
# Use CTRL+SHIFT+F10 to clear

########

library(picante)
library(car)
library(DHARMa)
library(MCMCglmm)
library(glmmTMB)
library(tidyverse)
library(adegenet)
library(emmeans)

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
table_1 = read.csv(here("data/RRA_table_1.csv"))
table_2 = read.csv(here("data/RRA_table_2.csv"))

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

png(here(paste("plots_temp/tree",threshold_used,".png",sep="")), res=300, width=8, height=8, units="in")
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
  labs(x="Gut Type", y="PD Z-Value")+
  scale_color_manual(values=c(animal_colors))+
  guides(col="none")
  pdz

pd = ggplot(PDData_sum, aes(x=GUT, y=pd.obs))+
  geom_boxplot(fill=NA, col="gray50", size=1)+
  geom_boxplot(data=PDData_sum, aes(x=GUT, y=pd.obs, col=Species), size=0.8, position=position_dodge(width=0.5), width=0.25, alpha=0.6)+
  theme_bw(base_size=16)+
  labs(x="Diet", y="Faith's PD")+
  scale_color_manual(values=c(animal_colors[-13]))+
  #facet_wrap(~Fermentation, scales="free_x")+
  guides(color="none")
pd

pdplot = gridExtra::grid.arrange(pd, pdz, ncol=2)

ggsave(here(paste("plots/pd_plots_sp",threshold_used,".png", sep="")), pdplot, device="png", dpi=300, width=10, height=5, units="in")


# species, ungrouped

pdplot_sp = PDData_sum %>% 
  mutate(Species = recode(Species, `Grevys zebra` = "Grevy's zebra", `Grants gazelle`="Grant's gazelle", `DikDik`="Dik-dik")) %>% 
  ggplot(aes(x=reorder(Species,-pd.obs,na.rm=T), y=pd.obs))+
  geom_jitter(width=0.01,height=0.01, alpha=0.5, aes(color=Species))+
  geom_boxplot(aes(fill=Species), color="gray50",outlier.shape=NA)+
  guides(fill="none",col="none")+
  theme_bw(base_size=16)+
  scale_color_manual(values=(animal_colors))+
  scale_fill_manual(values=(animal_colors))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0))+
  labs(x="Host Species",y="PD")
pdplot_sp

# total PD by species
PD_spmod = lm(pd.obs~ Species, data=PDData_sum)
summary(PD_spmod)
Anova(PD_spmod)
performance::r2(PD_spmod)[2]
testResiduals(simulateResiduals(PD_spmod)) # not fantastic

result_1 = data.frame(car::Anova(PD_spmod), R2 = performance::r2(PD_spmod)[2])

empd = emmeans(PD_spmod, ~Species)

# PD diff by species
PD_spmodz = lm(pd.obs.z ~ Species, data=PDData_sum)
summary(PD_spmodz)
DHARMa::testResiduals(simulateResiduals(PD_spmodz))

result_2 = data.frame(car::Anova(PD_spmodz), R2 = performance::r2(PD_spmodz)[2])

empdz= emmeans(PD_spmodz, ~Species)

pdzplot_sp = PDData_sum %>% 
  mutate(Species = recode(Species, `Grevys zebra` = "Grevy's zebra", `Grants gazelle`="Grant's gazelle", `DikDik`="Dik-dik")) %>% 
  ggplot(aes(x=reorder(Species,-pd.obs.z,na.rm=T), y=pd.obs.z))+
  geom_hline(yintercept = -2, linetype="dotted")+
  geom_jitter(width=0.01,height=0.01, alpha=0.5, aes(color=Species))+
  geom_boxplot(aes(fill=Species), color="gray50",outlier.shape=NA)+
  guides(fill="none",col="none")+
  theme_bw(base_size=16)+
  scale_color_manual(values=animal_colors)+
  scale_fill_manual(values=animal_colors)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0))+
  annotate(geom="text", x=4.5, y=-5.75, label="Lower than expected PD")+
  labs(x="Host Species",y="sesPD")


pdspplot = gridExtra::grid.arrange(pdplot_sp,pdzplot_sp,ncol=2)

ggsave(here(paste("plots/pdsplots",threshold_used,".png")),pdspplot, device="png", dpi=300, width=10, height=5, units="in")


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

prior_normal = list(G=list(G1=list(V=1,nu=0.02),
                           G2=list(V=1,nu=0.02),
                           G3=list(V=1, nu=0.02)),
                    R=list(V=1,nu=0.02))

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

# Calculate lambda
lambda_pdz = model_simple2pd$VCV[,'phylo']/
  (model_simple2pd$VCV[,'phylo']+model_simple2pd$VCV[,'units']+model_simple2pd$VCV[,'MSW93_Binomial']+model_simple2pd$VCV[,'Period'])
plot(density(lambda_pdz))
posterior.mode(lambda_pdz)
mean(lambda_pdz)

MCMCpdz = as.data.frame(summary(model_simple2pd)$solutions)
MCMCpdz_var = as.data.frame(summary(model_simple2pd)$Gcovariances)
lambda_pd_z = c(mean = mean(lambda_pdz), mode = posterior.mode(lambda_pdz))



PDZglmmtmb = glmmTMB(pd.obs.z ~ GUT+ UNDERSTORY_PROP+ GS+ log(BM_KG) + (1|MSW93_Binomial)+ (1|Period), data=PDMCMC2)

summary(PDZglmmtmb)
DHARMa::testResiduals(simulateResiduals(PDZglmmtmb))

# Save values
TMBpdz= as.data.frame(summary(PDZglmmtmb)$coefficients$cond)
TMBpdz_var = as.data.frame(unlist(summary(PDZglmmtmb)$varcor$cond))

pdz_em = as.data.frame(emmeans(PDZglmmtmb, ~ GUT, type="response"))

###
pdz_sp = PDMCMC2 %>% 
  ggplot(aes(x=GUT, y=pd.obs.z))+
  geom_boxplot(aes(col=Species), position=position_dodge(0.8), alpha=0.2)+
  geom_errorbar(data=pdz_em, aes(x=GUT, y=emmean, ymin=lower.CL, ymax=upper.CL), size=1)+
  geom_point(data=pdz_em, aes(x=GUT, y=emmean), size=2)+
  scale_color_manual(values=animal_colors[-3])+
  theme_bw(base_size=14)+
  labs(x="Gut Type", y="PD Z-Value")+guides(col="none")

pdz_sp


#### Just PD
# slow model
set.seed(123)
model_simple2pd.obs = MCMCglmm(pd.obs ~GUT+UNDERSTORY_PROP+ GS+log(BM_KG),
                               random=~phylo+MSW93_Binomial+Period,
                               family = "gaussian",
                               ginverse=list(phylo=inv.phylo$Ainv),
                               prior=prior_normal, data=PDMCMC2,
                               nitt=500000, burnin=10000, thin=250, verbose=F, pl=T)
summary(model_simple2pd.obs)

# Calculate lambda
lambda_pdobs = model_simple2pd.obs$VCV[,'phylo']/
  (model_simple2pd.obs$VCV[,'phylo']+model_simple2pd.obs$VCV[,'units']+model_simple2pd.obs$VCV[,'MSW93_Binomial']+model_simple2pd.obs$VCV[,'Period'])
plot(density(lambda_pdobs))
posterior.mode(lambda_pdobs)
mean(lambda_pdobs)

MCMCpd = as.data.frame(summary(model_simple2pd.obs)$solutions)
MCMCpd_var = as.data.frame(summary(model_simple2pd.obs)$Gcovariances)
lambda_pd = c(mean = mean(lambda_pdobs), mode = posterior.mode(lambda_pdobs))


PDglmmtmb = glmmTMB(pd.obs ~ GUT+ UNDERSTORY_PROP+ GS+ log(BM_KG)+ (1|MSW93_Binomial)+(1|Period), data=PDMCMC2)
summary(PDglmmtmb)

# Save values
TMBpd= as.data.frame(summary(PDglmmtmb)$coefficients$cond)
TMBpd_var = as.data.frame(unlist(summary(PDglmmtmb)$varcor$cond))

pd_em = as.data.frame(emmeans(PDglmmtmb, ~ GUT, type="response"))

###
pd_sp = PDMCMC2 %>% 
  ggplot(aes(x=GUT, y=pd.obs))+
  geom_boxplot(aes(col=Species), position=position_dodge(0.8), alpha=0.2)+
  geom_errorbar(data=pd_em, aes(x=GUT, y=emmean, ymin=lower.CL, ymax=upper.CL), size=1)+
  geom_point(data=pd_em, aes(x=GUT, y=emmean), size=2)+
  scale_color_manual(values=animal_colors[-3])+
  theme_bw(base_size=14)+
  labs(x="Gut Type", y="Faith's PD")+guides(col="none")

pd_sp

pdspplot = gridExtra::grid.arrange(pd_sp,pdz_sp,ncol=2)

ggsave(here(paste("plots/pd_pdz_sp_plots",threshold_used,".png",sep="")),pdspplot, device="png", dpi=300, width=10, height=5, units="in")



# write out information
result_2_3 = rbind(result_1, result_2)
write_delim(result_2_3, here(paste("docs/anova_tests_pd",threshold_used,".txt", sep="")))

MCMCpd$predictor = row.names(MCMCpd)
MCMCpdz$predictor = row.names(MCMCpdz)
TMBpd$predictor = row.names(TMBpd)
TMBpdz$predictor = row.names(TMBpdz)
MCMCpd_var$predictor = row.names(MCMCpd_var)
MCMCpdz_var$predictor = row.names(MCMCpdz_var)

MCMCmods = rbind(MCMCpd, MCMCpdz)
MCMCmods$Metric = c(rep("PD",5), rep("sesPD", 5))
MCMCvars = rbind(MCMCpd_var, MCMCpdz_var)
MCMCvars$Metric = c(rep("PD",3), rep("sesPD", 3))
TMBmods = rbind(TMBpd, TMBpdz)
TMBmods$Metric = c(rep("PD",5), rep("sesPD", 5))
TMBvars = cbind(TMBpd_var, TMBpdz_var)
names(TMBvars) = c("PD","sesPD")

## Models MCMC
MCMCmods2 = MCMCmods %>% 
  mutate_at(vars(post.mean:pMCMC), funs(round(.,3))) %>% 
  mutate(CI = paste(`l-95% CI`, `u-95% CI`, sep=", ")) %>% 
  mutate(M = paste(post.mean, " (",pMCMC, ")", sep="")) %>% 
  pivot_longer(cols=c(M,CI), names_to="Measure", values_to="Value") %>% 
  dplyr::select(predictor, Value, eff.samp, Metric)

## Variance
MCMCvars2 = MCMCvars %>% 
  mutate_at(vars(post.mean:eff.samp), funs(round(.,3))) %>% 
  mutate(CI=paste(`l-95% CI`, `u-95% CI`, sep=", ")) %>% 
  mutate_at(vars(post.mean), funs(as.character(.))) %>% 
  pivot_longer(cols=c(post.mean, CI), names_to="Measure", values_to="Value") %>% 
  dplyr::select(predictor,Value, eff.samp, Metric)
MCMCvars2

## TMB
TMBmods2 = TMBmods %>% 
  mutate_at(vars(Estimate:`Pr(>|z|)`), funs(round(.,3))) %>% 
  mutate(M_SE=paste(Estimate, `Std. Error`, sep=", ")) %>% 
  mutate(Z_p = paste(`z value`," (",`Pr(>|z|)`,")",sep="")) %>% 
  pivot_longer(cols=c(M_SE,Z_p), names_to="Measure", values_to="Value") %>% 
  dplyr::select(predictor, Value, Metric)

# Combine into table

allMods = cbind(MCMCmods2[1:10, ], TMBmods2[1:10,], MCMCmods2[11:20,], TMBmods2[11:20,])
# check all variables in correct order
names(allMods)=c("Predictor", "MCMC PD", "eff.samp.prev", "Metric", "P2", "TMB PD", "Metric", "P3",
                 "MCMC sesPD", "eff.samp.rich", "Metric", "P4", "TMB sesPD", "Metric")
allMods = allMods %>% 
  select(Predictor, `MCMC PD`, eff.samp.prev, `TMB PD`, `MCMC sesPD`, eff.samp.rich,`TMB sesPD`)

MCMCvars2
TMBvars

MCMCvars2 = data.frame(MCMCvars2[1:6,c(1:3)], TMB_Prev = c("","","",TMBvars[1,1],"",TMBvars[2,1]),
                       MCMCvars2[7:12,2:3], TMB_Rich =  c("","","",TMBvars[1,1],"",TMBvars[2,1]))
MCMCvars2 = MCMCvars2 %>% 
  mutate_at(vars(TMB_Prev, TMB_Rich), funs(round(as.numeric(.),3))) %>% 
  mutate_at(vars(TMB_Prev, TMB_Rich), funs(ifelse(is.na(.), "",.)))
names(MCMCvars2) = c("Predictor", "MCMC PD", "eff.samp.prev", "TMB PD", "MCMC sesPD", "eff.samp.rich",
                     "TMB sesPD")            

full_results_table = rbind(allMods, MCMCvars2)

write.csv(full_results_table, here(paste("docs/pd_sesPD_model_results",threshold_used,".csv", sep="")), row.names = F)

lambdas = as.data.frame(rbind(lambda_pd, lambda_pd_z))
write_delim(lambdas, here(paste("docs/lambda_pds",threshold_used,".txt",sep="")))
