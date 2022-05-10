
# 2. Summary Information and Alpha Diversity
# 6 October 2021
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

# A. read in RRA tables and combine with metadata

table_1 = read.csv(here("data/RRA_table_1r.csv"))
table_2 = read.csv(here("data/RRA_table_2r.csv"))
head(table_1)[,1:10]

# read in host information table
host_info = read.csv(here("data/host_metadata.csv"))
names(host_info)[1]="Sample"

# calculate richness
table_1$richness = specnumber(table_1[,-c(1,2)])
table_2$richness = specnumber(table_2[,-c(1,2)])


# join to host information
tab1 = full_join(host_info, table_1, by="Sample")
tab2 = full_join(host_info, table_2, by="Sample")

# sequenced yes, but not in table?
tab1 %>% 
  filter(rd2_sequenced=="Y") %>% 
  filter(is.na(mOTU_1)==T) %>% 
  nrow()
# these are due to filtering steps reads

# mOTU results but no host info -- should be zero!
tab1 %>% 
  filter(is.na(Location)==T)

dim(tab1);dim(tab2)
# add 27 sheep and goats for full sample size (1 hippo, 1 hybrid zebra, 2 waterbuck are still included)

# use foregut fermentation for hippos and camels
tab1 = tab1 %>% 
  mutate_at(vars(GUT), funs(ifelse(.=="FG*", "FG", .)))
tab2 = tab2 %>% 
  mutate_at(vars(GUT), funs(ifelse(.=="FG*", "FG", .)))

# correlation between table 1 and table 2 richness
cor.test(tab1$richness, tab2$richness)


##### Dataset-specific ############# 

# decide which data frame
tab = tab1
# type in the threshold as a character to use for saving plots correctly
threshold_used = "0.002" 

tab = tab %>% 
  filter(Species !="Hybrid zebra") %>% 
  filter(Species != "Waterbuck") %>% 
  filter(Sample != "MRC_17_HIP_93")

# sample sizes and facts


# 1. qPCR n
# plus 27 sheep and goats
fact_1 = tab %>% 
  filter(rd_2_qpcr_conducted=="Y") %>% 
  nrow() + 27

# 2. number with Ct < 35
fact_2 = tab %>% 
  filter(rd_2_qpcr_present=="Y") %>% 
  nrow()

# 3. number selected for metabarcoding from qPCR
fact_3 = tab %>% 
  filter(rd_2_qpcr_present=="Y")%>% 
  filter(rd2_sequenced == "Y") %>% 
  nrow()

# 4. number re-amplified and used in metabarcoding
fact_4 = tab %>% 
  filter(rd2_sequenced=="Y") %>% 
  nrow()

# 5. mean sample read depth after filters
fact_5 = mean(tab$Filtered_Reads, na.rm=T)

# 6. number of samples that dropped out
fact_6 = tab %>% 
  filter(rd2_sequenced=="Y") %>% 
  filter(is.na(mOTU_1)==T) %>% 
  nrow()

# 7 final mOTU table dimensions
fact_7 = tab %>% 
  filter(is.na(mOTU_1)==F) %>%
  select(mOTU_1:mOTU_94) %>% 
  dim(.)

fact_8 = tab %>% 
  filter(Species == "Cattle") %>% 
  nrow()

fact_9 = tab %>% 
  filter(is.na(mOTU_1)==F) %>% 
  filter(Species == "Cattle") %>% 
  nrow()

# fact summary:
fact_summary = data.frame(Fact = c("qPCR N", "Ct <35 N", "N selected from qPCR", "N metabarcoded", "Depth", "N Dropped", "N pos","N mOTU", "N cattle", "N cattle RRA"),
           Stat = round(c(fact_1, fact_2, fact_3, fact_4, fact_5, fact_6, fact_7, fact_8, fact_9),0))
write.table(fact_summary,here(paste("docs/2_sample_facts_", threshold_used,".txt", sep="")), sep="\t", row.names=F)




#### Richness and Prevalence ####
# Individual Richness #

# set animal colors
levels(as.factor(tab$Species))
animal_colors=c("darkorchid4","goldenrod1","blueviolet", "deepskyblue3", "hotpink",  "dodgerblue","green3",  "goldenrod","dodgerblue3", "maroon1",  "deepskyblue2", "greenyellow", "dodgerblue4", "lightskyblue","lightskyblue1", "maroon3", "green4", "cyan2")

tab_present = tab %>% 
  filter(is.na(mOTU_1)==F)

richgraph = tab_present %>% 
  mutate(Species = recode(Species, `Grevys zebra` = "Grevy's zebra", `Grants gazelle`="Grant's gazelle", `DikDik`="Dik-dik")) %>% 
  ggplot(aes(x=reorder(Species,-richness,na.rm=T), y=richness))+
  geom_jitter(width=0.01,height=0.3, alpha=0.5, aes(color=Species))+
  geom_boxplot(aes(fill=Species), color="gray50",outlier.shape=NA)+
  guides(fill="none",col="none")+
  theme_bw(base_size=16)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  scale_fill_manual(values=animal_colors)+
  scale_color_manual(values=animal_colors)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0))+
  scale_y_continuous(limits = c(0,(max(tab_present$richness)+1)))+
  labs(x="Host Species",y="Richness")
richgraph


# Prevalence #

# calculate prevalence with CI
tab_prev = tab %>% 
  mutate_at(vars(rd_2_qpcr_present), funs(ifelse(.=="Y",1,0))) %>% 
  group_by(Species) %>% 
  summarize_at(vars(rd_2_qpcr_present), funs(n=n(), pos=sum(.))) %>% 
  ungroup()

tab_prev = binom.confint(tab_prev$pos, tab_prev$n, methods="wilson")%>%
  cbind(tab_prev[,1])


prevgraph = tab_prev %>% 
  mutate(Species = recode(Species, `Grevys zebra` = "Grevy's zebra", `Grants gazelle`="Grant's gazelle", `DikDik`="Dik-dik")) %>% 
  ggplot(aes(x=reorder(Species, -mean), y=mean))+
  geom_point(aes(color=Species), size=2)+
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Species), size=1.2)+
  guides(col="none")+
  theme_bw(base_size = 16)+
  theme(panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank())+
  scale_color_manual(values=animal_colors)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0))+
  scale_y_continuous(limits=c(0,1))+
  labs(x="Host Species", y="Prevalence")
prevgraph



# Richness and Prevalence comparison
compare = tab_present %>% 
  group_by(Species) %>% 
  summarize_at(vars(richness), funs(mean)) %>% 
  arrange(desc(richness)) %>% 
  left_join(tab_prev)
temp = cor.test(compare$mean, compare$richness, method = "spearman")

result_1 = unlist(temp) %>% as.data.frame()

prev_rich_graph = gridExtra::grid.arrange(prevgraph, richgraph, ncol=2)
ggsave(here(paste("plots/2_prevalence_richness_",threshold_used,".png", sep="")), prev_rich_graph, width=10, height=5, dpi=300, device="png")


### qPCR value
head(tab)
ggplot(tab, aes(x=reorder(Species,rd_2_qpcr_mean_ct), y=rd_2_qpcr_mean_ct))+
  geom_boxplot(aes(fill=Species))+
  scale_fill_manual(values=animal_colors)+
  guides(fill="none")+
  theme_bw(base_size = 14)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0))

#### How much variation is explained by species?

## Prevalence #################
# Test prevalence
MP_prev = tab %>% 
  mutate_at(vars(rd_2_qpcr_present), funs(ifelse(.=="Y",1,0)))

# fit logistic with species
MPprev_mod = glmmTMB(rd_2_qpcr_present ~ Species, data=MP_prev, family="binomial")
# residuals
testResiduals(simulateResiduals(MPprev_mod))

# assess pct variation
temp = car::Anova(MPprev_mod)
as.data.frame(performance::r2_tjur(MPprev_mod))

# save results and sample size
result_2 = cbind(as.data.frame(temp), data.frame(R2 = performance::r2_tjur(MPprev_mod)))
dim(MP_prev)

## Richness ################
# fit model
dim(tab_present)
MPsp = glmmTMB(richness ~ Species, family="poisson", data=tab_present)

# diagnostics
testResiduals(simulateResiduals(MPsp))
temp = car::Anova(MPsp)

# results
y_hat = predict(MPsp, type="response")
y_bar = mean(tab_present$richness)
rtemp = 1-sum((tab_present$richness-y_hat)^2)/sum((tab_present$richness-y_bar)^2)

result_3 = cbind(as.data.frame(temp), data.frame(R2 = rtemp))

## Species-level analyses ##################################################

## Setup ##
# Fit a model with a species term
tab$phylo = tab$MSW93_Binomial
MP_prev$phylo = MP_prev$MSW93_Binomial
tab_present$phylo = tab_present$MSW93_Binomial


# format tree
tree = read.tree(here("data/new_mammal_tree_pruned.newick"))

# rename tips to match phylogeny
tree$tip.label[which(tree$tip.label=="Nanger_granti")]="Gazella_granti"
tree$tip.label[which(tree$tip.label=="Tragelaphus_oryx")]="Taurotragus_oryx"
tree$tip.label[which(tree$tip.label=="Equus_africanus")]="Equus_asinus"
tree$tip.label[which(tree$tip.label=="Equus_quagga")]="Equus_burchellii"

# remove unused tips
tree = drop.tip(tree, c("Tragelaphus_scriptus","Kobus_ellipsiprymnus"))

# format phylo for use with MCMCglmm
inv.phylo = inverseA(tree, nodes="TIPS", scale=T)

# set prior for poisson: weak inverse Wishart for residual variance and G structure

# Note: parameter expanded priors for G structure (as in Hadfield notes)
# gave unacceptable autocorrelation values for phylogeny and species

# with phylogeny
## inverse wishart
prior_pois = list(G=list(G1=list(V=1,nu=0.002),
                         G2=list(V=1,nu=0.002),
                         G3=list(V=1,nu=0.002)),
                  R=list(V=1,nu=0.002))

# without phylogeny
## inverse wishart
prior_pois2 = list(G=list(G1=list(V=1,nu=0.002),
                         G2=list(V=1,nu=0.002)),
                  R=list(V=1,nu=0.002))

#### Richness -MCMCglmm #####

# format data correctly
MPrichinfoMCMC = as.data.frame(tab_present)
MPrichinfoMCMC$animal = MPrichinfoMCMC$Sample_ID
MPrichinfoMCMC$phylo = MPrichinfoMCMC$MSW93_Binomial
MPrichinfoMCMC = as.data.frame(MPrichinfoMCMC)

# exclude cattle because they were treated with anthelminthics
MPrichinfoMCMC = MPrichinfoMCMC[-which(MPrichinfoMCMC$Species=="Cattle"),]

# use group mean for missing sample-level data for understory
MPrichinfoMCMC[which(is.na(MPrichinfoMCMC$UNDERSTORY_PROP)),]$UNDERSTORY_PROP=MPrichinfoMCMC[which(is.na(MPrichinfoMCMC$UNDERSTORY_PROP)),]$UNDERSTORY_SP_MEAN

dim(MPrichinfoMCMC) 
# 257 samples (281 - 24 cattle)

# double-check variable correlations
MPrichinfoMCMC %>% 
  ggplot(aes(x=log(BM_KG), y=log(RS_KM2)))+
  geom_point()+
  geom_smooth(method="lm", se=F);
cor.test(MPrichinfoMCMC$BM_KG, MPrichinfoMCMC$RS_KM2, method="spearman")
# Body mass and range size are closely related; will exclude range size

################
# Fit model

# with phylogeny
set.seed(123)
model_rich = MCMCglmm(richness ~ log(BM_KG)+ GS+ GUT +UNDERSTORY_PROP,
                      random= ~ phylo + MSW93_Binomial+Period,
                      family = "poisson",
                      ginverse=list(phylo=inv.phylo$Ainv),
                      prior=prior_pois,
                      data=MPrichinfoMCMC,
                      nitt=500000, burnin=10000, thin=250,
                      verbose=F, pl=T)

# Visually inspect mixing and posterior distributions
plot(model_rich)

# plot autocorrelation
plot.acfs = function(x) {
  n <- dim(x)[2]
  par(mfrow=c(ceiling(n/2),2), mar=c(3,2,3,0))
  for (i in 1:n) {
    acf(x[,i], lag.max=100, main=colnames(x)[i])
    grid()
  }
}

plot.acfs(model_rich$VCV)
plot.acfs(model_rich$Sol)

# view results
summary(model_rich)


###################
# without phylogeny
set.seed(123)
model_rich2 = MCMCglmm(richness ~ log(BM_KG)+ GS+ GUT +UNDERSTORY_PROP,
                       random= ~ MSW93_Binomial+Period,
                       family = "poisson",
                       prior=prior_pois2,
                       data=MPrichinfoMCMC,
                       nitt=500000, burnin=10000, thin=250,
                       verbose=F, pl=T)

# Visually inspect mixing and posterior
plot(model_rich2)

# Check autocorrelations
plot.acfs(model_rich2$VCV)
plot.acfs(model_rich2$Sol)

# view results
summary(model_rich2)

# store lambda
lambda_r = model_rich$VCV[,'phylo']/
  (model_rich$VCV[,'phylo']+model_rich$VCV[,'units']+model_rich$VCV[,'MSW93_Binomial']+model_rich$VCV[,'Period'])
posterior.mode(lambda_r)
mean(lambda_r)

# store results
MCMCrichness = as.data.frame(summary(model_rich)$solutions)
MCMCrichness_var = as.data.frame(summary(model_rich)$Gcovariances)
lambda_rich = c(mean = mean(lambda_r), mode = posterior.mode(lambda_r))

MCMCrichness2 = as.data.frame(summary(model_rich2)$solutions)
MCMCrichness_var2 = as.data.frame(summary(model_rich2)$Gcovariances)


### Compare to glmmTMB ###
# fit model
rmod = glmmTMB(richness ~  log(BM_KG)+ GS+ GUT +UNDERSTORY_PROP+ (1|Period)+(1|MSW93_Binomial), family="poisson", data=MPrichinfoMCMC)
summary(rmod)

drmod = glmmTMB(richness ~ GUT+(1|Period)+(1|MSW93_Binomial), family="poisson", data=MPrichinfoMCMC)
summary(drmod)
GUTests = as.data.frame(emmeans(rmod, ~GUT, type="response"))


# Plot
MPrichinfoMCMC %>% 
  group_by(Species, GS, BM_KG, RS_KM2, UNDERSTORY_SP_MEAN) %>% 
  summarize_at(vars(richness), funs(mean, sd, n())) %>% 
  mutate(se = sd/sqrt(n)) %>% 
  ungroup() %>% 
  mutate_at(vars(BM_KG, RS_KM2), funs(log)) %>% 
  pivot_longer(cols=GS:UNDERSTORY_SP_MEAN, names_to="Variable", values_to="Value") %>% 
  ggplot(aes(x=Value, y=mean))+
  geom_smooth(method="glm", method.args=list(family="poisson"),se=F, col="gray50")+
  geom_point(aes(col=Species), size=1.5)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, col=Species),size=1)+
  theme_bw(base_size=14)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  scale_color_manual(values=animal_colors[-3])+
  facet_wrap(~Variable, scales="free_x")+
  labs(y="Parasite mOTU Richness")

GUTests = MPrichinfoMCMC %>% 
  group_by(Species, GUT) %>% 
  summarize_at(vars(richness), funs(mean)) %>% 
  group_by(GUT) %>% 
  summarize_at(vars(richness), funs(mean, sd, n())) %>% 
  mutate(se=sd/sqrt(n))

rich_sp_plot = MPrichinfoMCMC %>% 
  ggplot(aes(x=GUT, y=richness))+
  geom_boxplot(aes(col=Species), position=position_dodge(0.8), alpha=0.2)+
  geom_errorbar(data=GUTests, aes(x=GUT, y=mean, ymin=mean-1.96*se, ymax=mean+1.96*se), size=1)+
  geom_point(data=GUTests, aes(x=GUT, y=mean), size=2)+
  scale_color_manual(values=animal_colors[-3])+
  theme_bw(base_size=14)+
  guides(col="none")+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  labs(x="Gut Type", y="Richness")


#### Prevalence MCMCglmm #####

# format data correctly
MPprevinfoMCMC = as.data.frame(MP_prev)
MPprevinfoMCMC$animal = MPprevinfoMCMC$Sample_ID
MPprevinfoMCMC$phylo = MPprevinfoMCMC$MSW93_Binomial
MPprevinfoMCMC = as.data.frame(MPprevinfoMCMC)
dim(MPprevinfoMCMC)

#exclude cattle
MPprevinfoMCMC=MPprevinfoMCMC[-which(MPprevinfoMCMC$Species=="Cattle"),]

# use group mean
MPprevinfoMCMC[which(is.na(MPprevinfoMCMC$UNDERSTORY_PROP)),]$UNDERSTORY_PROP=MPprevinfoMCMC[which(is.na(MPprevinfoMCMC$UNDERSTORY_PROP)),]$UNDERSTORY_SP_MEAN
dim(MPprevinfoMCMC)
# 519 - 49 cattle

MPprevinfoMCMC$L_BM_KG = log(MPprevinfoMCMC$BM_KG)

# Binomial priors -- parameter expanded for G component, fixed variance for residual variance
prior_binom = list(R = list(V = 1, fix = T),
                   G = list(G1 =list(V = 1, nu = 0.002),
                            G2 =list(V = 1, nu = 0.002),
                            G3 =list(V = 1, nu = 0.002)))

# no phylogeny version:
prior_binom2 = list(R = list(V = 1, fix = T),
                   G = list(G1 =list(V = 1, nu = 1, alpha.mu=0, alpha.V=1000),
                            G2 =list(V = 1, nu = 1, alpha.mu=0, alpha.V=1000))) 
             
# Fit full model
set.seed(123)
modt_1 = MCMCglmm(fixed = rd_2_qpcr_present ~ L_BM_KG + GS+ GUT +UNDERSTORY_PROP,
                 random = ~phylo+MSW93_Binomial+Period,
                 ginverse = list(phylo = inv.phylo$Ainv),
                 data = MPprevinfoMCMC,
                 family = "categorical",
                 prior = prior_binom,
                 nitt = 500000, burnin = 10000, thin = 250,
                 verbose = F, pl = T)

# plot model
plot(modt_1)

# check autocorrelation
plot.acfs(modt_1$VCV)
plot.acfs(modt_1$Sol)

# view results
summary(modt_1)


# without phylo
set.seed(123)
modt_2 = MCMCglmm(fixed = rd_2_qpcr_present ~ L_BM_KG + GS+ GUT +UNDERSTORY_PROP,
                 random = ~MSW93_Binomial+Period,
                 data = MPprevinfoMCMC,
                 family = "categorical",
                 prior = prior_binom2,
                 nitt = 500000, burnin = 10000, thin = 250,
                 verbose = F, pl = T)

# plot model
plot(modt_2)

# check autocorrelation
plot.acfs(modt_2$VCV)
plot.acfs(modt_2$Sol)

# view results
summary(modt_2)

# Calculate lambda
lambda_p = modt_1$VCV[,'phylo']/
  (modt_1$VCV[,'phylo']+modt_1$VCV[,'units']+modt_1$VCV[,'MSW93_Binomial']+modt_1$VCV[,"Period"])
posterior.mode(lambda_p,0.25)
mean(lambda_p)

MCMCprevalence = as.data.frame(summary(modt_1)$solutions)
MCMCprevalence_var = as.data.frame(summary(modt_1)$Gcovariances)
lambda_prev = c(mean = mean(lambda_p), mode=posterior.mode(lambda_p, 0.25))

MCMCprevalence2 = as.data.frame(summary(modt_2)$solutions)
MCMCprevalence_var2 = as.data.frame(summary(modt_2)$Gcovariances)

dim(MPprevinfoMCMC)
#470 samples


### Prevalence -glmmTMB ###
model_prevalence = glmmTMB(rd_2_qpcr_present ~ log(BM_KG)+ GS+ GUT +UNDERSTORY_PROP+ (1|MSW93_Binomial)+(1|Period), family="binomial", data=MPprevinfoMCMC)
summary(model_prevalence)
# Save values
TMBprevalence = as.data.frame(summary(model_prevalence)$coefficients$cond)
TMBprevalence_var = unlist(summary(model_prevalence)$varcor$cond)

MPdata_prev2 = left_join(tab_prev,unique(dplyr::select(MPrichinfoMCMC,Species,GUT,GS)))

MPdata_prev2[which(MPdata_prev2$Species=="Cattle"),]$GUT="FG"

# Plot gut type
prevplotgp = MPdata_prev2 %>% 
  group_by(GUT) %>% 
  summarize_at(vars(mean), funs(sd, mean, n())) %>% 
  mutate(se = sd/sqrt(n)) %>% 
  ggplot(aes(x=GUT, y=mean))+
  geom_errorbar(data=MPdata_prev2, aes(x=GUT, y=mean, ymin=lower, ymax=upper, col=Species), position = position_dodge(width=0.5), width=0.5)+
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), size=1)+
  geom_point(size=2)+
  scale_color_manual(values=animal_colors)+
  theme_bw()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  guides(col="none")+
  labs(x="Gut Type", y="Prevalence")


gpplots = gridExtra::grid.arrange(prevplotgp, rich_sp_plot, ncol=2)

ggsave(here(paste("plots/2_sp_prev_rich",threshold_used,".png", sep="")), gpplots, device="png", dpi=300, width=10, height=5, units="in")



# write out information
result_2_3 = rbind(result_2, result_3)

write.table(result_1, here(paste("docs/2_rich_prev_correlation",threshold_used,".txt", sep="")), sep="\t")
write.table(result_2_3, here(paste("docs/2_anova_tests",threshold_used,".txt", sep="")), sep="\t")


# species summary
prev_result = tab %>% 
  group_by(Species,MSW93_Binomial) %>% 
  mutate_at(vars(rd_2_qpcr_present), funs(ifelse(.=="Y",1,0))) %>% 
  summarize_at(vars(rd_2_qpcr_present), funs(mean,n())) %>% 
  rename(Prevalence = mean, Prev_N = n)
rich_result = tab %>% 
  filter(is.na(mOTU_1)==F) %>% 
  group_by(Species) %>% 
  summarize_at(vars(richness), funs(mean,n())) %>% 
  rename(Richness = mean, Rich_N = n)

sp_sum = left_join(prev_result, rich_result)

write.csv(sp_sum, here(paste("docs/2_prev_rich_summary",threshold_used,".csv", sep="")), row.names=F)
