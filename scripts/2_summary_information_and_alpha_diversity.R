
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

table_1 = read.csv(here("data/RRA_table_1.csv"))
table_2 = read.csv(here("data/RRA_table_2.csv"))
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
tab = tab2
# type in the threshold as a character to use for saving plots correctly
threshold_used = "0.02" 

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

# set animal colors based on alphabet
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

# calculate prevalence with correct errorbars
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



#### How much variation is explained by species?

## Prevalence #################
# Test prevalence
MP_prev = tab %>% 
  mutate_at(vars(rd_2_qpcr_present), funs(ifelse(.=="Y",1,0)))

MPprev_mod = glmmTMB(rd_2_qpcr_present ~ Species, data=MP_prev, family="binomial")
temp = car::Anova(MPprev_mod)
as.data.frame(performance::r2_tjur(MPprev_mod))

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

# rename some tips
tree$tip.label[which(tree$tip.label=="Nanger_granti")]="Gazella_granti"
tree$tip.label[which(tree$tip.label=="Tragelaphus_oryx")]="Taurotragus_oryx"
tree$tip.label[which(tree$tip.label=="Equus_africanus")]="Equus_asinus"
tree$tip.label[which(tree$tip.label=="Equus_quagga")]="Equus_burchellii"

#remove unused tips
tree = drop.tip(tree, c("Tragelaphus_scriptus","Kobus_ellipsiprymnus"))

inv.phylo = inverseA(tree, nodes="TIPS", scale=T) # decide all or tips

# set prior for poisson
prior_pois = list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02), G3=list(V=1,nu=0.02)), R=list(V=1,nu=0.02))


#### Richness -MCMCglmm #####

# format data correctly
MPrichinfoMCMC = as.data.frame(tab_present)
MPrichinfoMCMC$animal = MPrichinfoMCMC$Sample_ID
MPrichinfoMCMC$phylo = MPrichinfoMCMC$MSW93_Binomial
MPrichinfoMCMC = as.data.frame(MPrichinfoMCMC)

# exclude cattle because they were treated
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
MPrichinfoMCMC %>% 
  ggplot(aes(x=log(BM_KG), y=GS))+
  geom_point()+
  geom_smooth(method="lm", se=F)
cor.test(MPrichinfoMCMC$BM_KG, MPrichinfoMCMC$RS_KM2, method="spearman")

# Body mass and range size are closely related; will want to fit separately

# Fit model

set.seed(123)
model_rich = MCMCglmm(richness ~ log(BM_KG)+ GS+ GUT +UNDERSTORY_PROP, random= ~ phylo + MSW93_Binomial+Period, family = "poisson", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_pois, data=MPrichinfoMCMC, nitt=500000, burnin=10000, thin=250, verbose=F, pl=T)


# check autocorrelation
autocorr(model_rich$Sol)
autocorr(model_rich$VCV) 

summary(model_rich)


HPDinterval(mcmc(model_rich$Sol)) %>%
  as.data.frame() %>%
  mutate(predictor=rownames(.)) %>%
  filter(predictor != "(Intercept)") %>%
  ggplot(aes(x=predictor))+
  geom_errorbar(aes(ymin=lower, ymax=upper))+
  geom_hline(yintercept=0)


lambda_r = model_rich$VCV[,'phylo']/
  (model_rich$VCV[,'phylo']+model_rich$VCV[,'units']+model_rich$VCV[,'MSW93_Binomial']+model_rich$VCV[,'Period'])
plot(density(lambda_r))
posterior.mode(lambda_r)

# Save results
MCMCrichness = as.data.frame(summary(model_rich)$solutions)
MCMCrichness_var = as.data.frame(summary(model_rich)$Gcovariances)
lambda_rich = c(mean = mean(lambda_r), mode = posterior.mode(lambda_r))


## Richness - glmmTMB #####
# Full model
rmod = glmmTMB(richness ~  log(BM_KG)+ GS+ GUT +UNDERSTORY_PROP+ (1|Period)+(1|MSW93_Binomial), family="poisson", data=MPrichinfoMCMC)
summary(rmod)


# Save values
TMBrichness = as.data.frame(summary(rmod)$coefficients$cond)
TMBrichness_var = as.data.frame(unlist(summary(rmod)$varcor$cond))


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

rich_sp_plot = MPrichinfoMCMC %>% 
  ggplot(aes(x=GUT, y=richness))+
  geom_boxplot(aes(col=Species), position=position_dodge(0.8), alpha=0.2)+
  geom_errorbar(data=GUTests, aes(x=GUT, y=rate, ymin=lower.CL, ymax=upper.CL), size=1)+
  geom_point(data=GUTests, aes(x=GUT, y=rate), size=2)+
  scale_color_manual(values=animal_colors[-3])+
  theme_bw(base_size=14)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  labs(x="Gut Type", y="Richness")

rich_sp_plot = rich_sp_plot+guides(col="none")
rich_sp_plot



#### Prevalence -MCMCglmm #####

# advice from https://r-sig-mixed-models.r-project.narkive.com/I9vUTpvM/r-sig-me-rare-binary-outcome-mcmcglmm-and-priors
prior_binom = list(R = list(V = 1, fix = 1),
                   B = list(mu = c(rep(0,6)),
                            V = diag(6)*(4.3+pi^2/3)),
                   G = list(G1 =list(V = 1, nu = 0),
                            G1 =list(V = 1, nu = 0),
                            G1 =list(V = 1, nu = 0)))

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
# 520 - 49 cattle


a = 1000
prior_binom2 = list(R = list(V = 1, fix = 1),
                    B = list(mu = c(rep(0,5)),
                             V = diag(5)*(4.3+pi^2/3)),
                    G = list(G1 =list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                             G1 =list(V = diag(1), nu = 1,  alpha.mu = 0, alpha.V = diag(1)*a),
                             G1 =list(V = diag(1), nu = 1,  alpha.mu = 0, alpha.V = diag(1)*a)))

set.seed(123)
modt2 = MCMCglmm(fixed = rd_2_qpcr_present ~ log(BM_KG)+ GS+ GUT +UNDERSTORY_PROP,
                 random = ~phylo + MSW93_Binomial+Period,
                 ginverse = list(phylo = inv.phylo$Ainv),
                 data = MPprevinfoMCMC,
                 family = "categorical",
                 prior = prior_binom2,
                 nitt = 500000, burnin = 10000, thin = 250,
                 verbose = F, pl = T)


Vdam.2 = diag(colMeans(modt2$VCV)[1:3])
colnames(Vdam.2) = colnames(modt2$VCV)[1:3]
Vdam.2


# Calculate lambda
lambda_p = modt2$VCV[,'phylo']/
  (modt2$VCV[,'phylo']+modt2$VCV[,'units']+modt2$VCV[,'MSW93_Binomial']+modt2$VCV[,"Period"])
plot(density(lambda_p))
posterior.mode(lambda_p,0.25)
mean(lambda_p)

MCMCprevalence = as.data.frame(summary(modt2)$solutions)
MCMCprevalence_var = as.data.frame(summary(modt2)$Gcovariances)
lambda_prev = c(mean = mean(lambda_p), mode=posterior.mode(lambda_p, 0.25))

dim(MPprevinfoMCMC)
#471 samples


### Prevalence -glmmTMB ###
model_prevalence = glmmTMB(rd_2_qpcr_present ~ log(BM_KG)+ GS+ GUT +UNDERSTORY_PROP+ (1|MSW93_Binomial)+(1|Period), family="binomial", data=MPprevinfoMCMC)

# Save values
TMBprevalence = as.data.frame(summary(model_prevalence)$coefficients$cond)
TMBprevalence_var = unlist(summary(model_prevalence)$varcor$cond)


MPdata_prev2 = left_join(tab_prev,unique(dplyr::select(MPrichinfoMCMC,Species,GUT,GS)))

MPdata_prev2[which(MPdata_prev2$Species=="Cattle"),]$GUT="FG"

# Plot gut type
prevplotgp = emmeans(model_prevalence, ~GUT, type="response") %>% 
  as.data.frame() %>% 
  ggplot(aes(x=GUT, y=prob))+
  geom_errorbar(data=MPdata_prev2, aes(x=GUT, y=mean, ymin=lower, ymax=upper, col=Species), position=position_dodge(width=0.5), width=0.5)+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL ),size=1)+
  geom_point(size=1.5) +
  scale_color_manual(values=animal_colors)+
  theme_bw(base_size=14)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  guides(col="none")+
  labs(x="Gut Type", y="Prevalence")


preds = predict(model_prevalence, type="response", se=T)
MPprevinfoMCMC$pred = preds$fit
MPprevinfoMCMC$pred_se = preds$se.fit


# plot group size
ggplot(MPprevinfoMCMC, aes(x=GS, y=pred))+
  geom_point(data=MPdata_prev2, aes(x=GS, y=mean, col=Species), size=1.5)+
  geom_errorbar(data=MPdata_prev2, aes(x=GS, y=mean, ymin=lower, ymax=upper, col=Species), position = "dodge", width=0.5, alpha=0.8)+
  geom_smooth(method="glm", method.args=list(family="binomial"), col="gray20", se=F)+
  scale_color_manual(values=animal_colors)+
  facet_wrap(~GUT)+
  theme_bw(base_size=14)+guides(col=F)

gpplots = gridExtra::grid.arrange(prevplotgp, rich_sp_plot, ncol=2)

ggsave(here(paste("plots/2_sp_prev_rich",threshold_used,".png", sep="")), gpplots, device="png", dpi=300, width=10, height=5, units="in")



# write out information
result_2_3 = rbind(result_2, result_3)

write.table(result_1, here(paste("docs/2_rich_prev_correlation",threshold_used,".txt", sep="")), sep="\t")
write.table(result_2_3, here(paste("docs/2_anova_tests",threshold_used,".txt", sep="")), sep="\t")


MCMCprevalence$predictor = row.names(MCMCprevalence)
MCMCrichness$predictor = row.names(MCMCrichness)
TMBprevalence$predictor = row.names(TMBprevalence)
TMBrichness$predictor = row.names(TMBrichness)
MCMCprevalence_var$predictor = row.names(MCMCprevalence_var)
MCMCrichness_var$predictor = row.names(MCMCrichness_var)

MCMCmods = rbind(MCMCprevalence, MCMCrichness)
MCMCmods$Metric = c(rep("Prevalence",5), rep("Richness", 5))
MCMCvars = rbind(MCMCprevalence_var, MCMCrichness_var)
MCMCvars$Metric = c(rep("Prevalence",3), rep("Richness", 3))
TMBmods = rbind(TMBprevalence, TMBrichness)
TMBmods$Metric = c(rep("Prevalence",5), rep("Richness", 5))
TMBvars = cbind(TMBprevalence_var, TMBrichness_var)
names(TMBvars) = c("Prevalence","Richness")

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
names(allMods)=c("Predictor", "MCMC Prevalence", "eff.samp.prev", "Metric", "P2", "TMB Prevalence", "Metric", "P3",
                 "MCMC Richness", "eff.samp.rich", "Metric", "P4", "TMB Richness", "Metric")
allMods = allMods %>% 
  select(Predictor, `MCMC Prevalence`, eff.samp.prev, `TMB Prevalence`, `MCMC Richness`, eff.samp.rich,`TMB Richness`)

MCMCvars2
TMBvars

MCMCvars2 = data.frame(MCMCvars2[1:6,c(1:3)], TMB_Prev = c("","","",TMBvars[1,1],"",TMBvars[2,1]),
           MCMCvars2[7:12,2:3], TMB_Rich =  c("","","",TMBvars[1,1],"",TMBvars[2,1]))
MCMCvars2 = MCMCvars2 %>% 
  mutate_at(vars(TMB_Prev, TMB_Rich), funs(round(as.numeric(.),2))) %>% 
  mutate_at(vars(TMB_Prev, TMB_Rich), funs(ifelse(is.na(.), "",.)))
names(MCMCvars2) = c("Predictor", "MCMC Prevalence", "eff.samp.prev", "TMB Prevalence", "MCMC Richness", "eff.samp.rich",
                     "TMB Richness")            

full_results_table = rbind(allMods, MCMCvars2)

# cleaned results table
write.csv(full_results_table, here(paste("docs/2_prev_rich_model_results",threshold_used,".csv", sep="")), row.names = F)
# lambda
write.table(as.data.frame(rbind(lambda_prev,lambda_rich)), here(paste("docs/2_prev_rich_lambda",threshold_used,".txt")), sep="\t")


# species summary
prev_result = tab %>% 
  group_by(Species) %>% 
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
