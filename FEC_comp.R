nem <- read.csv("FEC_comparison.csv")

head(nem)

library(tidyverse)

ggplot(nem, aes(x=FEC, y=Reads))+
  geom_point(aes(color=Species))

cor(nem$Reads,nem$FEC)             

hist(nem$Reads)             
hist(log(nem$Reads+1))
hist(nem$FEC)
hist(log(nem$FEC+1))

nem$lRead <- log(nem$Reads + 1)
nem$lFEC <- log(nem$FEC +1)

ggplot(nem, aes(x=FEC, y=lRead))+
  geom_point(aes(color=Species))+
  geom_vline(xintercept=80)+
  geom_hline(yintercept=7.95)

head(nem)
p <- mutate(nem, ifelse(nem$FEC >0,1,0))
q <- mutate(nem, ifelse(nem$Reads >0,1,0))
nem$Rbin <- q[,16]
nem$Fbin <- p[,16]

library(ROCR)
preds <- prediction(predictions=nem$Rbin, labels=nem$Fbin)
perf <- performance(preds, "tpr","fpr")
plot(perf)

nem

mod <- glm(Fbin~lRead, data=nem, family="binomial")
summary(mod)

ps <- predict(mod)

preds <- prediction(predictions=ps, labels=nem$Fbin)
perf <- performance(preds, "tpr","fpr")
plot(perf, colorize=T)
abline(b=1, a=0)
# best threshold is where there are no false positives

cor.test(nem$lFEC, nem$lRead, method="spearman")
