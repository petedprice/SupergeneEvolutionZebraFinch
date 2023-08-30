#Dn_cont, Ds_cont, and Dtotal_cont are D_nonsyn, D_syn, D_total controlled for total exon length
#Reads are measured as log rpkm, to account for library size differences between samples and control for gene length
library(lmerTest)
library(lme4)
library(ggplot2)
#I've shortcut the model selection here by using lmerTest as it outputs 
#p values for variables from lmer

data <- read.table("scripts/work/Jon_send/JON_DEGDATA.txt")

hist(data$reads)
ggplot(data, aes(x = reads, y = D_nonsyn)) + geom_point()
#Data appears relatively normally distributed so using default gaussian dist
#data <- filter(data, Ds_cont > 0) #If included (and removing models with dif_or_not it tends to improve diagnostics and panterns remain the same)
data$logDn_cont <- log(data$Dn_cont)
data$logDn_cont[is.infinite(data$logDn_cont)] <- NA
data$logDs_cont <- log(data$Ds_cont)
data$logDs_cont[is.infinite(data$logDs_cont)] <- NA
data$logDtotal_cont <- log(data$Dtotal_cont)
data$logDtotal_cont[is.infinite(data$logDtotal_cont)] <- NA

(lm(reads ~ Dn_cont, data) %>% summary())$r.squared
(lm(reads ~ Dn_cont * genotype * gene, data) %>% summary())$r.squared

(lm(reads ~ exon_length, data) %>% summary())$r.squared


lm(reads ~ Dn_cont * genotype, data) %>% summary()
lm(reads ~ Ds_cont * genotype, filter(data, reads < 8)) %>% summary()
lm(reads ~ Dtotal_cont* genotype, filter(data, reads < 8)) %>% summary()

ggplot(data, aes(x = Dn_cont, y = reads)) + geom_point() + 
  geom_smooth(method = 'lm')
ggplot(data, aes(x = Ds_cont, y = reads)) + geom_point() + 
  geom_smooth(method = 'lm')
ggplot(data, aes(x = Dtotal_cont, y = reads)) + geom_point() + 
  geom_smooth(method = 'lm')

lm(reads ~ genotype * Ds_cont , data) %>% summary()
lm(reads ~ genotype * Ds_cont + genotype:Ds_cont, data) %>% summary()

model1 <- lm(reads ~ genotype * Ds_cont * gene, data)

model1_sum <- summary(model1)
model1_sum$coefficients[rownames(model1_sum$coefficients) 
                        == "genotypeAB:Ds_cont"]
model1_sum$coefficients[rownames(model1_sum$coefficients) 
                        == "genotypeAB"]
model1_sum$coefficients[rownames(model1_sum$coefficients) 
                        == "Ds_cont"]



model2 <- lm(reads ~ genotype * Dn_cont * gene, data)
model2_sum <- summary(model2)
model2_sum$coefficients[rownames(model2_sum$coefficients) 
                        == "genotypeAB:Dn_cont"]
model2_sum$coefficients[rownames(model2_sum$coefficients) 
                        == "genotypeAB"]
model2_sum$coefficients[rownames(model2_sum$coefficients) 
                        == "Dn_cont"]

model3 <- lm(reads ~ genotype * Dtotal_cont * gene, data)
model3_sum <- summary(model3)
model3_sum$coefficients[rownames(model3_sum$coefficients) 
                        == "genotypeAB:Dtotal_cont"]
model3_sum$coefficients[rownames(model3_sum$coefficients) 
                        == "genotypeAB"]
model3_sum$coefficients[rownames(model3_sum$coefficients) 
                        == "Dtotal_cont"]




lm(reads ~ Dtotal_cont, filter(data, reads < 8)) %>% summary()
lm(reads ~ gene, data) %>% summary()


ggplot(data, aes(x = logDtotal_cont, y = reads)) + geom_point() + 
   geom_smooth(method = 'lm')

model <- lmer(reads ~ (log(Dn_cont) * genotype * gene) + (1|sample), data)
summary(model)

model1 <- lmer(reads ~ Ds_cont + (1|sample) + (1|gene), data)
model2 <- lmer(reads ~ Dn_cont  + (1|sample) + (1|gene), data)
model3 <- lmer(reads ~ genotype + (1|sample) + (1|gene), data)
model4 <- lmer(reads ~ (Ds_cont * genotype * gene)+ (1|sample), data)
model5 <- lmer(reads ~ (Dn_cont * genotype)+ (1|sample) + (1|gene), data)

model6 <- lmer(reads ~ Ds_cont + dif_or_not + (1|sample) + (1|gene), data)
model7 <- lmer(reads ~ Dn_cont + dif_or_not + (1|sample) + (1|gene), data)
model8 <- lmer(reads ~ (Ds_cont * genotype)+ dif_or_not+ (1|sample) + (1|gene), data)
model9 <- lmer(reads ~ (Dn_cont * genotype)+ dif_or_not+ (1|sample) + (1|gene), data)

anova(model1, model2, model3, model4, model5)

anova(model1, model2, model3, model4, model5, model6, model7, model8, model9)

#FOR DEG data, it appears that model 6, 4, or 8 are the best fits. 
model4_final <- lmerTest::lmer(reads ~ (Ds_cont * genotype)+
                                 (1|sample) + (1|gene), data)
model6_final <- lmerTest::lmer(reads ~ Ds_cont + dif_or_not + 
                                 (1|sample) + (1|gene), data)
model8_final <- lmerTest::lmer(reads ~ (Ds_cont * genotype)+ dif_or_not+ 
                                 (1|sample) + (1|gene), data)

model4_final %>% summary()
model6_final %>% summary()
model8_final %>% summary()

simulationOutput4 <- simulateResiduals(model4_final)
plot(simulationOutput4)
simulationOutput6 <- simulateResiduals(model6_final)
plot(simulationOutput6)
simulationOutput8 <- simulateResiduals(model8_final)
plot(simulationOutput8)

############## ASE 
data2 <- read.table("JON_ASEDATA.txt")
#data2 <- data2 %>% filter(Dtotal > 0)
data2$gene <- data2$name

ASEmodel1 <- glmer(cbind(aCount, totalCount) ~ Ds_cont + 
                     (1|sample) + (1|gene), data2, family = "binomial")
ASEmodel2 <- glmer(cbind(aCount, totalCount) ~ Dn_cont  + 
                     (1|sample) + (1|gene), data2, family = "binomial")
ASEmodel3 <- glmer(cbind(aCount, totalCount) ~ genotype + 
                     (1|sample) + (1|gene), data2, family = "binomial")
ASEmodel4 <- glmer(cbind(aCount, totalCount) ~ (Ds_cont * genotype)+ 
                     (1|sample) + (1|gene), data2, family = "binomial")
ASEmodel5 <- glmer(cbind(aCount, totalCount) ~ (Dn_cont * genotype)+ 
                     (1|sample) + (1|gene), data2, family = "binomial")

ASEmodel6 <- glmer(cbind(aCount, totalCount) ~ Ds_cont + dif_or_not + 
                     (1|sample) + (1|gene), data2, family = "binomial")
ASEmodel7 <- glmer(cbind(aCount, totalCount) ~ Dn_cont + dif_or_not + 
                     (1|sample) + (1|gene), data2, family = "binomial")
ASEmodel8 <- glmer(cbind(aCount, totalCount) ~ (Ds_cont * genotype)+ dif_or_not+ 
                     (1|sample) + (1|gene), data2, family = "binomial")
ASEmodel9 <- glmer(cbind(aCount, totalCount) ~ (Dn_cont * genotype)+ dif_or_not+ 
                     (1|sample) + (1|gene), data2, family = "binomial")


anova(ASEmodel1, ASEmodel2, ASEmodel3, ASEmodel4, ASEmodel5, ASEmodel6, ASEmodel7, ASEmodel8, ASEmodel9)


#best fit is model4 so looking at variables
ASEmodel4 %>% summary()
#it appears that it's mostly genotype that contributes to ASE differences.

simulationOutput4 <- simulateResiduals(ASEmodel8)
plot(simulationOutput4)

