setwd("C:/Users/Elen/Desktop/Regressão")


require(pROC)
require(car)
require(lmtest)
require(effects)
require(statmod)
require(labestData)
require(hnp)
require(coefplot)
require(latticeExtra)
require(MASS)
require(glmnet)
require(ROCR)


### Selecionando as covariaveis para analise.
dados <- read.csv(file = "indian_liver_patient_dataset.csv",header = T)

# Tem que tratar os dados porque algum ser das trevas colocou
# 1 para doença e 2 para não doença
n <- length(dados$class)
for (i in 1:n) {
  if (dados$class[i]==2) {
    dados$class[i] <- 0
  }
}

## TRANSFORMANDO A VARIÁVEL RESPOSTA EM FACTOR
dados$class <- as.factor(dados$class)
dados$gender <- as.factor(dados$gender)
dados$TB <- log(dados$TB)
dados$DB <- log(dados$DB)
dados$alkphos <- log(dados$alkphos)
dados$sgpt <- log(dados$sgpt)
dados$sgot <- log(dados$sgot)
dados$A_G <- log(dados$A_G)

summary(dados)

x11()
par(mfrow = c(4,3))
barplot(table(dados$class), xlab = "CLASS")
hist(dados$age, xlab = "AGE", main = "")
barplot(table(dados$gender), xlab = "GENDER")
hist(dados$TB, xlab = "log(TB)", main = "")
hist(dados$DB, xlab = "log(DB)", main = "")
hist(dados$alkphos, xlab = "log(ALKPHOS)", main = "")
hist(dados$sgpt, xlab = "log(SGPT)", main = "")
hist(dados$sgot, xlab = "log(SGOT)", main = "")
hist(dados$TP, xlab = "TP", main = "")
hist(dados$ALB, xlab = "log(ALB)", main = "")
hist(dados$A_G, xlab = "Log(AG)", main = "")
par(mfrow = c(1,1))



### DEFINIR A QUANTIDADE DE DADOS PARA AJUSTE E VALIDACAO ###
set.seed(666)
indices <- sample(1:583, size = 408) 
### indices de um vetor com numeros de 1 a 583 numa sequencia aleatória.

dadosajuste <- dados[indices,]
### dataframe com as 583 linhas para ajuste (70%).

dadosvalid <- dados[-indices,]
### dataframe com 175 linhas, apenas para validacao (30%).

## REGRESSÃO VIA STEPWISE LOGITO
reglog <- step(glm(class ~ age + TB + sgpt + TP + A_G,
              family = binomial, data=dadosajuste),
               direction = 'backward')
summary(reglog)

reglog1 <- update(reglog, subset = -c(116,272,476,573,576))
compareCoefs(reglog, reglog1)

## REGRESSÃO VIA STEPWISE PROBITO
regprob <- step(glm(class ~ age + TB + sgpt + TP + A_G,
                family = binomial(link = "probit"), data=dadosajuste),
                direction = 'both')
summary(regprob)

regprob1 <- update(regprob, subset = -c(116,272,476,573,576))
compareCoefs(regprob, regprob1)

## REGRESSÃO VIA STEPWISE CLOGLOG
regclog <- step(glm(class ~ age + TB + sgpt + TP + A_G,
                family = binomial(link = "cloglog"), data=dadosajuste),
                direction = 'both')
summary(regclog)

regclog1 <- update(regclog, subset = -c(116,272,476,573,576))
compareCoefs(regclog, regclog1)

AIC(reglog, regprob, regclog)
BIC(reglog, regprob, regclog)

## ANÁLISE DE RESÍDUOS PARA TODOS OS MODELOS
set.seed(666)
residuoslog <- qresid(reglog)
ajustadoslog <- predict(reglog)
residuosprob <- qresid(regprob)
ajustadosprob <- predict(regprob)
residuosclog <- qresid(regclog)
ajustadosclog <- predict(regclog)

x11()
par(mfrow = c(3,2))
plot(residuoslog ~ ajustadoslog, pch = 20, cex = 1.4, col = 'blue', 
     main = "Link LOGITO")
qqnorm(residuoslog, pch = 20, cex = 1.4, col = 'blue')
qqline(residuoslog)
plot(residuosprob ~ ajustadosprob, pch = 20, cex = 1.4, col = 'blue',
     main = "Link PROBITO")
qqnorm(residuosprob, pch = 20, cex = 1.4, col = 'blue')
qqline(residuosprob)
plot(residuosclog ~ ajustadosclog, pch = 20, cex = 1.4, col = 'blue',
     main = "Link CLOGLOG")
qqnorm(residuosclog, pch = 20, cex = 1.4, col = 'blue')
qqline(residuosclog)
par(mfrow = c(1,1))

shapiro.test(residuoslog) # 0.2939
shapiro.test(residuosprob)# 0.9309
shapiro.test(residuosclog)# 0.5005

par(mfrow = c(1,3))
hist(residuoslog, main = "Link LOGITO")
hist(residuosprob, main = "Link PROBITO")
hist(residuosclog, main = "Link CLOGLOG")


## MATRIZ DE COVARIÂNcIA
vcov(reglog)
vcov(regprob)
vcov(regclog)


## TESTE VIF PARA MULTICOLINEARIDADE
round(vif(reglog),digits = 3)
round(vif(regprob),digits = 3)
round(vif(regclog),digits = 3)

## ENVELOPES SIMULADOS ESTÃO LINDÕES 
hnp(reglog, pch = 20, cex = 1.2, main = "Link LOGITO")
hnp(regprob, pch = 20, cex = 1.2, main = "Link PROBITO")
hnp(regclog, pch = 20, cex = 1.2, main = "Link CLOGLOG")
par(mfrow = c(1,1))

## GRÁFICOS PARA PONTOS CAGADOS
influenceIndexPlot(reglog, vars = c('Studentized','Cook','Hat'), 
                   id.n = 3, cex = 1.4)
influenceIndexPlot(regprob, vars = c('Studentized','Cook','Hat'), 
                   id.n = 3, cex = 1.4)
influenceIndexPlot(regclog, vars = c('Studentized','Cook','Hat'), 
                   id.n = 3, cex = 1.4)

## TUDO AS ANOVAS
anova(reglog, regprob, regclog, test='Chisq')

pajlog <- predict(reglog, type='response', newdata=dadosajuste)
pajprob <- predict(regprob, type='response', newdata=dadosajuste)
pajclog <- predict(regclog, type='response', newdata=dadosajuste)
### Probabilidades estimadas para os indivíduos da base de ajuste;

pvallog <-predict(reglog, type='response', newdata=dadosvalid) 
pvalprob <-predict(regprob, type='response', newdata=dadosvalid) 
pvalclog <-predict(regclog, type='response', newdata=dadosvalid) 
### Probabilidades estimadas para os indivíduos da base de validação.


################################################################################
### Usando o modelo ajustado para classificação da probabilidade de doença
### da base de validação.

predlog <- prediction(pvallog, dadosvalid$class)
predprob <- prediction(pvalprob, dadosvalid$class)
predclog <- prediction(pvalclog, dadosvalid$class)

### Vamos plotar a curva ROC
perflog <- performance(predlog, measure = "tpr" , x.measure = "fpr")
perfprob <- performance(predprob, measure = "tpr" , x.measure = "fpr")
perfclog <- performance(predclog, measure = "tpr" , x.measure = "fpr")
# tpr: True Positive Rate; fpr: False Positive Rate.

## CURVA ROC 
par(mfrow = c(1,3))
plot(perflog, colorize=TRUE, lwd = 2, main = "Link LOGITO")
abline(0,1, lty = 2)
plot(perfprob, colorize=TRUE, lwd = 2, main = "Link PROBITO")
abline(0,1, lty = 2)
plot(perfclog, colorize=TRUE, lwd = 2, main = "Link CLOGLOG")
abline(0,1, lty = 2)

## ÁREA EMBAIXO DA CURVA ROC
auclog <- performance(predlog, 'auc')
aucprob <- performance(predprob, 'auc')
aucclog <- performance(predclog, 'auc')

# CURVA ROC COM DADOS DE AJUSTE
predajlog <- prediction(pajlog, dadosajuste$class)
perfajlog <- performance(predajlog, measure = "tpr" , x.measure = "fpr") 
plot(perfajlog, col = 'red', lwd = 2, add = T)
abline(0,1, lty = 2)

predajprob <- prediction(pajprob, dadosajuste$class)
perfajprob <- performance(predajprob, measure = "tpr" , x.measure = "fpr") 
plot(perfajprob, col = 'red', lwd = 2, add = T)
abline(0,1, lty = 2)

predajclog <- prediction(pajclog, dadosajuste$class)
perfajclog <- performance(predajclog, measure = "tpr" , x.measure = "fpr") 
plot(perfajclog, col = 'red', lwd = 2, add = T)
abline(0,1, lty = 2)
par(mfrow = c(1,1))

# ÁREA EMBAIXO DA CURVA AJUSTADA
performance(predajlog, 'auc')
performance(predajprob, 'auc')
performance(predajclog, 'auc')

# KS DO MODELO AJUSTADO
# A estatística KS nos diz o quanto o modelo separa os pacientes 
# doentes dos não doentes. Entre 30% e 60% está bom demais.
kslog <- max(attr(perfajlog, "y.values")[[1]] - (attr(perfajlog, "x.values")[[1]]))
ksprob <- max(attr(perfajprob, "y.values")[[1]] - (attr(perfajprob, "x.values")[[1]]))
ksclog <- max(attr(perfajclog, "y.values")[[1]] - (attr(perfajclog, "x.values")[[1]]))

## EM NOSSA BASE TEMOS 71% DE DOENTES E 29% DE NÃO DOENTES

predroc1 <- roc(dadosvalid$class, pvallog, percent = TRUE)
predroc1$auc
plot(predroc1, print.thres = seq(0.1,0.95,0.05), print.thres.pattern.cex = 0.8)
c1 <- coords(predroc1, 'best'); c1

predroc2 <- roc(dadosvalid$class, pvalprob, percent = TRUE)
predroc2$auc
plot(predroc2, print.thres = seq(0.1,0.95,0.05), print.thres.pattern.cex = 0.8)
c2 <- coords(predroc2, 'best'); c2

predroc3 <- roc(dadosvalid$class, pvalclog, percent = TRUE)
predroc3$auc
plot(predroc3, print.thres = seq(0.1,0.95,0.05), print.thres.pattern.cex = 0.8)
c3 <- coords(predroc3, 'best'); c3

# TABELA DE CONTINGÊNCIA 
table(dadosvalid$class, pvallog > 0.6926)
table(dadosvalid$class, pvalprob > 0.6909)
table(dadosvalid$class, pvalclog > 0.6832)

total <- 175
a <- 83  # VERDADEIROS POSITIVOS
d <- 42  # VERDADEIROS NEGATIVOS
c <- 11 # FALSOS POSITIVOS
b <- 39  # FALSOS NEGATIVOS

sensibilidade <- a/(a+b)
especificidade <- d/(c+d)
acuracia <- (a+d)/total

reg <- data.frame(age = 44, TB = 0.4,
                       sgpt = 3, A_G = 0)
### Criando um data frame para os dados que vamos predizer.

final <- glm(class ~ age + TB + sgpt + A_G,
    family = binomial, data=dadosajuste)

predict(final, newdata=reg, type='link') 

?predict

-3.105+(0.016*44)+(0.556*0.40)+(0.882*3)-0.685*(-0.02)
