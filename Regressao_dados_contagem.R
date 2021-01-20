setwd("C:/Users/Elen/Desktop/Regressão")

require(corrplot)
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
require(gamlss)

dados <- read.csv2(file = "amazonas2.csv",header = T)

summary(dados)

corrplot(cor(dados), type = "lower")

# Análise Exploratória

x11()
par(mfrow = c(3,5))
boxplot(dados$obitos, xlab = "OBITOS")         # Total de óbitos por  municípios
boxplot(dados$dens, xlab = "DENS")             # % da população em domicílios com densidade > 2
boxplot(log(dados$total), xlab = "log(POP)")   # População total
boxplot(dados$fec, xlab = "FEC")               # Taxa de fecundidade total
boxplot(dados$est, xlab = "EST")               # Expectativa de anos de estudo
boxplot(log(dados$renda), xlab = "log(RENDA)") # Renda per capita
boxplot(dados$anal, xlab = "ANALF")            # Taxa de Abalfabetismo
boxplot(log(dados$pib), xlab = "log(PIB)")     # PIB per Capita
boxplot(dados$des, xlab = "DES")               # Taxa de desemprego
boxplot(dados$gini, xlab = "GINI")             # Índice de Gini
boxplot(dados$esp, xlab = "ESP")               # Esperança de vida ao nascer
boxplot(dados$mort, xlab = "MORT")             # Mortalidade infantil
boxplot(dados$dep, xlab = "DEP")               # Razão de dependência
boxplot(dados$env, xlab = "ENV")               # Taxa de envelhecimento
boxplot(dados$pob, xlab = "POB")               # % de pobres 2010
par(mfrow = c(1,1))

### Modelo 1: USANDO GLM.NB Binomial Negativa

reg1 <- step(glm.nb(obitos ~ dens + fec + est + anal + des +
                      gini + esp + mort + dep + env + pob +
                      offset(log(total)), data = dados),direction = 'back')
summary(reg1)
phi1 <- sum(rstandard(reg1, type='pearson')**2)/57
summary(reg1, dispersion = phi1)
par(mfrow = c(2,2))
plot(reg1)
par(mfrow = c(1,1))
vif(reg1)

#Cãoferindo os resíduos - Tão bunitinhos
set.seed(23)
residuos <- qresid(reg1)
ajustados <- predict(reg1)

par(mfrow = c(1,2))
plot(residuos ~ ajustados, pch = 20, cex = 1.4, col = 'blue')
qqnorm(residuos, pch = 20, cex = 1.4, col = 'blue')
qqline(residuos)

#Shapirão okay
shapiro.test(residuos)

#Checando pontos influentes
influenceIndexPlot(reg1, vars = c('Studentized','Cook','Hat'), 
                   id.n = 3, cex = 1.4)

#Modelando sem o ponto possivelmente infuente
reg11 <- update(reg1, subset = -c(20,53,38))
summary(reg11)
anova(reg1, reg11)
compareCoefs(reg1, reg11)

hnp(reg1)

### Modelo 6: USANDO GLM QuasiPoisson

reg6 <- glm(obitos ~ anal + dep + 
              offset(log(total)), 
                 family = quasipoisson(link = "log"), data = dados)
summary(reg6)
phi2 <- sum(rstandard(reg6, type='pearson')**2)/58
summary(reg6, dispersion = phi2)
par(mfrow = c(2,2)) 
plot(reg6)
vif(reg6)
par(mfrow = c(1,1))

#Cãoferindo os resíduos 
set.seed(23)
residuos6 <- qresid(reg6)
ajustados6 <- predict(reg6)

par(mfrow = c(1,2))
plot(residuos6 ~ ajustados6, pch = 20, cex = 1.4, col = 'blue')
qqnorm(residuos6, pch = 20, cex = 1.4, col = 'blue')
qqline(residuos6)

#Shapirão
shapiro.test(residuos6)

#Checando pontos influentes
influenceIndexPlot(reg6, vars = c('Studentized','Cook','Hat'), 
                   id.n = 3, cex = 1.4)

#Modelando sem o ponto possivelmente infuente
reg66 <- update(reg6, subset = -c(20,8))
summary(reg66, dispersion = phi2)
compareCoefs(reg66, reg6)
anova(reg6, reg66)
vif(reg66)
hnp(reg6)



### Gráficos de efeitos para a regressão Binomial Negativa
plot(allEffects(reg1), type = 'response')

x11()
bootreg1 <- Boot(reg1, R = 999)
plot(bootreg1)

bootreg6 <- Boot(reg6, R = 999)
plot(bootreg6)


par(mfrow = c(1,2))
hnp(reg1,main = "Poisson-Gama")
hnp(reg6, main = "quase-Poisson")
par(mfrow = c(1,1))

