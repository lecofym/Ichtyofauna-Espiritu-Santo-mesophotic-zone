#############################################################################################################################
#Manuscript Number: DSR1-D-23-00116 
#The mesophotic zone of the Marine Protected Area Espiritu Santo can function as a refuge for the ichthyofauna 

#Script to analysis
#use data from supplementary material
#############################################################################################################################

##############################------PERMANOVA-----##############################

library(vegan)
library(dplyr)
library(devtools)
library(pairwiseAdonis)

#--------------
#---Richness
data <- read.csv("richness.csv",header=TRUE,row.names=1,stringsAsFactors = T)
richness.dist <- vegdist(richness, method='jaccard')

#dispersion tests
dispersion.site=betadisper(richness.dist, group=data$Site)
permutest(dispersion.site,permutations = how(nperm = 10000))
plot(dispersion.site, hull=FALSE, ellipse=TRUE) ##sd ellipse

dispersion.zone=betadisper(richness.dist, group=date$zone)
permutest(dispersion.zone,permutations = how(nperm = 10000))
plot(dispersion.zone, hull=FALSE, ellipse=TRUE) ##sd ellipse

#permanova
richness.site.zone <-adonis2(richness.dist~Site*Zone, data=data, permutations = 10000, method="jaccard")
richness.site.zone
pairwise.adonis(richness.dist, data$Site, sim.method = 'jaccard', p.adjust.m = 'bonferroni',perm =10000)
pairwise.adonis(richness.dist, data$Zone, sim.method = 'jaccard', p.adjust.m = 'bonferroni',perm =10000)
densityplot(permustats(richness.site.zone))

#--------------
#---Abundance
data <- read.csv("abundance.csv",header=TRUE,row.names=1,stringsAsFactors = T)
abundance.dist <- vegdist(abundance, method='bray')

#dispersion tests
dispersion.site=betadisper(abundance.dist, group=data$Site)
permutest(dispersion.site,permutations = how(nperm = 10000))
plot(dispersion.site, hull=FALSE, ellipse=TRUE) ##sd ellipse

dispersion.zone=betadisper(abundance.dist, group=data$Zone)
permutest(dispersion.zone,permutations = how(nperm = 10000))
plot(dispersion.zone, hull=FALSE, ellipse=TRUE) ##sd ellipse

#permanova
abundancia.site.zone <-adonis2(abundance.dist~Site*Zone, data=data, permutations = 10000, method="bray")
abundancia.sitio.estrato
pairwise.adonis(abundance.dist, data$Site, sim.method = 'bray', p.adjust.m = 'bonferroni', perm =10000)
pairwise.adonis(abundance.dist, data$Zone, sim.method = 'bray', p.adjust.m = 'bonferroni', perm =10000)
densityplot(permustats(abundancia.site.zone))

#--------------
#---Biomass
data <- read.csv("biomass.csv",header=TRUE,row.names=1,stringsAsFactors = T)
biomass.dist <- vegdist(biomass, method='bray')

#dispersion tests
dispersion.site=betadisper(biomasa.dist, group=data$Site)
permutest(dispersion.site,permutations = how(nperm = 10000))
plot(dispersion.site, hull=FALSE, ellipse=TRUE) ##sd ellipse

dispersion.zone=betadisper(biomass.dist, group=data$Zone)
permutest(dispersion.zone,permutations = how(nperm = 10000))
plot(dispersion.zone, hull=FALSE, ellipse=TRUE) ##sd ellipse

#permanova
biomass.site.zone<-adonis2(biomasa.dist~Site*Zone, data=data, permutations = 10000, method="bray")
biomass.site.zone
pairwise.adonis(biomass.dist, datos$Site, sim.method = 'bray', p.adjust.m = 'bonferroni', perm =10000)
pairwise.adonis(biomass.dist, datos$Zone, sim.method = 'bray', p.adjust.m = 'bonferroni', perm =10000)
densityplot(permustats(biomass.site.zone))


##############################------FUNCTIONAL INDICES-----##############################

library("lme4")
library("jtools")
library("vegan")
library("arm")
library("ggplot2")
library("ade4")
library("diveRsity")
library("cluster")
library("clue")
library("FD")
library("ggrepel")
library("reshape2")
library("dplyr")

#--------------
#'multidimFD': function to compute and illustrate multidimensional
# functional diversity indices for a set of species assemblages. For details 
# about indices formulas see Mouillot et al. 2013, Trends in ecology and 
# Evolution (28:167-177) and references therein.

traits <- read.csv("functional matrix.csv",header=TRUE,row.names=1,stringsAsFactors = T)
#---with species per transect or zone
abundance <- read.csv("abundance.csv",header=TRUE,row.names=1)
biomass <- read.csv("biomass.csv",header=TRUE,row.names=1)

#transform to ln(x+1)
abundance.1 <- log1p(abundance.1)
biomass.1 <- log1p(biomass.1)

#construct dissimilarity matrix with Gower distance coefficient
dissimilarity <- daisy(traits,"gower")

#run Principal coordinate analysis (PCoA)
PCoA <- pcoa(dissimilarity)

#analyze the quality of the functional space using the method de Maire et al. 2015
source("quality_funct_space.R")
source("plot_funct_space.R")
source("multidimFbetaD.R")

quality <- quality_funct_space(traits, nbdim = 6, metric = "Gower", dendro = FALSE, plot = "quality_funct_space")
quality$meanSD

#select the first 4 coordinates of the PCoA matrix
four.coordinates <- PCoA$vectors[,1:4]
coordinates <- as.data.frame(four.coordinates)

#Calculate indices with abundance
source("multidimFD.R")
abundance.1 <- as.matrix(abundance.1)
indices <- multidimFD(four.coordinates, abundance.1)

#Calculate indices with biomass
source("multidimFD.R")
biomass.1 <- as.matrix(biomass.1)
indices <- multidimFD(four.coordinates, biomass.1)


##############################------RLM-----##############################

library(gvlma)
library(olsrr)
library(MASS)
library(tseries)
library(lmtest)
library(car)

#--------------
data.all <- read.csv("TableA1data.csv",header=TRUE,row.names=1,stringsAsFactors = T)
#standardized data (Mean =0, SD =1)

#---Richness
S.RLM.full=lm(S~.,data.all)
summary(S.RLM.full)
S.RLM1 <- stepAIC(S.RLM.full, trace=TRUE, direction="backward")
summary(S.RLM1) #
S.RLM1$anova

S.RLM.empty <- lm(S ~ 1, data=data.all)
horizon <- formula(S ~ T+OD+Ilum)
S.RLM2 <- stepAIC(S.RLM.empty, trace=FALSE, direction="forward", scope=horizon)
summary(S.RLM2)#
S.RLM2=lm(S ~ T,data.all)
S.RLM2$anova
vif(S.RLM2)

S.RLM3 <- stepAIC(S.RLM.empty, trace=FALSE, direction="both", scope=horizon)
S.RLM3$anova

S.RLM3=lm(S..spp.~ T -1,data.all)
summary(S.RLM3)
anova(S.RLM3)
summary(S.RLM3)$adj.r.squared#0.7038303 (determination coefficient)
sqrt(0.7038303)# 0.8389459 (correlation coefficient) 

#residuals
mean(S.RLM3$residuals)#zero mean
shapiro.test(resid(S.RLM3))#normal
durbinWatsonTest(S.RLM3)#independence
bptest(S.RLM3)#homoscedasticity
outoliers=S.RLM3$residuals[S.RLM3$residuals<mean(S.RLM3$residuals) - 2 * sd(S.RLM3$residuals)]
outoliers #

par(mfrow = c(1,2))
plot(S.RLM3)

#---Abundance
N.RLM.full=lm(N~.,data.all)
summary(N.RLM.full)
N.RLM1 <- stepAIC(N.RLM.full, trace=TRUE, direction="backward")
summary(N.RLM1)#
N.RLM1$anova

N.RLM.empty <- lm(N ~ 1, data=data.all)
horizon <- formula(N ~ T+OD+Ilum)
N.RLM2 <- stepAIC(N.RLM.empty, trace=FALSE, direction="forward", scope=horizon)
summary(N.RLM2)#
N.RLM2$anova
vif(N.RLM2)

N.RLM3 <- stepAIC(N.RLM.empty, trace=FALSE, direction="both", scope=horizon)
summary(N.RLM3)#
N.RLM3$anova

N.RLM3=lm(N~Ilum -1,data.all)
summary(N.RLM3)
anova(N.RLM3)
summary(N.RLM3)$adj.r.squared#0.4152127 (determination coefficient)
sqrt(0.4152127)# 0.64437 (correlation coefficient) 

#residuales
mean(N.RLM3$residuals)#zero mean
shapiro.test(resid(N.RLM3))#normal
durbinWatsonTest(N.RLM3)#independence
bptest(N.RLM3)#homoscedasticity
outoliers=N.RLM3$residuals[N.RLM3$residuals<mean(N.RLM3$residuals) - 2 * sd(N.RLM3$residuals)]
outoliers #

par(mfrow = c(1,2))
plot(N.RLM3)

#---Biomass
B.RLM.full=lm(B~.,data.all)
summary(B.RLM.full)
B.RLM1 <- stepAIC(B.RLM.full, trace=TRUE, direction="backward")
summary(B.RLM1)#
B.RLM1$anova

B.RLM.empty <- lm(B ~ 1, data=data.all)
horizon <- formula(B ~ T+OD+Ilum)
B.RLM2 <- stepAIC(B.RLM.empty, trace=FALSE, direction="forward", scope=horizon)
summary(B.RLM2)#
B.RLM2$anova
vif(B.RLM2)

B.RLM3 <- stepAIC(B.RLM.empty, trace=FALSE, direction="both", scope=horizon)
summary(B.RLM3)#B ~ 1
B.RLM3$anova

#NO HAY MODELO IDEAL
B.RLM3=lm(B ~ 1,data.all)
summary(B.RLM3)
anova(B.RLM3)
summary(B.RLM3)$adj.r.squared#0 (determination coefficient)
sqrt(0)# 0 (correlation coefficient) 

#residuales
mean(B.RLM3$residuals)#zero mean
shapiro.test(resid(B.RLM3))#normal
durbinWatsonTest(B.RLM3)#independence
bptest(B.RLM3)#homoscedasticity
outoliers=B.RLM3$residuals[B.RLM3$residuals<mean(B.RLM3$residuals) - 2 * sd(B.RLM3$residuals)]
outoliers #

par(mfrow = c(1,2))
plot(B.RLM3)

#---FRic
FRic.RLM.full=lm(FRic~.,data.all)
summary(FRic.RLM.full)
FRic.RLM1 <- stepAIC(FRic.RLM.full, trace=TRUE, direction="backward")
summary(FRic.RLM1)#
FRic.RLM1$anova

FRic.RLM.empty <- lm(FRic ~ 1, data=data.all)
horizon <- formula(FRic ~ T+OD+Ilum)
FRic.RLM2 <- stepAIC(FRic.RLM.empty, trace=FALSE, direction="forward", scope=horizon)
summary(FRic.RLM2)#
FRic.RLM2$anova
vif(FRic.RLM2)

FRic.RLM3 <- stepAIC(FRic.RLM.empty, trace=FALSE, direction="both", scope=horizon)
summary(FRic.RLM3)#FRic ~ T
FRic.RLM3$anova

FRic.RLM3=lm(FRic ~ T -1,data.all)
summary(FRic.RLM3)
anova(FRic.RLM3)
summary(FRic.RLM3)$adj.r.squared#0.5351573 (determination coefficient)
sqrt(0.5351573)# O.7315445 (correlation coefficient) 

#residuales
mean(FRic.RLM3$residuals)#zero mean
shapiro.test(resid(FRic.RLM3))#normal
durbinWatsonTest(FRic.RLM3)#independence
bptest(FRic.RLM3)#homoscedasticity
outoliers=FRic.RLM3$residuals[FRic.RLM3$residuals<mean(FRic.RLM3$residuals) - 2 * sd(FRic.RLM3$residuals)]
outoliers #

par(mfrow = c(1,2))
plot(FRic.RLM3)

#---FEve
FEve.RLM.full=lm(FEve~.,data.all)
summary(FEve.RLM.full)
FEve.RLM1 <- stepAIC(FEve.RLM.full, trace=TRUE, direction="backward")
summary(FEve.RLM1)#FEve ~ 1
FEve.RLM$anova

FEve.RLM.empty <- lm(FEve ~ 1, data=data.all)
horizon <- formula(FEve ~ T+OD+Ilum)
FEve.RLM2 <- stepAIC(FEve.RLM.empty, trace=FALSE, direction="forward", scope=horizon)
summary(FEve.RLM2)#FEve ~ 1, data = datos5
FEve.RLM$anova
vif(FEve.RLM)

FEve.RLM3 <- stepAIC(FEve.RLM.empty, trace=FALSE, direction="both", scope=horizon)
summary(FEve.RLM3)#FEve ~ 1
FEve.RLM3$anova

#NO HAY MODELO IDEAL
FEve.RLM3=lm(FEve ~ 1,data.all)
summary(FEve.RLM3)
anova(FEve.RLM3)
summary(FEve.RLM3)$adj.r.squared#0 (determination coefficient)
sqrt(0)# 0 (correlation coefficient) 

#residuales
mean(FEve.RLM3$residuals)#zero mean
shapiro.test(resid(FEve.RLM3))#normal
durbinWatsonTest(FEve.RLM3)#independence
bptest(FEve.RLM3)#homoscedasticity
outoliers=FEve.RLM3$residuals[FEve.RLM3$residuals<mean(FEve.RLM3$residuals) - 2 * sd(FEve.RLM3$residuals)]
outoliers #

par(mfrow = c(1,2))
plot(FEve.RLM3)

#---FDiv
FDiv.RLM.full=lm(FDiv~.,data.all)
summary(FDiv.RLM.full)
FDiv.RLM1 <- stepAIC(FDiv.RLM.full, trace=TRUE, direction="backward")
summary(FDiv.RLM1)#
FDiv.RLM1$anova

FDiv.RLM.empty <- lm(FDiv ~ 1, data=data.all)
horizon <- formula(FDiv ~ T+OD+Ilum)
FDiv.RLM2 <- stepAIC(FDiv.RLM.empty, trace=FALSE, direction="forward", scope=horizon)
summary(FDiv.RLM2)#FDiv ~ 
FDiv.RLM2$anova
vif(FDiv.RLM2)

FDiv.RLM3 <- stepAIC(FDiv.RLM.empty, trace=FALSE, direction="both", scope=horizon)
summary(FDiv.RLM3)#FDiv ~ 1
FDiv.RLM3$anova

#NO HAY MODELO IDEAL
FDiv.RLM3=lm(FDiv ~ 1,data.all)
summary(FDiv.RLM3)
anova(FDiv.RLM3)
summary(FDiv.RLM3)$adj.r.squared#0 (determination coefficient)
sqrt(0)# 0 (correlation coefficient) 

#residuales
mean(FDiv.RLM3$residuals)#zero mean 
shapiro.test(resid(FDiv.RLM3))#normal
durbinWatsonTest(FDiv.RLM3)#independence
bptest(FDiv.RLM3)#homoscedasticity
outoliers=FDiv.RLM3$residuals[FDiv.RLM3$residuals<mean(FDiv.RLM3$residuals) - 2 * sd(FDiv.RLM3$residuals)]
outoliers #

par(mfrow = c(1,2))
plot(FDiv.RLM3)

#---FSpe
FSpe.RLM.full=lm(FSpe~.,data.all)
summary(FSpe.RLM.full)
FSpe.RLM1 <- stepAIC(FSpe.RLM.full, trace=TRUE, direction="backward")
summary(FSpe.RLM1)#FSpe 
FSpe.RLM1$anova

FSpe.RLM.empty <- lm(FSpe ~ 1, data=data.all)
horizon <- formula(FSpe ~ T+OD+Ilum)
FSpe.RLM2 <- stepAIC(FSpe.RLM.empty, trace=FALSE, direction="forward", scope=horizon)
summary(FSpe.RLM2)#FSpe ~
FSpe.RLM2$anova
vif(FSpe.RLM2)

FSpe.RLM3 <- stepAIC(FSpe.RLM.empty, trace=FALSE, direction="both", scope=horizon)
summary(FSpe.RLM3)#FSpe ~ OD
FSpe.RLM3$anova

FSpe.RLM3=lm(FSpe ~ OD,data.all)
summary(FSpe.RLM3)
anova(FSpe.RLM3)
summary(FSpe.RLM3)$adj.r.squared#0.07855123 (determination coefficient)
sqrt(0.07855123)# 0.2802699 (correlation coefficient) 

#residuales
mean(FSpe.RLM3$residuals)#zero mean
shapiro.test(resid(FSpe.RLM3))#normal
durbinWatsonTest(FSpe.RLM3)#independence
bptest(FSpe.RLM3)#homoscedasticity
outoliers=FSpe.RLM3$residuals[FSpe.RLM3$residuals<mean(FSpe.RLM3$residuals) - 2 * sd(FSpe.RLM3$residuals)]
outoliers #

par(mfrow = c(1,2))
plot(FSpe.RLM3)

#---FOri
FOri.RLM.full=lm(FOri~.,data.all)
summary(FOri.RLM.full)
FOri.RLM1 <- stepAIC(FOri.RLM.full, trace=TRUE, direction="backward")
summary(FOri.RLM1)#FOri ~ 
FOri.RLM1$anova

FOri.RLM.empty <- lm(FOri ~ 1, data=data.all)
horizon <- formula(FOri ~ T+OD+Ilum)
FOri.RLM2 <- stepAIC(FOri.RLM.empty, trace=FALSE, direction="forward", scope=horizon)
summary(FOri.RLM2)#FOri ~ 
FOri.RLM2$anova
vif(FOri.RLM2)

FOri.RLM3 <- stepAIC(FOri.RLM.empty, trace=FALSE, direction="both", scope=horizon)
summary(FOri.RLM3)#FOri ~ 1
FOri.RLM3$anova

#NO HAY MODELO IDEAL
FOri.RLM3=lm(FOri ~ 1,data.all)
summary(FOri.RLM3)
anova(FOri.RLM3)
summary(FOri.RLM3)$adj.r.squared# (determination coefficient)
sqrt()# (correlation coefficient) 

#residuales
mean(FOri.RLM3$residuals)#zero mean
shapiro.test(resid(FOri.RLM3))#normal
durbinWatsonTest(FOri.RLM3)#independence
bptest(FOri.RLM3)#homoscedasticity
outoliers=FOri.RLM3$residuals[FOri.RLM3$residuals<mean(FOri.RLM3$residuals) - 2 * sd(FOri.RLM3$residuals)]
outoliers #

par(mfrow = c(1,2))
plot(FOri.RLM3)

##############################------GLM-----##############################

library(MASS)
library(ggplot2)
library(dplyr)
library(car)

#--------------
data.all <- read.csv("TableA1data.csv",header=TRUE,row.names=1,stringsAsFactors = T)

#---Richness
S.GLM=glm(S ~. , family = poisson(link = "log"),data = data.all )
summary(S.GLM)
S.GLM2 <- stepAIC(S.GLM, trace=TRUE, direction="backward")
summary(S.GLM2)#AIC: 74.046
drop1(S.GLM2, test="Chi")
vif(S.GLM2)

#---Abundance
N.GLM=glm(N ~., family = poisson(link = "log"),data = data.all )
summary(N.GLM)
N.GLM2 <- stepAIC(N.GLM, trace=TRUE, direction="backward")
summary(N.GLM2)#AIC: 1114.9
N.GLM3=glm(N ~T+Ilum, family = poisson(link = "log") )
summary(N.GLM3)#AIC: 1115.3

drop1(N.GLM3, test="Chi")
vif(N.GLM3)

#---Biomass
B.GLM=glm(B ~., family = Gamma(link = "log"),data = data.all )
summary(B.GLM)

B.GLM2 <- stepAIC(B.GLM, trace=TRUE, direction="backward")

B.GLM2=glm(B ~T-1, family = Gamma(link = log) )
summary(B.GLM2)#AIC: 341.99

drop1(B.GLM2, test="Chi")
vif(B.GLM2)

#---FRic
FRic.GLM=glm(FRic ~., family = Gamma(link = "inverse"),data = data.all )
summary(FRic.GLM)
FRic.GLM2 <- stepAIC(FRic.GLM, trace=TRUE, direction="backward")
summary(FRic.GLM2)#

FRic.GLM3=glm(FRic ~T, family = Gamma(link = "inverse") )
summary(FRic.GLM3)#AIC: -30.012
drop1(FRic.GLM3, test="Chi")
vif(FRic.GLM3)

#---FEve
FEve.GLM=glm(FEve ~., family = Gamma(link = "inverse"),data = data.all )
summary(FEve.GLM)#

FEve.GLM2 <- stepAIC(FEve.GLM, trace=TRUE, direction="backward")
summary(FEve.GLM2)#AIC: -25.243

FEve.GLM=glm(FEve ~, family = Gamma(link = "inverse") )
drop1(FEve.GLM2, test="Chi")
vif(FEve.GLM2)

#---FDiv
FDiv.GLM=glm(FDiv ~., family = Gamma(link = "inverse"),data = data.all )
summary(FDiv.GLM)
FDiv.GLM2 <- stepAIC(FDiv.GLM, trace=TRUE, direction="backward")
summary(FDiv.GLM2)#-53.316

FDiv.GLM3=glm(FDiv ~, family = Gamma(link = "inverse") )
summary(FDiv.GLM3)#

drop1(FDiv.GLM2, test="Chi")
vif(FDiv.GLM2)

#---FSpe
FSpe.GLM=glm(FSpe ~., family = Gamma(link = "inverse"),data = data.all )
summary(FSpe.GLM)#
FSpe.GLM2 <- stepAIC(FSpe.GLM, trace=TRUE, direction="backward")
summary(FSpe.GLM2)#AIC: -60.855
drop1(FSpe.GLM2, test="Chi")
vif(FSpe.GLM2)

#---FOri
FOri.GLM=glm(FOri ~., family = Gamma(link = "inverse"),data = data.all )
summary(FOri.GLM)#
FOri.GLM2 <- stepAIC(FOri.GLM, trace=TRUE, direction="backward")
summary(FOri.GLM2)#AIC: -54.127
drop1(FOri.GLM2, test="Chi")
vif(FOri.GLM2)