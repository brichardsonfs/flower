require(lme4)
require(lmerTest)
require(ggplot2)
require(xlsx)
require(doBy)
require(MuMIn)
require(lattice)
require(grid)
require(gridExtra)

##Dataset is Table S2
fl11<- read.csv(file="TableS2.csv", sep=",",head=TRUE, na.string="na")

###Dataframe with 21 variables and 20 interactions
data_a=with(fl11, data.frame(
  lat=lat,
  tdiff=mtwm-mtcm,
  adi = (dd5**0.5)/map,
  adimindd0= ((dd5**0.5)/map)*mmindd0,
  d100,
  dd0,
  dd5,
  fday,
  ffp,
  gsdd5,
  gsp,
  dd5mtcm = dd5*mtcm,
  pratio = gsp/map,
  gspdd5 =(gsp*dd5)/1000,
  gspmtcm =(gsp*mtcm)/1000,
  gsptd  =(gsp*(mtwm-mtcm))/100,
  map=map,
  mapdd5 =(map*dd5)/1000,
  mapmtcm =(map*mtcm)/1000,
  maptd  =(map*(mtwm-mtcm))/100,
  mat=mat,
  mmindd0=mmindd0,
  mmax=mmax,
  mmin=mmin,
  mtcm=mtcm,
  mtcmgsp =mtcm/gsp,
  mtcmmap =mtcm/map,
  mtwm=mtwm,
  sday=sday,
  sdi=(gsdd5**0.5)/gsp,
  sdimindd0=((gsdd5**0.5)/gsp)*mmindd0,
  tdgsp  =(mtwm-mtcm)/gsp,
  tdmap  =(mtwm-mtcm)/map,
  smrpb,
  smrsprpb,
  sprp,
  winp,
  smrp,
  sdimtcm=((gsdd5**0.5)/gsp)*mtcm,
  dd0map=dd0/map,
  dd0gsp=dd0/gsp,
  julian))

###Correlation among variables
a_cor <- cor(data_a) 
write.xlsx(x = a_cor, file = "correl.xlsx")

###Table S2
Climate_vars <- cbind(fl11$pop, data_a)
colnames(Climate_vars)[1] <- "Population"
write.xlsx(x = Climate_vars, file = "TableS2.xlsx")

###Variable elimination
flm1 <- lmer (julian ~ data_a$lat + data_a$d100 + data_a$dd0 + data_a$dd5 + data_a$fday + data_a$mapmtcm + data_a$mat + data_a$mtcm + data_a$mtcmgsp + data_a$mtcmmap + data_a$smrpb + data_a$sprp + data_a$sdimtcm + data_a$dd0gsp + (1 | garden) + (1|ssp:garden) + (1|pop:(ssp:garden)),REML=TRUE , data=fl11)
s_flm1 <- step(flm1) 
s_flm1

####mtcmgsp eliminated because of collinearity with dd0gsp

###FINAL MODELS
flm2 <- lmer (julian ~ data_a$lat+ data_a$d100 + data_a$sdimtcm + data_a$dd0gsp + (1 | garden) + (1|pop:garden),REML=FALSE , data=fl11)

flm3 <- lmer (julian ~ data_a$lat+ data_a$d100 + data_a$sdimtcm +  (1 | garden) + (1|pop:garden),REML=FALSE , data=fl11)

flm5 <- lmer (julian ~ data_a$lat+ data_a$d100 + (1 | garden) + (1|pop:garden),REML=FALSE , data=fl11)

anova(flm2,flm3,flm5)

###BEST MODEL flm5 based on BIC
flm5 <- lmer (julian ~ data_a$lat+ data_a$d100 + (1 | garden) + (1|pop:garden),REML=TRUE , data=fl11)
summary(flm5)
rand(flm5)

###EXAMINE RESIDUALS
residuals <- resid(flm5) 
summary(residuals)
hist(residuals)
plot(residuals)

####R2
r.squaredGLMM(flm5)

####CONFIDENCE INTERVALS
ci_boot <- confint(flm5,level=0.6,method= c("boot"), nsim = 1000)

###EXTRACT FIXED EFFECTS
y.hat4 <- model.matrix(flm5 , type = "fixed") %*% fixef(flm5)

####RANDOM EFFECTS
re_pop <- ranef(flm5, condVar=TRUE, whichel = "pop:garden")
re_garden <- ranef(flm5, condVar=TRUE, whichel = "garden")
dotplot(re_pop)
dotplot(re_garden)

re_pop1 <- unlist(re_pop)
re_pop2 <- as.vector(re_pop1)

###SUMMARY
fit <- with(fl11, data.frame(pop=pop, garden=garden, ssp=ssp, family=family,lat=lat, d100=d100, observed=julian,fitted=fitted(flm5),resid=resid(flm5)))
fit <- cbind(fit,y.hat4)

fit_pop <- summaryBy(observed + y.hat4 + fitted ~ pop + ssp + garden + lat + d100, data= fit, FUN = c(mean))
fit_pop <- cbind(fit_pop,re_pop2)
write.xlsx(x = fit_pop, file = "flowersummary.xlsx")

fit_pop2 <- summaryBy(observed + y.hat4 + fitted ~ pop + ssp + lat + d100, data= fit, FUN = c(mean))
###Remove 2 outlier with one garden represented
fit_pop2 <- fit_pop2[-6:-7,]

###FIXED EFF GRAPH FOR MANUSCRIPT

#q<-ggplot(fit_pop, aes(y=observed.mean,x=fitted.mean))+theme_bw() + ylim(240,320) + xlim(250,290)
#q+stat_smooth(method=lm,se=FALSE,linetype=4,color="gray",size=1) + geom_point(aes(shape=garden),size=3) + xlab("Predicted") + ylab("Observed") + labs(shape = "Gardens") +scale_shape(solid=FALSE) 

l <- ggplot(fit_pop2,aes(y=observed.mean,x=lat)) + theme_bw() + 
  stat_smooth(method=lm,se=FALSE,linetype=4,size=1,color="gray") + geom_point(size=2) + 
  xlab("Latitude") + ylab("Observed") + labs(color = "Gardens") + scale_shape_manual(values=c(21,22,24))

d <- ggplot(fit_pop2,aes(y=observed.mean,x=d100)) + theme_bw() + 
  stat_smooth(method=lm,se=FALSE,linetype=4,size=1,color="gray") + geom_point(size=2)+ xlab("D100") + 
  ylab("Observed") + labs(color = "Gardens") + scale_shape_manual(values=c(21,22,24))

p<-ggplot(fit_pop, aes(y=observed.mean,x=y.hat4.mean,shape=garden)) + theme_bw() + ylim(250,300) + 
  xlim(250,290)+ stat_smooth(method=lm,se=FALSE,linetype=4,size=1,color="gray") + 
  geom_point(size=3) + xlab("Predicted") + ylab("Observed") + labs(color = "Gardens") + scale_shape_manual(values=c(21,22,24))
p

grid.arrange(l,d,nrow=1)

cor.test(fit_pop2$observed.mean,fit_pop2$lat)
cor.test(fit_pop2$observed.mean,fit_pop2$d100)
####TEST FOR VARIATION IN GXE BY LATITUDE
###FIGURE 2
r <- ggplot(fit_pop,aes(lat,re_pop2,shape=garden))+theme_bw()+ scale_shape(solid=FALSE)
r + facet_grid(garden~.) + geom_point() + stat_smooth(method=lm,se=FALSE,linetype=4, color="gray")

cor.test(~lat + re_pop2, data=fit_pop, subset=(garden=="Eph"))
cor.test(~lat + re_pop2, data=fit_pop, subset=(garden=="Maj"))
cor.test(~lat + re_pop2, data=fit_pop, subset=(garden=="Orch"))

#LINEAR MODEL FOR GXE AND LATITUDE AT ORCHARD GARDEN
slopeO <- subset(fit_pop, fit_pop$garden=="Orch")
lm_orch <- lm(re_pop2~lat,data=slopeO)
summary(lm_orch)
slopeO <- cbind(slopeO,predict(lm_orch))

###plot climatypes
p<-ggplot(slopeO, aes(x=lat,y=observed.mean)) +  geom_point()
p+ geom_errorbar(aes(ymax = observed.mean + 13 + predict(lm_orch), ymin = observed.mean - 13 - predict(lm_orch))) + theme_bw()

###Fig S1
###flsum is fit_pop with population removed that are na from one or more gardens
flsum<- read.csv(file="flsum.csv", sep=",",head=TRUE, na.string="na")
flsum$pop <- factor(flsum$pop, levels = flsum$pop[order(flsum$lat)])
zy <- ggplot(flsum, aes(x=pop,y=observed.mean, label=pop))
zy +  geom_point(size=2) + geom_line(aes(group=garden,color=garden)) + theme_bw() + theme(axis.text.x=element_text(angle = 90,size = 7))

###TABLE S1

sumtable <- summaryBy(pop+julian~pop+garden+ssp+lat+long, data=fl11, FUN=c(length))
write.xlsx(x = sumtable, file = "suppl_table1.xlsx")

