# Statistical analyses for Dalziel et al. 2018 Science
#
#

rm(list=ls())
flu <- read.csv("influenza_cities.csv")



keep <- complete.cases(flu) 
flu <- flu[keep,]
rm(keep)



# Fit crowding as a function of population size
hcrowd_pop <- glm(log(h_crowding) ~ log(population_size), data = flu)  #residential crowding
wcrowd_pop <- glm(log(w_crowding) ~ log(population_size), data = flu)  #daytime crowding



# Fit base transmission potential (kappa) as a function of population size and crowding 
kappa_pop <- glm(log(kappa) ~ log(population_size), data = flu)
kappa_popXcrowd <- glm(log(kappa) ~ log(population_size)*log(h_crowding)*log(w_crowding), data = flu)



# Define predicted kappas from census data using both population size and crowding
# to use when we want a measure of kappa that derives from census data, not incidence data
flu$kpred <- kappa_popXcrowd$fitted.values



# Define excess kappa as the residual variation in base transmission not explained by population size
# and excess crowding as the residual variation in crowding not explained by population size
flu$exkappa <- kappa_pop$residuals
flu$exhcrowd <- hcrowd_pop$residuals
flu$exwcrowd <- wcrowd_pop$residuals



# Fit excess kappa as a function of excess crowding
exkappa_exhcrowd <- glm(exkappa ~ exhcrowd, data = flu)
exkappa_exwcrowd <- glm(exkappa ~ exwcrowd, data = flu)



# Fit epidemic intensity to population size, 
# transmission potential (estimated independently from census data), 
# and specific humidity
nu_pop <- glm(log(nu) ~ log(population_size), data = flu)
nu_sh <- glm(log(nu) ~ avgsh * sdsh + I(avgsh^2) + I(sdsh^2), data = flu)
nu_kappa <- glm(log(nu) ~ log(kpred), data = flu)
nu_popXsh <- glm(log(nu) ~ log(population_size) * avgsh * sdsh + I(avgsh^2) + I(sdsh^2), data = flu)
nu_kappaXsh <- glm(log(nu) ~ log(kpred) * avgsh * sdsh + I(avgsh^2) + I(sdsh^2), data = flu)




# Do model comparison for the intensity GLMs 
modcomp <- AIC(nu_pop,nu_sh,nu_kappa,nu_popXsh,nu_kappaXsh)

modcomp$Model <- c("Population size (P)", "Specific humidity (SH)", "Base transmission (K)", "P x SH", "K x SH")
modcomp <- modcomp[,c(3,1,2)]
rownames(modcomp) <- NULL

modcomp <- modcomp[order(modcomp$AIC),]
modcomp$DeltaAIC <- modcomp$AIC - modcomp$AIC[1]



# Bar plot of deltaAICs for intensity models
quartz(h=4,w=4)
pal <- 1 #brewer.pal(5,"Set1")
x <- barplot(modcomp$DeltaAIC,xlab = expression(paste(Delta,"AIC")),ylab="Model",srt=45,horiz=T,col=pal)
text(rep(30,5),x,modcomp$Model,pos=4,col='white', cex=0.95)
text(35,x[2],modcomp$Model[2],pos=4,font=1,cex=0.95)
text(35,x[1],modcomp$Model[1],pos=4,font=2,cex=0.95)

