#Pilot v Laboratory Efficiency for Medications Development
#Simulation Study of Efficiency

library(SpPack)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gmodels) #CrossTable Function
library(Hmisc) #rcorr
library(compute.es)
library(MBESS)

#parameters to test
#d = c(0.2, 0.5, 0.8) -- small, medium and large canonical Cohen's d
#npilot = seq(6,36,6) -- sample size of pilot test from 6 to 36
#LabMul = seq(2,10,2) -- Lab multiple for sample size associated with laboratory studies being considerably cheaper to run
#r.LabPilot = seq(0.3,0.9,0.2) -- correlation between laboratory-based effect size and the Pilot effect size

parameters <- expand.grid(d = c(0.2, 0.5, 0.8), 
                         npilot = seq(6,36,6),
                         LabMul = seq(2,10,2),
                         r.LabPilot = seq(0.3,0.9,0.2))
parameters <- arrange(parameters, d, npilot, LabMul, r.LabPilot)
parameters$nlab <- parameters$npilot * parameters$LabMul

#Want to run a simulation of 1000 meds for each parameter combination and compute the efficiency of pilot and lab based studies 
#efficiency defined in terms of the number of medications that would receive a positive signal justifying a large Pilot
nSims <- 100
for (i in 1:dim(parameters)[1]){
  #tracking simulations
  print(noquote(paste("parameter combination ", i, "out of ", dim(parameters)[1])))
  
  SimData <- SimCor(n = nSims, 
                    xmean = parameters$d[i], xsd = sqrt(des(d=parameters$d[i], n.1=parameters$npilot[i]/2, n.2=parameters$npilot[i]/2, verbose=F, dig=8)$var.d), 
                    ymean = parameters$d[i], ysd = sqrt(des(d=parameters$d[i], n.1=parameters$nlab[i]/2, n.2=parameters$nlab[i]/2, verbose=F, dig=8)$var.d),
                    rho = parameters$r.LabPilot[i])
  SimData <- rename(SimData, dPilot = x, dLab = y)

  #use lapply to get the sqrt of each dPilot and dLab
  SimData$dPilot.SD <- as.double(lapply(SimData$dPilot, function(x) sqrt(des(d=x, n.1=parameters$npilot[i]/2, n.2=parameters$npilot[i]/2, verbose=F, dig=8)$var.d)))
  SimData$dLab.SD <- as.double(lapply(SimData$dLab, function(x) sqrt(des(d=x, n.1=parameters$nlab[i]/2, n.2=parameters$nlab[i]/2, verbose=F, dig=8)$var.d)))
  
  #get p-values from d's
  SimData$dPilot.p <- as.double(lapply(SimData$dPilot, function(x) des(d=x, n.1=parameters$npilot[i]/2, n.2=parameters$npilot[i]/2, verbose=F, dig=8)$pval.d))
  SimData$dLab.p <- as.double(lapply(SimData$dLab, function(x) des(d=x, n.1=parameters$nlab[i]/2, n.2=parameters$nlab[i]/2, verbose=F, dig=8)$pval.d))

  #save the power of Pilot and associated lab studies
  parameters$EffPilot[i] <- mean(SimData$dPilot.p<0.05)
  parameters$EffLab[i] <- mean(SimData$dLab.p<0.05)
}

View(parameters)

#plotting----
#Dont need to plot the pilot repeated for each LabMul
parameters.pilot <- subset(parameters, LabMul == 2)
Eff.plot <- ggplot(data=parameters, aes(x=npilot, y=EffLab, colour=as.factor(LabMul))) +
  facet_grid(as.factor(r.LabPilot) ~ as.factor(d)) +
  geom_line() + 
  geom_line(data=parameters.pilot, aes(x=npilot, y=EffPilot), colour = "black", size=2) +
  DotRTheme()
Eff.plot


# SpHist(data=SimData, variable = "dPilot")
# SpHist(data=SimData, variable = "dLab")
# SpHist(data=SimData, variable = "dPilot.p")
# SpHist(data=SimData, variable = "dLab.p")


