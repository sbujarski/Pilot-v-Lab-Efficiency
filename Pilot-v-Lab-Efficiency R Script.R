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
#r.LabRCT = seq(0.3,0.9,0.2) -- correlation between laboratory-based effect size and the RCT effect size

parameters <- expand.grid(d = c(0.2, 0.5, 0.8), 
                         npilot = seq(6,36,6),
                         LabMul = seq(2,10,2),
                         r.LabRCT = seq(0.3,0.9,0.2))
parameters <- arrange(parameters, d, npilot, LabMul, r.LabRCT)

#compute SD of d
for (i in 1:dim(parameters)[1]){
  parameters$d.SD <- sqrt(des(d=parameters$d[i], n.1=parameters$npilot/2, n.2=parameters$npilot/2, verbose=F)$var.d)
}

#Want to run a simulation of 1000 meds for each parameter combination and compute the efficiency of pilot and lab based studies 
#efficiency defined in terms of the number of medications that would receive a positive signal justifying a large RCT

#testing with i = 1
nSims <- 1000
SimData <- SimCor(n = nSims, 
                  xmean = parameters$d[i], xsd = parameters$d.SD[i], 
                  ymean = parameters$d[i], ysd = parameters$d.SD[i],
                  rho = parameters$r.LabRCT[i])
SimData <- rename(SimData, dRCT = x, dLab = y)
# SpDesc(SimData)
# cor.test(~ dRCT + dLab, data=SimData)
#use lapply to get the sqrt of each dRCT and dLab
SimData$dRCT.SD <- lapply(SimData$dRCT, function(x) sqrt(des(d=x, n.1=parameters$npilot[i]/2, n.2=parameters$npilot[i]/2, verbose=F)$var.d))
SimData$dLab.SD <- lapply(SimData$dLab, function(x) sqrt(des(d=x, n.1=parameters$npilot[i]/2, n.2=parameters$npilot[i]/2, verbose=F)$var.d))


