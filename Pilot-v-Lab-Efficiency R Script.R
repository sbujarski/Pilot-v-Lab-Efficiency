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
library(pwr)

#parameters to test
#d = c(0.2, 0.5, 0.8) -- small, medium and large canonical Cohen's d
#npilot = seq(6,36,6) -- sample size of pilot test from 6 to 36
#LabMul = seq(2,10,2) -- Lab multiple for sample size associated with laboratory studies being considerably cheaper to run
#r.LabPilot = seq(0.3,0.9,0.2) -- correlation between laboratory-based effect size and the Pilot effect size

parameters <- expand.grid(npilot = seq(6,36,6),
                         LabMul = c(1,2,4),
                         r.LabPilot = seq(0.3,0.9,0.3))
parameters <- arrange(parameters, npilot, LabMul, r.LabPilot)
parameters$nlab <- parameters$npilot * parameters$LabMul

#Want to run a simulation of 1000 meds for each parameter combination and compute the efficiency of pilot and lab based studies 
#efficiency defined in terms of the number of medications that would receive a positive signal justifying a large Pilot
nMeds <- 10
for (i in 1:dim(parameters)[1]){
  #tracking simulations
  print(noquote(paste("parameter combination ", i, "out of ", dim(parameters)[1])))
  
  #Simulate nMeds of each d (0, 0.2, 0.5, 0.8) with SD = 0.2
  #Simulate simultaneously the medication lab efficacy based on parameters$r.LabPilot[i]
  SimMeds.d0 <- data.frame(dgroup=0, SimCor(n = nMeds, xmean=0, xsd=0.2, ymean=0, ysd=0.2, rho=parameters$r.LabPilot[i]))
  SimMeds.d0 <- rename(SimMeds.d0, dPilot = x, dLab = y)
  SimMeds.d2 <- data.frame(dgroup=0.2, SimCor(n = nMeds, xmean=0.2, xsd=0.2, ymean=0.2, ysd=0.2, rho=parameters$r.LabPilot[i]))
  SimMeds.d2 <- rename(SimMeds.d2, dPilot = x, dLab = y)
  SimMeds.d5 <- data.frame(dgroup=0.5, SimCor(n = nMeds, xmean=0.5, xsd=0.2, ymean=0.5, ysd=0.2, rho=parameters$r.LabPilot[i]))
  SimMeds.d5 <- rename(SimMeds.d5, dPilot = x, dLab = y)
  SimMeds.d8 <- data.frame(dgroup=0.8, SimCor(n = nMeds, xmean=0.8, xsd=0.2, ymean=0.8, ysd=0.2, rho=parameters$r.LabPilot[i]))
  SimMeds.d8 <- rename(SimMeds.d8, dPilot = x, dLab = y)
  
  #merge Simulated meds
  SimMeds <- rbind(SimMeds.d0, SimMeds.d2, SimMeds.d5, SimMeds.d8)
  
  #SimMeds %>% group_by(as.factor(dgroup)) %>% summarise(meandPilot = mean(dPilot), sddPilot = sd(dPilot), meandLab = mean(dLab), sddLab = sd(dLab))  
  
  #compute power for each medication in pilot study and in lab
  
  
  
  
  
  # SimData <- SimCor(n = nSims, 
  #                   xmean = parameters$d[i], xsd = sqrt(des(d=parameters$d[i], n.1=parameters$npilot[i]/2, n.2=parameters$npilot[i]/2, verbose=F, dig=8)$var.d), 
  #                   ymean = parameters$d[i], ysd = sqrt(des(d=parameters$d[i], n.1=parameters$nlab[i]/2, n.2=parameters$nlab[i]/2, verbose=F, dig=8)$var.d),
  #                   rho = parameters$r.LabPilot[i])
  # SimData <- rename(SimData, dPilot = x, dLab = y)
  # 
  # #use lapply to get the sqrt of each dPilot and dLab
  # SimData$dPilot.SD <- as.double(lapply(SimData$dPilot, function(x) sqrt(des(d=x, n.1=parameters$npilot[i]/2, n.2=parameters$npilot[i]/2, verbose=F, dig=8)$var.d)))
  # SimData$dLab.SD <- as.double(lapply(SimData$dLab, function(x) sqrt(des(d=x, n.1=parameters$nlab[i]/2, n.2=parameters$nlab[i]/2, verbose=F, dig=8)$var.d)))
  # 
  # #get p-values from d's
  # SimData$dPilot.p <- as.double(lapply(SimData$dPilot, function(x) des(d=x, n.1=parameters$npilot[i]/2, n.2=parameters$npilot[i]/2, verbose=F, dig=8)$pval.d))
  # SimData$dLab.p <- as.double(lapply(SimData$dLab, function(x) des(d=x, n.1=parameters$nlab[i]/2, n.2=parameters$nlab[i]/2, verbose=F, dig=8)$pval.d))
  # 
  # #save the power of Pilot and associated lab studies
  # parameters$EffPilot[i] <- mean(SimData$dPilot.p<0.05)
  # 
  # #How to account for the fact that a significant lab result might be a false positive
  # #Only count a significant Lab test as "efficient positive screen" if the SimData$dPilot which is the medication RCT effect size is >= 0.2
  # parameters$EffLab[i] <- sum(subset(SimData, dPilot >= 0.2)$dLab.p<0.05)/nSims
}

View(parameters)

#plotting----
#Dont need to plot the pilot repeated for each LabMul
parameters.pilot <- subset(parameters, LabMul == 2)

colorscale <- scales::seq_gradient_pal("lightblue", "navyblue", "Lab")(seq(0,1,length.out=4))
Eff.plot <- ggplot(data=parameters, aes(x=npilot, y=EffLab, colour=as.factor(LabMul))) +
  facet_grid(paste("Lab-RCT r = ", r.LabPilot) ~ paste("Cohen's d = ",d)) +
  geom_line() + 
  scale_colour_manual("Lab Sample\nMultiple", values=colorscale) +
  geom_line(data=parameters.pilot, aes(x=npilot, y=EffPilot), colour = "black", size=2) +
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", fill=NA))
Eff.plot


# SpHist(data=SimData, variable = "dPilot")
# SpHist(data=SimData, variable = "dLab")
# SpHist(data=SimData, variable = "dPilot.p")
# SpHist(data=SimData, variable = "dLab.p")


