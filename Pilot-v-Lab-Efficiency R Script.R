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
library(MASS) #mvrnorm simulation function

#parameters to test
#d = c(0.2, 0.5, 0.8) -- small, medium and large canonical Cohen's d
#nPilot = seq(6,36,6) -- sample size of pilot test from 6 to 36
#LabMul = seq(2,10,2) -- Lab multiple for sample size associated with laboratory studies being considerably cheaper to run
#r.LabPilot = seq(0.3,0.9,0.2) -- correlation between laboratory-based effect size and the Pilot effect size

parameters <- expand.grid(nPilot = seq(6,36,6),
                         LabMul = c(1,2,4),
                         r.LabPilot = seq(0.3,0.9,0.3))
parameters <- arrange(parameters, nPilot, LabMul, r.LabPilot)
parameters$nLab <- parameters$nPilot * parameters$LabMul

#testing mvrnorm
test <- data.frame(mvrnorm(n=1000, mu = c(0.8,0.8), Sigma = matrix(c(0.2^2, 0.9*0.2^2, 0.9*0.2^2, 0.2^2), nrow=2)))
SpDesc(test)
cor(test)

#Want to run a simulation of 1000 meds for each parameter combination and compute the efficiency of pilot and lab based studies 
#efficiency defined in terms of the number of medications that would receive a positive signal justifying a large Pilot
nMeds <- 1000
for (i in 1:dim(parameters)[1]){
  #tracking simulations
  print(noquote(paste("parameter combination ", i, "out of ", dim(parameters)[1])))
  
  #Simulate nMeds of each d (0, 0.2, 0.5, 0.8) with SD = 0.2
  #Simulate simultaneously the medication lab efficacy based on parameters$r.LabPilot[i]
  covariances <- matrix(c(0.2^2, parameters$r.LabPilot[i]*0.2^2, parameters$r.LabPilot[i]*0.2^2, 0.2^2), nrow=2)
  SimMeds.d0 <- data.frame(dgroup=0, data.frame(mvrnorm(n=nMeds, mu = c(0,0), Sigma = covariances)))
  SimMeds.d0 <- rename(SimMeds.d0, dPilot = X1, dLab = X2)
  SimMeds.d2 <- data.frame(dgroup=0.2, data.frame(mvrnorm(n=nMeds, mu = c(0.2,0.2), Sigma = covariances)))
  SimMeds.d2 <- rename(SimMeds.d2, dPilot = X1, dLab = X2)
  SimMeds.d5 <- data.frame(dgroup=0.5, data.frame(mvrnorm(n=nMeds, mu = c(0.5,0.5), Sigma = covariances)))
  SimMeds.d5 <- rename(SimMeds.d5, dPilot = X1, dLab = X2)
  SimMeds.d8 <- data.frame(dgroup=0.8, data.frame(mvrnorm(n=nMeds, mu = c(0.8,0.8), Sigma = covariances)))
  SimMeds.d8 <- rename(SimMeds.d8, dPilot = X1, dLab = X2)
  
  #merge Simulated meds
  SimMeds <- rbind(SimMeds.d0, SimMeds.d2, SimMeds.d5, SimMeds.d8)
  #testing
  SimMeds %>% group_by(as.factor(dgroup)) %>% summarise(meandPilot = mean(dPilot), sddPilot = sd(dPilot), 
                                                        meandLab = mean(dLab), sddLab = sd(dLab))
  
  #SimMeds %>% group_by(as.factor(dgroup)) %>% summarise(meandPilot = mean(dPilot), sddPilot = sd(dPilot), meandLab = mean(dLab), sddLab = sd(dLab))  
  
  # #use lapply to compute power for each medication in pilot study and in lab
  #need to use one sided because a negative effect size being detected isn't a positive screen
  SimMeds$power.Pilot <- as.double(lapply(SimMeds$dPilot, function(x) pwr.t.test(n=(parameters$nPilot[i]/2), d=x, sig.level=0.025, 
                                                                                 type="two.sample", alternative="greater")$power))
  SimMeds$power.Lab <- as.double(lapply(SimMeds$dLab, function(x) pwr.t.test(n=(parameters$nLab[i]/2), d=x, sig.level=0.025, 
                                                                                 type="two.sample", alternative="greater")$power))
  
  # SimMeds.d0$power.Pilot <- as.double(lapply(SimMeds.d0$dPilot, function(x) pwr.t.test(n=(parameters$nPilot[i]/2), d=x, sig.level=0.025,
  #                                                                                  type="two.sample", alternative="greater")$power))
  # SimMeds.d0$power.Lab <- as.double(lapply(SimMeds.d0$dLab, function(x) pwr.t.test(n=(parameters$nLab[i]/2), d=x, sig.level=0.025,
  #                                                                            type="two.sample", alternative="greater")$power))
  # SimMeds %>% group_by(as.factor(dgroup)) %>% summarise(meanpowerpilot = mean(power.Pilot), meanpowerlab = mean(power.Lab))
  # ggplot(data=SimMeds.d0, aes(x=dPilot, y=power.Pilot)) + geom_point() + scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(-1,1))
  # ggplot(data=SimMeds.d0, aes(x=dLab, y=power.Lab)) + geom_point() + scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(-1,1))

  #calculate expected number of positive findings by effect size group
  positives <- SimMeds %>% group_by(as.factor(dgroup)) %>% summarise(posPilot = sum(power.Pilot), posLab = sum(power.Lab))  
  
  #calculate the efficiency
  #defined as number of positive screens from pilot/lab studies as a multiple as those with 0 effect size
  #Pilot Efficiencies
  parameters$eff.Pilot.d2[i] <- positives$posPilot[2] / positives$posPilot[1]
  parameters$eff.Pilot.d5[i] <- positives$posPilot[3] / positives$posPilot[1]
  parameters$eff.Pilot.d8[i] <- positives$posPilot[4] / positives$posPilot[1]
  
  #Lab Efficiencies
  parameters$eff.Lab.d2[i] <- positives$posLab[2] / positives$posLab[1]
  parameters$eff.Lab.d5[i] <- positives$posLab[3] / positives$posLab[1]
  parameters$eff.Lab.d8[i] <- positives$posLab[4] / positives$posLab[1]
}

View(parameters)

#plotting----
#average pilot efficiency for each nPilot
PilotEfficiency <- parameters %>% group_by(nPilot) %>% summarise(meaneff.Pilot.d2 = mean(eff.Pilot.d2),
                                                                 meaneff.Pilot.d5 = mean(eff.Pilot.d5),
                                                                 meaneff.Pilot.d8 = mean(eff.Pilot.d8))

colorscale <- scales::seq_gradient_pal("lightblue", "navyblue", "Lab")(seq(0,1,length.out=3))
Eff.plot.d2 <- ggplot(data=parameters, aes(x=nPilot, y=eff.Lab.d2, colour=as.factor(r.LabPilot))) +
  facet_grid(LabMul ~ .) +
  geom_line() + 
  geom_line(data=PilotEfficiency, aes(x=nPilot, y=meaneff.Pilot.d2), colour = "black", size=2) +
  scale_colour_manual("Lab-Pilot\nCorrelation", values=colorscale) +
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", fill=NA))
Eff.plot.d2

colorscale <- scales::seq_gradient_pal("lightblue", "navyblue", "Lab")(seq(0,1,length.out=3))
Eff.plot.d5 <- ggplot(data=parameters, aes(x=nPilot, y=eff.Lab.d5, colour=as.factor(r.LabPilot))) +
  facet_grid(LabMul ~ .) +
  geom_line() + 
  geom_line(data=PilotEfficiency, aes(x=nPilot, y=meaneff.Pilot.d5), colour = "black", size=2) +
  scale_colour_manual("Lab-Pilot\nCorrelation", values=colorscale) +
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", fill=NA))
Eff.plot.d5

colorscale <- scales::seq_gradient_pal("lightblue", "navyblue", "Lab")(seq(0,1,length.out=3))
Eff.plot.d8 <- ggplot(data=parameters, aes(x=nPilot, y=eff.Lab.d8, colour=as.factor(r.LabPilot))) +
  facet_grid(LabMul ~ .) +
  geom_line() + 
  geom_line(data=PilotEfficiency, aes(x=nPilot, y=meaneff.Pilot.d8), colour = "black", size=2) +
  scale_colour_manual("Lab-Pilot\nCorrelation", values=colorscale) +
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", fill=NA))
Eff.plot.d8

# SpHist(data=SimData, variable = "dPilot")
# SpHist(data=SimData, variable = "dLab")
# SpHist(data=SimData, variable = "dPilot.p")
# SpHist(data=SimData, variable = "dLab.p")


