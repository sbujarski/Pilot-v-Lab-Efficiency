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
#d = c(0.2, 0.5, 0.8) -- medication effect size
#nPilot = seq(6, 36, 6) -- sample size of pilot test from 6 to 36
#LabMul = seq(1, 2, 4) -- Lab multiple for sample size associated with laboratory studies being considerably cheaper to run
#r.LabPilot = seq(0.3, 0.9, 0.3) -- correlation between laboratory-based effect size and the Pilot effect size

parameters <- expand.grid(d = c(0.2, 0.5, 0.8),
                          nPilot = seq(6, 36, 6),
                          LabMul = c(1, 2, 4),
                          r.LabPilot = seq(0.3,0.9,0.3))
parameters <- arrange(parameters, nPilot, LabMul, r.LabPilot)
parameters$nLab <- parameters$nPilot * parameters$LabMul
dim(parameters)

#testing mvrnorm
# test <- data.frame(mvrnorm(n=1000, mu = c(0.8,0.8), Sigma = matrix(c(0.2^2, 0.9*0.2^2, 0.9*0.2^2, 0.2^2), nrow=2)))
# SpDesc(test)
# cor(test)

#testing
# SpDesc(SimMeds)
# cor(SimMeds)
# SimMeds %>% group_by(as.factor(dgroup)) %>% summarise(meandPilot = mean(dPilot), sddPilot = sd(dPilot), meandLab = mean(dLab), sddLab = sd(dLab))
# #simulation testing the rescaling approach to counter the upward bias of lab SD
# meanPilot0 <- rep(NA,1000)
# sdPilot0 <- rep(NA,1000)
# meanLab0<- rep(NA,1000)
# sdLab0 <- rep(NA,1000)
# meanPilot2 <- rep(NA,1000)
# sdPilot2 <- rep(NA,1000)
# meanLab2<- rep(NA,1000)
# sdLab2 <- rep(NA,1000)
# r <- rep(NA,1000)
# for(k in 1:1000){
#   #simulate nMeds of 0 effect size and nMeds of parameters$d[i] effect size with SD = 0.2
#   SimPilot <- c(rnorm(n=nMeds, mean=0, sd=0.2), rnorm(n=nMeds, mean=0.2, sd=0.2))
#   #Simulate each medications medication lab effect size based on parameters$r.LabPilot[i]
#   SimLab <- SimCorX(x=SimPilot, ymean=mean(SimPilot), ysd=sd(SimPilot), rho=0.9)$y
#   
#   #SimCorX results in biased up SD, Attempting to renormalize
#   SimLab <- scale(SimLab)
#   SimLab <- SimLab * rnorm(1, sd(SimPilot), 0.02) + mean(SimPilot)
#   
#   #Combine into a SimMeds dataset
#   SimMeds <- data.frame(dgroup=c(rep(0,nMeds),rep(0.2,nMeds)),
#                         dPilot=SimPilot, 
#                         dLab=SimLab)
#   #testing
#   r[k] <- cor(SimMeds)[3,2]
#   Desc<- SimMeds %>% group_by(as.factor(dgroup)) %>% summarise(meandPilot = mean(dPilot), sddPilot = sd(dPilot), 
#                                                         meandLab = mean(dLab), sddLab = sd(dLab))
#   meanPilot0[k] <- Desc$meandPilot[1]
#   sdPilot0[k] <- Desc$sddPilot[1]
#   meanLab0[k] <- Desc$meandLab[1]
#   sdLab0[k] <- Desc$sddLab[1]
#   meanPilot2[k] <- Desc$meandPilot[2]
#   sdPilot2[k] <- Desc$sddPilot[2]
#   meanLab2[k] <- Desc$meandLab[2]
#   sdLab2[k] <- Desc$sddLab[2]
# }
# SpHist(r)
# SpHist(meanPilot0)
# SpHist(meanLab0)
# SpHist(sdPilot0)
# SpHist(sdLab0)
# SpHist(meanPilot2)
# SpHist(meanLab2)
# SpHist(sdPilot2)
# SpHist(sdLab2)
#SimMeds %>% group_by(as.factor(dgroup)) %>% summarise(meandPilot = mean(dPilot), sddPilot = sd(dPilot), meandLab = mean(dLab), sddLab = sd(dLab))  


#Want to run a simulation of 1000 meds for each parameter combination and compute the efficiency of pilot and lab based studies 
#efficiency defined in terms of the number of medications that would receive a positive signal justifying a large Pilot
nMeds <- 10
for (i in 1:dim(parameters)[1]){
  #tracking simulations
  print(noquote(paste("parameter combination ", i, "out of ", dim(parameters)[1])))
  
  #simulate nMeds of 0 effect size and nMeds of parameters$d[i] effect size with SD = 0.2
  SimPilot <- c(rnorm(n=nMeds, mean=0, sd=0.2), rnorm(n=nMeds, mean=parameters$d[i], sd=0.2))
  #Simulate each medications medication lab effect size based on parameters$r.LabPilot[i]
  SimLab <- SimCorX(x=SimPilot, ymean=mean(SimPilot), ysd=sd(SimPilot), rho=parameters$r.LabPilot[i])$y

  #SimCorX results in biased up SD, renormaling to not have dramatically upward biased sd
  meanSimLab <- mean(SimLab)
  SimLab <- scale(SimLab)
  SimLab <- SimLab * rnorm(1, sd(SimPilot), 0.02) + meanSimLab #0.02 based on simulation testing of sd of sd
  
  #Combine into a SimMeds dataset
  SimMeds <- data.frame(dgroup=c(rep(0,nMeds),rep(parameters$d[i],nMeds)),
                        dPilot=SimPilot, 
                        dLab=SimLab)
  
  # #use lapply to compute power for each medication in pilot study and in lab
  #need to use one sided because a negative effect size being detected isn't a positive screen
  SimMeds$power.Pilot <- as.double(lapply(SimMeds$dPilot, function(x) pwr.t.test(n=(parameters$nPilot[i]/2), d=x, sig.level=0.025, 
                                                                                 type="two.sample", alternative="greater")$power))
  SimMeds$power.Lab <- as.double(lapply(SimMeds$dLab, function(x) pwr.t.test(n=(parameters$nLab[i]/2), d=x, sig.level=0.025, 
                                                                                 type="two.sample", alternative="greater")$power))
  
  #calculate expected number of positive findings by effect size group
  positives <- SimMeds %>% group_by(as.factor(dgroup)) %>% summarise(posPilot = sum(power.Pilot), posLab = sum(power.Lab))  
  
  #What do we want to know?
  #1) what is the probability that an effective medication will screen positive (sensitivity)
  parameters$sens.Pilot[i] <- positives$posPilot[2] / nMeds
  parameters$sens.Lab[i] <- positives$posLab[2] / nMeds
  
  #2) What is the probability that an ineffective medication screens negative (specificity) - independent of medication effect size
  parameters$spec.Pilot[i] <- (nMeds-positives$posPilot[1]) / nMeds 
  parameters$spec.Lab[i] <- (nMeds-positives$posLab[1]) / nMeds 
  
  #3) what is the probability that positive screen is a true positive (positive predictive value)
  parameters$ppv.Pilot[i] <- positives$posPilot[2] / (positives$posPilot[2] + positives$posPilot[1]) 
  parameters$ppv.Lab[i] <- positives$posLab[2] / (positives$posLab[2] + positives$posLab[1]) 

  #4) what is the probability that negative screen is a true negative (negative predictive value)
  parameters$npv.Pilot[i] <- (nMeds - positives$posPilot[1]) / (2*nMeds - (positives$posPilot[2] + positives$posPilot[1])) 
  parameters$npv.Lab[i] <- (nMeds - positives$posLab[1]) / (2*nMeds - (positives$posLab[2] + positives$posLab[1])) 
}

View(parameters)

#plotting----
#average pilot efficiency for each nPilot
PilotEfficiency <- parameters %>% group_by(nPilot) %>% summarise(sens.Pilot = mean(sens.Pilot),
                                                                 spec.Pilot = mean(spec.Pilot),
                                                                 ppv.Pilot = mean(ppv.Pilot),
                                                                 npv.Pilot = mean(npv.Pilot))
#make strings for facet labels
parameters$dstr <- paste("d = ", parameters$d)
parameters$LabMulstr <- paste("Lab Multiple = ", parameters$LabMul)


#Sensitivity
colorscale <- scales::seq_gradient_pal("lightblue", "navyblue", "Lab")(seq(0,1,length.out=3))
sens.plot <- ggplot(data=parameters, aes(x=nPilot, y=sens.Lab, colour=as.factor(r.LabPilot))) +
  facet_grid(LabMulstr ~ dstr) +
  geom_line(size=1.5) + 
  geom_line(data=PilotEfficiency, aes(x=nPilot, y=sens.Pilot), colour = "black", size=2) + 
  scale_colour_manual("Correlation between\nLab and Clinic\nEffects", values=colorscale) +
  scale_x_continuous("Pilot Sample Size", limits=c(6,36), breaks=seq(6,36,6)) +
  scale_y_continuous("Sensitivity", limits=c(0,1)) +
  ggtitle("Sensitivity - Pilot versus Lab Screening") +
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", fill=NA))
sens.plot

#Specificity
colorscale <- scales::seq_gradient_pal("lightblue", "navyblue", "Lab")(seq(0,1,length.out=3))
spec.plot <- ggplot(data=parameters, aes(x=nPilot, y=spec.Lab, colour=as.factor(r.LabPilot))) +
  facet_grid(LabMulstr ~ dstr) +
  geom_line(size=1.5) + 
  geom_line(data=PilotEfficiency, aes(x=nPilot, y=spec.Pilot), colour = "black", size=2) + 
  scale_colour_manual("Correlation between\nLab and Clinic\nEffects", values=colorscale) +
  scale_x_continuous("Pilot Sample Size", limits=c(6,36), breaks=seq(6,36,6)) +
  scale_y_continuous("Specificity", limits=c(0,1)) +
  ggtitle("Specificity - Pilot versus Lab Screening") +
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", fill=NA))
spec.plot

#ppv
colorscale <- scales::seq_gradient_pal("lightblue", "navyblue", "Lab")(seq(0,1,length.out=3))
ppv.plot <- ggplot(data=parameters, aes(x=nPilot, y=ppv.Lab, colour=as.factor(r.LabPilot))) +
  facet_grid(LabMulstr ~ dstr) +
  geom_line(size=1.5) + 
  geom_line(data=PilotEfficiency, aes(x=nPilot, y=ppv.Pilot), colour = "black", size=2) + 
  scale_colour_manual("Correlation between\nLab and Clinic\nEffects", values=colorscale) +
  scale_x_continuous("Pilot Sample Size", limits=c(6,36), breaks=seq(6,36,6)) +
  scale_y_continuous("Positive Predictive Value", limits=c(0,1)) +
  ggtitle("Positive Predictive Value - Pilot versus Lab Screening") +
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", fill=NA))
ppv.plot


#npv
colorscale <- scales::seq_gradient_pal("lightblue", "navyblue", "Lab")(seq(0,1,length.out=3))
npv.plot <- ggplot(data=parameters, aes(x=nPilot, y=npv.Lab, colour=as.factor(r.LabPilot))) +
  facet_grid(LabMulstr ~ dstr) +
  geom_line(size=1.5) + 
  geom_line(data=PilotEfficiency, aes(x=nPilot, y=npv.Pilot), colour = "black", size=2) + 
  scale_colour_manual("Correlation between\nLab and Clinic\nEffects", values=colorscale) +
  scale_x_continuous("Pilot Sample Size", limits=c(6,36), breaks=seq(6,36,6)) +
  scale_y_continuous("Negative Predictive Value", limits=c(0,1)) +
  ggtitle("Negative Predictive Value - Pilot versus Lab Screening") +
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", fill=NA))
npv.plot




