## =========================================================================================================================================================
##
## R code to generate Figure 2 Appendix for manuscript titled "The Eco-Epidemiological Dynamics of Evolutionary Rescue in a Host Metapopulation"

## =========================================================================================================================================================


library(ggplot2)
library(reshape)
library(here)

timestep <- 5000

## set the current directory as default to save all the following results
#set_here()
setwd("C:\\Users\\lucyj\\Documents\\New folder\\Files for Avian Malaria Project\\SIR model\\Code")

onepatchdynamics_0.0005 <- read.csv("SIRmodel_one_patch_dynamics_TEST0.0005.csv",header=FALSE)
wildsize_0.0005 <- onepatchdynamics_0.0005[,1] + onepatchdynamics_0.0005[,2] + onepatchdynamics_0.0005[,3]
robustsize_0.0005 <- onepatchdynamics_0.0005[,4] + onepatchdynamics_0.0005[,5] + onepatchdynamics_0.0005[,6]
totalinfec_0.0005 <- onepatchdynamics_0.0005[,2] + onepatchdynamics_0.0005[,5]


onepatchdynamics_0.0001 <- read.csv("SIRmodel_one_patch_dynamics_TEST0.0001.csv",header=FALSE)
wildsize_0.0001 <- onepatchdynamics_0.0001[,1] + onepatchdynamics_0.0001[,2] + onepatchdynamics_0.0001[,3]
robustsize_0.0001 <- onepatchdynamics_0.0001[,4] + onepatchdynamics_0.0001[,5] + onepatchdynamics_0.0001[,6]
totalinfec_0.0001 <- onepatchdynamics_0.0001[,2] + onepatchdynamics_0.0001[,5]


onepatchdynamics_0.0009 <- read.csv("SIRmodel_one_patch_dynamics_TEST0.0009.csv",header=FALSE)
wildsize_0.0009 <- onepatchdynamics_0.0009[,1] + onepatchdynamics_0.0009[,2] + onepatchdynamics_0.0009[,3]
robustsize_0.0009 <- onepatchdynamics_0.0009[,4] + onepatchdynamics_0.0009[,5] + onepatchdynamics_0.0009[,6]
totalinfec_0.0009 <- onepatchdynamics_0.0009[,2] + onepatchdynamics_0.0009[,5]

# trans level 3 and growth level 1

totalinfec31_0.0005 <- totalinfec_0.0005[(4*timestep + 1): (timestep*5)]

totalinfec31_0.0001 <- totalinfec_0.0001[(4*timestep + 1): (timestep*5)]

totalinfec31_0.0009 <- totalinfec_0.0009[(4*timestep + 1): (timestep*5)]

x_axis <- seq(1, timestep, 1)


ylim_totinfec <- c(min(totalinfec31_0.0005,totalinfec31_0.0001,totalinfec31_0.0009), max(totalinfec31_0.0005,totalinfec31_0.0001,totalinfec31_0.0009))


theme_set(theme_bw(20))

tiff("Figure_2_Appendix.tiff", width = 10, height = 10, units = 'in', res = 600)

#par(mfrow = c(3, 2))

plot(x_axis, totalinfec31_0.0005, type = "l", lty = 1, xlim = c(0, 500), ylim = ylim_totinfec, xlab = '', ylab = '', lwd = 2)
points(x_axis, totalinfec31_0.0001, type = "l", lty = 3, lwd = 2)
points(x_axis, totalinfec31_0.0009, type = "l", lty = 2, lwd = 2)

legend(350, 250, c("Intial_Infecteds = 0.0001", '0.0005', "0.0009"), lty = c(3,1,2), bty ="n", lwd= c(2,2,2))

dev.off()
