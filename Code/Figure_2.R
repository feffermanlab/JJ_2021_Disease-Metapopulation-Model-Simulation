## =========================================================================================================================================================
##
## R code to generate Figure 3 and 4 for manuscript titled "The Eco-Epidemiological Dynamics of Evolutionary Rescue in a Host Metapopulation"

## =========================================================================================================================================================


library(ggplot2)
library(reshape)
library(here)

timestep <- 5000

## set the current directory as default to save all the following results
#set_here()
setwd("C:\\Users\\lucyj\\Documents\\New folder\\Files for Avian Malaria Project\\SIR model\\Code")

onepatchdynamics <- read.csv("SIRmodel_one_patch_dynamics_TEST3-14-21.csv",header=FALSE)
wildsize <- onepatchdynamics[,1] + onepatchdynamics[,2] + onepatchdynamics[,3]
robustsize <- onepatchdynamics[,4] + onepatchdynamics[,5] + onepatchdynamics[,6]
totalinfec <- onepatchdynamics[,2] + onepatchdynamics[,5]

## trans level 1 and growth level 1
wildsize11 <- wildsize[1:timestep]
robustsize11 <- robustsize[1:timestep]
totalinfec11 <- totalinfec[1:timestep]

# trans level 1 and growth level 2
wildsize12 <- wildsize[(timestep + 1): (timestep*2)]
robustsize12 <- robustsize[(timestep + 1): (timestep*2)]
totalinfec12 <- totalinfec[(timestep + 1): (timestep*2)]

# trans level 2 and growth level 1

wildsize21 <- wildsize[(2*timestep + 1): (timestep*3)]
robustsize21 <- robustsize[(2*timestep + 1): (timestep*3)]
totalinfec21 <- totalinfec[(2*timestep + 1): (timestep*3)]

# trans level 2 and growth level 2
wildsize22 <- wildsize[(3*timestep + 1): (timestep*4)]
robustsize22 <- robustsize[(3*timestep + 1): (timestep*4)]
totalinfec22 <- totalinfec[(3*timestep + 1): (timestep*4)]

# trans level 3 and growth level 1

wildsize31 <- wildsize[(4*timestep + 1): (timestep*5)]
robustsize31 <- robustsize[(4*timestep + 1): (timestep*5)]
totalinfec31 <- totalinfec[(4*timestep + 1): (timestep*5)]

# trans level 3 and growth level 2
wildsize32 <- wildsize[(5*timestep + 1): (timestep*6)]
robustsize32 <- robustsize[(5*timestep + 1): (timestep*6)]
totalinfec32 <- totalinfec[(5*timestep + 1): (timestep*6)]

x_axis <- seq(1, timestep, 1)

ylim_wild <- c(min(wildsize11, wildsize12, wildsize21, wildsize22, wildsize31, wildsize32), (max(wildsize11, wildsize12, wildsize21, wildsize22, wildsize31, wildsize32) + 500))

ylim_robust <- c(min(robustsize11, robustsize12, robustsize21, robustsize22, robustsize31, robustsize32), max(robustsize11, robustsize12, robustsize21, robustsize22, robustsize31, robustsize32))

ylim_totinfec <- c(min(totalinfec11, totalinfec12, totalinfec21, totalinfec22, totalinfec31, totalinfec32), max(totalinfec11, totalinfec12, totalinfec21, totalinfec22, totalinfec31, totalinfec32))
  

theme_set(theme_bw(20))

tiff("Figure_2-test_0.0005-3-14-21.tiff", width = 10, height = 10, units = 'in', res = 600)

par(mfrow = c(3, 2))

plot(x_axis, robustsize11, type = "l", lty = 1, ylim = c(0, 2000), xlab = '', ylab = '', lwd = 2, xlim=c(0,3000))
points(x_axis, robustsize21, type = "l", lty = 3, lwd = 2)
points(x_axis, robustsize31, type = "l", lty = 2, lwd = 2)


plot(x_axis, robustsize12, type = "l", lty = 1, ylim = c(0, 2000), xlab = '', ylab = '', lwd = 2,xlim=c(0,3000))
points(x_axis, robustsize22, type = "l", lty = 3, lwd = 2)
points(x_axis, robustsize32, type = "l", lty = 2, lwd = 2)

plot(x_axis, wildsize11, type = "l", lty = 1, ylim =  c(0, 2200), xlab = '', ylab = '', lwd = 2,xlim=c(0,3000))
points(x_axis, wildsize21, type = "l", lty = 3, lwd = 2)
points(x_axis, wildsize31, type = "l", lty = 2, lwd = 2)


plot(x_axis, wildsize12, type = "l", lty = 1, ylim =  c(0, 2200), xlab = '', ylab = '', lwd = 2, xlim=c(0,3000))
points(x_axis, wildsize22, type = "l", lty = 3, lwd = 2)
points(x_axis, wildsize32, type = "l", lty = 2, lwd = 2)

plot(x_axis, totalinfec11, type = "l", lty = 1, ylim = ylim_totinfec, xlab = '', ylab = '', lwd = 2,xlim=c(0,3000))
points(x_axis, totalinfec21, type = "l", lty = 3, lwd = 2)
points(x_axis, totalinfec31, type = "l", lty = 2, lwd = 2)


plot(x_axis, totalinfec12, type = "l", lty = 1, ylim = ylim_totinfec, xlab = '', ylab = '', lwd = 2,xlim=c(0,3000))
points(x_axis, totalinfec22, type = "l", lty = 3, lwd = 2)
points(x_axis, totalinfec32, type = "l", lty = 2, lwd = 2)

dev.off()
