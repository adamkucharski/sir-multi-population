# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# R Code for SIR epidemic in 2 coupled populations
# Author: AJ Kucharski (2017-)
#
# Main code
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Set up source functions
library(deSolve)
library(magrittr)
library(tidyverse)

setwd("~/Documents/GitHub/sir-multi-population/")

source("R/Functions_SIR.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define model parameters

popsize1 <- 10000 # population 1 size
popsize2 <- 5000 # population 2 size
popsize3 <- 9000 # population 2 size
dt <- 1  # time increment for simulation output (days)
time.vals <- seq(0,800,7) # time intervals

r0 <- 1.4 # basic reproduction number
  
# Set parameter vector
theta <- c(beta = NA,
          gamma = 1/10, # mean duration of infectiousness
          alpha1 = 0.0001, # relative contact rate between populations 1-2
          alpha2 = 0.0001, # relative contact rate between populations 2-3
          alpha3 = 0.000, # relative contact rate between populations 3-1
          npop1 = popsize1, # size of population 1
          npop2 = popsize2, # size of population 2
          npop3 = popsize3, # size of population 2
          rep = 0.25 # proportion of cases reported
          )

theta[["beta"]] <- r0*theta[["gamma"]] # define beta as function of R0 and gamma
  
# Set initial conditions - assume everyone susceptible
initial_inf <- 1 # initially infectious in population 1
theta_init <- c(s_init=NA,i1_init=initial_inf,r_init=0,s2_init=NA,i2_init=0,r2_init=0,s3_init=NA,i3_init=0,r3_init=0)
theta_init[["s_init"]] <- theta[["npop1"]]-theta_init[["i1_init"]]
theta_init[["s2_init"]] <- theta[["npop2"]]-theta_init[["i2_init"]]
theta_init[["s3_init"]] <- theta[["npop3"]]-theta_init[["i3_init"]]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run simulation
output_sim <- run_simulation(dt,theta,theta_init,time.vals)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot new reported cases over time
par(mfrow=c(1,1),mar = c(3.2,3.5,1,1),mgp=c(2.2,0.6,0),las=1)

# Total cases
plot(time.vals/7,output_sim$C1_trace,type="l",col='white',xlab="weeks",ylab="total cases",ylim=c(0,500))
lines(time.vals/7, output_sim$C1_trace ,col='blue') # population 1
lines(time.vals/7, output_sim$C2_trace,col='red') # population 2
lines(time.vals/7, output_sim$C3_trace,col='dark green') # population 3
lines(time.vals/7, output_sim$C1_trace + output_sim$C2_trace + output_sim$C3_trace ,col='black',lwd=2)  # total cases

# Reported cases
plot(time.vals/7, output_sim$C1_trace,type="l",col='white',xlab="weeks",ylab="reported cases",ylim=c(0,120))
lines(time.vals/7, cases_reported(0.8*output_sim$C1_trace,theta) ,col='blue') # population 1
lines(time.vals/7, cases_reported(0.8*output_sim$C2_trace,theta),col='red') # population 2
lines(time.vals/7, cases_reported(0.8*output_sim$C3_trace,theta),col='dark green') # population 3
lines(time.vals/7, cases_reported(1.3*( output_sim$C1_trace + output_sim$C2_trace + output_sim$C3_trace),theta) ,col='black',lwd=2)  # total cases


# Save plots
dev.copy(pdf,paste("plots/SIR_model.pdf",sep=""),width=8,height=4)
dev.off()
